#!/usr/bin/env python
# coding: utf-8

import pysam
import os
import pandas as pd
import numpy as np
import time
import argparse
import sys
from multiprocessing import Pool
import yaml
import pybedtools

map_bed = sys.argv[1]
samples = os.path.basename(map_bed).split('.')[0]
out_dir = sys.argv[2]
src_path=os.path.dirname(sys.argv[0])
tmp_dir = f"{out_dir}/tmp_{samples}/"

GC_counts_file = f"{out_dir}/GC_counts/{samples}.GC_counts.txt"
#params
sample_name = f"{samples}"
mappable_regions_path = f"{src_path}/../Ref/k100_minus_exclusion_lists.mappable_regions.hg38.bed" ##config
ref_seq = f"{src_path}/../Ref/hg38.fa" ##config
chrom_sizes = f"{src_path}/../Ref/hg38.standard.chrom.sizes" ##config
map_q = 30 ##config
size_range = [150, 210] ##config
CPU = 5


## rename 
bed_file         = map_bed # bed_file 
bed_file_name    = sample_name
ref_seq_path     = ref_seq
chrom_sizes_path = chrom_sizes
size_range       = size_range
out_file         = out_dir +'/GC_counts/'+ bed_file_name+'.GC_counts.txt'

print('out_file',out_file)
#create a directory for the GC data
for fd in [out_dir, tmp_dir, out_dir +'/GC_counts/']:
    if not os.path.exists(fd):
        os.makedirs(fd)

pybedtools.set_tempdir(tmp_dir)

###########################################输出配置
print('arguments provided:')

print('\tbed_file_name = "'+bed_file_name+'"')
print('\tmappable_regions_path = "'+mappable_regions_path+'"')

print('\tref_seq_path = "'+ref_seq_path+'"')
print('\tchrom_sizes_path = "'+chrom_sizes_path+'"')
print('\tout_dir = "'+out_dir+'"')

print('\tmap_q = '+str(map_q))
print('\tsize_range = '+str(size_range))
print('\tCPU = '+str(CPU))


###########################################筛选常染色体的region
#import filter
mappable_intervals = pd.read_csv(mappable_regions_path, sep='\t', header=None)

#remove non standard chromosomes and X and Y
chroms = ['chr'+str(m) for m in range(1,23)]
mappable_intervals = mappable_intervals[mappable_intervals[0].isin(chroms)]

print('chroms:', chroms)
print('number_of_intervals:',len(mappable_intervals))

sys.stdout.flush()

###################################################################################
def collect_reads(sub_list_file):
    #create a dict for holding the frequency of each read length and GC content
    GC_dict = {}
    for length in range(size_range[0],size_range[1]+1):
        GC_dict[length]={}
        for num_GC in range(0,length+1):
            GC_dict[length][num_GC]=0
    #import the bam file
    #this needs to be done within the loop otherwise it gives a truncated file warning
    bed_file_sub = pybedtools.BedTool(bed_file).intersect(pybedtools.BedTool(sub_list_file), wa = True)
    #this might also need to be in the loop
    #import the ref_seq
    ref_seq=pysam.FastaFile(ref_seq_path)
    n=0
    #fetch any read
    for read in bed_file_sub:
        try:
            fragment_mapq = int(read.name)
        except:
            fragment_mapq = 30
        fragment_chr    = read.chrom
        fragment_start  = read.start
        fragment_end    = read.end
        fragment_length = read.length
        if fragment_chr in chroms and fragment_mapq >= 30 and fragment_length >= size_range[0] and fragment_length <= size_range[1] :
            #count the GC content
            fragment_seq = ref_seq.fetch(fragment_chr,fragment_start,fragment_end)
            fragment_seq = np.array(list(fragment_seq.upper()))
            fragment_seq[np.isin(fragment_seq, ['A','T','W'])] = 0
            fragment_seq[np.isin(fragment_seq, ['C','G','S'])] = 1
            rng = np.random.default_rng(fragment_start)
            fragment_seq[np.isin(fragment_seq, ['N','R','Y','K','M','B','D','H','V'])] = rng.integers(2, size=len(fragment_seq[np.isin(fragment_seq, ['N','R','Y','K','M','B','D','H','V'])])) #random integer in range(2) (i.e. 0 or 1)
            fragment_seq = fragment_seq.astype(float)
            num_GC = int(fragment_seq.sum())
            GC_dict[fragment_length][num_GC]+=1
        n=n+1
        if n%200000==0:
            print('frgmt ',n,':',fragment_chr,fragment_start,fragment_end,'seconds:',np.round(time.time()-start_time))
            sys.stdout.flush()
    print('done')
    return(GC_dict)


###########################
########## main ###########
###########################
start_time = time.time()
p = Pool(processes=CPU) #use the available CPU
sublists = np.array_split(mappable_intervals,CPU) #split the list into sublists, one per CPU
sub_list_files = []
for sublist_n in range(len(sublists)):
    sublists[sublist_n].to_csv(f'{tmp_dir}/map_region.{sublist_n}.bed',sep='\t',index=False,header=None)
    sub_list_files.append(f'{tmp_dir}/map_region.{sublist_n}.bed')


GC_dict_list = p.map(collect_reads, sub_list_files, 1)


all_GC_df = pd.DataFrame()
for i,GC_dict in enumerate(GC_dict_list):
    GC_df = pd.DataFrame()
    for length in GC_dict.keys():
        current = pd.Series(GC_dict[length]).reset_index()
        current = current.rename(columns={'index':'num_GC',0:'number_of_fragments'})
        current['length']=length
        current = current[['length','num_GC','number_of_fragments']]
        GC_df = pd.concat([GC_df, current], ignore_index=True)
    GC_df = GC_df.set_index(['length','num_GC'])
    all_GC_df[i] = GC_df['number_of_fragments']
    del(GC_df,GC_dict)

all_GC_df = all_GC_df.sum(axis=1)
all_GC_df = pd.DataFrame(all_GC_df).rename(columns = {0:'number_of_fragments'})
all_GC_df = all_GC_df.reset_index()
all_GC_df.to_csv(out_file,sep='\t',index=False)
print('done')
