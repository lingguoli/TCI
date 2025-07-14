import os,sys
import pysam
import pandas as pd
import numpy as np
import time


bp_bed, GC_bias_wgs, out_dir = sys.argv[1:]
sample = os.path.basename(bp_bed).split('.')[0]
src_path=os.path.dirname(sys.argv[0])
TPM = f'{src_path}/../Ref/TPM_atlas.csv'
ref_seq = f"{src_path}/../Ref/hg38.fa"
fa = pysam.FastaFile(ref_seq)
size_range = [150,210]

def GC(fragment_chr,fragment_start,fragment_end):
    #count the GC content
    fragment_seq = fa.fetch(fragment_chr,fragment_start,fragment_end)
    fragment_seq = np.array(list(fragment_seq.upper()))
    fragment_seq[np.isin(fragment_seq, ['A','T','W'])] = 0
    fragment_seq[np.isin(fragment_seq, ['C','G','S'])] = 1
    rng = np.random.default_rng(fragment_start)
    fragment_seq[np.isin(fragment_seq, ['N','R','Y','K','M','B','D','H','V'])] = rng.integers(2, size=len(fragment_seq[np.isin(fragment_seq, ['N','R','Y','K','M','B','D','H','V'])])) #random integer in range(2) (i.e. 0 or 1)
    fragment_seq = fragment_seq.astype(float)
    num_GC = int(fragment_seq.sum())
    return num_GC

df_bias_wgs = pd.read_table(GC_bias_wgs)
dict_bias_wgs = df_bias_wgs.set_index(['length','num_GC']).to_dict()

df_bed = pd.read_table(bp_bed,names=['chr','start','end', 'TSS', 'bps'])
df_bed['length'] = df_bed.apply(lambda x: x['end'] - x['start'], axis=1)
df_bed['num_GC'] = df_bed.apply(lambda x: GC(x['chr'],x['start'],x['end']), axis=1)


df_bed['WGS_GC_bias'] = df_bed.apply(lambda x: dict_bias_wgs['smoothed_GC_bias'][(x['length'],x['num_GC'])], axis=1)
df_bed['bps_WGS_corrected'] = df_bed.apply(lambda x: x['bps']/x['WGS_GC_bias'] if x['WGS_GC_bias'] > 0.05 else 0, axis=1)

df_bed[['chr', 'start', 'end', 'TSS', 'WGS_GC_bias','bps', 'bps_WGS_corrected']].to_csv(f'{out_dir}/{sample}.bps.corrected.csv',index=False)

df_dep = df_bed.groupby('TSS').sum()[['bps', 'bps_WGS_corrected']]
df_dep = df_dep*1e6/df_dep.sum()

tpm = pd.read_table(TPM, sep = ',')
tci = df_dep.T[[]]

for tissue_, sub_dat in tpm.groupby('Specific high expression tissue'):
    tissue = tissue_.split(' -')[0]
    TSSs = list(tpm[tpm['Specific high expression tissue'] == tissue_]['TSS ID'])
    if tissue == 'Whole Blood':
        tci['Blood cell'] = 50 - (df_dep.T[list(set(TSSs).intersection(set(df_dep.index)))].sum(axis=1) / len(TSSs))
    else:
        tci[tissue] = 50 - (df_dep.T[list(set(TSSs).intersection(set(df_dep.index)))].sum(axis=1) / len(TSSs))

tci.to_csv(f'{out_dir}/{sample}.tci.csv', index_label = sample)
