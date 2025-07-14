import os,sys,pysam

bam, bed, min_map_qual, max_num_of_mismatch = sys.argv[1:]
min_map_qual = int(min_map_qual)
max_num_of_mismatch = int(max_num_of_mismatch)

reference_list=[]
for i in list(range(1,23))+['X','Y','M']:
    reference_list.append(f'chr{i}')
reference_list=set(reference_list)
sam=pysam.AlignmentFile(bam, 'rb')
with open(bed,'w') as bed_file:
    for r in sam:
        if r.reference_name not in reference_list:
            continue #  chr1-22 X Y M only
        elif r.cigarstring!=f'{r.rlen}M':
            continue # only 100% map is allowed
        elif r.mapping_quality < min_map_qual:
            continue #  mapping quality cutoff
        elif r.is_duplicate : 
            continue #  the read is either a PCR duplicate or an optical duplicate
        elif r.is_unmapped  : 
            continue #  new: filter unmapped seqment
        elif r.is_secondary : 
            continue #  filter secondary alignment
        elif r.get_tag('NM') > max_num_of_mismatch: 
            continue #  filter DNA fragment with many mismatches
        print(r.reference_name, r.reference_start, r.reference_end, r.mapping_quality, '-' if r.is_reverse else '+', sep = '\t', file = bed_file)

