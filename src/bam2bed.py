import os,sys,pysam

bam, bed, min_map_qual, max_num_of_mismatch = sys.argv[1:]
min_map_qual = int(min_map_qual)
max_num_of_mismatch = int(max_num_of_mismatch)

reference_list=[]
for i in list(range(1,23))+['X','Y','M']:
    reference_list.append(f'chr{i}') #  chr1-22 X Y M only


sam=pysam.AlignmentFile(bam, 'rb')
with open(bed,'w') as bed_file:
    for chrx in reference_list:
        filter_ID = set()
        for r in sam.fetch(chrx):
            if r.query_name in filter_ID : continue   #
            elif r.mapping_quality < min_map_qual:
                filter_ID.add(r.query_name); continue #  mapping quality cutoff
            elif r.reference_id != r.next_reference_id  :
                filter_ID.add(r.query_name); continue #  the two reads should map to the sample chr
            elif int(r.is_reverse)+int(r.mate_is_reverse) != 1  :
                filter_ID.add(r.query_name); continue #  keep "only one read is reverse complemented"
            elif r.mate_is_unmapped  :
                filter_ID.add(r.query_name); continue #  filter "next segment in the template unmapped"
            elif abs(r.template_length) > 600  :
                filter_ID.add(r.query_name); continue #  filter fragments with length more than 600bp
            elif len(r.cigar) != 1 :
                filter_ID.add(r.query_name); continue #  kepp fragments with 100% matche
            elif r.is_duplicate :
                filter_ID.add(r.query_name); continue #  the read is either a PCR duplicate or an optical duplicate
            elif r.is_secondary  :
                filter_ID.add(r.query_name); continue #  filter secondary alignment
            elif r.get_tag('NM') > max_num_of_mismatch:
                filter_ID.add(r.query_name); continue #  filter DNA fragment with many mismatches

        for r in sam.fetch(chrx):
            if r.is_read2: continue
            elif r.query_name not in filter_ID:
                r_len = abs(r.template_length)
                if r.is_reverse :
                    print(r.reference_name, r.reference_end - r_len, r.reference_end , r.mapping_quality, '-', r_len, sep = '\t', file = bed_file)
                else :
                    print(r.reference_name, r.reference_start, r.reference_start + r_len, r.mapping_quality, '+', r_len, sep = '\t', file = bed_file)

