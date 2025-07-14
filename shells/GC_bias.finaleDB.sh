bed=$1
dir=$2
sample=$(basename ${bed}|awk -v FS='.' '{print $1}')
scr_dir=$(dirname $0)
mkdir -p ${dir}

hg38_intervals=${scr_dir}/../Ref/k100_minus_exclusion_lists.mappable_regions.hg38.bed
tss_intervals=${scr_dir}/../Ref/ucsc.refgene.hg38.TSS.UD1210.merged.intervals.bed

if [ ! -e "${dir}/wgs/GC_bias/${sample}.GC_bias.txt" ]; then
    bedtools intersect -a ${bed} -b ${hg38_intervals} -wa |awk '$4>=30' |uniq > ${dir}/${sample}.mappable.hg38.bed
    bedtools intersect -a ${bed} -b ${tss_intervals}  -wa |awk '$4>=30' |uniq > ${dir}/${sample}.tss.hg38.bed
    
    python3 ${scr_dir}/../src/TCI.GC_count.wgs.py ${dir}/${sample}.mappable.hg38.bed ${dir}/wgs
    python3 ${scr_dir}/../src/TCI.GC_bias.wgs.py  ${dir}/${sample}.mappable.hg38.bed ${dir}/wgs
   
    rm ${dir}/${sample}.mappable.hg38.bed
    rm -rf ${dir}/tmp_${sample}  ${dir}/wgs/tmp_${sample}
else
    echo "${sample} GC bias finished."
fi

if [ ! -e "${dir}/${sample}.tci.txt" ]; then
    awk -v min=150 -v max=210 -v OFS='\t' '($3-$2)>=min && ($3-$2)<max{print $1,$2,$3}' ${dir}/${sample}.tss.hg38.bed | bedtools intersect -a - -b ${scr_dir}/../Ref/ucsc.refgene.hg38.TSS.UD1000.bed -wao | awk -v OFS='\t' '$8>0{print $1,$2,$3,$7,$8}' > ${dir}/${sample}.tss.hg38.bps.bed
    python3 ${scr_dir}/../src/TCI.GC_corrected.wgs.py ${dir}/${sample}.tss.hg38.bps.bed ${dir}/wgs/GC_bias/${sample}.GC_bias.txt ${dir}
else
    echo "${sample} TCI finished."
fi
