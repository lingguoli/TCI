### Required Arguments ###
fastq1=$1
adapter=$2

### Optional Arguments ###
outdir=$3
if [ ! -n "$outdir" ]; then
    outdir=`pwd`;
fi
mkdir -p ${outdir};
REF=$4
if [ ! -n "$REF" ]; then
    REF="../Ref/hg38.fa"
fi

THREAD=$5
if [ ! -n "$threads" ]; then
    THREAD=8;
fi
### Preperation and Software ###
prefix=$(echo `basename $fastq1` | awk -v FS='_' '{print $1}'| awk -v FS='.' '{print $1}')

############
### Main ###
############

if [ -s "${outdir}/${prefix}.fastp.finish" ]; then
    info1="${outdir}/${prefix}.fastp.finish exists"
    echo "[Skip] fastp: ${info1}."

else
    echo "[FASTP] Start: `date`"
    time fastp \
        --in1=$fastq1 --out1=${outdir}/${prefix}.clean_R1.fq.gz \
        --adapter_sequence=${adapter} \
        --html=${outdir}/${prefix}.report.html \
        --json=${outdir}/${prefix}.report.json \
        --qualified_quality_phred=5 \
        --unqualified_percent_limit=50 \
        --n_base_limit=10 \
        --disable_trim_poly_g \
        --thread=${THREAD}
    echo "[FASTP] End: `date`"
    echo "Practice_makes_perfect" > ${outdir}/${prefix}.fastp.finish
fi

if [ -s "${outdir}/${prefix}.mkdup.bam" ] && \
   [ -s "${outdir}/${prefix}.mkdup.bam.bai" ] && \
   [ -s "${outdir}/${prefix}.align.finish" ]; then
    info1="${outdir}/${prefix}.mkdup.bam exists"
    info2="${outdir}/${prefix}.mkdup.bam.bai exists"
    info3="${outdir}/${prefix}.align.finish exists"
    echo "[Skip] align: ${info1}, ${info2} and ${info3}."

else
    echo "[ALIGN] Start: `date`"
        minimap2 -ax sr ${REF} ${outdir}/${prefix}.clean_R1.fq.gz  -t ${THREAD} |samtools sort -@ ${THREAD} -T temp_sorted_${prefix} -o ${outdir}/${prefix}.sorted.bam -
        samtools markdup ${outdir}/${prefix}.sorted.bam ${outdir}/${prefix}.mkdup.bam
        samtools index ${outdir}/${prefix}.mkdup.bam
        #rm ${outdir}/${prefix}.clean_R1.fq.gz ${outdir}/${prefix}.sorted.bam
    echo "[ALIGN] End: `date`" 
    echo "Practice_makes_perfect" > ${outdir}/${prefix}.align.finish
fi


if [ -s "${outdir}/${prefix}.bed" ] && \
   [ -s "${outdir}/${prefix}.bed.finish" ]; then
    info1="${outdir}/${prefix}.bed exists"
    info2="${outdir}/${prefix}.bed.finish exists"
    echo "[Skip] bed: ${info1} and ${info2}."

else
    echo "[BED] Start: `date`"
        python3 "$(dirname $0)/../src/bam2bed.SE.py" ${outdir}/${prefix}.mkdup.bam ${outdir}/${prefix}.bed 30 5
    echo "[BED] End: `date`"
    echo "Practice_makes_perfect" > ${outdir}/${prefix}.bed.finish
fi

if [ -s "${outdir}/${prefix}.SE.tci.csv" ] && \
   [ -s "${outdir}/${prefix}.SE.tci.csv.finish" ]; then
    info1="${outdir}/${prefix}.SE.tci.csv exists"
    info2="${outdir}/${prefix}.SE.tci.csv.finish exists"
    echo "[Skip] bed: ${info1} and ${info2}."

else
    echo "[TCI] Start: `date`"
        in_bed=${outdir}/${prefix}.bed
        SE_tci=$(dirname $0)/../src/SE_cnt_to_TCI.py
        TPM=$(dirname $0)/../Ref/TPM_atlas.csv
        TSS_region=$(dirname $0)/../Ref/ucsc.refgene.hg38.TSS.UD1210.merged.intervals.bed
        TSS=$(dirname $0)/../Ref/ucsc.refgene.hg38.TSS.UD1000.bed
    # STEP1
        if [ "${in_bed##*.}" = "gz" ] || [ "${in_bed##*.}" = "bgz" ]; then ## gzip file
          zcat ${in_bed}| bedtools intersect -a - -b ${TSS_region} -wa|bedtools intersect -b -  -a ${TSS} -wa|uniq > ${outdir}/${prefix}.SE.cnt
        else ## ungzip file
          bedtools intersect -a ${in_bed} -b ${TSS_region} -wa|bedtools intersect -b -  -a ${TSS} -wa > ${outdir}/${prefix}.SE.cnt
        fi
    # STEP2
        python3 ${SE_tci} ${TPM} ${outdir}/${prefix}.SE.cnt ${outdir}/${prefix}.SE.depth ${outdir}/${prefix}.SE.tci.csv
    #rm ${outdir}/${prefix}.SE.cnt
    echo "[TCI] End: `date`"
    echo "Practice_makes_perfect" > ${outdir}/${prefix}.SE.tci.csv.finish
fi
