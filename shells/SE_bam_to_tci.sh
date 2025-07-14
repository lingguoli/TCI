### Required Arguments ###
bam=$1

### Optional Arguments ###
outdir=$2
if [ ! -n "$outdir" ]; then
    outdir=`pwd`;
fi
mkdir -p ${outdir};
REF=$3
if [ ! -n "$REF" ]; then
    REF="../Ref/hg38.fa"
fi

THREAD=$4
if [ ! -n "$threads" ]; then
    THREAD=8;
fi
### Preperation and Software ###
prefix=$(echo `basename $bam` | awk -v FS='_' '{print $1}'| awk -v FS='.' '{print $1}')

############
### Main ###
############

if [ -s "${outdir}/${prefix}.bed" ] && \
   [ -s "${outdir}/${prefix}.bed.finish" ]; then
    info1="${outdir}/${prefix}.bed exists"
    info2="${outdir}/${prefix}.bed.finish exists"
    echo "[Skip] bed: ${info1} and ${info2}."

else
    echo "[BED] Start: `date`"
        python3 "$(dirname $0)/../src/bam2bed.SE.py" ${bam} ${outdir}/${prefix}.bed 30 5
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
