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
        python3 "$(dirname $0)/../src/bam2bed.py" ${bam} ${outdir}/${prefix}.bed 30 5
    echo "[BED] End: `date`" 
    echo "Practice_makes_perfect" > ${outdir}/${prefix}.bed.finish
fi

if [ -s "${outdir}/${prefix}.tci.csv" ] && \
   [ -s "${outdir}/${prefix}.tci.csv.finish" ]; then
    info1="${outdir}/${prefix}.tci.csv exists"
    info2="${outdir}/${prefix}.tci.csv.finish exists"
    echo "[Skip] bed: ${info1} and ${info2}."

else
    echo "[TCI] Start: `date`"
        sh "$(dirname $0)/GC_bias.sh" ${outdir}/${prefix}.bed ${outdir}
    echo "[TCI] End: `date`" 
    echo "Practice_makes_perfect" > ${outdir}/${prefix}.tci.csv.finish
fi
