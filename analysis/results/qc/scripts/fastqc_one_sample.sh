#sample
s=$1
fq=$2
outdir="analysis/results/qc/fastqc_out/$s";
mkdir -p $outdir;
mkdir $TMPDIR/fastqc;
module add java fastqc
fastqc -o $outdir --extract -f fastq -d $TMPDIR/fastqc $fq;
