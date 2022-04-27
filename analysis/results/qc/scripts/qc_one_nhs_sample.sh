#sample
s=$1
bam=$2
gtf=$3
rnaseqc_out="analysis/results/qc/rnaseqc_out/nhs.$s";
mkdir -p $rnaseqc_out;
module add java fastqc
rnaseqc $gtf $bam $rnaseqc_out -u
