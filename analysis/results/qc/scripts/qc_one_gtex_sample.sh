#sample
s=$1
bam=$2
gtf=/srv/scratch/oateslab/share/data/hg38/collapsed_genes_gtf/gencode.v39.GRCh38.genes.gtf
rnaseqc_out="analysis/results/qc/rnaseqc_out/GTEx.$s";
mkdir -p $rnaseqc_out;
module add java fastqc
rnaseqc $gtf $bam $rnaseqc_out -u
