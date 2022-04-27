#sample
s=$1
indir=/srv/scratch/oateslab/rawData/2021/SRA_RNA_seq;
gtf=/srv/scratch/oateslab/share/data/hg38/collapsed_genes_gtf/gencode.v39.GRCh38.genes.gtf
bam_dir=analysis/results/SRA_samples/star/
p=`ls $indir/*/$s*.gz|head -1|cut -d "/" -f8`;
outdir="analysis/results/qc/fastqc_out/$p.$s";
rnaseqc_out="analysis/results/qc/rnaseqc_out/$p.$s";
mkdir -p $outdir $rnaseqc_out;
mkdir $TMPDIR/fastqc;
module add java fastqc
fastqc -o $outdir --extract -f fastq -d $TMPDIR/fastqc $indir/*/$s*.gz;
rnaseqc $gtf $bam_dir/$s.Aligned.sortedByCoord.out.md.bam $rnaseqc_out -u
