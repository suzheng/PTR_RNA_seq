in=$1
ba=${in%.bam}
# sort reads by identifier-name (-n)
module add bedtools
#samtools sort -n $in -o $ba.sortn.bam
bedtools bamtofastq -i $ba.sortn.bam -fq ${ba}_1.fastq -fq2 ${ba}_2.fastq && gzip -f ${ba}_1.fastq && gzip -f ${ba}_2.fastq
#bedtools bamtofastq -i $ba.sortn.bam -fq ${ba}_1.fastq -fq2 ${ba}_2.fastq && gzip ${ba}_1.fastq && gzip ${ba}_2.fastq && rm $ba.sortn.bam
