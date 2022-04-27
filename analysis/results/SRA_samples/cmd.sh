#download and analyze PE and SE SRA samples
grep _1.fastq.gz /srv/scratch/oateslab/rawData/2021/SRA_RNA_seq/*/meta.txt|cut -f12|perl -pe 's/_1.fastq.gz//'|l >pe_samples.txt
grep -vP "_1.fastq.gz|_2.fastq.gz" /srv/scratch/oateslab/rawData/2021/SRA_RNA_seq/*/meta.txt|grep fastq.gz|cut -f12|perl -pe 's/.fastq.gz//'|l >se_samples.txt

cat pe_samples.txt|awk '{print "sh analyze_one_PE_sample.sh "$1}'|sh
cat se_samples.txt|awk '{print "sh analyze_one_SE_sample.sh "$1}'|sh

#count reads at either end of the genes
for i in analysis/results/SRA_samples/star/*.Aligned.sortedByCoord.out.md.bam; do ba=`basename $i`; echo "sh analysis/src/stat_with_end_exons_reads//CountTwoEndsReads.sh $i analysis/results/SRA_samples/star/";done|split -l 10 - tmp/stat_ee.
for i in tmp/stat_ee.*;do ba=`basename $i`; cat analysis/src/stat_with_end_exons_reads/template.pbs $i >pbsee/$ba.pbs;qsub pbsee/$ba.pbs;done
mv star/*3prime* prime_star/&
ls prime_star/*intron|cut -d "/" -f2|cut -d "." -f1|sort|uniq -d|l >quanti_done_prime_end
for i in `ls analysis/results/SRA_samples/star/*.Aligned.sortedByCoord.out.md.bam|fgrep -v -f quanti_done_prime_end -`; do ba=`basename $i`; echo "sh analysis/src/stat_with_end_exons_reads//CountTwoEndsReads.sh $i analysis/results/SRA_samples/star/"|cat analysis/src/stat_with_end_exons_reads/template.pbs - >pbsee/$ba.pbs;qsub pbsee/$ba.pbs;done
mv prime_star count_reads_5p
mv count_reads_5p/*3prime* count_reads_3p/
