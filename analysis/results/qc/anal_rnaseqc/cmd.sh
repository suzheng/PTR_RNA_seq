#get the QC result header
l ../rnaseqc_out/E-MTAB-7247.human_brain_rna_1_ERR2812349/human_brain_rna_1_ERR2812349.Aligned.sortedByCoord.out.md.bam.metrics.tsv|cut -f1|perl -ne 's/\s$/\t/;print'|perl -pe 's/\t$/\n/'|l >header
#collect the QC result for all samples
for i in ../rnaseqc_out/*/*metrics.tsv;do awk -v i=$i 'BEGIN{print "project\t"i}NR>0' $i|cut -f2|paste -s -;done|cut -d "/" -f3-|perl -pe 's/\//\t/g'|cut -f1,4-|cat header - >metrics.tsv.pasted 
#clean up the QC result
l metrics.tsv.pasted|grep -v Pan_troglodytes|grep -v Macaca_mulatta|grep -v Mus_musculus|grep -v E-MTAB-7247|grep -v SRP127360|l >metrics.tsv.pasted.clean
#extract the QC-failed SRA human tissue samples
for i in ../rnaseqc_out/*/*metrics.tsv;do sh filter.sh $i;done >filter.out1
l filter.out1 |grep -v Pan_troglodytes|grep -v Macaca_mulatta|grep -v Mus_musculus|grep -v E-MTAB-7247|grep -v SRP127360|l >filter.out1.clean
l filter.out1.clean|grep FAILED|cut -d "/" -f3|cut -d "." -f2|l >failed_samples
#extract the QC-failed non-human species samples and GTEx samples
for i in analysis/results/qc/rnaseqc_out/nhs*/*metr*;do sh filter.sh $i;done >filter.out.nhs
for i in analysis/results/qc/rnaseqc_out/GTE*/*metr*;do sh filter.sh $i;done >filter.out.GTEx
l filter.out.GTEx |grep FAIL|cut -d "/" -f13|perl -pe 's/GTEx.//'|sort|uniq|l >filter.out.GTEx.failed
l filter.out.nhs|grep FAILED|cut -d "/" -f13|perl -pe 's/nhs.//'|sort|uniq|sort|uniq|l >filter.out.nhs.failed
