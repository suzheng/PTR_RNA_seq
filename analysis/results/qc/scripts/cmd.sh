#get the SRA sample list
ls analysis/results/SRA_samples/star/*bam|perl -pe 's/.Aligned.sortedByCoord.out.md.bam//'|cut -d "/" -f13|l >sample_list
for i in `cat sample_list`;do echo "sh qc_one_sample.sh $i";done|split -l 10 --additional-suffix=.txt - tmp/split
#Perform QC on the SRA samples
for i in tmp/*txt;do ba=`basename $i`;cat template.pbs $i >pbs/$ba.pbs;qsub pbs/$ba.pbs;done
#get the GTEx bam file list
ls analysis/results/stat_with_boundary_reads/out/*/*bam|grep -v Boundar >oateslab_file_list.GTEx
#get the GTEx sample list
cat oateslab_file_list.GTEx|awk -F "/" '{print $(NF-1)}'|sort|uniq|l >GTEx.samples
#Perform QC on the GTEx samples
for i in `cat GTEx.samples`;do bam=`grep $i oateslab_file_list.GTEx|cut -d ":" -f2-|head -1`;echo "sh qc_one_gtex_sample.sh $i $bam";done|split --additional-suffix=.cmd -l 10 - tmp/split.
for i in tmp/split.*cmd;do ba=`basename $i`;cat template.pbs $i >pbs/$ba.pbs;qsub pbs/$ba.pbs;done

#get the non-human species sample list and gtf file list
ls analysis/results/non_human_species/star/*bam|grep SRR|l >non_human_species1.bams
ls /srv/scratch/oateslab/share/data/non_human_species/collapsed_genes_gtf/*gtf >non_human_species_gtf_list
#perform QC on the non-human tissue samples
for i in `cat non_human_species1.bams`;do s=`echo $i|awk -F "/" '{print $NF}'|cut -d "." -f1`;gtf=`echo $i|cut -d "_" -f6-7|fgrep -i -f - non_human_species_gtf_list`;echo "sh qc_one_nhs_sample.sh $s $i $gtf"|cat template.pbs - >pbs/run.$s.pbs;qsub pbs/run.$s.pbs;done
