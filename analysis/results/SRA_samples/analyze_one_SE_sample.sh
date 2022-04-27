sample=$1

data_dir=/srv/scratch/oateslab/rawData/2021/SRA_RNA_seq
anal_dir=analysis/results/SRA_samples

#do the download, file structure has to be SRP_folder/meta.txt
cd $data_dir
SRP_dir=`grep $sample */meta.txt|cut -d "/" -f1`
cd $SRP_dir
grep $sample meta.txt|awk -F "\t" '{print "wget -c "$9" -O "$12}'|sh

#do the analysis
cd $anal_dir
mkdir pbs
echo "snakemake --cores all -s Snakefile.SE count_reads/$sample.Aligned.sortedByCoord.out.md.bam.intron" |cat template.pbs - >pbs/$sample.pbs
qsub pbs/$sample.pbs
