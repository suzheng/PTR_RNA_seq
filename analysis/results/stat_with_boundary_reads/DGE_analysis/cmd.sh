#get the sucessfully quantified files
ls analysis/results/stat_with_boundary_reads/out/*/*intron|grep -v prime|fgrep -f analysis/results/100_samples_each_tissue/DGE_analysis/all_Reads_QC_passed_samples - >all.intron_files
ln ../../100_samples_each_tissue/DGE_analysis/grouping_info ./
#create meta data file for each tissue
for i in `cut -f2 grouping_info|sort|uniq`;do grep -w $i grouping_info|cut -f1|fgrep -w -f - all.intron_files|perl -pe 's/.intron//'|l >meta/$i.prefix_list;done
wc -l meta/*|grep -v total|awk '$1>5{print $2}'|cut -d "/" -f2|cut -d "." -f1 >tissue_list
mkdir merged_counts/
#collect read count data for each tissue
for t in `cat tissue_list`;do echo "sh scripts/collectCounts.sh meta/$t.prefix_list merged_counts/$t&";done|l
#perform differential expression analysis
n=`cat tissue_list|wc -l`;for i in $(seq 1 $n);do i2=`expr $i + 1`; for j in $(seq $i2 $n );do g1=`head -$i tissue_list|tail -n1`;g2=`head -$j tissue_list|tail -n1`;echo "snakemake --cores all --use-singularity --singularity-args \"-B analysis/src/dge_functions_data:analysis/src/dge_functions_data\" -C out_dir=$PWD -s analysis/src/dge_functions_data/Snakefile_dge.GTEx.py DGE_out/figures_tables.${g1}.${g2}.out/DESeq2.exon_counts_subsampled.DESeq2.${g1}.vs.${g2}.txt"|cat template.pbs - >pbs/dge.$g1.$g2.pbs;qsub pbs/dge.$g1.$g2.pbs;done;done
