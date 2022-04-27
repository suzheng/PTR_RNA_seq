#collect the sucessfully quantified intron files
ls analysis/results/stat_with_boundary_reads/out/*/*5prime.intron >all.5pintron_files
ls analysis/results/stat_with_boundary_reads/out/*/*3prime.intron >all.3pintron_files
#create meta data file for each tissue
for p in 5p 3p;do for i in `cut -f2 grouping_info|sort|uniq`;do grep -w $i grouping_info|cut -f1|fgrep -w -f - all.${p}intron_files|perl -pe 's/.intron//'|l >${p}_meta/$i.prefix_list;done;done
l ../../100_samples_each_tissue/DGE_analysis/grouping_info|cut -f2|sort|uniq|l >tissue_list
#collect read count for each tissue
p=5p;for t in `cat tissue_list`;do echo "sh scripts/collectCounts.sh ${p}_meta/$t.prefix_list ${p}_merged_counts/$t";done|parallel -j 10
p=3p;for t in `cat tissue_list`;do echo "sh scripts/collectCounts.sh ${p}_meta/$t.prefix_list ${p}_merged_counts/$t";done|parallel -j 10
#perform differential expression analysis for 5'end reads
mv 5p_merged_counts merged_counts
n=`cat tissue_list|wc -l`;for i in $(seq 1 $n);do i2=`expr $i + 1`; for j in $(seq $i2 $n );do g1=`head -$i tissue_list|tail -n1`;g2=`head -$j tissue_list|tail -n1`;echo "snakemake --cores all --use-singularity --singularity-args \"-B analysis/src/dge_functions_data:analysis/src/dge_functions_data\" -C out_dir=$PWD -s analysis/src/dge_functions_data/Snakefile_dge.GTEx.py DGE_out/figures_tables.${g1}.${g2}.out/DESeq2.exon_counts_subsampled.DESeq2.${g1}.vs.${g2}.txt"|cat template.pbs - >pbs/dge.$g1.$g2.pbs;qsub pbs/dge.$g1.$g2.pbs;done;done
mv DGE_out 5p_DGE_out; 
mv merged_counts 5p_merged_counts;
#perform differential expression analysis for 3'end reads
mv 3p_merged_counts merged_counts
n=`cat tissue_list|wc -l`;for i in $(seq 1 $n);do i2=`expr $i + 1`; for j in $(seq $i2 $n );do g1=`head -$i tissue_list|tail -n1`;g2=`head -$j tissue_list|tail -n1`;echo "snakemake --cores all --use-singularity --singularity-args \"-B analysis/src/dge_functions_data:analysis/src/dge_functions_data\" -C out_dir=$PWD -s analysis/src/dge_functions_data/Snakefile_dge.GTEx.py DGE_out/figures_tables.${g1}.${g2}.out/DESeq2.exon_counts_subsampled.DESeq2.${g1}.vs.${g2}.txt"|cat template.pbs - >pbs/dge.$g1.$g2.pbs;qsub pbs/dge.$g1.$g2.pbs;done;done
mv DGE_out 3p_DGE_out; 
mv merged_counts 3p_merged_counts;
