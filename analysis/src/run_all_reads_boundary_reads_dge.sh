cd analysis/results/100_samples_each_tissue/DGE_analysis

n=`cat tissue_list|wc -l`;for i in $(seq 1 $n);do i2=`expr $i + 1`; for j in $(seq $i2 $n );do g1=`head -$i tissue_list|tail -n1`;g2=`head -$j tissue_list|tail -n1`;echo "snakemake --cores all --use-singularity --singularity-args \"-B analysis/src/dge_functions_data:analysis/src/dge_functions_data\" -C out_dir=$PWD -s analysis/src/dge_functions_data/Snakefile_dge.GTEx.py DGE_out/figures_tables.${g1}.${g2}.out/DESeq2.exon_counts_subsampled.DESeq2.${g1}.vs.${g2}.txt"|cat template.pbs - >pbs/dge.$g1.$g2.pbs;qsub pbs/dge.$g1.$g2.pbs;done;done

cd analysis/results/SRA_samples/DGE_analysis

for i in meta_ERP003613_SRP028336/*tsv;do p=`echo $i|cut -d "/" -f 2|perl -pe 's/.samples.tsv//'`;g1=`cat $i|cut -f2|grep -v condition|sort|uniq|head -1`;g2=`cat $i|cut -f2|grep -v condition|sort|uniq|tail -n1`;echo "snakemake --cores all -F --use-singularity --singularity-args \"-B analysis/src/dge_functions_data:analysis/src/dge_functions_data\" -C out_dir=$PWD meta_dir=meta_ERP003613_SRP028336 -s analysis/src/dge_functions_data/Snakefile_dge.SRA.py  $p/figures_tables.${g1}.${g2}.out/DESeq2.exon_counts_subsampled.DESeq2.${g1}.vs.${g2}.txt"|cat template.pbs - >pbs/$p.pbs;qsub pbs/$p.pbs;done
for i in meta/*tsv;do p=`echo $i|cut -d "/" -f 2|perl -pe 's/.samples.tsv//'`;g1=`cat $i|cut -f2|grep -v condition|sort|uniq|head -1`;g2=`cat $i|cut -f2|grep -v condition|sort|uniq|tail -n1`;echo "snakemake --cores all -F --use-singularity --singularity-args \"-B analysis/src/dge_functions_data:analysis/src/dge_functions_data\" -C out_dir=$PWD meta_dir=meta -s analysis/src/dge_functions_data/Snakefile_dge.SRA.py  $p/figures_tables.${g1}.${g2}.out/DESeq2.exon_counts_subsampled.DESeq2.${g1}.vs.${g2}.txt"|cat template.pbs - >pbs/$p.pbs;qsub pbs/$p.pbs;done

cd analysis/results/non_human_species/DGE_analysis
for i in meta_SRP028336/*tsv;do p=`echo $i|cut -d "/" -f 2|perl -pe 's/.samples.tsv//'`;g1=`cat $i|cut -f2|grep -v condition|sort|uniq|head -1`;g2=`cat $i|cut -f2|grep -v condition|sort|uniq|tail -n1`;echo "snakemake --cores all -F --use-singularity --singularity-args \"-B analysis/src/dge_functions_data:analysis/src/dge_functions_data\" -s Snakefile_dge.py $p/figures_tables.${g1}.${g2}.out/DESeq2.exon_counts_subsampled.DESeq2.${g1}.vs.${g2}.txt"|cat template.pbs - >pbs/$p.pbs;qsub pbs/$p.pbs;done

cd analysis/results/stat_with_boundary_reads/DGE_analysis
n=`cat tissue_list|wc -l`;for i in $(seq 1 $n);do i2=`expr $i + 1`; for j in $(seq $i2 $n );do g1=`head -$i tissue_list|tail -n1`;g2=`head -$j tissue_list|tail -n1`;echo "snakemake --cores all --use-singularity --singularity-args \"-B analysis/src/dge_functions_data:analysis/src/dge_functions_data\" -C out_dir=$PWD -s analysis/src/dge_functions_data/Snakefile_dge.GTEx.py DGE_out/figures_tables.${g1}.${g2}.out/DESeq2.exon_counts_subsampled.DESeq2.${g1}.vs.${g2}.txt"|cat template.pbs - >pbs/dge.$g1.$g2.pbs;qsub pbs/dge.$g1.$g2.pbs;done;done

cd analysis/results/boundary_reads_only_SRA/DGE_analysis
for i in meta/*tsv;do p=`echo $i|cut -d "/" -f 2|perl -pe 's/.samples.tsv//'`;g1=`cat $i|cut -f2|grep -v condition|sort|uniq|head -1`;g2=`cat $i|cut -f2|grep -v condition|sort|uniq|tail -n1`;echo "snakemake --cores all -F --use-singularity --singularity-args \"-B analysis/src/dge_functions_data:analysis/src/dge_functions_data\" -C out_dir=$PWD meta_dir=meta -s analysis/src/dge_functions_data/Snakefile_dge.SRA.py  $p/figures_tables.${g1}.${g2}.out/DESeq2.exon_counts_subsampled.DESeq2.${g1}.vs.${g2}.txt"|cat template.pbs - >pbs/$p.pbs;qsub pbs/$p.pbs;done
for i in meta_ERP003613_SRP028336/*tsv;do p=`echo $i|cut -d "/" -f 2|perl -pe 's/.samples.tsv//'`;g1=`cat $i|cut -f2|grep -v condition|sort|uniq|head -1`;g2=`cat $i|cut -f2|grep -v condition|sort|uniq|tail -n1`;echo "snakemake --cores all -F --use-singularity --singularity-args \"-B analysis/src/dge_functions_data:analysis/src/dge_functions_data\" -C out_dir=$PWD meta_dir=meta_ERP003613_SRP028336 -s analysis/src/dge_functions_data/Snakefile_dge.SRA.py  $p/figures_tables.${g1}.${g2}.out/DESeq2.exon_counts_subsampled.DESeq2.${g1}.vs.${g2}.txt"|cat template.pbs - >pbs/$p.pbs;qsub pbs/$p.pbs;done
