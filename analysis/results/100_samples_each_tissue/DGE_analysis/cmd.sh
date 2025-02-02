#get the list of all quantification files for intron
ls analysis/results/100_samples_each_tissue/out/*/*intron >all.intron_files
#extract the sample grouping info
cat ../info_for_samples_to_donwload.txt|cut -f2,15|perl -F"\t" -nale '$F[1]=~s/\s//g;$F[1]=~s/\W/_/g;$F[1]=~s/_$//;print "$F[0]\t$F[1]"'|sort|uniq >grouping_info
#generate the meta data file for sucessfully quantified samples
for i in `cut -f2 grouping_info|sort|uniq`;do grep -w $i grouping_info|cut -f1|fgrep -w -f - all.intron_files|perl -pe 's/.intron//'|l >meta/$i.prefix_list;done
l grouping_info |cut -f2|sort|uniq|l >tissue_list

#collect read count data for every tissue
for t in `cat tissue_list`;do echo "sh scripts/collectCounts.sh meta/$t.prefix_list merged_counts/$t";done|l
l all.intron_files|awk -F "/" '{print $NF}'|cut -d "." -f1|fgrep -v -f analysis/results/qc/anal_rnaseqc/filter.out.GTEx.failed - >all_Reads_QC_passed_samples

#perform differential expression analysis
n=`cat tissue_list|wc -l`;for i in $(seq 1 $n);do i2=`expr $i + 1`; for j in $(seq $i2 $n );do g1=`head -$i tissue_list|tail -n1`;g2=`head -$j tissue_list|tail -n1`;echo "snakemake --cores all --use-singularity --singularity-args \"-B analysis/src/dge_functions_data:analysis/src/dge_functions_data\" -C out_dir=$PWD -s analysis/src/dge_functions_data/Snakefile_dge.GTEx.py DGE_out/figures_tables.${g1}.${g2}.out/DESeq2.exon_counts_subsampled.DESeq2.${g1}.vs.${g2}.txt"|cat template.pbs - >pbs/dge.$g1.$g2.pbs;qsub pbs/dge.$g1.$g2.pbs;done;done

#select only 5 samples in each tissue to perform differential expression analysis
mv DGE_out DGE_out_all_samples
then edit the read_counts_file function to only read the first N samples analysis/src/dge_functions_data/shared_functions.R
mv DGE_out DGE_out_5samples
mv DGE_out_all_samples DGE_out
for i in pbs/*pbs;do ba=`basename $i`;echo "sh $i 1>logs/$ba.log 2>logs/$ba.err";done|parallel -j 12
ls pbs/*pbs|awk '{print "qsub "$1}'|l
