#get the directory list for all_reads (not only boundary reads or 5', 3' end reads) for SRA tissue and diseases (abnormal_condition).
dt=all_reads;g=sra_tissues;d=analysis/results/SRA_samples/DGE_analysis/;l analysis/results/SRA_samples/DGE_analysis/unfiltered_meta_ERP003613_SRP028336/good_comparison_list|perl -pe 's/.samples./\t/'|awk -v d=$d '{print "ls -d "d"/"$1"/*.out"}'|sh|l >$g.$dt
dt=all_reads;g=abnormal_conditions;d=analysis/results/SRA_samples/DGE_analysis/;l analysis/results/SRA_samples/DGE_analysis/unfiltered_meta/good_comparison_list|perl -pe 's/.samples./\t/'|awk -v d=$d '{print "ls -d "d"/"$1"/*.out"}'|sh|l >$g.$dt
#get the directory list for boundary reads for SRA tissue and diseases (abnormal_condition) datasets
dt=boundary_reads;g=sra_tissues;d=analysis/results/boundary_reads_only_SRA/DGE_analysis/;l analysis/results/SRA_samples/DGE_analysis/unfiltered_meta_ERP003613_SRP028336/good_comparison_list|perl -pe 's/.samples./\t/'|awk -v d=$d '{print "ls -d "d"/"$1"/*.out"}'|sh|l >$g.$dt
dt=boundary_reads;g=abnormal_conditions;d=analysis/results/boundary_reads_only_SRA/DGE_analysis/;l analysis/results/SRA_samples/DGE_analysis/unfiltered_meta/good_comparison_list|perl -pe 's/.samples./\t/'|awk -v d=$d '{print "ls -d "d"/"$1"/*.out"}'|sh|l >$g.$dt
#get the directory list for 5' and 3' end reads for SRA tissue and diseases (abnormal_condition) datasets
for prime in 5 3;do dt=${prime}prime_reads;g=sra_tissues;d=analysis/results/SRA_samples/${prime}eeDGE_analysis/;l analysis/results/SRA_samples/DGE_analysis/unfiltered_meta_ERP003613_SRP028336/good_comparison_list|perl -pe 's/.samples./\t/'|awk -v d=$d '{print "ls -d "d"/"$1"/*.out"}'|sh|l >$g.$dt;done
for prime in 5 3;do dt=${prime}prime_reads;g=abnormal_conditions;d=analysis/results/SRA_samples/${prime}eeDGE_analysis/;l analysis/results/SRA_samples/DGE_analysis/unfiltered_meta/good_comparison_list|perl -pe 's/.samples./\t/'|awk -v d=$d '{print "ls -d "d"/"$1"/*.out"}'|sh|l >$g.$dt;done
#get the directory list for non-human tissues
ls -d analysis/results/non_human_species/DGE_analysis/*/*.out|grep -v Homo_sapiens|l >non-human_species.all_reads
#get the directory list for GTEx tissues
ls -d analysis/results/100_samples_each_tissue/DGE_analysis/DGE_out/*.out >gtex_tissue.all_reads
ls -d analysis/results/stat_with_boundary_reads/DGE_analysis/DGE_out/*.out >gtex_tissue.boundary_reads
ls -d analysis/results/stat_with_boundary_reads/eeDGE_analysis/5p_DGE_out/*.out >gtex_tissue.5prime_reads
ls -d analysis/results/stat_with_boundary_reads/eeDGE_analysis/3p_DGE_out/*.out >gtex_tissue.3prime_reads
ls -d analysis/results/100_samples_each_tissue/DGE_analysis/DGE_out_5samples/*.out >gtex_tissue.5samples_reads
#get the directory list for 5 sample GTEx analysis
l gtex_tissue.all_reads|perl -pe 's/DGE_out/DGE_out_5samples/'|l >gtex_tissues_5_samples.all_reads
#collect the fold change stat data 
for f in *reads;do for t in AllReads.DESeq2CPM5 AllReads.edgeRCPM5 DESeq2.logCPM0 DESeq2.logCPM5 edgeR.logCPM0 edgeR.logCPM5;do for i in `cat $f`;do awk '{print FILENAME"\t"$0}' $i/$t.fc_stat.stat.txt|perl -pe 's/DGE_analysis\//\t/;'|cut -f2-|grep -v V1|perl -pe 's/\s/\t/g;s/\s+$/\n/;'|cut -f1,3-;done >fc_stat_stat/$f.$t.stat;done;done
#collect the DESeq2 differential expression analysis data for exonic reads sub-sampled analysis
for f in `ls *all_reads|grep -v 5_samples`;do for i in `cat $f`;do awk '{print FILENAME"\t"$0}' $i/DESeq2.exon_counts_subsampled.DESeq2*.txt|perl -pe 's/DGE_analysis[\/]+/\t/;'|cut -f2-|perl -pe 's/\/deseq2/\t/i'|cut -f1,3-|grep -v log2FoldChange;done >stat_collected/DEG_tables.DESeq2.exon_counts_subsampled.$f.table;done
#collect the DESeq2 differential expression analysis data for exonic reads NOT sub-sampled analysis
for f in `ls *all_reads|grep -v 5_samples`;do for i in `cat $f`;do awk '{print FILENAME"\t"$0}' $i/DESeq2.norm_intron_counts.DESeq2*.txt|perl -pe 's/DGE_analysis[\/]+/\t/;'|cut -f2-|perl -pe 's/\/deseq2/\t/i'|cut -f1,3-|grep -v log2FoldChange;done >stat_collected/DEG_tables.DESeq2.norm_intron_counts.$f.table;done
#perform filtering on the collected DESeq2 differential expression analysis data
for i in stat_collected/DEG_tables.DESeq2.*.table;do awk '$(NF-1)<0.05 && ($4>0.5 || $4<-0.5){print $1"\t"$2"\t"$4"\t"$3"\t"$6"\t"$7}' $i >$i.sig;done
#collect the exon intron cbind (column combined) data
for t in edgeR DESeq2;do for f in `ls *.all_reads|fgrep -v 5_sam`;do for i in `cat $f`;do awk '{print FILENAME"\t"$0}' $i/$t.logCPM5.exon_intron.cbind_tables.txt|perl -pe 's/DGE_analysis[\/]+/\t/;'|cut -f2-|perl -pe 's/\.out\//\t/'|cut -f1,3-|awk '$3!=""';done|awk 'NR==1 || $0!~"logFC"' >tables_cbind/$t.logCPM5.exon_intron.$f.cbind_tables;done;done
#collect the fc_stat data for all_reads analysis in 4 datasets
for f in abnormal_conditions.all_reads  gtex_tissue.all_reads non-human_species.all_reads  sra_tissues.all_reads;do for t in DESeq2.logCPM5;do for i in `cat $f`;do awk '{print FILENAME"\t"$0}' $i/$t.fc_stat.data.txt|perl -pe 's/DGE_analysis\//\t/;'|cut -f2-|grep -v e_fc|perl -pe 's/\s/\t/g;s/\s+$/\n/;s/\/figures_tables./-/;s/^\///;s/.fc_stat.data.txt//'|cut -f1-;done >stat_collected/$f.$t.fc_data;done;done
