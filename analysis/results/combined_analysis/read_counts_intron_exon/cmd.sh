#collect mean intronic and exonic read count by gene
for i in analysis/results/100_samples_each_tissue/DGE_analysis/merged_counts/*intron.merged.txt;do cat $i|cut -f7-|awk 'NR>1'|awk '{sum=0; for (i=1; i<=NF; i++) {sum=sum+$i;} m=sum/NF; print m; }';done >gene_intron_counts
for i in analysis/results/100_samples_each_tissue/DGE_analysis/merged_counts/*exon.merged.txt;do cat $i|cut -f7-|awk 'NR>1'|awk '{sum=0; for (i=1; i<=NF; i++) {sum=sum+$i;} m=sum/NF; print m; }';done >gene_exon_counts
#collect mean intronic and exonic read count by sample
for i in analysis/results/100_samples_each_tissue/out/*/*intron;do cat $i|awk 'NR>2{a+=$NF}END{print a/(NR-2)}';done >mean_intronic_count_by_sample
for i in analysis/results/100_samples_each_tissue/out/*/*exon;do cat $i|awk 'NR>2{a+=$NF}END{print a/(NR-2)}';done >mean_exonic_count_by_sample&
