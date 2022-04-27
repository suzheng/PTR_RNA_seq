#extract the
#1. differential expression magnitude amplified, up-regulated genes
#2. differential expression magnitude amplified, down-regulated genes
#3. differential expression magnitude reduced, up-regulated genes
#4. differential expression magnitude reduced, down-regulated genes
for i in DESeq2.logCPM5.exon_intron.*.all_reads.cbind_tables;do cat $i|perl -nale 'if((abs($F[3])>1 && $F[6]<0.05) || (abs($F[9])>1 && $F[12]<0.05)){if($F[3]*$F[9]>0 && abs($F[3]) > abs($F[9]) && $F[3] >0){print}}' |cut -f2|sort|uniq -c|sort -k1,1nr|awk '$1>1{print $2}'|cut -d "." -f1|fgrep -f - analysis/src/dge_functions_data/ENSG_Genename_mapping.txt|cut -f12|grep -P "\S"|l >$i.amp.up;done
for i in DESeq2.logCPM5.exon_intron.*.all_reads.cbind_tables;do cat $i|perl -nale 'if((abs($F[3])>1 && $F[6]<0.05) || (abs($F[9])>1 && $F[12]<0.05)){if($F[3]*$F[9]>0 && abs($F[3]) > abs($F[9]) && $F[3] <0){print}}' |cut -f2|sort|uniq -c|sort -k1,1nr|awk '$1>1{print $2}'|head -1000|cut -d "." -f1|fgrep -f - analysis/src/dge_functions_data/ENSG_Genename_mapping.txt|cut -f12|grep -P "\S"|l >$i.amp.down;done
for i in DESeq2.logCPM5.exon_intron.*.all_reads.cbind_tables;do cat $i|perl -nale 'if((abs($F[3])>1 && $F[6]<0.05) || (abs($F[9])>1 && $F[12]<0.05)){if($F[3]*$F[9]>0 && abs($F[3]) < abs($F[9]) && $F[3] >0){print}}' |cut -f2|sort|uniq -c|sort -k1,1nr|awk '$1>1{print $2}'|head -1000|cut -d "." -f1|fgrep -f - analysis/src/dge_functions_data/ENSG_Genename_mapping.txt|cut -f12|grep -P "\S"|l >$i.red.up;done
for i in DESeq2.logCPM5.exon_intron.*.all_reads.cbind_tables;do cat $i|perl -nale 'if((abs($F[3])>1 && $F[6]<0.05) || (abs($F[9])>1 && $F[12]<0.05)){if($F[3]*$F[9]>0 && abs($F[3]) < abs($F[9]) && $F[3] <0){print}}' |cut -f2|sort|uniq -c|sort -k1,1nr|awk '$1>1{print $2}'|head -1000|cut -d "." -f1|fgrep -f - analysis/src/dge_functions_data/ENSG_Genename_mapping.txt|cut -f12|grep -P "\S"|l >$i.red.down;done
#extract data for delta fold change and expression correlation analysis
i=DESeq2.logCPM5.exon_intron.sra_tissues.all_reads.cbind_tables;cut -f-2,4,8,10,14 $i|perl -pe 's/\/figures\S+//'|l >$i.4deltaFC_exp_cor
i=DESeq2.logCPM5.exon_intron.gtex_tissue.all_reads.cbind_tables;cut -f-2,4,8,10,14 $i|cut -d "/" -f2-|l >$i.4deltaFC_exp_cor
i=DESeq2.logCPM5.exon_intron.abnormal_conditions.all_reads.cbind_tables;cut -f-2,4,8,10,14 $i|perl -pe 's/\/figures\S+//'|l >$i.4deltaFC_exp_cor
