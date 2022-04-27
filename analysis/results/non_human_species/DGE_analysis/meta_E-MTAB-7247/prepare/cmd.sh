l /srv/scratch/oateslab/rawData/2021/SRA_RNA_seq/E-MTAB-7247/meta.txt|awk '{print $NF}'|perl -pe 's/.fastq.gz//'|grep -v filename|l> all_species_sample_names
l all_species_sample_names|cut -d "_" -f1|sort|uniq|l >species_list
for sp in `cat species_list`;do n=`cat tissue_list|wc -l`;for i in $(seq 1 $n);do i2=`expr $i + 1`; for j in $(seq $i2 $n );do g1=`head -$i tissue_list|tail -n1`;g2=`head -$j tissue_list|tail -n1`;grep -P "$g1|$g2" all_species_sample_names|grep $sp |awk -F "_" '{print $0"\t"$2}'|cat header - >../E-MTAB-7247.${g1}_$g2.$sp.samples.tsv;done;done;done
for i in ../*.tsv;do n=`cat $i|cut -f2|sort|uniq|wc -l`;echo -ne "$i\t$n\n";done|awk '$2<3{print "rm "$1}'|sh
