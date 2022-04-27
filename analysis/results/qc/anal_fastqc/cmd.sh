#collect fastq QC summary, extract failed samples
awk '{print FILENAME"\t"$0}' ../fastqc_out/*/*/summary.txt|cut -d "/" -f3,5-|perl -pe 's/\/summary.txt//' >summary.txt.cat
l summary.txt.cat |grep "Per sequence quality scores"|grep FAIL|cut -f1|cut -d "." -f2|l >failed_samples
