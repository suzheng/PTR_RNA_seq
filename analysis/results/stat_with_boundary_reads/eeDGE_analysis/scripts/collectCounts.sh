prefix_list=$1
out_prefix=$2
for feature_type in exon intron; do
		
	first=`cat $prefix_list |awk -v feature_type=$feature_type '{print $1"."feature_type}'|head -1`;
	chmod 750 $first
	cp $first $out_prefix.$feature_type.A;
	for i in `cat $prefix_list |awk -v feature_type=$feature_type 'NR>1{print $1"."feature_type}'`;do
		cut -f7 $i |paste $out_prefix.$feature_type.A - >$out_prefix.$feature_type.B;
		mv $out_prefix.$feature_type.B $out_prefix.$feature_type.A;
	done;
	cat $out_prefix.$feature_type.A |grep -v featureCounts|perl -pe 's/\S+\///g;s/.Aligned.\S+.bam//g;' >$out_prefix.$feature_type.merged.txt;
	#gzip -f $out_prefix.$group.$feature_type.merged.txt;
	rm $out_prefix.$feature_type.A
done
touch $out_prefix.DONE
