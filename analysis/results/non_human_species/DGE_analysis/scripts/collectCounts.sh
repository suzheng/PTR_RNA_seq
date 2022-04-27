dir=analysis/results/non_human_species/DGE_analysis;
p=$1;
fc_dir=analysis/results/non_human_species/count_reads/;
out_prefix=$dir/$p/$p
mkdir -p $dir/$p
for group in `cat $dir/meta*/$p.samples.tsv|grep -v -w condition|cut -f2|sort|uniq`; do
	for feature_type in exon intron; do
		grep -w $group $dir/meta*/$p.samples.tsv|cut -f1 >$out_prefix.$group.samples.tmp
		first=`ls $fc_dir/*.$feature_type|fgrep -w -f $out_prefix.$group.samples.tmp - |head -1`;
		chmod 750 $first
		cp $first $out_prefix.$feature_type.A;
		for i in `ls $fc_dir/*.$feature_type|fgrep -w -f $out_prefix.$group.samples.tmp - |awk 'NR>1'`;do
			cut -f7 $i |paste $out_prefix.$feature_type.A - >$out_prefix.$feature_type.B;
			mv $out_prefix.$feature_type.B $out_prefix.$feature_type.A;
		done;
		cat $out_prefix.$feature_type.A |grep -v featureCounts|perl -pe 's/\S+\///g;s/.Aligned.\S+.bam//g;' >$out_prefix.$group.$feature_type.merged.txt;
		#gzip -f $out_prefix.$group.$feature_type.merged.txt;
		rm $out_prefix.$feature_type.A $out_prefix.$group.samples.tmp
	done
done
touch $out_prefix.DONE
