#!/bin/bash
s=$1
dir=$PWD
mkdir -p $dir/out/$s 
if [[ -e $dir/out/$s/$s.download_done  &&  (! -f $dir/out/$s/$s.analysis_started) ]]
then
    echo "sh CountReads.sh $dir/out/$s/$s.Aligned.sortedByCoord.out.patched.md.bam $dir/out/$s/"|cat pbs_templ.pbs - >$dir/out/$s/$s.pbs && \
    qsub $dir/out/$s/$s.pbs && \
    touch $dir/out/$s/$s.analysis_started
else
    echo "$s not ready or analysis started" >log
fi

