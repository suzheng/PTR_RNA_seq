#!/usr/bin/perl -w
use strict;

my $in = shift;
my $in_gtf = shift;
my $out = shift||"$in.addedParent.gff3.gz";
my %m;
open I, "$in_gtf";
while(<I>){
        chomp;
        if(/gene_id \"(\w+)\";.*transcript_id \"(\w+)\";/){
                $m{$2}=$1;
        }
}
close I;
if($out =~ /\.gz$/){
	open O, "|gzip -c >$out";
}else{
	open O, ">$out";
}
open I, "$in";
while(<I>){
        chomp;
        my @F = split /\t/;
        if($F[$#F] =~ /Parent=transcript:(\w+)/){
                print O "$_;gene_id=$m{$1}\n";
        }else{
                print O "$_\n";
        }
}
close I;
close O;