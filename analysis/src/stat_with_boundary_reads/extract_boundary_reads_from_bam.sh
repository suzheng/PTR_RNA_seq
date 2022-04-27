#input bam
i=$1
#/srv/scratch/oateslab/share/data/hg38/Gencode
#l gencode.v38.GRCh38.exonsOnly.clean.gff3|perl -nale 'print "$F[0]\t$F[3]\t$F[3]\n$F[0]\t$F[4]\t$F[4]"'|sort -k1,1 -k2,2n|uniq >gencode.v38.GRCh38.exonsOnly.clean.gff3.exonBoundary.bed
exon_boundary_bed=/srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.exonsOnly.clean.gff3.exonBoundary.bed
#l gencode.v38.GRCh38.intronsOnly.exonsSubtracted.gff3.gz |perl -nale 'my $overhang=10;my $s=$F[3]-$overhang;my $e=$F[4]+$overhang;print "$F[0]\t$s\t$s\n$F[0]\t$e\t$e"'|l >gencode.v38.GRCh38.intronsOnly.exonsSubtracted.gff3.intronBoundary.oh10.bed
intron_boundary_bed=/srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.intronsOnly.exonsSubtracted.gff3.intronBoundary.oh10.bed

#get the basename
#ba=`basename $i .bam`;
#make sure there is 10bp overhang, with 10M or more at either end of xxxN in Cigar string
bedtools intersect -a $i -b $exon_boundary_bed|samtools view -h -| perl -F"\t" -nale 'if($F[0]=~/^\@/){print}elsif($F[5]=~/(\d+)M\d+N(\d+)M/ and $1>=10 and $2>=10){print}' | samtools view -S -b - >$i.exonBoundary.bam

#get intron boundary reads
bedtools intersect -a $i -b $intron_boundary_bed|samtools view -h -| \
awk -F "\t" '$1~"@" || $6!~"N"'|samtools view -S -b - >$i.intronBoundary.bam

samtools index $i.exonBoundary.bam
samtools index $i.intronBoundary.bam
