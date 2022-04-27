#add intron to the gff file
singularity exec -B `pwd` $SI/genometools.sif gt gff3 -retainids -addintrons gencode.v38.annotation.gff3.gz >gencode.v38.annotation.addedIntrons.gff3
#check the line that cause an error and get ride of that line
l gencode.v38.annotation.gff3.gz|head -1033189|tail -n 1 >err_line.1033189.txt
#re-run the add intron command
singularity exec -B `pwd` $SI/genometools.sif gt gff3 -retainids -addintrons -tidy gencode.v38.annotation.gff3.gz  >gencode.v38.annotation.addedIntrons.gff3&
#add gene IDs to the added introns
perl add_gene_id_to_introns.pl gencode.v38.annotation.addedIntrons.gff3 1>O 2>E&
#extract exon and intron separately
i=gencode.v38.annotation.addedIntrons.gff3.addedParent.gff3.gz;l $i|awk '/#/ || $3=="exon"'|gzip -c >gencode.v38.GRCh38.exonsOnly.gff3.gz
i=gencode.v38.annotation.addedIntrons.gff3.addedParent.gff3.gz;l $i|awk '/#/ || $3=="intron"'|gzip -c >gencode.v38.GRCh38.intronsOnly.gff3.gz
#extract overlapping exons from introns
bedtools subtract -a gencode.v38.GRCh38.intronsOnly.gff3.gz -b gencode.v38.GRCh38.exonsOnly.gff3.gz|l >gencode.v38.GRCh38.intronsOnly.exonsSubtracted.gff3.gz
#extract overlapping introns from exons
bedtools subtract -a gencode.v38.GRCh38.exonsOnly.gff3.gz -b gencode.v38.GRCh38.intronsOnly.gff3.gz|l >gencode.v38.GRCh38.exonsOnly.intronsSubtracted.gff3.gz
l gencode.v38.GRCh38.exonsOnly.gff3.gz|grep -vP "^\#" >gencode.v38.GRCh38.exonsOnly.clean.gff3

#extract the mapping information between gene symbols and EnsemblIDs
l gencode.v38.annotation.gff3.gz|grep -w gene|perl -pe 's/;/\t/g;s/\t\w+=/\t/g;'|l >ENSG_Genename_mapping.txt
#convert to bed format
l gencode.v38.GRCh38.exonsOnly.clean.gff3|perl -nale 'print "$F[0]\t$F[3]\t$F[3]\n$F[0]\t$F[4]\t$F[4]"'|sort -k1,1 -k2,2n|uniq >gencode.v38.GRCh38.exonsOnly.clean.gff3.exonBoundary.bed&
#generate bed file with 10bp overhang
l gencode.v38.GRCh38.intronsOnly.exonsSubtracted.gff3.gz |perl -nale 'my $overhang=10;my $s=$F[3]-$overhang;my $e=$F[4]+$overhang;print "$F[0]\t$s\t$s\n$F[0]\t$e\t$e"'|l >gencode.v38.GRCh38.intronsOnly.exonsSubtracted.gff3.intronBoundary.oh10.bed

#generate the gff files for the 5' end and 3' end exons and introns.
#some genes are in minus strand, so to multiple -1 to their coordinates for sorting together with plus strand genes
#sort exons by gene names, then by coordinates
#output the first two exons of the gene, note that different transcripts of a same gene may have same exons, so this is not equal to output first two lines of the gene.
l /srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.exonsOnly.clean.gff3|perl -F"\t" -nale '$F[8]=~/gene_id=([^\;]+);/ and my $gene=$1;my $pos=$F[3];if($F[6] eq "-"){$pos=(-1)*$pos;}print "$gene\t$pos\t$_"'|sort -k1,1 -k2,2n|awk 'BEGIN{n=0;preGene="";preExon=""}{if($1!=preGene){n=0;preGene=$1}if($2!=preExon){n=n+1;preExon=$2}if(n<=2){print $0}}'|cut -f3-|sort -k1,1 -k4,4n|l >/srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.exonsOnly.clean.5prime.gff3
l /srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.exonsOnly.clean.gff3|perl -F"\t" -nale '$F[8]=~/gene_id=([^\;]+);/ and my $gene=$1;my $pos=$F[3];if($F[6] eq "-"){$pos=(-1)*$pos;}print "$gene\t$pos\t$_"'|sort -k1,1 -k2,2nr|awk 'BEGIN{n=0;preGene="";preExon=""}{if($1!=preGene){n=0;preGene=$1}if($2!=preExon){n=n+1;preExon=$2}if(n<=2){print $0}}'|cut -f3-|sort -k1,1 -k4,4n|l >/srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.exonsOnly.clean.3prime.gff3&
i=gencode.v38.GRCh38.intronsOnly.exonsSubtracted.gff3.gz;l $i|perl -F"\t" -nale '$F[8]=~/gene_id=([^\;]+)/ and my $gene=$1;my $pos=$F[3];if($F[6] eq "-"){$pos=(-1)*$pos;}print "$gene\t$pos\t$_"'|sort -k1,1 -k2,2n|awk 'BEGIN{n=0;preGene="";preExon=""}{if($1!=preGene){n=0;preGene=$1}if($2!=preExon){n=n+1;preExon=$2}if(n<=2){print $0}}'|cut -f3-|sort -k1,1 -k4,4n >gencode.v38.GRCh38.intronsOnly.exonsSubtracted.5prime.gff3&
i=gencode.v38.GRCh38.intronsOnly.exonsSubtracted.gff3.gz;l $i|perl -F"\t" -nale '$F[8]=~/gene_id=([^\;]+)/ and my $gene=$1;my $pos=$F[3];if($F[6] eq "-"){$pos=(-1)*$pos;}print "$gene\t$pos\t$_"'|sort -k1,1 -k2,2nr|awk 'BEGIN{n=0;preGene="";preExon=""}{if($1!=preGene){n=0;preGene=$1}if($2!=preExon){n=n+1;preExon=$2}if(n<=2){print $0}}'|cut -f3-|sort -k1,1 -k4,4n >gencode.v38.GRCh38.intronsOnly.exonsSubtracted.3prime.gff3&
