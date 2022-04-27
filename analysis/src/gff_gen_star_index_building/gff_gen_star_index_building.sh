#input gff file
ingff=$1
sjdbGTFfile=$2
genomeFastaFiles=$3
#exon intron annot file output prefix
prefix=$4
#star index output directory
genomeDir=$5
genomeSAindexNbases=14

if [ $# -eq 6 ]
  then
    genomeSAindexNbases=$6
fi

module add singularity bedtools star

mkdir $genomeDir
#add introns to gff
singularity exec -B `pwd` -B /srv/scratch/oateslab/share/data/ $SI/genometools.sif gt gff3 -retainids -addintrons -tidy $ingff  >$prefix.addedIntrons.gff3
#add gene IDs to the added introns
perl /srv/scratch/oateslab/share/data/hg38/Gencode/add_gene_id_to_introns.pl $prefix.addedIntrons.gff3 $sjdbGTFfile $prefix.addedIntrons.gff3.addedParent.gff3.gz
#extract intron and exon separately
i=$prefix.addedIntrons.gff3.addedParent.gff3.gz;less $i|awk '/#/ || $3=="exon"'|gzip -c >$prefix.exonsOnly.gff3.gz
i=$prefix.addedIntrons.gff3.addedParent.gff3.gz;less $i|awk '/#/ || $3=="intron"'|gzip -c >$prefix.intronsOnly.gff3.gz
less $prefix.exonsOnly.gff3.gz >$prefix.exonsOnly.gff3
#subtract overlapping exon from introns
/srv/scratch/oateslab/softwares/bedtools/bedtools subtract -a $prefix.intronsOnly.gff3.gz -b $prefix.exonsOnly.gff3 >$prefix.intronsOnly.exonsSubtracted.gff3.gz

#add gene ID to extracted exons
less $prefix.exonsOnly.gff3.gz|grep -vP "^\#" >$prefix.exonsOnly.tmp.gff3
perl /srv/scratch/oateslab/share/data/hg38/Gencode/add_gene_id_to_introns.pl $prefix.exonsOnly.tmp.gff3 $sjdbGTFfile $prefix.exonsOnly.clean.gff3

#extracted the mapping info between gene symbols and Ensembl IDs
less $ingff|grep -P "\tgene\t"|perl -nale 'my $id="Unknown";/gene_id=([^\;]+);/ and $id=$1; my $name=$id;/Name=([^\;]+);/ and $name=$1;print "$id\t$name\t$_"'|cut -f-2 >$prefix.ENSG_Genename_mapping.txt

#generate STAR index
rm -r $TMPDIR/star
STAR   --runMode genomeGenerate   --outTmpDir $TMPDIR/star  --runThreadN 22   --genomeDir $genomeDir   --genomeFastaFiles $genomeFastaFiles   --sjdbGTFfile $sjdbGTFfile   --sjdbOverhang 100 --limitGenomeGenerateRAM=45000000000 --genomeSAindexNbases $genomeSAindexNbases 
