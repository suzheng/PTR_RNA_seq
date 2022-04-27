bam=$1
out_dir_full_path="$2"
ba=`basename $bam`
module add singularity
mkdir -p $out_dir_full_path
cd $out_dir_full_path 
indir=`dirname $bam`
singularity exec -B /srv/scratch/oateslab/share/data/hg38/ -B $indir -B $out_dir_full_path $SI/featurecounts.sif featureCounts -a /srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.exonsOnly.clean.gff3 -o $ba.exon -F GTF -t exon --ignoreDup -p -J --minOverlap 10 $bam;
singularity exec -B /srv/scratch/oateslab/share/data/hg38/ -B $indir -B $out_dir_full_path $SI/featurecounts.sif featureCounts -a /srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.intronsOnly.exonsSubtracted.gff3.gz -o $ba.intron -F GTF -t intron --ignoreDup -p -J --minOverlap 10 $bam;
#for feature level (i.e. each exon, not used)
singularity exec -B /srv/scratch/oateslab/share/data/hg38/ -B $indir -B $out_dir_full_path $SI/featurecounts.sif featureCounts -a /srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.exonsOnly.clean.gff3 -o $ba.exon.fl -O -T 4 -f -F GTF -t exon --ignoreDup -p -J --minOverlap 10 $bam;
singularity exec -B /srv/scratch/oateslab/share/data/hg38/ -B $indir -B $out_dir_full_path $SI/featurecounts.sif featureCounts -a /srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.intronsOnly.exonsSubtracted.gff3.gz -o $ba.intron.fl -O -T 4 -f -F GTF -t intron --ignoreDup -p -J --minOverlap 10 $bam;
