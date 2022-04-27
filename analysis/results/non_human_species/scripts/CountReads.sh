bam=$1
out_dir_full_path="$2"
exon_gff=$3
intron_gff=$4
mount_dir=$5
#mount_dir=/srv/scratch/oateslab/share/data/

ba=`basename $bam`
module add singularity
mkdir -p $out_dir_full_path
cd $out_dir_full_path 
indir=`dirname $bam`
singularity exec -B $mount_dir -B $indir -B $out_dir_full_path $SI/featurecounts.sif featureCounts -a $exon_gff -o $ba.exon -F GTF -t exon --ignoreDup -p -J --minOverlap 10 $bam;
singularity exec -B $mount_dir -B $indir -B $out_dir_full_path $SI/featurecounts.sif featureCounts -a $intron_gff -o $ba.intron -F GTF -t intron --ignoreDup -p -J --minOverlap 10 $bam;
singularity exec -B $mount_dir -B $indir -B $out_dir_full_path $SI/featurecounts.sif featureCounts -a $exon_gff -o $ba.exon.fl -O -T 4 -f -F GTF -t exon --ignoreDup -p -J --minOverlap 10 $bam;
singularity exec -B $mount_dir -B $indir -B $out_dir_full_path $SI/featurecounts.sif featureCounts -a $intron_gff -o $ba.intron.fl -O -T 4 -f -F GTF -t intron --ignoreDup -p -J --minOverlap 10 $bam;
