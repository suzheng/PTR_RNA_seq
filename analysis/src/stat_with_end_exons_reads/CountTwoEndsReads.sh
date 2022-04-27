bam=$1
out_dir_full_path=$2
SD="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
ba=`basename $bam`
module add singularity
export SINGULARITY_CACHEDIR=$TMPDIR/SINGULARITY_CACHEDIR/
module add samtools bedtools

cd $out_dir_full_path
#sh $SD/extract_boundary_reads_from_bam.sh $bam

singularity exec -B /srv/scratch/oateslab/share/data/hg38/ -B /srv/scratch/oateslab/rawData/2021/ -B $out_dir_full_path $SI/featurecounts.sif featureCounts -a /srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.exonsOnly.clean.5prime.gff3 -o $bam.5prime.exon -F GTF -t exon --ignoreDup -p -J --minOverlap 10 $bam;
singularity exec -B /srv/scratch/oateslab/share/data/hg38/ -B /srv/scratch/oateslab/rawData/2021/ -B $out_dir_full_path $SI/featurecounts.sif featureCounts -a /srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.exonsOnly.clean.3prime.gff3 -o $bam.3prime.exon -F GTF -t exon --ignoreDup -p -J --minOverlap 10 $bam;
singularity exec -B /srv/scratch/oateslab/share/data/hg38/ -B /srv/scratch/oateslab/rawData/2021/ -B $out_dir_full_path $SI/featurecounts.sif featureCounts -a /srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.intronsOnly.exonsSubtracted.5prime.gff3 -o $bam.5prime.intron -F GTF -t intron --ignoreDup -p -J --minOverlap 10 $bam;
singularity exec -B /srv/scratch/oateslab/share/data/hg38/ -B /srv/scratch/oateslab/rawData/2021/ -B $out_dir_full_path $SI/featurecounts.sif featureCounts -a /srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.intronsOnly.exonsSubtracted.3prime.gff3 -o $bam.3prime.intron -F GTF -t intron --ignoreDup -p -J --minOverlap 10 $bam;
