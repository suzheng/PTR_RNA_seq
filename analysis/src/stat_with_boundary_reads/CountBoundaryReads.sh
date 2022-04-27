bam=$1
out_dir_full_path=$2
SD="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
ba=`basename $bam`
module add singularity
export SINGULARITY_CACHEDIR=$TMPDIR/SINGULARITY_CACHEDIR/
module add samtools bedtools

cd $out_dir_full_path
sh $SD/extract_boundary_reads_from_bam.sh $bam

singularity exec -B /srv/scratch/oateslab/share/data/hg38/ -B /srv/scratch/oateslab/rawData/2021/ -B analysis/results/SRA_samples/star/ -B $out_dir_full_path $SI/featurecounts.sif featureCounts -a /srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.exonsOnly.clean.gff3 -o $bam.exon -F GTF -t exon --ignoreDup -p -J --splitOnly --minOverlap 10 $bam.exonBoundary.bam;
singularity exec -B /srv/scratch/oateslab/share/data/hg38/ -B /srv/scratch/oateslab/rawData/2021/ -B analysis/results/SRA_samples/star/ -B $out_dir_full_path $SI/featurecounts.sif featureCounts -a /srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.intronsOnly.exonsSubtracted.gff3.gz -o $bam.intron -F GTF -t intron --ignoreDup -p -J --nonSplitOnly --minOverlap 10 $bam.intronBoundary.bam;
