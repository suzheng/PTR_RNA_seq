bam=$1
out_dir_full_path=$2
ba=`basename $bam`
module add singularity
cd $out_dir_full_path 
dir=analysis/results/100_samples_each_tissue/
singularity exec -B $dir -B /srv/scratch/oateslab/share/data/hg38/ -B /srv/scratch/oateslab/rawData/2021/ -B $out_dir_full_path $SI/featurecounts.sif featureCounts -a /srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.exonsOnly.clean.gff3 -o $ba.exon -F GTF -t exon --ignoreDup -p -J --minOverlap 10 $bam && \
singularity exec -B $dir -B /srv/scratch/oateslab/share/data/hg38/ -B /srv/scratch/oateslab/rawData/2021/ -B $out_dir_full_path $SI/featurecounts.sif featureCounts -a /srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.intronsOnly.exonsSubtracted.gff3.gz -o $ba.intron -F GTF -t intron --ignoreDup -p -J --minOverlap 10 $bam && \
rm $bam $bam.bai 
