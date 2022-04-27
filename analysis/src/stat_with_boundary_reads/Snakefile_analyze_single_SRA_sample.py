import os

#sd = os.path.realpath(__file__)
sd = "analysis/src/stat_with_boundary_reads/"
#exon_boundary_bed="/srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.exonsOnly.clean.gff3.exonBoundary.bed"
#intron_boundary_bed="/srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.intronsOnly.exonsSubtracted.gff3.intronBoundary.oh10.bed"

rule extract_boundary_reads_count_reads:
    input:
        "analysis/results/SRA_samples/star/{sample}.Aligned.sortedByCoord.out.md.bam"
    output:
        "pbs/{sample}.pbs"
    log:
        "logs/extract_boundary_reads_count_reads/{sample}.log"
    threads:  # Samtools takes additional threads through its option -@
        1     # This value - 1 will be sent to -@
    shell:
        "mkdir pbs;echo sh " + sd + "/CountBoundaryReads_for_SRA.sh analysis/results/SRA_samples/star/{wildcards.sample}.Aligned.sortedByCoord.out.md.bam analysis/results/boundary_reads_only_SRA/star|cat " + sd + "/template.pbs - >pbs/{wildcards.sample}.pbs;"
