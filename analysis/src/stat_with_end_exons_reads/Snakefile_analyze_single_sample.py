import os

#sd = os.path.realpath(__file__)
sd = "analysis/src/stat_with_end_exons_reads/"
#exon_boundary_bed="/srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.exonsOnly.clean.gff3.exonBoundary.bed"
#intron_boundary_bed="/srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.intronsOnly.exonsSubtracted.gff3.intronBoundary.oh10.bed"
rule extract_end_exon_reads_count_reads:
    input:
        "out/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam"
    output:
        "pbsee/{sample}.pbs"
    log:
        "logs/extract_end_exon_reads_count_reads/{sample}.log"
    threads:  # Samtools takes additional threads through its option -@
        1     # This value - 1 will be sent to -@
    shell:
        "mkdir pbsee;echo sh " + sd + "/CountTwoEndsReads.sh out/{wildcards.sample}/{wildcards.sample}.Aligned.sortedByCoord.out.patched.md.bam analysis/results/stat_with_boundary_reads|cat " + sd + "/template.pbs - >pbsee/{wildcards.sample}.pbs;qsub pbsee/{wildcards.sample}.pbs"
