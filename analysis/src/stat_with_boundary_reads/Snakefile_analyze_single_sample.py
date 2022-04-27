import os

#sd = os.path.realpath(__file__)
sd = "analysis/src/stat_with_boundary_reads/"
#exon_boundary_bed="/srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.exonsOnly.clean.gff3.exonBoundary.bed"
#intron_boundary_bed="/srv/scratch/oateslab/share/data/hg38/Gencode/gencode.v38.GRCh38.intronsOnly.exonsSubtracted.gff3.intronBoundary.oh10.bed"
#print(sd)
#rule download_one_sample:
#    input:
#        sd + "/fakeInput"
#    output:
#        "out/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam"
#    log:
#        "logs/download_one_sample/{sample}.log"
#    params:
#        "" # optional params string
#    threads:  # Samtools takes additional threads through its option -@
#        1     # This value - 1 will be sent to -@
#    shell:
#        "sh " + sd + "/download_one_sample.sh {wildcards.sample}"

rule extract_boundary_reads_count_reads:
    input:
        "out/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam"
    output:
        "out/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam.intron"
    log:
        "logs/extract_boundary_reads_count_reads/{sample}.log"
    threads:  # Samtools takes additional threads through its option -@
        1     # This value - 1 will be sent to -@
    shell:
        "sh " + sd + "/CountBoundaryReads.sh out/{wildcards.sample}/{wildcards.sample}.Aligned.sortedByCoord.out.patched.md.bam analysis/results/stat_with_boundary_reads"
