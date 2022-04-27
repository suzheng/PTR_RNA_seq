
#out_dir = "analysis/results/SRA_samples/DGE_analysis/"
out_dir=config["out_dir"]
meta_dir=config["meta_dir"]
rule collect_count:
    input:
        meta_dir + "/{project}.samples.tsv"
    output:
        "{project}/{project}.DONE"
    log:
        "logs/collect_count/{project}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    shell:
        "sh scripts/collectCounts.sh {wildcards.project}"



rule dge_analysis:
    input:
        "{project}/{project}.DONE"
    output:
    	"{project}/figures_tables.{group1}.{group2}.out/DESeq2.exon_counts_subsampled.DESeq2.{group1}.vs.{group2}.txt"
    log:
        "logs/dge_analysis/{project}.{group1}.{group2}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    container:
        "$SI/exon_intron_dge_r_packages.sif"
    shell:
        "Rscript analysis/src/dge_functions_data/SRA_all_reads.DGE.R {wildcards.project} {wildcards.group1} {wildcards.group2} " + out_dir ##External files have to be in 'scripts' folder
