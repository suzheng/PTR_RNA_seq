out_dir="analysis/results/100_samples_each_tissue/DGE_analysis"
out_dir=config["out_dir"]
rule dge_analysis:
    input:
        "merged_counts/{group1}.DONE",
        "merged_counts/{group2}.DONE"
    output:
    	"DGE_out/figures_tables.{group1}.{group2}.out/DESeq2.exon_counts_subsampled.DESeq2.{group1}.vs.{group2}.txt"
    log:
        "logs/dge_analysis/DGE_out.{group1}.{group2}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    container:
        "$SI/exon_intron_dge_r_packages.sif"
    shell:
        "Rscript analysis/src/dge_functions_data/GTEx_all_reads.DGE.R {wildcards.group1} {wildcards.group2} " + out_dir ##External files have to be in 'scripts' folder
