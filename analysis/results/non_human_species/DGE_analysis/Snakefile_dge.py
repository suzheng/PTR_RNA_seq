import glob
#genename_mapping_file = config["genename_mapping_file"]
genename_mapping_file = "analysis/results/non_human_species/DGE_analysis/scripts/ENSG_Genename_mapping.txt" 
out_dir = "analysis/results/non_human_species/DGE_analysis/"

def get_meta_file(wildcards):
    # code that returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    return glob.glob(out_dir + "/meta/" + wildcards.project + '*.samples.tsv')

rule collect_count:
    input:
        get_meta_file
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
        "Rscript analysis/src/dge_functions_data/non_human_species_all_reads.DGE.R {wildcards.project} {wildcards.group1} {wildcards.group2} " + out_dir + " " + genename_mapping_file##External files have to be in 'scripts' folder
