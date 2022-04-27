mkdir pbs cat_pbs
for sample in `cat sample_names.txt`;do echo "snakemake --cores all --use-singularity -s analysis/src/stat_with_boundary_reads/Snakefile_analyze_single_SRA_sample.py pbs/${sample}.pbs";done|l
cat pbs/*|split -l 120 --additional-suffix=.pbs - cat_pbs/split 
for i in cat_pbs/*;do qsub $i;done
