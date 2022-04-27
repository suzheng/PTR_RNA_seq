#quantify using only boundary reads
for sample in `cat ../100_samples_each_tissue/sample_names.txt`;do echo "snakemake --cores all --use-singularity -s analysis/src/stat_with_boundary_reads/Snakefile_analyze_single_sample.py out/${sample}/${sample}.Aligned.sortedByCoord.out.patched.md.bam.intron"|cat template.pbs - >pbs/$sample.pbs;done
#quantify using only reads at two ends of gene
for i in `cat ../100_samples_each_tissue/sample_names.txt`;do echo "sh analysis/src/stat_with_end_exons_reads//CountTwoEndsReads.sh out/$i/$i.Aligned.sortedByCoord.out.patched.md.bam analysis/results/stat_with_boundary_reads"|cat analysis/src/stat_with_end_exons_reads/template.pbs - >pbsee/$i.pbs;done
for i in pbsee/*pbs;do echo "qsub $i";done|sh
