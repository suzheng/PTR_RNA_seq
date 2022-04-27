#The python script was downloaded from https://github.com/broadinstitute/gtex-pipeline
for i in /srv/scratch/oateslab/share/data/non_human_species/anno/*gtf;do ba=`basename $i`;echo "python3 collapse_annotation.py $i $ba";done|sh
