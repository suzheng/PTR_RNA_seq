#get the species list
ls /srv/scratch/oateslab/share/data/non_human_species/anno/*|cut -d "/" -f9|cut -d "." -f1-2|sort|uniq -c|awk '$1==3{print $2}'|l >species_list
#generate the annotation files, and STAR reference index files
for i in `cat species_list`;do echo "sh analysis/src/gff_gen_star_index_building/gff_gen_star_index_building.sh /srv/scratch/oateslab/share/data/non_human_species/anno/$i.*gff3.gz /srv/scratch/oateslab/share/data/non_human_species/anno/$i.*gtf /srv/scratch/oateslab/share/data/non_human_species/anno/$i.*.fa /srv/scratch/oateslab/share/data/non_human_species/$i/$i /srv/scratch/oateslab/share/data/non_human_species/$i"|cat template.pbs - >pbs/$i.pbs;qsub pbs/$i.pbs;done
#get the mapping info between gene symbols and Ensembl IDs
for i in `cat species_list`;do less /srv/scratch/oateslab/share/data/non_human_species/anno/$i.*gff3.gz |grep -P "\tgene\t"|perl -nale 'my $id="Unknown";/gene_id=([^\;]+);/ and $id=$1; my $name=$id;/Name=([^\;]+);/ and $name=$1;print "$id\t$name\t$_"'|cut -f-2 >/srv/scratch/oateslab/share/data/non_human_species/$i/$i.ENSG_Genename_mapping.txt;done
