for i in `cat /srv/scratch/oateslab/share/data/non_human_species/scripts/species_list`;do sed "s/SPECIES_NAME/$i/g" template.yaml.txt >$i.readCount.yaml;done
