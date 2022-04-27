cp -pr ../SRA_samples/scripts ./
cp -pr ../SRA_samples/template.pbs ./
#generate the config files for Snakemake template
l /srv/scratch/oateslab/rawData/2021/SRA_RNA_seq/E-MTAB-7247/meta.txt |cut -f12|grep -v FastQ|perl -pe 's/.fastq.gz//'|perl -nale 'my $yy="Unknown";/Human/i and $yy="Homo_sapiens.GRCh38.readCount.yaml";/mouse/i and $yy="Mus_musculus.GRCm39.readCount.yaml";/Pan_/i and $yy="Pan_troglodytes.Pan_tro_3.readCount.yaml";/macaque/i and $yy="Macaca_mulatta.Mmul_10.readCount.yaml";/chicken/i and $yy="Gallus_gallus.GRCg6a.readCount.yaml";/opossum/ and $yy="Monodelphis_domestica.ASM229v1.readCount.yaml";/platypus/ and $yy="Ornithorhynchus_anatinus.mOrnAna1.readCount.yaml";print "$_\tconfig/$yy"'|l >E-MTAB-7247.samples.configs.txt
l /srv/scratch/oateslab/rawData/2021/SRA_RNA_seq/SRP028336/meta.txt |cut -f12|grep -v FastQ|perl -pe 's/.fastq.gz//'|perl -nale 'my $yy="Unknown";/Homo/i and $yy="Homo_sapiens.GRCh38.readCount.yaml";/musculus/i and $yy="Mus_musculus.GRCm39.readCount.yaml";/Pan_/i and $yy="Pan_troglodytes.Pan_tro_3.readCount.yaml";/Macaca/i and $yy="Macaca_mulatta.Mmul_10.readCount.yaml";print "$_\tconfig/$yy"'|l >SRP028336.samples.configs.txt
#Perform data download, alignment and quantification
cat SRP028336.samples.configs.txt E-MTAB-7247.samples.configs.txt|awk '{print "sh analyze_one_SE_sample.sh "$1" "$2" 1>logs/"$1".log 2>logs/"$1".err"}'|l
