infile=$1
awk -F"\t" -v infile=$infile '{f="Mapping Rate";\
res="Unknown";val="NaN";if($1==f){val=$2;res="FAILED";\
if($2>0.7){res="PASS"}
	print f"\t"res"\t"val"\t"infile};}' $infile

awk -F"\t" -v infile=$infile '{f="Base Mismatch";\
res="Unknown";val="NaN";if($1==f){val=$2;res="FAILED";\
if($2<0.02){res="PASS"}
	print f"\t"res"\t"val"\t"infile};}' $infile

awk -F"\t" -v infile=$infile '{f="High Quality Rate";\
res="Unknown";val="NaN";if($1==f){val=$2;res="FAILED";\
if($2>0.6){res="PASS"}
	print f"\t"res"\t"val"\t"infile};}' $infile

awk -F"\t" -v infile=$infile '{f="Exonic Rate";\
res="Unknown";val="NaN";if($1==f){val=$2;res="FAILED";\
if($2>0.3){res="PASS"}
	print f"\t"res"\t"val"\t"infile};}' $infile

awk -F"\t" -v infile=$infile '{f="Ambiguous Alignment Rate";\
res="Unknown";val="NaN";if($1==f){val=$2;res="FAILED";\
if($2<0.1){res="PASS"}
	print f"\t"res"\t"val"\t"infile};}' $infile

awk -F"\t" -v infile=$infile '{f="rRNA Rate";\
res="Unknown";val="NaN";if($1==f){val=$2;res="FAILED";\
if($2<0.4){res="PASS"}
	print f"\t"res"\t"val"\t"infile};}' $infile


awk -F"\t" -v infile=$infile '{f="Avg. Splits per Read";\
res="Unknown";val="NaN";if($1==f){val=$2;res="FAILED";\
if($2<0.4){res="PASS"}
	print f"\t"res"\t"val"\t"infile};}' $infile

awk -F"\t" -v infile=$infile '{f="Genes Detected";\
res="Unknown";val="NaN";if($1==f){val=$2;res="FAILED";\
if($2<30000 && $2>5000){res="PASS"}
	print f"\t"res"\t"val"\t"infile};}' $infile
