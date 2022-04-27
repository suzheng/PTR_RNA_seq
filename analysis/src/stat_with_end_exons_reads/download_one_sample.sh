s=$1
dir=$PWD
mkdir -p $dir/out/$s 
gen3-client download-multiple --no-prompt --profile=AnVIL --manifest=$dir/jsons/GTEx.$s.json --download-path=$dir/out/$s --protocol=s3 --skip-completed --numparallel 1

