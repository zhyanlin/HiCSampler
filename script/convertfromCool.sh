if test "$#" -ne 4; then
    echo "Please provide all arguments:"
    echo "bash convertfromCool.sh coolfile region resol output"
    exit
fi
input=$1
region=$2
resol=$3
output=$4
chrom=$(echo $region | cut -d ":" -f 1)
startbp=$(echo $region | cut -d ":" -f 2 | cut -d '-' -f1)
endbp=$(echo $region | cut -d ":" -f 2 | cut -d '-' -f2)

cooler dump $input --join -r $region | awk -v startbp=$startbp -v resol=$resol '{print $1"\t"($2-startbp)/resol"\t"$1"\t"($5-startbp)/resol"\t"$NF}' > ${output}.RC.tsv
scale=$( python $( dirname -- "${BASH_SOURCE[0]}" )/getScale.py $input $chrom)
cooler  dump $input -t bins  -c 'chrom,start,weight' --na-rep 999999 --join -H | awk -v scale=$scale -v chrom=$chrom -v startbp=$startbp -v endbp=$endbp -v resol=$resol '{if($1==chrom && $2>=startbp-resol && $2<endbp-resol)print 1/$NF/scale}' >  ${output}.bias.tsv
