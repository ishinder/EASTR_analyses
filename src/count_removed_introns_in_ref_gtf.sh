#!/bin/bash

if [ $# -ne 2 ]; then
  echo "Usage: $0 bed_file jsize"
  exit 1
fi

bed_file=$1
jsize=$2

njunc=$(wc -l $bed_file | awk '{print $1}')
genes=($(awk '{split($NF,a,";"); print a[1]}' $bed_file | sort | uniq))
ngenes=${#genes[@]}
ntxns=$(awk '{split($NF,a,";"); print length(a)-1}' $bed_file | awk '{sum += $1} END {print sum}')


#count introns less than or equal to jsize
count=0
while read line; do
  start=$(echo $line | awk '{print $2}')
  end=$(echo $line | awk '{print $3}') 
  intron_length=$((end - start))
  if [ $intron_length -le $jsize ]; then
    count=$((count + 1))
  fi
done < $bed_file

proportion=$(echo "scale=4; $count / $njunc" | bc)

echo "Number of removed introns = $njunc"
echo "Number of removed transcripts = ${ntxns}"
echo "Number of genes with removed transcripts = ${ngenes}"
echo "Number of introns less than or equal to $jsize: $count"
echo "Proportion of introns less than or equal to $jsize: $proportion"
