#!/bin/bash
set -x

for i in 19 23 29 31 37 41  
do
   for j in -1 2 3 4 5 6 8 10 15   	
   do
	   python3 contig_merge.py -i ./data/lambda.r1.no_error.fq -k $i -c $j -o ./output/assembly$i.$j.fasta -r ./output/report$i.$j.txt >> ./output/output_merge1.txt
   done
done
