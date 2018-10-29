#!/bin/bash

for i in 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 
do
  python3 contig_merge.py -i ./data/lambda.r2.no_error.fq -c $i -k 31 -r ./reports/report_lambda/report_cutoff$i
done
