#!/bin/sh

## To Run batch jobs for multiple samples.
## R1_Files.txt contains full paths of Read 1 for Sample.
## R2_Files.txt contains full paths of Read 2 for Sample.

while IFS= read -r p && IFS= read -r q <&3; do
sbatch --ntasks=1 --time=72:00:00 --mem=50000 --nodes=1 --ntasks-per-node=16 BulkATAC_Processing.sh -a $p -b $q -o /path/to/Output
done <R1_Files.txt 3<R2_Files.txt


