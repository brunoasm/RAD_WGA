#!/bin/bash
#first, initialize output table with column names
echo -e "sample\tAnchylorhynchus\tAndranthobius\tCeletes_impar\tMicrostrates_bondari\tMicrostrates_ypsilon\tHomo\tnone" > counts.txt

#now, call slurm to:
#1 - use bowtie to make an alignment
#2 - use awk to get the third column of bowtie output (i. e. name of target sequence)
#3 - use sort and uniq to count the number of reads to each target taxon
#4 - use awk again to reformat this output to a single line and append to counts.txt
#5 - remove bowtie alignment

for infile in ../../all_weevils_consens/*.consens.gz
do
export infile
export sample=`basename $infile | cut -d . -f 1`

echo processing ${sample}, file $infile 
sbatch -N 1 -n 8 -p serial_requeue -o $sample.%A.out --export=ALL --mem=5000 -J r_$sample -t 7-00:00:00 count.sh
done
