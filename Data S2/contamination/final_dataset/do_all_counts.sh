#!/bin/bash
#first, initialize output table with column names
echo -e "sample\tAnchylorhynchus\tAndranthobius\tCeletes_impar\tMicrostrates_bondari\tMicrostrates_ypsilon\tHomo\tnone" > counts.txt

#now, call slurm to:
#1 - produce a fasta file for each sample from *.loci file
#2 - use bowtie to make an alignment
#3 - use awk to get the third column of bowtie output (i. e. name of target sequence)
#4 - use sort and uniq to count the number of reads to each target taxon
#5 - use awk again to reformat this output to a single line and append to counts.txt
#6 - remove bowtie alignment



find ../../*_outfiles -name '*.loci' -exec cat {} + | cut -d ' ' -f 1 | sort | uniq | grep BdM | while read sample
do
    export sample
    echo processing ${sample} 
    sbatch -N 1 -n 8 -p serial_requeue -o $sample.%A.out --export=ALL --mem=5000 -J r_$sample -t 0-02:00:00 count.sh
done
