#!/bin/bash
module load bowtie2/2.3.1-fasrc01
export BOWTIE2_INDEXES=../references/cont_index/
bowtie2 --very-sensitive --no-hd -p $SLURM_NTASKS -U $infile -x cont_index -S $sample.sam
cat $sample.sam | awk '{if ($3 ~ /Homo/) {print "Homo"} else if ($3 ~ "*" ) {print "none"} else {gsub(/_[0-9]+/,"",$3); print $3}}' | sort | uniq -c | awk -v output=$sample 'BEGIN {split("Anchylorhynchus Andranthobius Celetes_impar Microstrates_bondari Microstrates_ypsilon Homo none",names); for (i in names) {counts[names[i]] = 0}} {counts[$2] += $1} END {for (i=1; i<=length(names); i++) {output=output "\t" counts[names[i]]} print output}' >> counts.txt
rm $sample.sam
