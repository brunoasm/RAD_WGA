#!/bin/bash
#SBATCH -n 8 #Number of cores 
#SBATCH -N 1 #All cores in same machine
#SBATCH -p serial_requeue,interact #Partition to submit to 
#SBATCH --mem=8000
#SBATCH -J reference
#SBATCH -o make_reference.%A.out
#SBATCH -t 0-05:00:00


module load cd-hit/4.6.4-fasrc02 bowtie2/2.3.1-fasrc01

while read taxa;
do
taxon=`echo $taxa | cut -d , -f 1`
seqs=`echo $taxa | cut -d , -f 2`

find ../../wga_comparisons -name $taxon.loci | xargs -I {} awk -v seqs="$seqs" -v taxon="$taxon" 'BEGIN {split(seqs, valuesAsValues); for (i in valuesAsValues) valuesAsKeys[valuesAsValues[i]] = ""; OFS="\n"}($1 in valuesAsKeys) {gsub("-","",$2); print ">"taxon"_"NR,$2}' {} > $taxon.fasta

cd-hit -T $SLURM_NTASKS -i $taxon.fasta -o $taxon.derep.fasta -T 1 -c 0.95

done < samples_no_wga.txt 

rm *.clstr

cat *.derep.fasta > all_seqs
rm *.fasta
mv all_seqs all_seqs.fasta

#run once to get human genome:
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.25_GRCh38.p10/GCA_000001405.25_GRCh38.p10_genomic.fna.gz
#zcat GCA_000001405.25_GRCh38.p10_genomic.fna.gz | sed 's/ /_/g' | sed -r 's/[acgt]{1}/N/g' |gzip > human_genome.fasta.gz

zcat human_genome.fasta.gz >> all_seqs.fasta

mkdir -p cont_index
cd cont_index

bowtie2-build --threads $SLURM_NTASKS ../all_seqs.fasta cont_index

