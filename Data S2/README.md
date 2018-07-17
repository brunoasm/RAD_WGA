# Code for Medeiros & Farrell. 2018. PeerJ

This folder will be made available as a repository at the first author's [github page] (https://github.com/brunoasm) upon publication.

## Contents

1. **wga_comparisons/**

 This folder contains scripts that calculate a number of statistics from assemblies, as well as pairwise comparisons.
It also contains ipyrad output in the form of `*.loci` files.


2. **contamination/**

  + **references/**

  This folder contains scripts that generate a bowtie database from samples not subjected to MDA and the human genome.

  + **reads/**

  Scripts that use bowtie to match raw reads.

  + **clusters/**

  Scripts that use bowtie to match *clusters*, i. e., loci assembled for each library.

  + **final_dataset/**

  Scripts that use bowtie to match loci kept in the final dataset. 
