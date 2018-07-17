#!/usr/bin/env python
# -*- coding: utf-8 -*-

#this script calculates a number of statistics from pyrad results
#written by B. Medeiros Mar-2017


import re, pandas, numpy as np, argparse, os
from joblib import Parallel, delayed
#import multiprocessing

## this exception will be raised when samples have no loci in common
class OverlapError(Exception):
    pass


##function to parse a *.loci file
##*.alleles file not available anymore, rewrote everything to work with loci file!
def loci_file_parser(filename):
    def sample_parser(samplename): #this function obtains sample id, pool from sample name
        sampleid = samplename[0:7]
        pool = samplename[11:13]
        return(sampleid, pool)
    def seq_parser(line): #this function parses a line containing sequences, returning sample name and sequence
        line = re.sub('[>\n]','',line)
        (samplename, sequence) = re.split('\s*',line)
        sequence = sequence.upper()
        return(samplename, sequence)
    def separator_parser(line): #this parses separator lines. currently, only returns number of locus
        line = re.sub('[>\n]','',line)
        return(re.split('\|',line)[1])
    with open(filename, 'r') as infile:
        seq_dict = {}
        temp_dict = {}
        for line in infile:
            if '//' not in line: #while parsing a single RAD locus, create a temporary dictionary
                if line != '\n': #skip empty lines
                    (samplename, sequence) = seq_parser(line)
                    (sampleid, pool) = sample_parser(samplename)
                    temp_dict[(sampleid, pool)] = sequence
            else: #when finished locus, add locus id and merge to main dictionary
                locus_id = separator_parser(line)
                temp_dict = {(sampleid, pool, locus_id):sequence for [(sampleid, pool),sequence] in temp_dict.iteritems()}
                seq_dict.update(temp_dict)
                temp_dict = {}
    return(seq_dict)

##functions to recover statistics
#the table of iupac codes may be used by several functions
iupac_codes = {'A':set(['A']),
                'C':set(['C']),
                'G':set(['G']),
                'T':set(['T']),
                'R':set(['A', 'G']),
                'Y':set(['C', 'T']),
                'S':set(['G', 'C']),
                'W':set(['A', 'T']),
                'K':set(['G', 'T']),
                'M':set(['A', 'C']),
                'B':set(['C', 'G', 'T']),
                'D':set(['A', 'G', 'T']),
                'H':set(['A', 'C', 'T']),
                'V':set(['A', 'C', 'G']),
                'N':set(['A', 'C', 'G', 'T']),
                }
                
iupac_het = set(['R','Y','S','W','K','M'])

#counts total number of nucleotides and loci recovered for sample (including gaps)
def N_nuc(sample, pool, seq_dict):
    len_seqs = [len(seq) for [(samp_key, pool_key, locus_key),seq] in seq_dict.iteritems()
                     if samp_key == sample and pool_key == pool]
    return (len(len_seqs), sum(len_seqs))

#counts number of heterozygous nucleotide positions recovered for sample
def N_het(sample, pool, seq_dict):
    het_counter = 0
    seqs = {locus_key:seq
             for [(samp_key, pool_key, locus_key),seq] in seq_dict.iteritems()
             if samp_key == sample and pool_key == pool} #make a dict only with sequences for this sample
    loci = tuple(set(locus for locus in seqs.iterkeys())) #record loci names
    for locus in loci: #for each locus, compare sequences
        seq = seqs[locus]
        for i in xrange(len(seq)):
            if seq[i] in iupac_het: #count all nucleotides that are different within a single locus
                het_counter += 1
    return(het_counter)

#calculates the proportion of GC across all loci recovered for one sample
#ignores ambiguities and gaps
def GC_content(sample, pool, seq_dict):
    GC_counter = 0
    AT_counter = 0

    seqs = {locus_key:seq
             for [(samp_key, pool_key, locus_key),seq] in seq_dict.iteritems()
             if samp_key == sample and pool_key == pool} #make a dict only with sequences for this sample
    loci = tuple(set(locus for locus in seqs.iterkeys())) #record loci names
    for locus in loci: #for each locus, compare sequences
        seq = seqs[locus]
#        print 'seq'
#        print seq
        for i in xrange(len(seq)): #ignore gaps and ambiguities, count only GC
            if seq[i] in set(['C','G','S']):
                GC_counter += 2
            elif seq[i] in set(['A','T','W']):
                AT_counter += 2
            elif seq[i] in set(['M','K','Y','R']):
                AT_counter += 1
                GC_counter += 1
    GC_counter = float(GC_counter)
    return(GC_counter/(AT_counter + GC_counter))




#counts number of ambiguous nucleotide positions recovered for sample
def N_amb(sample, pool, seq_dict):
    iupac_ambiguous = set(['B', 'D','H','V','N'])
    amb_counter = 0
    seqs = {locus_key:seq
             for [(samp_key, pool_key, locus_key),seq] in seq_dict.iteritems()
             if samp_key == sample and pool_key == pool} #make a dict only with sequences for this sample
    loci = tuple(set(seqs.keys())) #record loci names
    for locus in loci: #for each locus, compare sequences
        seq = seqs[locus]
        for i in xrange(len(seq)):
            if seq[i] in iupac_ambiguous: #count all nucleotide positions with ambiguity
                amb_counter += 1
    return(amb_counter)

#counts number of common loci and nucleotides recovered for two samples (gaps included)
def N_common(sample_pool_1, sample_pool_2, seq_dict):
    #samples and pools have to be provided as lists of tuples. E.g. ('BdM1000', '10')
    locus_count = 0
    nuc_count = 0
    loci = tuple(set(locus_key for (samp_key, pool_key, locus_key) in seq_dict.iterkeys())) #record loci names
    for locus in loci:
        seq1_test = (sample_pool_1[0], sample_pool_1[1], locus) in seq_dict
        seq2_test = (sample_pool_2[0], sample_pool_2[1], locus) in seq_dict
        if seq1_test and seq2_test:
            locus_count += 1
            nuc_count += len(seq_dict[(sample_pool_1[0], sample_pool_1[1], locus)])
    if locus_count == 0:
        raise OverlapError
    else:
        return(locus_count, nuc_count)


#counts number of gaps between two samples (i. e. gaps found on the first sample in relation to the second)
def N_gap(sample_pool_1, sample_pool_2, seq_dict):
    #samples and pools have to be provided as lists of tuples. E.g. ('BdM1000', '10')
    locus_count = 0
    nuc_count = 0
    loci = tuple(set(locus_key for (samp_key, pool_key, locus_key) in seq_dict.iterkeys())) #record loci names
    for locus in loci:
        locus_has_gap = False
        seq1_test = (sample_pool_1[0], sample_pool_1[1], locus) in seq_dict
        seq2_test = (sample_pool_2[0], sample_pool_2[1], locus) in seq_dict
        if seq1_test and seq2_test:
            for i in xrange(len(seq_dict[(sample_pool_1[0], sample_pool_1[1], locus)])):
                #for each nucleotide position, will count if any gap in seq1 and no gap in seq2
                gap_test = seq_dict[(sample_pool_1[0], sample_pool_1[1], locus)][i] == '-' and \
                           seq_dict[(sample_pool_2[0], sample_pool_2[1], locus)][i] != '-' 
                if gap_test:
                    nuc_count += 1
                    locus_has_gap = True
        if locus_has_gap:
            locus_count += 1
    return(locus_count, nuc_count)

#calculates, per nucleotide position, the average number of shared alleles between two samples
#ignores gaps and ambiguous nucleotides
def N_matching_nucleotides(sample_pool_1, sample_pool_2, seq_dict):
    #samples and pools have to be provided as lists of tuples. E.g. ('BdM1000', '10')
    iupac_ambiguous = ['B', 'D','H','V','N']
    nuc_matching = []
    nuc_all = []
    loci = tuple(set(locus_key for (samp_key, pool_key, locus_key) in seq_dict.iterkeys())) #record loci names
    for locus in loci:
        temp_matching = []
        seq1_test = (sample_pool_1[0], sample_pool_1[1], locus) in seq_dict
        seq2_test = (sample_pool_2[0], sample_pool_2[1], locus) in seq_dict
        if seq1_test and seq2_test:
            for i in xrange(len(seq_dict[(sample_pool_1[0], sample_pool_1[1], locus)])):
                #retrieve alleles for this position
                nucs = (seq_dict[(sample_pool_1[0], sample_pool_1[1], locus)][i],
                        seq_dict[(sample_pool_2[0], sample_pool_2[1], locus)][i])
                if '-' in nucs: continue #ignore if any gap
                elif set(iupac_ambiguous) & set(nucs): continue #ignore if any ambiguous code
                elif set(iupac_codes[nucs[0]]) == set(iupac_codes[nucs[1]]): #if sets are the same, 2 nucleotides in common
                    temp_matching.append(2)
                else: #if not the same, count nucleotides in common
                    temp_matching.append(len(set(iupac_codes[nucs[0]]) & set(iupac_codes[nucs[1]])))
        if len(temp_matching):
            nuc_matching.append(float(sum(temp_matching)))
            nuc_all.append(len(temp_matching)) 
    return(sum(nuc_matching)/sum(nuc_all))


##functions to parallelize the code
def recover_stats(i): #i must be in sample_info.index
    #first, record sample stats
    try: #some samples have no locus in the end and have to be excluded
        taxon_dict = loci_dicts[sample_info.loc[i,'taxon']]
        sampleid = sample_info.loc[i,'sample']
        poolid = '{:0>2d}'.format(sample_info.loc[i, 'pool'])
        (N_loci, N_nucs) = N_nuc(sampleid, poolid, taxon_dict)
        N_hets = N_het(sampleid, poolid, taxon_dict)
        N_ambs = N_amb(sampleid, poolid, taxon_dict)
        GC = GC_content(sampleid, poolid, taxon_dict)
        return((i, N_loci, N_nucs, N_hets, N_ambs, GC))
    except ZeroDivisionError:
        return(i, np.nan, np.nan, np.nan, np.nan, np.nan)


def recover_pairwise_stats(i):
    pairwise_comp=dict()
    try:
        taxon_dict = loci_dicts[sample_info.loc[i,'taxon']]
    except KeyError:
        sample_info.loc[i,'taxon'] +  ' not in loci file, skipping (' + sample_info.loc[i,'sample'] + '{:0>2d})'.format(sample_info.loc[i, 'pool'])
        return pairwise_comp
    else:
        for j in sample_info.index:
            if sample_info.loc[i,'taxon'] == sample_info.loc[j,'taxon']: #only compute if same taxon
                sample_pool_1 = (sample_info.loc[i,'sample'], '{:0>2d}'.format(sample_info.loc[i, 'pool']))
                sample_pool_2 = (sample_info.loc[j,'sample'], '{:0>2d}'.format(sample_info.loc[j, 'pool']))

                try: #will fail if some of the samples not in loci file
                    (loci, nuc) = N_common(sample_pool_1, sample_pool_2, taxon_dict)
                except OverlapError:
                    pass
                else: #if the first function works, rest should work too
                    pairwise_comp[(sample_pool_1, sample_pool_2, 'N_loci_total')] = loci
                    pairwise_comp[(sample_pool_1, sample_pool_2, 'N_nuc_total')] = nuc

                    (loci, nuc) = N_gap(sample_pool_1, sample_pool_2, taxon_dict)
                    pairwise_comp[(sample_pool_1, sample_pool_2, 'N_loci_withgaps')] = loci
                    pairwise_comp[(sample_pool_1, sample_pool_2, 'N_nuc_gaps')] = nuc

                    pairwise_comp[(sample_pool_1, sample_pool_2, 'N_matching_nuc')] = N_matching_nucleotides(sample_pool_1, sample_pool_2, taxon_dict)

        return(pairwise_comp)

########Now, obtain stats for samples
if __name__ == "__main__":    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--info', default='./sample_info_final.csv', help = 'csv file with sample information')
    parser.add_argument('-t', '--threads', default=os.environ['SLURM_NTASKS'], type = int, help = 'number of threads for parallelization')
    args = parser.parse_args()

    
    print 'PROGRAM STARTING'
    ##First, read info_table
    sample_info = pandas.read_csv(args.info)
    
    ##Then, create new columns
    sample_info.loc[:,'N_loci_s7'] = np.nan
    sample_info.loc[:,'N_nucleotides_s7'] = np.nan
    sample_info.loc[:,'N_nucleotides_heteroz'] = np.nan
    sample_info.loc[:,'N_nucleotides_ambiguous'] = np.nan
    sample_info.loc[:,'GC_content'] = np.nan
    
    
    ##read loci files and save in a dictionary
    loci_dicts = {}
    print 'READING LOCI FILES'
    for taxon in sample_info.taxon.unique():
        try: #not all genera will have loci files
            loci_dicts[taxon] = loci_file_parser(taxon + '.loci')
            print 'loci file read for ' + taxon
        except:
            print 
            print taxon + " doesn't have outfiles"
    
    
    
    ##now, recover information for each sample
    results = Parallel(n_jobs=args.threads, verbose=50)(delayed(recover_stats)(i) for i in sample_info.index)
    print 'FINISHED RECOVERING SAMPLE STATS'
    
    for stats in results:
        sample_info.loc[stats[0],'N_loci_s7'] = stats[1]
        sample_info.loc[stats[0],'N_nucleotides_s7'] = stats[2]
        sample_info.loc[stats[0],'N_nucleotides_heteroz'] = stats[3]
        sample_info.loc[stats[0],'N_nucleotides_ambiguous'] = stats[4]
        sample_info.loc[stats[0],'GC_content'] = stats[5]
    
    print 'LOADED SAMPLE STATS TO TABLE'
    
    
    
    ##now, recover pairwise stats
    results = Parallel(n_jobs=args.threads, verbose=50)(delayed(recover_pairwise_stats)(i) for i in sample_info.index)
    print 'RECOVERED PAIRWISE STATS'
    
    pairwise_comp=dict()
    for i in results:
        #print str(len(i))
        pairwise_comp.update(i)
    print 'JOINED PAIRWISE STATS'
    
    
    ##finally, save info in tables and export
    
    sample_info.to_csv(os.path.basename(args.info).split('.')[0] + '_WGA.csv', index = False)
    
    stats_to_record = set([i for (a, b, i) in pairwise_comp.iterkeys()])
    
    for stat in stats_to_record:
        sample_pools = sorted(set([a for (a, b, c) in pairwise_comp.iterkeys() if c == stat]))
        sample_pools_str = [str(i[0]) + 'pool' + str(i[1]) for i in sample_pools]
        #print sample_pools_str
        temp_dataframe = pandas.DataFrame(data = np.nan, index = sample_pools_str, columns = sample_pools_str)
        for i in xrange(len(sample_pools)):
            for j in xrange(len(sample_pools)):
                try:
                    temp_dataframe.loc[sample_pools_str[i],sample_pools_str[j]] = pairwise_comp[(sample_pools[i], sample_pools[j], stat)]
                except KeyError:
                    pass
    
        temp_dataframe.to_csv('pairwise_{}.csv'.format(stat))
    
    print 'DONE'
