GenomeR
=======
This is, for the time, a directory containing scripts I use day-to-day for dealing with genomic data in R. It is updated often, but assume all scripts and functions, are works-in-progress unless otherwise stated

###VCFFunctions.r 
Functions for dealing with VCF files in R

###VarFunct.r
Various useful functions in R

###snipre.r 
Used for running [Bustamante's SnIPRE](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002806) to estimate the selection coefficient on genes in a give genome

script for running SNIPRE to estimate gamma on genes in a genome

###ReferencetoDataframe.r
Takes in a fasta file, loads it as a data frame

###manhattan.r
Manhattan plots for GWAS analysis

###Q_correct.r
FDR correction using [Storey's Q](http://www.genomine.org/papers/directfdr.pdf)

###creeping.r
For [Qanbari et al. (2012)](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0049525) creeping window analysis with Fst values. The file contains a general work flow for identifying significantly high Fst windows using permutation and corrects for FDR using [Storey's Q](http://www.genomine.org/papers/directfdr.pdf).  It does not yet collate significant windows together into larger regions, but I'll integrate "IRanges" from the Bioconductor suite later. It also doesn't yet plot the result.

