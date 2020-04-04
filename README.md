GenomeR
=======
This is, for the time, a directory containing hacky scripts I have used for dealing with genomic data in R. They are in various states of completion. When frequetly used, they git pushed into developing R packages or pipelines.

It also contains some data sets I use for honey bee data. I will add to this occasionally. 

# Data Sets
"Hunt_markers.txt" contains a list of markers from Greg Hunt. These markers are from his QTL papers and the set here contains those. Some of these have been published and some have not. I found this list on Greg's desktop computer when I inherited the lab.


"Solignac_x_Cornuet_2007_BeeLinkMap_Supp.xls" contains the OG honey bee recombination map. 



# VCFFunctions.r 

Functions for dealing with VCF files in R

# VarFunct.r
(outdated)
Various useful functions in R

# snipre.r 
(outdated)
Used for running [Bustamante's SnIPRE](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002806) to estimate the selection coefficient on genes in a give genome

script for running SNIPRE to estimate gamma on genes in a genome

# ReferencetoDataframe.r

Takes in a fasta file, loads it as a data frame

# manhattan.r

Manhattan plots for GWAS analysis

# Q_correct.r

FDR correction using [Storey's Q](http://www.genomine.org/papers/directfdr.pdf)

# creeping.r

(outdated)
For [Qanbari et al. (2012)](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0049525) creeping window analysis with Fst values. After removing windows that don't contain enough SNP sites, 


To come is a general work flow for identifying significantly high Fst windows. using perm.creeper() function and then correcting for multiple tests using [Storey's Q](http://www.genomine.org/papers/directfdr.pdf), found in Q_correct.r.  High Fst windwos can then be collated with IRanges" from the BioconductoR suite

An example of this would be:

<pre><code># Load Required  Packages ------------------------
source("Q_correct.r")
source("VarFunct.r") #for permutation function
	
#run Permutation --------------------------------- 	
	#"creeper" data frame is generated from creeper() function
head(creeper)
        fst num_snps start  end
2 0.1428318        3  3790 4296
3 0.1912891        2  3798 4296
5 0.1949204        3  5579 6429
num <- creeper$num_snps
lec <- perm.creeper(n=100000, num=num)
creeper$pnorm.fst <- pnorm(creeper$fstHvL, mean=mean(lec), sd=sd(lec),lower.tail=F)
creeper$q <- qvalue1(creeper$pnorm.fst)$q
</code></pre>

# ancestryPlots.R 

This function takes output from [Ancestry HMM](http://biorxiv.org/content/early/2016/07/15/064238) for many individuals and plots the data as an average +/- SD ancestry of the focal lineage. 

usage: 
<pre><code>
Rscript ancestryPlots.R <FILE1 .ahmm.maxpost>
</code></pre>




