######
#Qanbari Creeping Windows 
######

#BAH 1-Apr-16

creeper<-function(fst.frame, window.size = 1000 ){
	#Creeping Window Analyses based on Qanbari et al. (2012; PLOS ONE) 
		#http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0049525
		#This script takes in a data frame of Fst values for each position in A SINGLE chromosome/scaffold (henceforth "region")
		#The script divides the region into equal bp windows and places SNPs into their appropriate windows
		#For each window, it then estimates mean FST, skipping any window in which the distance between snps is > window.size
		#Outputs 2 plots to help trim windows based on the number of SNPs included in them 
			#NumSNPsbyBin.png - plots the variation in FST per bin (bin is a bin of window lengths)
			#SDbyBin.png - plots the average variation in FST per bin 
			
	#This script requires from the user 2 inputs:
		#1) Data must be in data frame format with 3 columns corresponding to region, position, and fst value. Order data frame on ascending on position. Do not include missing data (i.e. "NA")
			#	       chr1 156 0.00070403
			#	       chr1 203 0.00008735
			#					...
			#	       chr1 22552003 0.00007730
		#2) A window size (variable "window.size") in bp that will be the minimum size of a region in the scan 

	#after running creeper function, a data frame is output. 
		#Col 1 = mean Fst per window
		#Col 2 = number of SNPs per window
		#Col 3 = first SNP in the window
		#Col 4 = last SNP in the window (NOTE: this is a creeping window, so all SNPs will be used more than once)
		
		
		#uses base pair position and it's corresponding FST value to estimate mean fst in 	
		#creeping windows of basepair length (window.size arg). 
	
	names(fst.frame) = c("V2", "POS", "fst45")
	#Gather SNPs and FST as numeric variables
	snp = as.numeric(fst.frame$POS)
	fst = as.numeric(fst.frame$fst45)
	#estimate difference between SNPs	
	p = snp
	p = p[-1]
	p = c(p, NA)
	diff = abs(snp-p)
	diff = diff[!is.na(diff)]

	x=rep(0, (max(snp)+ window.size))
	x[snp]=1
	x=cumsum(x)
	endsnp=x[snp+ window.size]
	n=length(fst)
	vec=rep(0,n)
	len=rep(0,n)
	end=rep(0,n)

	for (i in 1:(n-1)){
			if (diff[i] < window.size ){
				snp_len=length(fst[i:endsnp[i]]) #number of SNPs per window
				fst_mean=mean(fst[i:endsnp[i]]) #mean FST per window
				vec[i]=fst_mean
				len[i]=snp_len
				end[i]=(snp[endsnp[i]])
				
			}else{
				vec[i]=NA
				len[i]=NA
				end[i]=(snp[endsnp[i]])
			}
	}
	creeper = data.frame(cbind(vec,len, snp, end))
	names(creeper)=c("fst","num_snps", "start", "end")
	creeper = creeper[which(complete.cases(creeper$num_snps)),] #remove windows with no SNPs 
	creeper = creeper[which(creeper$num_snps>0),]
	return(creeper)
}
