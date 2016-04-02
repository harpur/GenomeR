###
# ReftoDf
###

#Take the AMEL45new.fasta and outputs a df with 2 colums (chrom, seqs) scaffold (1c1g) and sequence vector 
	
	
	
	
	
	
	
	
	
	
	
	
	
	
Reftodf<-function(data.frame="/home/amel45/AM45/am45new.fasta",nms="Y"){	
#setwd("/media/data1/forty3/brock/AM45Regions")
r=readLines(data.frame);r=as.vector(r)
x=grep(">", r) #so, regions are from 2:19750
ma=x+1;mi=x-1;mi=mi[-1];mi=c(mi, length(r))
chrom=r[x];chrom=gsub(">", "", chrom)
chrom=gsub(" ","",chrom)
if(nms=="Y"){
	chrom=gsub("[.]", "c", chrom);g=rep("g", length(x));chrom=paste(chrom,g, sep="")
	seqs=c(0,1)
	for (i in 1:length(ma)){x=(r[ma[i]:mi[i]]);x=paste(x, collapse = '');seqs=c(seqs,x)}
	seqs=seqs[-c(1,2)]
	AMEL=data.frame(cbind(chrom, seqs))
	return(AMEL)
}else{
	seqs=c(0,1)
	for (i in 1:length(ma)){x=(r[ma[i]:mi[i]]);x=paste(x, collapse = '');seqs=c(seqs,x)}
	seqs=seqs[-c(1,2)]
	AMEL=data.frame(cbind(chrom, seqs))
	return(AMEL)
}
}



CDStodf<-function(data.frame="/home/amel45/AM45/amel_OGSv3.2_cds.fa"){	
#setwd("/media/data1/forty3/brock/AM45Regions")
r=readLines(data.frame);r=as.vector(r)
x=grep(">", r) #so, regions are from 2:19750
ma=x+1;mi=x-1;mi=mi[-1];mi=c(mi, length(r))
chrom=r[x];
chrom=gsub("gnl[|]Amel_4.5[|]", "", chrom);
chrom=gsub("-RA", "", chrom)
chrom=gsub(">", "", chrom)
seqs=c(0,1)
for (i in 1:length(ma)){x=(r[ma[i]:mi[i]]);x=paste(x, collapse = '');seqs=c(seqs,x)}
seqs=seqs[-c(1,2)]
AMEL=data.frame(cbind(chrom, seqs))
return(AMEL)
}





