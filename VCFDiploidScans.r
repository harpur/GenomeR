###
# Pi Across the Whole Genome in size W windows
###







  
#Functions
Read.VCF=function(file=i){
	#To load a VCF file set i="FILENAME"
	#then type "vcf=Read.VCF()
	x=readLines(i);has=grep(pattern="##", x)
	x=x[-has]
	x=strsplit(x,split="\t")
	return(x)	
	}
		


MinAllele=function(x){
	miss=length(grep("[.]/[.]",unlist(x)))
	if(miss>1){
		return("NA")}
	else{
		zero=2*length(grep("0/0",unlist(x)))+length(grep("0/1",unlist(x)))
		one=2*length(grep("1/1",unlist(x)))+length(grep("0/1",unlist(x)))
		if(zero>one){return(one)}else{return(zero)}
		}
	}

	
	
	
WhichMinAllele=function(x){
	miss=length(grep("[.]/[.]",unlist(x)))
	if(miss>1){
		return("NA")}
	else{
		zero=2*length(grep("0/0",unlist(x)))+length(grep("0/1",unlist(x)))
		one=2*length(grep("1/1",unlist(x)))+length(grep("0/1",unlist(x)))
		if(zero>one){return("1")}else{return("0")}
		}
	}	
	
	

		
VCF.Pi=function(file=i){
	x=Read.VCF(file=i)
	st=grep("FORMAT", unlist(x[1]))+1
	x=x[-1]
	en=length(unlist(x[1]))
	#open df and get the number of samples (after FORMAT)
	chrom=sapply(x,function(x) return(c(x[1],x[2])))
	ones=sapply(x,function(x) (length(x[st:en][grep("^1/1:", x[st:en])]))) #counts number of "1:'s"
	twos=sapply(x,function(x) (length(x[st:en][grep("^0/0:", x[st:en])]))) #counts number of "1:'s"
	hets=sapply(x,function(x) (length(x[st:en][grep("^0/1:", x[st:en])])))
	miss=sapply(x,function(x) (length(x[st:en][grep("[.]/[.]:", x[st:en])]))) 
	miss=(abs(twos+ones+hets));miss=abs(miss-((en-st)+1))
	miss=((en-st)+1)-miss #now, gets sample size
	numInds=(en-st)+1
	#Function assumes haploid, counts 0's, 1's and missing
	#1=ref, 0=alt, 
	#Now, need which SNP is least common 
	blah=as.character(ones>=twos)
	#If TRUE, 2*ones, if FALSE, 2*twos
	blah[which(blah=="TRUE")]	= (2 * twos[which(blah=="TRUE")] + hets[which(blah=="TRUE")])
	blah[which(blah=="FALSE")] =	(2 * ones[which(blah=="FALSE")] + hets[which(blah=="FALSE")])
	#blah that have miss==0 need to be removed.
	blah[miss=="0" | miss=="1"]=0 	
	blah=as.numeric(blah)
	blah=blah*(miss-(blah/2))
	blah=round(blah/(miss*(miss-1)),4)
	blah=cbind(chrom[1,],chrom[2,],miss,blah)
	return(blah)
	}	
	
	
	
Window.Theta=function(S=S, a1=a1, w=w){
			k = sum(S)/a1 ;k = k/w
			return(k)
	}
	
VCF.TDThetaPiWindows=function(file=i, w=1000){
	x=Read.VCF(file=i)
	st=grep("FORMAT", unlist(x[1]))+1
	x=x[-1]
	en=length(unlist(x[1]))
	#open df and get the number of samples (after FORMAT)
	chrom=sapply(x,function(x) return(c(x[1],x[2])))
	numInds=(en-st)+1

	##Pre-amble 1) Get Windows:
	#Logic from BAH's "Creepinv2.1.r" for making nucleotide windows
	#Strike that as of 8-Aug-14 I've changed this. Creepin jumps and makes windows from each SNP starting at 1
	#so, second window isn't 5001-10000 but that next SNP to 5K after
	#instead, I'm using ceiling function. Clever, I know.
	snps=as.numeric(sapply(x, function(x) return(x[2])))
	endsnp=ceiling(snps/w) #0-w=1; w-2w=2, etc.
	#So, I can use this logic to get out windows of w ncl's in length
	
	##Pre-amble 2) Get TD coefficients:
	#TD-this logic from Wright Lab's "slidingWindowStats.py" and Tajima 1989
	N=numInds
	a1=0
	a2=0
	for(u in c(1:(N-1))){
		a1 = c(a1,(1/u))
		a2 = c(a2, 1/(u*u))
		}
	a1=sum(a1);a2=sum(a2)
	b1 = (N+1)/(3*(N-1))#see Tajima 1989
    b2 = (2*(N*N+N+3))/(9*N*(N-1))
    c1 = b1 - 1/a1
    c2 = b2-(N+2)/(a1*N)+a2/(a1*a1)
    e1 = c1/a1
    e2 = c2/(a1*a1+a2)
	
	##Pre-amble 3) Theta coefficients:
	#Geting Theta--Follow's AD's logic 
	S=c()
	nS=c()
	for(site in 1:length(x)){
		#site=1
		hets = length(x[[site]][st:en][grep("^0/1:", x[[site]][st:en])])
		if(hets >= 1){
			S=c(S,1)
			nS=c(nS,0)
		}else{
			ones = length(x[[site]][st:en][grep("^1/1:", x[[site]][st:en])])
			twos = length(x[[site]][st:en][grep("^0/0:", x[[site]][st:en])])
			miss = sapply(x,function(x) (length(x[st:en][grep("[.]/[.]:", x[st:en])]))) 
			miss = (abs(twos+ones+hets));miss=abs(miss-((en-st)+1))
			miss = (abs(twos+ones+hets));miss=abs(miss-((en-st)+1))
			miss = ((en-st)+1)-miss
			ones = 2*ones+hets;twos=2*twos+hets
		if(miss == ones | miss == twos){ #Now, count segregating sites (S) and non-segregating (nS)
			S=c(S,0)
			nS=c(nS,1)
			}else{
			S=c(S,1)
			nS=c(nS,0)
			}
		}
		}
		#Above, miss becomes the samplesize (after correcting for missing values)
		#If the site is fixed for a polymorphism, miss==ones or miss==twos, so I
		#add on to the nS vector a 1, allowing me to sum the number of sites as sum(nS) and
		#sum(S). Theta can then be calculated as S/(windowSize * a1)
		
		##Theta Windows (using function Window.Theta):
		#I can just aggregate across "endsnp" for every S
		Theta.Win=aggregate(S,by=list(endsnp), Window.Theta,a1,w) #Theta/window
		PolySNPsinWindow=aggregate(S,by=list(endsnp), function(x) length(x[x>0])) #get out number of SNPs
		StaEnd=aggregate(snps,by=list(endsnp), range) #Get out the Start and End SNPs of the window
				
		##PiWindows
		pi=VCF.Pi(file=i);pi=as.numeric(pi[,4])
			#outputs Pi at each SNP (blah, the 4th column), and the number of sites (miss)
		S.Win=aggregate(S,by=list(endsnp), sum)
		Pi.Win=aggregate(pi,by=list(endsnp), function(x) sum(x)/w)

		
		##Tajima's D Windows
		#Tajima's D, Page 260:
		#https://encrypted.google.com/books?id=0nt-qaAfIbAC&printsec=frontcover&source=gbs_ViewAPI#v=onepage&q&f=false
			#d=pi-theta.....easy enough.
		TD=w*(Pi.Win$x) - w*(Theta.Win$x)
		TD[Pi.Win$x==0 | Theta.Win$x==0]=NA
			
		p1=aggregate(S,by=list(endsnp), function(x) sum(x)*e1)
		p2=aggregate(S,by=list(endsnp), function(x) sum(x)*e2*(sum(x)-1))
		p1p2=sqrt(p1$x+p2$x)
		TD=TD/p1p2
		#head(TD)
		output=cbind(start=StaEnd$x[,1],end=StaEnd$x[,2], pSNPs=PolySNPsinWindow$x,Theta=Theta.Win$x, Pi=Pi.Win$x, TD=TD)
		i=gsub(".vcf", ".divers", i)
		write.table(output, file=i, row.names=F, col.names=T, quote=F)
		}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

