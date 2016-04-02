######
# Various, Useful R functions

#BAH



		
perm.test <- function(x, samp.size, y,N){
samp=c()
for(i in 1:N){
	samp=c(samp, mean(sample(x, samp.size)))
	}
p=1-pnorm(mean(y),mean(samp),sd(samp) )
return(p)
}	
	

norm.interval = function(data, variance = var(data), conf.level = 0.95) {
 z = qnorm((1 - conf.level)/2, lower.tail = FALSE)
 xbar = mean(data)
 sdx = sqrt(variance/length(data))
 c(xbar - z * sdx, xbar + z * sdx)
 }





cohens_d <- function(x, y) {
    lx <- length(x)- 1
    ly <- length(y)- 1
    md  <- abs(mean(x) - mean(y))        ## mean difference (numerator)
    csd <- lx * var(x) + ly * var(y)
    csd <- csd/(lx + ly)
    csd <- sqrt(csd)                     ## common sd computation

    cd  <- md/csd                        ## cohen's d
}




MultiWayOverlapper = function(win.start,win.end,gene.start,gene.end,gene.list) {
  #this is a monster, but basically, looks within each row for genes  overlapping with whatever you want
  blah=outer(as.numeric(unlist(win.start)), as.numeric(unlist(gene.end)), "<=") 
  blah1=outer(as.numeric(unlist(win.end)), as.numeric(unlist(gene.start)), ">=") 
  blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
  if(!is.null(nrow(blah))){return(as.character(unlist(gene.list)[blah[,2]]))}	
  }

lenunique<-function(x){length(unique(x))}


write.list <- function(x,file,...){write.table(x,file,col.names=F,row.names=F,quote=F,...)}





Read.Fasta.DF <- function(data.frame){
#setwd("/media/data1/forty3/brock/AM45Regions")
r=readLines(data.frame);r=as.vector(r)
x=grep(">", r) #so, regions are from 2:19750
ma=x+1;mi=x-1;mi=mi[-1];mi=c(mi, length(r))
chrom=r[x];chrom=gsub(">", "", chrom)
chrom=gsub(" ","",chrom)
	seqs=c(0,1)
	for (i in 1:length(ma)){x=(r[ma[i]:mi[i]]);x=paste(x, collapse = '');seqs=c(seqs,x)}
	seqs=seqs[-c(1,2)]
	AMEL=data.frame(cbind(chrom, seqs))
	return(AMEL)
}


Write.Fasta.DF <- function(data.frame){
	#must have same format as Read.Fasta.DF output
		#col 1 = ID, col 2 = sequence
	seqs=as.character(unlist(data.frame[2]))
	ids=paste(">",as.character(unlist((data.frame[1]))),sep="")
	for (i in 1:nrow(data.frame)){
		sink(file="out.fa", append=T)
		cat(as.character(ids[i]));cat("\n")
		cat(seqs[i]);cat("\n")	
		sink()
	}
}



GC.content <- function(x){
#taxes in a string (fasta format), outputs GC
length(unlist(strsplit(x,""))[unlist(strsplit(x,"")) %in% c("G", "C")])/length(unlist(strsplit(x,"")))
}






#False Discovery Rate (B-H correction)
FDRate=function (pvec, rate) {
    pord <- sort(pvec)
    np <- length(pord)
    maxp <- tail(pord, 1)
    medp <- pord[round(np/2)]
    D <- (np/2)/(maxp - medp)
    ivec <- 1:np
    b <- rate < D * pord/ivec
    nb <- sum(b)
    if (nb == 0) {
        p <- maxp
    }
    else {
        if (nb == np) {
            p <- 0
        }
        else {
            p <- pord[match(TRUE, b) - 1]
        }
    }
    
}

FDRcontrol<-function(pvec,FDR=0.05){
# use the Simes/Benjamini/Hochberg method to choose which of a list of p values satifies a given FDR
s<-sort(pvec, index.return=TRUE); p<-s$x; n<-length(p); ivec<-1:n
b<-p<=FDR*ivec/n; i<-max(c(0,ivec[b]))
b<-ivec<=i
b<-b[(sort(s$ix, index.return=TRUE))$ix]
b
}





#Reverse complement a fasta file
rev.comp<-function(x,rev=TRUE)
{
x<-toupper(x)
y<-rep("N",nchar(x))
xx<-unlist(strsplit(x,NULL))
for (bbb in 1:nchar(x))
    {
        if(xx[bbb]=="A") y[bbb]<-"T"    
        if(xx[bbb]=="C") y[bbb]<-"G"    
        if(xx[bbb]=="G") y[bbb]<-"C"    
        if(xx[bbb]=="T") y[bbb]<-"A"
    }
if(rev==FALSE) 
    {
    for(ccc in (1:nchar(x)))
        {
        if(ccc==1) yy<-y[ccc] else yy<-paste(yy,y[ccc],sep="")
        }
    }
if(rev==T)
    {
    zz<-rep(NA,nchar(x))
    for(ccc in (1:nchar(x)))
        {
        zz[ccc]<-y[nchar(x)+1-ccc]
        if(ccc==1) yy<-zz[ccc] else yy<-paste(yy,zz[ccc],sep="")
        }
    }
    return(as.character(yy))  
} 
 
#Fabio Marroni's rev comp function (thanks) 






######
#Takes a vector of p-values and runs Storey's correction on it
######
# Data stored in in /home/brock/Fst_Scan/creepin:
#
#
#


#############################################
qvalue0<-function (p = NULL, lambda = seq(0, 0.9, 0.05), pi0.method = "smoother", 
    fdr.level = NULL, robust = FALSE, gui = FALSE, smooth.df = 3, 
    smooth.log.pi0 = FALSE) {
# abridgement of Dabney and Storey's qvalue to remove printed errors for non-gui uses
    if (is.null(p)) {
        qvalue.gui()
        return("Launching point-and-click...")
    }
    if (gui & !interactive()) 
        gui = FALSE
    if (min(p) < 0 || max(p) > 1) {
        if (gui) 
            eval(expression(postMsg(paste("ERROR: p-values not in valid range.", 
                "\n"))), parent.frame())
        else print("ERROR: p-values not in valid range.")
        return(0)
    }
    if (length(lambda) > 1 && length(lambda) < 4) {
        if (gui) 
            eval(expression(postMsg(paste("ERROR: If length of lambda greater than 1, you need at least 4 values.", 
                "\n"))), parent.frame())
        else print("ERROR: If length of lambda greater than 1, you need at least 4 values.")
        return(0)
    }
	
	
    if (length(lambda) > 1 && (min(lambda) < 0 || max(lambda) >= 
        1)) {
        if (gui) 
            eval(expression(postMsg(paste("ERROR: Lambda must be within [0, 1).", 
                "\n"))), parent.frame())
        else print("ERROR: Lambda must be within [0, 1).")
        return(0)
    }
    m <- length(p)
    if (length(lambda) == 1) {
        if (lambda < 0 || lambda >= 1) {
            if (gui) 
                eval(expression(postMsg(paste("ERROR: Lambda must be within [0, 1).", 
                  "\n"))), parent.frame())
            else print("ERROR: Lambda must be within [0, 1).")
            return(0)
        }
        pi0 <- mean(p >= lambda)/(1 - lambda)
        pi0 <- min(pi0, 1)
    }
    else {
        pi0 <- rep(0, length(lambda))
        for (i in 1:length(lambda)) {
            pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
        }
        if (pi0.method == "smoother") {
            if (smooth.log.pi0) 
                pi0 <- log(pi0)
            spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
            pi0 <- predict(spi0, x = max(lambda))$y
            if (smooth.log.pi0) 
                pi0 <- exp(pi0)
            pi0 <- min(pi0, 1)
        }
        else if (pi0.method == "bootstrap") {
            minpi0 <- min(pi0)
            mse <- rep(0, length(lambda))
            pi0.boot <- rep(0, length(lambda))
            for (i in 1:100) {
                p.boot <- sample(p, size = m, replace = TRUE)
                for (i in 1:length(lambda)) {
                  pi0.boot[i] <- mean(p.boot > lambda[i])/(1 - 
                    lambda[i])
                }
                mse <- mse + (pi0.boot - minpi0)^2
            }
            pi0 <- min(pi0[mse == min(mse)])
            pi0 <- min(pi0, 1)
        }
		
		
		
        else {
            print("ERROR: 'pi0.method' must be one of 'smoother' or 'bootstrap'.")
            return(0)
        }
    }
    if (pi0 <= 0) {
        if (gui) 
            eval(expression(postMsg(paste("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method.", 
                "\n"))), parent.frame())
        #else print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method.")
        return(0)
    }
    if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 
        1)) {
        if (gui) 
            eval(expression(postMsg(paste("ERROR: 'fdr.level' must be within (0, 1].", 
                "\n"))), parent.frame())
        #else print("ERROR: 'fdr.level' must be within (0, 1].")
        return(0)
    }
    u <- order(p)
    qvalue.rank <- function(x) {
        idx <- sort.list(x)
        fc <- factor(x)
        nl <- length(levels(fc))
        bin <- as.integer(fc)
        tbl <- tabulate(bin)
        cs <- cumsum(tbl)
        tbl <- rep(cs, tbl)
        tbl[idx] <- tbl
        return(tbl)
    }
    v <- qvalue.rank(p)
    qvalue <- pi0 * m * p/v
    if (robust) {
        qvalue <- pi0 * m * p/(v * (1 - (1 - p)^m))
    }
    qvalue[u[m]] <- min(qvalue[u[m]], 1)
    for (i in (m - 1):1) {
        qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 1)
    }
	
	
	
    if (!is.null(fdr.level)) {
        retval <- list(call = match.call(), pi0 = pi0, qvalues = qvalue, 
            pvalues = p, fdr.level = fdr.level, significant = (qvalue <= 
                fdr.level), lambda = lambda)
    }
    else {
        retval <- list(call = match.call(), pi0 = pi0, qvalues = qvalue, 
            pvalues = p, lambda = lambda)
    }
    class(retval) <- "qvalue"
    return(retval)
}
#############################################
qvalue1<-function(p){
# call Storey qvalue with defaults; if fails call with bootstrap; if fails call with lambda=0
x<-try(qvalue0(p),silent=TRUE)
if("qvalue"==class(x)) {q<-x$qvalues;pi0<-x$pi0}
else {
  x<-try(qvalue0(p,pi0.method="bootstrap"),silent=TRUE)
  if("qvalue"==class(x)) {q<-x$qvalues;pi0<-x$pi0}
  else {x<-qvalue0(p,lambda=0); q<-x$qvalues; pi0<-x$pi0}
  }
list(q=q,pi0=pi0)
}
#############################################








