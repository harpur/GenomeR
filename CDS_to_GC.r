###
# GC and CpGOE for genes
###






#Libraries --------------------------------------------------------------------
library(stringr)
library(plyr)
source("ReferencetoDataframe.r")




#data --------------------------------------------------------------------
bimp <- Reftodf(data.frame="Amel_HAv3.1_genomic.fna", nms="N")
gff <- read.table("AMEL_cds_t_trim.gtf")

#mung -------------------------------------------------------------------------
chrom <- bimp$chrom
chrom <- gsub('Apismellif.*','',chrom) #remove everything after the name of the chromosome so gff$V1 matches
bimp$chrom <- chrom
bimp <- bimp[which(bimp$chrom %in% gff$V1),]


#Estimate GC --------------------------
GC.df <- c()
for(i in 1:nrow(bimp)){
	temp <- gff[which(gff$V1 == as.character(bimp$chrom[i])),]
	seqs <- unlist(strsplit(as.character(bimp$seqs[i]),split=""))
	for(k in 1:nrow(temp)){
		seqs.temp <- seqs[c(temp$V4[k]:temp$V5[k] )]
		cs <- length(seqs.temp[which(seqs.temp=="C")])
		gs <- length(seqs.temp[which(seqs.temp=="G")])
		lens <- length(seqs.temp)
		temp.gc <- cbind(as.character(temp$V10[k]), cs, gs, lens)
		GC.df <- rbind(GC.df, temp.gc)
	}

	print(i)

}
temp <- GC.df
GC.df <- data.frame(temp)


#number of snps per scaffold
GC <- ddply(GC.df , c("V1"), summarize,
		GC = sum(as.numeric(as.character(cs))) + 
		sum(as.numeric(as.character(gs))), 
		len = sum(as.numeric(as.character(lens)))
		)
GC$GCperc <- GC$GC/GC$len
names(GC)[1] <- 'transcript'
write.table(GC,'/depot/bharpur/data/ref_genomes/AMEL/resources/GC/AMEL_CDS_GC', row.names=F, quote=F)


#GC.df$GC = apply(GC.df, 1, function(x) (as.numeric(as.character(x[2])) + 
#	as.numeric(as.character(x[3])))/ as.numeric(as.character(x[4])))

#GC.df = GC.df[which(as.numeric(GC.df$lens) >= 3),]
#names(GC.df)[1] = "Group"



#Estimate CpGOE for gene body ---------
#CpG
#Obs/Exp CpG = Number of CpG * N / (Number of C * Number of G)   
#where N = length of sequence.


gff.temp <- ddply(gff , c("V1", "V10"), summarize,
		st = min(V4),
		en = max(V5)
		)
bimp <- bimp[which(bimp$chrom %in% gff.temp$V1),]


CpG.df <- c()
for(i in 1:nrow(bimp)){
	temp <- gff.temp[which(gff.temp$V1 == as.character(bimp$chrom[i])),]
	seqs <- unlist(strsplit(as.character(bimp$seqs[i]),split=""))
	for(k in 1:nrow(temp)){
		seqs.temp <- paste(seqs[c(temp$st[k]:temp$en[k] )],collapse="")
		cpg <- str_count(seqs.temp, "CG")
		C <- str_count(seqs.temp, "C")
		G <- str_count(seqs.temp, "G")
		lens <- nchar(seqs.temp)
		temp.gc <- cbind(as.character(temp$V10[k]), cpg, C, G, lens)
		CpG.df <- rbind(CpG.df, temp.gc)
	}
	print(i)

}

CpG.df <- data.frame(CpG.df)


CpG.df$C <- as.numeric(as.character(CpG.df$C))
CpG.df$G <- as.numeric(as.character(CpG.df$G))
CpG.df$cpg <- as.numeric(as.character(CpG.df$cpg))
CpG.df$lens <- as.numeric(as.character(CpG.df$lens))



CpG.df$CpGOE = (CpG.df$cpg/CpG.df$lens)/((CpG.df$C/CpG.df$lens) * (CpG.df$G/CpG.df$lens))
CpG.df$GC = (CpG.df$C + CpG.df$G) / CpG.df$lens


write.table(CpG.df,'/depot/bharpur/data/ref_genomes/AMEL/resources/GC/AMEL_GENE_CPGOE', row.names=F, quote=F)




























#GC vs gamma ----------
GC.df = GC.bombus


#gamma
gm.bmel = read.csv("/data2/bombus/Snipre/snipre_run/bayesian_results_BmelToBimp.csv", header=T)
gm.bter = read.csv("/data2/bombus/Snipre/snipre_run/bayesian_results_BterToBimp.csv", header=T)
test = merge(gm.bter, gm.bmel, by="Group",suffix=c(".bter",".bmel"))
gamm = test[c(1,21,45)]


test = merge(gamm, GC.df, by = "Group")
cor.test(test$BSnIPRE.gamma.bter, test$GC,method="spearman")

#	Spearmans rank correlation rho
#
#data:  test$BSnIPRE.gamma.bter and test$GC
#S = 1.807e+11, p-value = 2.869e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#        rho 
#-0.08162812 

#conservative set
mk.cons = read.table("/data2/bombus/bam_Bter_to_Bimp/MKresults_NEW.bayesianresults_new",sep=",",header=T)

names(mk.cons)[1] = "Group"
test = merge(mk.cons, GC.df, by = "Group")

cor.test(test$BSnIPRE.gamma, test$GC,method="spearman")

cor.test(test$BSnIPRE.gamma, test$GC,method="spearman")
#	Spearman's rank correlation rho
#
#data:  test$BSnIPRE.gamma and test$GC
#S = 5101200000, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#     rho 
#-0.199529 



#negative relationship between GC and gamma



#CpG?


#conservative set
mk.cons = read.table("/data2/bombus/bam_Bter_to_Bimp/MKresults_NEW.bayesianresults_new",sep=",",header=T)

names(mk.cons)[1] = "Group"
test = merge(mk.cons, GC.df, by = "Group")

cor.test(test$BSnIPRE.gamma, test$CpGOE,method="spearman")
#	Spearman's rank correlation rho
#
#data:  test$BSnIPRE.gamma and test$CpGOE
#S = 3691100000, p-value = 6.342e-13
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.1320425 
#
#

#strong positive with CpGOE



#W Q, etc genes? ---------------------
degs = read.table(file="/data2/bombus/DEG/rawDEGs.txt",header=T)
degs = aggregate(degs[c(2:10)], by=list(degs$XP), mean)
names(degs)[1]="Group"
sp = rep("NDEG", nrow(degs))

sp[which(degs$FW<0.05 & degs$QG>0.05 & degs$FQ<0.05 & degs$WG>0.05)] = "F" #F
sp[which(degs$FW>0.05 & degs$QG<0.05 & degs$FQ<0.05 & degs$WG>0.05)] = "Q" #Q
sp[which(degs$FW<0.05 & degs$QG<0.05 & degs$FQ<0.05 & degs$WG<0.05)] = "W" #W
sp[which(degs$FW>0.05 & degs$QG<0.05 & degs$FQ>0.05 & degs$WG<0.05)] = "G"#G
table(sp)
degs$sp = sp

test = merge(GC.df, degs, by = "Group")

boxplot(test$GC~test$sp)
x = aov(test$GC ~ test$sp)
TukeyHSD(x)

#               diff          lwr          upr     p adj
#G-F    -0.033403381 -0.087077420 0.0202706582 0.4351157
#NDEG-F -0.023520232 -0.076999761 0.0299592963 0.7513879
#Q-F    -0.034728880 -0.089368247 0.0199104866 0.4128258
#W-F    -0.009608007 -0.065170966 0.0459549517 0.9898669
#NDEG-G  0.009883149  0.004031509 0.0157347885 0.0000407 #
#Q-G    -0.001325499 -0.013960455 0.0113094564 0.9985352
#W-G     0.023795374  0.007626760 0.0399639881 0.0005744 #
#Q-NDEG -0.011208648 -0.022989973 0.0005726772 0.0711900 #
#W-NDEG  0.013912225 -0.001598468 0.0294229184 0.1031028
#W-Q     0.025120873  0.005990217 0.0442515296 0.0031547 #
#worker genes have high GC....


test = merge(test, mk.cons, by = "Group")
test$resid = residuals(lm(test$GC~test$BSnIPRE.gamma))
boxplot(test$resid~test$sp)


#tested with residuals and the pattern vanishes :(





#CpG
#W Q, etc genes? ---------------------
boxplot(test$CpGOE~test$sp)

x = aov(test$CpGOE ~ test$sp)
#CpGOE no relationship

















