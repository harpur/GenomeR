snipre <- function(mydata, CPU=1){
## Main function call
# Package: SnIPRE. Author: Kirsten Eilertson
# Ref: Eilertson KE, Booth JG, Bustamante CD (2012) SnIPRE: Selection Inference Using a Poisson Random Effects Model. PLoS Comput Biol 8(12): e1002806
# Mod: Small changes for speed, Clement F. Kent, February 2013, plus addition of CPU argument; set CPU to the number of cores you wish to use
# also, automate check for required libraries arm and lme4

# load required libraries or exit/change CPU mode:
if (!require("arm",quietly=TRUE)) {print("arm package must be installed. Exiting SnIPRE"); return(NA)}
if (!require("lme4",quietly=TRUE)) {print("lme4 package must be installed. Exiting SnIPRE"); return(NA)}
if (CPU>1) { if (!require("multicore",quietly=TRUE)) {CPU<-1}   }

 data <- mydata
 PS <- data$PS
 PR <- data$PR
 FR <- data$FR
 FS <- data$FS

 n <- length(FS)
 TS <- data$Tsil
 TR <- data$Trepl
 nout = data$nout
 npop = data$npop

 Ivec <- matrix(1, nrow = n)  # makes one vector of appropriate size subset
 d.mu <-as.numeric( matrix(c(PS,PR,FS,FR), ncol =1))
 d.replacement <- as.numeric(c(0,1,0,1)%x%Ivec)
 d.fixed <- as.numeric(c(0,0,1,1)%x%Ivec)
 d.gene <- as.vector(rep(1,4)%x%c(1:n))
 d.TS <- as.vector(rep(1,4)%x%c(TS))
 d.TR <- as.vector(rep(1,4)%x%c(TR))

 count <- as.vector(d.mu)
 R <- as.vector(d.replacement)
 F <- as.vector(d.fixed)
 G <- as.vector(d.gene)
 RF <- as.vector(R*F)
 TR <- as.vector(d.TR)
 TS <- as.vector(d.TS)

 modGEN <- glmer(count~ 1+ R + F + RF +(1+R+F+RF|G),offset =  log(TS*(1-R)+TR*R), family = poisson)
                                        # sel effect
 se.RFG = se.ranef(modGEN)$G[,4]
 re.RFG = ranef(modGEN)$G$RF
 mydata$SnIPRE.est <- fixef(modGEN)[4]+re.RFG     # replacement effect
 clim <- 1.96*se.RFG
 lbound <- mydata$SnIPRE.est - clim
 ubound <- mydata$SnIPRE.est + clim

 beta <- fixef(modGEN)[1]
 betaG <- ranef(modGEN)$G[,1]
 betaF <- fixef(modGEN)[3]
 betaFG <- ranef(modGEN)$G[,3]

 mydata$SnIPRE.lbound <- lbound
 mydata$SnIPRE.ubound <- ubound

 negC <- which(ubound<=0)
 posC <- which(lbound>=0)

 mydata$SnIPRE.class <- "neut"
 mydata$SnIPRE.class[negC] <- "neg"
 mydata$SnIPRE.class[posC] <- "pos"

re.RG = ranef(modGEN)$G$R
se.RG = se.ranef(modGEN)$G[,2]
mydata$SnIPRE.Rest <- fixef(modGEN)[2]+ranef(modGEN)$G[,2]
m<-fixef(modGEN)[2]+re.RG; clim<-1.96*se.RG
Rlbound <- m-clim
Rubound <- m+clim

mydata$SnIPRE.Rlbound <- Rlbound
mydata$SnIPRE.Rubound <- Rubound

negR <- which(Rubound<=0)
posR <- which(Rlbound>=0)
params <- Snipre_tau.theta(beta, betaG,betaF,betaFG, npop, nout)
mydata$SnIPRE.tau = params$tau.est
mydata$SnIPRE.theta = params$theta.est
mydata$SnIPRE.Rclass <- "neut"
mydata$SnIPRE.Rclass[negR] <- "neg"
mydata$SnIPRE.Rclass[posR] <- "pos"
b<-mydata$SnIPRE.est > -4.4
mydata$SnIPRE.gamma[!b] <- -Inf
mydata$SnIPRE.gamma[b] <-  Snipre_LLest2gamm(mydata$SnIPRE.est[b], mydata$SnIPRE.tau[b], n = npop, m= nout)
na.set <- which(is.na(mydata$SnIPRE.gamma) == TRUE)
good.set <- which(b)
use <- setdiff(good.set, na.set)
mydata$SnIPRE.f <- NA
mydata$SnIPRE.f.lb <- NA
mydata$SnIPRE.f.ub <- NA
g.ests <- mydata$SnIPRE.gamma
g.zeros <- which(mydata$SnIPRE.class == "neut")
g.ests[g.zeros] <- .000001
mydata$SnIPRE.f[use] <- Snipre_frac.est(mydata$SnIPRE.Rest[use], g.ests[use], npop[use], nout[use])
mydata$SnIPRE.f.lb[use] <-Snipre_frac.est(mydata$SnIPRE.Rlbound[use], g.ests[use], npop[use],nout[use])
mydata$SnIPRE.f.ub[use] <-Snipre_frac.est(mydata$SnIPRE.Rubound[use], g.ests[use], npop[use],nout[use])
mydata$SnIPRE.f.class <- "neut"
mydata$SnIPRE.f.class[which(mydata$SnIPRE.f.ub <1)] <- "neg"
mydata$SnIPRE.f.class[which(mydata$SnIPRE.f.lb >1)] <- "pos"
return(list(new.dataset = mydata, model = modGEN))
}
