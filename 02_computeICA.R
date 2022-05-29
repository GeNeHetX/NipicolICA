source("/Users/remy.nicolle/Workspace/IMSI/SRC_load.R")
library(ggridges)
library(ggpubr)
library(pROC)
library(officer)
library(magrittr)
library(rvg)
library(GSVA)
library(qpathway) # devtools::install_github("GeNeHetX/qpathway")
library(stringr)
library(gtsummary)

newPATHS=qpathway::loadPath()




# NIPICOLEXPT EXPL ALLCENTROIDS
# load("/Users/remy.nicolle/Workspace/IMSI/data/historicalDatasets/tmsick.cittcga_20160407.RData")
# cittcgaa

#do/get ICA

k= 10 ; ng=5000 ; anorm="sc" ; normprojtype="sc"

XX=normalize(NIPICOLEXPT[which(rank(1/rowMeanNSd(NIPICOLEXPT))<=ng),],anorm)
SELICA=doica(XX,k=k,docor=T,maxiter = 1e+08, eps = 10^-6)
SELICA$corg=getUniqueGeneMat(SELICA$cor,ensannot[rownames(SELICA$cor),"GeneName"],rowMaxs(abs(SELICA$cor)))
SELICA$Sg=getUniqueGeneMat(SELICA$S,ensannot[rownames(SELICA$S),"GeneName"],rowMaxs(abs(SELICA$S)))

maxsign=apply(colRanges(SELICA$S),1,function(xr){sign(xr)[which.max(abs(xr))] })
SELICA$dirS=t(t(SELICA$S)*maxsign)
SELICA$dirSg=t(t(SELICA$Sg)*maxsign)
SELICA$dirA=t(t(SELICA$A)*maxsign)
saveRDS(SELICA,file="computedICA.rds")


SELICA=readRDS(file="SELICA.rds")



PROJL=lapply(EXPL,function(anexp){
  comensg=intersect(rownames(SELICA$dirS),rownames(anexp))
  comsymg=intersect(rownames(SELICA$dirSg),rownames(anexp))

  invs=NULL;comg=NULL
  if(length(comensg) >1000){
    comg=comensg;invs = MASS::ginv(SELICA$dirS[comg,])
  }else if(length(comsymg) >1000){
    comg=comsymg;invs = MASS::ginv(SELICA$dirSg[comg,])
  }else{
    return(NULL)
  }
  prj=t(invs%*%normalize(anexp[comg, ],normprojtype))
  colnames(prj)=colnames(SELICA$S)
  prj
})
