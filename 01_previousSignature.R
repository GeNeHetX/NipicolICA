dev.new()
source("/Users/remy.nicolle/Workspace/IMSI/SRC_load.R")
library(ggridges)
library(ggpubr)
library(pROC)
library(officer)
library(magrittr)
library(rvg)
library(GSVA)
library(qpathway)
library(stringr)
library(gtsummary)

newPATHS=qpathway::loadPath()

# nipidf  rckidf ideaIntradf ideaFrontdf TCGAsurvdf EXPL




WISPL=setNames(lapply(names(EXPL),function(anexpN){
  print(anexpN)
  expg=getUniqueGeneMat(EXPL[[anexpN]],PROBEL[[anexpN]],rowSds(EXPL[[anexpN]]))
  comg=intersect(rownames(expg),rownames(ALLCENTROIDS$TCGA.HiSeq))
  cmswisped=wisp(scale(expg[comg,],scale=T),scale(ALLCENTROIDS$TCGA.HiSeq[comg,],scale=T))[,c(1:4,16)]
  colnames(cmswisped)=sub("weight.TCGA.HiSeq.centered.","",colnames(cmswisped))
  cmswisped$topWeightedClass=sub("TCGA.HiSeq.centered.","",cmswisped$topWeightedClass)
  cmswisped
}),names(EXPL))


#

# immunoMSIUQg  nipicolUQg
# rickiresp nipiresp

# Direct gene/pop surv. figure 1 of RNA : predetermined rna signature and genes


bioppn=grep("BioPlanet",names(newPATHS),value=T)

selpathwayneg=c(
  "Elsevier_Pathway_Collection__Epithelial_to_Mesenchymal_Transition_in_Cancer:_Overview" ,
  "c2.cp.wikipathways__WP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_IN_COLORECTAL_CANCER",
  "Elsevier_Pathway_Collection__TGFB_Family_in_Epithelial_to_Mesenchymal_Transition_in_Cancer",
  "BioPlanet__TGF-beta_signaling_pathway",
  "Elsevier_Pathway_Collection__WNT_in_Epithelial_to_Mesenchymal_Transition_in_Cancer",
  "BioPlanet__Wnt_signaling_pathway" ,
  "BioPlanet__TGF-beta_receptor_signaling_in_EMT_(epithelial_to_mesenchymal_transition)",
  "BioPlanet__VEGF,_hypoxia,_and_angiogenesis"  ,
  "BioPlanet__Angiogenesis",
  "BioPlanet__mTOR_signaling_pathway",
  "BioPlanet__Tumor_necrosis_factor_(TNF)_pathway" ,
  "BioPlanet__TNF-alpha_effects_on_cytokine_activity,_cell_motility,_and_apoptosis",
  "BioPlanet__TNF-alpha_signaling_pathway",
  "Elsevier_Pathway_Collection__KRAS_Signaling",
  "c2.cp.biocarta__BIOCARTA_IFNA_PATHWAY",
  "c2.cp.biocarta__BIOCARTA_IFNG_PATHWAY"

)



selrickipathv=gsva(immunoMSIUQg, newPATHS[selpathwayneg], method="gsva",kcdf="Gaussian",min.sz=2,max.sz=20000)
selnipipathv=gsva(nipicolUQg, newPATHS[selpathwayneg], method="gsva",kcdf="Gaussian",min.sz=2,max.sz=20000)



pathpronodf=do.call(rbind,mclapply(compath,function(g){
  rickicox=summary(coxph(Surv(PFS,PD)~rickipathv[g,rownames(rickiresp)],rickiresp))
  nipicox=summary(coxph(Surv(PFS,PD)~nipipathv[g,rownames(nipiresp)],nipiresp))
  data.frame(gene=g,NIPIHR=coef(nipicox)[1,2],NIPIpv=nipicox$sctest[3],
  RICKIHR=coef(rickicox)[1,2],RICKIpv=rickicox$sctest[3],row.names=g)
},mc.cores=8))

pathpronodf[which(rowMaxs(as.matrix(pathpronodf[,c(3,5)]))<0.05),]



h=2;w=3
appt=read_pptx() %>%
add_slide(layout = "Blank", master = "Office Theme") %>%
ph_with(dml({
  print(shortunivforest(

    list(
      TCells= (coxph(Surv(PFS,PD)~rickisignatures[rownames(rickiresp),"mcpcount.Tcells"],rickiresp)),
      BCells= (coxph(Surv(PFS,PD)~rickisignatures[rownames(rickiresp),"mcpcount.B.lineage"],rickiresp)),
      Mono= (coxph(Surv(PFS,PD)~rickisignatures[rownames(rickiresp),"mcpcount.Mono.lineage"],rickiresp)),
      Neutro= (coxph(Surv(PFS,PD)~rickisignatures[rownames(rickiresp),"mcpcount.Neutrophils"],rickiresp)),
      Fibro= (coxph(Surv(PFS,PD)~rickisignatures[rownames(rickiresp),"mcpcount.Fibroblasts"],rickiresp)),
      PDL1=   (coxph(Surv(PFS,PD)~immunoMSIUQg["CD274",rownames(rickiresp)],rickiresp)),
      PD1=   (coxph(Surv(PFS,PD)~immunoMSIUQg["PDCD1",rownames(rickiresp)],rickiresp)),
      CTLA4=  (coxph(Surv(PFS,PD)~immunoMSIUQg["CTLA4",rownames(rickiresp)],rickiresp))
    )

    ,new_page=F,
    mar= unit(rep(1, times = 4), "mm"),
    colgap=unit(1, "mm"),xticks=c(0.5,1,2,4)
  ))
},pointsize=8),ph_location(left=w+0.2,top=0,width=w,height=h))%>%
ph_with(dml({
  print(shortunivforest(

    list(
      TCells= (coxph(Surv(PFS,PD)~signatures[rownames(nipiresp),"mcpcount.Tcells"],nipiresp)),
      BCells= (coxph(Surv(PFS,PD)~signatures[rownames(nipiresp),"mcpcount.B.lineage"],nipiresp)),
      Mono= (coxph(Surv(PFS,PD)~signatures[rownames(nipiresp),"mcpcount.Mono.lineage"],nipiresp)),
      Neutro= (coxph(Surv(PFS,PD)~signatures[rownames(nipiresp),"mcpcount.Neutrophils"],nipiresp)),
      Fibro= (coxph(Surv(PFS,PD)~signatures[rownames(nipiresp),"mcpcount.Fibroblasts"],nipiresp)),
      PDL1=   (coxph(Surv(PFS,PD)~nipicolUQg["CD274",rownames(nipiresp)],nipiresp)),
      PD1=   (coxph(Surv(PFS,PD)~nipicolUQg["PDCD1",rownames(nipiresp)],nipiresp)),
      CTLA4=  (coxph(Surv(PFS,PD)~nipicolUQg["CTLA4",rownames(nipiresp)],nipiresp))
    )

    ,new_page=F,
    mar= unit(rep(1, times = 4), "mm"),
    colgap=unit(1, "mm"),xticks=c(0.5,1,2,4)
  ))
},pointsize=8),ph_location(left=0,top=0,width=w,height=h))%>%
ph_with(dml({
  print(shortunivforest(
    setNames(lapply(paste0("CMS",1:4),\(x){
      y=(10*WISPL$nipicol[rownames(nipiresp),x])
      (coxph(Surv(PFS,PD)~y,nipiresp))
    }),paste0("CMS",1:4))

    ,new_page=F,addN=F,
    mar= unit(rep(1, times = 4), "mm"),
    colgap=unit(1, "mm"),xticks=c(0.2,0.5,1,2,4)
  ))
},pointsize=8),ph_location(left=0,top=w,width=w,height=h))%>%
ph_with(dml({
  print(shortunivforest(
    setNames(lapply(paste0("CMS",c(1,2,4)),\(x){
      y=10*(WISPL$ricki[rownames(rickiresp),x])
      (coxph(Surv(PFS,PD)~y,rickiresp,iter.max=100))
    }),paste0("CMS",c(1,2,4)))

    ,new_page=F,addN=F,
    mar= unit(rep(1, times = 4), "mm"),
    colgap=unit(1, "mm"),xticks=c(0.2,0.5,1,2,4)
  ))
},pointsize=8),ph_location(left=h,top=w,width=w,height=h)) %>%



add_slide(layout = "Blank", master = "Office Theme") %>%
ph_with(dml({
  print(shortunivforest(
    setNames(lapply(selpathwayneg,\(x){
      y=scale(selnipipathv[x,rownames(nipiresp)])
      (coxph(Surv(PFS,PD)~y,nipiresp))
    }),sub("^[^_]+","",selpathwayneg))

    ,new_page=F,addN=F,
    mar= unit(rep(1, times = 4), "mm"),
    colgap=unit(1, "mm"),xticks=c(0.2,0.5,1,2,4)
  ))
},pointsize=8),ph_location(left=0,top=0,width=w+4,height=h*1.5))%>%
ph_with(dml({
  print(shortunivforest(
    setNames(lapply(selpathwayneg,\(x){
      y=scale(selrickipathv[x,rownames(rickiresp)])
      (coxph(Surv(PFS,PD)~y,rickiresp,iter.max=100))
    }),sub("^[^_]+","",selpathwayneg))

    ,new_page=F,addN=F,
    mar= unit(rep(1, times = 4), "mm"),
    colgap=unit(1, "mm"),xticks=c(0.2,0.5,1,2,4)
  ))
},pointsize=8),ph_location(left=0,top=h*2,width=w+4,height=h*1.5))%>%



print(appt  , target = "RNAPREDEF.pptx")

}

# Figure showing all pathways, immune signatures and genes
# GENERAL pop and gene surv test. figure 1 of RNA : all rna signature and genes
{

  comborickiresp=rickiresp[which(rickiresp$DRUGD1=="Nivolumab | Ipilimumab"),]

  rickipathv=gsva(immunoMSIUQg, newPATHS, method="gsva",kcdf="Gaussian",min.sz=20,max.sz=200)
  nipipathv=gsva(nipicolUQg, newPATHS, method="gsva",kcdf="Gaussian",min.sz=20,max.sz=200)

  compath=intersect(rownames(rickipathv),rownames(nipipathv))
  pathpronodf=do.call(rbind,mclapply(compath,function(g){
    rickicox=summary(coxph(Surv(PFS,PD)~rickipathv[g,rownames(rickiresp)],rickiresp))
    ricki2cox=summary(coxph(Surv(PFS,PD)~rickipathv[g,rownames(comborickiresp)],comborickiresp))

    nipicox=summary(coxph(Surv(PFS,PD)~nipipathv[g,rownames(nipiresp)],nipiresp))
    data.frame(gene=g,NIPIHR=coef(nipicox)[1,2],NIPIpv=nipicox$sctest[3],
    RICKIHR=coef(rickicox)[1,2],RICKIpv=rickicox$sctest[3],
    RICKI2HR=coef(ricki2cox)[1,2],RICKI2pv=ricki2cox$sctest[3],row.names=g)
  },mc.cores=8))

  pathpronodf[which(rowMaxs(as.matrix(pathpronodf[,c(3,5)]))<0.05),]

  # a=cor(t(nipipathv[compath,rownames(PROJL$nipicol)]),PROJL$nipicol[,9],method="sp")
  # b=cor(t(rickipathv[compath,rownames(PROJL$ricki)]),PROJL$ricki[,9],method="sp")

  comsig=intersect(colnames(rickisignatures),colnames(signatures))
  comsig=comsig[grep("^xcell|^epic|^cibersort|^mcpcount|^Estimate",comsig)]
  sigpronodf=do.call(rbind,mclapply(comsig,function(g){
    rickicox=summary(coxph(Surv(PFS,PD)~rickisignatures[rownames(rickiresp),g],rickiresp))
    ricki2cox=summary(coxph(Surv(PFS,PD)~rickisignatures[rownames(comborickiresp),g],comborickiresp))

    nipicox=summary(coxph(Surv(PFS,PD)~signatures[rownames(nipiresp),g],nipiresp))
    data.frame(gene=g,NIPIHR=coef(nipicox)[1,2],NIPIpv=nipicox$sctest[3],
    RICKIHR=coef(rickicox)[1,2],RICKIpv=rickicox$sctest[3],
    RICKI2HR=coef(ricki2cox)[1,2],RICKI2pv=ricki2cox$sctest[3],row.names=g)
  },mc.cores=8))


  # openxlsx::write.xlsx(sigpronodf,file="table/TMEsignatureSurv.xlsx")

  a=cor((signatures[rownames(PROJL$nipicol),comsig]),PROJL$nipicol[,9],method="sp")
  b=cor((rickisignatures[rownames(PROJL$ricki),comsig]),PROJL$ricki[,9],method="sp")

  comsig[which(a>0.5 & b>0.5)]
  comsig[which(a < -0.5 & b < -0.5)]


  comg=intersect(rownames(nipicolUQg),rownames(immunoMSIUQg))
  comg=comg[which(rowMeans(nipicolUQg[comg,])>2 & rowMeans(immunoMSIUQg[comg,])>2)]

  gpronodf=do.call(rbind,mclapply(comg,function(g){
    rickicox=summary(coxph(Surv(PFS,PD)~immunoMSIUQg[g,rownames(rickiresp)],rickiresp))
    ricki2cox=summary(coxph(Surv(PFS,PD)~immunoMSIUQg[g,rownames(comborickiresp)],comborickiresp))

    nipicox=summary(coxph(Surv(PFS,PD)~nipicolUQg[g,rownames(nipiresp)],nipiresp))
    data.frame(gene=g,NIPIHR=coef(nipicox)[1,2],NIPIpv=nipicox$sctest[3],
    RICKIHR=coef(rickicox)[1,2],RICKIpv=rickicox$sctest[3],
    RICKI2HR=coef(ricki2cox)[1,2],RICKI2pv=ricki2cox$sctest[3],row.names=g)
  },mc.cores=8))

  gpronodf[which(gpronodf$NIPIpv<0.01 & gpronodf$RICKIpv<0.01),]

  gpronodf[which(p.adjust(gpronodf$RICKIpv,"fdr")<0.05 & gpronodf$NIPIpv<0.05),]

  alpha=-log10(0.05)




  plot2pvdf=function(df,apre="NIPI",bpre="RICKI2",laba="NIPICOL",labb="RICKI",minlogpv=6,...){

    xa=df[,paste0(apre,"pv")]
    xb=df[,paste0(bpre,"pv")]
    a=-log10(xa)*sign( ((df[,paste0(apre,"HR")] >1)*2)-1 )
    b=-log10(xb)*sign( ((df[,paste0(bpre,"HR")] >1)*2)-1 )


    if(max(abs(c(a,b)))>minlogpv){
      print("WARNIGN increased min log pv")
      # minlogpv=ceiling(max(c(-log10(a),-log10(b))))
    }


    pvticks=c(0.05,0.01,0.001,0.0001,10^-5,10^-6,10^-7,10^-10)
    ti=c(rev(log10(pvticks)),-log10(pvticks))
    plot(x=a,y=b,xlim=c(-minlogpv,minlogpv),ylim=c(-minlogpv,minlogpv),pch=16,axes=F,xlab=laba,ylab=labb,
    col=c("black","red")[(1+(p.adjust(xa,"fdr")<0.05|p.adjust(xb,"fdr")<0.05))],...)

    # plot(x=a,y=b,xlim=c(-minlogpv,minlogpv),ylim=c(-minlogpv,minlogpv),pch=16,axes=F,xlab=laba,ylab=labb,
    # col=c("black","red")[(1+(p.adjust(xa,"fdr")<0.05|p.adjust(xb,"fdr")<0.05))])


    axis(1,at=ti,labels=c(rev(pvticks),pvticks),las=2)
    axis(2,at=ti,labels=c(rev(pvticks),pvticks),las=2)

    abline(h=0,col="grey20",lty=3)
    abline(v=0,col="grey20",lty=3)
  }


  plot2pvlog=function(a,b,laba="NIPICOL",labb="RICKI",minlogpv=6,...){
    if(max(c(-log10(a),-log10(b)))>minlogpv){
      print("WARNIGN increased min log pv")
      # minlogpv=ceiling(max(c(-log10(a),-log10(b))))
    }

    pvticks=c(1,0.2,0.05,0.01,0.001,0.0001,0.00001)
    plot(x=-log10(a),y=-log10(b),xlim=c(0,minlogpv),ylim=c(0,minlogpv),pch=16,axes=F,xlab=laba,ylab=labb,
    col=c("black","red")[(1+(p.adjust(a,"fdr")<0.05|p.adjust(b,"fdr")<0.05))],...)
    axis(1,at=-log10(pvticks),labels=pvticks,las=2)
    axis(2,at=-log10(pvticks),labels=pvticks,las=2)
  }

  subpathpronodf=pathpronodf[which(pathpronodf$NIPIpv < 0.02 &pathpronodf$RICKIpv< 0.02 ),]

  # sapply(strsplit(sigpronodf$gene,"\\."),\(x)x[1])



  print(
    read_pptx() %>%
    add_slide(layout = "Blank", master = "Office Theme") %>%
    ph_with(dml({
      # plot2pvlog(sigpronodf$NIPIpv,sigpronodf$RICKIpv,minlogpv=5,main="Immune populations")
      plot2pvdf(sigpronodf,minlogpv=4,main="TME")
    },pointsize=8),ph_location(left=0,top=0,width=3,height=3)) %>%
    ph_with(dml({
      pgpronodf=gpronodf[which(gpronodf$gene %in%unlist(newPATHS)),]
      plot2pvdf(pgpronodf,minlogpv=6,main="Genes")
      subgpronodf=pgpronodf[which(pgpronodf$NIPIpv < 0.01 &pgpronodf$RICKIpv< 0.01 ),]
      # text(x=log10(subgpronodf$NIPIpv),y=log10(subgpronodf$RICKIpv),subgpronodf$gene,adj=1.1)
    },pointsize=8),ph_location(left=3.4,top=3.5,width=3,height=3)) %>%
    ph_with(dml({
      # plot2pvlog(pathpronodf$NIPIpv,pathpronodf$RICKIpv,minlogpv=5,main="Pathway")
      plot2pvdf(pathpronodf,minlogpv=5,main="Pathway")
      subppronodf=pathpronodf[which(pathpronodf$NIPIpv < 0.01 &pathpronodf$RICKIpv< 0.01 ),]
    },pointsize=8),ph_location(left=0,top=3.2,width=3,height=3))

    , target = "plot/CompareRickiNipicolCoxPv.pptx"
  )






  print(
    read_pptx() %>%
    add_slide(layout = "Blank", master = "Office Theme") %>%
    ph_with(dml({
      plot2pvdf(sigpronodf,apre="NIPI",bpre="RICKI2",minlogpv=4,main="TME")

    },pointsize=8),ph_location(left=0,top=0,width=3,height=3)) %>%
    ph_with(dml({
      pgpronodf=gpronodf[which(gpronodf$gene %in%unlist(newPATHS)),]
      plot2pvdf(pgpronodf,apre="NIPI",bpre="RICKI2",minlogpv=7,main="Genes")
    },pointsize=8),ph_location(left=3.4,top=3.5,width=3,height=3)) %>%
    ph_with(dml({
      plot2pvdf(pathpronodf,apre="NIPI",bpre="RICKI2",minlogpv=4,main="Pathway")
    },pointsize=8),ph_location(left=0,top=3.2,width=3,height=3))

    , target = "plot/CompareRickiNipicolCoxPvCOMBOONLY.pptx"
  )
