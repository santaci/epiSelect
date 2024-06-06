#!/usr/bin/env Rscript

cat("\nepiSelect:\nGet Top Candidates\n\n")

.libPaths(c("/home/lnd775/R/x86_64-pc-linux-gnu-library/4.2.2/",.libPaths()))

library(this.path)
options(warn=-1, keep.source=T)
args = commandArgs(trailingOnly=TRUE)

# Get source file for Manhattan plots
script.dir=this.dir()
source(paste0(script.dir,"/src/ManhattanQQ.R"))

# Path of working directory
i <-args[1]

# Relevant file in working directory depending on analysis
j <-args[2]

# Can choose from fst, gwas, ihs_before, ihs_after, nsl
analysis <- args[3]

# Known selected variant in simulations 
selected<- args[4]

# Generations being compared
gen1<-args[5]
gen2<-args[6]

thetop<-args[7]

if(length(args) < 1) {
  stop("Usage: Rscript get_top.R <path to scenario wd> <relevant file for analysis> <analysis> <selected_variant> <gen1> <gen2> [# of top candidates. Default = 10].\n Analysis can be: fst, gwas, ihs_before, ihs_after, nsl" )
}

if (missing(thetop) | is.null(thetop) | is.na(thetop)) {
  cat("You didn't provide a top number of candidates. Using top 10 as default.")
  thetop=c(10)
} else {
  thetop=as.integer(thetop)
}

library(utils)
library(ggplot2)
library(scales)
library(tidyr)
library(lattice)
library(latticeExtra)
library(ggpointdensity)

# Set the working directory to the path specified in i
setwd(paste0(i))

## Candidates cannot be within 1 MB of each other.

## FST FUNCTION ##
run_fst<-function(df_fst) {
  #Find top ten
  fst_val_index<-sort(df_fst$FST,index.return=TRUE,decreasing = TRUE)
  top_idx<-fst_val_index$ix[1]
  topTen<-df_fst[top_idx,]
  new_df<-df_fst[-top_idx,]
  fst_ind = sort(new_df$FST,index.return=TRUE, decreasing = TRUE)
  fst_ind = fst_ind$ix
  #fst_ind<-fst_val_index$ix
  prev<-topTen
  for (y in 1:length(fst_ind)){
    minflank=prev$POS - 1000000
    maxflank=prev$POS + 1000000
    new_df<-new_df[!is.na(new_df$FST),]
    y<-fst_ind[1]
    nex<-new_df[y,]
    new_df<-new_df[-y,]
    cat("Previous hit was: ",as.character(prev$SNP),". \n")
    cat("New hit is: ",as.character(nex$SNP),". \n")
    if (length(topTen$POS) < thetop) {
      if (nex$CHR==prev$CHR & nex$POS > minflank & nex$POS < maxflank) {
        new_df <- subset(new_df,new_df$CHR != prev$CHR | (new_df$CHR == prev$CHR & new_df$POS < minflank | new_df$POS > maxflank))
        fst_ind = sort(new_df$FST,index.return=TRUE, decreasing = TRUE)
        fst_ind = fst_ind$ix
        cat("Too close to SNP",as.character(paste0(prev$CHR,"_",prev$POS))," : ",abs(prev$POS-nex$POS)," bp away.\n")
      }

      else {
        prev<-nex
        topTen<-rbind(topTen,nex)
        fst_ind = sort(new_df$FST,index.return=TRUE, decreasing = TRUE)
        fst_ind = fst_ind$ix
      }
    }
    else
      break
  }
  return(topTen)
}

## JSFS FUNCTION ##
run_jsfs<-function(df_jsfs) {
  #Find top candidates
  jsfs_val_index<-df_jsfs[order(-df_jsfs$DIFF, df_jsfs$DEN), ][c("DEN","ix")]
  top_idx<-jsfs_val_index$ix[1]
  topTen<-df_jsfs[top_idx,]
  new_df<-df_jsfs[-top_idx,]
  jsfs_ind = sort(new_df$DEN,index.return=TRUE)
  jsfs_ind = jsfs_ind$ix
  #jsfs_ind<-jsfs_val_index$ix
  prev<-topTen
  for (y in 1:length(jsfs_ind)){
    minflank=prev$POS - 1000000
    maxflank=prev$POS + 1000000
    y<-jsfs_ind[1]
    nex<-new_df[y,]
    new_df<-new_df[-y,]
    cat("Previous hit was: ",as.character(prev$SNP),". \n")
    cat("New hit is: ",as.character(nex$SNP),". \n")
    if (length(topTen$POS) < thetop) {
      if (nex$CHR==prev$CHR & nex$POS > minflank & nex$POS < maxflank) {
        new_df <- subset(new_df,new_df$CHR != prev$CHR | (new_df$CHR == prev$CHR & new_df$POS < minflank | new_df$POS > maxflank))
        jsfs_ind = sort(new_df$DEN,index.return=TRUE)
        jsfs_ind = jsfs_ind$ix
        
        cat("Too close to SNP",as.character(paste0(prev$CHR,"_",prev$POS))," : ",abs(prev$POS-nex$POS)," bp away.\n")
      }
      
      else {
        prev<-nex
        topTen<-rbind(topTen,nex)
        jsfs_ind = sort(new_df$DEN,index.return=TRUE)
        jsfs_ind = jsfs_ind$ix
      }
    }
    else
      break
  }
  return(topTen)
}

## GWAS FUNCTION ##
run_gwas<- function(df_gwas) {
  #Find top ten
  gwas_val_index<-sort(df_gwas$P,index.return=TRUE)
  top_idx<-gwas_val_index$ix[1]
  topTen<-df_gwas[top_idx,]
  new_df<-df_gwas[-top_idx,]
  gwas_ind = sort(new_df$P,index.return=TRUE)
  gwas_ind = gwas_ind$ix
  #gwas_ind<-gwas_val_index$ix
  prev<-topTen
  for (y in 1:length(gwas_ind)){
    minflank=prev$BP - 1000000
    maxflank=prev$BP + 1000000
    new_df<-new_df[!is.na(new_df$P),]
    y<-gwas_ind[1]
    nex<-new_df[y,]
    new_df<-new_df[-y,]
    if (length(topTen$BP) < thetop) {
      if (nex$CHR==prev$CHR & nex$BP > minflank & nex$BP < maxflank) {
        new_df <- subset(new_df,new_df$CHR != prev$CHR | (new_df$CHR == prev$CHR & new_df$BP < minflank | new_df$BP > maxflank))
        gwas_ind = sort(new_df$P,index.return=TRUE)
        gwas_ind = gwas_ind$ix
        cat("Too close to SNP",as.character(prev$SNP)," : ",abs(prev$BP-nex$BP)," bp away.\n")
      }

      else {
        prev<-nex
        topTen<-rbind(topTen,nex)
        gwas_ind = sort(new_df$P,index.return=TRUE)
        gwas_ind = gwas_ind$ix
      }
    }
    else
      break
  }
  return(topTen)
}


## IHS FUNCTION ##
run_ihs <- function(df_ihs) {
  #Find top ten
  ihs_val_index<-sort(df_ihs$V7,index.return=TRUE,decreasing = TRUE)
  top_idx<-ihs_val_index$ix[1]
  topTen<-df_ihs[top_idx,]
  new_df<-df_ihs[-top_idx,]
  ihs_ind = sort(new_df$V7,index.return=TRUE,decreasing = TRUE)
  ihs_ind = ihs_ind$ix
  #ihs_ind<-ihs_val_index$ix
  prev<-topTen
  for (y in 1:length(ihs_ind)){
    minflank=prev$V2 - 1000000
    maxflank=prev$V2 + 1000000
    new_df<-new_df[!is.na(new_df$V7),]
    y<-ihs_ind[1]
    nex<-new_df[y,]
    new_df<-new_df[-y,]
    cat("y now",y,"\n")
    if (length(topTen$V2) < thetop) {
      if (nex$CHR==prev$CHR & nex$V2 > minflank & nex$V2 < maxflank) {
        new_df <- subset(new_df,new_df$CHR != prev$CHR | (new_df$CHR == prev$CHR & new_df$V2 < minflank | new_df$V2 > maxflank))
        ihs_ind = sort(new_df$V7,index.return=TRUE,decreasing = TRUE)
        ihs_ind = ihs_ind$ix
        cat("Too close to SNP",as.character(prev$V1)," : ",abs(prev$V2-nex$V2)," bp away.\n")
      }

      else {
        prev<-nex
        topTen<-rbind(topTen,nex)
        ihs_ind = sort(new_df$V7,index.return=TRUE,decreasing = TRUE)
        ihs_ind = ihs_ind$ix
      }
    }
    else
      break
  }
  return(topTen)
}
## Tajima's D FUNCTION ##
run_td<-function(df_td) {
  #Find top ten
  td_val_index<-sort(df_td$abTD,index.return=TRUE,decreasing = TRUE)
  top_idx<-td_val_index$ix[1]
  topTen<-df_td[top_idx,]
  new_df<-df_td[-top_idx,]
  td_ind=sort(new_df$abTD,index.return=TRUE,decreasing = TRUE)
  td_ind=td_ind$ix
  prev<-topTen
  for (y in 1:length(td_ind)){
    minflank=prev$POS - 1000000
    maxflank=prev$POS + 1000000
    y<-td_ind[1]
    nex<-new_df[y,]
    new_df<-new_df[-y,]
    print(length(new_df$SNP))
    cat("Previous hit was: ",as.character(prev$SNP),". \n")
    cat("New hit is: ",as.character(nex$SNP),". \n")
    if (length(topTen$POS) < 10) {
      if (nex$CHR==prev$CHR & nex$POS > minflank & nex$POS < maxflank) {
        new_df <- subset(new_df,new_df$CHR != prev$CHR | (new_df$CHR == prev$CHR & new_df$POS < minflank | new_df$POS > maxflank))
        td_ind = sort(new_df$abTD,index.return=TRUE, decreasing = TRUE)
        td_ind = td_ind$ix
        cat("Too close to SNP",as.character(paste0(prev$CHR,"_",prev$POS))," : ",abs(prev$POS-nex$POS)," bp away.\n")
      }
      
      else {
        prev<-nex
        topTen<-rbind(topTen,nex)
        td_ind = sort(new_df$abTD,index.return=TRUE, decreasing = TRUE)
        td_ind = td_ind$ix
      }
    }
    else
      break
  }
  return(topTen)
}

## Pi FUNCTION ##
run_pi<-function(df_pi) {
  #Find top ten
  pi_val_index<-sort(df_pi$PI,index.return=TRUE,decreasing = TRUE)
  top_idx<-pi_val_index$ix[1]
  topTen<-df_pi[top_idx,]
  new_df<-df_pi[-top_idx,]
  print(length(new_df$SNP))
  pi_ind=sort(new_df$PI,index.return=TRUE,decreasing = TRUE)
  pi_ind=pi_ind$ix
  print(length(pi_ind))
  
  prev<-topTen
  for (y in 1:length(pi_ind)){
    print(length(pi_ind))
    minflank=prev$POS - 1000000
    maxflank=prev$POS + 1000000
    y<-pi_ind[1]
    nex<-new_df[y,]
    new_df<-new_df[-y,]
    print(length(new_df$SNP))
    print(y)
    
    cat("Previous hit was: ",as.character(prev$SNP),". \n")
    cat("New hit is: ",as.character(nex$SNP),". \n")
    if (length(topTen$POS) < 10) {
      if (nex$CHR==prev$CHR & nex$POS > minflank & nex$POS < maxflank) {
        new_df <- subset(new_df,new_df$CHR != prev$CHR | (new_df$CHR == prev$CHR & new_df$POS < minflank | new_df$POS > maxflank))
        pi_ind = sort(new_df$PI,index.return=TRUE, decreasing = TRUE)
        pi_ind = pi_ind$ix
        cat("Too close to SNP",as.character(paste0(prev$CHR,"_",prev$POS))," : ",abs(prev$POS-nex$POS)," bp away.\n")
      }
      
      else {
        prev<-nex
        topTen<-rbind(topTen,nex)
        pi_ind = sort(new_df$PI,index.return=TRUE, decreasing = TRUE)
        pi_ind = pi_ind$ix
      }
    }
    else
      break
  }
  return(topTen)
}

if (analysis == "gwas") {
  df_gwas<-read.table(gzfile(j),header=TRUE)
  df_gwas<-df_gwas[!is.na(df_gwas$P),]
  topTen <- run_gwas(df_gwas = df_gwas)
  snippies10<-as.vector(topTen$SNP)
  if (is.na(match(selected,snippies10))==FALSE) {
    snippies10<-snippies10[-match(selected,snippies10)]
    snippies10<-c(snippies10,selected)
    snp_cols<-c(rep("blue4",thetop-1),"darkred")
    snp_shape<-c(rep(20,thetop-1),11)
  } else{snippies10<-c(snippies10,selected)
  snp_cols<-c(rep("blue4",thetop),"darkred")
  snp_shape<-c(rep(20,thetop),11) }
  write.table(topTen,paste0("gwas_",gen1,"_",gen2,"_top_cand.txt"),sep='\t', quote = FALSE, row.names=FALSE)
  ann<- annotateSNPRegions(df_gwas$SNP,df_gwas$CHR,df_gwas$BP,df_gwas$P,snippies10,labels = snippies10,col=snp_cols,kbaway=0,pch=snp_shape)
  bonfer<- 0.05/(length(df_gwas$SNP)*2)
  bitmap(file = paste0("df_gwas","_",gen1,"_",gen2,".png"), res=150,height=3,width=8)
  print(manhattan.plot(df_gwas$CHR,df_gwas$BP,df_gwas$P,sig.level=bonfer,xlab=list(label="Chromosome",cex=1.5),ylab=list(label=expression(-log[10](p-value)),cex=1.5),selected=selected,annotate=ann))
  dev.off()
  bitmap(file = "df_qq.png", res=150)
  print(qqunif.plot(df_gwas$GC,col='red2'))
  dev.off()
} else if (analysis == "fst") {
  df_fst <- read.table(gzfile(j),header = TRUE, sep = '\t',na.strings = "NaN")
  df_fst <- df_fst[!is.na(df_fst$FST),]
  df_fst <-subset(df_fst,df_fst$FST >= 0)
  topTen <- run_fst(df_fst = df_fst)
  snippies10<-as.vector(topTen$SNP)
  if (is.na(match(selected,snippies10))==FALSE) {
    snippies10<-snippies10[-match(selected,snippies10)]
    snippies10<-c(snippies10,selected)
    snp_cols<-c(rep("blue4",thetop-1),"darkred")
    snp_shape<-c(rep(20,thetop-1),11)
  } else{snippies10<-c(snippies10,selected)
  snp_cols<-c(rep("blue4",thetop),"darkred")
  snp_shape<-c(rep(20,thetop),11)}
  write.table(topTen,paste0(analysis,"_",gen1,"_",gen2,"_top_cand.txt"),sep='\t', quote = FALSE, row.names=FALSE)
  fst_snps<-annotateSNPRegions(df_fst$SNP,df_fst$CHR,df_fst$POS,df_fst$FST,snippies10,labels = snippies10,col=snp_cols,kbaway=0,pch=snp_shape)
  bitmap(file = paste0("df_fst","_",gen1,"_",gen2,".png"), res=150,height=3,width=8)
  print(ihs.plot(df_fst$CHR,df_fst$POS,df_fst$FST,xlab=list(label="Chromosome",cex=1.5),ylab=list(label=expression( italic('F'['ST']) ),cex=1.5),selected=selected, annotate = fst_snps))
  dev.off()
} else if (analysis == "jsfs") {
  jsfs_raw <- read.table(gzfile(j),header = TRUE, sep = '\t',na.strings = "NaN")
  df_jsfs <- jsfs_raw[!is.na(jsfs_raw$FST),]
  plot_jsfs<-ggplot(df_jsfs,aes(x=AC1,y=AC2)) + geom_pointdensity(adjust=2) + scale_color_gradientn(colors = c("grey",viridis_pal(option = "H")(20)))

## Note this adjustment of pointdensities is fixed but can be changed based on one's data. Look at data distribution and/or range.
  jsfs_adj<-max(ggplot_build(plot_jsfs)$data[[1]]$density)/0.001 
  plot_jsfs<-ggplot(df_jsfs,aes(x=AC1,y=AC2)) + geom_pointdensity(adjust=2+jsfs_adj) + scale_color_gradientn(colors = c("grey",viridis_pal(option = "H")(20)))
  df_jsfs <-cbind(df_jsfs,ggplot_build(plot_jsfs)$data[[1]]$density)

## We are specifically looking for increases in allele frequency when comparing time-points but this can be changed (i.e. AC1 < AC2).
  df_jsfs<-subset(df_jsfs,df_jsfs$AC1 < df_jsfs$AC2)
  lastCol=ncol(df_jsfs)
  colnames(df_jsfs)[lastCol] <-"DEN"
  df_jsfs$ix<-c(1:nrow(df_jsfs))
  df_jsfs$DIFF=abs(df_jsfs$AC1-df_jsfs$AC2)
  topTen <- run_jsfs(df_jsfs = df_jsfs)
  snippies10<-as.vector(topTen$SNP)
  if (is.na(match(selected,snippies10))==FALSE) {
  snippies10<-snippies10[-match(selected,snippies10)]
  snippies10<-c(snippies10,selected)
  snp_cols<-c(rep("purple1",thetop-1),"magenta")
  snp_shape<-c(rep(20,thetop-1),18)} else{snippies10<-c(snippies10,selected)
  snp_cols<-c(rep("purple1",thetop),"magenta")
  snp_shape<-c(rep(20,thetop),18)}
  # Clean up annotated SNPs
  jsfs_snps<-df_jsfs[df_jsfs$SNP %in% snippies10,]
  if(nrow(jsfs_snps[-(which(duplicated(jsfs_snps[jsfs_snps$SNP %in% snippies10,]$SNP))),])!=0) {
  jsfs_snps<-jsfs_snps[-(which(duplicated(jsfs_snps[jsfs_snps$SNP %in% snippies10,]$SNP))),]
  }
  jsfs_snps<-jsfs_snps[match(snippies10, jsfs_snps$SNP),]
  write.table(topTen,paste0(analysis,"_",gen1,"_",gen2,"_top_cand.txt"),sep='\t', quote = FALSE, row.names=FALSE)
  plot_jsfs<-plot_jsfs + theme_bw() + geom_point(size=2, data=jsfs_snps,aes(x=AC1,y=AC2), colour=snp_cols, shape=snp_shape) + labs(x=gen1,y=gen2,color="Point-Density") 
  ggsave(paste0("df_jsfs","_",gen1,"_",gen2,".png"),plot = plot_jsfs,dpi =150,height=3,width=8)
  } else if (analysis == "tajd"){
    df_td <- read.table(gzfile(j),header = TRUE, sep = '\t',na.strings = "NaN")
    df_td <- df_td[!is.na(df_td$TD),]
    df_td$abTD <- abs(df_td$TD)
    topTen <- run_td(df_td = df_td)
    snippies10<-as.vector(topTen$SNP)
    if (is.na(match(selected,snippies10))==FALSE) {
      snippies10<-snippies10[-match(selected,snippies10)]
      snippies10<-c(snippies10,selected)
      snp_cols<-c(rep("blue4",thetop-1),"darkred")
      snp_shape<-c(rep(20,thetop-1),11)
    } else if (is.na(match(selected,df_td$SNP))==TRUE){
      snp_cols<-rep("blue4",thetop)
      snp_shape<-rep(20,10) } else {snippies10<-c(snippies10,selected)
      snp_cols<-c(rep("blue4",thetop),"darkred")
      snp_shape<-c(rep(20,thetop),11) }
    print(snippies10)
    write.table(topTen,paste0(analysis,"_",gen1,"_",gen2,"_top_cand.txt"),sep='\t', quote = FALSE, row.names=FALSE)
    td_snps<-annotateSNPRegions(df_td$SNP,df_td$CHR,df_td$POS,df_td$TD,snippies10,labels = snippies10,col=snp_cols,kbaway=0,pch=snp_shape)
    print(df_td)
    bitmap(file = paste0("df_td_",gen2,".png"), res=150,height=3,width=8)
    print(ihs.plot(df_td$CHR,df_td$POS,df_td$TD,xlab=list(label="Chromosome",cex=1.5),ylab=list(label=expression(italic('T'['D'])),cex=1.5),selected=selected, annotate = td_snps))
    dev.off()
  } else if (analysis == "pi") {
    df_pi <- read.table(gzfile(j),header = TRUE, sep = '\t',na.strings = "NaN")
    df_pi <- df_pi[!is.na(df_pi$PI),]
    df_pi <-subset(df_pi,df_pi$PI >= 0)
    topTen <- run_pi(df_pi = df_pi)
    snippies10<-as.vector(topTen$SNP)
    if (is.na(match(selected,snippies10))==FALSE) {
      snippies10<-snippies10[-match(selected,snippies10)]
      snippies10<-c(snippies10,selected)
      snp_cols<-c(rep("blue4",thetop-1),"darkred")
      snp_shape<-c(rep(20,thetop-1),11)
    } else{snippies10<-c(snippies10,selected)
    snp_cols<-c(rep("blue4",thetop),"darkred")
    snp_shape<-c(rep(20,thetop),11)}
    write.table(topTen,paste0(analysis,"_",gen1,"_",gen2,"_top_cand.txt"),sep='\t', quote = FALSE, row.names=FALSE)
    pi_snps<-annotateSNPRegions(df_pi$SNP,df_pi$CHR,df_pi$POS,df_pi$PI,snippies10,labels = snippies10,col=snp_cols,kbaway=0,pch=snp_shape)
    bitmap(file = paste0("df_pi","_",gen2,".png"), res=150,height=3,width=8)
    print(ihs.plot(df_pi$CHR,df_pi$POS,df_pi$PI,xlab=list(label="Chromosome",cex=1.5),ylab=list(label=expression(pi),cex=1.5),selected=selected, annotate = pi_snps))
    dev.off()
  } else if (analysis == "ihs_before" | analysis == "ihs_after" | analysis == "nsl") {
  df_ihs <- read.table(j,header=FALSE, sep='\t')
  df_ihs <- df_ihs[!is.na(df_ihs$V7),]
  df_ihs <- separate(df_ihs,col="V1", sep="_",into = c("CHR"),remove = FALSE)
  df_ihs$V7 <-abs(df_ihs$V7)
  # Normalized iHS
  topTen<-run_ihs(df_ihs = df_ihs)
  snippies10<-as.vector(topTen$V1)
  if (is.na(match(selected,snippies10))==FALSE) {
    snippies10<-snippies10[-match(selected,snippies10)]
    snippies10<-c(snippies10,selected)
    snp_cols<-c(rep("blue4",thetop-1),"darkred")
    snp_shape<-c(rep(20,thetop-1),11)
  } else{snippies10<-c(snippies10,selected)
  snp_cols<-c(rep("blue4",thetop),"darkred")
  snp_shape<-c(rep(20,thetop),11) }
  gen<-strsplit(strsplit(j,split = "_")[[1]][3], split="i")[[1]][1]
  write.table(topTen,paste0(analysis,"_",gen,"top_cand.txt"),sep='\t', quote = FALSE, row.names=FALSE)
  ihs_snps<-annotateSNPRegions(df_ihs$V1,df_ihs$CHR,df_ihs$V2,df_ihs$V7,maxpvalue = max(df_ihs$V7),snippies10,labels = snippies10,col=snp_cols,kbaway=0,pch=snp_shape)
  # Not in windows
  bitmap(file = paste0("df_",gen,analysis,".png"), res=150,height=3,width=8)
  if (analysis == "nsl") { print(ihs.plot(df_ihs$CHR,df_ihs$V2,df_ihs$V7,xlab=list(label="Chromosome",cex=1.5),ylab=list(label="| Standardized nSL |",cex=1.5),selected=selected,annotate = ihs_snps)) } else { print(ihs.plot(df_ihs$CHR,df_ihs$V2,df_ihs$V7,xlab=list(label="Chromosome",cex=1.5),ylab=list(label="| Standardized iHS |",cex=1.5),selected=selected,annotate = ihs_snps)) }
  dev.off()
} else{ print("You can only perform: gwas, jsfs, fst, tajd, pi, ihs or nsl.") }
