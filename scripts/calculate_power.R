#!/usr/bin/env Rscript

cat("\nepiSelect:\nCalculate Power\n\n")

.libPaths(c("/home/lnd775/R/x86_64-pc-linux-gnu-library/4.2.2/",.libPaths()))
library(utils)
args = commandArgs(trailingOnly=TRUE)

# ESTABLISH WORKING DIRECTORY AND PARAMETERS.
#Working directory with key generations of interest in it.
workd<-args[1]
# Top Ten turns to Top 3 or whatever top
the_top<-args[2]

if(length(args) < 1) {
  stop("Usage: Rscript calculate_power.R <path to candidate tables> [# of top candidates. Default = 3].")
}

library(ggplot2)
library(scales)
library(tidyr)
library(lattice)
library(latticeExtra)
library(patchwork)
library(gtools)
library(data.table)
library(ggsci)
library(kableExtra)

if (missing(the_top) | is.null(the_top) | is.na(the_top)) {
  cat("You didn't provide a top number of candidates. Using top 3 as default.")
  the_top=3
}

setwd(paste0(workd))

# Single-level up recursive list of directories for generation comparisons
list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n-1)
    c(res, add)
  } else {
    res
  }
}

# Get top candidates files from directories
  files<-vector()
  filep <- list.files(pattern = c("topCand_mastertable.txt"),path=list.dirs.depth.n(".", n = 1),full.names = TRUE)
  filep <- filep[ grepl("topCand", filep) ]
  files<-append(filep,files)
  
  all_add <- data.frame()
  for (x in files) { 
    topTen<-read.table(x,sep="\t",header = FALSE,strip.white = TRUE)
    x<-gsub("^./","",x)
    y<-strsplit(x,split = "/")[[1]][1]
    x<-strsplit(x,split = "/")[[1]][2]
    
    model<-paste(strsplit(x,split = "_")[[1]][1:4],collapse="_")
    if (length(strsplit(x,split = "_")[[1]]) < 9) {
      Ninds<-strsplit(x,split = "_")[[1]][6]
    } else{
      Ninds<-paste(collapse = "_",strsplit(x,split = "_")[[1]][6:7])
    }
    colnames(topTen)[9]<-"V14"
    colnames(topTen)[10]<-"V15"
    colnames(topTen)[11]<-"V16"
    colnames(topTen)[12]<-"V9"
    topTen$V9<-factor(topTen$V9, levels = c("One","Five","Ten","None"))
    topTen<-separate(topTen,col=V7,sep="_", into = c("V10","V11"), remove= FALSE)
    topTen<-separate(topTen,col=V10,sep="chr", into = c("V12","V13"), remove=TRUE)
    topTen<-subset(topTen,select = -V12)
    topTen$V1<-gsub("chr","",topTen$V1)
    topTen$V1<-factor(topTen$V1)
    topTen$V5<-factor(topTen$V5)
    topTen<-subset(topTen,as.integer(as.character(topTen$V5))<=the_top)
    topTen$V6<-factor(topTen$V6)
    topTen<-data.table(topTen)
    
    ## These are the hard-coded conditions: 
    ## An r2 >= 0.8 (high LD) with selected variant.
    r2=0.8
    topTen$V9[as.character(topTen$V3)!=as.character(topTen$V7) & topTen$V1==topTen$V13 & topTen$V16 >= r2] <-"One"
    
    #Category 1
    topTen$V15<-round(topTen$V15,digits = 2)
    topTen$V17<-y
    topTen$V18<-Ninds
    all_add<-rbind(all_add,topTen)
  }
  
   # Two-colour scheme
  custom.col<-c("beige","coral")
  
  fst_cat1_cat2<-all_add[V8=="FST" & V9=="None",.N, by=.(V7,V15,V6,V14,V17,V18)]
  fst_cat1_cat2<-fst_cat1_cat2[order(V15)]
  fst_cat1_cat2$an<-"FST"
  
  jsfs_cat1_cat2<-all_add[V8=="JSFS" & V9=="None",.N, by=.(V7,V15,V6,V14,V17,V18)]
  jsfs_cat1_cat2<-jsfs_cat1_cat2[order(V15)]
  jsfs_cat1_cat2$an<-"JSFS"
  ans_cat1_cat2<-rbind(jsfs_cat1_cat2,fst_cat1_cat2)
  
  if ("iHS" %in% all_add$V8) {
    ihs_cat1_cat2<-all_add[V8=="iHS" & V9=="None",.N, by=.(V7,V15,V6,V14,V17,V18)]
    ihs_cat1_cat2<-ihs_cat1_cat2[order(V15)]
    ihs_cat1_cat2$an<-"iHS"
    ans_cat1_cat2<-rbind(ans_cat1_cat2,ihs_cat1_cat2)
  }
  if ("GWAS" %in% all_add$V8) {
    gwas_cat1_cat2<-all_add[V8=="GWAS" & V9=="None",.N, by=.(V7,V15,V6,V14,V17,V18)]
    gwas_cat1_cat2<-gwas_cat1_cat2[order(V15)]
    gwas_cat1_cat2$an<-"GWAS"
    ans_cat1_cat2<-rbind(ans_cat1_cat2,gwas_cat1_cat2)
  } 
  if ("TD" %in% all_add$V8) {
    td_cat1_cat2<-all_add[V8=="TD" & V9=="None",.N, by=.(V7,V15,V6,V14,V17,V18)]
    td_cat1_cat2<-td_cat1_cat2[order(V15)]
    td_cat1_cat2$an<-"TD"
    ans_cat1_cat2<-rbind(ans_cat1_cat2,td_cat1_cat2)
  }
  if ("PI" %in% all_add$V8) {
    pi_cat1_cat2<-all_add[V8=="PI" & V9=="None",.N, by=.(V7,V15,V6,V14,V17,V18)]
    pi_cat1_cat2<-pi_cat1_cat2[order(V15)]
    pi_cat1_cat2$an<-"PI"
    ans_cat1_cat2<-rbind(ans_cat1_cat2,pi_cat1_cat2)
  }
  mean_ans<-ans_cat1_cat2[, .(mean_sim = mean(N)), by =.(an,V17,V18,V15)]

## Needs to output a table as a tab delimited output and a plotted image.  
  ## CREATE LATEX TABLES ##
  mean_ans$V18<-factor(mean_ans$V18,levels=mixedsort(levels(factor(mean_ans$V18))))
  mean_ans$V17<-factor(mean_ans$V17,levels=mixedsort(levels(factor(mean_ans$V17))))
  mean_ans$perc<-(the_top-mean_ans$mean_sim)*100
  mean_ans=mean_ans[order(V15,V18),]
  mean_ans$an<-factor(mean_ans$an,levels = c("FST","JSFS","GWAS"))
  dcast(mean_ans,V18+V15~V17+an, value.var = "perc")
  power_table=dcast(mean_ans,V18+V15~V17+an, value.var = "perc")
  colnames(power_table)[1]<-c("samples")
  colnames(power_table)[2]<-c("fA")
  write.table(power_table,file = paste0(model,"_power.txt"),quote = FALSE,sep = "\t",row.names = FALSE)
  

  ## PLOTTING FOR SIMULATIONS
  #For topX - BIG FONT
  theme.size = 4 
  geom.text.size = theme.size/(14/5)  
  
  power_plot = ggplot(mean_ans,aes(x=factor(an),y=V17,fill=(the_top-mean_sim))) +
    geom_tile(alpha=0.8) + geom_text(size=geom.text.size,aes(x=an,y=V17,label=paste0(round((the_top-mean_sim)*100, digits=1),"%"))) + facet_grid(cols = vars(V15), rows = vars(V18))+ coord_equal() + theme_bw() + theme(axis.text.x = element_text(angle=90), panel.grid.major = element_blank(), text = element_text(size=theme.size), legend.text=element_text(size=3),strip.background=element_rect(fill="white"), legend.key.size = unit(0.25, 'cm')) +
    labs(x="",y="",fill=paste0("Found in Top ",the_top)) + scale_fill_gradientn(colours = custom.col, breaks = seq(0,1,0.1), limits=c(0,1))
  ggsave(plot = power_plot, filename = paste0(model,"_power.png"),device = "png",units = "px", height = 900, width=1200)
