.libPaths(c("/home/lnd775/R/x86_64-pc-linux-gnu-library/4.2.2/",.libPaths()))
library(utils)
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
args = commandArgs(trailingOnly=TRUE)

# ESTABLISH WORKING DIRECTORY AND PARAMETERS
workd<-args[1]
filep<-args[2]
nsim<-args[3]
# nsamp is resampling not sample size
nsamp<-args[4]
analysis<-args[5]

# Top Ten turns to Top 3 or whatever top
the_top<-3

setwd(paste0(workd))
model<-paste(strsplit(filep,split = "_")[[1]][1:4],collapse="_")
Ninds<-strsplit(filep,split = "_")[[1]][6]

  setwd("/home/lnd775/data/plague/results/plots/sweden/rec/VAA1/")
  files<-vector()
  # Specify "rec" or "add" model
  filep <- list.files(pattern = c("topTen_mastertable.txt"),path=c("58796_58800","58799_58799"),full.names = TRUE)
  filep <- filep[ grepl("topTen", filep) ]
  #filep <- list.files(pattern = c("add","topTen_mastertable.txt"),path=c("58284_58293","58292_58292","58292_58293"),full.names = TRUE)
  #filep <- list.files(pattern = c("add","topTen_mastertable.txt"),path=c("58292_58292"),full.names = TRUE)
  files<-append(filep,files)
  
  all_add <- data.frame()
  for (x in files) { 
    topTen<-read.table(x,sep="\t",header = FALSE,strip.white = TRUE)
    y<-strsplit(x,split = "/")[[1]][1]
    x<-strsplit(x,split = "/")[[1]][2]
    
    model<-paste(strsplit(x,split = "_")[[1]][1:4],collapse="_")
    Ninds<-strsplit(x,split = "_")[[1]][6]
    #Ninds2<-strsplit(x,split = "_")[[1]][6]
    colnames(topTen)[9]<-"V14"
    colnames(topTen)[10]<-"V15"
    colnames(topTen)[11]<-"V16"
    colnames(topTen)[12]<-"V9"
    #colnames(topTen)[13]<-"V17"
    topTen$V9<-factor(topTen$V9, levels = c("One","Five","Ten","None"))
    topTen<-separate(topTen,col=V7,sep="_", into = c("V10","V11"), remove= FALSE)
    topTen<-separate(topTen,col=V10,sep="chr", into = c("V12","V13"), remove=TRUE)
    topTen<-subset(topTen,select = -V12)
    topTen$V1[topTen$V1 == "chr21"] <-"21"
    topTen$V1[topTen$V1 == "chr22"] <-"22"
    topTen$V1<-factor(topTen$V1)
    topTen$V5<-factor(topTen$V5)
    
    topTen<-subset(topTen,as.integer(as.character(topTen$V5))<=the_top)
    topTen$V6<-factor(topTen$V6)
    topTen<-data.table(topTen)
    topTen$V9[as.character(topTen$V3)!=as.character(topTen$V7) & topTen$V1==topTen$V13 & topTen$V16 >= 0.8] <-"One"
    
    #Category 1
    topTen$V15<-round(topTen$V15,digits = 2)
    topTen$V17<-y
    topTen$V18<-Ninds
    all_add<-rbind(all_add,topTen)
    #all_add<-na.omit(all_add)
  }
  
   # Two-colour light
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
  mean_ans<-ans_cat1_cat2[, .(mean_sim = mean(N)), by =.(an,V17,V18,V15)]
  
  # CREATE LATEX TABLES
  tmp=unlist(strsplit(mean_ans$V18,"n"))
  tmp<-tmp[tmp!=""]
  mean_ans$V18<-as.integer(tmp)
  mean_ans$perc<-(the_top-mean_ans$mean_sim)*100
  mean_ans=mean_ans[order(V15,V18),]
  #kbl(mean_ans[,c(1,3:4,6)], format = "latex",caption = "Before versus Right after comparisons comparing methods, initial fAs, and sampling sizes.",col.names = c("Method","N","fA","% Detected"))
  #kbl(dcast(mean_ans,V18+V15~V17+an, value.var = "perc"), format = "latex",caption = "Before versus Right after comparisons comparing methods, initial fAs, and sampling sizes.",col.names = c("N","fA","% Detected\n (FST)","% Detected\n (JSFS)","iHS"))
  mean_ans$an<-factor(mean_ans$an,levels = c("FST","JSFS","GWAS"))
  dcast(mean_ans,V18+V15~V17+an, value.var = "perc")
  
  
    ### Change for sweden plots
mean_ans$V18<-factor(mean_ans$V18,levels = c("n20","n50","n100","n200","n500"))
  
mean_ans$V17[mean_ans$V17=="58796_58800"]<-"BeforeAfter"
mean_ans$V17[mean_ans$V17=="58799_58799"]<-"DeadAlive"

  
  ## LATEST PLOTTING FOR SIMULATIONS
  #For topX - BIG FONT
  ggplot(mean_ans,aes(x=factor(an),y=V17,fill=(the_top-mean_sim))) +
    geom_tile(alpha=0.8) + geom_text(aes(x=an,y=V17,label=paste0(round((the_top-mean_sim)*100, digits=1),"%"))) + facet_grid(cols = vars(V15), rows = vars(V18))+ coord_equal() + theme_bw() + theme(axis.text.x = element_text(angle=90), panel.grid.major = element_blank(), text = element_text(size=25), legend.text=element_text(size=12),strip.background=element_rect(fill="white"), legend.key.size = unit(1, 'cm')) +
    labs(x="",y="",fill=paste0("Found in Top ",the_top)) + scale_fill_gradientn(colours = custom.col, breaks = seq(0,1,0.1), limits=c(0,1))

  
  
  
  ### Change for london plots only
  mean_ans$V17[mean_ans$V17=="58284_58293"]<-"Before_RightAfter"
  mean_ans$V17[mean_ans$V17=="58292_58292"]<-"Dead_Survivor"
  mean_ans$V17[mean_ans$V17=="58292_58293"]<-"Dead_RightAfter"
  mean_ans$V18[mean_ans$V18=="n42"]<-"n42_n63"
  mean_ans$V18[mean_ans$V18=="n38"]<-"n38_n63"