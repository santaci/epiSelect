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
library(ggpubr)
args = commandArgs(trailingOnly=TRUE)

setwd("/home/lnd775/data/plague/results/plots/sweden/add/VAA08/")
files<-vector()
# Specify "rec" or "add" model
#filep <- list.files(pattern = c("0.txt"),path=c("58796_58800","58799_58799"),full.names = TRUE)

filep <- list.files(pattern = c("freq_"),path=c("58796_58800","58799_58799"),full.names = TRUE)
files<-append(filep,files)

all_freq <- data.frame()
for (x in files) { 
  the_freqs<-read.table(x,sep="\t",header = FALSE,col.names = c("GEN","N","fA","sim","recov","type"))
  y<-strsplit(x,split = "/")[[1]][1]
  x<-strsplit(x,split = "/")[[1]][2]
  
  model<-paste(strsplit(x,split = "_")[[1]][1:4],collapse="_")
  Ninds<-strsplit(x,split = "_")[[1]][7]
  aimfA<-paste(strsplit(x,split = "[_f]")[[1]][5])
  the_freqs$aimfA<-aimfA
  the_freqs<-separate(the_freqs,col=type,sep="_", into = c("type","rNum"), remove=TRUE)
  
  all_freq<-rbind(all_freq,the_freqs)
}

# Check if the True allele frequencies are actually there
if (any(all_freq$type=="Full")){
all_freq[all_freq$type=="Full",]$rNum<-0 
}

#Formatting for plague
if(any(all_freq$GEN=="58799")) {
  deadalive=unique(subset(all_freq,all_freq$GEN=="58799"))
  deadalive$GEN[deadalive$GEN=="58799"]<-"58799a"
  all_freq=rbind(subset(all_freq,all_freq$GEN!="58799"),deadalive) 
  
  ## Estimated fAs
  long_freqs=reshape(subset(all_freq,all_freq$type!="Full"),direction = "wide",idvar=c("sim","rNum","N","aimfA"),timevar="GEN",drop = as.vector(c("type","recov")))
  long_freqs$BADiff=long_freqs$fA.58800-long_freqs$fA.58796
  long_freqs$DADiff=long_freqs$fA.58799a-long_freqs$fA.58799d
  t1=dcast(subset(long_freqs,long_freqs$N!=1000),N+aimfA~.,value.var = "BADiff",fun.aggregate = mean)
  colnames(t1)=c("N","aimfA","BADelta")
  t3=dcast(subset(long_freqs,long_freqs$N!=1000),N+aimfA~.,value.var = "DADiff",fun.aggregate = mean)
  colnames(t3)=c("N","aimfA","DADelta")
  my_estimates=merge(t1,t3)
  
  ## True fAs
  if (any(all_freq$type=="Full")){
    long_freqs=reshape(subset(all_freq,all_freq$type=="Full"),direction = "wide",idvar=c("sim","aimfA"),timevar="GEN",drop = as.vector(c("type","recov","rNum","N")))
    long_freqs$BADeltaT=long_freqs$fA.58800-long_freqs$fA.58796
    #long_freqs$DADeltaT=long_freqs$fA.58799a-long_freqs$fA.58799d
    my_true=dcast(long_freqs,aimfA~.,value.var = "BADeltaT",fun.aggregate = mean)
   # my_true=rbind(my_true,dcast(long_freqs,aimfA~.,value.var = "DADeltaT",fun.aggregate = mean))
  }
 
  
} else {
  ## Estimated fAs
long_freqs=reshape(subset(all_freq,all_freq$type!="Full"),direction = "wide",idvar=c("sim","rNum","N","aimfA"),timevar="GEN",drop = as.vector(c("type","recov")))
long_freqs$BADiff=long_freqs$fA.58800-long_freqs$fA.58796
t1=dcast(subset(long_freqs,long_freqs$N!=1000),N+sim+aimfA~.,value.var = "BADiff",fun.aggregate = mean)
colnames(t1)=c("N","sim","aimfA","BADelta")
#t2=dcast(subset(long_freqs,long_freqs$N!=1000),N+sim+aimfA~.,value.var = "BPDiff",fun.aggregate = mean)
#colnames(t2)=c("N","sim","aimfA","BPDelta")
#my_estimates=merge(t1,t2)
my_estimates = t1

# This is used to copy into a spreadsheet along with other scenarios
kbl(my_estimates,format = 'simple')

## True fAs
if (any(all_freq$type=="Full")){
long_freqs=reshape(subset(all_freq,all_freq$type=="Full"),direction = "wide",idvar=c("sim","aimfA"),timevar="GEN",drop = as.vector(c("type","recov","rNum","N")))
long_freqs$BADeltaT=long_freqs$fA.58800-long_freqs$fA.58796
#long_freqs$BPDeltaT=long_freqs$fA.40017-long_freqs$fA.58796
my_trues=long_freqs[,c(1,2,6,7)]
my_trues$N="True"}
}

# Merge true and estimate fA datasets for plotting
deltas_plot=merge(my_trues,my_estimates,all = TRUE,by = c("sim","aimfA"))
deltas_plot$N.y<-as.factor(as.character(deltas_plot$N.y))

# Create tables 
# Viability plots using summary table
freq_power<-read.table("~/data/plague/results/plots/add_pwr_freq.out",header=TRUE,sep="\t")
freq_power=subset(freq_power,freq_power$N!=1000)
freq_power$aimfA<-as.factor(as.character(freq_power$aimfA))
freq_power$N<-as.factor(as.character(freq_power$N))
freq_power$N<-factor(freq_power$N,levels=c("20","50","100","200","500"))
freq_power$Scenario<-as.factor(as.character(freq_power$Scenario))
freq_power$Scenario<-factor(freq_power$Scenario,levels=c("DA","BA"))
freq_power$Scenario<-factor(freq_power$Scenario,levels=c("DA","BA"),labels = c("Dead v Survivors","Before v After"))
freq_power$Viability<-as.factor(as.character(freq_power$Viability))
freq_power$Viability<-factor(freq_power$Viability,levels = rev(levels(freq_power$Viability)))

ggline(freq_power,"FST","DeltaMean",color="N",shape = "aimfA",numeric.x.axis = TRUE,point.size=2,alpha=0.6,palette="aaas",facet.by=c("Viability","Scenario")) + labs(x=expression(paste("Power (",italic('F'['ST']),")" )),y=expression(paste("Estimated ",Delta,italic(bar(fA)))),shape="initial fA",color="N")
ggline(freq_power,"jSFS","DeltaMean",color="N",shape = "aimfA",numeric.x.axis = TRUE,point.size=2,alpha=0.6,palette="aaas",facet.by=c("Viability","Scenario")) + labs(x="Power (jSFS)",y=expression(paste("Estimated ",Delta,italic(bar(fA)))),shape="initial fA",color="N")
ggline(subset(freq_power,freq_power$Scenario=="DA"),"GWAS","DeltaMean",color="N",shape = "aimfA",palette="aaas",numeric.x.axis = TRUE,point.size=2,alpha=0.6,facet.by=c("Viability","Scenario")) + labs(x="Power (GWAS)",y=expression(paste("Estimated ",Delta,italic(bar(fA)))),shape="initial fA",color="N")

pp<-ggline(freq_power,"FST","aimfA",color="N",numeric.x.axis = TRUE,point.size=2,palette="aaas",facet.by=c("Viability","Scenario"), scales="free" ) + labs(color=expression(paste(italic('n'))),x=expression(paste("Power (",italic('F'['ST']),")" )),y=expression(paste("Starting ",italic('f'['A']))))
pp$layers[[1]]$aes_params$alpha <- 0.7
pp$layers[[2]]$aes_params$alpha <- 0.7


ggscatter(subset(freq_power,freq_power$Scenario=="BA"),"FST","DeltaMean",color="N",shape = "aimfA",alpha=0.6,facet.by="Viability") + labs(x="Power (FST)",y=expression(paste("Estimated ",Delta,italic(bar(fA)))),shape="initial fA",color="N")
ggscatter(subset(freq_power,freq_power$Scenario=="BP"),"FST","DeltaMean",color="N",shape = "aimfA",alpha=0.6,facet.by="Viability") + labs(x="Power (FST)",y=expression(paste("Estimated ",Delta,italic(bar(fA)))),shape="initial fA",color="N")
ggscatter(subset(freq_power,freq_power$Scenario=="DA"),"FST","DeltaMean",color="N",shape = "aimfA",alpha=0.6,facet.by="Viability") + labs(x="Power (FST)",y=expression(paste("Estimated ",Delta,italic(bar(fA)))),shape="initial fA",color="N")


ggscatter(subset(freq_power,freq_power$Scenario=="BA"),"FST","Delta.Mean",color="N",shape = "Viability",alpha=0.6,facet.by="aimfA") + labs(x="Power (FST)",y=expression(paste("Estimated ",Delta,italic(bar(x)))),shape="Viability",color="N")
ggscatter(subset(freq_power,freq_power$Scenario=="BP"),"FST","Delta.Mean",color="N",shape = "aimfA",alpha=0.6,facet.by="Viability") + labs(x="Power (FST)",y=expression(paste("Estimated ",Delta,italic(bar(x)))),shape="initial fA",color="N")
ggscatter(subset(freq_power,freq_power$Scenario=="DA"),"FST","Delta.Mean",color="N",shape = "Viability",alpha=0.6,facet.by="aimfA") + labs(x="Power (FST)",y=expression(paste("Estimated ",Delta,italic(bar(x)))),shape="Viability",color="N")



