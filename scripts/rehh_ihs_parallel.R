.libPaths(c("/home/lnd775/R/x86_64-pc-linux-gnu-library/4.2.2/",.libPaths()))
library(rehh)
library(data.table)
library(R.utils)
library(vcfR)
library(vcfppR)
library(parallel)
script.dir=getSrcDirectory(function(x) {x})
source(paste0(script.dir,"/src/rehh_vcf2haplo.R"))

options(warn=-1)
args = commandArgs(trailingOnly=TRUE)

# working directory
i <-args[1]
setwd(i)
#Generation performed on
gen <-args[2]
# Number of threads
t <-args[3]
print("Making cluster for parallelization...")
cl<-makeCluster(as.numeric(t))

# vcfs in working directory
j <-list.files(pattern = c("ihs_","vcf.gz"),full.names=TRUE)
j <-j[grep(gen,j)]
j <-j[grep("vcf.gz",j)]
j<-j[!grepl("csi",j)]

print(j)
chromList<-vector()
for (k in 1:length(j)){
  chromList<-append(x = chromList, strsplit(strsplit(j,split = "_")[[k]][5], split=".v")[[1]][1])
}

extract_data_from_file <- function(file) {
  k=strsplit(file,split = "[_.]")[[1]][6]
  # code to extract data from file
  hh <- vcf2haplohh(hap_file = file,min_maf=0.05,chr.name = k,polarize_vcf = FALSE,vcf_reader = 'vcfppR')
  scan <- scan_hh(hh, threads=t)
  return(scan)
}

print("Running REHH wgscan")
wgscan <- do.call(rbind,mclapply(j, extract_data_from_file,mc.cores = t))
stopCluster(cl)
print(colnames(wgscan))

# calculate genome-wide iHS values
wgscan.ihs <- ihh2ihs(wgscan, min_maf = 0.05,freqbin = 0.01)
df_ihs_rehh<-data.frame(V1=paste0(wgscan.ihs$ihs$CHR,"_",wgscan.ihs$ihs$POSITION),V2=wgscan.ihs$ihs$POSITION,V3=wgscan$FREQ_D,V4=wgscan$IHH_D,V5=wgscan$IHH_A,V6=10**(-wgscan.ihs$ihs$LOGPVALUE),V7=wgscan.ihs$ihs$IHS)
gz1 <-gzfile(paste0(strsplit(j,split = "_i")[[1]][1],".ihs.out.100bins.norm"), "w")
write.table(df_ihs_rehh,gz1,sep='\t', quote = FALSE, row.names=FALSE,col.names=FALSE)
close(gz1)

