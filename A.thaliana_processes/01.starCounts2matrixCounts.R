################################################################################
# RESULTS GENE COUNTS WITH STAR. BOTH PASSAGES 1 & 12.
# MJ Olmo-Uceda
# 2022/11/12
################################################################################
# The STAR -geneCounts output is a three column .tab where:
# column 1: gene id
# column 2: unstranded counts
# column 3: stranded-forward counts
# column 4: stranded-reverse counts
##################################################
library(tidyverse)
library(ggplot2)
library(readxl)
library(stringr)
library(DESeq2)

# path to results of star counts
path.reads <- "data/counts/star_counts/"
files <- list.files(path.reads) 

cts <- read.csv(paste0(path.reads, files[1]), 
                   header = F, 
                   skip = 4, 
                   sep = "\t")
sample <- unlist(strsplit(files[1], split = "_tair10_ReadsPerGene.out.tab"))

rownames(cts) <- cts$V1
cts.unstranded <- cts
cts.revstranded <- cts

cts.unstranded$V1 <- NULL
cts.unstranded$V3 <- NULL
cts.unstranded$V4 <- NULL
#colnames(cts.unstranded) <- sample

cts.revstranded$V1 <- NULL
cts.revstranded$V2 <- NULL
cts.revstranded$V3 <- NULL
#colnames(cts.revstranded) <- sample

samples <- c(sample)
for (i in 2:length(files)){
  sample <- unlist(strsplit(files[i], split = "_tair10_ReadsPerGene.out.tab"))
  #sample <- paste0(unlist(strsplit(sample, "_")), collapse = ".")
  samples[i] <- sample
  c <- read.csv(paste0(path.reads, files[i]), 
                header = F, 
                skip = 4, 
                sep = "\t")
  cts.unstranded <- cbind(cts.unstranded, c$V2)
  cts.revstranded <- cbind(cts.revstranded, c$V4)
}
colnames(cts.unstranded) <- samples
colnames(cts.revstranded) <- samples

# check order
all.equal(colnames(cts.unstranded), infoSamples_all$Sample) # no
all.equal(colnames(cts.revstranded), infoSamples_all$Sample) # no
# sort columns by vector
cts.unstranded <- cts.unstranded[infoSamples_all$Sample]
cts.revstranded <- cts.revstranded[infoSamples_all$Sample]
all.equal(colnames(cts.unstranded), infoSamples_all$Sample) # yes
all.equal(colnames(cts.revstranded), infoSamples_all$Sample) # yes

dim(cts.unstranded) # 32833    77
dim(cts.revstranded) # 32833    77

# write matrix cts of OrVProgression
write.csv(cts.unstranded, 
          file="data/rawCounts_wStarUnstranded_evolEpiMut221114.csv",
          quote = F,
          row.names = T)
write.csv(cts.revstranded, 
          file="data/rawCounts_wStarREVstranded_evolEpiMut221114.csv",
          quote = F,
          row.names = T)
