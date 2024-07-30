title: "compartment_analysis"
author: "Mingkee Achom"
date: "2024-04-17"


## load libraries
setwd("/path/to/cool_files_nf")

library(HiCool)
library(ggplot2)                                                                               
library(GenomicRanges)
library(InteractionSet)
library(HiCExperiment)
library(HiContactsData)
library(HiContacts)
library(rtracklayer)
library(BiocParallel, options(warn = 2))    
library(BSgenome.Hsapiens.UCSC.hg38)
library(vcfR)
library(tidyr)
library(patchwork)
library(WGCNA)
library(impute)
library(HiCool)

## for all sample loops, using the highest resolution 
## this can be zoomed later

gw_sample_list <- list.files("./", pattern = ".mcool")
gw_sample_list <- lapply(gw_sample_list, FUN = import, format = "mcool", resolution = 10240000) ## with highest resolution

## importing loops detected from chromisght run
## loop size used ranges from 100

chromo_loop_list <- list.files("./", pattern = "_50mb_loops.tsv")
chromo_loop_list <- lapply(chromo_loop_list, FUN = readr::read_tsv)
chromo_loop_list <- lapply(chromo_loop_list, FUN = plyinteractions::as_ginteractions, seqnames1 = chrom1, seqnames2 = chrom2)

## plot loops for all samples on chrX with refocus on 

library(gridExtra)

## for chrx first
## check on other chromosomes such as chr2, 5, 7, 9, 11, 12, 14, 17 based on DI analyses

plot_list <- list()
sample_list <- list()

## with score and q-value cut-off
for (i in 1:length(chromo_loop_list)){
  sample <- zoom(gw_sample_list[[i]], 40000) |> refocus('chr2')
  sample_list[[i]] <- (sample)
  loops <- chromo_loop_list[[i]][chromo_loop_list[[i]]$score >= 0.3 & chromo_loop_list[[i]]$qvalue <= 1e-6]
  plot_list[[i]] <- plotMatrix(sample_list[[i]], loops = loops, limits = c(-4, -1.2))
}
p <- grid.arrange(grobs=plot_list, nrow=2)

#ggsave(p, filename = "all_samples_50mb_loop_chr2.pdf", device = "pdf", height = 8.5, width = 12, units = "in")

#######################################
## for all chromosomes and sample-wise to visualize compartment tracks
## compartment tracks already exported before so this is generating plots for just visualization

plot_list <- list()
sample_list <- list()

## with score and q-value cut-off
for (i in 1:length(chromo_loop_list)){
  sample <- zoom(gw_sample_list[[i]], 40000) |> refocus('chr2')
  sample_list[[i]] <- (sample)
  loops <- chromo_loop_list[[i]][chromo_loop_list[[i]]$score >= 0.3 & chromo_loop_list[[i]]$qvalue <= 1e-6]
  plot_list[[i]] <- plotMatrix(sample_list[[i]], loops = loops, limits = c(-4, -1.2))
}
p <- grid.arrange(grobs=plot_list, nrow=2)

ggsave(p, filename = "all_samples_50mb_loop_chr2.pdf", device = "pdf",
       height = 8.5, width = 12, units = "in")

###########################################

## now analyze compartments for each of the samples for each chromosomes

#samples <- c("CT20_MicroC_1_S1", "CT20_MicroC_2_S2", "CT20_MicroC_3_S3", "CT20_MicroC_4_S4", "CT20_MicroC_5_S5", "CT20_MicroC_6_S6", "CT20_MicroC_7_S7", "CT20_MicroC_8_S8")
samples <- c('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8')

#test_list <- lapply(all_sample_list, FUN = getCompartments, genome = BSgenome.Hsapiens.UCSC.hg38, chromosomes = chromosomes[j])

compts_list <- vector('list', length(all_sample_list))
cov_list <- list()
for (j in 1:length(all_sample_list)){
  for (k in 1:length(samples)){
    compts <- getCompartments(all_sample_list[[j]][[k]], genome = BSgenome.Hsapiens.UCSC.hg38, chromosomes = chromosomes[j])
    cov <- coverage(metadata(compts)$eigens, weight = 'eigen') 
    cov <- cov[-c(25:195)] ## removing unwanted chromosomes
    #cov_list[[j]] <- cov
    export.bw(cov, paste("./result_files/", "eigen", "_", fileName(compts), "_", focus(compts), ".bw"))
    tf <- topologicalFeatures(compts, "compartments") 
    export.gff3(tf, paste("./result_files/", "compartments", "_", fileName(compts), "_", focus(compts), ".gff3"))
    compts_list[[j]][[k]] <- list(compts)
  }
}

#######################

## computing autocorrelatd contact maps
## first test on individual chromosomes with sample_list that contains hic exp of individual chromosomes only

sample_list <- list.files("./", pattern = ".mcool")

## get compartments only for those chromosomes with more DIs and chrX (containing AR regions)
## based on differential analysis, get compts for the following chromosomes

desired_chr <- c("chr2", "chr5", "chr7", "chr9", "chr11", "chr12", "chr14", "chr17", "chrX")

all_sample_list <- list()
for (i in 1:length(desired_chr)){
  all_sample <- lapply(sample_list, FUN = import, format = "mcool", resolution = 40000, focus= chromosomes[i])
  all_sample_list[[i]] <- all_sample
}

samples <- c('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8')

compts_list <- vector('list', length(all_sample_list))
cov_list <- list()
for (j in 1:length(all_sample_list)){
  for (k in 1:length(samples)){
    compts <- getCompartments(all_sample_list[[j]][[k]], genome = BSgenome.Hsapiens.UCSC.hg38, chromosomes = chromosomes[j])
    cov <- coverage(metadata(compts)$eigens, weight = 'eigen') 
    cov <- cov[-c(25:195)] ## removing unwanted chromosomes
    compts_list[[j]][[k]] <- list(compts)
  }
}


## now after calculating A/B compartments on desired chromosome for all samples
## plotting heatmap with phased coverage
  
plot_list_p1 <- list()
plot_list_p2 <- list()
sample_list <- list()
library(gridExtra)
for (i in 1:length(gw_sample_list)){
  #autocorrelated <- zoom(all_sample_list[[i]], 40000) |> refocus('chr1')
  autocorrelated <- autocorrelate(autocorrelated)
  autocorrelated_list[[i]] <- (autocorrelated)
   plot_list_p1[[i]] <- plotMatrix(autocorrelated_list[[i]], use.scores = 'autocorrelated', limits = c(-1, 1), scale = 'linear', cmap = bwrColors())
   eigen <- coverage(metadata(microC_compts)$eigens, weight = 'eigen')[[1]]
   eigen_df <- tibble(pos = cumsum(runLength(eigen)), eigen = runValue(eigen))
   plot_list_p2[[i]] <- ggplot(eigen_df) + geom_area(aes(x= pos, y = ifelse(eigen < 0, eigen, 0)), fill= "dark blue") + geom_area(aes(x= pos, y = ifelse(eigen > 0, eigen, 0)), fill= "dark red") + theme_void() + coord_cartesian(expand = FALSE) + labs(x = "Genomic position", y = "Eigenvector value") + theme(legend.position="none")
   
}
grid.arrange(grobs=plot_list, nrow=2)


plot_list <- list()
compts_list <- list()
## check on compts for chrX
for (i in 1:length(sample_list)){
  compts <- getCompartments(sample_list[[i]], genome = phasing_track, chromosomes = "chrX")
  compts_list[[i]] <- (compts)
  plot_list[[i]] <- plotMatrix(compts_list[[i]], use.scores = 'autocorrelated', limits = c(-1, 1), scale = 'linear')
}
