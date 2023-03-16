library(pairwiseAdonis)

myCol <- c(NAFLD = "#2c7bb6", Cirrhosis = "#fdae61", HCC = "#d7191c")
set.seed(24118)

#––––––––––––––––––––––––––#
# Beta-Diversity analysis  #
#––––––––––––––––––––––––––#
library(data.table)
library(egg)
library(vegan)
library(ggtext)
suppressMessages(library(DESeq2))
#source("analysis/v1/scripts/helper_gmpr.R") # I think I won't use this

# DNA samples
spl_meta <- fread("data/meta/samples_16S_DNA.csv")
spl_meta <- spl_meta[include == TRUE]
spl_meta <- spl_meta[sample != 576]
asv_tab <- read.table("data/dada/asv_tab.tsv", check.names = F)

# RNA samples
# asv_tab_rna <- read.table("../HCC_liver16S/data/dada/asv_tab.tsv")
# colnames(asv_tab_rna) <- gsub("HCC\\.","", colnames(asv_tab_rna))
# colnames(asv_tab_rna) <- gsub("\\.","-", colnames(asv_tab_rna))
# rna_si <- fread("data/rnaseq/sample_information_edit.csv")
# rna_si <- rna_si[external_name %in% colnames(asv_tab_rna)]


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Normalizing for sampling depth using deseq2 #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
cts_norm <- list()
cts_orig <- list()

# faeces
cts_tmp <- as.matrix(asv_tab[,spl_meta[Sample_type == "faeces", sample]])
cts_orig[["faeces (DNA)"]] <- cts_tmp[!apply(cts_tmp,1,function(x) sum(x > 0) < 3),]
rm(cts_tmp)
dds <- DESeqDataSetFromMatrix(countData = cts_orig[["faeces (DNA)"]],
                              colData = spl_meta[Sample_type == "faeces"],
                              design= ~ Condition)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- varianceStabilizingTransformation(dds)
cts_norm[["faeces (DNA)"]] <- assay(dds)

# blood
cts_tmp <- as.matrix(asv_tab[,spl_meta[Sample_type == "blood", sample]])
cts_orig[["blood (DNA)"]] <- cts_tmp[!apply(cts_tmp,1,function(x) sum(x > 0) < 3),]
rm(cts_tmp)
dds <- DESeqDataSetFromMatrix(countData = cts_orig[["blood (DNA)"]],
                              colData = spl_meta[Sample_type == "blood"],
                              design= ~ Condition)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- varianceStabilizingTransformation(dds)
cts_norm[["blood (DNA)"]] <- assay(dds)

# liver (DNA)
cts_tmp <- as.matrix(asv_tab[,spl_meta[Sample_type == "liver", sample]])
cts_orig[["liver (DNA)"]] <- cts_tmp[!apply(cts_tmp,1,function(x) sum(x > 0) < 3),]
rm(cts_tmp)
dds <- DESeqDataSetFromMatrix(countData = cts_orig[["liver (DNA)"]],
                              colData = spl_meta[Sample_type == "liver"],
                              design= ~ Condition)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- varianceStabilizingTransformation(dds)
cts_norm[["liver (DNA)"]] <- assay(dds)

# # liver (RNA)
# dds <- DESeqDataSetFromMatrix(countData = as.matrix(asv_tab_rna[,rna_si[, external_name]]),
#                               colData = rna_si,
#                               design= ~ Condition)
# dds <- estimateSizeFactors(dds, type = "poscounts")
# dds <- varianceStabilizingTransformation(dds)
# cts_norm[["liver (RNA)"]] <- assay(dds)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Claculate Bray-Curtis Dissimilarity and MDS axes  #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
BC_MDS <- lapply(cts_orig , FUN = function(mat) {
  dist.mat <- vegan::vegdist(t(mat), method = "bray")
  MDS <- cmdscale(dist.mat, k=2)
  
  nMDS_single <- metaMDS(t(mat), distance = "bray", k = 1, trymax = 350)
  nMDS <- metaMDS(t(mat), distance = "bray", k = 2, trymax = 350)
  
  dt_BCdist <- data.table(sample = rownames(MDS),
                          MDS1 = MDS[,1],
                          MDS2 = MDS[,2])
  
  dt_BCdist <- data.table(sample = rownames(nMDS$points),
                          MDS1 = nMDS$points[,1],
                          MDS2 = nMDS$points[,2])
  
  
  return(list(MDS = dt_BCdist,
              dist = dist.mat,
              expl_var = c(nMDS_single$stress, nMDS$stress)))
})

# add meta info
BC_MDS$`faeces (DNA)`$MDS <- merge(BC_MDS$`faeces (DNA)`$MDS, spl_meta)
BC_MDS$`blood (DNA)`$MDS <- merge(BC_MDS$`blood (DNA)`$MDS, spl_meta)
BC_MDS$`liver (DNA)`$MDS <- merge(BC_MDS$`liver (DNA)`$MDS, spl_meta)
# BC_MDS$`liver (RNA)`$MDS <- merge(BC_MDS$`liver (RNA)`$MDS, rna_si, by.x = "sample",
#                               by.y = "external_name")

names(BC_MDS) <- gsub(" \\(DNA\\)$","",names(BC_MDS))

# ~ ~ ~ ~ ~ #
# Plotting  #
# ~ ~ ~ ~ ~ #
beta_plots <- list()
posthoc_adonis <- list()
for(itype in names(BC_MDS)) {
  tmp <- copy(BC_MDS[[itype]]$MDS)
  tmp$Condition <- factor(tmp$Condition, levels = c("NAFLD","Cirrhosis","HCC"))
  tmp[["itype"]] <- itype
  
  # perform PERMANOVA
  p_test <- adonis(BC_MDS[[itype]]$dist ~ tmp$Condition,
                   method = "bray")
  p_test <- paste0("<br>(p=",round(p_test$aov.tab$`Pr(>F)`[1], digits = 3), ", R<sup>2</sup> =",
                   round(p_test$aov.tab$R2[1], digits = 3),")")
  tmp[["itype"]] <- paste0("<b>",tmp[["itype"]],"</b>",p_test)
  
  # if(itype == "liver")
  #   tmp <- tmp[MDS1 > -50]
  # Post-Hoc
  posthoc_adonis[[itype]] <- pairwise.adonis(BC_MDS[[itype]]$dist, factors = BC_MDS[[itype]]$MDS$Condition)
  
  p <- ggplot(tmp, aes(MDS1, MDS2, col = Condition, fill = Condition)) +
    stat_ellipse(geom = "polygon", alpha = 0.2) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = myCol) + 
    scale_fill_manual(values = myCol) + 
    labs(x = "NMDS 1", y = "NMDS 2") +
    facet_grid(. ~ itype) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(color = "black"),
          strip.background = element_blank(),
          strip.text = element_markdown(color = "black"),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(color = "black"))
  
  if(itype == "blood") {
    p <- p + theme(legend.position = "bottom")
  } else {
    p <- p + theme(legend.position = "none")
  }
  beta_plots[[itype]] <- p
}

p_beta_comb <- egg::ggarrange(plots = beta_plots, nrow = 1, draw = FALSE)
ggsave("output/plots/16S_betaDiversity.pdf", plot = p_beta_comb, 
       width = 7, height = 3.2)

p_comb_diversity <- ggpubr::ggarrange(p_comb_alpha, p_beta_comb, ncol = 1, labels = c("A","B"))
ggsave("output/plots/submission1/Fig1.pdf", plot = p_comb_diversity,
       height = 6.5, width = 7) 

fwrite(rbindlist(posthoc_adonis, idcol = "Sample_type"),
       "output/files/betaDiv_posthoc_permanova.csv")
