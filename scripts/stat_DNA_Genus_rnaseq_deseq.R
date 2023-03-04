# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Correlate Genera abundancies with gene            #
# expression                                        #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
library(data.table)
library(stringr)
library(ggplot2)
library(ggtext)
library(egg)
source("scripts/helper_gmpr.R")

pval_deg_cutoff <- 0.05

# get gene names
gene_names <- rsi[!duplicated(paste0(`Gene ID`,"$",`Gene Name`)),
                  .(GeneID = `Gene ID`, `Gene Name`)]
# gene_names <- fread("/mnt/nuuk/2020/HCC/RNA_SEQ_201215/stringtieFPKM/J22798-L1_S22_L002_R1_001Aligned.sortedByCoord.out.gene_abund.txt")
# gene_names <- gene_names[, .(GeneID = `Gene ID`, `Gene Name`)]
setkey(gene_names, "GeneID")

# get DEGS (based on Condition contrasts) (to extract directions of DEGs)
deg_stats <- readRDS("output/files/deg_stats.RDS")
#deg_stats <- deg_stats[GeneID != "ENSG00000273730"]
#deg_stats <- deg_stats[GeneID != "ENSG00000275757"]
#deg_stats <- deg_stats[GeneID != "ENSG00000277739"]
hcc_up <- deg_stats[pt.col %in% c("HCCtop","HCC") & contrast == "NAFLD_vs_HCC", GeneID]
hcc_dw <- deg_stats[pt.col %in% c("NAFLDtop","NAFLD") & contrast == "NAFLD_vs_HCC", GeneID]

# get Liver ASV counts
asv_tab <- read.table("data/dada/asv_tab.tsv", check.names = F)
spl_meta <- fread("data/meta/samples_16S_DNA.csv")
spl_meta <- spl_meta[include == TRUE]
spl_meta <- spl_meta[sample != 576]
spl_meta <- spl_meta[sample %in% colnames(asv_tab)]

# load taxonomic predictions
asv_tax <- read.table("data/dada/asv_tax.tsv", sep = "\t")
tmp <- rownames(asv_tax)
asv_tax <- data.table(asv_tax)
asv_tax$cluster <- tmp
asv_tax[Genus == "Unknown", Genus := paste(Genus, Family)]
asv_tax[, cluster_plot := paste0("(",cluster,")<br><i><b>",Genus,"</b></i>")]

# get gene expression counts
# tmp <- fread(paste0(rnas_data_path,"featureCounts/merged_gene_counts.txt.gz"))
tmp <- fread("data/rnaseq/salmon.merged.gene_counts.tsv")
tmp <- tmp[!grepl("_PAR_Y$", gene_id)]
tmp[, gene_id := gsub("\\..*$","",gene_id)]
expr_genes <- readLines("output/files/expressed_genes.txt")
tmp <- tmp[gene_id %in% expr_genes]
cts <- as.matrix(tmp[,-(1:2)])
rownames(cts) <- tmp$gene_id
colnames(cts) <- str_sub(colnames(cts), 1, 6)
#cts <- cts[!(rownames(cts) %in% c("ENSG00000277739","ENSG00000273730","ENSG00000275757")),] # pseudogene
cts <- round(cts)

# sample infor RNA-seq
si_ddeq <- copy(rseq_si)
si_ddeq[, name := gsub("-S1$","",name)]
si_ddeq <- si_ddeq[name %in% colnames(cts)]
si_ddeq[, age_sc := scale(age)]
si_ddeq[, Condition := factor(Condition)]
si_ddeq[, sex := factor(sex)]
si_ddeq[, genus_abun := NA_real_]
si_ddeq <- si_ddeq[sample_16S_data %in% colnames(tab_genus_counts)]
si_ddeq <- si_ddeq[Condition == "HCC"] # comment #2 by Reviewe #3 (Simpson's paradox)

# Normalizing Genus counts
liver_counts_norm <- tab_genus_counts[, spl_meta[Sample_type == "liver", sample]]
liver_counts_norm <- liver_counts_norm[, si_ddeq$sample_16S_data]
liver_counts_norm <- liver_counts_norm[apply(liver_counts_norm,1,function(x) sum(x > 0) > 0),]
#liver_counts_norm <- liver_counts_norm[rownames(liver_counts_norm) != "Unknown Unknown",]
dds_liver <- DESeqDataSetFromMatrix(countData = liver_counts_norm[, si_ddeq$sample_16S_data],
                                    colData = si_ddeq,
                                    design= ~ 1)
dds_liver <- estimateSizeFactors(dds_liver, type = "poscounts")
dds_liver <- estimateDispersions(dds_liver)
# liver_counts_norm <- assay(varianceStabilizingTransformation(dds_liver, blind = F))
# liver_counts_norm <- counts(dds_liver, normalize = TRUE)
liver_counts_norm <- counts(dds_liver, normalize = FALSE)
liver_counts_norm <- t(t(liver_counts_norm)/colSums(liver_counts_norm))

genus_keep <- names(which(apply(liver_counts_norm,1,
                                function(x) sum(x > 0.001) >= 8)))
genus_keep <- genus_keep[genus_keep != "Unknown Unknown"]


stat_transcr_micro <- list()
for(genusi in genus_keep) {
  cat(genusi,"\n")
  
  si_ddeq2 <- copy(si_ddeq)
  si_ddeq2$genus_abun <- scale(liver_counts_norm[genusi,])
  
  cts2 <- cts[,si_ddeq2$name]
  
  design <- ~ genus_abun + age_sc + sex
  
  dds <- DESeqDataSetFromMatrix(countData = cts2,
                                colData   = si_ddeq2,
                                design    = design)
  #dds <- nbinomWaldTest(dds, maxit = 500)
  dds <- DESeq(dds)
  
  res_abun <- results(dds, name = "genus_abun")
  tmp_names <- row.names(res_abun)
  res_abun <- as.data.table(res_abun)
  res_abun$GeneID <- tmp_names
  res_abun$genus <- genusi
  
  stat_transcr_micro[[genusi]] <- res_abun
}

stat_TM <- rbindlist(stat_transcr_micro)
stat_TM[!is.na(pvalue), padj := p.adjust(pvalue), by = GeneID]

# Heatmap
pval_cutoff <- 0.00001
rel_genus <- stat_TM[padj < pval_cutoff, unique(genus)]
rel_genes <- stat_TM[padj < pval_cutoff, unique(GeneID)]

stat_TM_small <- copy(stat_TM[genus %in% rel_genus & GeneID %in% rel_genes])

stat_TM_small[, association_dir := NA_character_]
stat_TM_small[padj < pval_cutoff & stat < 0, association_dir := "negative"]
stat_TM_small[padj < pval_cutoff & stat > 0, association_dir := "positive"]

stat_TM_small[, dir := "no diff. expr."]
stat_TM_small[GeneID %in% hcc_up, dir := "HCC upreg."]
stat_TM_small[GeneID %in% hcc_dw, dir := "NAFLD upreg."]

stat_TM_small[, GeneID2 := GeneID]
#stat_TM_small[dir != "no diff. expr.", GeneName := gene_names[GeneID, `Gene Name`]]
stat_TM_small[, GeneName := gene_names[GeneID, `Gene Name`]]
stat_TM_small[dir == "HCC upreg.", GeneID2 := paste0("<i style='color:#d7191c'>",GeneName,"</i>")]
stat_TM_small[dir == "NAFLD upreg.", GeneID2 := paste0("<i style='color:#2c7bb6'>",GeneName,"</i>")]
stat_TM_small[dir == "no diff. expr.", GeneID2 := paste0("<i style='color:#000000'>",GeneName,"</i>")]

labeled_genes <- stat_TM_small[dir != "no diff. expr.", unique(GeneID2)]

# stat_TM_small[, tmp := ifelse(is.na(association_dir),0,
#                               ifelse(association_dir == "negative",-1,1))]

stat_TM_small[, stat2 := ifelse(padj < pval_cutoff, stat, NA)]

genus_order <- cluster_factors(stat_TM_small, by.x = "genus", by.y = "GeneID",
                               value = "stat")
gene_order <- cluster_factors(stat_TM_small, by.x = "GeneID", by.y = "genus",
                              value = "stat")

dictio <- stat_TM_small[GeneID %in% gene_order,][!duplicated(GeneID)][, .(GeneID, GeneID2)]
setkey(dictio, "GeneID")

stat_TM_small$GeneID <- factor(stat_TM_small$GeneID, levels = gene_order, labels = dictio[gene_order, GeneID2])
stat_TM_small$genus <- factor(stat_TM_small$genus, levels = genus_order)

#lt_dt <- data.table(GeneID2 = labeled_genes)



p_TM <- ggplot(stat_TM_small[!is.na(stat2)],
               aes(GeneID, genus, fill = stat2,
                   col = stat2)) +
  #geom_vline(aes(xintercept = GeneID2), lt_dt) +
  geom_tile() +
  scale_fill_gradient2(high = "#b2182b", low = "#2166ac", mid = "#f7f7f7", na.value = NA) +
  scale_color_gradient2(high = "#b2182b", low = "#2166ac", mid = "#f7f7f7", na.value = NA) +
  # scale_x_discrete(labels = ifelse(gene_order %in% labeled_genes,
  #                                  gene_order, "")) +
  #scale_linetype_manual(values = c(NA,1,1)) +
  #facet_grid(. ~ dir, scales = "free", space = "free") +
  labs(x = "Gene", y = "Genus",
       fill = "Wald statistic",
       col = "Wald statistic") +
  theme_bw() +
  theme(axis.text.y = element_text(color = "black", face = "italic"),
        axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(color = "black"),
        #panel.grid.major.x = element_blank(),
        legend.position = "bottom")
p_TM

# p_TM <- ggplot(stat_TM_small,
#                aes(GeneID2, genus, fill = association_dir,
#                    col = association_dir)) +
#   #geom_vline(aes(xintercept = GeneID2), lt_dt) +
#   geom_tile() +
#   scale_fill_manual(values = c("#5e3c99","#e66101"), na.translate = F) +
#   scale_color_manual(values = c("#5e3c99","#e66101"), na.translate = F) +
#   # scale_x_discrete(labels = ifelse(gene_order %in% labeled_genes,
#   #                                  gene_order, "")) +
#   #scale_linetype_manual(values = c(NA,1,1)) +
#   #facet_grid(. ~ dir, scales = "free", space = "free") +
#   labs(x = "Gene", y = "Genus",
#        fill = "Association direction of\ngene expression\nand Genus abundance",
#        col = "Association direction of\ngene expression\nand Genus abundance") +
#   theme_bw() +
#   theme(axis.text.y = element_text(color = "black", face = "italic"),
#         axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust = 1),
#         axis.title.y = element_text(color = "black"),
#         #panel.grid.major.x = element_blank(),
#         legend.position = "right")
# p_TM

ggsave("output/plots/rna_16S_Genus-x-Expression-Associations.pdf",
       plot = p_TM, width = 16, height = 4.5)
