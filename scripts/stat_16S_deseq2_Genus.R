# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Differnetial Genus abundance analysis #
# via DESeq2                            #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
library(ggtext)
library(ggrepel)

comparison_names <- c("<i style='color:#fd8b1b'>Cirrhosis</i> vs. <i style='color:#d7191c'>HCC</i>",
                      "<i style='color:#2c7bb6'>NAFLD</i> vs. <i style='color:#fd8b1b'>Cirrhosis</i>",
                      "<i style='color:#2c7bb6'>NAFLD</i> vs. <i style='color:#d7191c'>HCC</i>")

myCol <- c(NAFLD = "#2c7bb6", Cirrhosis = "#fdae61", HCC = "#d7191c")

asv_tab <- read.table("data/dada/asv_tab.tsv", check.names = F)
spl_meta <- fread("data/meta/samples_16S_DNA.csv")
spl_meta <- spl_meta[include == TRUE]
spl_meta <- spl_meta[sample %in% colnames(asv_tab)]

#genus_tab <- read.table("data/dada/tax_tab_Genus.tsv", check.names = F, sep = "\t")
#family_tab <- read.table("data/dada/tax_tab_Family.tsv", check.names = F, sep = "\t")

# load taxonomic predictions
asv_tax <- read.table("data/dada/asv_tax.tsv", sep = "\t")
tmp <- rownames(asv_tax)
asv_tax <- data.table(asv_tax)
asv_tax$cluster <- tmp
asv_tax[Genus == "Unknown", Genus := paste(Genus, Family)]
asv_tax[, cluster_plot := paste0("(",cluster,")<br><i><b>",Genus,"</b></i>")]


#asv_tab_rel <- t(t(asv_tab)/colSums(asv_tab))
#asv_tab_rel_faecal <- asv_tab_rel[fecal_asvs,]

DT_asv_counts <- data.table(as.table(as.matrix(asv_tab)))
setnames(DT_asv_counts, c("cluster","sample","N"))
DT_asv_counts <- merge(DT_asv_counts, asv_tax, by = "cluster")
DT_genus_counts <- DT_asv_counts[, .(N = sum(N)), by = .(sample, Genus)]
tab_genus_counts <- dcast(DT_genus_counts, Genus ~ sample, value.var = "N")
tmp <- tab_genus_counts$Genus
tab_genus_counts <- as.matrix(tab_genus_counts[,-1])
rownames(tab_genus_counts) <- tmp
rel_genera <- names(which(apply(tab_genus_counts, 1, function(x) sum(x > 0) > 0)))
tab_genus_counts <- tab_genus_counts[rel_genera,]

# DESeq - faeces - Genus level counts
dds_faeces <- DESeqDataSetFromMatrix(countData = tab_genus_counts[, spl_meta[Sample_type == "faeces", sample]],
                                     colData = spl_meta[Sample_type == "faeces"],
                                     design= ~ Condition)
dds_faeces <- estimateSizeFactors(dds_faeces, type = "poscounts")
dds_faeces <- estimateDispersions(dds_faeces)
dds_faeces <- nbinomWaldTest(dds_faeces)

# DESeq - blood - Genus level counts
dds_blood <- DESeqDataSetFromMatrix(countData = tab_genus_counts[, spl_meta[Sample_type == "blood", sample]],
                                    colData = spl_meta[Sample_type == "blood"],
                                    design= ~ Condition)
dds_blood <- estimateSizeFactors(dds_blood, type = "poscounts")
dds_blood <- estimateDispersions(dds_blood)
dds_blood <- nbinomWaldTest(dds_blood)

# DESeq - liver - Genus level counts
dds_liver <- DESeqDataSetFromMatrix(countData = tab_genus_counts[, spl_meta[Sample_type == "liver", sample]],
                                    colData = spl_meta[Sample_type == "liver"],
                                    design= ~ Condition)
dds_liver <- estimateSizeFactors(dds_liver, type = "poscounts")
dds_liver <- estimateDispersions(dds_liver)
dds_liver <- nbinomWaldTest(dds_liver)

deseq_genus_cts <- list(faeces = dds_faeces, 
                        blood = dds_blood,
                        liver = dds_liver)
rm(dds_faeces)
rm(dds_blood)
rm(dds_liver)


stats_genus <- list()
k <- 1
for(itype in names(deseq_genus_cts)) {
  # N vs C
  zut <- results(deseq_genus_cts[[itype]], contrast=c("Condition","NAFLD","Cirrhosis"))
  tmp <- rownames(zut)
  zut <- data.table(data.frame(zut))
  zut$Genus <- tmp
  zut$Sample_type <- itype
  zut$grp1 <- "NAFLD"
  zut$grp2 <- "Cirrhosis"
  zut$contrast <- "<i style='color:#2c7bb6'>NAFLD</i> vs. <i style='color:#fd8b1b'>Cirrhosis</i>"
  
  stats_genus[[k]] <- copy(zut)
  k <- k + 1
  
  # N vs H
  zut <- results(deseq_genus_cts[[itype]], contrast=c("Condition","NAFLD","HCC"))
  tmp <- rownames(zut)
  zut <- data.table(data.frame(zut))
  zut$Genus <- tmp
  zut$Sample_type <- itype
  zut$grp1 <- "NAFLD"
  zut$grp2 <- "HCC"
  zut$contrast <- "<i style='color:#2c7bb6'>NAFLD</i> vs. <i style='color:#d7191c'>HCC</i>"
  stats_genus[[k]] <- copy(zut)
  k <- k + 1
  
  # N vs C
  zut <- results(deseq_genus_cts[[itype]], contrast=c("Condition","Cirrhosis","HCC"))
  tmp <- rownames(zut)
  zut <- data.table(data.frame(zut))
  zut$Genus <- tmp
  zut$Sample_type <- itype
  zut$grp1 <- "Cirrhosis"
  zut$grp2 <- "HCC"
  zut$contrast <- "<i style='color:#fd8b1b'>Cirrhosis</i> vs. <i style='color:#d7191c'>HCC</i>"
  
  stats_genus[[k]] <- copy(zut)
  k <- k + 1
  
  
  
}
stats_genus <- rbindlist(stats_genus)
stats_genus$contrast <- factor(stats_genus$contrast, levels = comparison_names)
stats_genus$Sample_type <- factor(stats_genus$Sample_type, levels = c("faeces",
                                                                      "blood",
                                                                      "liver"))
stats_genus <- stats_genus[!is.na(baseMean) & !is.na(log2FoldChange)]
stats_genus[, color := NA_character_]
stats_genus[padj < 0.05 & stat > 0, color := grp1]
stats_genus[padj < 0.05 & stat < 0, color := grp2]

stats_genus[, label := NA_character_]
stats_genus[padj < 0.05, label := Genus]


ggplot(stats_genus, aes(log10(baseMean), log2FoldChange, label = label,
                        col = color))+
  geom_point(size = 0.5) + 
  geom_text_repel(aes(angle = ifelse(log2FoldChange > 0, 25, -25)),
                  size = 2.25, color = "black") +
  theme_bw() +
  scale_color_manual(values = myCol, drop = F, na.value = "grey80") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        #panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text = element_markdown(face = "bold", colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none") + 
  facet_grid(`Sample_type`~`contrast`)

# heatmap solution
rel_genera <- stats_genus[padj < 0.05, unique(Genus)]
stats_genus[, pval.labs := pval_labs(padj)]
stats_genus[is.na(pval.labs), pval.labs := ""]
stats_genus[pval.labs == "****", pval.labs := "***"]

p_diff_genus <- ggplot(stats_genus[!(grp1 == "Cirrhosis" & grp2 == "HCC") & Genus %in% rel_genera],
       aes(grp2, Genus, fill = -log2FoldChange)) +
  geom_tile() +
  geom_text(aes(label = pval.labs)) + 
  scale_fill_gradient2(low = "#5e3c99", mid = "#ffffff", high = "#e66101") +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  facet_grid(.~Sample_type) +
  theme_bw() +
  labs(fill = expression(log[2]~FC), x = "Condition\n(compared to NAFLD)") +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(face = "italic", color = "black"),
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5,
                                   hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave(p_diff_genus, file = "output/plots/16S_diffGenus_deseq2.pdf",
       width = 4.8, height = 4.6)
ggsave(p_diff_genus, file = "output/plots/submission1/Fig3.pdf",
       width = 4.8, height = 4.6)
