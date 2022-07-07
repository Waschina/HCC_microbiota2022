library(data.table)
library(stringr)
library(ggplot2)
library(ggtext)
library(egg)
myCol <- c(NAFLD = "#2c7bb6", HCC = "#d7191c", Cirrhosis = "#fdae61")

#
# Collecting FPKM values for PCA
#
rseq_si <- fread("data/rnaseq/sample_information_edit.csv")
rseq_si[, name := gsub("-S1$","",name)]
rseq_si$Condition <- factor(rseq_si$Condition, levels = c("NAFLD", "Cirrhosis", "HCC"))

rsi <- list()
rnas_data_path <- "data/rnaseq/"
fpkm_files <- dir(paste0(rnas_data_path,"stringtieFPKM/"))
fpkm_files <- fpkm_files[grepl("gene_abund\\.txt\\.gz$", fpkm_files)]
for(i in fpkm_files)
  rsi[[i]] <- fread(paste0(rnas_data_path,"stringtieFPKM/",i))

names(rsi) <- str_sub(names(rsi), start = 1, end = 6)
rsi <- rbindlist(rsi, idcol = "name")
#rsi <- rsi[name != "J22812"] # bad in quality control

rsi <- merge(rsi, rseq_si, by = "name")
rel_genes <- rsi[TPM > 0.5, .N, by = "Gene ID"][N >= 3, `Gene ID`]
rsi <- rsi[`Gene ID` %in% rel_genes]
saveRDS(rsi, "output/files/TPM_FPKM.RDS")
writeLines(rel_genes, "output/files/expressed_genes.txt")


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# PCA on the basis of FPKM values   #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#T_S <- dcast(rsi, `Gene ID` + `Gene Name` ~ external_name, value.var="FPKM", fill = 0, fun.aggregate = sum)
T_S <- dcast(rsi, `Gene ID` + `Gene Name` ~ external_name, value.var="TPM", fill = 0, fun.aggregate = sum)
T_S <- T_S[apply(T_S[,-(1:2)],1,FUN = function(x) sum(x==0)) < 20]
saveRDS(T_S, "output/files/TPM_table.RDS")

#
# PCA
#
rna_pca <- prcomp(t(T_S[,-(1:2)]), scale. = T)
rp_dt      <- data.table(sample = colnames(T_S)[-(1:2)],
                         PC1    = rna_pca$x[,1],
                         PC2    = rna_pca$x[,2],
                         PC3    = rna_pca$x[,3])
#plot(rna_pca$sdev)
pc1_perc <- round(summary(rna_pca)[6]$importance[2,1]*100)
pc2_perc <- round(summary(rna_pca)[6]$importance[2,2]*100)
pc3_perc <- round(summary(rna_pca)[6]$importance[2,3]*100)

rp_dt <- merge(rp_dt, rseq_si, by.x = "sample", by.y = "external_name")

p <- ggplot(rp_dt, aes(PC1, PC2, col = Condition)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = c("#2c7bb6","#fdae61","#d7191c")) +
  scale_fill_manual(values = c("#2c7bb6","#fdae61","#d7191c")) +
  xlab(label = paste0("PC1 (",pc1_perc," %)")) +
  ylab(label = paste0("PC2 (",pc2_perc," %)")) +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none")
p
ggsave("output/plots/rna_PCA.pdf",p, width = 4.5, height = 3.15)
ggsave("output/plots/submission1/FigS2.pdf",p, width = 4.5, height = 3.15)
p_tpm_PCA <- p

# ~ ~ ~ ~ ~ ~ ~ ~ #
# Heatmap on TPM  #
# ~ ~ ~ ~ ~ ~ ~ ~ #
library(ComplexHeatmap)
library(circlize)
mat <- as.matrix(as.matrix(T_S[,3:ncol(T_S)]))
mat <- t(scale(t(mat)))

hc = hclust(dist(mat))
group = cutree(hc, k = 2)

#mat <- mat[apply(mat, 1, function(x) sum(x>10) >= 25),]
# mat <- log10(mat)
# z_inds <- which(mat == -Inf, arr.ind = T)
# mat[z_inds] <- Inf
# mat[z_inds] <- min(mat)
dim(mat)

# heat annotative (columns= samples)
#dt_dis <- rseq_si[name != "J22812"]
rseq_si
dt.fecal <- readRDS("output/files/dt.fecal.RDS")
rseq_si <- merge(rseq_si, dt.fecal[, .(sample, prop.fecal)], all.x = T, by.x = "sample_16S_data", by.y = "sample")
#dt_dis <- rseq_si[match(colnames(mat), external_name)][,.(Condition)]
dt_dis <- rseq_si[match(colnames(mat), external_name)][,.(Condition, `Proportion faecal bacteria` = prop.fecal)]

ha_column = HeatmapAnnotation(df = dt_dis,
                              col = list(Condition = myCol))

# ht1 <- Heatmap(mat, name = "scaled TPM", column_title = "", 
#                top_annotation = ha_column,
#                cluster_rows = cluster_within_group(t(mat), group),
#                row_split = 2, border = TRUE)
ht1 <- Heatmap(mat, name = "scaled TPM", column_title = "", 
               top_annotation = ha_column, border = TRUE, show_row_dend = FALSE)

pdf(file = "output/plots/rna_TPM_heatmap.pdf", width = 7.25, height = 7.5)
ht1
dev.off()


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# DEseq2 - differential expression  #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
library(DESeq2)

# ---- all combined ---- #
tmp <- fread(paste0(rnas_data_path,"featureCounts/merged_gene_counts.txt.gz"))
tmp <- tmp[Geneid %in% rel_genes]
cts <- as.matrix(tmp[,-(1:2)])
rownames(cts) <- tmp$Geneid
colnames(cts) <- str_sub(colnames(cts), 1, 6)

si_ddeq <- copy(rseq_si)
si_ddeq[, name := gsub("-S1$","",name)]
si_ddeq <- si_ddeq[name %in% colnames(cts)]
si_ddeq[, age_sc := scale(age)]
si_ddeq[, Condition := factor(Condition)]
si_ddeq[, sex := factor(sex)]
#si_ddeq <- si_ddeq[external_name != "FLC-LB8"]
cts <- cts[,si_ddeq$name]
design <- ~ Condition + age_sc + sex

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData   = si_ddeq,
                              design    = design)
dds <- DESeq(dds)
#resultsNames(dds)

res_N_H <- results(dds, contrast = c("Condition","NAFLD","HCC"))
tmp_names <- row.names(res_N_H)
res_N_H <- as.data.table(res_N_H)
res_N_H$GeneID <- tmp_names

res_N_C <- results(dds, contrast = c("Condition","NAFLD","Cirrhosis"))
tmp_names <- row.names(res_N_C)
res_N_C <- as.data.table(res_N_C)
res_N_C$GeneID <- tmp_names

res_C_H <- results(dds, contrast = c("Condition","Cirrhosis","HCC"))
tmp_names <- row.names(res_C_H)
res_C_H <- as.data.table(res_C_H)
res_C_H$GeneID <- tmp_names

#resTest <- res

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# MA-plots                      #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
library(ggrepel)
library(ggtext)
library(glue)
n_top <- 5

fc_cutoff <- 1
pv_cutoff <- 0.05

# C vs H
#tmp_C_H <- copy(res_C_H[!is.na(padj) & !is.na(log2FoldChange)])
tmp_C_H <- copy(res_C_H[!is.na(log2FoldChange)])
tmp_C_H$contrast <- "Cirr_vs_HCC" 

#tmp_N_H <- copy(res_N_H[!is.na(padj) & !is.na(log2FoldChange)])
tmp_N_H <- copy(res_N_H[!is.na(log2FoldChange)])
tmp_N_H$contrast <- "NAFLD_vs_HCC" 

#tmp_N_C <- copy(res_N_C[!is.na(padj) & !is.na(log2FoldChange)])
tmp_N_C <- copy(res_N_C[!is.na(log2FoldChange)])
tmp_N_C$contrast <- "NAFLD_vs_Cirr" 

tmp <- rbind(tmp_C_H, tmp_N_H, tmp_N_C)

fwrite(tmp, "output/files/HCC_deg_stats.tsv", quote = F, sep = "\t")

deg_stats <- copy(tmp)
deg_stats[is.na(padj), padj := 1]

gene_names <- copy(rsi[,.(GeneID = `Gene ID`,`Gene Name`)])
gene_names <- gene_names[!duplicated(GeneID)]

deg_stats <- merge(deg_stats, gene_names, by = "GeneID", all.x = T)

deg_stats <- deg_stats[order(contrast, log2FoldChange < 0, padj)]
deg_stats[log2FoldChange < 0 & padj < pv_cutoff, top_i := 1:.N, by = "contrast"]
deg_stats[log2FoldChange > 0 & padj < pv_cutoff, top_i := 1:.N, by = "contrast"]
deg_stats[, label := NA_character_]
deg_stats[top_i <= n_top, label := `Gene Name`]

# assign color factors
deg_stats[, pt.col := "0"]
deg_stats[contrast == "Cirr_vs_HCC" & log2FoldChange > 0 & padj < pv_cutoff, pt.col := "CIRR"]
deg_stats[contrast == "Cirr_vs_HCC" & log2FoldChange < 0 & padj < pv_cutoff, pt.col := "HCC"]
deg_stats[contrast == "NAFLD_vs_HCC" & log2FoldChange > 0 & padj < pv_cutoff, pt.col := "NAFLD"]
deg_stats[contrast == "NAFLD_vs_HCC" & log2FoldChange < 0 & padj < pv_cutoff, pt.col := "HCC"]
deg_stats[contrast == "NAFLD_vs_Cirr" & log2FoldChange > 0 & padj < pv_cutoff, pt.col := "NAFLD"]
deg_stats[contrast == "NAFLD_vs_Cirr" & log2FoldChange < 0 & padj < pv_cutoff, pt.col := "CIRR"]

deg_stats[pt.col == "NAFLD" & !is.na(label), pt.col := "NAFLDtop"]
deg_stats[pt.col == "CIRR" & !is.na(label), pt.col := "CIRRtop"]
deg_stats[pt.col == "HCC" & !is.na(label), pt.col := "HCCtop"]

deg_stats$pt.col <- factor(deg_stats$pt.col, levels = c("0","NAFLD", "NAFLDtop", "CIRR", "CIRRtop", "HCC", "HCCtop")) # ordering

#deg_stats[log2FoldChange < -10, log2FoldChange := -10]
#c("#2c7bb6","#fdae61","#d7191c")
saveRDS(deg_stats, file = "output/files/deg_stats.RDS")

# modify facet labels
comparison_names <- c("<i style='color:#fd8b1b'>Cirrhosis</i> vs. <i style='color:#d7191c'>HCC</i>",
                      "<i style='color:#2c7bb6'>NAFLD</i> vs. <i style='color:#fd8b1b'>Cirrhosis</i>",
                      "<i style='color:#2c7bb6'>NAFLD</i> vs. <i style='color:#d7191c'>HCC</i>")
names(comparison_names) <- c("Cirr_vs_HCC","NAFLD_vs_Cirr","NAFLD_vs_HCC")

p <- ggplot(deg_stats[order(-padj)], 
            aes(log10(baseMean), log2FoldChange, col = pt.col, label = label)) +
  geom_point(size = 0.5) + 
  geom_text_repel(aes(angle = ifelse(log2FoldChange > 0, 25, -25)),
                  size = 3.25, color = "black") +
  theme_bw() +
  scale_color_manual(values = c("lightgrey", # insig
                                "#96c4e6",   # NAFLD
                                "#2c7bb6",   # NAFLDtop
                                "#fecb9b",   # CIRR
                                "#fd8b1b",   # CIRRtopp
                                "#f28a8b",   # HCC
                                "#d7191c"    # HCCtop
  ), drop = F) +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text = element_markdown(face = "bold", colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none") + 
  facet_grid(.~`contrast`, labeller = labeller(contrast = comparison_names)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
p
p_deg <- p

ggsave("output/plots/rna_DEG_MAplot.png",p_deg, width = 14, height = 5)
ggsave("output/plots/rna_DEG_MAplot.pdf",p_deg, width = 14, height = 5)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Venn diagramme of DEG #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
library(VennDiagram)

CH_degs <- deg_stats[contrast == "Cirr_vs_HCC"   & padj < pv_cutoff, GeneID]
NH_degs <- deg_stats[contrast == "NAFLD_vs_HCC"  & padj < pv_cutoff, GeneID]
NC_degs <- deg_stats[contrast == "NAFLD_vs_Cirr" & padj < pv_cutoff, GeneID]

vd_contrast_degs <- venn.diagram(list(CH_degs, NH_degs, NC_degs),
                                 category.names = c("Cirrhosis vs. HCC","NAFLD vs. HCC","NAFLD vs. Cirrhosis"),
                                 filename = NULL, 
                                 ext.text = F,
                                 
                                 # circles
                                 #lty = 'blank',
                                 #fill = myCol[c("HCC","Cirrhosis")],
                                 
                                 # Numbers
                                 cex = 1,
                                 fontfamily = "sans",
                                 
                                 # catagories
                                 cat.default.pos = "outer",
                                 cat.fontfamily = "sans",
                                 cat.pos = c(0, 0, 180))

ggsave(vd_contrast_degs, file="output/plots/rna_Venn_DEG_contrast.svg", device = "svg", width = 3, heigh = 3)

# NAFDL as reference
H_ups <- deg_stats[contrast == "NAFLD_vs_HCC"   & padj < pv_cutoff & log2FoldChange < 0, GeneID]
C_ups <- deg_stats[contrast == "NAFLD_vs_Cirr"  & padj < pv_cutoff & log2FoldChange < 0, GeneID]
H_dws <- deg_stats[contrast == "NAFLD_vs_HCC"   & padj < pv_cutoff & log2FoldChange > 0, GeneID]
C_dws <- deg_stats[contrast == "NAFLD_vs_Cirr"  & padj < pv_cutoff & log2FoldChange > 0, GeneID]

vd_up <- venn.diagram(list(H_ups, C_ups), 
                      category.names = c("HCC upreg.","Cirrhosis upreg."), 
                      filename = NULL, 
                      ext.text = F,
                      
                      # circles
                      lty = 'blank',
                      fill = myCol[c("HCC","Cirrhosis")],
                      
                      # Numbers
                      cex = 1,
                      fontfamily = "sans",
                      
                      # catagories
                      cat.default.pos = "outer",
                      cat.fontfamily = "sans",
                      cat.pos = c(0, 0),
)
ggsave(vd_up, file="output/plots/rna_Venn_upreg.svg", device = "svg", width = 3, heigh = 3)

vd_dw <- venn.diagram(list(H_dws, C_dws), 
                      category.names = c("HCC downreg.","Cirrhosis downreg."), 
                      filename = NULL, 
                      ext.text = F,
                      
                      # circles
                      lty = 'blank',
                      fill = myCol[c("HCC","Cirrhosis")],
                      
                      # Numbers
                      cex = 1,
                      fontfamily = "sans",
                      
                      # catagories
                      cat.default.pos = "outer",
                      cat.fontfamily = "sans",
                      cat.pos = c(0, 0),
)
ggsave(vd_dw, file="output/plots/rna_Venn_downreg.svg", device = "svg", width = 3, heigh = 3)

# - - - - - - - - - - #
# GO-Term Enrichment  #
# - - - - - - - - - - #
library(topGO)
# reflist_genIDs <- c(res_N_H[baseMean > 0 & !is.na(pvalue), GeneID],
#                     res_C_H[baseMean > 0 & !is.na(pvalue), GeneID],
#                     res_N_C[baseMean > 0 & !is.na(pvalue), GeneID])
# reflist_genIDs <- unique(reflist_genIDs)
reflist_genIDs <- rel_genes

get_GOp_one_sided <- function(cond_a, cond_b, dir, go_level = "BP") {

  dt <- results(dds, contrast = c("Condition",cond_a, cond_b), altHypothesis = dir)
  tmp_names <- row.names(dt)
  dt <- as.data.table(dt)
  dt$GeneID <- tmp_names
  
  gene_list <- dt[!is.na(padj) & baseMean > 0 & GeneID %in% rel_genes, padj]
  names(gene_list) <- dt[!is.na(padj) & baseMean > 0 & GeneID %in% rel_genes, GeneID]
  
  
  GetGeneList <- function(input){
    input < 0.05 
  }
  
  
  # GOdata
  GOdata <-new ("topGOdata", 
                ontology = go_level, 
                allGenes = gene_list, 
                geneSel = GetGeneList,
                nodeSize = 10, 
                annot = annFUN.org,
                mapping="org.Hs.eg.db", 
                ID = "ensembl"
  )
  
  # run tests
  result.fisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  #result.fisher@score <- p.adjust(result.fisher@score, method = "BH")
  
  # run GenTable
  allRes <- data.table(GenTable(GOdata,
                                classicFisher = result.fisher,
                                topNodes = length(result.fisher@score)))
  allRes$classicFisher <- result.fisher@score[allRes$GO.ID]
  allRes[, padj := p.adjust(classicFisher, method = "BH")]
  
  allRes[, Condition := cond_b]
  allRes[, Ontology  := go_level]
  #allRes$classicFisher <- as.numeric(allRes$classicFisher)
  
  return(allRes)
}

# NAFLD as reference
goBP_stat_N_H_up <- get_GOp_one_sided("NAFLD", "HCC", "less", "BP")
goMF_stat_N_H_up <- get_GOp_one_sided("NAFLD", "HCC", "less", "MF")
goCC_stat_N_H_up <- get_GOp_one_sided("NAFLD", "HCC", "less", "CC")

goBP_stat_N_H_dw <- get_GOp_one_sided("NAFLD", "HCC", "greater", "BP")
goMF_stat_N_H_dw <- get_GOp_one_sided("NAFLD", "HCC", "greater", "MF")
goCC_stat_N_H_dw <- get_GOp_one_sided("NAFLD", "HCC", "greater", "CC")

goBP_stat_N_C_up <- get_GOp_one_sided("NAFLD", "Cirrhosis", "less", "BP")
goMF_stat_N_C_up <- get_GOp_one_sided("NAFLD", "Cirrhosis", "less", "MF")
goCC_stat_N_C_up <- get_GOp_one_sided("NAFLD", "Cirrhosis", "less", "CC")

goBP_stat_N_C_dw <- get_GOp_one_sided("NAFLD", "Cirrhosis", "greater", "BP")
goMF_stat_N_C_dw <- get_GOp_one_sided("NAFLD", "Cirrhosis", "greater", "MF")
goCC_stat_N_C_dw <- get_GOp_one_sided("NAFLD", "Cirrhosis", "greater", "CC")

GOdt_up <- rbindlist(list(goBP_stat_N_H_up,
                          goMF_stat_N_H_up,
                          goCC_stat_N_H_up,
                          goBP_stat_N_C_up,
                          goMF_stat_N_C_up,
                          goCC_stat_N_C_up))

GOdt_dw <- rbindlist(list(goBP_stat_N_H_dw,
                          goMF_stat_N_H_dw,
                          goCC_stat_N_H_dw,
                          goBP_stat_N_C_dw,
                          goMF_stat_N_C_dw,
                          goCC_stat_N_C_dw))

GOdt_up[, rnk := 1:.N, by = c("Condition","Ontology")]
GOdt_up[, top_go := min(rnk), by = c("GO.ID","Ontology")]
n_go <- 10
GOdt_up <- GOdt_up[top_go <= n_go]
GOdt_up[, direction := "Up-regulated"]
GOdt_up <- GOdt_up[order(classicFisher)]
GOdt_up$Term <- factor(GOdt_up$Term, levels = rev(GOdt_up[!duplicated(Term), Term]))


GOdt_dw[, rnk := 1:.N, by = c("Condition","Ontology")]
GOdt_dw[, top_go := min(rnk), by = c("GO.ID","Ontology")]
n_go <- 10
GOdt_dw <- GOdt_dw[top_go <= n_go]
GOdt_dw[, direction := "Down-regulated"]
GOdt_dw <- GOdt_dw[order(padj)]
GOdt_dw$Term <- factor(GOdt_dw$Term, levels = rev(GOdt_dw[!duplicated(Term), Term]))

min_p <- min(GOdt_up$padj, GOdt_dw$padj)

facet_labs <- c("Biological process",
                "Cellular component",
                "Molecular function")
names(facet_labs) <- c("BP","CC","MF")

p_go_up <- ggplot(GOdt_up, aes(Condition, Term, fill = -log10(padj))) +
  geom_tile(col = "black") +
  geom_point(aes(shape = ifelse(padj <= 0.05, "B","A")), fill = "white",
             color = "black", size = 1) +
  scale_fill_viridis_c(option = "inferno", limits = c(0, -log10(min_p))) +
  scale_shape_manual(values = c(A = NA, B = 21), drop = F) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  facet_grid(Ontology~direction, scales = "free_y", space = "free",
             labeller = labeller(Ontology = facet_labs)) +
  labs(x = "Contrast", y = "GO-Term", fill = "-log10(p-value)") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.background = element_rect(fill = "black"),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "bottom",
        legend.justification = 1
  ) +
  guides(shape = "none")
p_go_up

p_go_dw <- ggplot(GOdt_dw, aes(Condition, Term, fill = -log10(padj))) +
  geom_tile(col = "black") +
  geom_point(aes(shape = ifelse(padj <= 0.05, "B","A")), fill = "white",
             color = "black", size = 1) +
  scale_fill_viridis_c(option = "inferno", limits = c(0, -log10(min_p))) +
  scale_shape_manual(values = c(A = NA, B = 21), drop = F) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  facet_grid(Ontology~direction, scales = "free_y", space = "free",
             labeller = labeller(Ontology = facet_labs)) +
  labs(x = "Condition", y = "GO-Term", fill = "-log10(p^-value)") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.background = element_rect(fill = "black"),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none"
  ) +
  guides(shape = "none")
p_go_dw

#p_go_comb <- ggarrange(p_go_up, p_go_dw, nrow = 1)
p_go_comb <- ggarrange(p_go_up, nrow = 1) # no significant GO terms in downregs

ggsave("output/plots/rna_GOenrichments_compared_to_NAFLD.pdf",p_go_comb,
       width = 3.75, height = 11)
ggsave("output/plots/submission1/FigS4.pdf",p_go_comb,
       width = 3.75, height = 11)

