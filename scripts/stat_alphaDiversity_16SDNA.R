myCol <- c(NAFLD = "#2c7bb6", Cirrhosis = "#fdae61", HCC = "#d7191c")

alpha_metrics <- c("Shannon","J.eveness")
my_comparisons <- list(c("NAFLD","HCC"),
                       c("NAFLD","Cirrhosis"),
                       c("Cirrhosis","HCC"))
cst_fs <- 3

#––––––––––––––––––––––––––#
# Alpha-Diversity analysis #
#––––––––––––––––––––––––––#
library(data.table)
library(egg)
library(ggpubr)

spl_meta <- fread("data/meta/samples_16S_DNA.csv")
spl_meta <- spl_meta[include == TRUE]
spl_meta <- spl_meta[sample != 576]
dt_adiv  <- fread("data/dada/DT_alpha_diversity.tsv")

dt <- merge(spl_meta, dt_adiv, by = "sample")
# dt$Condition <- factor(dt$Condition, levels = c("NAFLD","Cirrhosis","HCC"))
# dt$Sample_type <- factor(dt$Sample_type, levels = c("faeces","blood","liver"))

# alpha_metrics <- c("Shannon","Simpson","Simpson_inv","Fisher.alpha","S.richness","J.eveness")
# p_list <- list()
# for(am in alpha_metrics) {
#   dt$y <- dt[[am]]
#   
#   legend_pos <- ifelse(am == "InvSimpson","bottom","none")
#   
#   cst_fs <- 3
#   
#   p <- ggplot(dt, aes(Condition, y, fill = Condition)) +
#     geom_boxplot(outlier.shape = NA) + 
#     geom_jitter(alpha = 0.3, width = 0.25, size = 0.75) +
#     scale_fill_manual(values = myCol) +
#     labs(y = am) +
#     facet_grid(.~Sample_type) + 
#     theme_bw() +
#     theme(panel.grid = element_blank(),
#           strip.background = element_blank(),
#           strip.text = element_text(color = "black", face = "bold"),
#           axis.text = element_text(color = "black"),
#           axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
#           legend.position = legend_pos) +
#     stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
#                        symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
#                        tip.length = 0, #hide.ns = T,
#                        size = cst_fs)
#   p$layers[[3]]$aes_params$textsize <- cst_fs
#   p_list[[am]] <- p
#   
# }
# 
# p_alpha <- egg::ggarrange(plots = p_list, ncol = 3)
# ggsave("analysis/v1/plots/16S_alphaDiversity.pdf", plot = p_alpha, width = 14, height = 8)

#––––––––––––––––––––––––––––––––––––––––––––––#
# Alpha diversity for 16S data from biopsy RNA #
#––––––––––––––––––––––––––––––––––––––––––––––#

# rna_adiv <- fread("../HCC_liver16S/data/dada/DT_alpha_diversity.tsv")
# rna_si <- fread("data/rnaseq/sample_information_edit.csv")
# rna_adiv[, sample := gsub("^HCC\\.","",sample)]
# rna_adiv[, sample := gsub("\\.","-",sample)]
# rna_si <- rna_si[external_name %in% rna_adiv$sample]
# 
# rna_adiv <- merge(rna_adiv, rna_si, by.x = "sample", by.y = "external_name", all.x = T)

# prna_list <- list()
# for(am in alpha_metrics) {
#   rna_adiv$y <- rna_adiv[[am]]
#   
#   legend_pos <- ifelse(am == "InvSimpson","bottom","none")
#   
#   cst_fs <- 3
#   
#   p <- ggplot(rna_adiv, aes(Condition, y, fill = Condition)) +
#     geom_boxplot(outlier.shape = NA) + 
#     geom_jitter(alpha = 0.3, width = 0.25, size = 0.75) +
#     scale_fill_manual(values = myCol) +
#     labs(y = am) +
#     theme_bw() +
#     theme(panel.grid = element_blank(),
#           strip.background = element_blank(),
#           strip.text = element_text(color = "black", face = "bold"),
#           axis.text = element_text(color = "black"),
#           axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
#           legend.position = legend_pos) +
#     stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
#                        symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
#                        tip.length = 0, #hide.ns = T,
#                        size = cst_fs)
#   
#   p$layers[[3]]$aes_params$textsize <- cst_fs
#   prna_list[[am]] <- p
# }
# 
# prna_alpha <- egg::ggarrange(plots = prna_list, ncol = 3)
# ggsave("analysis/v1/plots/16S_alphaDiversity_liver16S_rna.pdf", plot = prna_alpha, width = 5, height = 6)



# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Combine alpha div metric for all  #
# matrices and DNA/RNA              #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
dt_comb <- copy(dt)
dt_comb[, sample_type2 := paste(Sample_type)]
#dt_comb <- rbind(dt_comb, rna_adiv, fill = TRUE)
#dt_comb[is.na(sample_type2), sample_type2 := "liver (RNA)"]

dt_comb$Condition <- factor(dt_comb$Condition, levels = c("NAFLD","Cirrhosis","HCC"))
dt_comb$sample_type2 <- factor(dt_comb$sample_type2, levels = c("faeces",
                                                               "blood",
                                                               "liver"))

# Shannon
p_shannon <- ggplot(dt_comb, aes(Condition, Shannon, fill = Condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.3, width = 0.25, size = 0.75) +
  scale_fill_manual(values = myCol) +
  labs(y = "Shannon index") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        legend.position = "none") +
  facet_grid(. ~ sample_type2) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("****", "***", "**", "*", "ns")),
                     tip.length = 0.01, #hide.ns = T,
                     size = cst_fs)
p_shannon

# Evenness
p_even <- ggplot(dt_comb, aes(Condition, J.eveness, fill = Condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.3, width = 0.25, size = 0.75) +
  scale_fill_manual(values = myCol) +
  labs(y = "Species evenness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        legend.position = "right") +
  facet_grid(. ~ sample_type2) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("****", "***", "**", "*", "ns")),
                     tip.length = 0.01, #hide.ns = T,
                     size = cst_fs)
p_even

p_comb <- egg::ggarrange(p_shannon, p_even, nrow = 1, draw = FALSE)
ggsave("output/plots/16S_alpha_diversity.pdf", plot = p_comb,
       width = 4.7, height = 7)
p_comb_alpha <- p_comb

