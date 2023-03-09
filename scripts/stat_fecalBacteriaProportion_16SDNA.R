myCol <- c(NAFLD = "#2c7bb6", Cirrhosis = "#fdae61", HCC = "#d7191c")

#–––––––––––––––––––––––––––––––#
# Fecal bacteria in blood/liver #
#–––––––––––––––––––––––––––––––#
library(ggplot2)
library(ggpubr)
library(ggtext)
library(data.table)


asv_tab <- read.table("data/dada/asv_tab.tsv", check.names = F)
spl_meta <- fread("data/meta/samples_16S_DNA.csv")
spl_meta <- spl_meta[include == TRUE]
spl_meta <- spl_meta[sample != 576]
#spl_meta <- spl_meta[underlying_condition == "NAFLD"]
#spl_meta <- spl_meta[Condition == "NAFLD" | underlying_condition == "HCV"]
spl_meta[, table(Condition, Sample_type)]

# n_cts <- colSums(asv_tab)[1]

# Definition of faecal ASVs: The ASV has an relative abundance of >= 1% in at least 1 faecal sample.
#fecal_asvs <- names(which(apply(asv_tab[,spl_meta[Sample_type == "faeces",sample]],1, function(x) sum(x/n_cts >= 0.01)) >= 1))

# Definition of faecal ASVs: The ASV has an relative abundance of >= 1% in at least 1 faecal sample from the same etiology group.
n_cts <- colSums(asv_tab)
# abs > 1
fecal_asvs_N <- names(which(apply(asv_tab[,spl_meta[Sample_type == "faeces" & Condition == "NAFLD",sample]],1, function(x) sum(x/n_cts[names(x)] >= 0.001)) >= 1))
fecal_asvs_C <- names(which(apply(asv_tab[,spl_meta[Sample_type == "faeces" & Condition == "Cirrhosis",sample]],1, function(x) sum(x/n_cts[names(x)] >= 0.001)) >= 1))
fecal_asvs_H <- names(which(apply(asv_tab[,spl_meta[Sample_type == "faeces" & Condition == "HCC",sample]],1, function(x) sum(x/n_cts[names(x)] >= 0.001)) >= 1))
# # rel > 5%
# fecal_asvs_N <- names(which(apply(asv_tab[,spl_meta[Sample_type == "faeces" & Condition == "NAFLD",sample]],1, function(x) sum(x/n_cts[names(x)] >= 0.001)) >= spl_meta[Sample_type == "faeces" & Condition == "NAFLD",.N] * 0.05))
# fecal_asvs_C <- names(which(apply(asv_tab[,spl_meta[Sample_type == "faeces" & Condition == "Cirrhosis",sample]],1, function(x) sum(x/n_cts[names(x)] >= 0.001)) >= spl_meta[Sample_type == "faeces" & Condition == "Cirrhosis",.N] * 0.05))
# fecal_asvs_H <- names(which(apply(asv_tab[,spl_meta[Sample_type == "faeces" & Condition == "HCC",sample]],1, function(x) sum(x/n_cts[names(x)] >= 0.001)) >= spl_meta[Sample_type == "faeces" & Condition == "HCC",.N] * 0.05))


# dt.fecal <- data.table(sample = colnames(asv_tab),
#                        prop.fecal = apply(asv_tab[fecal_asvs,],2,sum)/n_cts)
# dt.fecal <- merge(dt.fecal, spl_meta)
# dt.fecal[Condition == "cirrhosis", Condition := "Cirrhosis"]

dt.fecal_new <- list()
# NAFLD
dt.fecal_new[["N"]] <- data.table(sample = spl_meta[Condition == "NAFLD", sample])
dt.fecal_new[["N"]][, prop.fecal := apply(asv_tab[,sample], 2, function(x) sum(x[fecal_asvs_N])/sum(x))]
# NAFLD
dt.fecal_new[["C"]] <- data.table(sample = spl_meta[Condition == "Cirrhosis", sample])
dt.fecal_new[["C"]][, prop.fecal := apply(asv_tab[,sample], 2, function(x) sum(x[fecal_asvs_C])/sum(x))]
# HCC
dt.fecal_new[["H"]] <- data.table(sample = spl_meta[Condition == "HCC", sample])
dt.fecal_new[["H"]][, prop.fecal := apply(asv_tab[,sample], 2, function(x) sum(x[fecal_asvs_H])/sum(x))]

dt.fecal <- rbindlist(dt.fecal_new)
dt.fecal <- merge(dt.fecal, spl_meta)

dt.fecal$Condition <- factor(dt.fecal$Condition, levels = c("NAFLD","Cirrhosis","HCC"))

saveRDS(dt.fecal, file = "output/files/dt.fecal.RDS")

# use these stats
wilcox.test(dt.fecal[Sample_type == "blood" & Condition == "NAFLD",     prop.fecal],
            dt.fecal[Sample_type == "blood" & Condition %in% c("Cirrhosis","HCC"), prop.fecal], alternative = "less")
wilcox.test(dt.fecal[Sample_type == "liver" & Condition == "NAFLD",     prop.fecal],
            dt.fecal[Sample_type == "liver" & Condition %in% c("Cirrhosis","HCC"), prop.fecal], alternative = "less")
wilcox.test(dt.fecal[Sample_type == "faeces" & Condition == "NAFLD",     prop.fecal],
            dt.fecal[Sample_type == "faeces" & Condition %in% c("Cirrhosis","HCC"), prop.fecal], alternative = "less")


wilcox.test(dt.fecal[Sample_type == "blood" & Condition == "NAFLD",     prop.fecal],
            dt.fecal[Sample_type == "blood" & Condition %in% c("Cirrhosis"), prop.fecal], alternative = "less")
wilcox.test(dt.fecal[Sample_type == "blood" & Condition == "NAFLD",     prop.fecal],
            dt.fecal[Sample_type == "blood" & Condition %in% c("HCC"), prop.fecal], alternative = "less")
wilcox.test(dt.fecal[Sample_type == "liver" & Condition == "NAFLD",     prop.fecal],
            dt.fecal[Sample_type == "liver" & Condition %in% c("Cirrhosis"), prop.fecal], alternative = "less")
wilcox.test(dt.fecal[Sample_type == "liver" & Condition == "NAFLD",     prop.fecal],
            dt.fecal[Sample_type == "liver" & Condition %in% c("HCC"), prop.fecal], alternative = "less")

my_comparisons <- list(c("NAFLD","HCC"),
                       c("NAFLD","Cirrhosis"),
                       c("HCC","Cirrhosis"))

p <- ggplot(dt.fecal[Sample_type != "faeces"], aes(Condition, prop.fecal, fill = Condition)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(alpha = 0.3, width = 0.25, size = 0.75) +
  scale_fill_manual(values = myCol) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  facet_grid(.~Sample_type) +
  labs(y = "Proportion of faecal bacteria") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     #method.args = list(alternative = "less"),
                     label.y = c(0.5, 0.54, 0.46)+0.1,
                     label.y.npc = c("middle", "middle", "middle"),
                     tip.length = 0.01)
p$layers[[3]]$aes_params$textsize <- 2.5
p

ggsave("output/plots/16S_fecalASVs_liver_blood.pdf", plot = p,
       width = 4, height =  4)
ggsave("output/plots/submission1/Fig2.pdf", plot = p,
       width = 4, height =  4)
