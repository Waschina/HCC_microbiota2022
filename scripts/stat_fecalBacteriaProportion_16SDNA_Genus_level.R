myCol <- c(NAFLD = "#2c7bb6", Cirrhosis = "#fdae61", HCC = "#d7191c")

#––––––––––––––––––––––––––––––––––#
# Fecal bacteria in blood/liver    #
# based on Genus aggregated counts #
#––––––––––––––––––––––––––––––––––#

# Bsicially here, we identify faecal genera instead of faecal ASVs
library(ggplot2)
library(ggpubr)
library(ggtext)
library(data.table)

# sample invo and asv counts
asv_tab <- read.table("data/dada/asv_tab.tsv", check.names = F)
spl_meta <- fread("data/meta/samples_16S_DNA.csv")
spl_meta <- spl_meta[include == TRUE]
spl_meta <- spl_meta[sample != 576]
spl_meta <- spl_meta[sample %in% colnames(asv_tab)]
#spl_meta <- spl_meta[underlying_condition == "NAFLD"]
asv_tab <- asv_tab[, colnames(asv_tab) %in% spl_meta$sample]


# load taxonomic predictions
asv_tax <- read.table("data/dada/asv_tax.tsv", sep = "\t")
tmp <- rownames(asv_tax)
asv_tax <- data.table(asv_tax)
asv_tax$cluster <- tmp
#asv_tax[Genus == "Unknown", Genus := paste(Genus, Family)]

asv_tax_tmp <- copy(asv_tax)
asv_tax_tmp <- asv_tax_tmp[cluster %in% rownames(asv_tab)]
# asv_tax_tmp[genus == "", genus := paste0("unknown ",family)]
# asv_tax_tmp[species == "", species := paste0("unknown ",genus)]
asv_tax_tmp[, Kingdom2 := paste0("k__",Kingdom)]
asv_tax_tmp[, Phylum2 := paste0("k__",Kingdom,"|p__",Phylum)]
asv_tax_tmp[, Class2 := paste0("k__",Kingdom,"|p__",Phylum,"|c__",Class)]
asv_tax_tmp[, Order2 := paste0("k__",Kingdom,"|p__",Phylum,"|c__",Class,"|o__",Order)]
asv_tax_tmp[, Family2 := paste0("k__",Kingdom,"|p__",Phylum,"|c__",Class,"|o__",Order,"|f__",Family)]
asv_tax_tmp[, Genus2 := paste0("k__",Kingdom,"|p__",Phylum,"|c__",Class,"|o__",Order,"|f__",Family,"|g__",Genus)]

# make aggregated count table
ilvl <- "Genus2"

lvls <- unique(asv_tax_tmp[[ilvl]])
cat(ilvl, "-", length(lvls),"\n")

tmpcts <- matrix(NA_integer_, nrow = length(lvls), ncol = ncol(asv_tab))
colnames(tmpcts) <- colnames(asv_tab)
rownames(tmpcts) <- lvls

for(itax in lvls) {
  tmpclust <- asv_tax_tmp[get(ilvl) == itax, cluster]
  tmpcts[itax, ] <- colSums(asv_tab[tmpclust,, drop = FALSE])
}


# define faecal genera
n_cts <- colSums(tmpcts)

# fecal_tax_N <- names(which(apply(tmpcts[,spl_meta[Sample_type == "faeces" & Condition == "NAFLD",sample]],1, function(x) sum(x/n_cts[names(x)] >= 0.001)) >= 1))
# fecal_tax_C <- names(which(apply(tmpcts[,spl_meta[Sample_type == "faeces" & Condition == "Cirrhosis",sample]],1, function(x) sum(x/n_cts[names(x)] >= 0.001)) >= 1))
# fecal_tax_H <- names(which(apply(tmpcts[,spl_meta[Sample_type == "faeces" & Condition == "HCC",sample]],1, function(x) sum(x/n_cts[names(x)] >= 0.001)) >= 1))

# rel > 5%
fecal_tax_N <- names(which(apply(tmpcts[,spl_meta[Sample_type == "faeces" & Condition == "NAFLD",sample]],1, function(x) sum(x/n_cts[names(x)] >= 0.001)) >= spl_meta[Sample_type == "faeces" & Condition == "NAFLD",.N] * 0.05))
fecal_tax_C <- names(which(apply(tmpcts[,spl_meta[Sample_type == "faeces" & Condition == "Cirrhosis",sample]],1, function(x) sum(x/n_cts[names(x)] >= 0.001)) >= spl_meta[Sample_type == "faeces" & Condition == "Cirrhosis",.N] * 0.05))
fecal_tax_H <- names(which(apply(tmpcts[,spl_meta[Sample_type == "faeces" & Condition == "HCC",sample]],1, function(x) sum(x/n_cts[names(x)] >= 0.001)) >= spl_meta[Sample_type == "faeces" & Condition == "HCC",.N] * 0.05))

dt.fecal_new <- list()
# NAFLD
dt.fecal_new[["N"]] <- data.table(sample = spl_meta[Condition == "NAFLD", sample])
dt.fecal_new[["N"]][, prop.fecal := apply(tmpcts[,sample], 2, function(x) sum(x[fecal_tax_N])/sum(x))]
# NAFLD
dt.fecal_new[["C"]] <- data.table(sample = spl_meta[Condition == "Cirrhosis", sample])
dt.fecal_new[["C"]][, prop.fecal := apply(tmpcts[,sample], 2, function(x) sum(x[fecal_tax_C])/sum(x))]
# HCC
dt.fecal_new[["H"]] <- data.table(sample = spl_meta[Condition == "HCC", sample])
dt.fecal_new[["H"]][, prop.fecal := apply(tmpcts[,sample], 2, function(x) sum(x[fecal_tax_H])/sum(x))]

dt.fecal <- rbindlist(dt.fecal_new)
dt.fecal <- merge(dt.fecal, spl_meta)

dt.fecal$Condition <- factor(dt.fecal$Condition, levels = c("NAFLD","Cirrhosis","HCC"))

saveRDS(dt.fecal, file = "output/files/dt.fecal.RDS")

# Plotting
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
                     label.y = c(0.5, 0.54, 0.46)+0.25,
                     label.y.npc = c("middle", "middle", "middle"),
                     tip.length = 0.01)
p$layers[[3]]$aes_params$textsize <- 2.5
p

ggsave("output/plots/16S_fecalASVs_liver_blood_.pdf", plot = p,
       width = 4, height =  4)
ggsave("output/plots/submission1/Fig2.pdf", plot = p,
       width = 4, height =  4)
