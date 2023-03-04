# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Differnetial Genus abundance analysis #
# via DESeq2                            #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
library(ggtext)
library(ggrepel)
library(lefser)

myCol <- c(NAFLD = "#2c7bb6", Cirrhosis = "#fdae61", HCC = "#d7191c")
comparison_names <- c("<i style='color:#fd8b1b'>Cirrhosis</i> vs. <i style='color:#d7191c'>HCC</i>",
                      "<i style='color:#2c7bb6'>NAFLD</i> vs. <i style='color:#fd8b1b'>Cirrhosis</i>",
                      "<i style='color:#2c7bb6'>NAFLD</i> vs. <i style='color:#d7191c'>HCC</i>")

asv_tab <- read.table("data/dada/asv_tab.tsv", check.names = F)
spl_meta <- fread("data/meta/samples_16S_DNA.csv")
spl_meta <- spl_meta[include == TRUE]
spl_meta <- spl_meta[sample != 576]
spl_meta <- spl_meta[sample %in% colnames(asv_tab)]
asv_tab <- asv_tab[, colnames(asv_tab) %in% spl_meta$sample]


#genus_tab <- read.table("data/dada/tax_tab_Genus.tsv", check.names = F, sep = "\t")
#family_tab <- read.table("data/dada/tax_tab_Family.tsv", check.names = F, sep = "\t")

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
lvlcts <- list()
for(ilvl in c("Phylum2","Class2","Order2","Family2","Genus2")) {
  lvls <- unique(asv_tax_tmp[[ilvl]])
  cat(ilvl, "-", length(lvls),"\n")
  
  lvlcts[[ilvl]] <- matrix(NA_integer_, nrow = length(lvls), ncol = ncol(asv_tab))
  colnames(lvlcts[[ilvl]]) <- colnames(asv_tab)
  rownames(lvlcts[[ilvl]]) <- lvls
  
  for(itax in lvls) {
    tmpclust <- asv_tax_tmp[get(ilvl) == itax, cluster]
    lvlcts[[ilvl]][itax, ] <- colSums(asv_tab[tmpclust,, drop = FALSE])
  }
}
lvlcts <- do.call(rbind, lvlcts)

# split by sample type
taxcts <- list()
s_exp <- list()
for(stype in unique(spl_meta$Sample_type)) {
  meta2 <- copy(spl_meta[sample %in% colnames(lvlcts)])
  meta2 <- meta2[Sample_type == stype]
  meta2[Condition == "NAFLD", Condition := "1_NAFLD"]
  meta2[Condition == "Cirrhosis", Condition := "2_Cirrhosis"]
  meta2[Condition == "HCC", Condition := "3_HCC"]
  #meta2$Condition <- factor(meta2$Condition, levels = c("NAFLD","Cirrhosis","HCC"))
  spls <- spl_meta[, sample]
  taxcts[[stype]] <- lvlcts[, spls]
  reltax <- which(apply(taxcts[[stype]], 1, function(x) any(x > 0)))
  taxcts[[stype]] <- taxcts[[stype]][reltax,]
  taxcts[[stype]] <- t(t(taxcts[[stype]])/colSums(taxcts[[stype]])) * 1e6
  taxcts[[stype]] <- taxcts[[stype]][,meta2[Sample_type==stype, sample]]
  s_exp[[stype]] <- SummarizedExperiment(assays = SimpleList(counts = taxcts[[stype]]),
                                       colData = meta2[,.(sample, Condition)])
}

set.seed(24118)
lefse_res_NH <- lapply(s_exp, function(x) {
  s_exp_tmp <- x[,x$Condition %in% c("1_NAFLD","3_HCC")]
  res <- lefser(s_exp_tmp, groupCol = "Condition",kruskal.threshold = 0.05, wilcox.threshold = 0.05)
  return(res)
})
lefse_res_NC <- lapply(s_exp, function(x) {
  s_exp_tmp <- x[,x$Condition %in% c("1_NAFLD","2_Cirrhosis")]
  res <- lefser(s_exp_tmp, groupCol = "Condition",kruskal.threshold = 0.05, wilcox.threshold = 0.05)
  return(res)
})
lefse_res_HC <- lapply(s_exp, function(x) {
  s_exp_tmp <- x[,x$Condition %in% c("3_HCC","2_Cirrhosis")]
  res <- lefser(s_exp_tmp, groupCol = "Condition",kruskal.threshold = 0.05, wilcox.threshold = 0.05)
  return(res)
})

lefse_res_NH <- rbindlist(lefse_res_NH, idcol = "Sample_type")
lefse_res_NH$contrast <- comparison_names[3]
lefse_res_NH[scores < 0, Condition := "NAFLD"]
lefse_res_NH[scores > 0, Condition := "HCC"]

lefse_res_NC <- rbindlist(lefse_res_NC, idcol = "Sample_type")
lefse_res_NC$contrast <- comparison_names[2]
lefse_res_NC[scores < 0, Condition := "NAFLD"]
lefse_res_NC[scores > 0, Condition := "Cirrhosis"]

lefse_res_HC <- rbindlist(lefse_res_HC, idcol = "Sample_type")
lefse_res_HC$contrast <- comparison_names[1]
lefse_res_HC[scores < 0, Condition := "Cirrhosis"]
lefse_res_HC[scores > 0, Condition := "HCC"]

lefse_res <- rbindlist(list(lefse_res_NH,lefse_res_NC,lefse_res_HC))
lefse_res$Sample_type <- factor(lefse_res$Sample_type, levels = c("faeces","blood","liver"))
lefse_res$Condition <- factor(lefse_res$Condition, levels = c("NAFLD","Cirrhosis","HCC"))

lefse_res[, Names := gsub("^k__Bacteria\\|","", Names)]

lefse_res[, Names := gsub("p__","<b style='color:#bbbbbb'>[p]</b>", Names)]
lefse_res[, Names := gsub("c__","<b style='color:#bbbbbb'>[c]</b>", Names)]
lefse_res[, Names := gsub("o__","<b style='color:#bbbbbb'>[o]</b>", Names)]
lefse_res[, Names := gsub("f__","<b style='color:#bbbbbb'>[f]</b>", Names)]
lefse_res[, Names := gsub("g__","<b style='color:#bbbbbb'>[g]</b>", Names)]
lefse_res[, Names := gsub("\\|"," ", Names)]
rel_tax <- lefse_res[abs(scores) > 3.5, Names]

p_lefse <- ggplot(lefse_res[Names %in% rel_tax], aes(scores, Names, fill = Condition)) +
  geom_col() +
  geom_vline(xintercept = 0, lty = 3) +
  facet_grid(contrast~Sample_type, space = "free_y", scales = "free_y") +
  scale_x_continuous(breaks = c(-3.5,0,3.5)) +
  scale_fill_manual(values = myCol) +
  labs(x = "LDA score (log10)") +
  theme_bw() +
  theme(strip.text.y = element_markdown(),
        strip.text.x = element_text(face = "bold"),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_markdown(),
        legend.position = "bottom")

ggsave("output/plots/16S_LEFSe.pdf", plot = p_lefse,
       width = 10, height = 9)
ggsave("output/plots/submission1/Fig3.pdf", plot = p_lefse,
       width = 10, height = 9)
#-------------------#
# Composition plots #
#-------------------#
ilvl <- "Class2"

lvls <- unique(asv_tax_tmp[[ilvl]])
cat(ilvl, "-", length(lvls),"\n")

tmpcts <- matrix(NA_integer_, nrow = length(lvls), ncol = ncol(asv_tab))
colnames(tmpcts) <- colnames(asv_tab)
rownames(tmpcts) <- lvls

for(itax in lvls) {
  tmpclust <- asv_tax_tmp[get(ilvl) == itax, cluster]
  tmpcts[itax, ] <- colSums(asv_tab[tmpclust,, drop = FALSE])
}

rownames(tmpcts) <- gsub("^.*__","",rownames(tmpcts))
tmpcts <- data.table(as.table(tmpcts))
colnames(tmpcts) <- c("Taxon","sample","count")
tmpcts <- tmpcts[, .(count = sum(count)), by = c("Taxon","sample")]
tmpcts[, rel.abun := count/sum(count), by = sample]

rel_taxons <- tmpcts[, sum(sqrt(rel.abun)), by = Taxon][order(-V1)][1:12, Taxon]
tmpcts[!(Taxon %in% rel_taxons), Taxon := "Other"]
rel_taxons <- c(rel_taxons, "Other")
tmpcts <- tmpcts[, .(rel.abun = sum(rel.abun)), by = .(Taxon, sample)]
tmpcts <- merge(tmpcts, spl_meta)
colcode <- c(as.character(MetBrewer::met.brewer("Redon",12)),"#151515")
names(colcode) <- rel_taxons

tmpcts$Taxon <- factor(tmpcts$Taxon, levels = rel_taxons)
tmpcts$Condition <- factor(tmpcts$Condition, levels = c("NAFLD","Cirrhosis","HCC"))

p_f <- ggplot(tmpcts[Sample_type=="faeces"], aes(rel.abun, sample, fill = Taxon)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = colcode, drop = FALSE) +
  facet_grid(Condition~Sample_type, scales = "free", space = "free") +
  scale_x_continuous(expand = c(0,0), labels = scales::percent) +
  labs(x = "Rel. Abundance",
       y = "Sample",
       fill = "Class") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"),
        axis.text.y = element_blank(),
        legend.position = "none")

p_b <- ggplot(tmpcts[Sample_type=="blood"], aes(rel.abun, sample, fill = Taxon)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = colcode, drop = FALSE) +
  facet_grid(Condition~Sample_type, scales = "free", space = "free") +
  scale_x_continuous(expand = c(0,0), labels = scales::percent) +
  labs(x = "Rel. Abundance",
       y = "Sample",
       fill = "Class") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"),
        axis.text.y = element_blank(),
        legend.position = "bottom")
  
p_l <- ggplot(tmpcts[Sample_type=="liver"], aes(rel.abun, sample, fill = Taxon)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = colcode, drop = FALSE) +
  facet_grid(Condition~Sample_type, scales = "free", space = "free") +
  scale_x_continuous(expand = c(0,0), labels = scales::percent) +
  labs(x = "Rel. Abundance",
       y = "Sample",
       fill = "Class") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"),
        axis.text.y = element_blank(),
        legend.position = "none")

p_compositions <- egg::ggarrange(p_f, p_b, p_l, nrow = 1,
                                 labels = c("A","B","C"),
                                 draw = FALSE)

ggsave("output/plots/16S_compositions.pdf", plot = p_compositions,
       width = 9, height = 11)

