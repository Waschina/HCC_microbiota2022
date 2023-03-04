# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Make plot for cohort overview #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(bit64))
suppressMessages(library(ggtext))

myCol <- c(NAFLD = "#2c7bb6", HCC = "#d7191c", Cirrhosis = "#fdae61")

# load 16S DNA sample info
spls_16S_DNA <- fread("data/meta/samples_16S_DNA.csv")[include == 1]
spls_16S_DNA[, case_ID := paste(Patient_ID, Condition, sep = "_")]
spls_16S_DNA <- spls_16S_DNA[sample != 576]

# load Kiel rna-seq sample info
spls_rnaseq <- fread("data/rnaseq/sample_information_edit.csv")
spls_rnaseq <- merge(spls_rnaseq, spls_16S_DNA[,.(Patient_ID, sample)],
                  by.x = "sample_16S_data", by.y = "sample",
                  all.x = T)
spls_rnaseq[is.na(Patient_ID), Patient_ID := .I]
spls_rnaseq[, case_ID := paste(Patient_ID, Condition, sep = "_")]

# # load 16S RNA sample info for biospies
# spls_16S_RNA <- fread("../HCC_liver16S/data/dada/DT_alpha_diversity.tsv")
# spls_16S_RNA <- spls_16S_RNA[,.(sample)]
# spls_16S_RNA[, sample := gsub("HCC\\.","",sample)]
# spls_16S_RNA[, sample := gsub("\\.","-",sample)]
# spls_16S_RNA <- merge(spls_16S_RNA, spls_rnaseq[,.(sample = external_name, Patient_ID, Condition, case_ID)])

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Construct a table with all cases as single rows/entries #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
caseDT <- rbindlist(list(spls_16S_DNA[,.(case_ID, Patient_ID, Condition)],
                         spls_rnaseq[,.(case_ID, Patient_ID, Condition)]))
caseDT <- caseDT[!duplicated(case_ID)]
caseDT[,c("fecal_16S","blood_16S","biopsy_16S","biopsy_rna") := F]



# fecal 16S samples
tmp <- spls_16S_DNA[Sample_type == "faeces", case_ID]
caseDT[, fecal_16S := case_ID %in% tmp]

# blood 16S samples
tmp <- spls_16S_DNA[Sample_type == "blood", case_ID]
caseDT[, blood_16S := case_ID %in% tmp]

# biopsy 16S samples
tmp <- spls_16S_DNA[Sample_type == "liver", case_ID]
caseDT[, biopsy_16S := case_ID %in% tmp]

# # biopsy 16S-Kiel samples
# tmp <- spls_16S_RNA$case_ID
# caseDT[, biopsy_16SKiel := case_ID %in% tmp]

# biopsy RNAseq samples
tmp <- spls_rnaseq$case_ID
caseDT[, biopsy_rna := case_ID %in% tmp]

n <- nrow(caseDT)

# write output as cohort overview table
fwrite(caseDT, "output/files/#cohort_overview.csv")

# Some PAtients are in the list with different Conditions. This is because the patients
# came to the clinic at different times, when the disease has progressed e.g from cirr-
# hosis to HCC. The following two lines are only to get a quick picture on the sample
# distribution of those cases.
dupl_pat <- caseDT[duplicated(Patient_ID), Patient_ID]
caseDT[Patient_ID %in% dupl_pat][order(Patient_ID)]

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# make plot as overview #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
caseDT[, spl_type_cov := fecal_16S + blood_16S + biopsy_16S + biopsy_rna, by = "case_ID"]
caseDT <- caseDT[order(Condition, -spl_type_cov)]
case_order <- caseDT$case_ID

caseDT_long <- melt(caseDT, id.vars = c("case_ID","Patient_ID","Condition"),
                    measure.vars = c("fecal_16S","blood_16S","biopsy_16S","biopsy_rna"),
                    value.name = "sample_presence", variable.name = "sample_type")

n_cond <- c(HCC = 0, Cirrhosis = 0, NAFLD = 0)
n_cond["HCC"] <- caseDT[Condition=="HCC",.N]
n_cond["Cirrhosis"] <- caseDT[Condition=="Cirrhosis",.N]
n_cond["NAFLD"] <- caseDT[Condition=="NAFLD",.N]

caseDT_long$case_ID <- factor(caseDT_long$case_ID, levels = case_order)
caseDT_long$Condition <- factor(caseDT_long$Condition, levels = c("NAFLD","Cirrhosis","HCC"),
                                labels = paste0("<b><i style='color:",myCol[c("NAFLD","Cirrhosis","HCC")],"'>",c("NAFLD","Cirrhosis","HCC"),"</i></b> (",n_cond[c("NAFLD","Cirrhosis","HCC")],")"))



caseDT_long[sample_type == "fecal_16S", sample_type := "16S from DNA (feces)"]
caseDT_long[sample_type == "blood_16S", sample_type := "16S from DNA (blood)"]
caseDT_long[sample_type == "biopsy_16S", sample_type := "16S from DNA (liver)"]
#caseDT_long[sample_type == "biopsy_16SKiel", sample_type := "16S from RNA (liver)"]
caseDT_long[sample_type == "biopsy_rna", sample_type := "Transcriptome (liver)"]

caseDT_long[, strip_col := myCol[Condition]]
#caseDT_long <- caseDT_long[sample_type != "16S from RNA (liver)"]

caseDT_long$sample_type <- factor(caseDT_long$sample_type, levels = c("16S from DNA (feces)",
                                                                      "16S from DNA (blood)",
                                                                      "16S from DNA (liver)",
                                                                      #"16S from RNA (liver)",
                                                                      "Transcriptome (liver)"))

p <- ggplot(caseDT_long, aes(case_ID, sample_type, fill = sample_presence, col = sample_presence)) +
  #geom_tile() +
  geom_point(aes(shape = sample_presence)) +
  facet_grid(.~Condition, space = "free", scales = "free") +
  scale_fill_manual(values = c(NA,"black")) +
  scale_color_manual(values = c(NA,"black")) +
  scale_shape_manual(values = c(NA, 21)) +
  #scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  labs(x = paste0("Cases (n = ",n,")"), y = "Sample type") +
  theme(axis.text.x = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text = element_markdown(),
        strip.background = element_blank(),
        axis.text.y = element_text(color = "black"),
        legend.position = "none")
p

ggsave("output/plots/cohort_overview.pdf", plot = p, width = 14, height = 1.4)
ggsave("output/plots/submission1/FigS1.pdf", plot = p, width = 14, height = 1.4)


pval_labs <- function(pval, ns.lab = "") {
  res <- ifelse(pval < 0.0001, "****",
                ifelse(pval < 0.001, "***",
                       ifelse(pval < 0.01, "**",
                              ifelse(pval < 0.05, "*",ns.lab))))
}
