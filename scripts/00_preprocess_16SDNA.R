min_seq_count <- 2500

cat("Loading packages ...\n")
suppressMessages(library(dada2))
suppressMessages(library(Biostrings))
suppressMessages(library(data.table))
suppressMessages(library(vegan))
cat("\tDADA2:", as.character(packageVersion("dada2")),"\n")
cat("\tBiostrings:", as.character(packageVersion("Biostrings")),"\n")
cat("\tdata.table:", as.character(packageVersion("data.table")),"\n")
cat("\tvegan:", as.character(packageVersion("vegan")),"\n")

# create directory structure if not already there
if(!dir.exists("data/dada/")) {
  dir.create("data/dada/")
}
if(!dir.exists("data/silva_refs/")) {
  dir.create("data/silva_refs/")
}

# get SILVA DB if needed
if(!file.exists("data/silva_refs/silva_nr99_v138.1_train_set.fa.gz")) {
  cat("Downloading silva reference 16S database for taxonomy assignments ...\n")
  download.file("https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1", 
                destfile = "data/silva_refs/silva_nr99_v138.1_train_set.fa.gz")
  download.file("https://zenodo.org/record/4587955/files/SILVA_LICENSE.txt?download=1", 
                destfile = "data/silva_refs/SILVA_LICENSE.txt")
}


# paths to individual files
fastq <- dir("data/fastq/DNA/", full.names = T)
fastq.names <- sapply(strsplit(basename(fastq), "_"), `[`, 1)

# Quality filtering
cat("Quality filtering ...\n")
fq.filt  <- file.path("data/fastq/DNA/", "../DNA-filtered/",
                      fastq.names)
names(fq.filt) <- fastq.names
flt.out <- filterAndTrim(fastq, fq.filt,
                         maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE) 

# dereplicate
cat("Dereplicating sequences ...\n")
fq.derep <- derepFastq(fq.filt)

# learn errors
cat("Learning error rates from sequences ...\n")
fq.err <- learnErrors(fq.derep, multithread=TRUE)
#plotErrors(fq.err, nominalQ=TRUE)

# infer sample composition
cat("Sample inference from amplicon data (dada2) ...\n")
fq.dada <- dada(fq.derep, err=fq.err, multithread=TRUE)

# # remove chimeras
# cat("Removing Chimeras ...\n")
# fq.nochim <- removeBimeraDenovo(fq.dada, multithread=TRUE)

# count table
cat("Contructing ASV-table and representative sequences ...\n")
seqtab <- makeSequenceTable(fq.dada)
seqtab.nochim <- removeBimeraDenovo(seqtab, multithread=TRUE)

sum(seqtab.nochim)/sum(seqtab)

# Assign numbers to ASVs
ASVs <- paste0("asv",1:ncol(seqtab.nochim))
asv.tab <- seqtab.nochim
colnames(asv.tab) <- ASVs
rownames(asv.tab) <- gsub("\\.fastq\\.gz|\\.fq\\.gz","", rownames(asv.tab))
#write.csv(asv.tab, "test.csv")


# write ASV representative sequences (temporary!)
asv_tmp_fna <- tempfile()
asv_fasta <- DNAStringSet(colnames(seqtab.nochim))
names(asv_fasta) <- ASVs
writeXStringSet(asv_fasta, filepath = asv_tmp_fna)

# assign taxonomy
cat("Predicting taxonomic assignments ...\n")
asv.tax <- assignTaxonomy(asv_tmp_fna, refFasta = "data/silva_refs/silva_nr99_v138.1_train_set.fa.gz")
#all(colnames(seqtab.nochim) == rownames(asv.tax))
rownames(asv.tax) <- ASVs
#write.table(asv.tax, "analysis/v1/files/asv_tax_prefilter.txt")

# Removing non-bacterial ASVs
cat("Removing non-bacterial ASVs ...\n")
ind_NA <- which(is.na(asv.tax[,"Kingdom"]))
ind_noBac <- which(asv.tax[,"Kingdom"] != "Bacteria" | asv.tax[,"Order"] == "Chloroplast" | asv.tax[,"Family"] == "Mitochondria")
ind_noBac <- c(ind_NA, ind_noBac)
asv.tax <- asv.tax[-ind_noBac,]
asv.tab <- asv.tab[,-ind_noBac]
asv_fasta <- asv_fasta[-ind_noBac]
cat("\tremoved",length(ind_noBac),"ASVs\n")
cat("\tremaining ASVs:",length(asv_fasta),"\n")

# Downsampling
cat("Filter samples by counts ...\n")
asv.tab <- t(asv.tab)

ind_spl_rm <- which(colSums(asv.tab) < min_seq_count)
if(length(ind_spl_rm) != 0) {
  cat("Removing", length(ind_spl_rm), "samples, which have less than",min_seq_count, "reads:")
  cat("\n\t",paste(colnames(asv.tab)[ind_spl_rm], collapse = "\n\t"),"\n")
}

# source("analysis/v1/scripts/rarefyTable.R")
# if(length(ind_spl_rm) != 0) {
#   asv.tab_rarefied <- rarefyTable(asv.tab[,-ind_spl_rm])
# } else {
#   asv.tab_rarefied <- rarefyTable(asv.tab)
# }
# 
# cat("\tRarefied to",sum(asv.tab_rarefied[,1]),"Counts per sample\n")

# aggregating counts by taxonomy predictions
ind_unknown <- which(is.na(asv.tax) | grepl("Unknown|unknown", asv.tax),arr.ind = T)
asv.tax[ind_unknown] <- "Unknown"

tax.tab <- list()
tax.levels <- colnames(asv.tax)[-1]
for(taxlvl in 1:length(tax.levels)) {
  if(taxlvl == 1) {
    lvls_tmp <- asv.tax[,2]
  } else {
    lvls_tmp <- apply(asv.tax[,2:(taxlvl+1)],1,function(x) paste(x, collapse = ";"))
  }
  
  grps <- unique(lvls_tmp)
  mat_tmp <- matrix(0, ncol = ncol(asv.tab), nrow = length(grps))
  rownames(mat_tmp) <- grps
  colnames(mat_tmp) <- colnames(asv.tab)
  
  for(grp_i in grps) {
    if(sum(lvls_tmp == grp_i) > 1) {
      mat_tmp[grp_i,] <- colSums(asv.tab[lvls_tmp == grp_i,])
    } else {
      mat_tmp[grp_i,] <- asv.tab[lvls_tmp == grp_i,]
    }
  }
  tax.tab[[tax.levels[taxlvl]]] <- mat_tmp
}

# Alpha Diversity
cat("Calculating alpha diversities, eveness, and richness ...\n")
alpha_div <- data.table(sample       = colnames(asv.tab),
                        Shannon      = diversity(asv.tab, index = "shannon", MARGIN = 2),
                        Simpson      = diversity(asv.tab, index = "simpson", MARGIN = 2),
                        Simpson_inv  = diversity(asv.tab, index = "invsimpson", MARGIN = 2),
                        Fisher.alpha = fisher.alpha(asv.tab, MARGIN = 2),
                        S.richness   = colSums(asv.tab > 0))
alpha_div[, J.eveness := Shannon/log(S.richness)]

# Export processed data
cat("Exporting processed data tables ...\n")
write.table(asv.tab, file = "data/dada/asv_tab.tsv", sep = "\t", quote = F)
#write.table(asv.tab_rarefied, file = "data/dada/asv_tab_rarefied.tsv", sep = "\t", quote = F)
write.table(asv.tax, file = "data/dada/asv_tax.tsv", sep = "\t", quote = F)
writeXStringSet(asv_fasta, filepath = "data/dada/asv_seqs.fna")
fwrite(alpha_div, file = "data/dada/DT_alpha_diversity.tsv", sep = "\t", quote = F)
for(i in names(tax.tab)) {
  write.table(tax.tab[[i]], file = paste0("data/dada/tax_tab_",i,".tsv"), sep = "\t", quote = F)
}


