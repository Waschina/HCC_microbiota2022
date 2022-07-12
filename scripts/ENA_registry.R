#––––––––––––––––––––––––––––––––––#
# Generate sample registry for ENA #
#––––––––––––––––––––––––––––––––––#
library(data.table)

study_acc <- "PRJEB54571"

samples_meta <- fread("data/meta/samples_16S_DNA.csv")
samples_meta <- samples_meta[include == 1]

# taxonomies:
# gut: human feces metagenome (2705415)
# blood: human blood metagenome (1504969)
# liver: liver metagenome (1685930)

# broad-scale environmental context:
# all: An environmental system determined by an animal. (ENVO_01001002)

# local environmental context
# gut: feces (UBERON:0001555)
# blood: blood (UBERON:0000178)
# liver: liver (UBERON:0002107)


# environment medium
# gut: feces (UBERON:0001988)
# blood: blood (UBERON:0000178)
# liver: liver (UBERON:0002107)

# 
seq_files <- dir("data/fastq/DNA/")

#–––––––––––––––––#
# Sample registry #
#–––––––––––––––––#
splreg <- data.table(tax_id = NA_integer_,
                     scientific_name = NA_character_,
                     sample_alias = unique(gsub("\\.fastq\\.gz$",
                                                "",seq_files)))
splreg <- merge(splreg, samples_meta[, .(sample, Sample_type)],
                by.x = "sample_alias", by.y = "sample")
splreg <- splreg[,.(tax_id, scientific_name, sample_alias, Sample_type)]
splreg[, sample_title := sample_alias]
splreg[, sample_description := sample_alias]
splreg[, `project name` := study_acc]
splreg[, `sequencing method` := "MiSeq"]
splreg[, `collection date` := "not provided"]
splreg[, `geographic location (country and/or sea)` := "Austria"]
splreg[, `geographic location (latitude)` := 47.259659]
splreg[, `geographic location (longitude)` := 11.400375]

splreg[, `broad-scale environmental context` := "ENVO:01001002"]

splreg[Sample_type == "faeces", `local environmental context` := "UBERON:0001555"]
splreg[Sample_type == "blood", `local environmental context` := "UBERON:0000178"]
splreg[Sample_type == "liver", `local environmental context` := "UBERON:0002107"]

splreg[Sample_type == "faeces", `environmental medium` := "UBERON:0001988"]
splreg[Sample_type == "blood", `environmental medium` := "UBERON:0000178"]
splreg[Sample_type == "liver", `environmental medium` := "UBERON:0002107"]

splreg[Sample_type == "faeces", tax_id := 2705415]
splreg[Sample_type == "blood", tax_id := 1504969]
splreg[Sample_type == "liver", tax_id := 1685930]

splreg[Sample_type == "faeces", scientific_name := "human feces metagenome"]
splreg[Sample_type == "blood", scientific_name := "human blood metagenome"]
splreg[Sample_type == "liver", scientific_name := "liver metagenome"]

splreg[, Sample_type := NULL]

file.copy("output/ENA/Checklist_GSC-MIxS human associated_1657624790854.tsv",
          "output/ENA/Checklist_GSC-MIxS human associated_1657624790854_filled.tsv", overwrite = TRUE)

fwrite(splreg, file = "output/ENA/Checklist_GSC-MIxS human associated_1657624790854_filled.tsv",
       append = TRUE, sep = "\t", quote = FALSE)


#–––––––––––––––––––#
# PE-Reads registry #
#–––––––––––––––––––#
splreg <- fread("output/ENA/samples-2022-07-12T13 25 40.csv")

if(dir.exists("output/ENA/manifest_files/"))
  file.remove(dir("output/ENA/manifest_files/", full.names = TRUE))
dir.create("output/ENA/manifest_files/", showWarnings = FALSE)

# write manifest files for webin_cli
for(i in 1:nrow(splreg)) {
  olin <- rep(NA_character_, 8)
  olin[1] <- paste0("STUDY ", study_acc)
  olin[2] <- paste0("SAMPLE ", splreg[i, id])
  olin[3] <- paste0("NAME ", splreg[i, alias])
  olin[4] <- "INSTRUMENT Illumina MiSeq"
  olin[5] <- "LIBRARY_SOURCE METAGENOMIC"
  olin[6] <- "LIBRARY_SELECTION PCR"
  olin[7] <- "LIBRARY_STRATEGY AMPLICON"
  olin[8] <- paste0("FASTQ ",splreg[i, alias],".fastq.gz")
  
  writeLines(olin, paste0("output/ENA/manifest_files/", splreg[i, id],
                          ".txt"))
}


