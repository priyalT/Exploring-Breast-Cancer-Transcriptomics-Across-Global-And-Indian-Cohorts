if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
library(dplyr)
library(magrittr)

query <- GDCquery(
  project = "TCGA-GBM", 
  data.category = "Clinical",
  data.type = "Clinical Supplement", 
  data.format = "BCR Biotab"
)

GDCdownload(query)

clinical.BCRtab.all <- GDCprepare(query)

names(clinical.BCRtab.all)

if (!requireNamespace("DT", quietly = TRUE)) {
  install.packages("DT")
}

library(DT)

clinical.BCRtab.all$clinical_patient_brca  %>% 
  head() %>% 
  DT::datatable(options = list(scrollX = TRUE, keys = TRUE))

data = clinical.BCRtab.all$clinical_patient_brca
#write.csv(data,"/Users/priyaltripathi/Documents/iit/TCGA-BRCA_clinical.csv", row.names = FALSE)


mol_subtypes <- TCGAquery_subtype(tumor = "brca")
  #write.csv(mol_subtypes,"/Users/priyaltripathi/Documents/iit/TCGA-BRCA_molsubtypes.csv", row.names = FALSE)


data <- read.csv("/Users/priyaltripathi/Documents/iit/TCGA-BRCA_molsubtypes.csv")
df <- data[1]
df <- df %>%
  mutate(age = data$age_at_initial_pathologic_diagnosis, stage = data$pathologic_stage, pathology = data$BRCA_Pathology,
         PAM50_subtype = data$BRCA_Subtype_PAM50)
missing_rows <- data[!data$patient %in% data2$bcr_patient_barcode, ]
data <- data %>% filter(!patient %in% missing_rows$patient)
df <- df %>% filter(!patient %in% missing_rows$patient)


#Genetic data
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
data <- GDCprepare(query)
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("SummarizedExperiment")
}
library(SummarizedExperiment)
assayNames(data)
fpkm <- assay(data, "fpkm_unstrand")
log_fpkm <- log2(fpkm + 1)
scaled_fpkm <- scale(log_fpkm)
gene_dataset <- t(scaled_fpkm)

#mapping
genes <- colnames(gene_dataset)
genes <- gsub("\\..*", "", genes)

library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                 filters = "ensembl_gene_id", 
                 values = genes, 
                 mart = mart)

ensembl_to_symbol <- setNames(mapping$hgnc_symbol, mapping$ensembl_gene_id)

# Step 4: Map the column names of gene_dataset using the Ensembl to HGNC symbol mapping
mapped_colnames <- sapply(genes, function(id) {
  if (id %in% names(ensembl_to_symbol)) {
    return(ensembl_to_symbol[id])  # Return the mapped HGNC symbol
  } else {
    return(id)  # Return Ensembl ID if no mapping found
  }
})

colnames(gene_dataset) <- mapped_colnames

