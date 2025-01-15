install.packages("janitor")

library(TCGAbiolinks)
library(GEOquery)
library(tidyverse)
library(janitor)

cancerTypesamples <- getGEOSuppFiles("GSE62944", makeDirectory = F, baseDir = getwd(), filter_regex = "GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt.gz")
CancerType <- rownames(cancerTypesamples) %>%
  read_tsv(col_names = F)
colnames(CancerType) <- c("Sample_ID", "cancer_type")
CancerType <- CancerType %>%
  dplyr::filter(cancer_type == "BRCA")

normalTypesamples <- getGEOSuppFiles("GSE62944", makeDirectory = F, baseDir = getwd(), filter_regex = "GSE62944_06_01_15_TCGA_24_Normal_CancerType_Samples.txt.gz")
NormalType <- rownames(normalTypesamples) %>%
  read_tsv(col_names = F)
colnames(NormalType) <- c("Sample_ID", "cancer_type")
NormalType <- NormalType %>%
  dplyr::filter(cancer_type == "BRCA")%>%
  mutate(bcr_patient_barcode = substr(Sample_ID, 1, 12))

na_strings <- c("NA", "[Unknown]", "[Not available]", "[Not Evaluated]", "[Not Applicable]")
Sys.setenv("VROOM_CONNECTION_SIZE" = 500000)
clinicalVariables <- getGEOSuppFiles("GSE62944", makeDirectory = F, baseDir = getwd(), filter_regex = "GSE62944_06_01_15_TCGA_24_548_Clinical_Variables_9264_Samples.txt.gz")
Clinical_Variables <- rownames(clinicalVariables) %>%
  read_tsv(col_names = F, na = na_strings) %>%
  dplyr::select(-(X2:X3))

Transposed_df <- as_tibble(t(Clinical_Variables), stringsAsFactors = F)
Transposed_df[1,1] <- "Sample_ID"
Transposed_df <- row_to_names(Transposed_df, 1, remove_row = TRUE, remove_rows_above = TRUE)

Merged_tumor_df <- Transposed_df %>%
  inner_join(CancerType, by = "Sample_ID") %>%
  dplyr::select(-cancer_type)


Merged_normal_df <- Transposed_df %>%
  inner_join(NormalType, by ="bcr_patient_barcode") %>%
  rename(Sample_ID = Sample_ID.x) %>%
  dplyr::select(-c(Sample_ID.y, cancer_type))


query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Clinical",
  data.type = "Clinical Supplement", 
  data.format = "BCR Biotab"
)


#Download the data
GDCdownload(query)

#Prepare data
clinical.BCRtab.all <- GDCprepare(query)

#List the names in clinical.BCRtab.all
names(clinical.BCRtab.all)

#Install DT library and load it
if (!requireNamespace("DT", quietly = TRUE)) {
  install.packages("DT")
}

library(DT)

#Load the head() or top 6 rows of clinical_patient_brca which is the clinical information
clinical.BCRtab.all$clinical_patient_brca  %>% 
  head() %>% 
  DT::datatable(options = list(scrollX = TRUE, keys = TRUE))

#query molecular subtypes and store them under mol_subtypes and write it into a csv
mol_subtypes <- TCGAquery_subtype(tumor = "brca")
write.csv(mol_subtypes,"/Users/priyaltripathi/Documents/iit/TCGA-BRCA_molsubtypes.csv", row.names = FALSE)

#Manipulating the molecular subtypes data
data <- read.csv("/Users/priyaltripathi/Documents/iit/TCGA-BRCA_molsubtypes.csv")
df <- data[1]
df <- df %>%
  mutate(age = data$age_at_initial_pathologic_diagnosis, stage = data$pathologic_stage, pathology = data$BRCA_Pathology,
         PAM50_subtype = data$BRCA_Subtype_PAM50)
missing_rows <- data[!data$patient %in% data2$bcr_patient_barcode, ]
data <- data %>% filter(!patient %in% missing_rows$patient)
df <- df %>% filter(!patient %in% missing_rows$patient)




#Genetic data curation
#Query the genetic data
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
#Download it
GDCdownload(query)
#store it in a data variable
data <- GDCprepare(query)

#If SummarizedExperiment is not present install it, and load it
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("SummarizedExperiment")
}
library(SummarizedExperiment)

#Obtain the names of assays or data objects associated with an object, typically a SummarizedExperiment object in R.
assayNames(data)
#Extract the “fpkm_unstrand” assay from the data object and store it in fpkm variable
fpkm <- assay(data, "fpkm_unstrand")
#Apply a log2 transformation to the FPKM values.
log_fpkm <- log2(fpkm + 1)
#Scale the log-transformed FPKM values so that each gene has a mean of 0 and standard deviation of 1
scaled_fpkm <- scale(log_fpkm)
#Transpose the dataset to have genes as colnames and samples as rownames
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

#Map the column names of gene_dataset using the Ensembl to HGNC symbol mapping
mapped_colnames <- sapply(genes, function(id) {
  if (id %in% names(ensembl_to_symbol)) {
    return(ensembl_to_symbol[id])  # Return the mapped HGNC symbol
  } else {
    return(id)  # Return Ensembl ID if no mapping found
  }
})

colnames(gene_dataset) <- mapped_colnames

