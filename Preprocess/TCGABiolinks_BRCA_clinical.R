install.packages("janitor")

#load the libraries
library(TCGAbiolinks)
library(GEOquery)
library(tidyverse)
library(janitor)

#get samples for BRCA via GEO
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

#Load the patient metadata columns
na_strings <- c("NA", "[Unknown]", "[Not available]", "[Not Evaluated]", "[Not Applicable]")
Sys.setenv("VROOM_CONNECTION_SIZE" = 500000)
clinicalVariables <- getGEOSuppFiles("GSE62944", makeDirectory = F, baseDir = getwd(), filter_regex = "GSE62944_06_01_15_TCGA_24_548_Clinical_Variables_9264_Samples.txt.gz")
Clinical_Variables <- rownames(clinicalVariables) %>%
  read_tsv(col_names = F, na = na_strings) %>%
  dplyr::select(-(X2:X3))

#Just basically converting clinical_variables into a proper dataframe by transposing and doing other processing
Transposed_df <- as_tibble(t(Clinical_Variables), stringsAsFactors = F)
Transposed_df[1,1] <- "Sample_ID"
Transposed_df <- row_to_names(Transposed_df, 1, remove_row = TRUE, remove_rows_above = TRUE)
write.csv(Transposed_df,"/Users/priyaltripathi/Multicluster Analysis/Preprocess/TCGA/TCGA-BRCA_clinical.csv", row.names = FALSE)
