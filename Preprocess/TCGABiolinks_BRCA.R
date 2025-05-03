#libraries
library(TCGAbiolinks)
library(dplyr)
library(magrittr)

#quering the TCGA database for clinical data related to TCGA-BRCA
query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Clinical",
  data.type = "Clinical Supplement", 
  data.format = "BCR Biotab"
)

#downloading the data
GDCdownload(query)

#preparing the data
clinical.BCRtab.all <- GDCprepare(query)

#displaying the names of the datasets available
names(clinical.BCRtab.all)

#installing the DT package if not already installed
if (!requireNamespace("DT", quietly = TRUE)) {
  install.packages("DT")
}

#loading the DT package
library(DT)

# Display the first few rows of the clinical_drug_acc dataset in an interactive table
clinical.BCRtab.all$clinical_patient_brca  %>% 
  head() %>% 
  DT::datatable(options = list(scrollX = TRUE, keys = TRUE))

data = clinical.BCRtab.all$clinical_patient_brca
write.csv(data,"/Users/priyaltripathi/Documents/iit/TCGA-BRCA_clinical.csv", row.names = FALSE)

data2 <- read.csv("/Users/priyaltripathi/Documents/iit/TCGA-BRCA_clinical.csv")
data2 <- data2[data2$bcr_patient_barcode %in% data$patient, ]
df <- df %>%
  mutate(gender = data2$gender, Nodal_status = data2$ajcc_nodes_pathologic_pn, ER_status = data2$er_status_by_ihc, PR_status = data2$pr_status_by_ihc,
         HER_status = data2$her2_status_by_ihc, histological_type = data2$histological_type)

write.csv(data,"/Users/priyaltripathi/Documents/iit/TCGA-BRCA_patient_characteristics.csv", row.names = FALSE)
