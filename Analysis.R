#libraries
library(dplyr)
library(tibble)

#Reaading the dataset into variables
df_indian <- read.csv("/Users/priyaltripathi/iit/Intersection/Indian_intersect.csv")
df_tcga <- read.csv("/Users/priyaltripathi/iit/Intersection/TCGA-BRCA_intersect.csv")

#variance-based filtering
#indian
ind_var <- apply(df_indian, 2, var)
ind_var <- as.data.frame(ind_var)
ind_genes <- arrange(ind_var, desc(ind_var))
ind_genes <- rownames_to_column(ind_genes, var = "gene_names")
ind_genes <- ind_genes[1:2000,]

#tcga
tcga_var <- apply(df_tcga, 2, var)
tcga_var <- as.data.frame(tcga_var)
tcga_genes <- arrange(tcga_var, desc(tcga_var))
tcga_genes <- rownames_to_column(tcga_genes, var = "gene_names")
tcga_genes <- tcga_genes[1:2000,]

#UMAP
#subsetting the original dataset with the top 2000 variable genes
#indian
ind_gene_names <- ind_genes$gene_names
selected_columns <- c("X", colnames(df_indian)[colnames(df_indian) %in% ind_gene_names])
ind_df <- df_indian[, selected_columns]

#tcga
tcga_gene_names <- tcga_genes$gene_names
tcga_df <- df_tcga[, colnames(df_tcga) %in% tcga_gene_names]

#including the patient characters in either both datasets
#indian
ind_char <- read.csv("/Users/priyaltripathi/iit/Indian/BRCA_Indian_patient_characteristics.csv")
ind_char <- ind_char[, c(3,4,5)]
ind_char <- as.data.frame(ind_char)
df_ind <- cbind(ind_df, ind_char)

#tcga
tcga_char <- read.csv("/Users/priyaltripathi/iit/TCGA/TCGA-BRCA_clinical.csv")
tcga_char <- tcga_char[, c(44,50,56)]


