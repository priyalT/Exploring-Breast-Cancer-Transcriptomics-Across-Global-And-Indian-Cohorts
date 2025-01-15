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
ind_var <- arrange(ind_var, desc(ind_var))
ind_var <- rownames_to_column(ind_var, var = "gene_names")
ind_var <- ind_var[1:2000,]

#tcga
tcga_var <- apply(df_tcga, 2, var)
tcga_var <- as.data.frame(tcga_var)
tcga_var <- arrange(tcga_var, desc(tcga_var))
tcga_var <- rownames_to_column(tcga_var, var = "gene_names")
tcga_var <- tcga_var[1:2000,]

#subsetting the original dataset with the top 2000 variable genes
#indian
ind_gene_names <- ind_var$gene_names
selected_columns <- c("X", colnames(df_indian)[colnames(df_indian) %in% ind_gene_names])
ind_df <- df_indian[, selected_columns]

#tcga
tcga_gene_names <- tcga_var$gene_names
selected_columns <- c("X", colnames(df_tcga)[colnames(df_tcga) %in% tcga_gene_names])
tcga_df <- df_tcga[, selected_columns]

#including the patient characters in either both datasets
#indian
ind_char <- read.csv("/Users/priyaltripathi/iit/Indian/BRCA_Indian_patient_characteristics.csv")
ind_char <- ind_char[, c(4,5,6)]
ind_char <- as.data.frame(ind_char)
df_ind <- cbind(ind_df, ind_char)
write.csv(df_ind,"/Users/priyaltripathi/Multicluster Analysis/Analysis/Indian_BRCA.csv")

#tcga
tcga_char <- read.csv("/Users/priyaltripathi/Multicluster Analysis/Preprocess/TCGA/TCGA-BRCA_clinical.csv")
#figuring out which columns to include
tcga_col <- colnames(tcga_char)
class(tcga_col)
matches <- grepl("ihc", tcga_col, ignore.case = TRUE)
matching_columns <- tcga_col[matches]
matching_columns <- as.data.frame(matching_columns)
matching_columns_indices <- which(matches)
matching_columns_indices <- as.data.frame(matching_columns_indices)

#including those columns
filtered_tcga_char <- filtered_tcga_char[, c(1,142,148,154)]
filtered_tcga_char <- as.data.frame(filtered_tcga_char)
sum(duplicated(filtered_tcga_char$Sample_ID))
filtered_tcga_char <- filtered_tcga_char[!duplicated(filtered_tcga_char$Sample_ID), ]
selected_rows <- rownames(tcga_df) %in% rownames(filtered_tcga_char)
filtered_tcga <- tcga_df[selected_rows, ]

merged_df <- dplyr::left_join(filtered_tcga, filtered_tcga_char, by = c("X" = "Sample_ID"))
write.csv(merged_df,"/Users/priyaltripathi/Multicluster Analysis/Analysis/TCGA_BRCA.csv")

