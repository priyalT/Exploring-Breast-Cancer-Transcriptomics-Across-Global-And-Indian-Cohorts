#installing the necessary libraries
install.packages("umap")
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("patchwork")
install.packages("cluster")
install.packages("BiocManager")
BiocManager::install("edgeR")


#loading the libraries
library(umap)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(patchwork)
library(cluster)
library(limma)
library(edgeR)


#reading the tcga dataframe, converting positive negative values to binary, and omitting NA values
tcga_df <- read.csv("/Users/priyaltripathi/Multicluster Analysis/Analysis/TCGA_BRCA.csv")
tcga_df <- tcga_df %>%
  mutate(
    er_status_by_ihc = ifelse(er_status_by_ihc == "Positive", 1, 0),
    pr_status_by_ihc = ifelse(pr_status_by_ihc == "Positive", 1, 0),
    her2_status_by_ihc = ifelse(her2_status_by_ihc == "Positive", 1, 0)
  )
tcga_df <- tcga_df %>%
  na.omit()

#UMAP clustering
#tcga
tcga_gene_data <- tcga_df[, -c(1, 2, 2003, 2004, 2005)]
tcga_umap_result <- umap(tcga_gene_data)
tcga_umap_data <- as.data.frame(tcga_umap_result$layout)
colnames(tcga_umap_data) <- c("UMAP1", "UMAP2")
tcga_umap_data$Sample_ID <- tcga_df$X
tcga_umap_data$ER <- tcga_df$er_status_by_ihc
tcga_umap_data$PR <- tcga_df$pr_status_by_ihc
tcga_umap_data$HER2 <- tcga_df$her2_status_by_ihc

tcga_df$U1 <- tcga_umap_data$UMAP1
tcga_df$U2 <- tcga_umap_data$UMAP2
tcga_df$coarse.cluster <- ifelse(tcga_df$U2 < -5, 'Cluster1', 'Cluster2')

plot_er <- ggplot(tcga_umap_data, aes(x = UMAP1, y = UMAP2, color = as.factor(ER))) +
  geom_point() +
  labs(title = "ER Status", x = "UMAP1", y = "UMAP2") +
  theme_minimal()
plot_er

plot_pr <- ggplot(tcga_umap_data, aes(x = UMAP1, y = UMAP2, color = as.factor(PR))) +
  geom_point() +
  labs(title = "PR Status", x = "UMAP1", y = "UMAP2") +
  theme_minimal()
plot_pr

plot_her2 <- ggplot(tcga_umap_data, aes(x = UMAP1, y = UMAP2, color = as.factor(HER2))) +
  geom_point() +
  labs(title = "HER2 Status", x = "UMAP1", y = "UMAP2") +
  theme_minimal()
plot_her2

grid.arrange(plot_er, plot_pr, plot_her2, ncol = 3)


#reading the indian dataframe, converting positive negative values to binary, and omitting NA values
ind_df <- read.csv("/Users/priyaltripathi/Multicluster Analysis/Analysis/Indian_BRCA.csv")
ind_df <- ind_df %>%
  mutate(
    ER.status = ifelse(ER.status == "ER+", 1, 0),
    PR.status = ifelse(PR.status == "PR+", 1, 0),
    HER.status = ifelse(HER.status == "HER+", 1, 0)
  )
ind_df <- ind_df %>%
  na.omit()
ind_gene_data <- ind_df[, -c(1, 2, 2003, 2004, 2005)]
ind_umap_result <- umap(ind_gene_data)
ind_umap_data <- as.data.frame(ind_umap_result$layout)
colnames(ind_umap_data) <- c("UMAP1", "UMAP2")
ind_umap_data$Sample_ID <- ind_df$X
ind_umap_data$ER <- ind_df$ER.status
ind_umap_data$PR <- ind_df$PR.status
ind_umap_data$HER2 <- ind_df$HER.status
ind_df$U1 <- ind_umap_data$UMAP1
ind_df$U2 <- ind_umap_data$UMAP2
ind_df$coarse.cluster <- ifelse(ind_df$U2 > 0, 'Cluster1', 'Cluster2')

plot_er_ind <- ggplot(ind_umap_data, aes(x = UMAP1, y = UMAP2, color = as.factor(ER))) +
  geom_point() +
  labs(title = "ER Status", x = "UMAP1", y = "UMAP2") +
  theme_minimal()
plot_er_ind

plot_pr_ind <- ggplot(ind_umap_data, aes(x = UMAP1, y = UMAP2, color = as.factor(PR))) +
  geom_point() +
  labs(title = "PR Status", x = "UMAP1", y = "UMAP2") +
  theme_minimal()
plot_pr_ind

plot_her2_ind <- ggplot(ind_umap_data, aes(x = UMAP1, y = UMAP2, color = as.factor(HER2))) +
  geom_point() +
  labs(title = "HER2 Status", x = "UMAP1", y = "UMAP2") +
  theme_minimal()
plot_her2_ind

grid.arrange(plot_er_ind, plot_pr_ind, plot_her2_ind, ncol = 3)

plot_er
plot_er_ind

plot_pr
plot_pr_ind

plot_her2
plot_her2_ind

#checking which metadata is significantly associated w the clusters
table_cluster_her <- table(tcga_df$coarse.cluster, tcga_df$her2_status_by_ihc)
chisq.test(table_cluster_her)

table_cluster_her <- table(ind_df$coarse.cluster, ind_df$HER.status)
chisq.test(table_cluster_her)

table_cluster_pr <- table(tcga_df$coarse.cluster, tcga_df$pr_status_by_ihc)
chisq.test(table_cluster_pr)

table_cluster_pr <- table(ind_df$coarse.cluster, ind_df$PR.status)
chisq.test(table_cluster_pr)

table_cluster_er <- table(tcga_df$coarse.cluster, tcga_df$er_status_by_ihc)
chisq.test(table_cluster_er)

table_cluster_er <- table(ind_df$coarse.cluster, ind_df$ER.status)
chisq.test(table_cluster_er)
