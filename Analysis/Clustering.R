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

  
#trying out dea 
#tcga
tcga_gene <- t(tcga_gene_data)
tcga_umap_data$Cluster <- ifelse(tcga_df$U2 < -5, 'Cluster1', 'Cluster2')
group <- factor(tcga_umap_data$Cluster)
design <- model.matrix(~ 0 + group)
fit <- lmFit(tcga_gene, design)
contrast.matrix <- makeContrasts(group1_vs_group2 = groupCluster1 - groupCluster2, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit)
results <- topTable(fit2, adjust="fdr", number=Inf)
head(results)
write.csv(results, "DEA_results_TCGA.csv", row.names = TRUE)
tcga_umap_data$CFH <- as.numeric(tcga_gene["CFH", ])
tcga_umap_data$TMEM176A <- as.numeric(tcga_gene["TMEM176A", ])
tcga_umap_data$DBNDD1 <- as.numeric(tcga_gene["DBNDD1", ])
tcga_umap_data$COPZ2 <- as.numeric(tcga_gene["COPZ2", ])
tcga_umap_data$PRKAR2B <- as.numeric(tcga_gene["PRKAR2B", ])
tcga_umap_data$CROT <- as.numeric(tcga_gene["CROT", ])

 # Plot UMAP with FOXA1 expression
CFH <- ggplot(tcga_umap_data, aes(x = UMAP1, y = UMAP2, color = CFH, shape = Cluster)) +
  geom_point() +
  scale_color_viridis_c(option = "viridis")+
labs(title = "CFH Expression by Cluster", color = "CFH Expression") +
  theme_minimal()
CFH

TMEM176A <- ggplot(tcga_umap_data, aes(x = UMAP1, y = UMAP2, color = TMEM176A, shape = Cluster)) +
  geom_point() +
  scale_color_viridis_c(option = "magma")+
  labs(title = "TMEM176A Expression by Cluster", color = "TMEM176A Expression") +
  theme_minimal()
TMEM176A

DBNDD1 <- ggplot(tcga_umap_data, aes(x = UMAP1, y = UMAP2, color = DBNDD1, shape = Cluster)) +
  geom_point() +
  scale_color_viridis_c(option = "plasma")+
  labs(title = "DBNDD1 Expression by Cluster", color = "DBNDD1 Expression") +
  theme_minimal()
DBNDD1

COPZ2 <- ggplot(tcga_umap_data, aes(x = UMAP1, y = UMAP2, color = COPZ2, shape = Cluster)) +
  geom_point() +
  scale_color_viridis_c(option = "cividis") +
  labs(title = "COPZ2 Expression by Cluster", color = "COPZ2 Expression") +
  theme_minimal()
COPZ2

PRKAR2B <- ggplot(tcga_umap_data, aes(x = UMAP1, y = UMAP2, color = PRKAR2B, shape = Cluster)) +
  geom_point() +
  scale_color_distiller(palette = "Spectral") +
  labs(title = "PRKAR2B Expression by Cluster", color = "PRKAR2B Expression") +
  theme_minimal()
PRKAR2B

CROT <- ggplot(tcga_umap_data, aes(x = UMAP1, y = UMAP2, color = CROT, shape = Cluster)) +
  geom_point() +
  scale_color_viridis_c(option = "plasma")+
  labs(title = "CROT Expression by Cluster", color = "CROT Expression") +
  theme_minimal()
CROT

#indian
ind_gene <- t(ind_gene_data)
ind_umap_data$Cluster <- ifelse(ind_df$U2 > 0, 'Cluster1', 'Cluster2')
group <- factor(ind_umap_data$Cluster)
design <- model.matrix(~ 0 + group)
fit <- lmFit(ind_gene, design)
contrast.matrix <- makeContrasts(group1_vs_group2 = groupCluster1 - groupCluster2, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit)
results <- topTable(fit2, adjust="fdr", number=Inf)
head(results)
results["CROT", ]
sorted_results <- results[order(results$adj.P.Val), ]
head(sorted_results, 5)
write.csv(sorted_results, "DEA_results_ind.csv", row.names = TRUE)


ind_umap_data$KIT <- as.numeric(ind_gene["KIT", ])
ind_umap_data$SCARA5 <- as.numeric(ind_gene["SCARA5", ])
ind_umap_data$COL10A1 <- as.numeric(ind_gene["COL10A1", ])
ind_umap_data$DST <- as.numeric(ind_gene["DST", ])
ind_umap_data$SFRP1 <- as.numeric(ind_gene["SFRP1", ])
ind_umap_data$KCNJ16 <- as.numeric(ind_gene["KCNJ16", ])

KIT <- ggplot(ind_umap_data, aes(x = UMAP1, y = UMAP2, color = KIT, shape = Cluster)) +
  geom_point() +
  scale_color_viridis_c(option = "viridis")+
  labs(title = "KIT Expression by Cluster", color = "KIT Expression") +
  theme_minimal()
KIT

SCARA5 <- ggplot(ind_umap_data, aes(x = UMAP1, y = UMAP2, color = SCARA5, shape = Cluster)) +
  geom_point() +
  scale_color_viridis_c(option = "magma")+
  labs(title = "SCARA5 Expression by Cluster", color = "SCARA5 Expression") +
  theme_minimal()
SCARA5

COL10A1 <- ggplot(ind_umap_data, aes(x = UMAP1, y = UMAP2, color = COL10A1, shape = Cluster)) +
  geom_point() +
  scale_color_viridis_c(option = "plasma")+
  labs(title = "COL10A1 Expression by Cluster", color = "COL10A1 Expression") +
  theme_minimal()
COL10A1

DST <- ggplot(ind_umap_data, aes(x = UMAP1, y = UMAP2, color = DST, shape = Cluster)) +
  geom_point() +
  scale_color_viridis_c(option = "cividis") +
  labs(title = "DST Expression by Cluster", color = "DST Expression") +
  theme_minimal()
DST

SFRP1 <- ggplot(ind_umap_data, aes(x = UMAP1, y = UMAP2, color = SFRP1, shape = Cluster)) +
  geom_point() +
  scale_color_viridis_c(option = "cividis") +
  labs(title = "SFRP1 Expression by Cluster", color = "SFRP1 Expression") +
  theme_minimal()
SFRP1

KCNJ16 <- ggplot(ind_umap_data, aes(x = UMAP1, y = UMAP2, color = KCNJ16, shape = Cluster)) +
  geom_point() +
  scale_color_viridis_c(option = "cividis") +
  labs(title = "KCNJ16 Expression by Cluster", color = "KCNJ16 Expression") +
  theme_minimal()
KCNJ16

ind_umap_data$CROT <- as.numeric(ind_gene["CROT", ])
CROT <- ggplot(ind_umap_data, aes(x = UMAP1, y = UMAP2, color = CROT, shape = Cluster)) +
  geom_point() +
  scale_color_viridis_c(option = "cividis") +
  labs(title = "CROT Expression by Cluster", color = "CROT Expression") +
  theme_minimal()
CROT
