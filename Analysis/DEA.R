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

# Plot UMAP
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
