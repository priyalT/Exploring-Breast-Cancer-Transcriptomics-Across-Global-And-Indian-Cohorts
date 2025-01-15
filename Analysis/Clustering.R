#installing the necessary libraries
install.packages("umap")
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("patchwork")
install.packages("cluster")

#loading the libraries
library(umap)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(patchwork)
library(cluster)

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
tcga_gene_data <- tcga_df[, -c(1, 2, 2003, 2004, 2005)]
tcga_umap_result <- umap(tcga_gene_data)
tcga_umap_data <- as.data.frame(tcga_umap_result$layout)

colnames(tcga_umap_data) <- c("UMAP1", "UMAP2")
tcga_umap_data$Sample_ID <- tcga_df$X
tcga_umap_data$ER <- tcga_df$er_status_by_ihc
tcga_umap_data$PR <- tcga_df$pr_status_by_ihc
tcga_umap_data$HER2 <- tcga_df$her2_status_by_ihc

plot_er <- ggplot(tcga_umap_data, aes(x = UMAP1, y = UMAP2, color = as.factor(ER))) +
  geom_point() +
  labs(title = "ER Status", x = "UMAP1", y = "UMAP2") +
  theme_minimal()

plot_pr <- ggplot(tcga_umap_data, aes(x = UMAP1, y = UMAP2, color = as.factor(PR))) +
  geom_point() +
  labs(title = "PR Status", x = "UMAP1", y = "UMAP2") +
  theme_minimal()

plot_her2 <- ggplot(tcga_umap_data, aes(x = UMAP1, y = UMAP2, color = as.factor(HER2))) +
  geom_point() +
  labs(title = "HER2 Status", x = "UMAP1", y = "UMAP2") +
  theme_minimal()

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

plot_er_ind <- ggplot(ind_umap_data, aes(x = UMAP1, y = UMAP2, color = as.factor(ER))) +
  geom_point() +
  labs(title = "ER Status", x = "UMAP1", y = "UMAP2") +
  theme_minimal()

plot_pr_ind <- ggplot(ind_umap_data, aes(x = UMAP1, y = UMAP2, color = as.factor(PR))) +
  geom_point() +
  labs(title = "PR Status", x = "UMAP1", y = "UMAP2") +
  theme_minimal()

plot_her2_ind <- ggplot(ind_umap_data, aes(x = UMAP1, y = UMAP2, color = as.factor(HER2))) +
  geom_point() +
  labs(title = "HER2 Status", x = "UMAP1", y = "UMAP2") +
  theme_minimal()

grid.arrange(plot_er_ind, plot_pr_ind, plot_her2_ind, ncol = 3)

plot_er
plot_er_ind

plot_pr
plot_pr_ind

plot_her2
plot_her2_ind


# Assuming `ind_umap_data` contains your UMAP coordinates and ER status
umap_coords <- ind_umap_data[, c("UMAP1", "UMAP2")]  # Extract UMAP coordinates
er_status <- as.factor(ind_umap_data$ER)  # ER status as a factor (cluster labels)

# Calculate silhouette scores
sil <- silhouette(as.numeric(er_status), dist(umap_coords))

# Visualize average silhouette width
plot(sil, main = "Silhouette Plot")

# Average silhouette score
avg_sil <- mean(sil[, "sil_width"])
print(paste("Average Silhouette Score:", avg_sil))
