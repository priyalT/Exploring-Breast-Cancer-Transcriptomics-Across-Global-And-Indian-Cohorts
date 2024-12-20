#Intersection of gene_dataset and expression_data
a <- gene_dataset
b <- expression_data
tcga <- colnames(a)
indian <- colnames(b)
class(tcga)
intersection_col <- intersect(tcga,indian)
a <- a[, intersection_col]
b <- b[, intersection_col]
write.csv(a,"/Users/priyaltripathi/Documents/iit/TCGA-BRCA_intersect.csv", row.names = TRUE)
write.csv(b,"/Users/priyaltripathi/Documents/iit/Indian_intersect.csv", row.names = TRUE)

TC <- read.csv("TCGA-BRCA_intersect.csv")
save(TC, file = "TCGA-BRCA_intersect.RData")

India <- read.csv("Indian_intersect.csv")
save(India, file = "Indian_intersect.RData")
