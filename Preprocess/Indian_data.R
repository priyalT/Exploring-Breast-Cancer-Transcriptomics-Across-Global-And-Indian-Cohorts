#Setting up the dataset
#libraries
library(GEOquery)
library(dplyr)
library (stringr)

#Procuring and cleaning data
gse <- getGEO("GSE40206", GSEMatrix = TRUE)
expression_data <- exprs(gse[[1]])
expression_data <- as.data.frame(expression_data)
gene_info <- fData(gse[[1]])
gene_symbols <- gene_info[[10]]
expression_data$gene_symbols <- gene_symbols
expression_data <- expression_data[,c(ncol(expression_data), 1:ncol(expression_data)-1)]
colnames(expression_data)[1] <- "Gene Symbols"
expression_data$`Gene Symbols`[expression_data$`Gene Symbols` == ""] <- NA
expression_data <- na.omit(expression_data)
expression_data <- t(expression_data)
expression_data <- as.data.frame(expression_data)
colnames(expression_data) <- expression_data[1,]
expression_data <- expression_data[-1,]
row_names <- rownames(expression_data)
View(expression_data)
expression_data <- lapply(expression_data, as.numeric)
sapply(expression_data, class)
range(expression_data)
expression_data<-as.data.frame(expression_data)
rownames(expression_data) <- row_names
#write.csv(expression_data,"/Users/priyaltripathi/Documents/iit/BRCA_Indian_Gene_data.csv")

#Setting up the characteristics dataset
tissue <- pData(gse[[1]])
characteristics <- tissue[1]
characteristics[c("Patient ID", "NOD Status", "ER status", "PR status", "HER status")] <- do.call(rbind, strsplit(characteristics$title, "_"))
characteristics <- characteristics %>% 
  select(-title)
characteristics <- characteristics %>%
  mutate(follow_up_period = tissue$characteristics_ch2.18, diagnosis = tissue$characteristics_ch2.11, related_investigation = tissue$characteristics_ch2.16,
         related_investigation_2 = tissue$characteristics_ch2.17, age = tissue$characteristics_ch2.9)
colnames(characteristics) <- c("Patient_ID", "NOD_status", "ER_status", "PR_status", "HER_status", "follow_up_period", "diagnosis", "related_investigation", "related_investigation_2", "age")
characteristics <- subset(characteristics, !(age == "size: 2.2x2x2.5" | age == "size: 5x4x2.5"))
#write.csv(characteristics,"/Users/priyaltripathi/Documents/iit/BRCA_Indian_patient_characteristics.csv")


#Making Table1
characteristics$diagnosis <- replace(characteristics$diagnosis, characteristics$diagnosis=="diag.: IDC-III", "IDC-III")
characteristics$diagnosis <- replace(characteristics$diagnosis, grepl("IDC.*DCIS", characteristics$diagnosis), "IDC-III with DCIS")
characteristics$diagnosis <- replace(characteristics$diagnosis, grepl("IDC.*infiltrative margins", characteristics$diagnosis), "IDC-III with Necrosis/Infiltrative Margins")
characteristics$diagnosis <- replace(characteristics$diagnosis, grepl("IDC.*necrosis", characteristics$diagnosis), "IDC-III with Necrosis/Infiltrative Margins")
characteristics$diagnosis <- replace(characteristics$diagnosis, grepl("IDC.*infiltrating margins", characteristics$diagnosis), "IDC-III with Necrosis/Infiltrative Margins")
characteristics$diagnosis <- replace(characteristics$diagnosis, 
                                     grepl("DCIS", characteristics$diagnosis) & 
                                       !grepl("IDC-III with DCIS", characteristics$diagnosis), 
                                     "DCIS")
characteristics$diagnosis <- replace(characteristics$diagnosis, 
                                     !characteristics$diagnosis %in% c("IDC-III", "IDC-III with DCIS", 
                                                                       "IDC-III with Necrosis/Infiltrative Margins", 
                                                                       "DCIS"), 
                                     "Other")
characteristics$age <- gsub("age: ", "", characteristics$age)
characteristics$age <- as.numeric(characteristics$age)
characteristics$age_group <- cut(characteristics$age, breaks = seq(20, 90, by = 10), right = FALSE)

characteristics$related_investigation <- replace(characteristics$related_investigation, grepl("Passed away | other relevant investigations: -------- | follow-up period: 3yr.5m", characteristics$related_investigation), "Deceased")
characteristics$related_investigation <- replace(characteristics$related_investigation, grepl("didn't | refused | not willing | didnt | No follow up", characteristics$related_investigation), "Lost to follow-up")
characteristics$related_investigation <- replace(characteristics$related_investigation, 
                                     !characteristics$related_investigation %in% c("Deceased", "Lost to follow-up"), 
                                     "Alive")


characteristics$status <- 
  factor(characteristics$diagnosis,
         levels = c("IDC-III", "IDC-III with DCIS", "IDC-III with Necrosis/Infiltrative Margins", "DCIS", "Other"),
         labels = c("IDC-III", "IDC-III with DCIS", "IDC-III with Necrosis/Infiltrative Margins", "DCIS", "Other"))
characteristics$investigation <- 
  factor(characteristics$related_investigation,
         levels = c("Alive", "Deceased", "Lost to follow-up"),
         labels = c("Alive", "Deceased", "Lost to follow-up"))


library(table1)
table1(~ age_group + factor(NOD_status) + factor(ER_status) + factor(PR_status) * factor(HER_status) | factor(status) + factor(investigation), data = characteristics)

