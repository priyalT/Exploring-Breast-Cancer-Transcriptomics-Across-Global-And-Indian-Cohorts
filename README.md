# Exploring Breast Cancer Transcriptomics Across Global And Indian Cohorts: A Comparative Study

1. [Introduction](https://github.com/priyalT/Multicluster-and-Survival-Analysis/edit/main/README.md#introduction)
2. [Processing, Clustering and Differential Expression Analysis of Indian and Global Breast Cancer Transcriptomics data](https://github.com/priyalT/Multicluster-and-Survival-Analysis/edit/main/README.md#processing-clustering-and-differential-expression-analysis-of-indian-and-global-breast-cancer-transcriptomics-data)
3. [Understanding the significance of UMAP clusters and Differential Expression Analysis of both cohorts](https://github.com/priyalT/Multicluster-and-Survival-Analysis/edit/main/README.md#understanding-the-significance-of-umap-clusters-and-differential-expression-analysis-of-both-cohorts)
4. [Cross-population analysis of Breast Cancer Transcriptomics](https://github.com/priyalT/Multicluster-and-Survival-Analysis/edit/main/README.md#cross-population-analysis-of-breast-cancer-transcriptomics)
5. [Conclusion](https://github.com/priyalT/Multicluster-and-Survival-Analysis/edit/main/README.md#conclusion)
6. [References](https://github.com/priyalT/Multicluster-and-Survival-Analysis/edit/main/README.md#references)



Introduction 
=======
Heterogeneity is seen widely in breast cancer, which is also one of the leading causes of malignancy-related deaths in women [1]. According to Global Cancer Observatory (GLOBOCAN) 2022, 2.3 million new cases and 67,000 deaths had been reported from female breast cancer [2]. In the same report, India ranked third in terms of cancer diagnoses and second in terms of cancer-related mortality [3]. Female breast cancer has been reported to be the most prevalent. Hence, it becomes imperative to study breast cancer because of its variation in epidemiology and differences in its molecular subtypes. 
A report from Indian Journal of Cancer analysed the differences in rural and urban breast cancers in India [4]. The study concluded that the lifestyle in rural areas is protective against risk of developing breast cancer. Extrapolating these results, it is safe to say that there may also be differences seen between breast cancer between India and the rest of the world. For example, TP53 mutations are more pronounced in the Indian breast cancer population (approximately 55% of the patients carry this mutation) compared to other populations worldwide [5]. While these differences exist and have aided in prediction, diagnosis, drug design and a lot more, there is lack thereof the similarities that occur between global and Indian breast cancer patients. 
Comparative analysis of global and Indian breast cancer patients can lead to better population-specific oncological research, as seen in the study led by Paul et. al [6]. The study focused on oral squamous cell-carcinoma and identified biomarkers and mutations found in both global and Indian populations. Similarly, in yet another study, it was found that in non-small cell lung cancer (NSCLC), mutations in EGFR and ALK are prevalent in both Indian and global populations [7], [8].
There is an evident lack in comparative oncology research between India and the global circumstance. Considering the number of studies that have been carried out focusing on biomarkers, shared and different, among the Indian and global cohort, it becomes absolutely imperative for intended and focussed research in this direction. We aim to bridge this lack, by proposing a study which not only aims to identify recurring and different clustering or biomarker patterns in breast cancer datasets from India and the global TCGA cohort, highlighting similarities that may indicate shared biological significance, but also investigate unique cluster-specific gene expression patterns, which may help characterize distinct molecular subtypes.

![image](https://github.com/user-attachments/assets/13be807d-424f-4987-a737-2a2c94e2b78b)


Processing, Clustering and Differential Expression Analysis of Indian and Global Breast Cancer Transcriptomics data
=========

All analyses were conducted in R (version 4.4.0) using packages as mentioned per subsection.


* Data curation and processing of both Indian and Global datasets:

  The Indian data was obtained via the respective GEO website (GSE40206) [9], whereas TCGA data was obtained using getGEOSuppFiles() from R (GSE62944) [10], [11]. Different forms of the same cohort's data were obtained, including gene expression, molecular subtypes, and clinical characteristics. The respective patient metadata was also acquired. NA values were omitted from both the datasets. The Indian dataset was already normalized. To prepare TCGA's genetic expression data for analysis, the data had to be log-normalized and scaled. Transposing the dataframe was also required. All of this was achieved with TCGABiolinks and SummaryExperiments. EMBL IDs were later mapped from the TCGA dataset to HGNC symbols, which were also seen in the Indian dataset. This prepared the TCGA dataset for further comparative analysis. Using intersect(), both TCGA and Indian dataset were subsetted to include common genes. This solidified both our datasets for comparison. The summary of both the datasets is mentioned in Table 1.

   | Dataset | Samples      | Total Genes |
   | ------------- | ------------- |----------|
   | TCGA | 815 |15243 |
   | Indian (GSE40206) | 81 | 15243 |

   Table 1. Summary of datasets used

* Variance-based gene filtering:
   To both the datasets, the statistical method of variance was applied to genes to filter them. Using apply() [12], sample variance for both the datasets was calculated. Top 2000 genes were selected based on the variance score and the original datasets were then subsetted to only include these top 2000 genes. Further processing included adding a ‘patient_characteristics’ column to the dataframe which detailed the diagnosis of each patient. 

* Clustering analysis:
   Once the final dataset was prepared with variance filtered genes, UMAP clustering was performed. Clusters were manually assigned as they could be visually segregated. No algorithm was applied. Segregation was done via Cluster1 and Cluster2. The UMAP plots were graphed in terms of ER status, PR status and HER status. The two values of the statuses were set to 0 (absent) and 1 (present). To check how significant the ER, PR and HER status is in both the cohorts, chi square tests were performed for all three in both the cohorts, and the results obtained are mentioned in Table 2. 

* Differential Expression Analysis (DEA):
   Differential expression analysis was further performed to identify the statistically significant genes present in both datasets that have a correlation with the type of breast cancer. Multiple hypothesis testing was also performed on the result to identify the most significant genes. While the processed datasets shared many common genes, the sets of differentially expressed genes (DEGs) identified in each cohort were largely distinct. This suggests that although the same genes were analyzed, their expression patterns and statistical significance varied across the two populations, potentially reflecting underlying biological or environmental differences. The significant genes found are included in Table 3.

Understanding the significance of UMAP clusters and Differential Expression Analysis of both cohorts
==========
* Results of clustering
  -------------
   The final results of chi-square tests validated the clustering observed via UMAP. The tests showed that the association between cluster assignments and ER, PR, HER2 status was strong and statistically significant due to the good p-values (less than 0.05) obtained. This strengthened the results from the UMAP plots performed, with the factor being presence or absence of estrogen receptor, progesterone receptor, and human epidermal growth factor 2. Clear clustering was observed in both cohorts. However, subtype separation (e.g., ER+ vs ER−) appeared more distinct in the TCGA dataset, possibly due to larger sample size.

   | Dataset | Gene | chi-square | p-value | degree of freedom |
   | ------------- | ------------- | ---------- | ---------- | ------ |
   | Indian | HER | 9.512 | 0.002 | 1 |
   |        | PR | 8.4301 | 0.003 | 1 |
   |        | ER | 12.188 | 0.0004 | 1 |
   | Global | HER | 14.899 | 0.002 | 1 |
   |        | PR | 285.221 | < 2.2e-16 | 1 |
   |        | ER | 439.96 | < 2.2e-16 | 1 |

     Table 2. Summary of statistical results of clustering analysis
  
   ![image](https://github.com/user-attachments/assets/effbdaff-9773-497e-9e09-35e907736852)
   ![image](https://github.com/user-attachments/assets/41bcb10e-adf5-4501-9ef7-b8ca96e1cc72)
   ![image](https://github.com/user-attachments/assets/b89622f6-a553-4f63-b019-aa862ec18749)
     Figure 1. UMAP plots of HER2, PR, ER status of the global cohort, respectively. ‘0’ represents absence and ‘1’ represents presence of the receptor

  ![image](https://github.com/user-attachments/assets/cd3f9433-4782-4ee0-a889-9c90bf649707)
  ![image](https://github.com/user-attachments/assets/f5a0924d-6d40-4c7b-ad66-9cc5d65ee44d)
  ![ind_UMAP_pr](https://github.com/user-attachments/assets/722d7950-f2ac-44d9-8069-eca8b2398e48)
     Figure 2. UMAP plots of HER2, PR, ER status of the indian cohort, respectively. ‘0’ represents absence and ‘1’ represents presence of the receptor

* Results of Differential Expression Analysis (DEA)
  --------------------
  DEA identified several notable genes in the cohorts. The list differed in the order of the genes in both the cohorts. This gives credence to the fact that, although receptor-based clustering is similar, the gene expression underlying those subtypes can be different in Indian and worldwide cohorts. Table 3 presents the top five most notable genes identified to be differentially expressed in both datasets. Interestingly, cell signaling and tumorigenesis-related genes KIT and SFRP1 were found in the Indian cohort and are breast cancer-associated genes [13], [14]. The presence of common genes such as CFH, DST, and COPZ2 in the sets can be indicative of core molecular characteristics common to breast cancer tumors independent of population origin.

  These common DEGs can also be used as potential universal biomarkers, and uniquely expressed genes in each of the cohorts identify population-specific molecular profiles. Such duality emphasizes the need for population-conscious methodologies in cancer genomics, in which both common and distinct attributes guide diagnosis, prognosis, and treatment strategies.

   | Global | Indian | 
   | ------ | ------| 
   | CFH | KIT | 
   | TMEM176A | SCARA5 | 
   | DBNDD1 | COL10A1 | 
   | COPZ2 | DST | 
   | PRKAR2B | SFRP1 |

  Table 3. Most significant genes present in both global and indian cohorts.

  ![image](https://github.com/user-attachments/assets/76b43eec-97ed-4f20-adde-7e0b6128bc53)
  ![image](https://github.com/user-attachments/assets/c470bdae-afdf-4b24-bb19-3afe01e4b41c)
  ![image](https://github.com/user-attachments/assets/4b1236a5-2e4e-4a8e-95d9-26d085c18727)
  ![image](https://github.com/user-attachments/assets/1f5d070f-dda7-4c01-b38f-2f2fcdfc2cf9)
  ![image](https://github.com/user-attachments/assets/da236df5-a42c-49a4-afbc-630fa67cb58a)

  *** Figure 3. UMAP plots of CFH, TMEM176A, DBNDD1, COPZ2, PRKAR2B expression of the global cohort, respectively. Color grading has been performed for the expression values. ***

  ![image](https://github.com/user-attachments/assets/2d03b595-cd83-40a6-830a-f0e56819e601)
  ![image](https://github.com/user-attachments/assets/d8b04393-c908-446d-adf0-98bd04855c5b)
  ![image](https://github.com/user-attachments/assets/c3f2ac18-a68f-40d1-974d-ca216facd83b)
  ![image](https://github.com/user-attachments/assets/b1b378fb-c14a-41dc-8b7f-99d801843c5e)
  ![image](https://github.com/user-attachments/assets/5bbe81df-8e3a-43b1-82c2-85448182dd0c)

  Figure 4. UMAP plots of KIT, SCARA5, COL10A1, DST, SFRP1 expression of the Indian cohort, respectively. Color grading has been performed for the expression values.


Cross-population analysis of Breast Cancer Transcriptomics
============================================================
This research sought to contrast and also find similar gene expression profiles and molecular subtypes of breast cancer in Indian and global (TCGA) cohorts, with emphasis on common and unique clustering patterns and differentially expressed genes (DEGs). Our findings present evidence that, in spite of population-specific variations, some molecular characteristics of breast cancer are geographically conserved.
UMAP-based clustering identified and evidenced unique molecular subtypes in both cohorts, highly correlated with ER, PR, and HER2 status. The statistical significance of these clusters, confirmed by chi-square tests, highlights the importance of these receptors as major classifiers of breast cancer subtypes worldwide. The more distinct separation of subtypes in the TCGA dataset is probably due to the greater sample size, which increases resolution and statistical power. Nevertheless, the occurrence of distinct clusters within the small Indian cohort suggests that receptor status-based molecular subtyping is resilient across populations.
To further recognise the unique molecular subtypes, DEGs were found to have substantial divergence between the two populations based on differential gene expression analysis, even though they had analyzed the same group of 2000 high-variance genes. This discovery indicates that although the general molecular architecture of the disease is similar, local environmental, genetic, or epigenetic determinants could modulate gene expression patterns in population-specific manners. This is consistent with earlier findings indicating variability in mutation frequencies—such as the Indian breast cancer patient population's high frequency of TP53 mutations—reflecting population-specific oncogenic processes.
Notably, a group of genes, such as CFH, DST, and COPZ2, were differentially expressed in both cohorts, albeit not in the same priority. These genes can be considered a part of central elements of breast tumor biology and are of interest as universal biomarkers or targets for therapy. On the other hand, genes such as KIT and SFRP1, which were mostly present in the Indian dataset, could indicate regional specificity, highlighting the value in incorporating diverse groups of people in genomic analysis. These findings of overlapping biomarkers present opportunities for broad-spectrum diagnostic and therapeutic application, whereas population-specific observations inform the establishment of localized, precision medicine strategies. 

Limitations:
--------------
While yielding promising findings, this investigation has some limitations. The Indian sample set was relatively small in size (n=81), potentially limiting the generalizability of certain findings and diminishing statistical power. The second limitation is heterogeneity of platforms and preprocessing pipelines between GEO and TCGA datasets, which, although effort has been put into normalization, can cause batch effects or noise. More inclusion of Indian datasets or utilization of single-platform data could alleviate this.

Future Directions: 
--------------------
Subsequent studies must build on this comparative platform by integrating multi-omics information (e.g., proteomics, methylomics), clinical endpoints, and survival data. Longitudinal samples from Indian patients might be used to confirm the prognostic utility of identified biomarkers. Lastly, functional work must explore the roles of common DEGs such as DST and CFH in breast cancer biology.


Conclusion
=======
The current research offers comparative gene expression profile analysis of Indian and global (TCGA) breast cancer cohorts. In contrast to variations in sample size and population-specific genetic backgrounds, identical clustering according to ER, PR, and HER2 status was noted in both cohorts. Although the majority of DEGs were population-specific, a small subset of common DEGs could possibly reflect conserved molecular characteristics of breast cancer. These results underscore the need to enroll underrepresented populations in cancer genomics studies to identify both universal and region-specific biomarkers. Subsequent research with larger, more diverse Indian patient series and inclusion of further clinical and molecular information is necessary to better understand and inform more inclusive precision oncology approaches.

References
===========
1. Xiong, X., Zheng, L., Ding, Y., Chen, Y., Cai, Y., Wang, L., Huang, L., Liu, C., Shao, Z., & Yu, K. (2025). Breast cancer: pathogenesis and treatments. Signal Transduction and Targeted Therapy, 10(1). https://doi.org/10.1038/s41392-024-02108-4
2. Kim, J., Harper, A., McCormack, V., Sung, H., Houssami, N., Morgan, E., Mutebi, M., Garvey, G., Soerjomataram, I., & Fidler-Benaoudia, M. M. (2025). Global patterns and trends in breast cancer incidence and mortality across 185 countries. Nature Medicine. https://doi.org/10.1038/s41591-025-03502-3
3. Singh, K., Grover, A., & Dhanasekaran, K. (2025). Unveiling the cancer epidemic in India: a glimpse into GLOBOCAN 2022 and past patterns. The Lancet Regional Health - Southeast Asia, 34, 100546. https://doi.org/10.1016/j.lansea.2025.100546
4. Nagrani, R., Budukh, A., Koyande, S., Panse, N., Mhatre, S., & Badwe, R. (2014). Rural urban differences in breast cancer in India. Indian Journal of Cancer, 51(3), 277. https://doi.org/10.4103/0019-509x.146793
5. Rohatgi, N., Limaye, S. A., Bahl, A., Sirohi, B., Sahni, S., Patil, D., Datta, V., Limited, D., Crook, T., Kurzrock, R., & Cristofanilli, M. (2023). Comprehensive molecular profiling from liquid and tissue for analysing the somatic mutational landscape in an Indian breast cancer cohort. Journal of Clinical Oncology, 41(16_suppl), e13080. https://doi.org/10.1200/jco.2023.41.16_suppl.e13080
6. Paul, A. D., Prabhu, A., Nidhi, S., Thomas, M. R., Shetty, R., Shenoy, P. U., & Das, R. (2024). Identification of novel genetic variants associated with oral squamous cell carcinoma (OSCC) in South-West coast of India using targeted exome sequencing. Gene, 148947. https://doi.org/10.1016/j.gene.2024.148947
7. Roy, M., Bal, A., Gupta, N., Prasad, K. T., Wakelee, H. A., & Singh, N. (2022). A brief report on the mutational landscape in non-small cell lung cancer of South Asian patients: Comparison at a US and an Indian Institution. Lung India, 39(4), 315–318. https://doi.org/10.4103/lungindia.lungindia_428_21
8. Sisoudiya, S. D., Houle, A. A., Fernando, T. M., Wilson, T. R., Schutzman, J. L., Lee, J. K., Schrock, A. B., Sokol, E. S., Sivakumar, S., Shi, Z., & Pathria, G. (2024). Abstract 6120: Genomic profiling of KRAS and EGFR-altered non-squamous non-small cell lung cancer reveal ancestry-specific co-alterations with therapeutic implications. Cancer Research, 84(6_Supplement), 6120. https://doi.org/10.1158/1538-7445.am2024-6120
9. GEO Accession viewer. (n.d.). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40206
10. TCGAbiolinks: Clinical data. (2025, March 6). https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/clinical.html
11. GEO Accession viewer. (n.d.). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62944
12. apply function - RDocumentation. (n.d.). https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/apply
13. Rahimi, M., Kakroodi, S. T., & Tajvidi, M. (2022). The Importance of RTK Signaling Genes and their Inhibitors in Breast Cancer. Journal of Obstetrics Gynecology and Cancer Research, 7(4), 258–271. https://doi.org/10.30699/jogcr.7.4.258
14. Wu, Z., Zhang, Y., Yue, J., & Zhou, T. (2020). Comprehensive analysis of the expression and prognosis for SFRPs in breast carcinoma. Cell Transplantation, 29, 096368972096247. https://doi.org/10.1177/0963689720962479





  


  





