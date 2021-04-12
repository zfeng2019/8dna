# 8dna

#### Code Summary:


​	**Code_S1.R**

​		Data cleansing (TCGA)

​		Filter samples with Recur/ no_Recur

​		Differential analysis with volcano plot 

​		Negative correlation between overlapped Methylation and RNAseq data

​		Heatmap

​		SVM classifier & SVM Score

​		T-test and Cox analysis




​	**Code_S2.R**

​		Data cleansing (GEO)

​		Calculate SVM Score

​		Plot t-SNE visualization

​		Calculate AUC with TCGA dataset (self-validation)

​		Calculate AUC with GEO dataset

​	

​	**Code_S3.R**

​		Discussion on histological types 

​		Univariate and Cox regression analysis for squamous carcinomas only



#### Required data files:

​	**(1)** *'gene_beta.RData'*

​		*'gene_beta.RData'* includes data with 24,278 DNA methylation genes of 309 patients. 



​	**(2)** *'Clinical_Methylation.RData'*
​		Binding the clinical information;

​		It includes data with 31 clinical features and 24,278 DNA methylation genes 

​		of 248 patients.



​	**(3)** *'rnas.RData'*
​		*'rnas.RData'* includes data with 20,501 expressed genes of 306 patients.

* Due to the limited max upload size, these three data files in the current study are available from the corresponding author on reasonable request.

#### R packages:

​	**Required:**

​		BiocManager (limma, GEOquery); 

​		dplyr;



​	**Optional:** 

​		e1071 (for SVM)

​		pheatmap (for heatmap plot)

​		survival (for survival analysis)

​		survminer (for drawing survival curves)

​		Rtsne (for t-SNE visualization)

​		ggplot2 (for ROC)

​		plotROC (for AUC) 

​		vcd (for fisher test)



#### Contact address: 

​	zfeng2019@foxmail.com

​	Zhen Feng (First Affiliated Hospital of Wenzhou Medical University)
