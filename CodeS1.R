rm(list = ls())  
options(stringsAsFactors = F)


load('gene_beta.RData')
load('Clinical_Methylation.RData')
load('rnas.RData')


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager") ## installation required 


# ============================
# DMG in TCGA:

# install.packages("dplyr")
library(dplyr)  ## installation required 

clin_methy_prog = filter(clin_methy, patient.survival_days != -1 & patient.survival_conditions != -1 & patient.recurrency != -1)  ## delete items with label â€?-1â€?
rownames(clin_methy_prog)<-clin_methy_prog[,1]

clin_methy_prog = clin_methy_prog[,-c(1:27)]
clin_methy_prog[,3] = 1-as.numeric(clin_methy_prog[,4])
colnames(clin_methy_prog)[3]<-"patient.no_recurrency"
# 3rd col 1: no_recurrencyï¼?0: recurrencyã€?
# 4th col 1: recurrencyï¼?0: no_recurrencyã€?

# create group_list
clin_methy_copy<-clin_methy_prog

group_list<-clin_methy_copy[,-c(5:24282)]
group_list<-group_list[,-c(1:3)]
for(i in 1:135){
  if(group_list[i]=='1'){
    group_list[i]<-'recurrency'
  }
  else{
    group_list[i]<-'no_recurrency'
  }
}
#

group <- factor(clin_methy_copy$patient.recurrency, levels = c('no_recurrency','recurrency'))

 # clin_methy_prog 
 # rownamesï¼šgene_name ï¼›colnamesï¼špatient
clin_methy_prog = clin_methy_prog[,-c(1:4)]
clin_methy_prog <- t(clin_methy_prog)

design=model.matrix(~factor( group_list ))

for (i in 1:135){
  design[i,1] <- 1-design[i,2]
}

colnames(design) <- c('no_recurrency','recurrency')
rownames(design) <- colnames(clin_methy_prog)

clin_methy_data <- mapply(clin_methy_prog, FUN=as.numeric)
clin_methy_data <- matrix(data=clin_methy_data, nrow=24278, byrow=FALSE)

rownames(clin_methy_data)<-rownames(clin_methy_prog)
colnames(clin_methy_data)<-colnames(clin_methy_prog)

# ============================
# clin_methy_data

# BiocManager::install("limma")
library(limma)  ## installation required 

contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)

deg = function(clin_methy_data1,design,contrast.matrix){
  fit <- lmFit(clin_methy_data1,design)
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2)
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  head(nrDEG)
  return(nrDEG)
}

deg_Met = deg(clin_methy_data,design,contrast.matrix)
head(deg_Met)

save(deg_Met,file = 'degTCGA_Methyl.Rdata') 
 # degTCGA_Methyl.Rdata : used to extract P.values

# ============================
# DEG in TCGA:

# load('rnas.RData')
rnas<-t(rnas)
rnas_prog=rnas[rownames(rnas) %in% rownames(clin_methy_copy),]
rnas_prog<-t(rnas_prog)
deg = function(rnas_prog,design,contrast.matrix){
  fit <- lmFit(rnas_prog,design)
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2)
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  head(nrDEG)
  return(nrDEG)
}

deg_rna = deg(rnas_prog,design,contrast.matrix)
head(deg_rna)
save(deg_rna,file = 'degTCGA_RNA.Rdata')

# ============================
# Methylation and RNAseq data Opverlapping.

deg_Metp <- deg_Met[deg_Met$P.Value<=0.01,]
deg_Metp <- cbind(rownames(deg_Metp), deg_Metp)
colnames(deg_Metp)[1] <- 'gene_names'

deg_rnap <- deg_rna[deg_rna$P.Value<=0.01,]
deg_rnap <- cbind(rownames(deg_rnap), deg_rnap)
colnames(deg_rnap)[1] <- 'gene_names'

met_rna <- merge(deg_Metp,deg_rnap,by = 'gene_names')

# ============================
# Negative Correlation.

# Methylation
clin_methy_data1 <- cbind(rownames(clin_methy_data), clin_methy_data)
clin_methy_data1=clin_methy_data1[clin_methy_data1[,1] %in% met_rna[,1],]
clin_methy_data1 <- t(clin_methy_data1)
clin_methy_data1 <- clin_methy_data1[-1,]

 # RNAseq
rnas_prog1 <- cbind(rownames(rnas_prog), rnas_prog)
rnas_prog1=rnas_prog1[rnas_prog1[,1] %in% met_rna[,1],]
rnas_prog1 <- t(rnas_prog1)
rnas_prog1 <- rnas_prog1[-1,]

row.names <- c('r.squared','adj.r.squared','p-value','slope')
column.names <- met_rna[,1]
r2<-array(0,dim=c(4,16),dimnames = list(row.names,column.names))

for (i in 1:16){
  lm.reg <- lm(as.numeric(clin_methy_data1[,i])~ as.numeric(rnas_prog1[,i]))
  r2[1,i] <- summary(lm.reg)[8][[1]] # r.squared
  r2[2,i] <- summary(lm.reg)[9][[1]] # adj.r.squared
  r2[3,i] <- summary(lm.reg)$coefficients[,4][[2]] # p-value
  r2[4,i] <- summary(lm.reg)[4][[1]][2,1] # æ–œçŽ‡
}

# ============================
# 8 DNA Methylation generation

mylist <- list()
for (i in 1:16){
  if (r2[3,i]<0.05 & r2[4,i]<0 ){
    mylist <- c(mylist, colnames(r2)[i])
  }
}
# load('mylist.RData') # "mylist.RData" file is the selected eight genes.
clin_methy_data1<-t(clin_methy_data1)
clin_methy_prog <- clin_methy_data1[rownames(clin_methy_data1) %in% mylist,]
rnas_prog1<-t(rnas_prog1)
rnas_prog <- rnas_prog1[rownames(rnas_prog1) %in% mylist,]

# ============================
# plot code :

  # ===========
  # vocano plot
  # Methylation
nrDEG=deg_Met
head(nrDEG)
attach(nrDEG)
plot(logFC,-log10(P.Value))

  # RNA seq
nrDEG=deg_rna
head(nrDEG)
attach(nrDEG)
plot(logFC,-log10(P.Value))

  # ===========
  # pheatmap

# install.packages("pheatmap")
library(pheatmap) ## if plot heatmap, required.

  # ===========
  # 200 genes pheatmap
clin_methy_data1[1:4,1:4]
table(group_list)
x=deg_Met$logFC
names(x)=rownames(deg_Met) 
cg=c(names(head(sort(x),100)),
     names(tail(sort(x),100)))

pheatmap(clin_methy_data1[cg,],show_colnames =F,show_rownames = F) 
n=t(scale(t(clin_methy_data1[cg,])))
n[n>2]=2
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
ac=data.frame(g=group_list)
rownames(ac)=colnames(n)
pheatmap(nn,show_colnames =F,
         show_rownames = F,
         #cluster_cols = TRUE, 
         #annotation_col=ac, 
         cluster_rows = TRUE, 
         gaps_row = 100,
         cutree_rows = 2,
         filename = 'heatmap_top200_DEG.png')

  # ===========
  # 8 gene pheatmap
  # Methylation
met8 <- clin_methy_data1[rownames(clin_methy_data1) %in% mylist,]
rna8 <- rnas_data[rownames(rnas_data) %in% mylist,]
g<-as.matrix(group_list)
g<-t(g)
met8n<-rbind(g,met8)
met8n<-met8n[,order(met8n[1,])]
g<-data.frame(met8n[1,])
met8n<-met8n[-1,]
met8nn<-met8n
met8n <- mapply(met8n, FUN=as.numeric)
met8n <- matrix(data=met8n, nrow=8, byrow=FALSE)
rownames(met8n)<-rownames(met8nn)
colnames(met8n)<-colnames(met8nn)
pheatmap(met8n,show_colnames =F,show_rownames = F)
colnames(g)[1]<-'g'
pheatmap(met8n,show_colnames =F,
         show_rownames = F,
         cluster_cols = F, 
         annotation_col=g, 
         cluster_rows = F,
         cellwidth = 10, 
         cellheight = 100,
         #gaps_row = 100,
         #cutree_rows = 2,
         filename = 'heatmap_met8n_DEG.png')

  # 8 gene pheatmap
  # RNAseq
rna8 <- rnas_prog1[rownames(rnas_prog1) %in% mylist,]
g<-as.matrix(group_list)
g<-t(g)
rna8n<-rbind(g,rna8)
rna8n<-rna8n[,order(rna8n[1,])]
g<-data.frame(rna8n[1,])
rna8n<-rna8n[-1,]
rna8nn<-rna8n
rna8n <- mapply(rna8n, FUN=as.numeric)
rna8n <- matrix(data=rna8n, nrow=8, byrow=FALSE)
rownames(rna8n)<-rownames(rna8nn)
colnames(rna8n)<-colnames(rna8nn)
pheatmap(rna8n,show_colnames =F,show_rownames = F)
colnames(g)[1]<-'g'
pheatmap(rna8n,show_colnames =F,
         show_rownames = F,
         cluster_cols = F, 
         annotation_col=g, 
         cluster_rows = F,
         cellwidth = 10, 
         cellheight = 100,
         #gaps_row = 100,
         #cutree_rows = 2,
         filename = 'heatmap_rna8n_DEG.png')

# ============================
# SVM classifier & SVM Score :

# load('mylist.RData')
clin_methy_train <- clin_methy_data[rownames(clin_methy_data) %in% mylist,]
clin_methy_train <- t(clin_methy_train)
clin_methy_train <- cbind(rownames(clin_methy_train),clin_methy_train)
out<-design
colnames(clin_methy_train)[1] <-'ID'
out <- cbind(rownames(out),out)
colnames(out)[1] <-'ID'
clin_methy_train <- merge(clin_methy_train, out, by = 'ID')
clin_methy_train <- clin_methy_train[,-1]
clin_methy_tr <- mapply(clin_methy_train, FUN=as.numeric)
clin_methy_tr <- matrix(data=clin_methy_tr, nrow=135, byrow=FALSE)
rownames(clin_methy_tr)<-rownames(clin_methy_train)
colnames(clin_methy_tr)<-colnames(clin_methy_train)

clin_methy_tr<-clin_methy_tr[,-9]
###
save(clin_methy_tr,file = 'clin_methy_train.RData')

# install.packages("e1071")
library(e1071)  ## installation required.

#svmfit = svm(patient.recurrency ~ C17orf62 + CD3E + CUL7 + LOC100132215 + MPHOSPH10 + PDIA6 + TIGIT + ZFYVE21, data = clin_methy_tr, kernel = "linear", cost = 10, scale = FALSE, C-classification)
x <- clin_methy_tr[,-9]
y <- as.factor(clin_methy_tr[,-c(1:8)])
svmfit <- svm(y ~ x, data = clin_methy_tr,  kernel = "linear")
pred <- predict(svmfit,x)
table(pred,y)

print(svmfit)
#plot(svmfit, clin_methy_tr)
svmfit <- svm(as.factor(recurrency) ~ .,data = clin_methy_tr, kernel = "linear")

clin_methy_svm <- clin_methy_data1[rownames(clin_methy_data1) %in% mylist,]
clin_methy_svm<-t(clin_methy_svm)
clin_methy_svm1<-clin_methy_svm
clin_methy_svm <- mapply(clin_methy_svm, FUN=as.numeric)
clin_methy_svm <- matrix(data=clin_methy_svm, nrow=135, byrow=FALSE)
rownames(clin_methy_svm)<-rownames(clin_methy_svm1)
colnames(clin_methy_svm)<-colnames(clin_methy_svm1)

#a <- 0.48*clin_methy_svm[,1]-0.61*clin_methy_svm[,2]+0.17*clin_methy_svm[,3]+0.37*clin_methy_svm[,4]+0.12*clin_methy_svm[,5]+0.25*clin_methy_svm[,6]-0.49*clin_methy_svm[,7]+0.49*clin_methy_svm[,8]-0.74
svmfit[29]#svm decision values
decision <- as.data.frame(svmfit[29])

 # Univariate analysis and multivariate analysis
 # T-test and cox
clin <- clin_methy[,-c(32:24309)]
clin = filter(clin, patient.survival_days != -1 & patient.survival_conditions != -1 & patient.recurrency != -1)

clin_t <- clin[,c(2,28,31)]
clin_t <- mapply(clin_t, FUN=as.numeric)
clin_t <- matrix(data=clin_t, nrow=135, byrow=FALSE)
clin_t[,1] <- -clin_t[,1]/365 #age
colnames(clin_t)<-c('age','survival_days','recurrency')
a<-t.test(clin_t[,1] ~ clin_t[,3],clin_t, var.equal=TRUE)
b<-t.test(clin_t[,1] ~ clin_t[,3],clin_t)
boxplot(clin_t[,1] ~ clin_t[,3],data=clin_t)


#library(vcd)
clin_f <- clin[,c(23:25,31)]
clin_f[,4] <- mapply(clin_f[,4], FUN=as.numeric)

mytable <- xtabs(~clin$patient.survival_conditions + clin$patient.recurrency, data = clin)
fisher.test(mytable)

clin_1 = filter(clin, patient.recurrency == 1)
clin_0 = filter(clin, patient.recurrency == 0)

# install.packages("survival")
library("survival")  ## installation required

# install.packages("survminer")
library("survminer")  ## if drawing survival curves using 'ggsurvplot', required.

clin_c <- clin[,c(2,12,17,20,23:26,28,29,31,14)] # [14] "patient.histological_type" 
clin_c[,1] <- mapply(clin_c[,1], FUN=as.numeric)
clin_c$patient.days_to_birth <- -clin_c[,1]/365 # age
clin_c$age = ifelse(clin_c$patient.days_to_birth<40,1,ifelse(clin_c$patient.days_to_birth<60,2,ifelse(clin_c$patient.days_to_birth>60,3,0))) # 1:<40; 2:<60; 3:>60

colnames(clin_c)[9]<-'time' # disease-free survival time
clin_c[,9] <- mapply(clin_c[,9], FUN=as.numeric)

clin$figo = ifelse(clin$patient.stage_event.clinical_stage=='stage ivb',4,ifelse(clin$patient.stage_event.clinical_stage=='stage iva',4,ifelse(clin$patient.stage_event.clinical_stage=='stage iiib',3,ifelse(clin$patient.stage_event.clinical_stage=='stage iiia',3,ifelse(clin$patient.stage_event.clinical_stage=='stage iii',3,ifelse(clin$patient.stage_event.clinical_stage=='stage iib',2,ifelse(clin$patient.stage_event.clinical_stage=='stage iia2',2,ifelse(clin$patient.stage_event.clinical_stage=='stage iia1',2,ifelse(clin$patient.stage_event.clinical_stage=='stage iia',2,ifelse(clin$patient.stage_event.clinical_stage=='stage ii',2,ifelse(clin$patient.stage_event.clinical_stage=='stage ib2',1,ifelse(clin$patient.stage_event.clinical_stage=='stage ib1',1,ifelse(clin$patient.stage_event.clinical_stage=='stage ib',1,ifelse(clin$patient.stage_event.clinical_stage=='stage ia1',1,ifelse(clin$patient.stage_event.clinical_stage=='stage i',1,0)))))))))))))))

clin_c$patient.histological_type=ifelse(clin_c$patient.histological_type == 'cervical squamous cell carcinoma','squamous',ifelse(clin_c$patient.histological_type == 'adenosquamous','adenosquamous','adenos'))

colnames(clin_c)[10] <- "status" # "recurrence status"
clin_c$status=ifelse(clin_c$status == 1,2,ifelse(clin_c$status==0,1,0))

clin_c<-cbind(clin_c,clin$figo)

#load('SVMScoreTCGA.RData') # decision is the SVM Score in TCGA. Generated in the previous parts.
clin_c<-cbind(clin_c,decision) 

fit <- survfit(Surv(time, status) ~ age, data = clin_c)

clin_c<-clin_c[,-c(1,2,3)]
clin_cox<-clin_c
colnames(clin_cox)[11]<-'figo'
colnames(clin_cox)[12]<-'svm'
res.cox <- coxph(Surv(time, status) ~ patient.neoplasm_histologic_grade 
                 + patient.stage_event.tnm_categories.pathologic_categories.pathologic_m 
                 + patient.stage_event.tnm_categories.pathologic_categories.pathologic_n 
                 + patient.stage_event.tnm_categories.pathologic_categories.pathologic_t 
                 + patient.vital_status 
                 + patient.recurrency 
                 + figo 
                 + svm, data = clin_cox)


res.cox <- coxph(Surv(time, status) ~ #
                   age
                 + patient.vital_status
                 + patient.stage_event.tnm_categories.pathologic_categories.pathologic_m 
                 + patient.stage_event.tnm_categories.pathologic_categories.pathologic_n
                 + figo 
                 + svm, data = clin_cox)
print(res.cox)

summary(res.cox)
ggsurvplot(survfit(res.cox), data = clin_c, palette = "#2E9FDF", ggtheme = theme_minimal(), legend = "none")

ggsurvplot(survfit(res.cox), data = clin_c$age, conf.int = TRUE, legend.labs=c("age=1", "age=2", "age=3"), ggtheme = theme_minimal())

# "svmscoreTCGA" is further used for ROC plot and AUC calculation in CodeS2.R.

svmscoreTCGA<-cbind(clin_methy_tr,decision) 
colnames(svmscoreTCGA)[10] <- 'svm'
