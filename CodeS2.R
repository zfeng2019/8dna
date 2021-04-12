f='GSE30759_eSet.Rdata'
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30759
library(GEOquery)
if(!file.exists(f)){
  gset <- getGEO('GSE30759', destdir=".",
                 AnnotGPL = F,     
                 getGPL = F)  
  save(gset,file=f)   
}
load('GSE30759_eSet.Rdata')  
a=gset[[1]] #
dat=exprs(a) #
dim(dat)#
# GPL6244
dat[1:4,1:4] 
#dat=log2(dat+0.1)
pd=pData(a) 
phe=pd[,c(43:46,48)]
gpl <- getGEO('GPL8490', destdir=".")
colnames(Table(gpl))
probe2gene=Table(gpl)[,c(1,23)] 
head(probe2gene)
# library(stringr)
save(probe2gene,file='probe2gene.Rdata')
ids=probe2gene 
colnames(ids)=c('probe_id','symbol')  
dat<-cbind(rownames(dat),dat)
colnames(dat)[1] <- 'probe_id'

dat1 <- merge(ids,dat,by='probe_id')
dat1 <- dat1[,-1]

library(dplyr) # loaded in CodeS1.R
dim(dat1)

dat1[,c(2:64)] %>% mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.numeric) -> dat2
dat2 <- cbind(dat1$symbol,dat2)
colnames(dat2)[1] <- 'gene_name'

dat3<-aggregate(dat2[,c(2:64)],by=list(gene_name=dat2$gene_name),FUN = "max")
rownames(dat3)<-dat3[,1]
dat3$gene_name<-NULL

dat4=dat3[rownames(dat3) %in% mylist,]

save(dat4,file = 'GEOMethy64.RData')

GEOMethy64<-dat4
dim(GEOMethy64)
class(GEOMethy64)

GEOMethy64<-GEOMethy64[complete.cases(GEOMethy64),]
boxplot(GEOMethy64,las=2)

# ============================
# dealing with clinical information
phe=pd[,c(2,12,13,14,15,16,20,21,22,23)]

GEOMethy64<-t(GEOMethy64)
GEOMethy64<-cbind(rownames(GEOMethy64),GEOMethy64)
colnames(GEOMethy64)[1]<-'geo_accession'
clin<-merge(phe,GEOMethy64,by = 'geo_accession')
rownames(clin)<-clin[,1]

clin = filter(clin, characteristics_ch1.2 == 'status(0=normal,1=cancer): 1')
## clin = filter(clin, characteristics_ch1.10 != 'survival time (years): NA')
rownames(clin)<-clin[,1]
#clin_methy_prog = filter(clin_methy, patient.survival_days != -1 & patient.survival_conditions != -1 & patient.recurrency != -1)

colnames(clin)[2]<-'status'
colnames(clin)[3]<-'age'
colnames(clin)[4]<-'stage'
colnames(clin)[5]<-'histology'
colnames(clin)[6]<-'grade'
colnames(clin)[7]<-'survival time'
colnames(clin)[8]<-'death status'
colnames(clin)[9]<-'recurrence free survival'
colnames(clin)[10]<-'relapse'

#
k=48
for (i in 1:k) {
  b=unlist(strsplit(clin$status[i], ":", fixed = TRUE))
  clin$status[i]=1
  clin$status[i]=b[2]
}
for (i in 1:k) {
  b=unlist(strsplit(clin$age[i], ":", fixed = TRUE))
  clin$age[i]=1
  clin$age[i]=b[2]
}
for (i in 1:k) {
  b=unlist(strsplit(clin$stage[i], ":", fixed = TRUE))
  clin$stage[i]=1
  clin$stage[i]=b[2]
}
for (i in 1:k) {
  b=unlist(strsplit(clin$histology[i], ":", fixed = TRUE))
  clin$histology[i]=1
  clin$histology[i]=b[2]
}
for (i in 1:k) {
  b=unlist(strsplit(clin$grade[i], ":", fixed = TRUE))
  clin$grade[i]=1
  clin$grade[i]=b[2]
}
for (i in 1:k) {
  b=unlist(strsplit(clin$`survival time`[i], ":", fixed = TRUE))
  clin$`survival time`[i]=1
  clin$`survival time`[i]=b[2]
}
for (i in 1:k) {
  b=unlist(strsplit(clin$`death status`[i], ":", fixed = TRUE))
  clin$`death status`[i]=1
  clin$`death status`[i]=b[2]
}
for (i in 1:k) {
  b=unlist(strsplit(clin$`recurrence free survival`[i], ":", fixed = TRUE))
  clin$`recurrence free survival`[i]=1
  clin$`recurrence free survival`[i]=b[2]
}
for (i in 1:k) {
  b=unlist(strsplit(clin$relapse[i], ":", fixed = TRUE))
  clin$relapse[i]=1
  clin$relapse[i]=b[2]
}

clin<-clin[,-c(1,4,6)] #delete NA

clin<-clin[complete.cases(clin),]
clin1 <- mapply(clin, FUN=as.numeric)
clin1 <- matrix(data=clin1, nrow=48, byrow=FALSE)
rownames(clin1)<-rownames(clin)
colnames(clin1)<-colnames(clin)

prog<- array(-1:-1, c(1,48))
rownames(prog) <- 'prog'
colnames(prog)<-rownames(clin1)

for (i in 1:48) {
  
  if (clin1[i,4] == clin1[i,6] & clin1[i,4]>2){# good outcomes，porg=0
    prog[i] = 0
  }
  
  else if (clin1[i,4]-clin1[i,6]>0){ # bad outcomes，prog=1
    prog[i] = 1
  }
}
table(prog)
prog<-t(prog)
clin2<-cbind(clin1,prog)

clin2 = filter(as.data.frame(clin2), prog != -1)

save(clin2,file = 'GEOdata2.RData')

##============================
#svm score

svmscore<- array(-1:-1, c(1,47))
rownames(svmscore) <- 'svm_score'
colnames(svmscore)<-rownames(clin2)

svmscore <- 0.48*clin2[,8]-0.61*clin2[,9]+0.17*clin2[,10]+0.12*clin2[,11]+0.25*clin2[,12]
svmscore<-data.frame(svmscore)
colnames(svmscore) <- 'svm_score'
rownames(svmscore)<-rownames(clin2)
colnames(clin2)[13]<-'prog'
clin3<-cbind(clin2,svmscore)

# sort clin3 by prog
save(clin3,file = 'clin3.RData')
clin4 <- clin3[order(clin3$prog,clin3$svm_score),]

clin4[,4]<-clin4[,4]*365
clin4[,6]<-clin4[,6]*365
svmscoreTCGA<-as.data.frame(clin_methy_tr)
clin4$N_C17orf62 <- mean(svmscoreTCGA$C17orf62)/mean(clin4$C17orf62)*clin4$C17orf62
clin4$N_CD3E <- mean(svmscoreTCGA$CD3E)/mean(clin4$CD3E)*clin4$CD3E
clin4$N_CUL7 <- mean(svmscoreTCGA$CUL7)/mean(clin4$CUL7)*clin4$CUL7
clin4$N_MPHOSPH10 <- mean(svmscoreTCGA$MPHOSPH10)/mean(clin4$MPHOSPH10)*clin4$MPHOSPH10
clin4$N_PDIA6 <- mean(svmscoreTCGA$PDIA6)/mean(clin4$PDIA6)*clin4$PDIA6

N_svmscore<- array(-1:-1, c(1,47))
rownames(N_svmscore) <- 'N_svm_score'
colnames(N_svmscore)<-rownames(clin4)
N_svmscore <- 0.48*clin4$N_C17orf62-0.61*clin4$N_CD3E+0.17*clin4$N_CUL7+0.25*clin4$N_PDIA6

N_svmscore<-data.frame(N_svmscore)
colnames(N_svmscore) <- 'N_svmscore'
rownames(N_svmscore)<-rownames(clin4)
clin_N<-cbind(clin4,N_svmscore)

## sort clin_N by prog
clin_N <- clin_N[order(clin_N$svm,clin_N$N_svmscore),]  ##
#clin_N1<-clin_N[clin_N$dfs!=0,]
#clin_N1<-clin_N1[clin_N$dfs!=83,]
#clin_N1<-clin_N1[-1,]
#clin_N<-clin_N1
colnames(clin_N)[4]<-"time"
colnames(clin_N)[6]<-"dfs"

save(clin_N,file = 'clin_N.RData')
write.csv(clin_N,file = 'clin_N.csv')
# thres: reduced SVM score threshold
#thres <- 0.74-0.49*mean(svmscoreTCGA$ZFYVE21) +0.49*mean(svmscoreTCGA$TIGIT)-0.37*mean(svmscoreTCGA$LOC100132215)
#thres <- thres - 0.12*mean(svmscoreTCGA$MPHOSPH10)

clin0<-filter(svmscoreTCGA, svmscoreTCGA$recur == 0)
clin1<-filter(svmscoreTCGA, svmscoreTCGA$recur == 1)

thres0<-0.74-0.49*mean(clin0$ZFYVE21) +0.49*mean(clin0$TIGIT)-0.37*mean(clin0$LOC100132215)
thres1<-0.74-0.49*mean(clin1$ZFYVE21) +0.49*mean(clin1$TIGIT)-0.37*mean(clin1$LOC100132215)

thres0<-thres0-0.12*mean(clin0$MPHOSPH10)
thres1<-thres1-0.12*mean(clin1$MPHOSPH10)
clin_N$svm = ifelse(clin_N$N_svmscore<=thres1,0,ifelse(clin_N$N_svmscore>thres1,1,-1))

library(survival)  ## loaded in CodeS1.R
library(survminer)  ## loaded in CodeS1.R
clin_N$svm = ifelse(clin_N$N_svmscore<=thres,0,ifelse(clin_N$N_svmscore>thres,1,-1))
clin_N=filter(clin_N,clin_N$svm!=-1)
dim(clin_N)
clin_N <- clin_N[order(clin_N$recurrence,clin_N$svm),]
clin_N<-clin_N[,-1]

clin_N$recurrence <- ifelse(clin_N$time>clin_N$dfs,1,ifelse(clin_N$time==clin_N$dfs,0,-1))
table(clin_N$recurrence)

save(clin_N,file = 'Surv_clin_N.RData')
#fit <- survfit(Surv(dfs2,relapse) ~ svm, data = clin_N)
fit <- survfit(Surv(dfs,recurrence) ~ svm, data = clin_N)
ggsurvplot(fit)
surv_pvalue(fit)

# ============================
# t-SNE

# install.packages("Rtsne")
library(Rtsne)  ## if plot t-SNE, required.

colors = rainbow(length(unique(y)))
tsne <- Rtsne(x, dims = 2, perplexity=20, verbose=TRUE, max_iter = 10000)
plot(tsne$Y, main="tsne")
text(tsne$Y, labels=y, col=colors[y])

## ROC
#
install.packages("ggplot2")
library(ggplot2)
install.packages("plotROC")
library(plotROC)  ## if using ROC and AUC, required.

# calculate AUC with TCGA dataset (self-validation)
tcgaROC <- ggplot(svmscoreTCGA[,c(9,10)], aes(d = 1-svmscoreTCGA$recurrency, m = svmscoreTCGA$svm)) + geom_roc()
tcgaROC
calc_auc(tcgaROC)

# calculate AUC with GEO dataset
basicplot <- ggplot(clin_N[,c(19,22)], aes(d = 1-clin_N$recurrence, m = clin_N$`N_svmscore`)) + geom_roc()

basicplot
calc_auc(basicplot)
