## T-test and cox

# rm(list = ls())
# options(stringsAsFactors = F)

# load('Clinical_Methylation.RData')
library(dplyr)  ## loaded in previous R scripts
library(survival)  ## loaded in previous R scripts
library(survminer)  ## loaded in previous R scripts
library(vcd)  ## if performing fisher test, required.

clin <- clin_methy[,-c(32:24309)]

clin = filter(clin, patient.survival_days != -1 & patient.survival_conditions != -1 & patient.recurrency != -1)

clin_t <- clin[,c(2,28,31)] #t-test
clin_t <- mapply(clin_t, FUN=as.numeric)
clin_t <- matrix(data=clin_t, nrow=135, byrow=FALSE)
clin_t[,1] <- -clin_t[,1]/365 # age
colnames(clin_t)<-c('age','survival_days','recurrency')
a<-t.test(clin_t[,1] ~ clin_t[,3],clin_t, var.equal=TRUE)
print(a) #Two Sample t-test: age
boxplot(clin_t[,1] ~ clin_t[,3],data=clin_t)


clin_f <- clin[,c(23:25,31)] # m n t
clin_f[,4] <- mapply(clin_f[,4], FUN=as.numeric)

mytable <- xtabs(~clin$patient.survival_conditions + clin$patient.recurrency, data = clin)
fisher.test(mytable)

clin_1 = filter(clin, patient.recurrency == 1)
clin_0 = filter(clin, patient.recurrency == 0)

clin_cox<-clin_c
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

load('SVMScoreTCGA.RData') #decision is the SVM Score in TCGA
clin_c<-cbind(clin_c,decision) 

clin_c<-clin_c[,-c(1,2,3)]
clin_cox<-clin_c
colnames(clin_cox)[12]<-'svm'
colnames(clin_cox)[11]<-'figo'
save(clin_cox,file = 'clin_cox.RData')
res.cox <- coxph(Surv(time, status) ~ #
                 + patient.stage_event.tnm_categories.pathologic_categories.pathologic_m 
                 + patient.stage_event.tnm_categories.pathologic_categories.pathologic_n
                 + figo 
                 + svm, data = clin_cox)
print(res.cox)

summary(res.cox)

# univariate analysis

mytable <- xtabs(~clin_cox$patient.stage_event.tnm_categories.pathologic_categories.pathologic_m + clin_cox$patient.recurrency, data = clin_cox)
fisher.test(mytable)

mytable <- xtabs(~clin_cox$patient.stage_event.tnm_categories.pathologic_categories.pathologic_n + clin_cox$patient.recurrency, data = clin_cox)
fisher.test(mytable)

mytable <- xtabs(~clin_cox$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t + clin_cox$patient.recurrency, data = clin_cox)
fisher.test(mytable)

#test for histological_type
table(clin_cox$patient.histological_type)

mytable <- xtabs(~clin_cox$patient.histological_type + clin_cox$patient.recurrency, data = clin_cox)
fisher.test(mytable)

a<-t.test(clin_cox$svm ~ clin_cox$patient.histological_type, data = clin_cox)
print(a)

b<-chisq.test(clin_cox$patient.histological_type,clin_cox$svm)
print(b)


## univariate and Cox regression analysis for squamous carcinomas only
clin_cox<-clin_cox[,-10]
clin_cox<-cbind(clin_t[,1],clin_cox)
colnames(clin_cox)[1]<-'age'
table(clin_cox$patient.histological_type)
clin_cox_hist<-filter(clin_cox, patient.histological_type == 'squamous')
 # univariate analysis
clin_1 = filter(clin_cox_hist, patient.recurrency == 1)
clin_0 = filter(clin_cox_hist, patient.recurrency == 0)

histcox <- xtabs(~clin_cox_hist$patient.histological_type + clin_cox_hist$patient.recurrency, data = clin_cox_hist)
fisher.test(histcox)

table(clin_0$figo)
table(clin_1$figo)
#age
a<-t.test(clin_cox_hist$age ~ clin_cox_hist$patient.recurrency, data = clin_cox_hist)
sd(clin_1$age)
sd(clin_0$age[1:71]) #standard diviation
print(a)
#svm
a<-t.test(clin_cox_hist$svm ~ clin_cox_hist$patient.recurrency, data = clin_cox_hist)
sd(clin_1$svm)
sd(clin_0$svm) #standard diviation
print(a)

# keep the squamous carcinomas only
# 
clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t = ifelse(clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t=='t1a1'|clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t=='t1b'|clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t=='t1b1'|clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t=='t1b2','t1',
                                                                                             ifelse(clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t=='t2' | clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t=='t2a' | clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t== 't2a1' | clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t== 't2a2' ,'t2',
                                                                                                    ifelse(clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t=='t3','t3',
                                                                                                           ifelse(clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t=='t4','t4','tx')))) #
###############################
clin_cox_hist1<-clin_cox_hist
clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t = ifelse(clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t=='t1',1,
                                                                                             ifelse(clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t=='t2',2,
                                                                                                    ifelse(clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t=='t3',3,
                                                                                                           ifelse(clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t=='t4',4,0)))) #

clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_m = ifelse(clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_m=='m0',0,
                                                                                             ifelse(clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_m=='m1',1,2)) #mx:2

clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_n = ifelse(clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_n=='n1',1,
                                                                                             ifelse(clin_cox_hist$patient.stage_event.tnm_categories.pathologic_categories.pathologic_n=='n0',0,2)) #nx:2

###############################



res.cox <- coxph(Surv(time, status) ~ #
                   + patient.stage_event.tnm_categories.pathologic_categories.pathologic_m 
                 + patient.stage_event.tnm_categories.pathologic_categories.pathologic_n
                 + patient.stage_event.tnm_categories.pathologic_categories.pathologic_t
                 #+ figo 
                 + svm, data = clin_cox_hist)
print(res.cox)
summary(res.cox)

res.cox <- coxph(Surv(time, status) ~ patient.stage_event.tnm_categories.pathologic_categories.pathologic_m 
                 + patient.stage_event.tnm_categories.pathologic_categories.pathologic_n
                 + patient.stage_event.tnm_categories.pathologic_categories.pathologic_t
                 + svm, data = clin_cox_hist)
print(res.cox)
summary(res.cox)

# screen Nx items
# delete items with label “-1”
clin_cox_N <- filter(clin_cox, patient.stage_event.tnm_categories.pathologic_categories.pathologic_n != 'nx')
clin_1 = filter(clin_cox_N, patient.recurrency == 1)
clin_0 = filter(clin_cox_N, patient.recurrency == 0)

table(clin_cox_N$patient.stage_event.tnm_categories.pathologic_categories.pathologic_n)
nxcox <- xtabs(~clin_cox_N$patient.neoplasm_histologic_grade + clin_cox_N$patient.recurrency, data = clin_cox_N)
fisher.test(nxcox)

a<-t.test(clin_cox_N$svm ~ clin_cox_N$patient.stage_event.tnm_categories.pathologic_categories.pathologic_n, data = clin_cox_N)
print(a)
sd(clin_1$svm)
sd(clin_0$svm)

##============================
# delete Nx only, and 81 left.

clin_cox_N1<-clin_cox_N

clin_cox_N$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t = ifelse(clin_cox_N$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t=='t1a1'|clin_cox_N$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t=='t1b'|clin_cox_N$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t=='t1b1'|clin_cox_N$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t=='t1b2',1,
                                                                                             ifelse(clin_cox_N$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t=='t2' | clin_cox_N$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t=='t2a' | clin_cox_N$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t== 't2a1' | clin_cox_N$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t== 't2a2' ,2,
                                                                                                    ifelse(clin_cox_N$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t=='t3',3,
                                                                                                           ifelse(clin_cox_N$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t=='t4',4,0)))) #

clin_cox_N$patient.stage_event.tnm_categories.pathologic_categories.pathologic_n = ifelse(clin_cox_N$patient.stage_event.tnm_categories.pathologic_categories.pathologic_n=='n1',1,
                                                                                             ifelse(clin_cox_N$patient.stage_event.tnm_categories.pathologic_categories.pathologic_n=='n0',0,2)) #nx:2
##

res.cox <- coxph(Surv(time, status) ~ #
                 + patient.stage_event.tnm_categories.pathologic_categories.pathologic_n
                 + patient.stage_event.tnm_categories.pathologic_categories.pathologic_t 
                 + svm, data = clin_cox_N)
print(res.cox)
summary(res.cox)

##============================
## delete all N stage

clin_cox$patient.stage_event.tnm_categories.pathologic_categories.pathologic_m = ifelse(clin_cox$patient.stage_event.tnm_categories.pathologic_categories.pathologic_m=='m0',0,
                                                                                             ifelse(clin_cox$patient.stage_event.tnm_categories.pathologic_categories.pathologic_m=='m1',1,2)) #mx:2

res.cox <- coxph(Surv(time, status) ~ #
                   + patient.stage_event.tnm_categories.pathologic_categories.pathologic_m
                 + figo
                 + svm, data = clin_cox)
print(res.cox)
summary(res.cox)