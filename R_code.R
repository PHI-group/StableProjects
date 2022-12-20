

##subcortical volume comparison
alldata<-read.csv('/')
#alldata<-subset(alldata,alldata$age)  ##change age range
names(alldata)
nrow(alldata)
subcort<-alldata[,c(7:20)]
com<- matrix(0, nrow=14,ncol=8)

for (i in 1:14) {
  tem<-lm(subcort[,i]~scale(alldata$age)*alldata$diagnosis+alldata$sex*alldata$diagnosis+scale(alldata$TIV))
  tem2<-summary(step(tem))
  print(tem2)
}

for (i in 1:14) {
  tem<-lm(subcort[,i]~scale(alldata$age)+alldata$diagnosis+alldata$sex+scale(alldata$TIV))
  tem2<-summary(tem)
  com[i,1]<-tem2$coefficients[2,1]
  com[i,2]<-tem2$coefficients[2,4]
  com[i,3]<-tem2$coefficients[3,1]
  com[i,4]<-tem2$coefficients[3,4]
  com[i,5]<-tem2$coefficients[4,1]
  com[i,6]<-tem2$coefficients[4,4]
  com[i,7]<-tem2$coefficients[5,1]
  com[i,8]<-tem2$coefficients[5,4]
}
p<-com[,c(2,4,6,8)]
p<-as.data.frame(p)

p_fdr<-matrix(0, nrow=14,ncol=4)
for (i in 1:4) {
  p_fdr[,i]<-p.adjust(p[,i], method = 'fdr')
}


##correlation
for (i in 1:14) {
  tem<-lm(subcort[,i]~alldata$TIV)
  tem2<-residuals(tem)
  subcort[,i]<-tem2
}
names(alldata)
all<-cbind(subcort,alldata[,c(2:5)])
hcdat <- subset(all,all$diagnosis=='HC')
padat <- subset(all,all$diagnosis=='SAD')

library(ppcor)
t<-pcor.test(alldata$Left.Thalamus.Proper,alldata$DSRSC,c(alldata$sex,alldata$age)) ##SASC or SAD-scale or FNE-scale
t$p.value
for (i in 1:14) {
  tem<-pcor.test(padat[,i],padat$DSRSC,c(padat$age,padat$sex))
  com[i,1]<-tem$estimate
  com[i,2]<-tem$p.value
}

p<-com[,c(2)]
p<-as.data.frame(p)

p_fdr<-p.adjust(p, method = 'fdr')


##subcortical covairance
#alldata<-subset(alldata,alldata$age)  ##change age range
hcdat <- subset(alldata,alldata$diagnosis=='HC')
hcres <- NULL
names(hcdat)
for (i in 4:17){
  tmp_res <- resid(lm(hcdat[[i]]~hcdat$age+hcdat$sex+hcdat$TIV))##+hcdat$SASC总分/+hcdat$抑郁总分DSRSC
  hcres <- cbind(hcres, tmp_res)
}
names(hcdat)
colnames(hcres) <- names(hcdat)[4:17]
hccor <- cor(hcres)
diag(hccor) <- 0
## compute covariance in PA group
padat <- subset(alldata,alldata$diagnosis=='SAD')
pares <- NULL
for (i in 4:17){
  tmp_res <- resid(lm(padat[[i]]~padat$age+padat$sex+padat$TIV))##+padat$SASC总分/+padat$抑郁总分DSRSC
  pares <- cbind(pares, tmp_res)
}
colnames(pares) <- names(padat)[4:17]
pacor <- cor(pares)
diag(pacor) <- 0
table(alldata$diagnosis)

table(alldata$diagnosis)
diffcor <- hccor - pacor
Z <- (atanh(pacor)-atanh(hccor))/sqrt(1/(67-3)+1/(76-3))
P <- 2*(1-pnorm(abs(Z)))
P
Pfdr <- matrix(0, nrow=nrow(P), ncol=nrow(P), dimnames=list(rownames(P),colnames(P)))

Pfdr[lower.tri(P)] <- p.adjust(P[lower.tri(P)], method = 'fdr')


##subcortical cortical covariance centrality
#alldata<-subset(dat,dat$age>=15) ##change age range
table(alldata$group)
names(alldata)

hcdat <- subset(alldata,alldata$diagnosis=='HC')
hcres <- NULL
names(hcdat)
for (i in 10:91){
  tmp_res <- resid(lm(hcdat[[i]]~hcdat$age+hcdat$sex+hcdat$TIV))##+hcdat$SASC总分/+hcdat$抑郁总分DSRSC
  hcres <- cbind(hcres, tmp_res)
}
names(hcdat)
colnames(hcres) <- names(hcdat)[10:91]

padat <- subset(alldata,alldata$diagnosis=='SAD')
pares <- NULL
for (i in 10:91){
  tmp_res <- resid(lm(padat[[i]]~padat$age+padat$sex+padat$TIV))##+padat$SASC总分/+padat$抑郁总分DSRSC
  pares <- cbind(pares, tmp_res)
}
colnames(pares) <- names(padat)[10:91]

library(tidyverse)
library(psych)
library(progress)
hcres<-as.data.frame(hcres)
head(hcres)
pares<-as.data.frame(pares)
data<-rbind(hcres,pares)
nrow(data)
data$group <- as.factor(c(rep('HC',36),rep('PA',36)))  # participants number in HC and SAD group
head(data)
ncol(data)
table(data$group)
PA_origin <- data %>% 
  dplyr::filter(group == 'PA')
HC_origin <- data %>% 
  dplyr::filter(group == 'HC')

corr_diff <- function(data1,data2){
  data1_sub <- data1 %>% 
    select(-83) %>% 
    select(1:14)
  data1_cort <- data1 %>% 
    select(-83) %>% 
    select(-c(1:14))
  
  data2_sub <- data2 %>% 
    select(-83) %>% 
    select(1:14)
  data2_cort <- data2 %>% 
    select(-83) %>% 
    select(-c(1:14))
  
  data1_cor <- corr.test(data1_sub,data1_cort)
  data2_cor <- corr.test(data2_sub,data2_cort)
  
  return(data1_cor$r - data2_cor$r)
}

origin_diff <- corr_diff(PA_origin,HC_origin)
origin_diff_mean<-rowMeans(origin_diff)
origin_diff_mean<-as.data.frame(origin_diff_mean)

n_perm <- 10000 ### the number of permutation test 
n_PA <- 36 ### the participants number in PA

perm_index_PA <- list()

set.seed(123)
for(i in 1:n_perm) {
  perm_index_PA[[i]] <- sample(x = seq_len(nrow(data)),
                               size = n_PA,
                               replace = F)
}

r_diff_perm <- list()
mean_diff<-list()
pb <- progress_bar$new(total = n_perm, 
                       format = '[:bar] :percent eta: :eta',
                       clear = FALSE)


for (i in seq_along(perm_index_PA)) {
  pb$tick(0)
  PA <- data[perm_index_PA[[i]], ]
  HC <- data[-perm_index_PA[[i]], ]
  r_diff_perm[[i]] <- corr_diff(PA,HC)
  mean_diff[[i]]<-rowMeans(r_diff_perm[[i]])
  pb$tick()
}

p_matrix <- matrix(nrow = nrow(origin_diff_mean), ncol = ncol(origin_diff_mean))


for (i in seq_len(nrow(origin_diff_mean))) {
  
  diff_vector <- sapply(seq_len(n_perm), function(x) mean_diff[[x]][i])
  
  p <- rank(c(origin_diff_mean[i,], diff_vector))[1]/(n_perm + 1)
  
  if (origin_diff_mean[i,] < 0) {
    p <- 1-p
  }
  
  p_matrix[i,] <- p
}

p_matrix<-1-p_matrix
p_matrix
rownames(p_matrix) <- rownames(origin_diff_mean)

p_matrix %>%
  as.data.frame() %>% 
  write.csv('/')





