## NESSE SCRIPT MONTAMOS A ALL.TABLE COM O MÉTODO DO FUJITA. NESSA VERSÃO DO CÓDIGO GARANTIMOS QUE OS GENES SEJAM ABSOLUTAMENTE IGUAIS, VARIANDO UM POUCO ( 10%) OS SNPS. 



#My dir.
setwd("~/Dropbox/laboratorio/scripts/fujita/")

#My packages
library(data.table)
library(splitstackshape)
library(ggplot2)
library(GGally)
library(knitr)
library(dplyr)
library(dtplyr)
library(ggthemes)
library(tidyverse)



# Up main function
source("~/Dropbox/laboratorio/scripts/european_pop/load_functions.R")


############################################################
############################################################

# PREPARANDO OS GENES SOB SELEÇÃO NO CONTINENTE EUROPEU 

############################################################
############################################################

# Preparing positive windows with eur.pop and maf
maf <-fread('~/Dropbox/laboratorio/data/all.pop.maf.tsv')
#distance.eur.win(maf, dist_max = 1125000, dist_min = 1000000)

eur.win(maf, n=150000, m=160000, p=350000)
#eur.win(maf, n=125000, m=130000, p=350000)

rm(maf)  
# Analysis #

## Read fulldata  ##
system.time(
  coding_data <- fread("~/Dropbox/laboratorio/data/coding_data.2.tsv")
)

system.time(
  # conding snps in continental pop
  maf_selinfo <- 
    prepare_coding_data(coding_data)[POP == GRP]
)

maf_selinfo<-as.data.table(mutate(maf_selinfo, htz = 2*(MAF)*(1-MAF)))
maf_selinfo <- mutate(maf_selinfo, Polyphen_Sift = paste(PolyPhenCat, SIFTcat , sep = "_"))

### %in% nos GenesNames
EUR <- filter(maf_selinfo, POP == "EUR")
EAS <- filter(maf_selinfo, POP == "EAS")
AFR <- filter(maf_selinfo, POP == "AFR")
SAS <- filter(maf_selinfo, POP == "SAS")
EUR_EAS <- EUR[EUR$GeneName%in%EAS$GeneName]
AFR_SAS <- AFR[AFR$GeneName%in%SAS$GeneName]
EUR_EAS_AFR_SAS <- EUR_EAS[EUR_EAS$GeneName%in%AFR_SAS$GeneName]
maf_selinfo <- maf_selinfo[maf_selinfo$GeneName%in%EUR_EAS_AFR_SAS$GeneName]
rm(EUR,EAS,AFR,SAS,EUR_EAS,AFR_SAS,EUR_EAS_AFR_SAS)

#rm(maf)
#genes_load, cada linha é um gene, e as colunas sao as informacoes desse gene (pn, ps,...)

# Apenas com europeus
#genes_load <- eval_load(maf_selinfo[POP == "EUR"], groups = c("selection", "GeneName"))[order(selection, pdbpn, n_snp)]

#Com todas pop
genes_load <- eval_load(maf_selinfo, groups = c("selection", "GeneName", "POP"))[order(selection, pdbpn, n_snp)]

#Bad_genes é pra ver quais genes nao estao descritos pelo PolyPhen, pdbpn = Snps_Pol / Pn, nesse script, todos sao bons 
#bad_genes <- genes_load[ selection != 'control' & selection !=  "balanced" & selection !=  "MHC" & selection !=  "peri" & selection !=  "xMHC" & (pdbpn < 0.95 | pdbpn == "NaN"   ), GeneName]
#bad_genes2 <- genes_load[ selection != 'positive'  & selection !=  "balanced" & selection !=  "MHC" & selection !=  "peri" & selection !=  "xMHC" & (pdbpn < 0.95 | pdbpn == "NaN" ), GeneName]


#bad_genes <- genes_load[ selection != 'control' & selection !=  "balanced" & selection !=  "MHC" & selection !=  "peri" & selection !=  "xMHC"  & (pdbpn < 0.95 | sift_pn < 0.95 |pdbpn == "NaN" | sift_pn  == "NaN"), GeneName]
bad_genes <- genes_load[selection != 'control' & (pdbpn < 0.95), GeneName]


# well annotated genes data (good genes)
maf_selinfo <- 
  maf_selinfo[!(GeneName %in% bad_genes)]
maf_selinfo <- filter(maf_selinfo, !POP == "ALL")

############################################################
############################################################

data=maf_selinfo
data=filter(data,POP=="EUR")
data_ori=data


#Arrumando
delta_estat<-as.data.frame(matrix(NA,1,60))
data <- as.data.frame(data)
#parametros     
region = c('positive')
n_genes<- length(unique(data$GeneName[data$selection==region]))
n_SNPs<-  length(data$GeneName[data$selection==region])
n_samples<- 1000L
around<- 0.1
i=1
vetor_SNPs<-c(NA)
vetor_SNPs2<-c(NA)
ptm=proc.time()
#for(i in 1:n_samples){
while (sum(is.na(vetor_SNPs)==F) < n_samples) {
  data<-as.data.frame(data)
  #O sorteio para uma POP:
  genes_sort<-sample(unique(data$GeneName),size = n_genes)
  control_SNP<-sum(data$GeneName%in%genes_sort)
  if(((1-around)*n_SNPs)<=control_SNP & control_SNP<=((1+around)*n_SNPs)){
    vetor_SNPs[i]<-control_SNP
    vetor_SNPs2[i]<-NA
    
    data$selection2 <- data$GeneName%in%genes_sort
    data$selection2[data$selection2==TRUE]<-"positive"
    data$selection2[data$selection2==FALSE]<-"control"
    
    #Calculando estatiscas das duas novas regioes
    data <- as.data.table(data)
    data_r <- data[data$selection2%in%region, ]
    estat_r<-eval_load(data_r, groups = c("POP"))
    data_g <- data[data$selection2%in%"control", ]
    estat_g<-eval_load(data_g, groups = c("POP"))
    
    delta_estat[i,]<-as.numeric(estat_r)-as.numeric(estat_g)
    
  }
  else{vetor_SNPs2[i]<-control_SNP
  vetor_SNPs[i]<-NA}
  i=i+1
}
proc.time()-ptm
colnames(delta_estat)<-colnames(eval_load(data_g, groups = c("POP")))

delta_estat<-delta_estat[is.na(delta_estat$pn)==F,]

#### DADO REAL POSITIVO 

data_ori<-as.data.table(data_ori)
real<-eval_load(data_ori, groups = c("POP","selection"))
real<-real[,3:length(real)]
delta_real<-as.numeric(real[2,])-as.numeric(real[1,])

#teste <- select(eval_load(data_g, groups = c("POP")), -POP)
#teste <- eval_load(data_g, groups = c("POP"))
teste <-delta_estat[,-1]
delta_real<-t(as.data.frame(delta_real))
colnames(delta_real)<-colnames(teste)

teste <- melt(teste)
teste$treatment <- rep("control")
teste$POP <- rep("EUR")
teste$POP_US <- rep("EUR")
teste$boot <- rep("std_gene")

delta_real <- melt(delta_real)
delta_real$treatment <- rep("positive")
delta_real$POP <- rep("EUR")
delta_real$POP_US <- rep("EUR")
delta_real <- select(delta_real, - Var1)
delta_real <- rename(delta_real, variable = Var2 )
delta_real$boot <- rep("std_gene")



Delta_eur <- rbind(delta_real, teste)
Delta_eur <- as.data.table(Delta_eur)
Delta_eur[, pvalue := (.N - frank(value) + 1)/ .N, 
          by = .(POP, variable, POP_US)]
Delta_eur[, plabelpos := max(value) + 0.15 * (max(value) - min(value)), by = .(variable)]
Delta_eur <- mutate(Delta_eur, is_significant = c(pvalue <= 0.05 | pvalue >= 1 - 0.05))

############################################################
############################################################

data=maf_selinfo
data=filter(data,POP=="EAS")
data_ori=data


#Arrumando
delta_estat<-as.data.frame(matrix(NA,1,60))
data <- as.data.frame(data)
#parametros     
region = c('positive')
n_genes<- length(unique(data$GeneName[data$selection==region]))
n_SNPs<-  length(data$GeneName[data$selection==region])
n_samples<- 1000L
around<- 0.1
i=1
vetor_SNPs<-c(NA)
vetor_SNPs2<-c(NA)
ptm=proc.time()
#for(i in 1:n_samples){
while (sum(is.na(vetor_SNPs)==F) < n_samples) {
  data<-as.data.frame(data)
  #O sorteio para uma POP:
  genes_sort<-sample(unique(data$GeneName),size = n_genes)
  control_SNP<-sum(data$GeneName%in%genes_sort)
  if(((1-around)*n_SNPs)<=control_SNP & control_SNP<=((1+around)*n_SNPs)){
    vetor_SNPs[i]<-control_SNP
    vetor_SNPs2[i]<-NA
    
    data$selection2 <- data$GeneName%in%genes_sort
    data$selection2[data$selection2==TRUE]<-"positive"
    data$selection2[data$selection2==FALSE]<-"control"
    
    #Calculando estatiscas das duas novas regioes
    data <- as.data.table(data)
    data_r <- data[data$selection2%in%region, ]
    estat_r<-eval_load(data_r, groups = c("POP"))
    data_g <- data[data$selection2%in%"control", ]
    estat_g<-eval_load(data_g, groups = c("POP"))
    
    delta_estat[i,]<-as.numeric(estat_r)-as.numeric(estat_g)
    
  }
  else{vetor_SNPs2[i]<-control_SNP
  vetor_SNPs[i]<-NA}
  i=i+1
}
proc.time()-ptm
colnames(delta_estat)<-colnames(eval_load(data_g, groups = c("POP")))

delta_estat<-delta_estat[is.na(delta_estat$pn)==F,]

#### DADO REAL POSITIVO 

data_ori<-as.data.table(data_ori)
real<-eval_load(data_ori, groups = c("POP","selection"))
real<-real[,3:length(real)]
delta_real<-as.numeric(real[2,])-as.numeric(real[1,])

#teste <- select(eval_load(data_g, groups = c("POP")), -POP)
#teste <- eval_load(data_g, groups = c("POP"))
teste <-delta_estat[,-1]
delta_real<-t(as.data.frame(delta_real))
colnames(delta_real)<-colnames(teste)

teste <- melt(teste)
teste$treatment <- rep("control")
teste$POP <- rep("EAS")
teste$POP_US <- rep("EUR")
teste$boot <- rep("std_gene")

delta_real <- melt(delta_real)
delta_real$treatment <- rep("positive")
delta_real$POP <- rep("EAS")
delta_real$POP_US <- rep("EUR")
delta_real <- select(delta_real, - Var1)
delta_real <- rename(delta_real, variable = Var2 )
delta_real$boot <- rep("std_gene")



Delta_eas <- rbind(delta_real, teste)
Delta_eas <- as.data.table(Delta_eas)
Delta_eas[, pvalue := (.N - frank(value) + 1)/ .N, 
          by = .(POP, variable, POP_US)]
Delta_eas[, plabelpos := max(value) + 0.15 * (max(value) - min(value)), by = .(variable)]
Delta_eas <- mutate(Delta_eas, is_significant = c(pvalue <= 0.05 | pvalue >= 1 - 0.05))

############################################################
############################################################

data=maf_selinfo
data=filter(data,POP=="SAS")
data_ori=data


#Arrumando
delta_estat<-as.data.frame(matrix(NA,1,60))
data <- as.data.frame(data)
#parametros     
region = c('positive')
n_genes<- length(unique(data$GeneName[data$selection==region]))
n_SNPs<-  length(data$GeneName[data$selection==region])
n_samples<- 1000L
around<- 0.1
i=1
vetor_SNPs<-c(NA)
vetor_SNPs2<-c(NA)
ptm=proc.time()
#for(i in 1:n_samples){
while (sum(is.na(vetor_SNPs)==F) < n_samples) {
  data<-as.data.frame(data)
  #O sorteio para uma POP:
  genes_sort<-sample(unique(data$GeneName),size = n_genes)
  control_SNP<-sum(data$GeneName%in%genes_sort)
  if(((1-around)*n_SNPs)<=control_SNP & control_SNP<=((1+around)*n_SNPs)){
    vetor_SNPs[i]<-control_SNP
    vetor_SNPs2[i]<-NA
    
    data$selection2 <- data$GeneName%in%genes_sort
    data$selection2[data$selection2==TRUE]<-"positive"
    data$selection2[data$selection2==FALSE]<-"control"
    
    #Calculando estatiscas das duas novas regioes
    data <- as.data.table(data)
    data_r <- data[data$selection2%in%region, ]
    estat_r<-eval_load(data_r, groups = c("POP"))
    data_g <- data[data$selection2%in%"control", ]
    estat_g<-eval_load(data_g, groups = c("POP"))
    
    delta_estat[i,]<-as.numeric(estat_r)-as.numeric(estat_g)
    
  }
  else{vetor_SNPs2[i]<-control_SNP
  vetor_SNPs[i]<-NA}
  i=i+1
}
proc.time()-ptm
colnames(delta_estat)<-colnames(eval_load(data_g, groups = c("POP")))

delta_estat<-delta_estat[is.na(delta_estat$pn)==F,]

#### DADO REAL POSITIVO 

data_ori<-as.data.table(data_ori)
real<-eval_load(data_ori, groups = c("POP","selection"))
real<-real[,3:length(real)]
delta_real<-as.numeric(real[2,])-as.numeric(real[1,])

#teste <- select(eval_load(data_g, groups = c("POP")), -POP)
#teste <- eval_load(data_g, groups = c("POP"))
teste <-delta_estat[,-1]
delta_real<-t(as.data.frame(delta_real))
colnames(delta_real)<-colnames(teste)

teste <- melt(teste)
teste$treatment <- rep("control")
teste$POP <- rep("SAS")
teste$POP_US <- rep("EUR")
teste$boot <- rep("std_gene")

delta_real <- melt(delta_real)
delta_real$treatment <- rep("positive")
delta_real$POP <- rep("SAS")
delta_real$POP_US <- rep("EUR")
delta_real <- select(delta_real, - Var1)
delta_real <- rename(delta_real, variable = Var2 )
delta_real$boot <- rep("std_gene")



Delta_sas <- rbind(delta_real, teste)
Delta_sas <- as.data.table(Delta_sas)
Delta_sas[, pvalue := (.N - frank(value) + 1)/ .N, 
          by = .(POP, variable, POP_US)]
Delta_sas[, plabelpos := max(value) + 0.15 * (max(value) - min(value)), by = .(variable)]
Delta_sas <- mutate(Delta_sas, is_significant = c(pvalue <= 0.05 | pvalue >= 1 - 0.05))

############################################################
############################################################

data=maf_selinfo
data=filter(data,POP=="AFR")
data_ori=data


#Arrumando
delta_estat<-as.data.frame(matrix(NA,1,60))
data <- as.data.frame(data)
#parametros     
region = c('positive')
n_genes<- length(unique(data$GeneName[data$selection==region]))
n_SNPs<-  length(data$GeneName[data$selection==region])
n_samples<- 1000L
around<- 0.1
i=1
vetor_SNPs<-c(NA)
vetor_SNPs2<-c(NA)
ptm=proc.time()
#for(i in 1:n_samples){
while (sum(is.na(vetor_SNPs)==F) < n_samples) {
  data<-as.data.frame(data)
  #O sorteio para uma POP:
  genes_sort<-sample(unique(data$GeneName),size = n_genes)
  control_SNP<-sum(data$GeneName%in%genes_sort)
  if(((1-around)*n_SNPs)<=control_SNP & control_SNP<=((1+around)*n_SNPs)){
    vetor_SNPs[i]<-control_SNP
    vetor_SNPs2[i]<-NA
    
    data$selection2 <- data$GeneName%in%genes_sort
    data$selection2[data$selection2==TRUE]<-"positive"
    data$selection2[data$selection2==FALSE]<-"control"
    
    #Calculando estatiscas das duas novas regioes
    data <- as.data.table(data)
    data_r <- data[data$selection2%in%region, ]
    estat_r<-eval_load(data_r, groups = c("POP"))
    data_g <- data[data$selection2%in%"control", ]
    estat_g<-eval_load(data_g, groups = c("POP"))
    
    delta_estat[i,]<-as.numeric(estat_r)-as.numeric(estat_g)
    
  }
  else{vetor_SNPs2[i]<-control_SNP
  vetor_SNPs[i]<-NA}
  i=i+1
}
proc.time()-ptm
colnames(delta_estat)<-colnames(eval_load(data_g, groups = c("POP")))

delta_estat<-delta_estat[is.na(delta_estat$pn)==F,]

#### DADO REAL POSITIVO 

data_ori<-as.data.table(data_ori)
real<-eval_load(data_ori, groups = c("POP","selection"))
real<-real[,3:length(real)]
delta_real<-as.numeric(real[2,])-as.numeric(real[1,])

#teste <- select(eval_load(data_g, groups = c("POP")), -POP)
#teste <- eval_load(data_g, groups = c("POP"))
teste <-delta_estat[,-1]
delta_real<-t(as.data.frame(delta_real))
colnames(delta_real)<-colnames(teste)

teste <- melt(teste)
teste$treatment <- rep("control")
teste$POP <- rep("AFR")
teste$POP_US <- rep("EUR")
teste$boot <- rep("std_gene")

delta_real <- melt(delta_real)
delta_real$treatment <- rep("positive")
delta_real$POP <- rep("AFR")
delta_real$POP_US <- rep("EUR")
delta_real <- select(delta_real, - Var1)
delta_real <- rename(delta_real, variable = Var2 )
delta_real$boot <- rep("std_gene")



Delta_afr <- rbind(delta_real, teste)
Delta_afr <- as.data.table(Delta_afr)
Delta_afr[, pvalue := (.N - frank(value) + 1)/ .N, 
          by = .(POP, variable, POP_US)]
Delta_afr[, plabelpos := max(value) + 0.15 * (max(value) - min(value)), by = .(variable)]
Delta_afr <- mutate(Delta_afr, is_significant = c(pvalue <= 0.05 | pvalue >= 1 - 0.05))

############################################################
############################################################

eur.all.table <- rbind(Delta_afr, Delta_eur, Delta_eas, Delta_sas)

############################################################
############################################################

# PREPARANDO OS GENES SOB SELEÇÃO NO CONTINENTE ASIÁTICO 

############################################################
############################################################



# Preparing positive windows with eur.pop and maf
maf <-fread('~/Dropbox/laboratorio/data/all.pop.maf.tsv')
#asn.win(maf, n=150000, m=170000, p=350000)
asn.win(maf, n=235000, m=245000, p=350000)

rm(maf)
# Analysis #

## Read fulldata  ##
system.time(
  coding_data <- fread("~/Dropbox/laboratorio/data/coding_data.2.tsv")
)

system.time(
  # conding snps in continental pop
  maf_selinfo <- 
    prepare_coding_data(coding_data)[POP == GRP]
)

#maf_selinfo <- filter(maf_selinfo, POP == "EAS")

maf_selinfo<-as.data.table(mutate(maf_selinfo, htz = 2*(MAF)*(1-MAF)))

### %in% nos GenesNames
EUR <- filter(maf_selinfo, POP == "EUR")
EAS <- filter(maf_selinfo, POP == "EAS")
AFR <- filter(maf_selinfo, POP == "AFR")
SAS <- filter(maf_selinfo, POP == "SAS")
EUR_EAS <- EUR[EUR$GeneName%in%EAS$GeneName]
AFR_SAS <- AFR[AFR$GeneName%in%SAS$GeneName]
EUR_EAS_AFR_SAS <- EUR_EAS[EUR_EAS$GeneName%in%AFR_SAS$GeneName]
maf_selinfo <- maf_selinfo[maf_selinfo$GeneName%in%EUR_EAS_AFR_SAS$GeneName]
rm(EUR,EAS,AFR,SAS,EUR_EAS,AFR_SAS,EUR_EAS_AFR_SAS)

#Com todas pop
genes_load <- eval_load(maf_selinfo, groups = c("selection", "GeneName", "POP"))[order(selection, pdbpn, n_snp)]

bad_genes <- genes_load[selection != 'control' & (pdbpn < 0.9 |pdbpn == "NaN"), GeneName]

# well annotated genes data (good genes)
maf_selinfo <- 
  maf_selinfo[!(GeneName %in% bad_genes)]
table(maf_selinfo$selection[maf_selinfo$POP=='EAS'])
maf_selinfo <- filter(maf_selinfo, !POP == "ALL")


############################################################
############################################################

data=maf_selinfo
data=filter(data,POP=="EAS")
data_ori=data


#Arrumando
delta_estat<-as.data.frame(matrix(NA,1,60))
data <- as.data.frame(data)
#parametros     
region = c('positive')
n_genes<- length(unique(data$GeneName[data$selection==region]))
n_SNPs<-  length(data$GeneName[data$selection==region])
n_samples<- 1000L
around<- 0.1
i=1
vetor_SNPs<-c(NA)
vetor_SNPs2<-c(NA)
ptm=proc.time()
#for(i in 1:n_samples){
while (sum(is.na(vetor_SNPs)==F) < n_samples) {
  data<-as.data.frame(data)
  #O sorteio para uma POP:
  genes_sort<-sample(unique(data$GeneName),size = n_genes)
  control_SNP<-sum(data$GeneName%in%genes_sort)
  if(((1-around)*n_SNPs)<=control_SNP & control_SNP<=((1+around)*n_SNPs)){
    vetor_SNPs[i]<-control_SNP
    vetor_SNPs2[i]<-NA
    
    data$selection2 <- data$GeneName%in%genes_sort
    data$selection2[data$selection2==TRUE]<-"positive"
    data$selection2[data$selection2==FALSE]<-"control"
    
    #Calculando estatiscas das duas novas regioes
    data <- as.data.table(data)
    data_r <- data[data$selection2%in%region, ]
    estat_r<-eval_load(data_r, groups = c("POP"))
    data_g <- data[data$selection2%in%"control", ]
    estat_g<-eval_load(data_g, groups = c("POP"))
    
    delta_estat[i,]<-as.numeric(estat_r)-as.numeric(estat_g)
    
  }
  else{vetor_SNPs2[i]<-control_SNP
  vetor_SNPs[i]<-NA}
  i=i+1
}
proc.time()-ptm
colnames(delta_estat)<-colnames(eval_load(data_g, groups = c("POP")))

delta_estat<-delta_estat[is.na(delta_estat$pn)==F,]

#### DADO REAL POSITIVO 

data_ori<-as.data.table(data_ori)
real<-eval_load(data_ori, groups = c("POP","selection"))
real<-real[,3:length(real)]
delta_real<-as.numeric(real[2,])-as.numeric(real[1,])

#teste <- select(eval_load(data_g, groups = c("POP")), -POP)
#teste <- eval_load(data_g, groups = c("POP"))
teste <-delta_estat[,-1]
delta_real<-t(as.data.frame(delta_real))
colnames(delta_real)<-colnames(teste)

teste <- melt(teste)
teste$treatment <- rep("control")
teste$POP <- rep("EAS")
teste$POP_US <- rep("EAS")
teste$boot <- rep("std_gene")

delta_real <- melt(delta_real)
delta_real$treatment <- rep("positive")
delta_real$POP <- rep("EAS")
delta_real$POP_US <- rep("EAS")
delta_real <- select(delta_real, - Var1)
delta_real <- rename(delta_real, variable = Var2 )
delta_real$boot <- rep("std_gene")



Delta_eur <- rbind(delta_real, teste)
Delta_eur <- as.data.table(Delta_eur)
Delta_eur[, pvalue := (.N - frank(value) + 1)/ .N, 
          by = .(POP, variable, POP_US)]
Delta_eur[, plabelpos := max(value) + 0.15 * (max(value) - min(value)), by = .(variable)]
Delta_eur <- mutate(Delta_eur, is_significant = c(pvalue <= 0.05 | pvalue >= 1 - 0.05))

############################################################
############################################################

data=maf_selinfo
data=filter(data,POP=="EUR")
data_ori=data


#Arrumando
delta_estat<-as.data.frame(matrix(NA,1,60))
data <- as.data.frame(data)
#parametros     
region = c('positive')
n_genes<- length(unique(data$GeneName[data$selection==region]))
n_SNPs<-  length(data$GeneName[data$selection==region])
n_samples<- 1000L
around<- 0.1
i=1
vetor_SNPs<-c(NA)
vetor_SNPs2<-c(NA)
ptm=proc.time()
#for(i in 1:n_samples){
while (sum(is.na(vetor_SNPs)==F) < n_samples) {
  data<-as.data.frame(data)
  #O sorteio para uma POP:
  genes_sort<-sample(unique(data$GeneName),size = n_genes)
  control_SNP<-sum(data$GeneName%in%genes_sort)
  if(((1-around)*n_SNPs)<=control_SNP & control_SNP<=((1+around)*n_SNPs)){
    vetor_SNPs[i]<-control_SNP
    vetor_SNPs2[i]<-NA
    
    data$selection2 <- data$GeneName%in%genes_sort
    data$selection2[data$selection2==TRUE]<-"positive"
    data$selection2[data$selection2==FALSE]<-"control"
    
    #Calculando estatiscas das duas novas regioes
    data <- as.data.table(data)
    data_r <- data[data$selection2%in%region, ]
    estat_r<-eval_load(data_r, groups = c("POP"))
    data_g <- data[data$selection2%in%"control", ]
    estat_g<-eval_load(data_g, groups = c("POP"))
    
    delta_estat[i,]<-as.numeric(estat_r)-as.numeric(estat_g)
    
  }
  else{vetor_SNPs2[i]<-control_SNP
  vetor_SNPs[i]<-NA}
  i=i+1
}
proc.time()-ptm
colnames(delta_estat)<-colnames(eval_load(data_g, groups = c("POP")))

delta_estat<-delta_estat[is.na(delta_estat$pn)==F,]

#### DADO REAL POSITIVO 

data_ori<-as.data.table(data_ori)
real<-eval_load(data_ori, groups = c("POP","selection"))
real<-real[,3:length(real)]
delta_real<-as.numeric(real[2,])-as.numeric(real[1,])

#teste <- select(eval_load(data_g, groups = c("POP")), -POP)
#teste <- eval_load(data_g, groups = c("POP"))
teste <-delta_estat[,-1]
delta_real<-t(as.data.frame(delta_real))
colnames(delta_real)<-colnames(teste)

teste <- melt(teste)
teste$treatment <- rep("control")
teste$POP <- rep("EUR")
teste$POP_US <- rep("EAS")
teste$boot <- rep("std_gene")

delta_real <- melt(delta_real)
delta_real$treatment <- rep("positive")
delta_real$POP <- rep("EUR")
delta_real$POP_US <- rep("EAS")
delta_real <- select(delta_real, - Var1)
delta_real <- rename(delta_real, variable = Var2 )
delta_real$boot <- rep("std_gene")



Delta_eas <- rbind(delta_real, teste)
Delta_eas <- as.data.table(Delta_eas)
Delta_eas[, pvalue := (.N - frank(value) + 1)/ .N, 
          by = .(POP, variable, POP_US)]
Delta_eas[, plabelpos := max(value) + 0.15 * (max(value) - min(value)), by = .(variable)]
Delta_eas <- mutate(Delta_eas, is_significant = c(pvalue <= 0.05 | pvalue >= 1 - 0.05))

############################################################
############################################################

data=maf_selinfo
data=filter(data,POP=="SAS")
data_ori=data


#Arrumando
delta_estat<-as.data.frame(matrix(NA,1,60))
data <- as.data.frame(data)
#parametros     
region = c('positive')
n_genes<- length(unique(data$GeneName[data$selection==region]))
n_SNPs<-  length(data$GeneName[data$selection==region])
n_samples<- 1000L
around<- 0.1
i=1
vetor_SNPs<-c(NA)
vetor_SNPs2<-c(NA)
ptm=proc.time()
#for(i in 1:n_samples){
while (sum(is.na(vetor_SNPs)==F) < n_samples) {
  data<-as.data.frame(data)
  #O sorteio para uma POP:
  genes_sort<-sample(unique(data$GeneName),size = n_genes)
  control_SNP<-sum(data$GeneName%in%genes_sort)
  if(((1-around)*n_SNPs)<=control_SNP & control_SNP<=((1+around)*n_SNPs)){
    vetor_SNPs[i]<-control_SNP
    vetor_SNPs2[i]<-NA
    
    data$selection2 <- data$GeneName%in%genes_sort
    data$selection2[data$selection2==TRUE]<-"positive"
    data$selection2[data$selection2==FALSE]<-"control"
    
    #Calculando estatiscas das duas novas regioes
    data <- as.data.table(data)
    data_r <- data[data$selection2%in%region, ]
    estat_r<-eval_load(data_r, groups = c("POP"))
    data_g <- data[data$selection2%in%"control", ]
    estat_g<-eval_load(data_g, groups = c("POP"))
    
    delta_estat[i,]<-as.numeric(estat_r)-as.numeric(estat_g)
    
  }
  else{vetor_SNPs2[i]<-control_SNP
  vetor_SNPs[i]<-NA}
  i=i+1
}
proc.time()-ptm
colnames(delta_estat)<-colnames(eval_load(data_g, groups = c("POP")))

delta_estat<-delta_estat[is.na(delta_estat$pn)==F,]

#### DADO REAL POSITIVO 

data_ori<-as.data.table(data_ori)
real<-eval_load(data_ori, groups = c("POP","selection"))
real<-real[,3:length(real)]
delta_real<-as.numeric(real[2,])-as.numeric(real[1,])

#teste <- select(eval_load(data_g, groups = c("POP")), -POP)
#teste <- eval_load(data_g, groups = c("POP"))
teste <-delta_estat[,-1]
delta_real<-t(as.data.frame(delta_real))
colnames(delta_real)<-colnames(teste)

teste <- melt(teste)
teste$treatment <- rep("control")
teste$POP <- rep("SAS")
teste$POP_US <- rep("EAS")
teste$boot <- rep("std_gene")

delta_real <- melt(delta_real)
delta_real$treatment <- rep("positive")
delta_real$POP <- rep("SAS")
delta_real$POP_US <- rep("EAS")
delta_real <- select(delta_real, - Var1)
delta_real <- rename(delta_real, variable = Var2 )
delta_real$boot <- rep("std_gene")



Delta_sas <- rbind(delta_real, teste)
Delta_sas <- as.data.table(Delta_sas)
Delta_sas[, pvalue := (.N - frank(value) + 1)/ .N, 
          by = .(POP, variable, POP_US)]
Delta_sas[, plabelpos := max(value) + 0.15 * (max(value) - min(value)), by = .(variable)]
Delta_sas <- mutate(Delta_sas, is_significant = c(pvalue <= 0.05 | pvalue >= 1 - 0.05))

############################################################
############################################################

data=maf_selinfo
data=filter(data,POP=="AFR")
data_ori=data


#Arrumando
delta_estat<-as.data.frame(matrix(NA,1,60))
data <- as.data.frame(data)
#parametros     
region = c('positive')
n_genes<- length(unique(data$GeneName[data$selection==region]))
n_SNPs<-  length(data$GeneName[data$selection==region])
n_samples<- 1000L
around<- 0.1
i=1
vetor_SNPs<-c(NA)
vetor_SNPs2<-c(NA)
ptm=proc.time()
#for(i in 1:n_samples){
while (sum(is.na(vetor_SNPs)==F) < n_samples) {
  data<-as.data.frame(data)
  #O sorteio para uma POP:
  genes_sort<-sample(unique(data$GeneName),size = n_genes)
  control_SNP<-sum(data$GeneName%in%genes_sort)
  if(((1-around)*n_SNPs)<=control_SNP & control_SNP<=((1+around)*n_SNPs)){
    vetor_SNPs[i]<-control_SNP
    vetor_SNPs2[i]<-NA
    
    data$selection2 <- data$GeneName%in%genes_sort
    data$selection2[data$selection2==TRUE]<-"positive"
    data$selection2[data$selection2==FALSE]<-"control"
    
    #Calculando estatiscas das duas novas regioes
    data <- as.data.table(data)
    data_r <- data[data$selection2%in%region, ]
    estat_r<-eval_load(data_r, groups = c("POP"))
    data_g <- data[data$selection2%in%"control", ]
    estat_g<-eval_load(data_g, groups = c("POP"))
    
    delta_estat[i,]<-as.numeric(estat_r)-as.numeric(estat_g)
    
  }
  else{vetor_SNPs2[i]<-control_SNP
  vetor_SNPs[i]<-NA}
  i=i+1
}
proc.time()-ptm
colnames(delta_estat)<-colnames(eval_load(data_g, groups = c("POP")))

delta_estat<-delta_estat[is.na(delta_estat$pn)==F,]

#### DADO REAL POSITIVO 

data_ori<-as.data.table(data_ori)
real<-eval_load(data_ori, groups = c("POP","selection"))
real<-real[,3:length(real)]
delta_real<-as.numeric(real[2,])-as.numeric(real[1,])

#teste <- select(eval_load(data_g, groups = c("POP")), -POP)
#teste <- eval_load(data_g, groups = c("POP"))
teste <-delta_estat[,-1]
delta_real<-t(as.data.frame(delta_real))
colnames(delta_real)<-colnames(teste)

teste <- melt(teste)
teste$treatment <- rep("control")
teste$POP <- rep("AFR")
teste$POP_US <- rep("EAS")
teste$boot <- rep("std_gene")

delta_real <- melt(delta_real)
delta_real$treatment <- rep("positive")
delta_real$POP <- rep("AFR")
delta_real$POP_US <- rep("EAS")
delta_real <- select(delta_real, - Var1)
delta_real <- rename(delta_real, variable = Var2 )
delta_real$boot <- rep("std_gene")



Delta_afr <- rbind(delta_real, teste)
Delta_afr <- as.data.table(Delta_afr)
Delta_afr[, pvalue := (.N - frank(value) + 1)/ .N, 
          by = .(POP, variable, POP_US)]
Delta_afr[, plabelpos := max(value) + 0.15 * (max(value) - min(value)), by = .(variable)]
Delta_afr <- mutate(Delta_afr, is_significant = c(pvalue <= 0.05 | pvalue >= 1 - 0.05))

############################################################
############################################################

eas.all.table <- rbind(Delta_afr, Delta_eur, Delta_eas, Delta_sas)

############################################################
############################################################