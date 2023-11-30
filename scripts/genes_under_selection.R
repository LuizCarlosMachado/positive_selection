### EUROPA 
EUR_genes <- fread("~/Dropbox/laboratorio/data/Colonna_HighD.csv")
EUR_genes <- filter(EUR_genes , VT == 'SNP')
EUR_genes <- filter(EUR_genes , PAIR == 'AFR-EUR')
EUR_genes <- filter(EUR_genes , POPwithHighestDAF == 'EUR')
EUR_genes <- filter(EUR_genes, GeneBiotype == "protein_coding")
EUR_genes <- select(EUR_genes, c(CHR, POPwithHighestDAF, HGNCsymbol, POS, GeneStart_bp, GeneEnd_bp ))
EUR_genes <- mutate(EUR_genes, GeneName = HGNCsymbol )

filter(maf_selinfo,  CHR == 14 & (POS > 62037258 & POS <62125414))
filter(maf_selinfo,  CHR == 1 & (POS > 116915290 & POS <116952883))

### ÁSIA 

EAS_genes <- fread("~/Dropbox/laboratorio/data/colonna/aba2_EAS.csv")
EAS_genes <- filter(EAS_genes , VT == 'SNP')
EAS_genes <- filter(EAS_genes, GeneBiotype == "protein_coding")
EAS_genes <- filter(EAS_genes , PAIR == 'AFR-EAS')
EAS_genes <- filter(EAS_genes, highestDAF >  925)
EAS_genes <- filter(EAS_genes , POPwithHighestDAF == 'EAS')
EAS_genes <- select(EAS_genes, c(CHR, POPwithHighestDAF, HGNCsymbol, POS, GeneStart_bp, GeneEnd_bp ))
EAS_genes <- mutate(EAS_genes, GeneName = HGNCsymbol )




positive_genes <- rbind(EUR_genes, EAS_genes)
### AKEY 2009 
AKEY <-  fread("~/Dropbox/laboratorio/data/akey_2009_genes.csv")
AKEY<-as.data.frame(AKEY)
#####GABI

for(i in 1:nrow(positive_genes)){
   if(nrow(AKEY[AKEY$CHR == positive_genes$CHR[i] & (AKEY$GeneStart_bp < positive_genes$POS[i]) & (positive_genes$POS[i] < AKEY$GeneEnd_bp),])>0){
      positive_genes$Akey[i]<-
           AKEY[AKEY$CHR == positive_genes$CHR[i] & (AKEY$GeneStart_bp < positive_genes$POS[i]) & (positive_genes$POS[i] < AKEY$GeneEnd_bp),5]
   }
   else{
      positive_genes$Akey[i]<-NA
   }
}

sum(is.na(positive_genes$Akey)==F)


GROSSMAN.2014 <- fread("~/Dropbox/laboratorio/data/grossman_2014_genes.csv")
GROSSMAN.2014<-as.data.frame(GROSSMAN.2014)

for(i in 1:nrow(positive_genes)){
   if(nrow(GROSSMAN.2014[GROSSMAN.2014$CHR == positive_genes$CHR[i] & (GROSSMAN.2014$POS_start_hg18 < positive_genes$POS[i]) & (positive_genes$POS[i] < GROSSMAN.2014$POS_end_hg18),])>0){
      positive_genes$grossman2014[i]<-
           GROSSMAN.2014[GROSSMAN.2014$CHR == positive_genes$CHR[i] & (GROSSMAN.2014$POS_start_hg18 < positive_genes$POS[i]) & (positive_genes$POS[i] < GROSSMAN.2014$POS_end_hg18),4]
   }
   else{
      positive_genes$grossman2014[i]<-NA
   }
}


####gene ao inves de SNP

#Grossman
for(i in 1:nrow(positive_genes)){
   if(nrow(GROSSMAN.2014[(GROSSMAN.2014$CHR == positive_genes$CHR[i] & (GROSSMAN.2014$POS_start_hg18 > positive_genes$GeneStart_bp[i]) & (positive_genes$GeneEnd_bp[i] > GROSSMAN.2014$POS_start_hg18))|
                         (GROSSMAN.2014$CHR == positive_genes$CHR[i] & (GROSSMAN.2014$POS_end_hg18 > positive_genes$GeneStart_bp[i]) & (positive_genes$GeneEnd_bp[i] > GROSSMAN.2014$POS_end_hg18))                  ,])>0){
      positive_genes$grossman2014genes[i]<-
           GROSSMAN.2014[(GROSSMAN.2014$CHR == positive_genes$CHR[i] & (GROSSMAN.2014$POS_start_hg18 > positive_genes$GeneStart_bp[i]) & (positive_genes$GeneEnd_bp[i] > GROSSMAN.2014$POS_start_hg18))|
                         (GROSSMAN.2014$CHR == positive_genes$CHR[i] & (GROSSMAN.2014$POS_end_hg18 > positive_genes$GeneStart_bp[i]) & (positive_genes$GeneEnd_bp[i] > GROSSMAN.2014$POS_end_hg18))                  ,4]
   }
   #else{
   #   positive_genes$grossman2014genes[i]<-NA
   #}
}



#vetor<-c(NA)
#for(i in 1:nrow(positive_genes)){
#   vetor[i]<-nrow(AKEY[(AKEY$CHR == positive_genes$CHR[i] & (AKEY$GeneStart_bp > positive_genes$GeneStart_bp[i]) & (positive_genes$GeneEnd_bp[i] > AKEY$GeneStart_bp))|
#                          (AKEY$CHR == positive_genes$CHR[i] & (AKEY$GeneEnd_bp > positive_genes$GeneStart_bp[i]) & (positive_genes$GeneEnd_bp[i] > AKEY$GeneEnd_bp))                  ,])}
#which(vetor==2)
#Akey
for(i in 1:nrow(positive_genes)){
   if(nrow(AKEY[(AKEY$CHR == positive_genes$CHR[i] & (AKEY$GeneStart_bp > positive_genes$GeneStart_bp[i]) & (positive_genes$GeneEnd_bp[i] > AKEY$GeneStart_bp))|
                         (AKEY$CHR == positive_genes$CHR[i] & (AKEY$GeneEnd_bp > positive_genes$GeneStart_bp[i]) & (positive_genes$GeneEnd_bp[i] > AKEY$GeneEnd_bp))                  ,])>0){
      positive_genes$Akeygenes[i]<-
           AKEY[(AKEY$CHR == positive_genes$CHR[i] & (AKEY$GeneStart_bp > positive_genes$GeneStart_bp[i]) & (positive_genes$GeneEnd_bp[i] > AKEY$GeneStart_bp))|
                         (AKEY$CHR == positive_genes$CHR[i] & (AKEY$GeneEnd_bp > positive_genes$GeneStart_bp[i]) & (positive_genes$GeneEnd_bp[i] > AKEY$GeneEnd_bp))                  ,5]
   }
   #else{
   #   positive_genes$Akeygenes[i]<-NA
   #}
}
#i = 17 e 24 1 da problema!
sum(is.na(positive_genes$Akeygenes)==F)
