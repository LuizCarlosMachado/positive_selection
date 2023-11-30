#My dir.
setwd("~/Dropbox/laboratorio/scripts/european_pop/")

#My packages
library(data.table)
library(splitstackshape)
library(ggplot2)
library(GGally)
library(knitr)
library(dplyr)
library(dtplyr)
library(ggthemes)

# Up main function
source("~/Dropbox/laboratorio/scripts/european_pop/load_functions.R")

# Preparing positive windows with eur.pop and maf
maf <-fread('~/Dropbox/laboratorio/data/all.pop.maf.tsv')
maf <- as.data.frame(maf)
#distance.eur.win(maf, dist_max = 1125000, dist_min = 1000000)

#eur.win(maf, n=160000, m=170000, p=350000)
eur.win(maf, n=125000, m=130000, p=350000)

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

maf_selinfo <- as.data.frame(maf_selinfo)
maf_selinfo <- mutate(maf_selinfo, htz = 2*(MAF)*(1-MAF))
maf_selinfo <- mutate(maf_selinfo, Polyphen_Sift = paste(PolyPhenCat, SIFTcat , sep = "_"))

### %in% nos GenesNames
EUR <- filter(maf_selinfo, POP == "EUR")
EAS <- filter(maf_selinfo, POP == "EAS")
AFR <- filter(maf_selinfo, POP == "AFR")
SAS <- filter(maf_selinfo, POP == "SAS")

EUR_EAS <- EUR %>% filter(GeneName %in% pull(EAS,GeneName))
AFR_SAS <- AFR %>% filter(GeneName %in% pull(SAS,GeneName))
EUR_EAS_AFR_SAS <- EUR_EAS %>% filter(GeneName %in% pull(AFR_SAS,GeneName))
maf_selinfo <- maf_selinfo %>% filter(GeneName %in% pull(EUR_EAS_AFR_SAS,GeneName))
rm(EUR,EAS,AFR,SAS,EUR_EAS,AFR_SAS,EUR_EAS_AFR_SAS)

#rm(maf)
#genes_load, cada linha é um gene, e as colunas sao as informacoes desse gene (pn, ps,...)

# Apenas com europeus
#genes_load <- eval_load(maf_selinfo[POP == "EUR"], groups = c("selection", "GeneName"))[order(selection, pdbpn, n_snp)]

#Com todas pop
maf_selinfo <- as.data.table(maf_selinfo)
genes_load <- eval_load(maf_selinfo, groups = c("selection", "GeneName", "POP"))[order(selection, pdbpn, n_snp)]
maf_selinfo <- as.data.frame(maf_selinfo)

#Bad_genes é pra ver quais genes nao estao descritos pelo PolyPhen, pdbpn = Snps_Pol / Pn, nesse script, todos sao bons 
#bad_genes <- genes_load[ selection != 'control' & selection !=  "balanced" & selection !=  "MHC" & selection !=  "peri" & selection !=  "xMHC" & (pdbpn < 0.95 | pdbpn == "NaN"   ), GeneName]
#bad_genes2 <- genes_load[ selection != 'positive'  & selection !=  "balanced" & selection !=  "MHC" & selection !=  "peri" & selection !=  "xMHC" & (pdbpn < 0.95 | pdbpn == "NaN" ), GeneName]


#bad_genes <- genes_load[ selection != 'control' & selection !=  "balanced" & selection !=  "MHC" & selection !=  "peri" & selection !=  "xMHC"  & (pdbpn < 0.95 | sift_pn < 0.95 |pdbpn == "NaN" | sift_pn  == "NaN"), GeneName]
bad_genes <- genes_load[selection != 'control' & (pdbpn < 0.95), GeneName]


# well annotated genes data (good genes)

maf_selinfo <- maf_selinfo %>% 
  filter(!(GeneName %in% bad_genes))

maf_selinfo <- filter(maf_selinfo, !POP == "ALL")

table(maf_selinfo$selection[maf_selinfo$POP=='EUR'])


maf_selinfo <- as.data.table(maf_selinfo)
genes_load <- eval_load(maf_selinfo, groups = c("selection", "GeneName", "POP"))[order(selection, pdbpn, n_snp)]
maf_selinfo <- as.data.frame(maf_selinfo)
# N. de genes
#length(unique(maf_selinfo$GeneName[maf_selinfo$selection=='positive']))
kable(genes_load[ 
  selection == 'positive',
  .(GeneName, n_snp, ps, pn, 
    pb, pd2, pd1)], # tirei pd1pneu, vê isso
  format = 'markdown', digits = 2)


#positivo <- filter(select( genes_load, selection, GeneName, ps, pn, pd1, pd2, pb), selection == 'positive')
#control <- filter(select( genes_load, selection, GeneName, ps, pn, pd1, pd2, pb), selection == 'control')

#EUR <- filter(maf_selinfo, POP == "EUR")
#genes_load <- eval_load(EUR, groups = c("selection", "GeneName"))[order(n_snp, decreasing = T)]
#kable(genes_load[ 
#   selection == 'positive',
#   .(GeneName, n_snp, Pdel, PdelPneu)], # tirei pd1pneu, vê isso
#   format = 'markdown', digits = 2)
#unique(EUR$GeneName[EUR$selection=="positive"])





































##########################################################################################
## Vizinhos 
##########################################################################################

aba2 <- fread("~/Dropbox/laboratorio/data/colonna/aba2.csv")
aba2 <- filter(aba2 , POPwithHighestDAF == 'EUR')
aba2 <- filter(aba2 , VT == 'SNP')
aba2 <- filter(aba2 , PAIR == 'AFR-EUR')
aba2 <- filter(aba2, GeneBiotype == "protein_coding")


maf_selinfo_vizinhos <- 
  maf_selinfo[!(GeneName %in% aba2$HGNCsymbol)]

genes_load <- eval_load(maf_selinfo_vizinhos[POP == "EUR"], groups = c("selection", "GeneName"))[order(selection, pdbpn, n_snp)]
kable(genes_load[ 
  selection == 'positive',
  .(GeneName, n_snp, pdbpn, pn, ps, 
    pb, pd1, pd2,   
    PNPS, pdpb, PdelPneu)][pdbpn > 0.95], # tirei pd1pneu, vê isso
  format = 'markdown', digits = 2)







##########################################################################################
## HighD 
##########################################################################################

vizinhos <-unique(maf_selinfo_vizinhos$GeneName[maf_selinfo_vizinhos$selection=='positive'])

maf_selinfo_HighD <- 
  maf_selinfo[!(GeneName %in% vizinhos)]


genes_load <- eval_load(maf_selinfo_HighD[POP == "EUR"], groups = c("selection", "GeneName"))[order(selection, pdbpn, n_snp)]
kable(genes_load[ 
  selection == 'positive',
  .(GeneName, n_snp, pdbpn, pn, ps, 
    pb, pd1, pd2,   
    pnps, pdpb, pd1pneu)][pdbpn > 0.95], # tirei pd1pneu, vê isso
  format = 'markdown', digits = 2)




# Atualzação de análises # 




##########################################################################################
##  SIFT, Grantham e Polyphen 
##########################################################################################


#SIFT
## Sobreposição SIFT e Polyphen

maf_selinfo <- mutate(maf_selinfo, Polyphen_Sift = paste(PolyPhenCat, SIFTcat , sep = "_"))
table(maf_selinfo$Polyphen_Sift[maf_selinfo$selection == "positive"])
table(maf_selinfo$Polyphen_Sift[maf_selinfo$selection == "control"])
table(maf_selinfo$Polyphen_Sift)

###### Testes de deleterius ######

maf_selinfo <- mutate(maf_selinfo, teste = Polyphen_Sift )
maf_selinfo[maf_selinfo$teste == "benign_deleterious", "teste"]<-"pb"
maf_selinfo[maf_selinfo$teste == " benign_NA ", "teste"]<- NA
maf_selinfo[maf_selinfo$teste == "NA_deleterious", "teste"]<- NA
maf_selinfo[maf_selinfo$teste == "possibly_damaging_deleterious", "teste"]<- "pd1"
maf_selinfo[maf_selinfo$teste == "possibly_damaging_NA", "teste"]<- NA
maf_selinfo[maf_selinfo$teste == "possibly_damaging_tolerated", "teste"]<- "pb"
maf_selinfo[maf_selinfo$teste == "probably_damaging_deleterious", "teste"]<- "pd1"
maf_selinfo[maf_selinfo$teste == "probably_damaging_NA", "teste"]<- NA
maf_selinfo[maf_selinfo$teste == "probably_damaging_tolerated", "teste"]<- "pb"
maf_selinfo[maf_selinfo$teste == "NA_NA", "teste"]<- NA
maf_selinfo[maf_selinfo$teste == "NA_tolerated", "teste"]<- NA
maf_selinfo[maf_selinfo$teste == "benign_NA", "teste"]<- NA
maf_selinfo[maf_selinfo$teste == "benign_tolerated", "teste"]<- "pb"




table(maf_selinfo$Polyphen_Sift)

benign_deleterious                     benign_NA              benign_tolerated 
12064                                     34                         58805 
NA_deleterious                         NA_NA (PS)                NA_tolerated 
101                                      76946                       245 
possibly_damaging_deleterious          possibly_damaging_NA   possibly_damaging_tolerated 
10296                                     14                         8489 
probably_damaging_deleterious          probably_damaging_NA   probably_damaging_tolerated 
20135                                     4                          6105



#Total dos pN em Polyphen_SIFT 
12064 + 34 + 58805 + 101  + 245 + 10296 + 14 + 8489 + 20135 + 4 + 6105
116292
# Calculo da sobreposicao total
benign_tolerated
58805
possibly_damaging_deleterious
10296
probably_damaging_deleterious
20135

(58805 + 10296 + 20135)/116292
0.7673443



# pb em polyphen e del em sift - no genoma
12064/116292
0.1037389
# pb em polyphen e neutro em sift - no genoma
58805/116292
0.5056668
# pd em polyphen e del em sift - no genoma
(10296 + 20135)/116292
0.2616775
# pd em polyphen e neutro em sift - no genoma
(8489+6105)/116292
0.1254944

#SIFT
## Sobreposição SIFT e Polyphen em  POSITIVOS

table(maf_selinfo$Polyphen_Sift[maf_selinfo$selection=='positive'])

benign_deleterious                     benign_NA 
48                                     0 
benign_tolerated                       NA_deleterious 
202                                    1 
NA_NA                                  NA_tolerated 
372                                    1 
possibly_damaging_deleterious          possibly_damaging_NA 
37                                      0 
possibly_damaging_tolerated            probably_damaging_deleterious 
43                                     113 
probably_damaging_tolerated 
29 


#Total dos pN em Polyphen_SIFT  GRUPO POSITIVO 
48 + 202 + 2 + 37  + 43 + 29 + 113
474
# Calculo da sobreposicao total
benign_tolerated
202
possibly_damaging_deleterious
37
probably_damaging_deleterious
113

(202 + 37 + 113)/474
0.742616



# pb em polyphen e del em sift  no positivo
48/474
0.1012658
# pb em polyphen e neutro em sift - no positivo
202/474
0.4261603
# pd em polyphen e del em sift - no positivo
(113 + 37)/474
0.3164557
# pd em polyphen e neutro em sift - no positivo
(43+29)/474
0.1518987




