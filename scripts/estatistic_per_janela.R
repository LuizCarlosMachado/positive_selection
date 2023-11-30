##Script inteiro para populacoes Europeias


#Packages
library(data.table)
library(splitstackshape)
library(ggplot2)
library(GGally)
library(knitr)
library(dplyr)
library(ggthemes)
library(dtplyr)

setwd("~/Dropbox/laboratorio/scripts/european_pop/")
#maf <-fread('~/Dropbox/laboratorio/data/all.pop.maf.tsv')


###Funcoes necessárias

opa <- function (base_size = 12, base_family = "") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(legend.background = element_blank(), legend.key = element_blank(), 
          panel.background = element_blank(), legend.position = "right", panel.border = element_blank(), 
          strip.background = element_rect(fill = "gray92", colour = NA), plot.background = element_blank(), 
          axis.ticks = element_line(), axis.ticks.x = element_blank(), 
          axis.ticks.y = element_blank(), axis.ticks.length = unit(1, 
                                                                   "lines"))
  
}

set_snv_selection <- 
  function(data, n = 50000, selected_genes, exons = FALSE) {
    if(exons == FALSE){
      z <-as.data.table(unique(filter(select(filter(maf, 
                                                      (CHR==1 & POS>116935068-n & POS<116935068+n) | (CHR==1 & POS>246074181-n &POS<246074181+n)| 
                                                      (CHR==2 & POS>215975232-n&POS<215975232+n)|(CHR==2 & POS>72368190-n &POS<72368190+n)|(CHR==2 & POS>29980408-n &POS<29980408+n)| (CHR==2 & POS>158126458-n &POS<158126458+n)| (CHR==2 & POS>236886196-n &POS<236886196+n)|
                                                      (CHR==3 & POS>188703960-n&POS<188703960+n)|(CHR==3 & POS>188675732-n &POS<188675732+n)|(CHR==3 & POS>71088268-n &POS<71088268+n)| (CHR==3 & POS>96610897-n &POS<96610897+n)|
                                                      (CHR==4 & POS>38803255-n&POS<38803255+n)|(CHR==4 & POS>148820303-n &POS<148820303+n)|(CHR==4 & POS>75036044-n &POS<75036044+n)|
                                                      (CHR==6 & POS>136516257-n&POS<136516257+n)|
                                                      (CHR==7 & POS>75238617-n&POS<75238617+n)|(CHR==7 & POS>157197114-n &POS<157197114+n)|(CHR==7 & POS>146401434-n &POS<146401434+n)|                
                                                      (CHR==8 & POS>54704553-n&POS<54704553+n)|
                                                      (CHR==10 & POS>31786137-n&POS<31786137+n)|
                                                      (CHR==11 & POS>19620227-n & POS<19620227+n) | (CHR==11 & POS>64532579-n &POS<64532579+n)| 
                                                      (CHR==13 & POS>72272535-n & POS<72272535+n) | (CHR==13 & POS>49067103-n &POS<49067103+n)| 
                                                      (CHR==14 & POS>62038307-n&POS<62038307+n)|
                                                      (CHR==15 & POS>34258834-n&POS<34258834+n)|
                                                      (CHR==16 & POS>76523792-n & POS<76523792+n) | (CHR==16 & POS>84007565-n &POS<84007565+n)| 
                                                      (CHR==17 & POS>73349497-n&POS<73349497+n)|
                                                      (CHR==18 & POS>67624554-n&POS<67624554+n)|
                                                      (CHR==19 & POS>51062940-n&POS<51062940+n)|
                                                      (CHR==22 & POS>46771921-n&POS<46771921+n)) , c(CHR, POP, GeneName)), POP=="EUR")))
      
      
      z <- z[!is.na(GeneName),]
      z[,selection:="positive"]
      
      df <- merge(data, 
                  z,
                  all.x = TRUE, 
                  by = c("CHR", "GeneName", "POP"))
      
      df[is.na(selection), selection := 'control']
    }
  }

eval_load <- 
  function(data, groups = c("")){
    data[,.(
      pn = as.double(sum(Consequence == 'NON_SYNONYMOUS', na.rm = TRUE)), 
      ps = as.double(sum(Consequence == 'SYNONYMOUS', na.rm = TRUE)), 
      pd1 = as.double(sum(PolyPhenCat == "probably_damaging", na.rm = TRUE)), 
      pd2 = as.double(sum(PolyPhenCat  == "possibly_damaging", na.rm = TRUE)), 
      pb = as.double(sum(PolyPhenCat == "benign", na.rm = TRUE)),  
      sift_del = as.double(sum(SIFTcat == "deleterious", na.rm = TRUE)),
      sift_neu = as.double(sum(SIFTcat == "tolerated", na.rm = TRUE)),
      Grantham = mean(Grantham, na.rm = TRUE), 
      htz = mean(htz, na.rm = TRUE), 
      PHRED = mean(PHRED, na.rm = TRUE), 
      RawScore = mean(RawScore, na.rm = TRUE), 
      n_PHRED = sum(!is.na(PHRED)), 
      wn = sum((Consequence == 'NON_SYNONYMOUS') * MAF, na.rm = TRUE), 
      ws = sum((Consequence == 'SYNONYMOUS') * MAF, na.rm = TRUE), 
      wd1 = sum((PolyPhenCat == "probably_damaging") * MAF, na.rm = TRUE), 
      wd2 = sum((PolyPhenCat  == "possibly_damaging") * MAF, na.rm = TRUE), 
      wb = sum((PolyPhenCat == "benign") * MAF, na.rm = TRUE),  
      w_PHRED = mean(PHRED * MAF, na.rm = TRUE),
      w_RawScore = mean(RawScore * MAF, na.rm = TRUE),
      mMAF = mean(MAF, na.rm = TRUE),
      mhtz = 2 * mean(MAF * (1 - MAF), na.rm = TRUE)
    ), by = groups][,`:=`(
      n_snp = pn + ps,  
      n_db = pd1 + pd2 + pb,
      wd = wd1 + wd2, 
      pdpb = (pd1 + pd2) / pb,
      pdps = (pd1 + pd2) / ps,
      # Jonatas Pdelpneu = pd1 / (pb + ps + pd2),
      pd1ps = pd1 / ps,
      pdpn = (pd1 + pd2) / pn,
      pdbps = (pd1 + pd2 + pb) / ps, 
      
      ## LOAD STATS
      PdelPneu = (pd1)/(pb + ps + pd2),
      PdelPS = (pd1)/ps,
      Pd1 = pd1,
      Pdel = pd1,
      Pdel_Poly = pd1 + pd2,
      Pd = (pd1 + pd2),
      Pdel_Sift = sift_del,
      PNPS = pn / ps,
      PdelPneu_SIFT = sift_del/(sift_neu+ps),
      PNPT = pn/(pn+ps),
      PdelPT = (pd1)/(pb + ps + pd2 + pd1),
      #J_Pdelpneu = pd1 / (pb + ps + pd2),
      
      pdbpn = (pd1 + pd2 + pb) / pn, 
      sift_pn = (sift_del + sift_neu)/ pn,
      SNPs_Neu = pb + ps,
      pd1pb = pd1/pb,
      pn_nsnp = pn/(pn+ps),
      pNpS_Sift = (sift_del + sift_neu)/ps,
      pdelPneu_SIFT = sift_del/(sift_neu+ps),
      #pdelPneu_SIFT = sift_del/ps,
      pbps = pb / ps,
      pbpn = pb / pn,
      pd1pn = pd1 / pn,
      pd2pn = pd2 / pn,
      w_pd1pneu = wd1 / (wb + ws + wd2),
      w_pd1ps = wd1 / ws,
      w_pnps = wn / ws,
      w_pdpb = (wd1 + wd2) / wb,
      pd1_n_snp = pd1/(pn+ps),
      w_pdps = (wd1 + wd2) / ws)]
  }
stat_per_win<- function(data,valores_n) {
  #criando objetos intermediarios
  lista_result<-list(NA)
  for(i in 1:length(valores_n)){
    x <- set_snv_selection(data, valores_n[i], selected_genes = pop_info)
    x <- as.data.table(mutate(x, htz = 2*(MAF)*(1-MAF)))
    y <-select(eval_load(x, c("selection")),  c(htz, PdelPneu, PNPS, PdelPT, PNPT, pdelPneu_SIFT)) # ou qual mais quiser do eval_loa
    y$selection[1]<-"control"
    y$selection[2]<-"positive"
    lista_result[[i]]<- y
  }
  return(do.call("rbind",lista_result))
}

#Tem que gerar o Bad Genes
##Código
maf2 <- filter(maf, POP == "EUR")
maf2 <- maf2[!(GeneName %in% bad_genes)]
lista <-stat_per_win(maf2,(c( seq(10,1000, 5)*1000)))
lista <- mutate(lista, janela = rep(c(seq(10,1500, 5))*2, each= 2))
lista <- rename(lista, heterozigosidade=htz)

write.table(lista, "~/Dropbox/laboratorio/data/estatisticas_por_janela.POS_EUR_POP_EUR2.tsv", row.names=F)
#write.table(lista, "~/Dropbox/laboratorio/data/estatisticas_por_janela.POS_EUR_POP_EAS2.tsv", row.names=F)
#write.table(lista, "~/Dropbox/laboratorio/data/estatisticas_por_janela.POS_EUR_POP_AFR2.tsv", row.names=F)
#write.table(lista, "~/Dropbox/laboratorio/data/estatisticas_por_janela.POS_EUR_POP_SAS2.tsv", row.names=F)

# POS ASIATICAS NO MUNDO 

eur <- fread("~/Dropbox/laboratorio/data/estatisticas_por_janela.POS_EAS_POP_EUR.tsv")
eur <- melt(eur, measure=1:6, value.factor=TRUE)
eur <- mutate(eur, POP = "EUR")
eur <- filter(eur, !janela < 70)

eas <- fread("~/Dropbox/laboratorio/data/estatisticas_por_janela.POS_EAS_POP_EAS.tsv")
eas <- melt(eas, measure=1:6, value.factor=TRUE)
eas <- mutate(eas, POP = "EAS")
eas <- filter(eas, !janela < 70)

afr <- fread("~/Dropbox/laboratorio/data/estatisticas_por_janela.POS_EAS_POP_AFR.tsv")
afr <- melt(afr, measure=1:6, value.factor=TRUE)
afr <- mutate(afr, POP = "AFR")
afr <- filter(afr, !janela < 70)

sas <- fread("~/Dropbox/laboratorio/data/estatisticas_por_janela.POS_EAS_POP_SAS.tsv")
sas <- melt(sas, measure=1:6, value.factor=TRUE)
sas <- mutate(sas, POP = "SAS")
sas <- filter(sas, !janela < 70)

data1 <- rbind(eur,eas,afr,sas)
data1 <- mutate(data1, under_sel = "EAS")
data1$selection<- factor(data1$selection,levels = c("positive","control"))
# As posicoes sao as mesmas do EAS em diferentes regioes do mundo, e so EAS esta sob selecao

################################
#### POS EUROPEIAS NO MUNDO ####
################################


eur <- fread("~/Dropbox/laboratorio/data/estatisticas_por_janela.POS_EUR_POP_EUR.tsv")
eur <- melt(eur, measure=1:6, value.factor=TRUE)
eur <- mutate(eur, POP = "EUR")

eas <- fread("~/Dropbox/laboratorio/data/estatisticas_por_janela.POS_EUR_POP_EAS.tsv")
eas <- melt(eas, measure=1:6, value.factor=TRUE)
eas <- mutate(eas, POP = "EAS")

afr <- fread("~/Dropbox/laboratorio/data/estatisticas_por_janela.POS_EUR_POP_AFR.tsv")
afr <- melt(afr, measure=1:6, value.factor=TRUE)
afr <- mutate(afr, POP = "AFR")

sas <- fread("~/Dropbox/laboratorio/data/estatisticas_por_janela.POS_EUR_POP_SAS.tsv")
sas <- melt(sas, measure=1:6, value.factor=TRUE)
sas <- mutate(sas, POP = "SAS")


data2 <- rbind(eur,eas,afr,sas)
data2 <- mutate(data2, under_sel = "EUR")
data2$selection<- factor(data2$selection,levels = c("positive","control"))
# As posicoes sao as mesmas do EUR em diferentes regioes do mundo, e so EAS esta sob selecao

#write.table(data, "~/Dropbox/laboratorio/data/estatisticas_por_janela.completa.tsv", row.names=F)
#data<- fread("~/Dropbox/laboratorio/data/estatisticas_por_janela.completa.tsv")


#Gráficos e testes
# PdelPneu
  data <- rbind(data1, data2)
  teste <- filter(data, variable == "PdelPneu" )
  teste$selection<- factor(teste$selection,levels = c("positive","control"))
  m <- ggplot(teste, aes(x=janela, y=value, fill = POP))
  m + geom_line(aes(colour = POP, fill = selection, linetype = selection), size = 1) + facet_grid(~  variable ~ under_sel,  scales = "free") + opa() + xlim(0,1000) +  labs( x = "Distance of HighD SNP (kb) ", y = "") +  scale_colour_manual(values = c("#CC0000", "#009966", "#006699", "#FF9933")) +
    theme(axis.text = element_text(size = rel(1.2)), legend.position = "right", plot.title = element_text(size = rel(1.5)) ,axis.title.y = element_text(size = rel(1.5)), axis.title.x = element_text(size = rel(1.5)),legend.title = element_text(size=16, face="bold")
          , legend.text=element_text(size=16), legend.key.size = unit(1.2, "cm"))  + geom_vline(xintercept = 240)

 #ggsave("./figures/load_win.png")



  # Heterozigosidade
data3 <- filter(data1, POP == "EAS")
data4 <- filter(data2, POP == "EUR")
data <- rbind(data3, data4)
teste <- filter(data, variable == "heterozigosidade" )
teste$selection<- factor(teste$selection,levels = c("positive","control"))

m <- ggplot(teste, aes(x=janela, y=value, fill = POP))
m + geom_line(aes(colour = POP, fill = selection, linetype = selection), size = 1.5) + facet_grid(~  variable ~ under_sel,  scales = "free") + opa() + xlim(0,2250) +  labs(x = "Distance of HighD SNP (kb) ", y = "") +  
  scale_colour_manual( values = c("#009966", "#006699")) +
  theme(axis.text = element_text(size = rel(1.2)), legend.position = "right", plot.title = element_text(size = rel(1.5)) ,axis.title.y = element_text(size = rel(1.5)), axis.title.x = element_text(size = rel(2.0)),legend.title = element_text(size=16, face="bold")
        , legend.text=element_text(size=16), legend.key.size = unit(1.2, "cm")) + geom_vline(xintercept = 200)

#ggsave("./figures/heterozigosidade.png")



# Diferenca entre valor e media genomica no PdelPneu 

data<-as.data.frame(data)
positivo<-data[as.character(data$selection)=="positive",]
controle<-data[as.character(data$selection)=="control",]

positivo$dif<-positivo$value - controle$value
head(positivo)
positivo <- as.data.table(positivo)

teste <- filter(positivo, variable == "PdelPneu" )
m <- ggplot(teste, aes(x=janela, y=dif, fill = POP))
m + geom_line(aes(colour = POP), size = 1.5) + facet_grid(~  variable ~ under_sel,  scales = "free") + opa() + xlim(0,1000) +  
  labs( x = "Distance of HighD SNP (kb) ", y = "") +  scale_colour_manual(values = c("#CC0000", "#009966", "#006699", "#FF9933")) +
  theme(axis.text = element_text(size = rel(1.2)), legend.position = "right", plot.title = element_text(size = rel(1.5)) ,axis.title.y = element_text(size = rel(1.5)), axis.title.x = element_text(size = rel(1.5)),legend.title = element_text(size=16, face="bold")
        , legend.text=element_text(size=16), legend.key.size = unit(1.2, "cm")) + geom_hline(linetype = 2, yintercept = 0)


+ geom_vline(xintercept = c(125, 1500))
#gsave("./figures/load_win_mean_diff.png")








library(ggplot2)
library(data.table)
library(dplyr)
z_score <- fread("~/Dropbox/laboratorio/data/ZScore3000EUR.tsv")

z_score <- fread("~/Dropbox/laboratorio/data/ZScore2000EAS.tsv")
z_score <- fread("~/Dropbox/laboratorio/data/ZScore2000EUR.tsv")
z_score <- mutate(z_score, janela = seq(10,2000,10))
z_score <- mutate(z_score, estatistica = "zScore")
z_score <- melt(z_score,  measure=1:4, value.factor=TRUE)
z_score <- rename(z_score, POP = variable)

m <- ggplot(z_score, aes(x=janela, y=value, fill = POP))
m + geom_line(aes(colour = POP, fill = POP), size = 0.5) + opa()  +  xlim(0, 800) + labs( x = "Distance of HighD SNP (kb) ", y = "") +  scale_colour_manual(values = c("#CC0000", "#009966", "#006699", "#FF9933")) + geom_vline(xintercept = c(120, 510, 490, 530), color = c("#006699", "#CC0000", "#009966", "#FF9933" )) +
  theme(axis.text = element_text(size = rel(1.2)), legend.position = "right", plot.title = element_text(size = rel(1.5)) ,axis.title.y = element_text(size = rel(1.5)), axis.title.x = element_text(size = rel(1.5)),legend.title = element_text(size=16, face="bold")
        , legend.text=element_text(size=16), legend.key.size = unit(1.2, "cm"))  



###########################################



#positivo <- fread("~/Dropbox/laboratorio/data/dif_estatisticas_por_janela.completa.tsv")

#Gráficos e testes
# PdelPneu
teste <- filter(positivo, variable == "PdelPneu")
m <- ggplot(teste, aes(x=janela, y=dif, fill = POP))
m + geom_line(aes(colour = POP), size = 1.5) + facet_grid(~  variable ~ under_sel ~ .,  scales = "free") + opa() + xlim(20,1000) +  
  labs( x = "", y = "") +  scale_colour_manual(values = c("#CC0000", "#009966", "#006699", "#FF9933")) +
  theme(axis.text = element_text(size = rel(1.2)), legend.position = "none", plot.title = element_text(size = rel(1.5)) ,axis.title.y = element_text(size = rel(1.5)), axis.title.x = element_text(size = rel(1.5)),legend.title = element_text(size=16, face="bold")
        , legend.text=element_text(size=16), legend.key.size = unit(1.2, "cm")) + geom_hline(linetype = 2, yintercept = 0)


+ geom_vline(xintercept = c(125, 1500))
#gsave("./figures/load_win_mean_diff.png")