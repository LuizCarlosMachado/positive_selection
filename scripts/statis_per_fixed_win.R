##Script inteiro para populacoes Europeias


setwd("~/Dropbox/laboratorio/scripts/european_pop/")



library(data.table)
library(splitstackshape)
library(ggplot2)
library(GGally)
library(knitr)
library(dplyr)
library(dtplyr)
library(ggthemes)

set_snv_selection <- 
  function(data, i = 2, n = 39000, selected_genes, exons = FALSE) {
    if(exons == FALSE){
      z <-as.data.table(unique(filter(select(filter(data,  (CHR==1 &  POS > (116935068-i*n)+n &
                                     POS < (116935068-(i-1)*n)+n) |
                      (CHR==1 &
                         POS > (116935068+(i-1)*n)-n &
                         POS < (116935068+i*n)-n)
                    |
                      
                      (CHR==1 &  POS > (246074181-i*n)+n &
                         POS < (246074181-(i-1)*n)+n) |
                      (CHR==1 &
                         POS > (246074181+(i-1)*n)-n &
                         POS < (246074181+i*n)-n)
                    |
                      
                      (CHR==2 &  POS > (215975232-i*n)+n &
                         POS < (215975232-(i-1)*n)+n) |
                      (CHR==2 &
                         POS > (215975232+(i-1)*n)-n &
                         POS < (215975232+i*n)-n)
                    |
                      
                      (CHR==2 &  POS > (72368190-i*n)+n &
                         POS < (72368190-(i-1)*n)+n) |
                      (CHR==2 &
                         POS > (72368190+(i-1)*n)-n &
                         POS < (72368190+i*n)-n)
                    
                    |
                      
                      (CHR==2 &  POS > (29980408-i*n)+n &
                         POS < (29980408-(i-1)*n)+n) |
                      (CHR==2 &
                         POS > (29980408+(i-1)*n)-n &
                         POS < (29980408+i*n)-n)
                    |
                      
                      (CHR==2 &  POS > (158126458-i*n)+n &
                         POS < (158126458-(i-1)*n)+n) |
                      (CHR==2 &
                         POS > (158126458+(i-1)*n)-n &
                         POS < (158126458+i*n)-n)
                    |
                      
                      (CHR==2 &  POS > (116935068-i*n)+n &
                         POS < (116935068-(i-1)*n)+n) |
                      (CHR==2 &
                         POS > (116935068+(i-1)*n)-n &
                         POS < (116935068+i*n)-n)
                    |
                      
                      (CHR==2 &  POS > (236886196-i*n)+n &
                         POS < (236886196-(i-1)*n)+n) |
                      (CHR==2 &
                         POS > (236886196+(i-1)*n)-n &
                         POS < (236886196+i*n)-n)
                    
                    |
                      
                      (CHR==3 &  POS > (188703960-i*n)+n &
                         POS < (188703960-(i-1)*n)+n) |
                      (CHR==3 &
                         POS > (188703960+(i-1)*n)-n &
                         POS < (188703960+i*n)-n)
                    |
                      
                      (CHR==3 &  POS > (188675732-i*n)+n &
                         POS < (188675732-(i-1)*n)+n) |
                      (CHR==3 &
                         POS > (188675732+(i-1)*n)-n &
                         POS < (188675732+i*n)-n)
                    |
                      
                      (CHR==3 &  POS > (71088268-i*n)+n &
                         POS < (71088268-(i-1)*n)+n) |
                      (CHR==3 &
                         POS > (71088268+(i-1)*n)-n &
                         POS < (71088268+i*n)-n)
                    |
                      
                      (CHR==3 &  POS > (96610897-i*n)+n &
                         POS < (96610897-(i-1)*n)+n) |
                      (CHR==3 &
                         POS > (96610897+(i-1)*n)-n &
                         POS < (96610897+i*n)-n)
                    
                    |
                      
                      (CHR==4 &  POS > (38803255-i*n)+n &
                         POS < (38803255-(i-1)*n)+n) |
                      (CHR==4 &
                         POS > (38803255+(i-1)*n)-n &
                         POS < (38803255+i*n)-n)
                    
                    |
                      
                      (CHR==4 &  POS > (148820303-i*n)+n &
                         POS < (148820303-(i-1)*n)+n) |
                      (CHR==4 &
                         POS > (148820303+(i-1)*n)-n &
                         POS < (148820303+i*n)-n)
                    
                    |
                      
                      (CHR==4 &  POS > (75036044-i*n)+n &
                         POS < (75036044-(i-1)*n)+n) |
                      (CHR==4 &
                         POS > (75036044+(i-1)*n)-n &
                         POS < (75036044+i*n)-n)
                    
                    |
                      
                      (CHR==6 &  POS > (136516257-i*n)+n &
                         POS < (136516257-(i-1)*n)+n) |
                      (CHR==6 &
                         POS > (136516257+(i-1)*n)-n &
                         POS < (136516257+i*n)-n)
                    
                    |
                      
                      (CHR==7 &  POS > (75238617-i*n)+n &
                         POS < (75238617-(i-1)*n)+n) |
                      (CHR==7 &
                         POS > (75238617+(i-1)*n)-n &
                         POS < (75238617+i*n)-n)
                    |
                      
                      (CHR==7 &  POS > (157197114-i*n)+n &
                         POS < (157197114-(i-1)*n)+n) |
                      (CHR==7 &
                         POS > (157197114+(i-1)*n)-n &
                         POS < (157197114+i*n)-n)
                    |
                      
                      (CHR==7 &  POS > (146401434-i*n)+n &
                         POS < (146401434-(i-1)*n)+n) |
                      (CHR==7 &
                         POS > (146401434+(i-1)*n)-n &
                         POS < (146401434+i*n)-n)
                    |
                      
                      (CHR==8 &  POS > (54704553-i*n)+n &
                         POS < (54704553-(i-1)*n)+n) |
                      (CHR==8 &
                         POS > (54704553+(i-1)*n)-n &
                         POS < (54704553+i*n)-n)
                    |
                      
                      (CHR==10 &  POS > (31786137-i*n)+n &
                         POS < (31786137-(i-1)*n)+n) |
                      (CHR==10 &
                         POS > (31786137+(i-1)*n)-n &
                         POS < (31786137+i*n)-n)
                    
                    |
                      
                      (CHR==11 &  POS > (64532579-i*n)+n &
                         POS < (64532579-(i-1)*n)+n) |
                      (CHR==11 &
                         POS > (64532579+(i-1)*n)-n &
                         POS < (64532579+i*n)-n)
                    
                    |
                      
                      (CHR==11 &  POS > (19620227-i*n)+n &
                         POS < (19620227-(i-1)*n)+n) |
                      (CHR==11 &
                         POS > (19620227+(i-1)*n)-n &
                         POS < (19620227+i*n)-n)
                    |
                      
                      (CHR==13 &  POS > (72272535-i*n)+n &
                         POS < (72272535-(i-1)*n)+n) |
                      (CHR==13 &
                         POS > (72272535+(i-1)*n)-n &
                         POS < (72272535+i*n)-n)
                    
                    |
                      
                      (CHR==13 &  POS > (49067103-i*n)+n &
                         POS < (49067103-(i-1)*n)+n) |
                      (CHR==13 &
                         POS > (49067103+(i-1)*n)-n &
                         POS < (49067103+i*n)-n)
                    |
                      
                      (CHR==14 &  POS > (62038307-i*n)+n &
                         POS < (62038307-(i-1)*n)+n) |
                      (CHR==14 &
                         POS > (62038307+(i-1)*n)-n &
                         POS < (62038307+i*n)-n)
                    |
                      
                      (CHR==15 &  POS > (34258834-i*n)+n &
                         POS < (34258834-(i-1)*n)+n) |
                      (CHR==15 &
                         POS > (34258834+(i-1)*n)-n &
                         POS < (34258834+i*n)-n)
                    |
                      
                      (CHR==16 &  POS > (76523792-i*n)+n &
                         POS < (76523792-(i-1)*n)+n) |
                      (CHR==16 &
                         POS > (76523792+(i-1)*n)-n &
                         POS < (76523792+i*n)-n)
                    |
                      
                      (CHR==16 &  POS > (84007565-i*n)+n &
                         POS < (84007565-(i-1)*n)+n) |
                      (CHR==16 &
                         POS > (84007565+(i-1)*n)-n &
                         POS < (84007565+i*n)-n)
                    |
                      
                      (CHR==17 &  POS > (73349497-i*n)+n &
                         POS < (73349497-(i-1)*n)+n) |
                      (CHR==17 &
                         POS > (73349497+(i-1)*n)-n &
                         POS < (73349497+i*n)-n)
                    |
                      
                      (CHR==18 &  POS > (67624554-i*n)+n &
                         POS < (67624554-(i-1)*n)+n) |
                      (CHR==18 &
                         POS > (67624554+(i-1)*n)-n &
                         POS < (67624554+i*n)-n)
                    |
                      
                      (CHR==19 &  POS > (51062940-i*n)+n &
                         POS < (51062940-(i-1)*n)+n) |
                      (CHR==19 &
                         POS > (51062940+(i-1)*n)-n &
                         POS < (51062940+i*n)-n)
                    
                    |
                      
                      (CHR==22 &  POS > (46771921-i*n)+n &
                         POS < (46771921-(i-1)*n)+n) |
                      (CHR==22 &
                         POS > (46771921+(i-1)*n)-n &
                         POS < (46771921+i*n)-n))
        , c(CHR, POP, GeneName)), POP=="EUR")))
      
      
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
      PHRED = mean(PHRED, na.rm = TRUE), 
      htz = mean(htz, na.rm = TRUE), 
      RawScore = mean(RawScore, na.rm = TRUE), 
      n_PHRED = sum(!is.na(PHRED)), 
      wn = sum((Consequence == 'NON_SYNONYMOUS') * MAF, na.rm = TRUE), 
      ws = sum((Consequence == 'SYNONYMOUS') * MAF, na.rm = TRUE), 
      wd1 = sum((PolyPhenCat == "probably_damaging") * MAF, na.rm = TRUE), 
      wd2 = sum((PolyPhenCat  == "possibly_damaging") * MAF, na.rm = TRUE), 
      wb = sum((PolyPhenCat == "benign") * MAF, na.rm = TRUE),  
      w_PHRED = mean(PHRED * MAF, na.rm = TRUE),
      w_RawScore = mean(RawScore * MAF, na.rm = TRUE)
    ), by = groups][,`:=`(
      n_snp = pn + ps,  
      n_db = pd1 + pd2 + pb,
      pd = pd1 + pd2, 
      wd = wd1 + wd2, 
      pnps = pn / ps,
      pdpb = (pd1 + pd2) / pb,
      pdps = (pd1 + pd2) / ps,
      pd1pneu = pd1 / (pb + ps + pd2),
      pd1ps = pd1 / ps,
      pdpn = (pd1 + pd2) / pn,
      pdbpn = (pd1 + pd2 + pb) / pn, 
      pbps = pb / ps,
      w_pd1pneu = wd1 / (wb + ws + wd2),
      w_pd1ps = wd1 / ws,
      w_pnps = wn / ws,
      w_pdpb = (wd1 + wd2) / wb,
      w_pdps = (wd1 + wd2) / ws)]
  }


stat_per_win<- function(data,valores_n) {
  #criando objetos intermediarios
  lista_result<-list(NA)
  for(i in 1:length(valores_n)){
    x <- set_snv_selection(data, valores_n[i], n = 39000, selected_genes = pop_info)
    x <- as.data.table(mutate(x, htz = 2*(MAF)*(1-MAF)))
    y <-select(eval_load(x, c("selection")), c(htz, pd1pneu, pnps, PHRED, n_snp, pn, ps, pd1, pd2, pb)) # ou qual mais quiser do eval_loa
    y$treatment[1]<-"control"
    y$treatment[2]<-"positive"
    lista_result[[i]]<- y
  }
  return(do.call("rbind",lista_result))
}


###all data### 

sample(100000, 5, replace = T)
maf <-fread('~/Dropbox/laboratorio/data/eur.maf.tsv')

lista <-stat_per_win(maf,( sample(100000, 5, replace = T)))
lista <- as.data.table(filter(lista, treatment == "positive")) 
lista <- as.data.table(mutate(lista, n_da_janela_39kb = seq(1, 4996, 1)))
lista <- rename(lista, heterozigosidade=htz)
lista <- rename(lista, pdelpneu = pd1pneu )
#write.table(lista, "~/Dropbox/laboratorio/data/estatisticas_por_janela_deslizante_25kb.tsv", row.names=F)
lista <- fread("~/Dropbox/laboratorio/data/estatisticas_por_janela_deslizante_39kb.tsv")
#lista <- mutate(lista, treatment = rep("control"))
#list <- edit(lista)
#lista <- as.data.table(list)


measures_to_plot <- c("n_snp")
q<- melt(lista, measure.vars = measures_to_plot, 
         variable.name = "measure", 
         value.name = "value", 
         varialble.factor = TRUE)

high_sd <- mean((q$value)[q$treatment=="control"]) + sd((q$value)[q$treatment=="control"])
my_mean <-mean((q$value)[q$treatment=="control"])
low_sd <- mean((q$value)[q$treatment=="control"]) - sd((q$value)[q$treatment=="control"])
positive_value <- mean((q$value)[q$treatment=="positive"])
x <- ggplot(q, aes(n_da_janela_39kb, value)) +
  geom_line(aes(colour = treatment), size = 0.5) + 
  facet_wrap(~measure, nrow = 4, scales="free_y") + geom_hline(yintercept = my_mean, colour="seagreen") + geom_hline(yintercept = high_sd,  linetype = 2, colour="seagreen")+ geom_hline(yintercept = low_sd, linetype = 2, colour="seagreen") + geom_hline(yintercept = mean((q$value)[q$treatment=="positive"]))  + geom_hline(yintercept = mean((q$value)[q$treatment=="positive"])) + geom_hline(yintercept = mean((q$value)[q$treatment=="positive"]), colour="orange") +
  labs(title = "Janelas de 39kb", x = "Distância em janelas de 39kb do HighD ", y = "Valor") +
  annotate("text", min(q$value), x = -15,  y = 1.05*high_sd , label = "sd", colour="seagreen") +
  annotate("text", min(q$value), x = -15,  y = 1.05*positive_value, label = "Positive", colour="orange3") +  
  annotate("text", min(q$value), x = -15,  y = 1.02*positive_value, label = "Value", colour="orange3") +
  annotate("text", min(q$value), x = -15,  y = 1.05*my_mean, label = "mean", colour="seagreen") +
  annotate("text", min(q$value), x = -15,  y = 1.05*low_sd, label = "sd", colour="seagreen") 
x + theme_minimal() + scale_colour_brewer(palette = "Set2")
#x + opa() + scale_colour_brewer(palette = "Set2")
#x + opa() +   scale_colour_manual(values = c("royalblue4", "steelblue3"))
ggsave("./figures/janelas_fixas_n_snp.png")




dist_mean <-mean((q$value)[q$treatment=="positive"])
dist_quant <- quantile(q$value, 0.05)
ggplot(q, aes(value, colour = measure)) +
geom_density() + geom_vline(xintercept =  mean((q$value)[q$treatment=="positive"]), linetype = 2) +  geom_vline(xintercept = quantile(q$value, 0.05)) +
theme_hc() +
scale_colour_brewer(palette = "Set2")   +
labs(title = "Distribuição de janelas de 39kb") +
annotate("text", min(q$value), x = 0.94*dist_mean,  y = 0.04, label = "positive line", size = 3) +
annotate("text", min(q$value), x = dist_quant+ 0.04*dist_quant,  y = 0.04, label = "0.05 line", size = 3) 
ggsave("./figures/dist_janelas_fixas_htz.png")




#numero de SNP

p <- ggplot(q, aes(n_da_janela_39kb, value)) + geom_point()
p


p <- ggplot(q, aes(n_da_janela_39kb, value, colour = treatment))
geom_line(aes(ymin = 100, ymax = 200))

control_mean <- mean(q$value[q$treatment=="control"])
value <- c((q$value[q$treatment=="positive"]), mean(q$value[q$treatment=="control"]))
sd(q$value[q$treatment=="control"])
teste <- data.frame(value, names )
p <- ggplot(teste, aes(names, value, colour = names))
p +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = 100, ymax = 200), position = "dodge", width = 0.25)
