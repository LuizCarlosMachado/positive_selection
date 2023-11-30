#Algumas estatísticas por bins, diversidade nao muda nada entre o controle e a janela positiva. 

## Packages
library()
library(data.table)
library(dplyr)
library(ggplot2)
library(ggthemes)

##Function
set_snv_selection <- 
  function(data, selected_genes, exons = FALSE) {
    if(exons == FALSE){
      z <-   
        unique(fread('eur.pos.cadd.tsv'))
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
      Heterozygosity = mean(Pi, na.rm = TRUE),
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
      pdelpneu = (pd1 + pd2) / (pb + ps),
      pd1ps = pd1 / ps,
      pdpn = (pd1 + pd2) / pn,
      pdbpn = (pd1 + pd2 + pb) / pn, 
      pbps = pb / ps,
      w_pdelpneu = (wd1 + wd2) / (wb + ws ),
      w_pd1ps = wd1 / ws,
      w_pnps = wn / ws,
      w_pdpb = (wd1 + wd2) / wb,
      w_pdps = (wd1 + wd2) / ws) ]
  }

load_per_bin<-
  function(data, 
           figure_name = 'european_bin_load.png',
           n_bins  = 10,
           pop = "EUR", 
           measures_to_plot = c(
             "pbps", 
             "pdpn", 
             "pdbpn",
             "pnps", 
             "pdelpneu",
             "pdpb",
             "PHRED",
             "Heterozygosity"),
           ncols = 2,
           leg_position = 'right'){
    
    data[, bins := cut(logMAF, n_bins, labels=FALSE), keyby = .(POP)][
      , maxlogMAF := max(logMAF), by = .(POP, bins)]
    
    bin_load <-
      eval_load(data, groups = c("POP", "selection", "bins", "maxlogMAF"))
    
    set_data_plot <- 
      function(data){
        melt(data, measure.vars = measures_to_plot, 
             variable.name = "measure", 
             value.name = "value", 
             varialble.factor = TRUE)
      }
    
    load_bin_plot <-
      set_data_plot(bin_load)
    
    ggplot(data=transform(load_bin_plot[POP == pop], measure = factor(measure, levels = measures_to_plot)), aes(x= maxlogMAF, y=value, color = selection)) + 
      geom_point(size=2) +
      geom_smooth(se = FALSE, method=glm) + 
      facet_wrap(~measure, ncol=ncols, scales="free_y") +
      labs(title = "Paired Statistics per Bins", x = "MaxlogMAF", y = "Value")  +
      theme_hc() + scale_colour_hc() + scale_fill_brewer(palette = "Paired", direction = -1)  
    ggsave(paste0('./figures/', figure_name))
    
    data[, bins := NULL]
  }

## Full Data 

maf <-fread('/Users/luiz/Dropbox/2016/mestrado/data/eur.maf.tsv')
n=30000 
eur.pos.cadd <- filter(maf, 
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
                         (CHR==22 & POS>46771921-n&POS<46771921+n))

eur.pos.cadd <- select(eur.pos.cadd, c(CHR, POP, GeneName))
eur.pos.cadd <- filter(eur.pos.cadd, POP=="EUR")
write.table(eur.pos.cadd, "eur.pos.cadd.tsv", row.names=F)

## Contruindo os Pi 
# Pi = 2*(MAF)*(1-MAF) - "2pq"
# Formando a coluna coluna Pi
maf <- mutate(maf,
              Pi = 2*(MAF)*(1-MAF))
# Formando a coluna 'selection' com a janela positiva
maf <- set_snv_selection(maf)
load_per_bin(maf)


