library(ggplot2)
library(ggthemes)
aba2 <- fread("~/Dropbox/laboratorio/data/colonna/aba2.csv")
aba2 <- filter(aba2 , POPwithHighestDAF == 'EUR')
aba2 <- filter(aba2 , VT == 'SNP')
aba2 <- filter(aba2 , PAIR == 'AFR-EUR')
aba2 <- filter(aba2, GeneBiotype == "protein_coding")

table(maf_selinfo$selection)

genes.colonna <- select(aba2, HGNCsymbol)
vizinhos <- 
  select(filter(maf_selinfo[!(GeneName %in% genes.colonna$HGNCsymbol)], selection == "positive"), GeneName)


maf <-fread('~/Dropbox/laboratorio/data/eur.maf.tsv')





# Preparando janela
eur.win(maf, 100000)
genes_win <-select(eur.win(maf, 100000), GeneName)

##maf_selinfo
system.time(
  maf_selinfo <- set_snv_selection(maf, selected_genes = pop_info)
)
maf_selinfo<-as.data.table(mutate(maf_selinfo, htz = 2*(MAF)*(1-MAF)))

#rm(maf)
#genes_load, cada linha é um gene, e as colunas sao as informacoes desse gene (pn, ps,...)
genes_load <- eval_load(maf_selinfo[POP == "EUR"], groups = c("selection", "GeneName"))[order(selection, pdbpn, n_snp)]

#Bad_genes é pra ver quais genes nao estao descritos pelo PolyPhen, pdbpn = Snps_Pol / Pn, nesse script, todos sao bons 
bad_genes <- genes_load[ selection != 'control' & pdbpn < 0.98, GeneName]

maf_selinfo <- 
  maf_selinfo[!(GeneName %in% bad_genes)]


kable(genes_load[ 
  selection != 'control',
  .(selection, GeneName, n_snp, pdbpn, pn, ps, 
    pb, pd1, pd2,   
    pnps, pdpb, pd1pneu, PHRED)][pdbpn > 0.95], # tirei pd1pneu, vê isso
  format = 'markdown', digits = 2)


#########################################################################################################################################################################################################################################################################################################################################################################################################

maf_selinfo <- 
  maf_selinfo[!(GeneName %in% genes.colonna$HGNCsymbol)]
#maf_selinfo <- 
  #maf_selinfo[!(GeneName %in% vizinhos$GeneName)]

#rm(maf)
#genes_load, cada linha é um gene, e as colunas sao as informacoes desse gene (pn, ps,...)
genes_load <- eval_load(maf_selinfo[POP == "EUR"], groups = c("selection", "GeneName"))[order(selection, pdbpn, n_snp)]

#Bad_genes é pra ver quais genes nao estao descritos pelo PolyPhen, pdbpn = Snps_Pol / Pn, nesse script, todos sao bons 
bad_genes <- genes_load[ selection != 'control' & pdbpn < 0.98, GeneName]

maf_selinfo <- 
  maf_selinfo[!(GeneName %in% bad_genes)]

kable(genes_load[ 
  selection != 'control',
  .(selection, GeneName, n_snp, pdbpn, pn, ps, 
    pb, pd1, pd2,   
    pnps, pdpb, pd1pneu, PHRED)][pdbpn > 0.95 & ps > 0], 
  format = 'markdown', digits = 2)


#########################################################################################################################################################################################################################################################################################################################################################################################################


violin_plots <- 
  function(load_region, bootstrap_control,
           measures_to_plot = c(
                                "htz"
                                
                                
           ),
           pops = c("ALL", "AFR", "EUR", "EAS", "SAS"), 
           figure_name = 'violin.teste.png',
           th = 0.05,
           point_size = 2, 
           label_size = 3, 
           leg_position = 'none',
           ncols = 2){
    
    set_data_plot <- function(data){
      melt(data, measure.vars = measures_to_plot, 
           variable.name = "measure", 
           value.name = "value", 
           varialble.factor = TRUE)[
             ,.(POP, replica, measure, value)]
    }
    
    control_data_plot <-
      set_data_plot(
        eval_load(bootstrap_control, 
                  groups = c("POP", "replica"))
      )
    
    load_region_plot <-
      set_data_plot(load_region[, replica := NA])
    
    data_plot <- rbindlist(list(load_region_plot, control_data_plot))
    
    data_plot[, pvalue := (.N - frank(value) + 1)/ .N, 
              by = .(POP, measure)]
    
    data_plot[, `:=`(measure = factor(measure, levels = measures_to_plot),
                     POP = factor(POP, levels = pops))]
    
    data_plot[, plabelpos := max(value) + 0.2 * (max(value) - min(value)), by = .(measure)]
    data_plot[is.na(replica), is_significant := (pvalue <= th | pvalue >= 1 - th)]
    
    if(figure_name == 'none'){
      return(data_plot[, plabelpos := NULL])
    }else{
      ggplot(data = data_plot[!is.na(replica)], 
             aes(x=POP, y = value, fill=POP)) + 
        geom_violin() + 
        labs(x = "") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        geom_boxplot(width=0.1, outlier.size=0)  +
        geom_point(data = data_plot[is.na(replica)],
                   aes(x=POP, y=value), size= point_size, color='red') +
        geom_label(data = data_plot[is.na(replica)],
                   aes(x=POP, 
                       y=plabelpos, 
                       label = format(round(log10(pvalue),2), nsmall = 2), 
                       fill=is_significant), 
                   size=label_size) +
        facet_wrap(~measure, ncol=ncols, scales="free_y", switch = "y") +
        labs(title = "European Population Statistics")  +
        opa() + scale_colour_hc() +
        scale_fill_brewer(palette = "Paired")   
      ggsave(paste0('./figures/', figure_name))
      return(data_plot[, plabelpos := NULL])
    }
  }

#control bootstrap based on peri sfs
cbs_positive <- std_bootstrap_control(maf_selinfo, region = c("positive"))

#peri load 
positive_load <- 
  eval_load(maf_selinfo[selection %in% c("positive")], 
            groups = c("POP"))

violin_plots(positive_load, cbs_positive,
             figure_name = 'htz.png')
