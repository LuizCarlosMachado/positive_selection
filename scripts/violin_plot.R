library(ggthemes)
library(ggplot2)



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


violin_plots <- 
  function(load_region, bootstrap_control,
           measures_to_plot = c("pn",
                                "pd1",
                                "pd2",
                                "pb",
                                "ps",
                                "pd1pneu",
                                "pnps"),
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
             aes(x=POP, y = value)) + 
        geom_violin(fill='grey70') + 
        labs(x = "") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        geom_boxplot(fill = "grey70", width=0.1, outlier.size=0)  +
        geom_point(data = data_plot[is.na(replica)],
                   aes(x=POP, y=value), size= point_size, color='grey20') +
        geom_label(data = data_plot[is.na(replica)],
                   aes(x=POP, 
                       y=plabelpos, 
                       label = format(round(log10(pvalue),2), nsmall = 2), 
                       fill=is_significant), 
                   size=label_size) +
        facet_wrap(~measure, ncol=ncols, scales="free_y", switch = "y") +
        labs(title = "European Population Statistics")  +
        opa() + scale_colour_hc() +
        #scale_color_hue("groups") +
        #scale_fill_manual("População", values= c("grey70", "grey40", "grey90")) 
        scale_fill_manual(
          "", values=c( "grey40", "grey90"),
          guide=guide_legend(override.aes = list(colour=c("grey40", "grey90"))),
          labels=c("Non-Significant", "Significant")
        )
       # scale_fill_brewer(palette = "Paired")   
      ggsave(paste0('./figures/', figure_name))
      return(data_plot[, plabelpos := NULL])
    }
  }

#control bootstrap based on peri sfs
cbs_positive <- gene_sfs_bootstrap_control(maf_selinfo, region = c("positive"), n_samples = 3000)

#peri load 
positive_load <- 
  eval_load(maf_selinfo[selection %in% c("positive")], 
            groups = c("POP"))

violin_plots(positive_load, cbs_positive,
             figure_name = 'genic_violin.png')


###################################################################################################################################################################################################################################################   portugues    ################################################################################################################################################################


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


violin_plots <- 
  function(load_region, bootstrap_control,
           measures_to_plot = c(
             
              #"PdelPneu",
             # "PdelPS",
             #"Pd1Pneu",
              #"Pd1PS",
             #"Pd1",
             # "Pd",
             #"Pdel_Sift",
             #"PNPS",
             #"PdelPneu_SIFT",
              #"Pd1"
            
           ),
          pops = c("ALL", "AFR", "EUR", "EAS", "SAS"), 
           figure_name = 'violin.teste.png',
           th = 0.05,
           point_size = 5, 
           label_size = 5, 
           leg_position = 'none',
           ncols = 1){
    
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
             aes(x=POP, y = value)) + 
        geom_violin(fill='lightsteelblue3') + 
        #geom_violin(fill='#88AABB') + 
        #geom_violin(fill='#99AABB') + 
        # #99CCCC
        # darkseagreen
        labs(x = "") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        geom_boxplot(fill = "lightsteelblue3", width=0.1, outlier.size=0)  +
        #geom_boxplot(fill = "#88AABB", width=0.1, outlier.size=0)  +
        #geom_boxplot(fill = "#99AABB", width=0.1, outlier.size=0)  +
        
        geom_point(data = data_plot[is.na(replica)],
                   aes(x=POP, y=value), size= point_size, color='grey20') +
        geom_label(data = data_plot[is.na(replica)],
                   aes(x=POP, 
                       y=plabelpos, 
                       label = format(round(log10(pvalue),2), nsmall = 2), 
                       fill=is_significant), 
                   size=label_size) +
        facet_wrap(~measure, ncol=ncols, scales="free_y", switch = "y") +
        labs(title = "Hitchhiking Load")  +
        opa() + scale_colour_hc() +
        
        
        
        
        theme(axis.text = element_text(size = rel(1.5)),legend.position = "none", plot.title = element_text(size = rel(2.0)) , strip.text.y = element_text(size = rel(2.5)) ,axis.title.y = element_blank(), axis.title.x = element_text(size = rel(2.0)), legend.title = element_text(size=17, face="bold")
              , legend.text=element_text(size=17), legend.key.size = unit(1, "cm"))   +
        
        
        #scale_color_hue("groups") +
        #scale_fill_manual("População", values= c("grey70", "grey40", "grey90")) 
        scale_fill_manual(
          "Grupos", values=c( "#FF6633",  "seagreen3"),
          guide=guide_legend(override.aes = list(colour=c("#FF6633", "seagreen3"))),
          #   "", values=c( "grey40", "grey90"),
          #    guide=guide_legend(override.aes = list(colour=c("grey40", "grey90"))),
          labels=c("Não-Significante", "Significante")
        )
      # scale_fill_brewer(palette = "Paired")   
      
    }
  }

#control bootstrap based on peri sfs
cbs_positive <-  std_bootstrap_control(maf_selinfo, region = c("positive"), n_samples = 1000)
#cbs_positive <-  gene_sfs_bootstrap_control(maf_selinfo, region = c("positive"), n_samples = 1000)
#peri load 
positive_load <- 
  eval_load(maf_selinfo[selection %in% c("positive")], 
            groups = c("POP"))

violin_plots(positive_load, cbs_positive,
             measures_to_plot = c(
               "ps"
             #   "pn",
               #"pb",
               #"PNPS",
              # "Pdel_Sift",
               #"PdelPneu_SIFT"
               
               ))



ggsave("./figures/SIFT_gene_std_boot.png")
