std.boot <-  std_bootstrap_control(maf_selinfo, region = c("positive"), n_samples = 1000)
sfs.boot <-  sfs_bootstrap_control(maf_selinfo, region = c("positive"), n_samples = 1000)
gene_std.boot <- gene_std_bootstrap_control(maf_selinfo, region = c("positive"), n_samples = 1000)
gene_sfs.boot <- gene_sfs_bootstrap_control(maf_selinfo, region = c("positive"), n_samples = 1000)


    measures_to_plot = c("PdelPneu"
               # "PdelPS",
               #"Pd1Pneu",
               #"Pd1PS",
                # "Pd1"
                # "Pd",
               #"Pdel_Sift",
               #"PNPS",
               #"PdelPneu_SIFT",
                 #"htz"
            )
          
      set_data_plot <- function(data){
         melt(data, measure.vars = measures_to_plot, 
              variable.name = "measure", 
              value.name = "value", 
              varialble.factor = TRUE)[
                 ,.(POP, replica, measure, value)]
      }
      
      

      std <-
         set_data_plot(
            eval_load(std.boot, 
                      groups = c("POP", "replica"))
         )
      

      sfs <-
         set_data_plot(
            eval_load(sfs.boot, 
                      groups = c("POP", "replica"))
         )
      

      gene_std <-
         set_data_plot(
            eval_load(gene_std.boot, 
                      groups = c("POP", "replica"))
         )
      

      gene_sfs <-
         set_data_plot(
            eval_load(gene_sfs.boot, 
                      groups = c("POP", "replica"))
         )
      
      
std
std <- mutate(std, Boots = "STD-Boot")
std <- mutate(std, class.Boots = "SNP")
sfs
sfs <- mutate(sfs, Boots = "SFS-Boot")
sfs <- mutate(sfs, class.Boots = "SNP")

gene_std
gene_std <- mutate(gene_std, Boots = "STD-Boot-Gene")
gene_std <- mutate(gene_std, class.Boots = "Gene")

gene_sfs
gene_sfs <- mutate(gene_sfs, Boots = "SFS-Boot-Gene")
gene_sfs <- mutate(gene_sfs, class.Boots = "Gene")

boots_totalpop <- rbind(std, sfs, gene_std, gene_sfs)
boots <- filter(boots_totalpop, POP == "EUR")

ggplot(boots, aes(value,  fill = Boots, colour = Boots)) +
   geom_density(alpha = 0.2) + opa()


ggplot(boots, aes(value,  fill = class.Boots, colour = class.Boots)) +
   geom_density(alpha = 0.2, size = 1) + opa() +
   labs(title = "", x = "", y = "")  +
   theme(axis.text = element_text(size = rel(1.0)), legend.position = "right", plot.title = element_text(size = rel(0.5)) ,axis.title.y = element_text(size = rel(0.5)), axis.title.x = element_text(size = rel(0.5)),legend.title = element_text(size=16, face="bold")
         , legend.text=element_text(size=16), legend.key.size = unit(1.2, "cm"))
ggsave("./figures/bootstraps_dist/PdelPneu.png")



## HTZ
boots <- filter(boots, Boots ==  "STD-Boot" |Boots ==  "STD-Boot-Gene" )
ggplot(boots, aes(value,  fill = Boots, colour = Boots)) +
   geom_density(alpha = 0.2) + opa()

ggplot(boots, aes(value,  fill = class.Boots, colour = class.Boots)) +
   geom_density(alpha = 0.2, size = 1) + opa() +
   labs(title = "", x = "", y = "")  +
   theme(axis.text = element_text(size = rel(1.0)), legend.position = "right", plot.title = element_text(size = rel(0.5)) ,axis.title.y = element_text(size = rel(0.5)), axis.title.x = element_text(size = rel(0.5)),legend.title = element_text(size=16, face="bold")
         , legend.text=element_text(size=16), legend.key.size = unit(1.2, "cm"))
ggsave("./figures/bootstraps_dist/htz.png")
