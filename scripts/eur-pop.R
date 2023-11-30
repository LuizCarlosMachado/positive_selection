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

## Main functions ##

aba2 <- fread("~/Dropbox/laboratorio/data/Colonna_HighD.csv")
aba2 <- filter(aba2 , POPwithHighestDAF == 'EUR')
aba2 <- filter(aba2 , VT == 'SNP')
aba2 <- filter(aba2 , PAIR == 'AFR-EUR')
aba2 <- filter(aba2, GeneBiotype == "protein_coding")
HighD_POS <- aba2$POS
HighD_Cros <- aba2$CHR
#aba2_gene <- filter(aba2, GeneBiotype == "protein_coding")
colonna_genes <- select(aba2, c(CHR, POPwithHighestDAF, HGNCsymbol))
colonna_genes <- rename(colonna_genes,  POP = POPwithHighestDAF, GeneName = HGNCsymbol)

eur.win <- 
  function( data, 
            n=100000,
            m=120000,
            p=200000, 
            posicoes=HighD_POS,
            cromossomos=HighD_Cros
  ) {
    
    
    
    #A
    lista_a<-list(NA)
    for(i in 1:length(posicoes)) 
    { lista_a[[i]]  <-filter(data, 
                             (CHR==cromossomos[i] & POS>(posicoes[i]-n) & POS<(posicoes[i]+n))) 
    }             
    a<-do.call(rbind,lista_a)
    Genes_a <-unique(filter(select(a, c(CHR, POP, GeneName)), POP=="EUR")) 
    
    #C
    lista_b<-list(NA)
    for(i in 1:length(posicoes)) 
    { lista_b[[i]]  <-filter(data, 
                             ((CHR==cromossomos[i] & POS>(posicoes[i]+m) & POS<(posicoes[i]+p) )| (CHR==cromossomos[i] & POS>(posicoes[i]-p) & POS<(posicoes[i]-m))))
    }
    b<-do.call(rbind,lista_b)
    Genes_b <-unique(filter(select(b, c(CHR, POP, GeneName)), POP=="EUR")) 
    
    #A^C
    genes_interesse <- Genes_a[Genes_a$GeneName%in%Genes_b$GeneName==F]
    genes <- unique(rbind(colonna_genes, genes_interesse))
    
    write.table(genes, "~/Dropbox/laboratorio/data/genes_interesse.tsv", row.names=F)
    
    return(genes)
  }
ntile <- 
  function(x, n = 4){
    # Just like dplyr::ntile, if it  used data.table::frank it
    # shold be faster, but it is not.....
    floor((n * (rank(x, ties.method = "first", na.last = "keep") - 1)/ length(x)) + 1)
  }


eval_maf <- 
  function(data){
    data[nAL > 2, MAF := AF[which.min(AF)], by= .(POP, CHR, POS)][
      is.na(MAF), MAF := AF][
        MAF != 1 & MAF == AF][
          MAF > 0.5, MAF := 1 - MAF][
            , AF := NULL][
              , logMAF := log10(MAF)]
  }

set_snv_selection2 <- 
  function(data, pop_info) {
    
    # Analysis removing bitarello genes 
    bitarello_db <- 
      fread('~/Dropbox/laboratorio/data/bitarello.genes2.txt')
    bitarello_genes <- 
      bitarello_db[,V1]
    
    classics <- 
      union(
        paste0("HLA-", c("A", "B", "C", "DQB1", "DQA1", "DQB1", "DPA1", "DPB1", "DRA", "DRB1")),
        grep("HLA-", bitarello_genes, value = TRUE))
    
    peri_db  <- 
      fread('~/Dropbox/laboratorio/data/peri_transcripts_25kb_20sd.tsv')
    
    peri_genes <- peri_db[POP == 'ALL' & 
                            !(GeneName %in% classics) &
                            !(GeneName %in% bitarello_genes), 
                          GeneName]
    z <-
      unique(fread('~/Dropbox/laboratorio/data/genes_interesse.tsv'))
    
    
    positive_Luiz <- 
      z[,GeneName]
    
    data[GeneName %in% bitarello_genes, selection := 'balanced'][
      GeneName %in% positive_Luiz, selection := 'positive'][
        GeneName %in% peri_genes, selection := 'peri'][
          GeneName %in% classics, selection := 'classics'][
            CHR == 6  & 
              POS >= 29640169 & 
              POS <= 33115544 & 
             is.na(selection), 
            selection := 'MHC'][
              CHR == 6  & 
                POS >= data[GeneName == 'SCGN', min(POS)] & 
                POS <= data[GeneName == "SYNGAP1", max(POS)] & 
                is.na(selection), 
              selection := 'xMHC'][
                is.na(selection), selection := 'control']
    
    data <- merge(data, pop_info[, .(POP, GRP)], 
                  by = c('POP'),
                  allow.cartesian = TRUE, 
                  all.x = TRUE)
    
    data[is.na(GRP), GRP := POP]
    
    return(data)
  } 

set_snv_selection <- 
function(data, selected_genes, exons = FALSE) {
  if(exons == FALSE){
    z <-
      unique(fread('~/Dropbox/laboratorio/data/genes_interesse.tsv'))
    
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

sfs_bootstrap_control<- 
  function(data, 
           n_samples = 3000L, 
           region = c('positive'), 
           n_bins = 12 
  ){# booststrap respecting site freq. spectrum
    
    data[, bins := cut(logMAF, n_bins, labels=FALSE), 
         keyby = POP][, sampling_info := paste(POP, bins)]
    
    sfs_region <-
      data[selection %in% region,.(
        n = .N), 
        #min = min(MAF), 
        #max = max(MAF)), 
        keyby = .(POP, bins, sampling_info)]
    
    # Some time there are no snps in bins: make ns = 0 for theses
    ns <- n_samples * sfs_region[,n]
    names(ns) <- sfs_region[,sampling_info]
    skeys <- data[, unique(sampling_info)]
    miss_keys <- setdiff(skeys, sfs_region[,sampling_info])
    miss <- rep(0, length(miss_keys))
    names(miss) <- miss_keys
    ns <- c(ns, miss)
    
    bootstrap_control <- 
      stratified(data[selection == "control"], 
                 "sampling_info", 
                 replace = TRUE,
                 ns)[, 
                     replica := 1:n_samples, 
                     by = POP]
    
    data[, `:=`(bins = NULL, sampling_info = NULL)]
    return(bootstrap_control)
  }

gene_sfs_bootstrap_control<- 
  function(data, 
           n_samples = 1000, 
           region = c('case'), 
           n_bins = 12 
  ){# booststrap respecting site freq. spectrum
    
    data[, bins := cut(logMAF, n_bins, labels=FALSE), 
         keyby = POP][, sampling_info := paste(POP, bins)]
    
    sfs_region <-
      data[selection %in% region,.(
        n = .N), 
        keyby = .(POP, bins, sampling_info)]
    
    # Some time there are no snps in bins: make ns = 0 for theses
    ns <- sfs_region[,n]
    names(ns) <- sfs_region[,sampling_info]
    skeys <- data[, unique(sampling_info)]
    miss_keys <- setdiff(skeys, sfs_region[,sampling_info])
    miss <- rep(0, length(miss_keys))
    names(miss) <- miss_keys
    ns <- c(ns, miss)
    
    if(length(miss) > 0){
      miss_sfs <- data.table(sampling_info = miss_keys, n = 0)
      miss_sfs[, c("POP", "bins") := tstrsplit(sampling_info, " ")] 
      miss_sfs[, bins := as.numeric(bins)]
      
      sfs_region <-
        rbindlist(list(sfs_region, 
                       miss_sfs[, (names(sfs_region)), with = FALSE]))
    }
    #sfs_region[order(POP, bins)]
    
    gene_names <- 
      sample(data[selection == 'control', unique(GeneName)])
    #just to take in the begining of the list
    mg <- as.integer(length(gene_names) / 3)  
    
    num_g <- length(data[selection %in% region, unique(GeneName)])
    
    setkey(data, GeneName)
    
    get_genes <- 
      function(id){
        gn <- sample(gene_names[c(1:mg)], 1)
        gn_id  <- which(gene_names ==  gn)
        sfs_ok <- TRUE
        i <- 0
        while(sfs_ok){
          g_l <- gene_names[gn_id:(gn_id + as.integer(0.9 * num_g) + i)]
          data_bs <- data[g_l]
          
          msfs <- 
            merge(sfs_region, data_bs[, .N, keyby = .(POP, bins)], 
                  by = c("POP", "bins"), 
                  all.x = TRUE)
          
          sfs_ok <- 
            msfs[,any(is.na(N)) || 
                   any((n - N) > 0, na.rm = TRUE)]
          
          i <- i + 1
        }
        data_bs[, replica := id]
        return(data_bs)
      }
    
    bootstrap_control <- 
      rbindlist(
        lapply(1:n_samples, 
               function(x)
                 stratified(get_genes(x), 
                            "sampling_info", 
                            replace = FALSE,
                            ns)
        ))
    
    data[, `:=`(bins = NULL, sampling_info = NULL)]
    return(bootstrap_control)
  }

std_bootstrap_control<- 
  function(data, 
           n_samples = 1000L, 
           region = c('positive')
  ){# standard booststrap
    
    data_r <- data[selection %in% region]
    ns <- data_r[,.N * n_samples, by = POP][, V1]
    names(ns) <- data_r[, unique(POP)]
    
    control_b <- 
      stratified(data[selection == "control"], 
                 "POP", 
                 replace = TRUE,
                 ns)[, 
                     replica := 1:n_samples, 
                     by = POP]
    
    return(control_b)
  }

gene_std_bootstrap_control<- 
  function(data, 
           n_samples = 1000L, 
           region = c('positive')
  ){# standard booststrap
    
    data_r <- data[selection %in% region]
    ns <- data_r[, .N, by = POP]
    l_ns <-  ns[, N]
    names(l_ns) <- ns[, POP]
    gene_names <- sample(data[selection == 'control', unique(GeneName)])
    
    #just to take in the begining of the list
    mg <- as.integer(length(gene_names) / 3)  
    num_g <- length(data[selection %in% region, unique(GeneName)])
    setkey(data, GeneName)
    
    get_genes <- 
      function(id){
        gn <- sample(gene_names[c(1:mg)], 1)
        gn_id  <- which(gene_names ==  gn)
        
        nsnp_ok <- TRUE
        i <- 0
        while(nsnp_ok){
          g_l <- gene_names[gn_id:(gn_id + as.integer(0.9 * num_g) + i)]
          data_bs <- data[g_l]
          
          m_ns <- 
            merge(ns, data_bs[, .N, keyby = POP], 
                  by = c("POP"), )
          
          nsnp_ok <- 
            m_ns[, any((N.x - N.y) > 0)]
          
          i <- i + 1
        }
        data_bs[, replica := id]
        return(data_bs)
      }
    
    control_b <- 
      rbindlist(
        lapply(1:n_samples, 
               function(x)
                 stratified(get_genes(x), 
                            "POP", 
                            replace = FALSE,
                            l_ns)
        ))
    
    return(control_b)
  }

violin_plots <- 
  function(load_region, bootstrap_control,
           measures_to_plot = c("pnps", "w_pnps", 
                                "pd1pneu", "w_pd1pneu", 
                                "pd1", "wd1", 
                                "PHRED", "w_PHRED"),
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
        theme_hc() + scale_colour_hc() +
        scale_fill_brewer(palette = "Paired")   
      ggsave(paste0('./figures/', figure_name))
      return(data_plot[, plabelpos := NULL])
    }
  }

full_sfs <- 
  function(data,
           figure_name = 'sfs_teste.png', 
           n_bins = 12, 
           leg_position = 'right',
           pd = 0.05, 
           sel = c("control", "positive"), 
           ncols = 2){
    
    data_plot <- copy(data[selection %in% sel])
    data_plot[, bins := cut(logMAF, n_bins, label = FALSE)][
      , logMAF := max(logMAF), by = .(bins, POP)]
    
    ggplot(transform(data_plot, selection = factor(selection, levels = sel)), 
           aes(x=logMAF, fill = selection)) + 
      geom_bar(aes(y = ..prop..),
               position = 'dodge') +
      facet_wrap(~POP, ncol=ncols, scales="free_y") +
      theme(legend.position=leg_position)
    
    ggsave(paste0('./figures/', figure_name))
  }

load_per_bin<-
  function(data, 
           figure_name = 'teste_bin_load.png',
           n_bins  = 5,
           pop = "ALL", 
           measures_to_plot = c(
             "pbps", 
             "pdpn", 
             "pdbpn",
             "pnps", 
             "pdelpneu",
             "pdpb",
             "PHRED"),
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
      theme(legend.position=leg_position)
    ggsave(paste0('./figures/', figure_name))
    
    data[, bins := NULL]
  }

rank_freq_fold <- 
  function(data, 
           region = c('postive'),
           n_samples = 3000L 
  ){
    
    data_r <- data[selection %in% region][, replica := NA_integer_]
    
    control_b <- 
      std_bootstrap_control(data, 
                            region = region, 
                            n_samples = n_samples)
    data_p <- 
      rbindlist(list(data_r, control_b))
    
    data_p[, rank := frank(MAF, ties.method = 'first'), by = .(replica, POP)]
    
    lom <- dcast(data_p,  replica ~ POP + rank, value.var = c("MAF"))
    
    M <- lom[, 2:dim(lom)[2], with = FALSE]
    ground <- unlist(lom[is.na(replica), 2:dim(lom)[2], with = FALSE])
    
    lom[, 2:dim(lom)[2]] <- sweep(1 / M, 2, ground, FUN = "*")
    
    drank_MAF <- melt(lom, measure.vars = 2:dim(lom)[2], 
                      value.name = 'diff_MAF') 
    
    drank_MAF[, c("POP", "rank") := tstrsplit(variable, "_")]
    drank_MAF[, rank := as.numeric(rank)][, variable := NULL]
    
    merge(
      data_r[, .(n_snp = .N), by = POP],
      drank_MAF[!is.na(replica), .(mean(diff_MAF)), by = .(replica, POP)][, 
                                                                          .(mean_prop = mean(V1), 
                                                                            sd_prop = sd(V1), 
                                                                            median = median(V1),
                                                                            lower = quantile(V1, c(0.025)),
                                                                            upper = quantile(V1, c(0.975))
                                                                          ), by = POP],
      all = TRUE, 
      by = c("POP"))
  }


mean_freq_fold <- 
  function(data, 
           region = c('postive'),
           n_samples = 3000L 
  ){
    
    data_r <- data[selection %in% region][, replica := NA_integer_]
    
    control_b <- 
      std_bootstrap_control(data, 
                            region = region, 
                            n_samples = n_samples)
    data_p <- 
      rbindlist(list(data_r, control_b))
    
    m <- data_p[, .(n_snp = .N, mean_MAF = mean(MAF)), by = .(replica, POP)]
    
    m[, pvalue := (.N - frank(mean_MAF) + 1)/ .N, by = .(POP)]
    
    x <- dcast(m, replica ~ POP, value.var = c("mean_MAF"))
    
    ground <- m[is.na(replica), mean_MAF]
    
    x[, 2:dim(x)[2]] <- 
      sweep( 1 / x[, 2:dim(x)[2], with  = FALSE], 2, 
             ground, FUN = "*")
    
    y <- melt(x, measure.vars = 2:dim(x)[2],
              variable.name = 'POP',
              value.name = 'fold')
    
    z <- merge(m[is.na(replica)],
               m[!is.na(replica),
                 .(mean_MAF_control = mean(mean_MAF),
                   err_MAF_control = sd(mean_MAF)),
                 by = POP],
               by = c("POP"), all = TRUE)
    
    zz <- merge( z, 
                 y[!is.na(replica), 
                   .(mean_fold = mean(fold), err_fold  = sd(fold)), 
                   by = POP], 
                 by = c("POP"), all = TRUE)
    zz[, replica := NULL]
    
    return(list(mean_ff = zz, bs_data_ff = y))
  }


ff_bar_plots <- 
  function(data,
           pops = c("EUR"), 
           rgs = c('positive'),
           figure_name = 'freq_fold_teste.png',
           th = 0.05,
           type = 'del', 
           point_size = 2, 
           label_size = 3, 
           leg_position = 'top',
           ncols = 3){
    
    data_ffold <- copy(data[type == 'del' & region %in% rgs])
    
    data_ffold[, `:=`(region = factor(region, levels = rgs),
                      POP = factor(POP, levels = pops))]
    f <- 
      function(x){
        r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
        names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
        return(r)
      }
    
    box_markers <- 
      data_ffold[, .(variable = c("ymin", "lower", "middle", 
                                  "upper", "ymax"), 
                     value = f(fold)),
                 by = .(POP, region)]
    
    bx <- dcast(box_markers, POP + region ~ variable)
    
    ggplot(data =  bx,
           aes(x=POP, fill=POP)) + 
      geom_boxplot(stat = "identity",
                   aes(ymin = ymin, 
                       lower = lower, 
                       middle = middle, 
                       upper = upper, 
                       ymax = ymax)
      ) +
      geom_hline(yintercept = 1) +
      theme(legend.position=leg_position,
            axis.text.x = element_blank(),
            panel.margin = unit(2, "lines")) +
      facet_wrap(~region, ncol=ncols, switch = "x") 
    ggsave(paste0('./figures/', figure_name))
  }


eval_freq_fold <- 
  function(region = c('positive'), data){
    
    del_ffold <- 
      mean_freq_fold(data[PolyPhenCat == 'probably_damaging'],
                     region = region)
    del_ffold[[1]][, type := 'del']
    del_ffold[[2]][, type := 'del']
    
    neu_ffold <- 
      mean_freq_fold(data[PolyPhenCat != 'probably_damaging'],
                     region = region)
    neu_ffold[[1]][, type := 'neu']
    neu_ffold[[2]][, type := 'neu']
    
    mean_ffold <- rbindlist(list(del_ffold[[1]], neu_ffold[[1]]))
    bs_ffold <- rbindlist(list(del_ffold[[2]], neu_ffold[[2]]))
    
    return(list(mean_ffold = mean_ffold,
                bs_ffold = bs_ffold))
  }

pairs_wc_tests <-
  function(data){
    x <- data[POP == "ALL" & !(selection %in% c("positive"))][, logMAF]
    g <- data[POP == "ALL" & !(selection %in% c("positive"))][, selection]
    pairwise.wilcox.test(x,g, paired = FALSE)
  }


set_regions <- 
  function(){
    regions <- list(c('positive'))
    names(regions) <- 
      unlist(lapply(regions, function(X) paste0(X, collapse = "_")))
    
    return(regions)
  }

get_load_and_plot <-
  function(data,
           n_samples = 3000L, 
           region = c('positive'), 
           n_bins = 12,
           measures_to_plot = c("pnps", "w_pnps", 
                                "pd1pneu", "w_pd1pneu", 
                                "pd1", "wd1", 
                                "PHRED", "w_PHRED"),
           pops = c("ALL", "AFR", "EUR", "EAS", "SAS"), 
           figure_name = 'violin.teste.png',
           th = 0.05,
           point_size = 2, 
           label_size = 3, 
           leg_position = 'none',
           ncols = 2
  ){
    
    if( figure_name != 'none' | figure_name != 'violin.teste.png'){
      figure_name <-  paste0(c('load_', 
                               paste0(region, collapse = '-'), 
                               '.png'), 
                             collapse = "")
    }
    
    load_region <-
      eval_load(data[selection %in% region], groups = c('POP'))
    
    control_b <- 
      sfs_bootstrap_control(data, 
                            region = region,
                            n_samples = n_samples, 
                            n_bins = n_bins)
    
    data_plot <- 
      violin_plots(load_region, control_b, 
                   measures_to_plot = measures_to_plot, 
                   pops = pops,
                   figure_name = figure_name, 
                   th = th, 
                   point_size = point_size,
                   label_size = label_size, 
                   leg_position = leg_position,
                   ncols = ncols)
    
    x <- dcast(data_plot, 
               replica ~ POP + measure, 
               value.var = c("value"),
               sep = '-')
    
    ground <- 
      unlist(x[is.na(replica), 2:dim(x)[2], with = FALSE])
    
    x[, 2:dim(x)[2]] <- 
      sweep(1 / x[, 2:dim(x)[2], with  = FALSE], 2, 
            ground, FUN = "*")
    
    y <- melt(x, measure.vars = 2:dim(x)[2],
              value.name = 'fold')
    
    y[, c("POP", "measure") := tstrsplit(variable, '-')][,
                                                         variable := NULL]
    
    z <- merge(data_plot[is.na(replica)],
               data_plot[!is.na(replica),
                         .(mean_control = mean(value),
                           err_control = sd(value)),
                         by = .(POP, measure)],
               by = c("POP", "measure"), 
               all = TRUE)
    
    zz <- merge(z, 
                y[!is.na(replica), 
                  .(mean_fold = mean(fold), 
                    err_fold = sd(fold)),
                  by = .(POP, measure)],
                by = c("POP", "measure"), 
                all = TRUE)
    zz[, replica := NULL]
    
    return(zz)
  }



scatter_control <- 
  function(bootstrap_control, 
           figure_name = 'teste_scatter.png',
           leg_position = 'left'){
    
    pop_groups <- 
      readr::read_delim("pops.1kg.tsv", delim=";") %>%
      dplyr::select(POP, GRP)
    
    control_data_plot <-
      bootstrap_control %>%
      #dplyr::group_by(POP, replica, bins = cut(MAF, seq(0, 0.5, 0.05), labels=FALSE)) %>%
      dplyr::group_by(POP, replica) %>%
      eval_load() %>%
      dplyr::left_join(pop_groups)
    
    #ggplot(control_data_plot) +
    #geom_point(aes(x = pnps, y = PHRED, color=GRP)) 
    
    ggpairs(control_data_plot, 
            columns = grep("(GRP)|(PHRED)|(pnps)|(pdpb)", colnames))
    
    ggsave(paste0('./figures/', figure_name))
  }


##########################################################################################
#################################### Read fulldata  ######################################
##########################################################################################


maf <-fread('~/Dropbox/laboratorio/data/all.pop.maf.tsv')
maf <- filter(maf, POP == "EUR")
## HighD + Vizinhos
# Preparando janela
eur.win(maf, n=145000, m=160000, p=350000)

##maf_selinfo
system.time(
  maf_selinfo <- set_snv_selection(maf, selected_genes = pop_info)
)
maf_selinfo<-as.data.table(mutate(maf_selinfo, htz = 2*(MAF)*(1-MAF)))
maf_selinfo <- filter(maf_selinfo, POP == "EUR")

#rm(maf)
#genes_load, cada linha é um gene, e as colunas sao as informacoes desse gene (pn, ps,...)
genes_load <- eval_load(maf_selinfo, groups = c("selection", "GeneName", "POP"))[order(selection, pdbpn, n_snp)]


#Bad_genes é pra ver quais genes nao estao descritos pelo PolyPhen, pdbpn = Snps_Pol / Pn, nesse script, todos sao bons 
#bad_genes <- genes_load[ selection != 'control' & (pdbpn < 0.95 |pdbpn == "NaN"), GeneName]
#bad_genes <- genes_load[ selection != 'control' & pnps == Inf , GeneName]
#bad_genes <- genes_load[ selection != 'control' &  GeneName == "CELSR1", GeneName]


# well annotated genes data (good genes)
maf_selinfo <- 
  maf_selinfo[!(GeneName %in% bad_genes)]
table(maf_selinfo$selection)
