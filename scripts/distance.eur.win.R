# Script para comparar regioes adjacentes com regioes muito distantes. 

distance.eur.win <- 
  function( data, 
            posicoes=HighD_POS,
            dist_max = 1000000,
            dist_min = 800000,
            cromossomos=HighD_Cros
  ) {
    
    
    #Preparing_data
    aba2 <- fread("~/Dropbox/laboratorio/data/Colonna_HighD.csv")
    aba2 <- filter(aba2 , VT == 'SNP')
    aba2 <- filter(aba2 , PAIR == 'AFR-EUR')
    aba2 <- filter(aba2 , POPwithHighestDAF == 'EUR')
    aba2 <- filter(aba2, GeneBiotype == "protein_coding")
    HighD_POS <- aba2$POS
    HighD_Cros <- aba2$CHR
    colonna_genes <- select(aba2, c(CHR, POPwithHighestDAF, HGNCsymbol))
    colonna_genes <- rename(colonna_genes,  POP = POPwithHighestDAF, GeneName = HGNCsymbol)
    
    #A
    lista_a<-list(NA)
    for(i in 1:length(posicoes)) 
    { lista_a[[i]]  <-filter(data, 
                             (CHR==cromossomos[i] & (POS>(posicoes[i]-dist_max)
                                                    & POS<(posicoes[i]-dist_min))) |
                               (CHR==cromossomos[i] & (POS>(posicoes[i]+dist_max)
                                                      & POS<(posicoes[i]+dist_min)))
    )   }
    
    a<-do.call(rbind,lista_a)
    genes <-unique(filter(select(a, c(CHR, POP, GeneName)), POP=="EUR")) 
    write.table(genes, "~/Dropbox/laboratorio/data/genes_interesse.tsv", row.names=F)
    
    return(genes)
  }


# Script para comparar regioes adjacentes com regioes muito distantes. 

distance.asian.win <- 
  function( data, 
            posicoes=HighD_POS,
            dist_max = 1000000,
            dist_min = 800000,
            cromossomos=HighD_Cros
  ) {
    
    #Preparing_data
    aba2 <- fread("~/Dropbox/laboratorio/data/colonna/aba2_EAS.csv")
    aba2 <- filter(aba2 , VT == 'SNP')
    aba2 <- filter(aba2, GeneBiotype == "protein_coding")
    aba2 <- filter(aba2 , PAIR == 'AFR-EAS')
    aba2 <- filter(aba2, highestDAF >  925)
    aba2 <- filter(aba2 , POPwithHighestDAF == 'EAS')
    HighD_POS <- aba2$POS
    HighD_Cros <- aba2$CHR
    colonna_genes <- select(aba2, c(CHR, POPwithHighestDAF, HGNCsymbol))
    colonna_genes <- rename(colonna_genes,  POP = POPwithHighestDAF, GeneName = HGNCsymbol)
    
    
    #A
    lista_a<-list(NA)
    for(i in 1:length(posicoes)) 
    { lista_a[[i]]  <-filter(data, 
                             (CHR==cromossomos[i] & (POS>(posicoes[i]-dist_max)
                                                     & POS<(posicoes[i]-dist_min))) |
                               (CHR==cromossomos[i] & (POS>(posicoes[i]+dist_max)
                                                       & POS<(posicoes[i]+dist_min)))
    )   }          
    a<-do.call(rbind,lista_a)
    Genes_a <-unique(filter(select(a, c(CHR, POP, GeneName)), POP=="EAS")) 
    
    write.table(genes, "~/Dropbox/laboratorio/data/asian_genes_interesse.tsv", row.names=F)
    
    return(genes)
  }



##### ##### ##### ##### ##### ##### 
##### ##### Analisando  ##### ##### 
##### ##### ##### ##### ##### ##### 

# Aqui vou gerar uma janela com eur.win e uma janela de mesmo tamanho com o distance.eur.win... A ideia eh comparalas, ou seja, rodar o script pra uma e depois pra outra. 
numero_genes <- length(unique(maf_selinfo$GeneName[maf_selinfo$selection=="positive"]))
a <- 
  eval_load(maf_selinfo[selection %in% c("positive")], 
            groups = c("POP"))
a$select<-rep("positive",nrow(a))

tab<-as.data.frame(matrix(NA,1,4))
result<-NULL
for(i in 1:5000){
  genes<- maf_selinfo[maf_selinfo$GeneName%in%sample(maf_selinfo$GeneName, numero_genes)]
  b <- 
    eval_load(genes[selection %in% c("control")], 
              groups = c("POP"))
  b$select<-rep("control",nrow(b))
  result<-rbind(result,b)
}
result

tab[1,1:4]<-(a$PdelPneu-as.numeric(tapply(result$PdelPneu,result$POP,mean)))/as.numeric(tapply(result$PdelPneu,result$POP,sd))

colnames(tab)<-unique(a$POP)
rownames(tab)<-paste0("z",c(1:nrow(tab)),sep="")




#eur.win(maf, n=125000, m=130000, p=350000)
positive_win <- tab
positive_win <- melt(positive_win, measure=1:4, value.factor=TRUE)
positive_win <- mutate(positive_win, under_selection = "EUR")
positive_win <- rename(positive_win, POP = variable)

#distance.eur.win(maf, dist_max = 350000, dist_min = 200000)
control_win <- tab
control_win  <- melt(control_win , measure=1:4, value.factor=TRUE)
control_win  <- mutate(control_win , under_selection = "EUR")
control_win  <- rename(control_win , POP = variable)

#distance.eur.win(maf, dist_max = 350000, dist_min = 200000)
control_win2 <- tab
control_win2  <- melt(control_win2 , measure=1:4, value.factor=TRUE)
control_win2  <- mutate(control_win2 , under_selection = "EUR")
control_win2  <- rename(control_win2 , POP = variable)







control_teste <- c(rep(control_win$value[1], 65), rep(control_win$value[2], 150), rep(control_win$value[3], 138), rep(control_win$value[4], 16))
hist(control_teste)




p<-ggplot(control_win, aes(x=POP, y=number)) 
  p + geom_histogram()

  
  
  
  control_win <- melt(control_win, measure=1:4, value.factor=TRUE)
  p<-ggplot(control_win, aes(x=variable, y=value)) 
  p + geom_point()
  
  
  distance.eur.win(maf, dist_max = 350000, dist_min = 200000)
  
  #eur.win(maf, n=125000, m=130000, p=350000)
  