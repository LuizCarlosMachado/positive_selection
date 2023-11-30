##############
##############
# DIFERENTES FREQUENCIAS 

data=maf_selinfo
data=filter(data,POP=="EUR")
data <- filter(data, MAF > 0)
data_ori=data
#parametros     
n_samples = 1000L    
region = c('positive')
delta_estat<-as.data.frame(matrix(NA,1,59))

ids_df <- data %>%
  distinct(GeneName, selection)

rle_lens <- rle(data$GeneName)
data <- as.data.frame(data)
for(i in 1:n_samples){
  # tpy<-tapply(data$selection,INDEX = data$POP,"sample")
  #data$selection<-tapply(data$selection,INDEX = data$POP,"sample")
  bla<-sample(c(1:length(rle_lens$values)),replace=F)
  gene2 <- rle_lens$values[bla]
  nreps <- rle_lens$lengths[bla]
  
  data2 <- data %>%
    mutate(gene2 = rep(gene2, nreps)) %>%
    left_join(ids_df, by = c("gene2" = "GeneName")) %>%
    rename(selection = selection.x, selection2 = selection.y)
  #sobre regiao positiva
  data2 <- as.data.table(data2)
  data_r <- data2[selection2 %in% region]
  estat_r<-eval_load(data_r, groups = c("POP"))
  data_g <- data2[selection2 %in% "control"]
  estat_g<-eval_load(data_g, groups = c("POP"))
  
  delta_estat[i,]<-as.numeric(estat_r)-as.numeric(estat_g)
  
} 
colnames(delta_estat)<-colnames(eval_load(data_g, groups = c("POP")))
#PdelPneu
#hist(delta_estat$PdelPneu,xlim = c(-0.06,0.06))
data_ori<-as.data.table(data_ori)
real<-eval_load(data_ori, groups = c("POP","selection"))
real<-real[,3:length(real)]
delta_real<-as.numeric(real[2,])-as.numeric(real[1,])

#teste <- select(eval_load(data_g, groups = c("POP")), -POP)
#teste <- eval_load(data_g, groups = c("POP"))
teste<-delta_estat[,-1]
delta_real<-t(as.data.frame(delta_real))
colnames(delta_real)<-colnames(teste)

teste <- melt(teste)
teste$tratamento <- rep("control")
teste$POP <- rep("EUR")
teste$POP_US <- rep("EUR")
teste$boot <- rep("std_gene")
teste$MAF <- rep("0")

delta_real <- melt(delta_real)
delta_real$tratamento <- rep("positive")
delta_real$POP <- rep("EUR")
delta_real$POP_US <- rep("EUR")
delta_real <- select(delta_real, - Var1)
delta_real <- rename(delta_real, variable = Var2 )
delta_real$boot <- rep("std_gene")
delta_real$MAF <- rep("0")


Delta_eur_0 <- rbind(delta_real, teste)

#write.table(Delta_total, "~/Dropbox/laboratorio/data/delta_total_easUS.csv")



data=maf_selinfo
data=filter(data,POP=="EUR")
data <- filter(data, MAF > 0.001)
data_ori=data
#parametros     
n_samples = 1000L    
region = c('positive')
delta_estat<-as.data.frame(matrix(NA,1,59))

ids_df <- data %>%
  distinct(GeneName, selection)

rle_lens <- rle(data$GeneName)
data <- as.data.frame(data)
for(i in 1:n_samples){
  # tpy<-tapply(data$selection,INDEX = data$POP,"sample")
  #data$selection<-tapply(data$selection,INDEX = data$POP,"sample")
  bla<-sample(c(1:length(rle_lens$values)),replace=F)
  gene2 <- rle_lens$values[bla]
  nreps <- rle_lens$lengths[bla]
  
  data2 <- data %>%
    mutate(gene2 = rep(gene2, nreps)) %>%
    left_join(ids_df, by = c("gene2" = "GeneName")) %>%
    rename(selection = selection.x, selection2 = selection.y)
  #sobre regiao positiva
  data2 <- as.data.table(data2)
  data_r <- data2[selection2 %in% region]
  estat_r<-eval_load(data_r, groups = c("POP"))
  data_g <- data2[selection2 %in% "control"]
  estat_g<-eval_load(data_g, groups = c("POP"))
  
  delta_estat[i,]<-as.numeric(estat_r)-as.numeric(estat_g)
  
} 
colnames(delta_estat)<-colnames(eval_load(data_g, groups = c("POP")))
#PdelPneu
#hist(delta_estat$PdelPneu,xlim = c(-0.06,0.06))
data_ori<-as.data.table(data_ori)
real<-eval_load(data_ori, groups = c("POP","selection"))
real<-real[,3:length(real)]
delta_real<-as.numeric(real[2,])-as.numeric(real[1,])

#teste <- select(eval_load(data_g, groups = c("POP")), -POP)
#teste <- eval_load(data_g, groups = c("POP"))
teste<-delta_estat[,-1]
delta_real<-t(as.data.frame(delta_real))
colnames(delta_real)<-colnames(teste)

teste <- melt(teste)
teste$tratamento <- rep("control")
teste$POP <- rep("EUR")
teste$POP_US <- rep("EUR")
teste$boot <- rep("std_gene")
teste$MAF <- rep("0.001")

delta_real <- melt(delta_real)
delta_real$tratamento <- rep("positive")
delta_real$POP <- rep("EUR")
delta_real$POP_US <- rep("EUR")
delta_real <- select(delta_real, - Var1)
delta_real <- rename(delta_real, variable = Var2 )
delta_real$boot <- rep("std_gene")
delta_real$MAF <- rep("0.001")


Delta_eur_01 <- rbind(delta_real, teste)

#write.table(Delta_total, "~/Dropbox/laboratorio/data/delta_total_easUS.csv")

data=maf_selinfo
data=filter(data,POP=="EUR")
data <- filter(data, MAF > 0.0011)
data_ori=data
#parametros     
n_samples = 1000L    
region = c('positive')
delta_estat<-as.data.frame(matrix(NA,1,59))

ids_df <- data %>%
  distinct(GeneName, selection)

rle_lens <- rle(data$GeneName)
data <- as.data.frame(data)
for(i in 1:n_samples){
  # tpy<-tapply(data$selection,INDEX = data$POP,"sample")
  #data$selection<-tapply(data$selection,INDEX = data$POP,"sample")
  bla<-sample(c(1:length(rle_lens$values)),replace=F)
  gene2 <- rle_lens$values[bla]
  nreps <- rle_lens$lengths[bla]
  
  data2 <- data %>%
    mutate(gene2 = rep(gene2, nreps)) %>%
    left_join(ids_df, by = c("gene2" = "GeneName")) %>%
    rename(selection = selection.x, selection2 = selection.y)
  #sobre regiao positiva
  data2 <- as.data.table(data2)
  data_r <- data2[selection2 %in% region]
  estat_r<-eval_load(data_r, groups = c("POP"))
  data_g <- data2[selection2 %in% "control"]
  estat_g<-eval_load(data_g, groups = c("POP"))
  
  delta_estat[i,]<-as.numeric(estat_r)-as.numeric(estat_g)
  
} 
colnames(delta_estat)<-colnames(eval_load(data_g, groups = c("POP")))
#PdelPneu
#hist(delta_estat$PdelPneu,xlim = c(-0.06,0.06))
data_ori<-as.data.table(data_ori)
real<-eval_load(data_ori, groups = c("POP","selection"))
real<-real[,3:length(real)]
delta_real<-as.numeric(real[2,])-as.numeric(real[1,])

#teste <- select(eval_load(data_g, groups = c("POP")), -POP)
#teste <- eval_load(data_g, groups = c("POP"))
teste<-delta_estat[,-1]
delta_real<-t(as.data.frame(delta_real))
colnames(delta_real)<-colnames(teste)

teste <- melt(teste)
teste$tratamento <- rep("control")
teste$POP <- rep("EUR")
teste$POP_US <- rep("EUR")
teste$boot <- rep("std_gene")
teste$MAF <- rep("0.0011")

delta_real <- melt(delta_real)
delta_real$tratamento <- rep("positive")
delta_real$POP <- rep("EUR")
delta_real$POP_US <- rep("EUR")
delta_real <- select(delta_real, - Var1)
delta_real <- rename(delta_real, variable = Var2 )
delta_real$boot <- rep("std_gene")
delta_real$MAF <- rep("0.0011")


Delta_eur_011 <- rbind(delta_real, teste)

#write.table(Delta_total, "~/Dropbox/laboratorio/data/delta_total_easUS.csv")


data=maf_selinfo
data=filter(data,POP=="EUR")
data <- filter(data, MAF > 0.002)
data_ori=data
#parametros     
n_samples = 1000L    
region = c('positive')
delta_estat<-as.data.frame(matrix(NA,1,59))

ids_df <- data %>%
  distinct(GeneName, selection)

rle_lens <- rle(data$GeneName)
data <- as.data.frame(data)
for(i in 1:n_samples){
  # tpy<-tapply(data$selection,INDEX = data$POP,"sample")
  #data$selection<-tapply(data$selection,INDEX = data$POP,"sample")
  bla<-sample(c(1:length(rle_lens$values)),replace=F)
  gene2 <- rle_lens$values[bla]
  nreps <- rle_lens$lengths[bla]
  
  data2 <- data %>%
    mutate(gene2 = rep(gene2, nreps)) %>%
    left_join(ids_df, by = c("gene2" = "GeneName")) %>%
    rename(selection = selection.x, selection2 = selection.y)
  #sobre regiao positiva
  data2 <- as.data.table(data2)
  data_r <- data2[selection2 %in% region]
  estat_r<-eval_load(data_r, groups = c("POP"))
  data_g <- data2[selection2 %in% "control"]
  estat_g<-eval_load(data_g, groups = c("POP"))
  
  delta_estat[i,]<-as.numeric(estat_r)-as.numeric(estat_g)
  
} 
colnames(delta_estat)<-colnames(eval_load(data_g, groups = c("POP")))
#PdelPneu
#hist(delta_estat$PdelPneu,xlim = c(-0.06,0.06))
data_ori<-as.data.table(data_ori)
real<-eval_load(data_ori, groups = c("POP","selection"))
real<-real[,3:length(real)]
delta_real<-as.numeric(real[2,])-as.numeric(real[1,])

#teste <- select(eval_load(data_g, groups = c("POP")), -POP)
#teste <- eval_load(data_g, groups = c("POP"))
teste<-delta_estat[,-1]
delta_real<-t(as.data.frame(delta_real))
colnames(delta_real)<-colnames(teste)

teste <- melt(teste)
teste$tratamento <- rep("control")
teste$POP <- rep("EUR")
teste$POP_US <- rep("EUR")
teste$boot <- rep("std_gene")
teste$MAF <- rep("0.002")

delta_real <- melt(delta_real)
delta_real$tratamento <- rep("positive")
delta_real$POP <- rep("EUR")
delta_real$POP_US <- rep("EUR")
delta_real <- select(delta_real, - Var1)
delta_real <- rename(delta_real, variable = Var2 )
delta_real$boot <- rep("std_gene")
delta_real$MAF <- rep("0.002")


Delta_eur_02 <- rbind(delta_real, teste)

#write.table(Delta_total, "~/Dropbox/laboratorio/data/delta_total_easUS.csv")




data=maf_selinfo
data=filter(data,POP=="EUR")
data <- filter(data, MAF > 0.003)
data_ori=data
#parametros     
n_samples = 1000L    
region = c('positive')
delta_estat<-as.data.frame(matrix(NA,1,59))

ids_df <- data %>%
  distinct(GeneName, selection)

rle_lens <- rle(data$GeneName)
data <- as.data.frame(data)
for(i in 1:n_samples){
  # tpy<-tapply(data$selection,INDEX = data$POP,"sample")
  #data$selection<-tapply(data$selection,INDEX = data$POP,"sample")
  bla<-sample(c(1:length(rle_lens$values)),replace=F)
  gene2 <- rle_lens$values[bla]
  nreps <- rle_lens$lengths[bla]
  
  data2 <- data %>%
    mutate(gene2 = rep(gene2, nreps)) %>%
    left_join(ids_df, by = c("gene2" = "GeneName")) %>%
    rename(selection = selection.x, selection2 = selection.y)
  #sobre regiao positiva
  data2 <- as.data.table(data2)
  data_r <- data2[selection2 %in% region]
  estat_r<-eval_load(data_r, groups = c("POP"))
  data_g <- data2[selection2 %in% "control"]
  estat_g<-eval_load(data_g, groups = c("POP"))
  
  delta_estat[i,]<-as.numeric(estat_r)-as.numeric(estat_g)
  
} 
colnames(delta_estat)<-colnames(eval_load(data_g, groups = c("POP")))
#PdelPneu
#hist(delta_estat$PdelPneu,xlim = c(-0.06,0.06))
data_ori<-as.data.table(data_ori)
real<-eval_load(data_ori, groups = c("POP","selection"))
real<-real[,3:length(real)]
delta_real<-as.numeric(real[2,])-as.numeric(real[1,])

#teste <- select(eval_load(data_g, groups = c("POP")), -POP)
#teste <- eval_load(data_g, groups = c("POP"))
teste<-delta_estat[,-1]
delta_real<-t(as.data.frame(delta_real))
colnames(delta_real)<-colnames(teste)

teste <- melt(teste)
teste$tratamento <- rep("control")
teste$POP <- rep("EUR")
teste$POP_US <- rep("EUR")
teste$boot <- rep("std_gene")
teste$MAF <- rep("0.003")

delta_real <- melt(delta_real)
delta_real$tratamento <- rep("positive")
delta_real$POP <- rep("EUR")
delta_real$POP_US <- rep("EUR")
delta_real <- select(delta_real, - Var1)
delta_real <- rename(delta_real, variable = Var2 )
delta_real$boot <- rep("std_gene")
delta_real$MAF <- rep("0.003")


Delta_eur_03 <- rbind(delta_real, teste)

#write.table(Delta_total, "~/Dropbox/laboratorio/data/delta_total_easUS.csv")




data=maf_selinfo
data=filter(data,POP=="EUR")
data <- filter(data, MAF > 0.004)
data_ori=data
#parametros     
n_samples = 1000L    
region = c('positive')
delta_estat<-as.data.frame(matrix(NA,1,59))

ids_df <- data %>%
  distinct(GeneName, selection)

rle_lens <- rle(data$GeneName)
data <- as.data.frame(data)
for(i in 1:n_samples){
  # tpy<-tapply(data$selection,INDEX = data$POP,"sample")
  #data$selection<-tapply(data$selection,INDEX = data$POP,"sample")
  bla<-sample(c(1:length(rle_lens$values)),replace=F)
  gene2 <- rle_lens$values[bla]
  nreps <- rle_lens$lengths[bla]
  
  data2 <- data %>%
    mutate(gene2 = rep(gene2, nreps)) %>%
    left_join(ids_df, by = c("gene2" = "GeneName")) %>%
    rename(selection = selection.x, selection2 = selection.y)
  #sobre regiao positiva
  data2 <- as.data.table(data2)
  data_r <- data2[selection2 %in% region]
  estat_r<-eval_load(data_r, groups = c("POP"))
  data_g <- data2[selection2 %in% "control"]
  estat_g<-eval_load(data_g, groups = c("POP"))
  
  delta_estat[i,]<-as.numeric(estat_r)-as.numeric(estat_g)
  
} 
colnames(delta_estat)<-colnames(eval_load(data_g, groups = c("POP")))
#PdelPneu
#hist(delta_estat$PdelPneu,xlim = c(-0.06,0.06))
data_ori<-as.data.table(data_ori)
real<-eval_load(data_ori, groups = c("POP","selection"))
real<-real[,3:length(real)]
delta_real<-as.numeric(real[2,])-as.numeric(real[1,])

#teste <- select(eval_load(data_g, groups = c("POP")), -POP)
#teste <- eval_load(data_g, groups = c("POP"))
teste<-delta_estat[,-1]
delta_real<-t(as.data.frame(delta_real))
colnames(delta_real)<-colnames(teste)

teste <- melt(teste)
teste$tratamento <- rep("control")
teste$POP <- rep("EUR")
teste$POP_US <- rep("EUR")
teste$boot <- rep("std_gene")
teste$MAF <- rep("0.004")

delta_real <- melt(delta_real)
delta_real$tratamento <- rep("positive")
delta_real$POP <- rep("EUR")
delta_real$POP_US <- rep("EUR")
delta_real <- select(delta_real, - Var1)
delta_real <- rename(delta_real, variable = Var2 )
delta_real$boot <- rep("std_gene")
delta_real$MAF <- rep("0.004")


Delta_eur_04 <- rbind(delta_real, teste)

#write.table(Delta_total, "~/Dropbox/laboratorio/data/delta_total_easUS.csv")




data=maf_selinfo
data=filter(data,POP=="EUR")
data <- filter(data, MAF > 0.005)
data_ori=data
#parametros     
n_samples = 1000L    
region = c('positive')
delta_estat<-as.data.frame(matrix(NA,1,59))

ids_df <- data %>%
  distinct(GeneName, selection)

rle_lens <- rle(data$GeneName)
data <- as.data.frame(data)
for(i in 1:n_samples){
  # tpy<-tapply(data$selection,INDEX = data$POP,"sample")
  #data$selection<-tapply(data$selection,INDEX = data$POP,"sample")
  bla<-sample(c(1:length(rle_lens$values)),replace=F)
  gene2 <- rle_lens$values[bla]
  nreps <- rle_lens$lengths[bla]
  
  data2 <- data %>%
    mutate(gene2 = rep(gene2, nreps)) %>%
    left_join(ids_df, by = c("gene2" = "GeneName")) %>%
    rename(selection = selection.x, selection2 = selection.y)
  #sobre regiao positiva
  data2 <- as.data.table(data2)
  data_r <- data2[selection2 %in% region]
  estat_r<-eval_load(data_r, groups = c("POP"))
  data_g <- data2[selection2 %in% "control"]
  estat_g<-eval_load(data_g, groups = c("POP"))
  
  delta_estat[i,]<-as.numeric(estat_r)-as.numeric(estat_g)
  
} 
colnames(delta_estat)<-colnames(eval_load(data_g, groups = c("POP")))
#PdelPneu
#hist(delta_estat$PdelPneu,xlim = c(-0.06,0.06))
data_ori<-as.data.table(data_ori)
real<-eval_load(data_ori, groups = c("POP","selection"))
real<-real[,3:length(real)]
delta_real<-as.numeric(real[2,])-as.numeric(real[1,])

#teste <- select(eval_load(data_g, groups = c("POP")), -POP)
#teste <- eval_load(data_g, groups = c("POP"))
teste<-delta_estat[,-1]
delta_real<-t(as.data.frame(delta_real))
colnames(delta_real)<-colnames(teste)

teste <- melt(teste)
teste$tratamento <- rep("control")
teste$POP <- rep("EUR")
teste$POP_US <- rep("EUR")
teste$boot <- rep("std_gene")
teste$MAF <- rep("0.005")

delta_real <- melt(delta_real)
delta_real$tratamento <- rep("positive")
delta_real$POP <- rep("EUR")
delta_real$POP_US <- rep("EUR")
delta_real <- select(delta_real, - Var1)
delta_real <- rename(delta_real, variable = Var2 )
delta_real$boot <- rep("std_gene")
delta_real$MAF <- rep("0.005")


Delta_eur_05 <- rbind(delta_real, teste)

#write.table(Delta_total, "~/Dropbox/laboratorio/data/delta_total_easUS.csv")

Delta_total_freq <- rbind(Delta_eur_0, Delta_eur_01,Delta_eur_011, Delta_eur_02, Delta_eur_021, Delta_eur_03,Delta_eur_04,Delta_eur_05)
Delta_total_freq <- as.data.table(Delta_total_freq)
Delta_total_freq <- rename(Delta_total_freq, treatment = tratamento)
Delta_total_freq[, pvalue := (.N - frank(value) + 1)/ .N, 
            by = .(POP, variable)]
Delta_total_freq[, plabelpos := max(value) + 0.2 * (max(value) - min(value)), by = .(variable)]
Delta_total_freq <- mutate(Delta_total_freq, is_significant = c(pvalue <= 0.05 | pvalue >= 1 - 0.05))

write.table(Delta_total_freq, "~/Dropbox/laboratorio/data/Delta_total_freq_eur.csv")


ggplot(filter(Delta_total_freq, variable == "PdelPneu"), aes(x=MAF, y = value)) + 
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
  #facet_wrap(~variable, scales="free_y", switch = "y")
  facet_wrap(~variable, scales="free_y") +
  geom_label(data= filter(Delta_total_freq, treatment=="positive", variable == "PdelPneu", POP_US == "EUR"),
             aes(x=MAF, 
                 y=plabelpos, 
                 label = format(round(log10(pvalue),2), nsmall = 2), 
                 fill=is_significant), size=label_size) +
  geom_point(data= filter(Delta_total_freq, treatment=="positive", variable == "PdelPneu", POP_US == "EUR"),
             aes(x=MAF, y=value), size= point_size, color='grey20') + opa()

###


Delta_total_freq


p<-ggplot(data = filter(Delta_total_freq, variable == "PdelPneu", treatment == "positive"), aes(x=MAF, y=value, group=treatment)) +
  geom_line(aes(color=treatment))+
  geom_point(aes(color=treatment))
p