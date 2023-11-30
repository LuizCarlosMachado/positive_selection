# Teste de Fisher para pnps e para pdelpneu
# PnPs


table((maf_selinfo$Consequence)[maf_selinfo$selection == "positive"])
table((maf_selinfo$Consequence)[maf_selinfo$selection == "control"])
Ps <- c(381,  75669)
Pn <- c(483 , 115422)
cont_table <- as.matrix(data.frame(Ps, Pn, row.names = c("Positive", "Control")))
fisher.test(cont_table, alternative = "greater")



#MAF > 0.001
maf_selinfo2 <- filter(maf_selinfo, MAF > 0.001)
table((maf_selinfo2$Consequence)[maf_selinfo2$selection == "positive"])
table((maf_selinfo2$Consequence)[maf_selinfo2$selection == "control"])


Ps <- c(202, 48697)
Pn <- c(238, 59689)
cont_table2 <- as.matrix(data.frame(Ps, Pn, row.names = c("Positive", "Control")))
fisher.test(cont_table2, alternative = "two.sided")


##############################
########################### Para classificações do Polyphen 
#PdelPneu = Pd1+ Pd2 / Pb +Ps
# Pd1_Pd2 = Pd1 + Pd2
# Pb_Ps = Pb +Ps
table((maf_selinfo$PolyPhenCat)[maf_selinfo$selection == "positive"])
table((maf_selinfo$PolyPhenCat)[maf_selinfo$selection == "control"])
table((maf_selinfo$Consequence)[maf_selinfo$selection == "positive"])
table((maf_selinfo$Consequence)[maf_selinfo$selection == "control"])
Pd1 <- c(144 , 25917)
Ps_Pb <- c(637, 146085)
cont_table2 <- as.matrix(data.frame(Pd1, Ps_Pb, row.names = c("Positive", "Control")))
fisher.test(cont_table2, alternative = "greater")





#MAF > 0.0019
maf_selinfo2 <- filter(maf_selinfo, MAF > 0.0019)
table((maf_selinfo2$PolyPhenCat)[maf_selinfo2$selection == "positive"])
table((maf_selinfo2$PolyPhenCat)[maf_selinfo2$selection == "control"])
table((maf_selinfo2$Consequence)[maf_selinfo2$selection == "positive"])
table((maf_selinfo2$Consequence)[maf_selinfo2$selection == "control"])
Pd1 <- c(53, 8226)
Ps_Pb <- c(312, 73270)
cont_table2 <- as.matrix(data.frame(Pd1, Ps_Pb, row.names = c("Positive", "Control")))
fisher.test(cont_table2, alternative = "greater")










