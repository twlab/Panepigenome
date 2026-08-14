con<-read.table("/scratch/zdong/Projects/PanEpiG/V1-9/Model_Ind/Refmeth/continent.log")

library(dplyr)

con2 <- con %>%
  mutate(value = 1) %>%
  tidyr::pivot_wider(
    names_from  = V2,   # population column
    values_from = value,
    values_fill = 0
  )

pop<-read.table("/scratch/zdong/Projects/PanEpiG/V1-9/Model_Ind/Refmeth/pop.log")
pop<-pop[-which(pop$V2 %in% c("CH","EUR")),]

library(dplyr)

pop2 <- pop %>%
  mutate(value = 1) %>%
  tidyr::pivot_wider(
    names_from  = V2,   # population column
    values_from = value,
    values_fill = 0
  )



df<-read.table("/scratch/zdong/Projects/PanEpiG/V1-9/Model_Ind/NN/Saifur/part_00.txt",
               header=F)
he<-read.table("/scratch/zdong/Projects/PanEpiG/V1-9/Model_Ind/NN/Saifur/header.txt",header=T)
names(df)<-names(he)

df<-df[,c("gene",con$V1)]

### missing rate cutoff of 10%
df <- df[rowMeans(is.na(df[,2:ncol(df)])) <= 0.10, ]

### remove HG002 and HG005 as the undefined pop
df1 <- df[, !(colnames(df) %in% c("HG002","HG005"))]

r2     <- numeric(nrow(df))
adj_r2 <- numeric(nrow(df))
r2_pop     <- numeric(nrow(df))
adj_r2_pop <- numeric(nrow(df))

for (i in seq_len(nrow(df))) {
  
  y <- as.numeric(df[i, 2:ncol(df)])
  
  if (length(y) == 0 || all(is.na(y))) {
    r2[i] <- NA
    adj_r2[i] <- NA
  } else {
    fit <- lm(y ~ con2$EUR + con2$EAS + con2$AMR + con2$AFR + con2$SAS)
    r2[i] <- summary(fit)$r.squared
    adj_r2[i] <- summary(fit)$adj.r.squared
  }
  
  y <- as.numeric(df1[i, 2:ncol(df1)])
  
  if (length(y) == 0 || all(is.na(y))) {
    r2_pop[i] <- NA
    adj_r2_pop[i] <- NA
  } else {
    fit <- lm(y ~ pop2$GBR + pop2$FIN + pop2$CHS + pop2$PUR + pop2$CLM + pop2$IBS +
                pop2$ACB + pop2$PEL + pop2$KHV + pop2$CDX + pop2$GWD + pop2$PJL +
                pop2$ESN + pop2$MSL + pop2$STU + pop2$ITU + pop2$BEB + pop2$YRI +
                pop2$CHB + pop2$JPT + pop2$LWK + pop2$MXL + pop2$ASW +
                pop2$TSI + pop2$GIH + pop2$MKK)
    r2_pop[i] <- summary(fit)$r.squared
    adj_r2_pop[i] <- summary(fit)$adj.r.squared
  }
}

re<-data.frame(df[,1],r2,adj_r2,r2_pop,adj_r2_pop)
write.table(re,file="r2_00.txt",sep="\t",row.names = F,quote = F)




