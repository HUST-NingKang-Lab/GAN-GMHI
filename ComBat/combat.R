output <- "./metadata.csv"
df_output<-read.table(output, sep = ",", header = TRUE,row.names = 1)
eset_pca <- df_output
edata1<-eset_pca[,c(-1,-2,-3)]
phen01 <-eset_pca[,c(1,2,3)]
batch1 <- phen01$batch
mod1=model.matrix(~as.factor(Phenotype), data=phen01)
corrected_metadata=ComBat(dat=edata1,batch=phen01$batch1,mod=mod1)