setwd("D:/Study/Graduation/GMHI_2020-master/")

library(vegan)
library(ape)
library(ade4)
library(ggplot2)
library(ggpubr)
library(easyGgplot2)
library(dplyr)

relative_abundance_file <- "./4347_final_relative_abundances.txt"
meta_data_file <- "./Final_metadata_4347.csv"
output <- "./final.csv"

df_output<-read.table(output, sep = ",", header = TRUE,row.names = 1)
Fig1b_input <- data.frame( t(df_output[,c(4:906)]))

Fig1c_d_input <- read.csv(meta_data_file, sep = ",", header = T,row.names = 1,check.names = F)
Fig1b_classified <- df_output[,-grep('unclassified', names(df_output))] # unclassified species removed

#colnames(Fig1b_classified)[2] <- "Phenotype_all"

Fig1b_classified$Phenotype_all <- gsub(x = Fig1b_classified$Type, 
                                  pattern = "[^Healthy].+|advanced adenoma", replacement = "Nonhealthy")
Fig1b_classified<-Fig1b_classified[c(789,1:788)]
trs <- function(x) asin(sqrt(x*0.01)) 

Fig1c_d_dataset1 <- data.frame(t(Fig1c_d_input),check.rows = F,check.names = F)
Fig1c_d_dataset <- Fig1c_d_dataset1[,c(1,7,15,15,33:ncol(Fig1c_d_dataset1))]
Fig1c_d_dataset[,-c(1:4)] <- lapply(Fig1c_d_dataset[,-c(1:4)], function(x) as.numeric(as.character(x)))

permanova1 <- Fig1b_classified %>% mutate_each(funs(trs), colnames(Fig1b_classified[,-c(1:4)]))
permanova1$study <- gsub(x = Fig1c_d_dataset$study, pattern = "_.+", replacement = "")


BC.dist=vegdist(permanova1[,c(5:789)], method = 'bray')
adonis2(BC.dist ~ Phenotype_all, data = permanova1, permutations = 999)

pcoa_all <- dudi.pco(BC.dist, scannf=FALSE, nf=3)
evals <- eigenvals(pcoa_all)
Variance <- evals / sum(evals)
Variance1 <- 100 * signif(Variance[1], 2)
Variance2 <- 100 * signif(Variance[2], 2)
Variance3 <- 100 * signif(Variance[3], 2)

pc_plot_data<-data.frame(pcoa_all$li$A1, pcoa_all$li$A2, permanova1$Phenotype_all,
                         Fig1c_d_dataset$Phenotype.1, Fig1c_d_dataset$`Sample Accession or Sample ID`)
colnames(pc_plot_data)<-c("x","y","Phenotype","Phenotype_all","Sample_ID")
pc_plot_data1 <- merge(pc_plot_data,aggregate(cbind(mean.x=x,mean.y=y)~Phenotype, 
                                              pc_plot_data,mean), by="Phenotype")
library(randomcoloR)
pdf(file = 'fig1a_out1.pdf')
Fig1a <- ggplot(pc_plot_data1, aes(x,y,color=factor(Phenotype)))+geom_point(size=1,alpha=0.9)+
  geom_point(aes(x=mean.x,y=mean.y),size=0)+stat_ellipse(level = 0.95)+
  labs(x = paste("PCoA1 (", Variance1, "%)", sep=""), y = paste("PCoA2 (", Variance2, "%)", sep=""),
       title = 'Fig1a')+theme_bw()+theme(legend.title=element_blank())+
  scale_color_manual(values=c("Healthy"="steelblue","Nonhealthy"="orange2"))
Fig1a+annotate("text", label = paste0(c("R-squared = ",0.00582,'\n',"p-value = ",0.001),collapse = ''), 
               x = -0.25, y = 0.35, color = "black")
dev.off()
pdf(file = 'fig1b_out1.pdf',height = 7,width = 10)
Fig1b <- ggplot(pc_plot_data1, aes(x,y,color=factor(Phenotype_all)))+geom_point(size=1,alpha=0.5)+
  geom_point(aes(x=mean.x,y=mean.y),size=0)+stat_ellipse(level = 0.95)+
  labs(x = paste("PCoA1 (", Variance1, "%)", sep=""), y = paste("PCoA2 (", Variance2, "%)", sep=""),
       title = 'Fig1b')+theme_bw()+theme(legend.title=element_blank())+
  scale_color_manual(values = c("Healthy"="steelblue",distinctColorPalette(12)))
Fig1b 
dev.off()
