library(ggplot2)
library(ggpubr)
library(vegan)
library(rcompanion)
final <- read.csv('GMHI.csv')

fig4a <- ggscatterhist(
  final, x = "R_alpha", y = "R_GMHI",
  color = "Phenotype",
  palette = c("steelblue", "orange2"),
  margin.plot = "histogram",
  ggtheme = theme_bw(),
  show.legend = FALSE,xlab = "Shannon(RAW)",ylab = "GMHI(RAW)",
  margin.params = list(fill = "Phenotype", color = "black", size = 0.1)
)

fig4b <- ggscatterhist(
  final, x = "G_alpha", y = "G_GMHI",
  color = "Phenotype",
  palette = c("steelblue", "orange2"),
  margin.plot = "histogram",
  ggtheme = theme_bw(),
  show.legend = FALSE,xlab = "Shannon(GAN)",ylab = "GMHI(GAN)",
  margin.params = list(fill = "Phenotype", color = "black", size = 0.1)
)

fig4c<-ggscatterhist(
  final, x = "R_GMHI", y = "G_GMHI",
  color = "Phenotype",
  palette = c("steelblue", "orange2"),
  margin.plot = "histogram",
  ggtheme = theme_bw(),
  show.legend = FALSE,xlab = "RAW_GMHI",ylab = "GAN_GMHI",
  margin.params = list(fill = "Phenotype", color = "black", size = 0.1)
)

fig4d <- ggscatterhist(
  final, x = "R_alpha", y = "G_alpha",
  color = "Phenotype",
  palette = c("steelblue", "orange2"),
  margin.plot = "histogram",
  ggtheme = theme_bw(),
  show.legend = FALSE,xlab = "RAW_Shannon",ylab = "GAN_Shannon",
  margin.params = list(fill = "Phenotype", color = "black", size = 0.1)
)

pdf("图3-7.pdf",width=8,height=8)
fig4a
fig4b
fig4c
fig4d
dev.off()

pdf("图3-3.pdf")
fig3a <- ggplot(final, aes(x=Phenotype, y=G_GMHI, fill=Phenotype)) + 
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")+
  theme_classic()           
fig3a + rremove("legend") + theme(axis.text=element_text(size=14,face="bold"),
                                  axis.title=element_text(size=14,face="bold"))+
  scale_colour_manual(values=c("Healthy"="steelblue","Nonhealthy"="orange2"))+
  stat_compare_means(label = "p.format",method = "wilcox.test",label.x.npc = "middle")+
  scale_fill_manual(values=c('steelblue','orange2'))+
  labs(x = "",y="GAN_GMHI",title = "Fig2a")+
  geom_hline(aes(yintercept=0),color="grey",linetype='dashed')+
  coord_cartesian(ylim = c(-6,6))
cliffDelta(data = final,G_GMHI~Phenotype) #0.917

fig3b <- ggplot(final, aes(x=Phenotype, y=G_alpha, fill=Phenotype)) + 
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")+
  theme_classic()           
fig3b + rremove("legend") + theme(axis.text=element_text(size=14,face="bold"),
                                  axis.title=element_text(size=14,face="bold"))+
  scale_colour_manual(values=c("Healthy"="steelblue","Nonhealthy"="orange2"))+
  stat_compare_means(label = "p.format",method = "wilcox.test",label.x.npc = "middle")+
  scale_fill_manual(values=c('steelblue','orange2'))+
  labs(x = "",y="Shannon Diversity(GAN)",title = "Fig2b")+
  coord_cartesian(ylim = c(0,5))
cliffDelta(data = final,G_alpha~Phenotype) #-0.0506

fig3c <- ggplot(final, aes(x=Phenotype, y=R_GMHI, fill=Phenotype)) + 
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")+
  theme_classic()           
fig3c + rremove("legend") + theme(axis.text=element_text(size=14,face="bold"),
                                  axis.title=element_text(size=14,face="bold"))+
  scale_colour_manual(values=c("Healthy"="steelblue","Nonhealthy"="orange2"))+
  stat_compare_means(label = "p.format",method = "wilcox.test",label.x.npc = "middle")+
  scale_fill_manual(values=c('steelblue','orange2'))+
  labs(x = "",y="RAW_GMHI",title = "Fig2c")+
  geom_hline(aes(yintercept=0),color="grey",linetype='dashed')+
  coord_cartesian(ylim = c(-6,6))
cliffDelta(data = final,R_GMHI~Phenotype) #0.557

fig3d <- ggplot(final, aes(x=Phenotype, y=R_alpha, fill=Phenotype)) + 
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")+
  theme_classic()           
fig3d + rremove("legend") + theme(axis.text=element_text(size=14,face="bold"),
                                  axis.title=element_text(size=14,face="bold"))+
  scale_colour_manual(values=c("Healthy"="steelblue","Nonhealthy"="orange2"))+
  stat_compare_means(label = "p.format",method = "wilcox.test",label.x.npc = "middle")+
  scale_fill_manual(values=c('steelblue','orange2'))+
  labs(x = "",y="Shannon Diversity(RAW)",title = "Fig2d")+
  coord_cartesian(ylim = c(0,5))
cliffDelta(data = final,R_alpha~Phenotype) #0.103
dev.off()


#1+12
c<-c(distinctColorPalette(13))

pdf("图3-4.pdf",width=8,height=6)
fig3e<-ggplot(final, aes(x=PA, y=G_GMHI, fill=Phenotype_all)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1,fill="white")
fig3e+scale_fill_manual(values=c)+
  geom_hline(aes(yintercept=0),color="grey",linetype='dashed')+
  theme_classic()+
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5))+
  labs(x = "",y="GAN_GMHI",title = "Fig2e")+rremove("legend")+
  coord_flip(ylim = c(-6,6))

fig3f<-ggplot(final, aes(x=PA, y=R_GMHI, fill=Phenotype_all)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1,fill="white")
fig3f+scale_fill_manual(values=c)+
  geom_hline(aes(yintercept=0),color="grey",linetype='dashed')+
  theme_classic()+
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5))+
  labs(x = "",y="RAW_GMHI",title = "Fig2f")+rremove("legend")+
  coord_flip(ylim = c(-6,6))

fig3g<-ggplot(final, aes(x=PA, y=G_alpha, fill=Phenotype_all)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1,fill="white")
fig3g+scale_fill_manual(values=c)+
  theme_classic()+
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5))+
  labs(x = "",y="GAN_Shannon",title = "Fig2g")+rremove("legend")+
  coord_flip(ylim = c(0,5))

fig3h<-ggplot(final, aes(x=PA, y=R_alpha, fill=Phenotype_all)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1,fill="white")
fig3h+scale_fill_manual(values=c)+
  theme_classic()+
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5))+
  labs(x = "",y="RAW_Shannon",title = "Fig2h")+rremove("legend")+
  coord_flip(ylim = c(0,5))
dev.off()


library(tidyverse)
library( tidytidbits )
final_healthy <- final%>%filter(.$PA=='Healthy')
final_healthy$Batch = as.character(final_healthy$Batch)
final_healthy$Batch[which(final_healthy$Batch%in%c('Sankaranarayanan (2015)','Tanca (2017)','Hall (2017)','Karlsson (2012)'))] <-'Other'
final_healthy$Batch[which(final_healthy$Batch=='HMP1: Huttenhower (2012) and HMP2: Lloyd-Price (2017)')] <-'HMP1&2'
final_healthy$Phenotype_Batch <- paste(final_healthy$Batch,final_healthy$PA,sep = ' ')

final_CD <- final%>%filter(.$PA=='CD')
final_CD$Batch = as.character(final_CD$Batch)
final_CD$Batch[which(final_CD$Batch%in%c('Nielsen (2014)','Hall (2017)'))] <-'Other'
final_CD$Phenotype_Batch <- paste(final_CD$Batch,final_CD$PA,sep = ' ')

final_CRC <- final%>%filter(.$PA=='CRC')
final_CRC$Batch = as.character(final_CRC$Batch)
final_CRC$Phenotype_Batch <- paste(final_CRC$Batch,final_CRC$PA,sep = ' ')

final_T2D <- final%>%filter(.$PA=='T2D')
final_T2D$Batch = as.character(final_T2D$Batch)
final_T2D$Phenotype_Batch <- paste(final_T2D$Batch,final_T2D$PA,sep = ' ')

final_UC <- final%>%filter(.$PA=='UC',.$Batch!='Hall (2017)')
final_UC$Batch = as.character(final_UC$Batch)
final_UC$Phenotype_Batch <- paste(final_UC$Batch,final_UC$PA,sep = ' ')

final_batch = rbind.data.frame(final_CRC,final_CD,final_healthy,final_T2D,final_UC)
final_batch$Phenotype_Batch<-factor(final_batch$Phenotype_Batch,levels=unique(final_batch$Phenotype_Batch))

library(randomcoloR)
c1= c(distinctColorPalette(39))

pdf('图3-8.pdf',width=9,height=12)
fig5a<-ggplot(final_batch, aes(x=Phenotype_Batch, y=G_GMHI, fill=Phenotype_Batch)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1,fill="white")
fig5a+scale_fill_manual(values=c1)+
  geom_hline(aes(yintercept=0),color="grey",linetype='dashed')+
  theme_bw()+
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5))+
  labs(x = "",y="GAN_GMHI",title = "Fig5a")+rremove("legend")+
  coord_flip(ylim = c(-6,6))

fig5b<-ggplot(final_batch, aes(x=Phenotype_Batch, y=R_GMHI, fill=Phenotype_Batch)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1,fill="white")
fig5b+scale_fill_manual(values=c1)+
  geom_hline(aes(yintercept=0),color="grey",linetype='dashed')+
  theme_bw()+
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5))+
  labs(x = "",y="RAW_GMHI",title = "Fig5b")+rremove("legend")+
  coord_flip(ylim = c(-6,6))

fig5c<-ggplot(final_batch, aes(x=Phenotype_Batch, y=G_alpha, fill=Phenotype_Batch)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1,fill="white")
fig5c+scale_fill_manual(values=c1)+
  theme_bw()+
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5))+
  labs(x = "",y="GAN_Shannon",title = "Fig5c")+rremove("legend")+
  coord_flip(ylim = c(0,5))

fig5d<-ggplot(final_batch, aes(x=Phenotype_Batch, y=R_alpha, fill=Phenotype_Batch)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1,fill="white")
fig5d+scale_fill_manual(values=c1)+
  theme_bw()+
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5))+
  labs(x = "",y="RAW_Shannon",title = "Fig5c")+rremove("legend")+
  coord_flip(ylim = c(0,5))
dev.off()

test_final <- read.csv('test_GMHI2.csv')
test_final$Phenotype_Batch <- paste(test_final$Batch,test_final$PA,sep = ' ')

pdf('图3-9.pdf',width=6)
fig6a <-ggplot(test_final, aes(x=Phenotype, y=R_GMHI, fill=Phenotype)) + 
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")+
  theme_bw()       
fig6a + rremove("legend") + theme(axis.text=element_text(size=14,face="bold"),
                                  axis.title=element_text(size=14,face="bold"))+
  scale_colour_manual(values=c("Healthy"="steelblue","Nonhealthy"="orange2"))+
  stat_compare_means(label = "p.format",method = "wilcox.test",label.x.npc = "middle")+
  scale_fill_manual(values=c('steelblue','orange2'))+
  labs(x = "",y="RAW_GMHI",title = "Fig6a")+
  geom_hline(aes(yintercept=0),color="grey",linetype='dashed')+
  coord_cartesian(ylim = c(-6,6))
cliffDelta(data = test_final,R_GMHI~Phenotype) #0.643

fig6b <-ggplot(test_final, aes(x=Phenotype, y=G_GMHI1, fill=Phenotype)) + 
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")+
  theme_bw()        
fig6b + rremove("legend") + theme(axis.text=element_text(size=14,face="bold"),
                                  axis.title=element_text(size=14,face="bold"))+
  scale_colour_manual(values=c("Healthy"="steelblue","Nonhealthy"="orange2"))+
  stat_compare_means(label = "p.format",method = "wilcox.test",label.x.npc = "middle")+
  scale_fill_manual(values=c('steelblue','orange2'))+
  labs(x = "",y="GAN_GMHI",title = "Fig6b")+
  geom_hline(aes(yintercept=0),color="grey",linetype='dashed')+
  coord_cartesian(ylim = c(-6,6))
cliffDelta(data = test_final,GMHI_new~Phenotype) #0.443
dev.off()


test_final$Phenotype_Batch<-factor(test_final$Phenotype_Batch,levels = c("Wirbel (2019) CRC","Thomas (2019) CRC","Dhakan (2019) Healthy",
                                              "Wirbel (2019) Healthy","Bedarf (2017) Healthy",
                                              "Thomas (2019) CA","Vaughn (2016) CD","Gupta (2020) RA",
                                              "Loomba (2017) NAFLD","Wen (2017) AS","Qin (2014) LC"))

pdf('图3-10.pdf',width=8)
c11<-c(distinctColorPalette(11))
fig6c<-ggplot(test_final, aes(x=Phenotype_Batch, y=G_GMHI1, fill=Phenotype_Batch)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1,fill="white")
fig6c+scale_fill_manual(values=c11)+
  geom_hline(aes(yintercept=0),color="grey",linetype='dashed')+
  theme_bw()+
  theme(axis.text=element_text(size=11,face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "",y="GAN_GMHI",title = "Fig6c")+rremove("legend")+
  coord_cartesian(ylim = c(-6,6))

fig6d<-ggplot(test_final, aes(x=Phenotype_Batch, y=R_GMHI, fill=Phenotype_Batch)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1,fill="white")
fig6d+scale_fill_manual(values=c11)+
  geom_hline(aes(yintercept=0),color="grey",linetype='dashed')+
  theme_bw()+
  theme(axis.text=element_text(size=11,face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "",y="RAW_GMHI",title = "Fig6d")+rremove("legend")+
  coord_cartesian(ylim = c(-6,6))

dev.off()