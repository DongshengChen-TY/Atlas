#####################################
########calculate virus for TR
#####################################
library(Seurat)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)

dir.create("VirusScan")
setwd("VirusScan")

virus<-read.csv("viralReceptor.csv")

virus<-data.frame(Family=virus$Family,Virus=virus$Virus,Gene=virus$Gene)
gene<-virus$Gene
gene<-capitalize(tolower(gene))
virus$Gene<-gene

TR$Tissue <- TR$Tissue

TR$fi<-paste(TR$Tissue,TR$Celltype,sep=".")
p6<-DotPlot(TR,assay="RNA",features=unique(virus$Gene),cols=c("blue","red"),group.by="fi")
write.table(p6$data,"dotplotdata_virus_TR.txt",sep="\t",quote=F)

e<-data.frame("avg.exp"=p6$data$avg.exp,"pct.exp"=p6$data$pct.exp,"Gene"=p6$data$features.plot,"tissue"=p6$data$id)
virus$Gene <- virus$Gene
data<-merge(virus,e,by="Gene")

data<-data%>%separate(.,col=tissue,into=c("Tissue","Celltype"),sep="\\.",remove=F)
write.table(data,"virus_celltype_dotplot_TR.txt",sep="\t",quote=F,row.names=F)

data$Celltype<-as.character(data$Celltype)
data<-data[which(data$Celltype!="Unknown" & data$Celltype!="Others"),]

labels<-c("Lung"="Lung","Kidney"="Kidney","Liver"="Liver","Spleen"="Spleen","Heart"="Heart")
data$Tissue<-factor(data$Tissue,c("Lung","Kidney","Liver","Spleen","Heart"))

p7 <- ggplot(data, mapping=aes_string(x="Celltype",y="Gene",fill="pct.exp")) +
  geom_tile()+
  scale_fill_gradient(low="white",high = "red")+
  scale_y_discrete(position = "right")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) +
  theme_classic()+
  facet_grid(Family+Virus~Tissue,shrink = T,scales = "free",space = "free",switch = "y",as.table = T,labeller = labeller(Tissue=labels))+
  theme(strip.text.x = element_text(angle=0,face="bold",size = 16),strip.text.y=element_text(angle=180,face="bold",size = 16))+
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust = 1,face = "bold",color = "black",size = 16),axis.text.y.right = element_text(face = "bold.italic",color="black",size=16))

pdf("Virus_heatmap_TR.pdf",width=50,height=80)
print(p7)
dev.off()


data$VirusTissueCelltype <- paste(data$Virus,data$tissue,sep = ".")
EntryTop1 <- data %>% group_by(VirusTissueCelltype) %>% top_n(n = 1,wt = pct.exp)

SelecetVirus <- c("African swine fever virus","Classical swine fever virus","Influenza A virus","Japanese encephalitis virus",
                  "Ebola virus","Pseudorabies virus","Porcine reproductive and respiratory syndrome virus","Rabies lyssavirus",
                  "Porcine circovirus","Transmissible gastroenteritis virus","Mammalian orthoreovirus")

Number <- c()

EntryTop1_Select <- EntryTop1[which(EntryTop1$Virus %in% SelecetVirus),]
write.table(EntryTop1_Select, "TR_Selected11Virus_data.tsv",sep = "\t",quote = F)

p7 <- ggplot(EntryTop1_Select, mapping=aes_string(x="Celltype",y="Gene",fill="pct.exp")) +
  geom_tile()+
  scale_fill_gradient(low="white",high = "red")+
  scale_y_discrete(position = "right")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) +
  theme_classic()+
  facet_grid(Family+Virus~Tissue,shrink = T,scales = "free",space = "free",switch = "y",as.table = T,labeller = labeller(Tissue=labels))+
  theme(strip.text.x = element_text(angle=0,face="bold",size = 16),strip.text.y=element_text(angle=180,face="bold",size = 16))+
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust = 1,face = "bold",color = "black",size = 16),
        axis.text.y.right = element_text(face = "bold.italic",color="black",size=16))

pdf("Virus_heatmap_TR_Selected.pdf",width=50,height=20)
print(p7)
dev.off()

TRScan1<-EntryTop1_Select[,c(3,4,8,9,6)]
library(reshape2)
TRScanwide<-dcast(TRScan1,Virus~Tissue+Celltype,fun.aggregate = max)
rownames(TRScanwide) <- TRScanwide$Virus

rownames(TRScanwide)<-TRScanwide[,1]
TRScanwide<-TRScanwide[,-1]
TRScan<-TRScanwide
for (i in 1:ncol(TRScan)) {
  TRScan[,i]<-as.character(TRScan[,i])
  TRScan[,i][which(TRScanwide[,i]<1)]<-"Hard"
  TRScan[,i][which(10>TRScanwide[,i] & TRScanwide[,i]>=1)]<-"Mild"
  TRScan[,i][which(TRScanwide[,i]>=10)]<-"Easy"
  
}
TRScan <- TRScan[SelecetVirus,]
write.csv(TRScan,"TR_virus.csv",row.names = T)



TRScan2 <- TRScan1
TRScan2$signal <- 2
for (i in row.names(TRScan2)) {
  TRScan2[i,"signal"][which(TRScan1[i,"pct.exp"]<1)]<-0
  TRScan2[i,"signal"][which(10>TRScan1[i,"pct.exp"] & TRScan1[i,"pct.exp"]>=1)]<-1
  TRScan2[i,"signal"][which(TRScan1[i,"pct.exp"]>=10)]<-2
}

TRScan2 <- TRScan1
TRScan2$signal <- "green"
for (i in row.names(TRScan2)) {
  TRScan2[i,"signal"][which(TRScan1[i,"pct.exp"]<1)]<-"red"
  TRScan2[i,"signal"][which(10>TRScan1[i,"pct.exp"] & TRScan1[i,"pct.exp"]>=1)]<-"yellow"
  TRScan2[i,"signal"][which(TRScan1[i,"pct.exp"]>=10)]<-"green"
}
TRScan2$signal <- factor(TRScan2$signal, levels=c("red","yellow","green"))

p7 <- ggplot(TRScan2, mapping=aes_string(x="Celltype",y="Virus",fill="signal")) +
  geom_tile()+
  scale_y_discrete(position = "right")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) +
  theme_classic()+
  facet_grid(Virus~Tissue,shrink = T,scales = "free",space = "free",switch = "y",
             as.table = T,labeller = labeller(Tissue=labels))+
  scale_fill_manual(values=c("red","yellow","green"))+
  theme(strip.text.x = element_text(angle=0,face="bold",size = 16),strip.text.y=element_text(angle=180,face="bold",size = 16))+
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust = 1,face = "bold",color = "black",size = 16),
        axis.text.y.right = element_text(face = "bold.italic",color="black",size=16)) 

pdf("Virus_heatmap_TR_Selected_Signal.pdf",width=40,height=10)
print(p7)
dev.off()


