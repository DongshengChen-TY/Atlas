###########Plot
d<-rbind(a$data,b$data,c$data)
d$id<-as.character(d$id)
d<-d%>%separate(.,col="id",into=c("NewCelltype","Tissue"),sep="\\.",remove=F)
d<-d[which(d$NewCelltype!="Unknown"),]
##指定绘图顺序
d$features.plot<-factor(d$features.plot,levels = c("TMPRSS4","TMPRSS2","ACE2"))
d$Tissue<-factor(d$Tissue,levels = c("Lung","Kidney","Liver","Spleen","Heart"))
d$NewCelltype<-as.factor(d$NewCelltype)
d$NewCelltype <- factor(d$NewCelltype,levels(d$NewCelltype)[c(length(levels(d$NewCelltype)):1)])

ratio<-read.csv("ACE2_TMPRSS2_ratio_fi.txt",sep="\t")
ratio1<-read.csv("ACE2_TMPRSS4_ratio_fi.txt",sep="\t")
ratio<-rbind(ratio,ratio1)
co<-ratio[which(ratio$X1=="ACE2_TMPRSS2" | ratio$X1=="ACE2_TMPRSS4"),]
co<-data.frame("features.plot"=co$X1,"fi"=co$fi,per=co$X2)
co$per<-as.numeric(co$per)
co$fi<-as.character(co$fi)

co<-co%>%separate(.,col="fi",into=c("NewCelltype","Tissue"),sep="\\.",remove=F)
co<-co[which(co$NewCelltype!="Others"),]
co<-co[which(co$NewCelltype!="Unknown"),]
co<-co[which(co$NewCelltype!="Clara cell"),]
##指定绘图顺序
co$Tissue<-factor(co$Tissue,c("LG","Kidney","Liver","Spleen","Heart"))
co$NewCelltype<-as.factor(co$NewCelltype)
co$NewCelltype <- factor(co$NewCelltype,levels(co$NewCelltype)[c(length(levels(co$NewCelltype)):1)])

write.table(co,"co_exp_ratio_fi_.txt",sep = "\t",quote = F,row.names = F)
labels<-c("Lung"="Lung","Kidney"="Kidney","Liver"="Liver","Spleen"="Spleen","Heart"="Heart")
p3 <- ggplot(data=d,mapping=aes_string(x="features.plot",y="NewCelltype"))+
  geom_point(mapping=aes_string(size="pct.exp",color="avg.exp"),show.legend = TRUE)+
  scale_color_gradient(low="blue",high = "red")+
  scale_size(limits = c(min(d$pct.exp[which(d$pct.exp>0)]),max(d$pct.exp)))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) +
  guides(size = guide_legend(title = "Percent Expressed")) +
  labs(x = "", y = "",color="Average Expression",size="Percent Expressed")+
  theme_bw()+
  facet_grid(Tissue~features.plot,shrink = F,scales = "free",space = "free",switch = "y",labeller = labeller(Tissue=labels))+
  theme(strip.text.x = element_text(angle=0,face="bold"),strip.text.y=element_text(angle=90,face="bold"))+
  theme(axis.text.x = element_text(face = "bold",colour = "black"),axis.text.y = element_text(face="bold",colour = "black"))

p4 <- ggplot(co,aes(x=NewCelltype,y=per,fill=features.plot))+
  geom_bar(stat="identity",width = 0.5,position = position_dodge())+
  scale_color_brewer(palette = "Set3")+
  coord_flip()+
  theme_bw()+
  labs(x = "", y = "",title="Co-expression")+
  facet_grid(Tissue~.,shrink = F,scales = "free",space = "free",switch = "y",labeller = labeller(Tissue=labels))+
  theme(strip.text.x = element_text(angle=0,face="bold"),strip.text.y=element_text(angle=90,face="bold"))+
  theme(axis.text.x = element_text(face = "bold",colour = "black"),axis.text.y = element_text(face="bold",colour = "black"))

###########Save pdf
pdf("_plot_A_T_NewCelltype.pdf",width=13,height=27)
plot_grid(p4,p3)
dev.off()

p3 <- ggplot(data=d,mapping=aes_string(x="features.plot",y="NewCelltype"))+
  geom_point(mapping=aes_string(size="pct.exp",color="avg.exp"),show.legend = TRUE)+
  scale_color_gradient(low="blue",high = "red")+
  scale_size(limits = c(min(d$pct.exp[which(d$pct.exp>0)]),max(d$pct.exp)))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) +
  guides(size = guide_legend(title = "Percent Expressed")) +
  labs(x = "", y = "",color="Average Expression",size="Percent Expressed")+
  theme_bw()+
  facet_grid(Tissue~features.plot,shrink = F,scales = "free",space = "free",switch = "y",labeller = labeller(Tissue=labels))+
  theme(strip.text.x = element_text(angle=0,face="bold"),strip.text.y=element_text(angle=90,face="bold"))+
  theme(axis.text.x = element_text(face = "bold",colour = "black"),axis.text.y = element_text(face="bold",colour = "black"))+NoLegend()

p4 <- ggplot(co,aes(x=NewCelltype,y=per,fill=features.plot))+
  geom_bar(stat="identity",width = 0.5,position = position_dodge())+
  scale_color_brewer(palette = "Set3")+
  coord_flip()+
  theme_bw()+
  labs(x = "", y = "",title="Co-expression")+
  facet_grid(Tissue~.,shrink = F,scales = "free",space = "free",switch = "y",labeller = labeller(Tissue=labels))+
  theme(strip.text.x = element_text(angle=0,face="bold"),strip.text.y=element_text(angle=90,face="bold"))+
  theme(axis.text.x = element_text(face = "bold",colour = "black"),axis.text.y = element_text(face="bold",colour = "black"))+NoLegend()

###########Save pdf
pdf("_plot_A_T_NewCelltype_Nolegend.pdf",width=10,height=20)
plot_grid(p4,p3)
dev.off()



