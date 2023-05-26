
library(phytools)
library(ape)
library(ggplot2)
library(reshape2)
library(phangorn)
library(cowplot)
library(ggridges)
library(ggnewscale)
library(ggExtra)
library(ggtern)
theme_set(theme_cowplot())


fd = read.table("OOA1_OOA2_West_localFstats__20_20.txt",header=TRUE)
dxy2 = read.table("true_pops_dxy.txt",header=TRUE)
fst2 = read.table("true_pops_fst.txt",header=TRUE)
chrs = c("chr2L","chr2R","chr3L","chr3R","chrX")
fd$recomb = NA
fd$dxy_west_ooa2 = NA
fd$dxy_west_cen = NA
fd$dxy_cen_ooa2 = NA

fd$fst_west_ooa2 = NA
fd$fst_west_cen = NA
fd$fst_cen_ooa2 = NA

for(c in chrs){
  chrtable = read.table(paste("genetic_map_comeron2012_dm6_",c,".txt",sep=""),header=TRUE)
  idc = which(fd$chr==c)
  wndid = sapply(fd[idc,]$windowStart,function(x) which(chrtable$Position.bp>x)[1])
  fd$recomb[idc] = chrtable$Rate.cM.Mb.[wndid]
  
  #And then annotate fst/dxy:
  
  #fst - West OOA2
  wndid = unlist(sapply(fd[idc,]$windowStart,function(x) which(fst2$window_pos_1 == x & fst2$pop1 %in% c("West","OOA2") & fst2$pop2 %in% c("West","OOA2") & fst2$chromosome==c)))
  fd$fst_west_ooa2[idc] = fst2$avg_wc_fst[wndid]
  #fst - West Cen
  wndid = unlist(sapply(fd[idc,]$windowStart,function(x) which(fst2$window_pos_1 == x & fst2$pop1 %in% c("West","CEN") & fst2$pop2 %in% c("West","CEN") & fst2$chromosome==c)))
  fd$fst_west_cen[idc] = fst2$avg_wc_fst[wndid]
  #fst - Cen OOA2
  wndid = unlist(sapply(fd[idc,]$windowStart,function(x) which(fst2$window_pos_1 == x & fst2$pop1 %in% c("CEN","OOA2") & fst2$pop2 %in% c("CEN","OOA2") & fst2$chromosome==c)))
  fd$fst_cen_ooa2[idc] = fst2$avg_wc_fst[wndid]
  #dxy - West OOA2
  wndid = unlist(sapply(fd[idc,]$windowStart,function(x) which(dxy2$window_pos_1 == x & dxy2$pop1 %in% c("West","OOA2") & dxy2$pop2 %in% c("West","OOA2") & dxy2$chromosome==c)))
  fd$dxy_west_ooa2[idc] = dxy2$avg_dxy[wndid]
  #dxy - West Cen
  wndid = unlist(sapply(fd[idc,]$windowStart,function(x) which(dxy2$window_pos_1 == x & dxy2$pop1 %in% c("West","CEN") & dxy2$pop2 %in% c("West","CEN") & dxy2$chromosome==c)))
  fd$dxy_west_cen[idc] = dxy2$avg_dxy[wndid]
  #dxy - Cen OOA2
  wndid = unlist(sapply(fd[idc,]$windowStart,function(x) which(dxy2$window_pos_1 == x & dxy2$pop1 %in% c("CEN","OOA2") & dxy2$pop2 %in% c("CEN","OOA2") & dxy2$chromosome==c)))
  fd$dxy_cen_ooa2[idc] = dxy2$avg_dxy[wndid]
}


fd$west_branch = (fd$dxy_west_cen+fd$dxy_west_ooa2-fd$dxy_cen_ooa2)/(2*(fd$dxy_west_cen+fd$dxy_west_ooa2+fd$dxy_cen_ooa2))
fd$ooa_branch = (-1*fd$dxy_west_cen+fd$dxy_west_ooa2+fd$dxy_cen_ooa2)/(2*(fd$dxy_west_cen+fd$dxy_west_ooa2+fd$dxy_cen_ooa2))
fd$cen_branch = (fd$dxy_west_cen-fd$dxy_west_ooa2+fd$dxy_cen_ooa2)/(2*(fd$dxy_west_cen+fd$dxy_west_ooa2+fd$dxy_cen_ooa2))

fd$fst_cen_ooa2[which(fd$fst_cen_ooa2<0)] = 0
fd$fst_west_ooa2[which(fd$fst_west_ooa2<0)] = 0
fd$fst_west_cen[which(fd$fst_west_cen<0)] = 0

fd$t_west_cen = -log(1-fd$fst_west_cen)
fd$t_west_ooa = -log(1-fd$fst_west_ooa2)
fd$t_ooa_cen = -log(1-fd$fst_cen_ooa2) 


fd$pbs_west = (fd$t_west_cen+fd$t_west_ooa-fd$t_ooa_cen)/2
fd$pbs_ooa = (fd$t_west_ooa+fd$t_ooa_cen-fd$t_west_cen)/2
fd$pbs_cen = (fd$t_west_cen-fd$t_west_ooa+fd$t_ooa_cen)/2

fd$pbs_west2 = sapply(fd$pbs_west,function(x) max(0,x,na.rm=TRUE))
fd$pbs_ooa2 = sapply(fd$pbs_ooa,function(x) max(0,x,na.rm=TRUE))
fd$pbs_cen2 = sapply(fd$pbs_cen,function(x) max(0,x,na.rm=TRUE))


fd$pbs_west3 = sapply(1:46320,function(x) fd$pbs_west2[x]/(fd$pbs_cen2[x]+fd$pbs_ooa2[x]+fd$pbs_west2[x]))
fd$pbs_ooa3 = sapply(1:46320,function(x) fd$pbs_ooa2[x]/(fd$pbs_cen2[x]+fd$pbs_ooa2[x]+fd$pbs_west2[x]))
fd$pbs_cen3 = sapply(1:46320,function(x) fd$pbs_cen2[x]/(fd$pbs_cen2[x]+fd$pbs_ooa2[x]+fd$pbs_west2[x]))



fd2 = subset(fd, f_d>=0 & f_d <=1) 

p1 = ggplot(fd2,aes(x=recomb,y=f_d))+scale_fill_gradient(low = "#3FA9F5",high = "#B72467")+
  geom_hex(bins=10)+geom_smooth(method="lm",fill=NA,col="white",size=2)+geom_smooth(method="lm",fill=NA,col="black")+
  labs(y=expression(italic(f)[D]),x="Recombination (cM/Mb)",col="")

p2 = ggplot(fd2,aes(x=pbs_ooa3,y=f_d))+scale_fill_gradient(low = "#3FA9F5",high = "#B72467")+
  geom_hex(bins=10)+geom_smooth(method="lm",fill=NA,col="white",size=2)+geom_smooth(method="lm",fill=NA,col="black")+
  labs(x="OOA2 Branch",y=expression(italic(f)[D]),col="")+coord_cartesian(ylim=c(0,1))

legend_1 = get_legend(p2)

p3 = ggplot(fd2,aes(x=pbs_west3,y=f_d))+scale_fill_gradient(low = "#3FA9F5",high = "#B72467")+
  geom_hex(bins=10)+geom_smooth(method="lm",fill=NA,col="white",size=2)+geom_smooth(method="lm",fill=NA,col="black")+
  labs(x="West Branch",y=expression(italic(f)[D]),col="")+coord_cartesian(ylim=c(0,1))


plot_grid(p1,plot_grid(p2+theme(legend.position="none"),p3+theme(legend.position="none"),legend_1,nrow=1,rel_widths = c(0.9,0.9,0.2)),nrow=2)

plot_grid(p1,p2,p3,nrow=1,rel_widths = c(1,1,1))



library(ggtern)
library(colorspace)

fd2 = subset(fd,abs(f_d)<1 & f_d>=0 & pbs_ooa2>0)
fd3 = subset(fd2,f_d>0.5)

all=ggplotGrob(ggplot(fd2,aes(y=pbs_cen3,z=pbs_ooa3,x=pbs_west3))+
                 stat_density_tern(bdl = 0.01,geom="polygon",n=100,aes(fill= ..level..,alpha=..level..))+
                 theme_bw()+
                 coord_tern()+
                 scale_fill_gradient(low = "#3FA9F5",high = "#B72467")+
                 labs(x="West",y="South1",z="OOA2",title="Genome-wide windows")+theme(legend.position="none"))

high=ggplotGrob(ggplot(fd3,aes(y=pbs_cen3,z=pbs_ooa3,x=pbs_west3))+
                  theme_bw()+
                  coord_tern()+
                  stat_density_tern(bdl = 0.01,geom="polygon",n=100,aes(fill= ..level..,alpha=..level..))+
                  scale_fill_gradient(low = "#3FA9F5",high = "#B72467")+
                  labs(x="West",y="South1",z="OOA2",title="High introgression windows")+theme(legend.position="none"))

points = ggplotGrob(ggplot(fd2,aes(y=pbs_cen,z=pbs_ooa2,x=pbs_west,col=f_d))+
                      geom_point(size=0.1)+
                      theme_bw()+
                      coord_tern()+
                      scale_color_stepsn(breaks=c(0,0.1,0.25,0.5,0.75,1),colors = c(rgb(0.6,0.6,0.6,0.2,maxColorValue = 1),rgb(0.6,0.6,0.6,0.2,maxColorValue = 1),rgb(0.6,0.6,0.6,0.2,maxColorValue = 1),lighten("#3FA9F5"),"#3FA9F5"))+
                      labs(x="West",y="South1",z="OOA2",col=expression(italic(f)[D])))

row1 = plot_grid(NULL,points,nrow=1,labels=c("","A"))
row2 = plot_grid(all,high,nrow=1,labels=c("B",""))
row3 = plot_grid(p1+theme_cowplot(),p2+theme_cowplot(),nrow=1,labels=c("C","D"))

plot_grid(row1,row2,row3,nrow=3)


