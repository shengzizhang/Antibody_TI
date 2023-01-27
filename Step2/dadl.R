#!/usr/bin/Rscript
library(ggplot2);

dat<-read.table('dadl.txt',sep=' ',nrow=length(readLines('dadl.txt'))-1);

p<-ggplot(dat,aes(x=V1,y=V2))+geom_point()+geom_line()+geom_errorbar(aes(ymin=V2-V4*2,ymax=V2+V4*2),width=0.04)+xlab('lambda')+ylab('dV/dL')+theme_classic()

pdf('dadl_overview.pdf',width=4,height=3)
print(p)

dev.off()
t<-getwd()
if(grepl('protein',t)){
dat$g='protein'
dat1<-read.table('../complex/dadl.txt',sep=' ',nrow=length(readLines('../complex/dadl.txt'))-1);
dat1$g='complex'
dat2<-rbind(dat,dat1)
p<-ggplot(dat2,aes(x=V1,y=V2, group=g, color=g))+geom_point()+geom_line()+geom_errorbar(aes(ymin=V2-V4,ymax=V2+V4),width=0.04)+xlab('lambda')+ylab('dV/dL')+theme_classic()
pdf('dadl_compare.pdf',width=4,height=3)
print(p)

dev.off()
}
dat<-read.table('combined_energy.txt',sep='\t');

p<-ggplot(dat,aes(x=V3,y=V4))+geom_line()+geom_smooth(method='gam',se=T)+xlab('ps')+ylab('Total energy')+theme_classic()

pdf('Total_energy.pdf',width=5,height=3)
print(p)

dev.off()

dat<-read.table('combined_dvdl.txt',sep='\t');
dat$cumsum<-ave(dat$V4,dat$V2,FUN=cumsum)
dat$cumaverage<-dat$cumsum/dat$V3
p<-ggplot(dat,aes(x=V3,y=V4))+geom_line()+geom_smooth(method='gam',se=T)+facet_wrap(~V2,ncol=4,scales='free')+xlab('ps')+ylab('Dvdl')+theme_classic()

pdf('Dvdl_lambda.pdf',width=9,height=6)
print(p)

dev.off()

p<-ggplot(dat,aes(x=V3,y=cumaverage))+geom_line()+facet_wrap(~V2,ncol=4,scales='free')+xlab('ps')+ylab('Dvdl')+theme_classic()

