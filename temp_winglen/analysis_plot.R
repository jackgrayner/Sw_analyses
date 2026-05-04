library(ggplot2)
library(lme4)
library(car)
library(MuMIn)
library(ggbeeswarm)
library(gridExtra)
library(lmerTest)
library(dplyr)

mes<-read.csv("mes_nums.csv")

#what affects pronotum length?
lme.rpl<-lmer(scale((pronotum)) ~ sex + factor(temp) + sw.yn +
                scale(rearing_density) + (1|f0_box/rep),data=mes)
Anova(lme.rpl)
summary(lme.rpl)
qqnorm(resid(lme.rpl))
qqline(resid(lme.rpl))

#does sw frequency differ between temps? (note small sample of 15 Sw)
fisher.test(table(mes$sw.yn,mes$temp)) #no

#what affects wing length? 
lme.rwl<-lmer(scale(log10(rwing_length)) ~ sex + scale(log10(pronotum)) + factor(temp) + 
                scale(rearing_density) + (1|f0_box/rep),data=mes[mes$sw.yn=="Lw",])
#lme.rwl<-lmer(scale(log10(rwing_length)) ~ sex + scale(log10(pronotum)) + factor(temp) + sw.yn +
#                scale(rearing_density) + (1|f0_box/rep),data=mes)
Anova(lme.rwl)
summary(lme.rwl)
qqnorm(resid(lme.rwl))
qqline(resid(lme.rwl))

g.len<-ggplot(mes,aes(colour=temp,group=sw.yn,x=log10(pronotum),y=log10(as.numeric(rwing_length))))+
  geom_point(data=mes[mes$sw.yn=="Lw",],alpha=0.75,size=0.75)+
  geom_point(data=mes[mes$sw.yn=="Sw",],shape='triangle',alpha=1)+
  ylab('Log10 wing length (mm)')+xlab('Log10 pronotum length (mm)')+
  theme_minimal()+scale_colour_manual(values=c('#4789b5', '#e6b327', '#c73838'))+
  scale_fill_manual(values=c('#4789b5', '#e6b327', '#c73838'))+
  geom_smooth(data=mes[!mes$sw.yn=="Sw",],
              aes(colour=temp,fill=temp,group=paste0(temp,sex),x=log10(pronotum),y=log10(as.numeric(rwing_length))),
              method='lm')+
  facet_grid(.~sex,scales='free')+theme(panel.border=element_rect(colour='black'))
g.len

ggsave('wing-pronotum.png',plot=g.len,dpi=600,height=5,width=7)
