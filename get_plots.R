#script used to generate plots and analyze data

#load in required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(vegan)
library(RColorBrewer)

#function to get nice (I think) colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#function to pull pvalues from linear model objects
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#read data
df<-read.table(file = "meta_data.txt", sep = '\t', stringsAsFactors = F, header = T)

#filter entries without genus names
df<-subset(df,is.na(df$Genus)==F) 
#filter entries without species names 
df<-subset(df,is.na(df$Species)==F)
#filter data Jan 2010 - March 2017
df<-df[df$Date_Count>16,]
df<-df[df$Date_Count<104,]
#handle NAs in house
df[is.na(df)]<-"NA"

#count up each species in each month
count <- aggregate(ID ~ Date_Count + Date + Species + Genus + Family + Order + Class + Phylum + Kingdom,df,length)
colnames(count)[colnames(count)=="ID"] <- "Count"

#get eveness
E <- aggregate(Count ~ Date_Count,count,function(x){diversity(x)/log(length(x))}) 
#get unique species
R <- aggregate(Count ~ Date_Count+Date,count,function(x) length(x)) 

#Figure 1
ggplot(E, aes(Date_Count,Count)) + xlab(NULL) + ylab(label = "Species Eveness") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.minor.x = element_line(color = "white", size=0.5), panel.grid.major.x = element_line(color = "gray", size=1)) + 
  scale_x_continuous(minor_breaks = seq(0:200),breaks = c(17, 29, 41, 53, 65, 77, 89, 101)) +
  scale_y_continuous(sec.axis = sec_axis(~.*2000, name = "Species Richness")) +
  geom_smooth(aes(E$Date_Count,E$Count,fill=gg_color_hue(2)[2]), color=gg_color_hue(2)[1],size=2,show.legend = F) + 
  geom_point(color=gg_color_hue(2)[1],size=3) +
  geom_smooth(aes(R$Date_Count,R$Count/2000,fill=gg_color_hue(2)[1]), color=gg_color_hue(2)[2],size=2,show.legend = F) +
  geom_point(aes(R$Date_Count,R$Count/2000), size=3, color=gg_color_hue(2)[2]) +
  coord_cartesian(xlim = c(17, 103),ylim=c(0,1),expand = F)

#Figure 1 a, absolute one percenters
total_p <- aggregate(Count ~ Species, count, sum)
total_p <- total_p %>%
  mutate(percentile = percent_rank(Count))
list99 <- total_p[total_p$percentile>=0.99,]$Species
list90 <- total_p[total_p$percentile<0.99&total_p$percentile>=0.90,]$Species
abs <- count %>%
  group_by(Date_Count) %>% 
  mutate(percentile = percent_rank(Count))
abs$MO <- ifelse(abs$Species %in% list99,"99",ifelse(abs$Species %in% list90,"90","NMO"))
abs_ag<-aggregate(Count ~ Date_Count + MO, data = abs, sum, na.rm = TRUE)
abs_ag$MO <- factor(abs_ag$MO, levels=c("99","90","NMO"))
pal <- brewer.pal(3,name="YlGnBu")
g1 <- ggplot(abs_ag, aes(x=(Date_Count),y=Count,group=MO,fill=MO))  + geom_area(position='fill') +  scale_fill_manual(values = c(pal[3],pal[2],pal[1])) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(plot.margin = unit(c(1, 1, 4, 1), "lines"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
g1 + coord_cartesian(xlim = c(17, 103),expand = F) + geom_vline(xintercept=17,color="gray") + geom_vline(xintercept=29,color="gray") + geom_vline(xintercept=41,color="gray") + geom_vline(xintercept=53,color="gray") + geom_vline(xintercept=65,color="gray") + geom_vline(xintercept=77,color="gray") + geom_vline(xintercept=89,color="gray") + geom_vline(xintercept=101,color="gray") 

#Figure 1 b, relative one percenters by month
rel <- count %>% 
  group_by(Date_Count) %>% 
  mutate(percentile = percent_rank(Count))
rel$MO <- ifelse(rel$percentile < .90,"NMO",ifelse(rel$percentile < .99,"90","99"))
rel_ag<-aggregate(Count ~ Date_Count + MO, data = rel, sum, na.rm = TRUE)
rel_ag$MO <- factor(rel_ag$MO, levels=c("99","90","NMO"))
g1 <- ggplot(rel_ag, aes(x=(Date_Count),y=Count,group=MO,fill=MO))  + geom_area(position='fill') +  scale_fill_manual(values = c(pal[3],pal[2],pal[1])) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(plot.margin = unit(c(1, 1, 4, 1), "lines"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
g1 + coord_cartesian(xlim = c(17, 103),expand = F) + geom_vline(xintercept=17,color="gray") + geom_vline(xintercept=29,color="gray") + geom_vline(xintercept=41,color="gray") + geom_vline(xintercept=53,color="gray") + geom_vline(xintercept=65,color="gray") + geom_vline(xintercept=77,color="gray") + geom_vline(xintercept=89,color="gray") + geom_vline(xintercept=101,color="gray") 

#Figure 3
M <- count(df,Species,Date_Count)
M <- M %>%
  group_by(Date_Count,Species) %>%
  summarise (N = sum(n)) %>%
  mutate(freq = N / sum(N))
M <- M[M$Species %in% list99,]
M <- complete(M,Species,Date_Count=0:131)
M <- unique(M)
M[is.na(M)]<-0
M <- M[M$Date_Count<104 & M$Date_Count>16,]
NIH_M<-M[M$Species %in% c("Mus musculus","Danio rerio","Saccharomyces cerevisiae","Arabidopsis thaliana","Drosophila melanogaster","Caenorhabditis elegans","Rattus norvegicus"),] 
ggplot(NIH_M,aes(Date_Count, freq, color=Species)) + geom_smooth(se=F,size=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(plot.margin = unit(c(1, 1, 4, 1), "lines"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.minor.x = element_line(color = "white", size=0.5), panel.grid.major.x = element_line(color = "gray", size=1)) + 
  scale_x_continuous(minor_breaks = seq(0:200),breaks = c(17, 29, 41, 53, 65, 77, 89, 101)) 

#get slopes p values and r squared for absolute and relative lines for each species
slopes <- M %>%
  group_by(Species) %>%
  mutate(slope = coef(lm(N ~ Date_Count))[2]) %>%
  mutate(r = summary(lm(N ~ Date_Count))$r.squared) %>%
  mutate(p = lmp(lm(N ~ Date_Count))) %>%
  mutate(freqslope = coef(lm(freq ~ Date_Count))[2]) %>%
  mutate(freqr = summary(lm(freq ~ Date_Count))$r.squared) %>%
  mutate(freqp = lmp(lm(freq ~ Date_Count))) %>%
  mutate(count = sum(N)) %>%
  mutate(freq = mean(freq)) 
speciesums <- data.frame(Species=slopes$Species,Count=slopes$count,Freq=slopes$freq,Slope=slopes$slope,R=slopes$r,P=slopes$p,FreqSlope=slopes$freqslope,RFreq=slopes$freqr,PFreq=slopes$freqr,stringsAsFactors = F)
speciesums <- unique(speciesums)

#actual meta-analysis
library(metafor)
logitP = function(x) 1/(1+exp(-x))
meta <- rel
meta$MO <- ifelse(meta$Species %in% list99,"MO","NMO")
meta <- aggregate(Count ~ Date_Count+MO,meta,sum)
metaA <- meta[meta$MO=="MO",]
metaB <- meta[meta$MO=="NMO",]
meta <- merge(metaA,metaB,by = "Date_Count")
newmeta <- data.frame(meta$Date_Count)
newmeta$MO <- meta$Count.x
newmeta$NMO <- meta$Count.y
meta <- newmeta
meta$ni <- meta$MO + meta$NMO
meta = mutate(meta, ni=ni, prct = 100 * MO/ni)
meta = escalc(measure = "PLO",xi=MO,ni=ni,data=meta)
m0_w = rma(yi, vi,data=meta,knha=T)
logitP( m0_w$b[1] )       #mean
logitP( m0_w$ci.lb[1] )   #CI lower bound
logitP( m0_w$ci.ub[1] )   #CI upper bound
m0_w$tau2                 #heterogeneity
m0_w$QEp                  #p val for heterogeneity Q test
m0_w$I2                   #I^2 heterogeneity statistic


