#script used to generate plots and analyze data for sequencing disparity ms

##NOTES##
#before filtering 55 entries removes because did not have the expected number of rows (entry error, pooled samples)

library(dplyr)
library(vegan)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

#function to get nice (I think) colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#read sra input data
raw_data <- read.table(file = "data_out_filtered.txt", sep = '\t', stringsAsFactors = F, header = F, na.strings = "NA")
names(raw_data) <- c("accession","taxid","strategy","bases","date")
#filter duplicate entries
data_filtered <- subset(raw_data,is.na(raw_data$accession)==F) 
#filter entries without accession (58 removed)
data_filtered <- subset(raw_data,is.na(raw_data$accession)==F) 
#filter entries without date (1599 removed)
data_filtered <- subset(data_filtered,is.na(data_filtered$date)==F) 
#read lineage input data (line 6253 and 8082 have special characters "#" and "'" respectively which were manually removed)
lineage <- read.table(file = "lineage_out.txt", sep = '\t', stringsAsFactors = F, header = T, na.strings = c("NA ","NA"))
#add lineage data to sra data
data_lineage <- merge(data_filtered,lineage)
#filter entries without ncbi genus (14034 removed)
data_lineage_filtered <- subset(data_lineage,is.na(data_lineage$Genus)==F) 
#filter entries without ncbi species (9042 removed)
data <- subset(data_lineage_filtered,is.na(data_lineage_filtered$Species)==F) 
#just keep year and month in date column
data$date <- gsub('-\\d\\d \\d\\d:\\d\\d:\\d\\d','',data$date,perl=T)
#add column datecount which converts months into relative number value (first month = 0, second month = 1, etc)
datecount_list <- data.frame(date=sort(unique(data$date)),datecount=0:124)
data <- merge(data,datecount_list)
#filter to experiments published Jan 2010 - Jul 2018
data <- data[21<=data$datecount & data$datecount<=123,]

#get number of experiments for each species for each month
species_month_total <- data %>%
  group_by(datecount,Species) %>%
  mutate(experiments = length(Species)) %>%
  mutate(bases = sum(bases))  %>%
  select(-taxid, -accession, -strategy) %>%
  distinct()

#get list of top 1%ers
species_percentiles <- species_month_total %>%
  group_by(.dots=c("Species","Genus","Family","Order","Class","Phylum","Kingdom")) %>%
  mutate(experiments = sum(experiments)) %>%
  mutate(bases = sum(bases,na.rm = T)) %>%
  select(-date, -datecount) %>%
  distinct() %>%
  ungroup() %>%
  mutate(percentile = percent_rank(experiments)) %>%
  mutate(percentof = experiments/sum(experiments)) %>%
  mutate(rank = rank(-experiments))

one_percent <- species_percentiles[species_percentiles$percentile>=0.99,]
#what percent of total experiments to 1%ers represent?
sum(one_percent$experiments)/sum(species_percentiles$experiments)
#how many animals in the 1%?
length(one_percent[which(one_percent$Kingdom=="Metazoa"),]$Kingdom)
#plants?
length(one_percent[which(one_percent$Kingdom=="Viridiplantae"),]$Kingdom)
#fungi?
length(one_percent[which(one_percent$Kingdom=="Fungi"),]$Kingdom)
#protists (other?)
length(one_percent[is.na(one_percent$Kingdom),]$Kingdom)
#how many phyla
length(unique(one_percent$Phylum))
#lets look at them
one_percent_phyla <- aggregate(Species ~ Phylum,one_per,length)
#lets look at NIH models (when just genus is listed take the most popular species)
model_list = c("Mus musculus ","Rattus norvegicus ","Saccharomyces cerevisiae ","Schizosaccharomyces pombe ","Neurospora crassa ",
                   "Dictyostelium discoideum ","Caenorhabditis elegans ","Daphnia magna ","Drosophila melanogaster ","Danio rerio ",
                   "Xenopus tropicalis ","Gallus gallus ","Arabidopsis thaliana ")
#table 1
one_percent_models <- one_percent[one_percent$Species %in% model_list,]
#how many models?
sum(one_percent_models$experiments)/sum(species_percentiles$experiments)
#how many just mice?
sum(species_percentiles[species_percentiles$Species=="Mus musculus ",]$experiments)/sum(species_percentiles$experiments)

#get top 1% 10% and rest for each month
species_month_relative <- species_month_total %>% 
  group_by(datecount) %>% 
  mutate(percentile = percent_rank(experiments))
species_month_relative$category <- ifelse(species_month_relative$percentile < .90,"other",ifelse(species_month_relative$percentile < .99,"90","99"))
category_month_relative <- species_month_relative %>%
  group_by(category,datecount) %>%
  summarise(experiments=sum(experiments)) %>%
  ungroup() %>%
  group_by(datecount) %>%
  mutate(freq=experiments/sum(experiments))
category_month_relative$category <- factor(category_month_relative$category, levels=c("99","90","other"))
#linear model of percentage of 99th percentile over time
summary(lm(category_month_relative[category_month_relative$category=="99",]$freq ~ 
     category_month_relative[category_month_relative$category=="99",]$datecount))
#top 1% represented what percent of experiments in 2010
twenty_ten <- species_month_total[species_month_total$datecount >= 21 & species_month_total$datecount <= 32 ,]
twenty_ten <- aggregate(experiments ~ Species,twenty_ten,sum)
twenty_ten$percentile <- percent_rank(twenty_ten$experiments)
sum(twenty_ten[twenty_ten$percentile>=0.99,]$experiments)/sum(twenty_ten$experiments)
#top 1% represented what percent of experiments in 2017?
twenty_seventeen <- species_month_total[species_month_total$datecount >= 105 & species_month_total$datecount <= 116 ,]
twenty_seventeen <- aggregate(experiments ~ Species,twenty_seventeen,sum)
twenty_seventeen$percentile <- percent_rank(twenty_seventeen$experiments)
sum(twenty_seventeen[twenty_seventeen$percentile>=0.99,]$experiments)/sum(twenty_seventeen$experiments)
#figure 1
pal <- brewer.pal(3,name="YlGnBu")
g1 <- ggplot(category_month_relative, aes(x=datecount,y=experiments,group=category,fill=category)) + geom_area(position='fill') +
  scale_fill_manual(values = c(pal[3],pal[2],pal[1])) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), legend.position="none",axis.title.y=element_blank(),axis.text.y=element_blank()) +
  scale_x_continuous(breaks=c(21,33,45,57,69,81,93,105,117))
g1 + coord_cartesian(xlim = c(21, 123),expand = F) 

#calculate species richness and evenness for each month
month_richness <- setNames(aggregate(experiments ~ datecount,species_month_total,function(x) length(x)),c("datecount","richness"))
month_evenness <- setNames(aggregate(experiments ~ datecount,species_month_total,function(x){diversity(x)/log(length(x))}),c("datecount","evenness"))
#linear model of richness over time
lm_richness <- lm(month_richness$richness ~ month_richness$datecount)
summary(lm_richness)
#richness in 2010
length(twenty_ten$Species)
#in 2017
length(twenty_seventeen$Species)
#linear model of evenness over time
summary(lm(month_evenness$evenness ~ month_evenness$datecount))
#figure 2
g2 <- ggplot(month_evenness, aes(datecount,evenness)) +
  #geom_smooth(aes(month_evenness$datecount,month_evenness$evenness),method = 'lm',se=F,size=0.75,color="black") +
  geom_smooth(aes(month_evenness$datecount,month_evenness$evenness,fill=gg_color_hue(2)[2]), color=gg_color_hue(2)[1],size=2,show.legend = F) +
  geom_point(color=gg_color_hue(2)[1],size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(1, 1, 4, 1), "lines"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.minor.x = element_line(color = "white", size=0.5), panel.grid.major.x = element_line(color = "gray", size=1), axis.title.y=element_blank()) + 
  scale_x_continuous(minor_breaks = seq(0:200), breaks = c(21,33,45,57,69,81,93,105,117)) +
  scale_y_continuous(sec.axis = sec_axis(~.*2000)) +
  #geom_smooth(aes(month_richness$datecount,month_richness$richness/2000),method = 'lm',se=F,size=0.75,color="black") +
  geom_smooth(aes(month_richness$datecount,month_richness$richness/2000,fill=gg_color_hue(2)[1]), color=gg_color_hue(2)[2],size=2,show.legend = F) +
  geom_point(aes(month_richness$datecount,month_richness$richness/2000), size=3, color=gg_color_hue(2)[2]) 
g2 + coord_cartesian(xlim = c(21, 123),expand = F) 

#fill data frame for every possible species month combo, 
species_month_full <- unique(complete(species_month_relative[,c(9,2,10)],Species,datecount=21:123))
species_month_full[is.na(species_month_full$experiments),]$experiments = 0
species_month_full <- species_month_full %>%
  group_by(datecount) %>%
  mutate(freq = experiments/sum(experiments))
#get linear slopes and p values of # of experiments and % of dataset for each species over time
species_full <- species_month_full %>%
    group_by(Species) %>%
    mutate(exp_slope = coef(lm(experiments ~ datecount))[2]) %>%
    mutate(exp_p = summary(lm(experiments ~ datecount))$coefficients[2,"Pr(>|t|)"]) %>%
    mutate(freq_slope = coef(lm(freq ~ datecount))[2]) %>%
    mutate(freq_p = summary(lm(freq ~ datecount))$coefficients[2,"Pr(>|t|)"]) %>%
    select(-datecount, -experiments, -freq) %>%
    distinct()
#significant experiment change
sign_exp_full <- species_full[species_full$exp_p<0.05,]
#positive?
sign_pos_exp <- sign_exp_full[which(sign_exp_full$exp_slope>0),]
#negative?
sign_neg_exp <- sign_exp_full[which(sign_exp_full$exp_slope<0),]

#significant frequency change
sign_freq_full <- species_full[species_full$freq_p<0.05,]
#positive?
sign_pos_freq <- sign_freq_full[which(sign_freq_full$freq_slope>0),]
#negative?
sign_neg_freq <- sign_freq_full[which(sign_freq_full$freq_slope<0),]

sign_exp_full[which(sign_exp_full$Species %in% model_list),]
sign_freq_full[which(sign_freq_full$Species %in% model_list),]




species_freq_plotter <- function(x) {
  plot(species_month_full[species_month_full$Species==x,]$datecount,species_month_full[species_month_full$Species==x,]$freq)
  abline(lm(species_month_full[species_month_full$Species==x,]$freq ~ species_month_full[species_month_full$Species==x,]$datecount))
}

#figure 3
NIH_species_month <-species_month_full[species_month_full$Species %in% c("Mus musculus ","Saccharomyces cerevisiae ",
                                                        "Caenorhabditis elegans ","Drosophila melanogaster ","Danio rerio ",
                                                        "Rattus norvegicus ","Arabidopsis thaliana ","Gallus gallus "),] 

g3 <- ggplot(NIH_species_month,aes(datecount, freq, color=Species)) + geom_smooth(se=F,size=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(plot.margin = unit(c(1, 1, 4, 1), "lines"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.minor.x = element_line(color = "white", size=0.5), panel.grid.major.x = element_line(color = "gray", size=1)) + 
  scale_x_continuous(minor_breaks = seq(0:200),breaks = c(21,33,45,57,69,81,93,105,117)) #+ theme(legend.position="none",axis.title.y=element_blank())
g3 + coord_cartesian(xlim = c(21, 123), ylim = c(0,.40),expand = F)


new_model_month <-species_month_full[species_month_full$Species %in% sign_pos_freq[order(-sign_pos_freq$freq_slope),]$Species[2:9],] 
g4 <- ggplot(new_model_month,aes(datecount, freq, color=Species)) + geom_smooth(se=F,size=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(plot.margin = unit(c(1, 1, 4, 1), "lines"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.minor.x = element_line(color = "white", size=0.5), panel.grid.major.x = element_line(color = "gray", size=1)) + 
  scale_x_continuous(minor_breaks = seq(0:200),breaks = c(21,33,45,57,69,81,93,105,117)) #+ theme(legend.position="none",axis.title.y=element_blank())
g4 + coord_cartesian(xlim = c(21, 123),expand = F)


recentdanio=NIH_species_month[NIH_species_month$Species=="Danio rerio " & NIH_species_month$datecount>=81,]
summary(lm(recentdanio$freq ~ recentdanio$datecount))


for (species in sign_pos_freq$Species) {
  species_freq_plotter(species)
}

#how many species with whole genomes?
length(unique(data[data$strategy=="WGS",]$Species))

#how many iucn redlist species?
iucn <- read.table("iucn_species.txt",sep = "\t",stringsAsFactors = F)
length(which(iucn$V1 %in% data$Species))

#bases instead of experiments
bases_data <- subset(data,is.na(data$bases)==F) 

bases_data_species <- bases_data %>%
  group_by(datecount,Species) %>%
  mutate(bases = sum(bases)) %>%
  select(-taxid, -accession, -strategy) %>%
  distinct()

#supplemental I
species_merged <- merge(species_percentiles,species_full)
write.table(species_merged[c("Species","Genus","Family","Order","Class","Phylum","Kingdom","experiments","rank","percentile","percentof","exp_slope","exp_p","freq_slope","freq_p","bases")],"SupplementalI.txt",
            sep = "\t",quote=F,row.names = F,col.names = c("Species","Genus","Family","Order","Class","Phylum","Kingdom","Total Number of Experiments","Rank (Experiments)",
                                                           "Percentile (Experiments)","Percent of Database (Experiments)","Change in # of Experiments/Month","p","Change in % of Dataset/Month","p","Total Number of Bases"))


#supplemental II
write.table(species_month_total[c("date","Species","Genus","Family","Order","Class","Phylum","Kingdom","experiments","bases")],"SupplementalII.txt",sep = "\t",
            quote=F,row.names = F,col.names = c("Month","Species","Genus","Family","Order","Class","Phylum","Kingdom","Total Number of Experiments","Total Number of Bases"))


#figure S1
bases_month_relative <- bases_data_species %>% 
  group_by(datecount) %>% 
  mutate(percentile = percent_rank(bases))
bases_month_relative$category <- ifelse(bases_month_relative$percentile < .90,"other",ifelse(bases_month_relative$percentile < .99,"90","99"))
category_bases_relative <- bases_month_relative %>%
  group_by(category,datecount) %>%
  summarise(bases=sum(bases)) %>%
  ungroup() %>%
  group_by(datecount) %>%
  mutate(freq=bases/sum(bases))
category_bases_relative$category <- factor(category_bases_relative$category, levels=c("99","90","other"))

pal <- brewer.pal(3,name="YlGnBu")
g1bases <- ggplot(category_bases_relative, aes(x=datecount,y=bases,group=category,fill=category)) + geom_area(position='fill') +
  scale_fill_manual(values = c(pal[3],pal[2],pal[1])) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(),axis.text.y=element_blank(),legend.position="none") +
  scale_x_continuous(breaks=c(21,33,45,57,69,81,93,105,117))
g1bases + coord_cartesian(xlim = c(21, 123),expand = F) 

#linear model of percentage of 99th percentile over time
summary(lm(category_bases_relative[category_bases_relative$category=="99",]$freq ~ 
             category_bases_relative[category_bases_relative$category=="99",]$datecount))

#figure S3
bases_month_full <- unique(complete(bases_data_species[,c(2,3,10)],Species,datecount=21:123))
bases_month_full[is.na(bases_month_full$bases),]$bases = 0
bases_month_full <- bases_month_full %>%
  group_by(datecount) %>%
  mutate(freq = bases/sum(bases))
NIH_bases_month <-bases_month_full[bases_month_full$Species %in% c("Mus musculus ","Saccharomyces cerevisiae ",
                                                                         "Caenorhabditis elegans ","Drosophila melanogaster ","Danio rerio ",
                                                                         "Rattus norvegicus ","Arabidopsis thaliana ","Gallus gallus "),] 

g3 <- ggplot(NIH_bases_month,aes(datecount, freq, color=Species)) + geom_smooth(se=F,size=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(plot.margin = unit(c(1, 1, 4, 1), "lines"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.minor.x = element_line(color = "white", size=0.5), panel.grid.major.x = element_line(color = "gray", size=1)) + 
  scale_x_continuous(minor_breaks = seq(0:200),breaks = c(21,33,45,57,69,81,93,105,117)) + theme(legend.position="none",axis.title.y=element_blank(),axis.text.y=element_blank())
g3 + coord_cartesian(xlim = c(21, 123), ylim = c(0,.40),expand = F)







one_per_lm <- lm(tmp[tmp$MO=="99",]$freq ~ tmp[tmp$MO=="99",]$datecount)
coef(one_per_lm)[2]
lmp(one_per_lm)

plot(tmp[tmp$MO=="99",]$datecount,tmp[tmp$MO=="99",]$freq)


one_per_Phyla$Phylum

exp_R <- aggregate(exp_count ~ datecount,exp,function(x) length(x))

exp_E <- aggregate(exp_count ~ datecount,exp,function(x){diversity(x)/log(length(x))})

ggplot(exp_E, aes(datecount,exp_count)) +
  #geom_smooth(aes(exp_E$datecount,exp_E$exp_count),method = 'lm',se=F,size=0.75,color="black") +
  geom_smooth(aes(exp_E$datecount,exp_E$exp_count,fill=gg_color_hue(2)[2]), color=gg_color_hue(2)[1],size=2,show.legend = F) +
  geom_point(color=gg_color_hue(2)[1],size=3) +
  scale_y_continuous(sec.axis = sec_axis(~.*2000, name = "Species Richness")) +
  #geom_smooth(aes(exp_R$datecount,exp_R$exp_count/2000),method = 'lm',se=F,size=0.75,color="black") +
  geom_smooth(aes(exp_R$datecount,exp_R$exp_count/2000,fill=gg_color_hue(2)[1]), color=gg_color_hue(2)[2],size=2,show.legend = F) +
  geom_point(aes(exp_R$datecount,exp_R$exp_count/2000), size=3, color=gg_color_hue(2)[2]) 






bases_R <- aggregate(bases_count ~ datecount,bases,function(x) length(x))

bases_E <- aggregate(bases_count ~ datecount,bases,function(x){diversity(x)/log(length(x))})

ggplot(bases_E, aes(datecount,bases_count)) +
  geom_smooth(aes(bases_E$datecount,bases_E$bases_count,fill=gg_color_hue(2)[2]), color=gg_color_hue(2)[1],size=2,show.legend = F) + 
  geom_point(color=gg_color_hue(2)[1],size=3) +
  scale_y_continuous(sec.axis = sec_axis(~.*2000, name = "Species Richness")) +
  geom_smooth(aes(bases_R$datecount,bases_R$bases_count/2000,fill=gg_color_hue(2)[1]), color=gg_color_hue(2)[2],size=2,show.legend = F) +
  geom_point(aes(bases_R$datecount,bases_R$bases_count/2000), size=3, color=gg_color_hue(2)[2]) 




exp <- aggregate(taxid ~ datecount + date + Species + Genus + Family + Order + Class + Phylum + Kingdom,df,length)

colnames(exp)[colnames(exp)=="Species"] <- "exp_count"



bases <- aggregate(bases ~ datecount + date + Species + Genus + Family + Order + Class + Phylum + Kingdom,df,sum)
count <- merge(exp,bases)
colnames(count)[colnames(count)=="bases"] <- "bases_count"

#richness
exp_R <- aggregate(exp_count ~ datecount,count,function(x) length(x))
bases_R <- aggregate(bases_count ~ datecount,count,function(x) sum(x))

#eveness
exp_E <- aggregate(exp_count ~ datecount,count,function(x){diversity(x)/log(length(x))})
bases_E <- aggregate(bases_count ~ datecount,count,function(x){diversity(x)/log(sum(x))})

ggplot(exp_E, aes(datecount,exp_count)) +
  geom_smooth(aes(exp_E$datecount,exp_E$exp_count,fill=gg_color_hue(2)[2]), color=gg_color_hue(2)[1],size=2,show.legend = F) + 
  geom_point(color=gg_color_hue(2)[1],size=3) +
  scale_y_continuous(sec.axis = sec_axis(~.*2000, name = "Species Richness")) +
  geom_smooth(aes(exp_R$datecount,exp_R$exp_count/2000,fill=gg_color_hue(2)[1]), color=gg_color_hue(2)[2],size=2,show.legend = F) +
  geom_point(aes(exp_R$datecount,exp_R$exp_count/2000), size=3, color=gg_color_hue(2)[2]) 

ggplot(bases_E, aes(datecount,bases_count)) +
  geom_smooth(aes(bases_E$datecount,bases_E$bases_count,fill=gg_color_hue(2)[2]), color=gg_color_hue(2)[1],size=2,show.legend = F) + 
  geom_point(color=gg_color_hue(2)[1],size=3) +
  scale_y_continuous(sec.axis = sec_axis(~.*1000000000000000, name = "Species Richness")) +
  geom_smooth(aes(bases_R$datecount,bases_R$bases_count/1000000000000000,fill=gg_color_hue(2)[1]), color=gg_color_hue(2)[2],size=2,show.legend = F) +
  geom_point(aes(bases_R$datecount,bases_R$bases_count/1000000000000000), size=3, color=gg_color_hue(2)[2]) 


  


pal <- brewer.pal(3,name="YlGnBu")


exp_per <- aggregate(exp_count ~ Species, exp, sum) %>%
  mutate(percentile = percent_rank(exp_count))
exp99 <- exp_per[exp_per$percentile>=0.99,]$Species
exp90 <- exp_per[exp_per$percentile<0.99&exp_per$percentile>=0.90,]$Species
exp_month_per <- exp %>%
  group_by(datecount) %>% 
  mutate(percentile = percent_rank(exp_count))
exp_month_per$MO <- ifelse(exp_month_per$Species %in% exp99,"99",ifelse(exp_month_per$Species %in% exp90,"90","NMO"))
exp_month_per_ag <- aggregate(exp_count ~ datecount + MO, data = exp_month_per, sum, na.rm = TRUE)
exp_month_per_ag$MO <- factor(exp_month_per_ag$MO, levels=c("99","90","NMO"))
ggplot(exp_month_per_ag, aes(x=datecount,y=exp_count,group=MO,fill=MO)) + geom_area(position='fill')


bases_per <- aggregate(bases_count ~ Species, count, sum) %>%
  mutate(percentile = percent_rank(bases_count))
bases99 <- bases_per[bases_per$percentile>=0.99,]$Species
bases90 <- bases_per[bases_per$percentile<0.99&bases_per$percentile>=0.90,]$Species
bases_month_per <- count %>%
  group_by(datecount) %>% 
  mutate(percentile = percent_rank(bases_count))
bases_month_per$MO <- ifelse(bases_month_per$Species %in% bases99,"99",ifelse(bases_month_per$Species %in% bases90,"90","NMO"))
bases_month_per_ag <- aggregate(bases_count ~ datecount + MO, data = bases_month_per, sum, na.rm = TRUE)
bases_month_per_ag$MO <- factor(bases_month_per_ag$MO, levels=c("99","90","NMO"))
ggplot(bases_month_per_ag, aes(x=datecount,y=bases_count,group=MO,fill=MO)) + geom_area(position='fill')



bases_rel <- count %>% 
  group_by(datecount) %>% 
  mutate(percentile = percent_rank(bases_count))
bases_rel$MO <- ifelse(bases_rel$percentile < .90,"NMO",ifelse(bases_rel$percentile < .99,"90","99"))
bases_rel_ag<-aggregate(bases_count ~ datecount + MO, data = bases_rel, sum, na.rm = TRUE)
bases_rel_ag$MO <- factor(bases_rel_ag$MO, levels=c("99","90","NMO"))
ggplot(bases_rel_ag, aes(x=datecount,y=bases_count,group=MO,fill=MO)) + geom_area(position='fill')




M <- count(df,Species,datecount)
M <- M %>%
  group_by(datecount,Species) %>%
  summarise (N = sum(n)) %>%
  mutate(freq = N / sum(N))

NIH_M <- M[M$Species %in% c("Mus musculus ","Danio rerio ","Saccharomyces cerevisiae ","Arabidopsis thaliana ","Drosophila melanogaster ","Caenorhabditis elegans ","Rattus norvegicus "),] 
NIH_M <- complete(NIH_M,Species,datecount=21:123)
NIH_M[is.na(NIH_M)]<-0


ggplot(NIH_M,aes(datecount, freq, color=Species)) + geom_smooth(method='loess') 

NIH_M$Species <- factor(NIH_M$Species, levels=c("Mus musculus ","Saccharomyces cerevisiae ","Danio rerio ","Arabidopsis thaliana ","Drosophila melanogaster ","Caenorhabditis elegans ","Rattus norvegicus "))
ggplot() + geom_bar(aes(y = freq, x = datecount, fill = Species), data = NIH_M, stat="identity") + scale_fill_manual(values=brewer.pal(7,name="YlGnBu"))



M <- aggregate(bases ~ datecount+Species,df,sum) %>%
  group_by(datecount,Species) %>%
  summarise(N = sum(bases)) %>%
  mutate(freq = N / sum(N))
M <- complete(M,Species,datecount=21:111)
M <- unique(M)
M[is.na(M)]<-0

NIH_M <- M[M$Species %in% c("Mus musculus ","Danio rerio ","Saccharomyces cerevisiae ","Arabidopsis thaliana ","Drosophila melanogaster ","Caenorhabditis elegans ","Rattus norvegicus "),] 
NIH_M <- complete(NIH_M,Species,datecount=21:111)
M <- unique(M)
M[is.na(M)]<-0
ggplot(NIH_M,aes(datecount, freq, color=Species)) + geom_smooth(method='lm') 


tmp<-aggregate(N ~ Species,NIH_M,sum)

M <- count(df,Species,datecount)
M <- M %>%
  group_by(datecount,Species) %>%
  summarise (N = sum(n)) %>%
  mutate(freq = N / sum(N))

complete(M,Species,datecount=21:123)
M[is.na(M)]<-0


SUM = sum(count$exp_count)
slopes <- M %>%
  group_by(Species) %>%
  mutate(slope = coef(lm(N ~ datecount))[2]) %>%
  mutate(r = summary(lm(N ~ datecount))$r.squared) %>%
  mutate(p = lmp(lm(N ~ datecount))) %>%
  mutate(freqslope = coef(lm(freq ~ datecount))[2]) %>%
  mutate(freqr = summary(lm(freq ~ datecount))$r.squared) %>%
  mutate(freqp = lmp(lm(freq ~ datecount))) %>%
  mutate(count = sum(N)) %>%
  mutate(freq = (sum(N)/SUM))

M %>%
  mutate(p = lmp(lm(N ~ datecount)))

lmp(lm(M$N ~ M$datecount))

speciesums <- data.frame(Species=slopes$Species,Count=slopes$count,Freq=slopes$freq,Slope=slopes$slope,R=slopes$r,P=slopes$p,FreqSlope=slopes$freqslope,RFreq=slopes$freqr,PFreq=slopes$freqr,stringsAsFactors = F)
speciesums <- unique(speciesums)

x <- speciesums[speciesums$FreqSlope > 0 & speciesums$PFreq <= 0.05,]


x <- aggregate(bases ~ datecount,df,sum)

x <- aggregate(bases~datecount,df,sum)
plot(x,type="b")
