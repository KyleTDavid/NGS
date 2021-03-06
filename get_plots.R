#script used to generate plots and analyze data for sequencing disparity ms
#written to execute in RStudio

##NOTES##
#before filtering 56 entries removes because did not have the expected number of rows (entry error, pooled samples)

#install or load packages
library(dplyr)
library(vegan)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(ape)

#function to get nice (I think) colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
#function to darken colors
darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

#read sra input data (available https://figshare.com/articles/SRA_search_raw_output/7158659)
raw_data <- read.table(file = "data_out_filtered.txt", sep = '\t', stringsAsFactors = F, header = F, na.strings = "NA")
names(raw_data) <- c("accession","taxid","strategy","bases","date")
#filter duplicate entries (7 removed)
data_filtered <- unique(raw_data) 
#filter entries without accession (51 removed)
data_filtered <- subset(raw_data,is.na(raw_data$accession)==F) 
#filter entries without date (1019 removed)
data_filtered <- subset(data_filtered,is.na(data_filtered$date)==F) 
#read lineage output data (available https://figshare.com/articles/Lineage_output/7158692) (lines 7164 & 9244 & 11967 were edited manually to remove special characters)
lineage <- read.table(file = "lineage_out.txt", sep = '\t', stringsAsFactors = F, header = T, na.strings = c("NA ","NA"))
#add lineage data to sra data
data_lineage <- merge(data_filtered,lineage)
#filter entries without ncbi genus (17798 removed)
data_lineage_filtered <- subset(data_lineage,is.na(data_lineage$Genus)==F) 
#filter entries without ncbi species (10647 removed)
data <- subset(data_lineage_filtered,is.na(data_lineage_filtered$Species)==F) 
#just keep year and month in date column
data$date <- gsub('-\\d\\d \\d\\d:\\d\\d:\\d\\d','',data$date,perl=T)
#add column datecount which converts months into relative number value (first month = 0, second month = 1, etc)
datecount_list <- data.frame(date=sort(unique(data$date)),datecount=0:129)
data <- merge(data,datecount_list)
#filter to experiments published Jan 2010 - Dec 2018
data <- data[21<=data$datecount & data$datecount<=128,]

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
one_percent_phyla <- aggregate(Species ~ Phylum,one_percent,length)
#lets look at NIH models (when just genus is listed take the most popular species)
model_list = c("Mus musculus ","Rattus norvegicus ","Saccharomyces cerevisiae ","Schizosaccharomyces pombe ","Neurospora crassa ",
               "Dictyostelium discoideum ","Caenorhabditis elegans ","Daphnia magna ","Drosophila melanogaster ","Danio rerio ",
               "Xenopus tropicalis ","Gallus gallus ","Arabidopsis thaliana ")
#table of NIH models
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
#top 1% represented what percent of experiments in 2018?
twenty_eighteen <- species_month_total[species_month_total$datecount >= 117 & species_month_total$datecount <= 128 ,]
twenty_eighteen <- aggregate(experiments ~ Species,twenty_eighteen,sum)
twenty_eighteen$percentile <- percent_rank(twenty_eighteen$experiments)
sum(twenty_eighteen[twenty_eighteen$percentile>=0.99,]$experiments)/sum(twenty_eighteen$experiments)

#calculate species richness and evenness for each month
month_richness <- setNames(aggregate(experiments ~ datecount,species_month_total,function(x) length(x)),c("datecount","richness"))
month_evenness <- setNames(aggregate(experiments ~ datecount,species_month_total,function(x){diversity(x)/log(length(x))}),c("datecount","evenness"))
#linear model of richness over time
lm_richness <- lm(month_richness$richness ~ month_richness$datecount)
summary(lm_richness)
#richness in 2010
length(twenty_ten$Species)
#in 2018
length(twenty_eighteen$Species)
#linear model of evenness over time
summary(lm(month_evenness$evenness ~ month_evenness$datecount))
#figure 1a
g1a <- ggplot(month_evenness, aes(datecount,evenness)) +
  #geom_smooth(aes(month_evenness$datecount,month_evenness$evenness),method = 'lm',se=F,size=0.75,color="black") +
  geom_smooth(aes(month_evenness$datecount,month_evenness$evenness,fill=gg_color_hue(2)[2]), color=gg_color_hue(2)[1],size=2,show.legend = F) +
  geom_point(color=gg_color_hue(2)[1],size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(1, 1, 4, 1), "lines"), axis.title.x=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.minor.x = element_line(color = "white", size=0.5),
        panel.grid.major.x = element_line(color = "gray", size=1), axis.title.y=element_blank()) + 
  scale_x_continuous(minor_breaks = seq(0:200), breaks = seq(21,128,12)) +
  scale_y_continuous(sec.axis = sec_axis(~.*2000)) +
  #geom_smooth(aes(month_richness$datecount,month_richness$richness/2000),method = 'lm',se=F,size=0.75,color="black") +
  geom_smooth(aes(month_richness$datecount,month_richness$richness/2000,fill=gg_color_hue(2)[1]), color=gg_color_hue(2)[2],size=2,show.legend = F) +
  geom_point(aes(month_richness$datecount,month_richness$richness/2000), size=3, color=gg_color_hue(2)[2])
g1a + coord_cartesian(xlim = c(21, 128),ylim=c(0,1),expand = F) 

#figure 1b
pal <- brewer.pal(3,name="YlGnBu")
g1b <- ggplot(category_month_relative, aes(x=datecount,y=experiments,group=category,fill=category)) + geom_area(position='fill') +
  scale_fill_manual(values = c(pal[3],pal[2],pal[1])) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), legend.position="none",axis.title.y=element_blank(),axis.text.y=element_blank()) +
  scale_x_continuous(breaks = seq(21,128,12))
g1b + coord_cartesian(xlim = c(21, 128),expand = F) 

#fill data frame for every possible species month combo, 
species_month_full <- unique(complete(species_month_relative[,c(10,3,11)],Species,datecount=21:128))
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

#how many species with whole genomes?
length(unique(data[data$strategy=="WGS",]$Species))

#iucn data, file downloaded from IUCN website
iucn <- read.csv("taxonomy.csv")
iucn_species <- paste(iucn$genusName,iucn$speciesName)
sra_species <- trimws(species_data$Species)
length(intersect(iucn_species,sra_species))
iucn_species[1]

#look at seq strats
strat = aggregate(date ~ strategy,data,length)


#Supplementary experiments data table
write.table(species_month_total[c("date","Species","Genus","Family","Order","Class","Phylum","Kingdom","experiments","bases")],"SupplementalII.txt",sep = "\t",
            quote=F,row.names = F,col.names = c("Month","Species","Genus","Family","Order","Class","Phylum","Kingdom","Total Number of Experiments","Total Number of Bases"))

#Supplementary species data table 
species_merged <- merge(species_percentiles,species_full)
write.table(species_merged[c("Species","Genus","Family","Order","Class","Phylum","Kingdom","experiments","rank","percentile","percentof","exp_slope","exp_p","freq_slope","freq_p","bases")],"SupplementalI.txt",
            sep = "\t",quote=F,row.names = F,col.names = c("Species","Genus","Family","Order","Class","Phylum","Kingdom","Total Number of Experiments","Rank (Experiments)",
                                                           "Percentile (Experiments)","Percent of Database (Experiments)","Change in # of Experiments/Month","p","Change in % of Dataset/Month","p","Total Number of Bases"))


#old figure not included in final ms, part of code edited from https://bitbucket.org/caseywdunn/animal-genomes/src
dataclade=data
dataclade$clade=ifelse(is.na(data$Class)==F,data$Class,ifelse(is.na(data$Order)==F,data$Order,data$Genus))
dataclade$clade=trimws(dataclade$clade)
cladeprop = function(x) {
  length(x)/length(dataclade$clade)
}
clade_data=aggregate(date ~ clade, dataclade, cladeprop)
clade_data_filt=clade_data[clade_data$date>0.01,]
species_data=aggregate(date ~ Species + clade, dataclade, cladeprop)
topspecies_data <- species_data[species_data$clade %in% clade_data_filt$clade,] %>%
  group_by(clade) %>%
  top_n(n = 1, wt = date)
tre = read.tree("tree.txt")
plot( tre, edge.width=5, show.tip.label = F, x.lim=c(0,50),node.depth = 2)
x_offset = 7 # x position of start of bars
one_size = 110 # how big the bar should be for a value of 1
abline(v=( x_offset + .05*one_size), col="gray")
abline(v=( x_offset + .1*one_size), col="gray")
abline(v=( x_offset + .15*one_size), col="gray")
abline(v=( x_offset + .2*one_size), col="gray")
abline(v=( x_offset + .25*one_size), col="gray")
abline(v=( x_offset + .3*one_size), col="gray")
abline(v=( x_offset + .35*one_size), col="gray")
abline(v=( x_offset + .4*one_size), col="gray")
rownames(clade_data_filt) <- clade_data_filt$clade
clade_data_filt = clade_data_filt[match(tre$tip.label, clade_data_filt$clade),]
segments(x_offset,1:nrow(clade_data_filt), x_offset + pmax(clade_data_filt$date*one_size, 0 ), 1:nrow(clade_data_filt), lwd=9, col=gg_color_hue(12), lend=1 ) 
rownames(topspecies_data) <- topspecies_data$clade
topspecies_data = topspecies_data[match(tre$tip.label, topspecies_data$clade),]
segments(x_offset,1:nrow(clade_data_filt), x_offset + pmax(topspecies_data$date*one_size, 0 ), 1:nrow(topspecies_data), lwd=8, col=darken(gg_color_hue(12),1.4), lend=1 ) 

#old figure not included
NIH_species_month <-species_month_full[species_month_full$Species %in% c("Mus musculus ","Saccharomyces cerevisiae ",
                                                                         "Caenorhabditis elegans ","Drosophila melanogaster ","Danio rerio ",
                                                                         "Rattus norvegicus ","Arabidopsis thaliana ","Gallus gallus "),] 
g3 <- ggplot(NIH_species_month,aes(datecount, freq, color=Species)) + geom_smooth(se=F,size=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(plot.margin = unit(c(1, 1, 4, 1), "lines"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.minor.x = element_line(color = "white", size=0.5), panel.grid.major.x = element_line(color = "gray", size=1)) + 
  scale_x_continuous(minor_breaks = seq(0:200),breaks = seq(21,128,12)) #+ theme(legend.position="none",axis.title.y=element_blank())
g3 + coord_cartesian(xlim = c(21, 128), ylim = c(0,.40),expand = F)

#bases instead of experiments
bases_data <- subset(data,is.na(data$bases)==F) 

bases_data_species <- bases_data %>%
  group_by(datecount,Species) %>%
  mutate(bases = sum(bases)) %>%
  select(-taxid, -accession, -strategy) %>%
  distinct()

#old figure with bases
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

#old with bases 
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



