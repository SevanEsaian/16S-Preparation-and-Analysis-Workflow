#Load Necessary Libraries --------
library(ecodist)
library(vegan)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(devtools)
library(introdataviz)
library(ggsignif)
library(ggpubr)
library(corncob)

#Split Violin Plots -----------
env.data <- read.csv("environmental data ready for plotting.csv")
View(env.data)

ggplot(env.data, aes(x=Depth, y=SurfaceAreacm, fill=Age)) +
  geom_violin(position=position_dodge(0.9)) +
  geom_boxplot(binaxis='y', width=.05, 
               position=position_dodge(0.9), color="white", outlier.shape = NA) +
  theme_bw() +
  coord_flip() +
  theme(axis.text = element_text(size = 18)) +
  theme(axis.title = element_text(size = 20)) +
  scale_fill_manual(values = c('yellow3', 'darkblue')) +
  ylab(bquote('Surface Area '(cm^2))) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  annotate("text", x=3.4, y=8, label = "(A)", size=12) +
  annotate("text", x=3.2, y=512, label = "A", size=10) +
  annotate("text", x=2.75, y=512, label = "A", size=10) +
  annotate("text", x=2.2, y=512, label = "A", size=10) +
  annotate("text", x=1.8, y=512, label = "B", size=10) +
  annotate("text", x=1.2, y=510, label = "A", size=10) +
  annotate("text", x=0.8, y=510, label = "A", size=10) +
  theme(legend.position = "none")

ggplot(env.data, aes(x=Depth, y=FvFm, fill=Age)) +
  geom_violin(position=position_dodge(0.9)) +
  geom_boxplot(binaxis='y', width=.05, 
               position=position_dodge(0.9), color="white", outlier.shape = NA) +
  theme_bw() +
  coord_flip() +
  theme(axis.text = element_text(size = 18)) +
  theme(axis.title = element_text(size = 20)) +
  scale_fill_manual(values = c('yellow3', 'darkblue')) +
  ylab("Maximum Potential Quantfum Efficiency (Fv/Fm)") +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  ylim(0.2, 0.7) +
  annotate("text", x=3.4, y=0.2, label = "(B)", size=12) +
  annotate("text", x=3.2, y=0.7, label = "A", size=10) +
  annotate("text", x=2.8, y=0.7, label = "B", size=10) +
  annotate("text", x=2.2, y=0.7, label = "B", size=10) +
  annotate("text", x=1.8, y=0.7, label = "C", size=10) +
  annotate("text", x=1.2, y=0.7, label = "A", size=10) +
  annotate("text", x=0.8, y=0.7, label = "D", size=10) +
  theme(legend.position = "none")

ggplot(env.data, aes(x=Depth, y=ChlperC, fill=Age)) +
  geom_violin(position=position_dodge(0.9)) +
  geom_boxplot(binaxis='y', width=.05, 
               position=position_dodge(0.9), color="white", outlier.shape = NA) +
  theme_bw() +
  coord_flip() +
  theme(axis.text = element_text(size = 18)) +
  theme(axis.title = element_text(size = 20)) +
  scale_fill_manual(values = c('yellow3', 'darkblue')) +
  ylab("Chla:C (mg)") +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  ylim(0.001, 0.011) +
  annotate("text", x=3.4, y=0.001, label = "(C)", size=12) +
  annotate("text", x=3.2, y=0.011, label = "A", size=10) +
  annotate("text", x=2.8, y=0.011, label = "B", size=10) +
  annotate("text", x=2.2, y=0.011, label = "A", size=10) +
  annotate("text", x=1.8, y=0.011, label = "B", size=10) +
  annotate("text", x=1.2, y=0.011, label = "B", size=10) +
  annotate("text", x=0.8, y=0.011, label = "A", size=10) +
  theme(legend.position = "none")

ggplot(env.data, aes(x=Depth, y=CN, fill=Age)) +
  geom_violin(position=position_dodge(0.9)) +
  geom_boxplot(binaxis='y', width=.05, 
               position=position_dodge(0.9), color="white", outlier.shape = NA) +
  theme_bw() +
  coord_flip() +
  theme(axis.text = element_text(size = 18)) +
  theme(axis.title = element_text(size = 20)) +
  scale_fill_manual(values = c('yellow3', 'darkblue')) +
  ylab("C:N (%)") +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  ylim(6, 29) +
  annotate("text", x=3.4, y=6, label = "(D)", size=12) +
  annotate("text", x=3.25, y=29, label = "A", size=10) +
  annotate("text", x=2.8, y=29, label = "A", size=10) +
  annotate("text", x=2.25, y=29, label = "B", size=10) +
  annotate("text", x=1.8, y=29, label = "B", size=10) +
  annotate("text", x=1.25, y=29, label = "B", size=10) +
  annotate("text", x=0.8, y=29, label = "B", size=10) +
  theme(legend.text=element_text(size=20),
        legend.title = element_text(size=20),
        legend.position = "none")
##
##Redoing the diversity metrics calculations - 060723------------
DivMet <- read.csv("060723-DiversityMetrics-Barplot.csv")
View(DivMet)
DivMet %>%
  ggplot( aes(x=Age, y=AlphaDiversity, fill=Age)) +
  geom_boxplot(alpha=0.7, outlier.shape = NA, color="gray10") +
  scale_fill_manual(values = c('yellow3', 'darkblue')) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, alpha=0.5) +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()) 



DivMet %>%
  ggplot( aes(x=Depth, y=AlphaDiversity, fill=Depth)) +
  geom_boxplot(alpha=0.7, outlier.shape = NA, color="gray35") +
  scale_fill_manual(values = c('purple3', 'forestgreen')) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, alpha=0.5) +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()) 
#
#Run an NMDS of just blade samples no water ----------
blade.reads <- read.csv("031923 - Blade Rarefied Reads - For PCoA.csv")
View(blade.reads)
blade.reads.edit <- blade.reads %>%
  tibble::column_to_rownames("UniqueID")
View(blade.reads.edit)

samp.coords <- read.csv("031923 - Sample Coords - For NMDS.csv")
View(samp.coords)
samp.coords$UniqueID <- as.character(samp.coords$UniqueID)

asv.dist<-vegdist(blade.reads.edit, method='bray')
asv.dist
asv.div<-adonis2(asv.dist~as.factor(samp.coords$Depth), data=samp.coords, permutations=9999)
asv.div

age.div<-adonis2(asv.dist~as.factor(samp.coords$Age), data=samp.coords, permutations=9999)
age.div

BC.NMDS.shell<- function(dataframe) {
  metaMDS(dataframe, distance="bray", k=2, trymax=1000)
} 
BC.NMDS <- BC.NMDS.shell(blade.reads.edit)  
plot(BC.NMDS)
BC.NMDS$stress
View(BC.NMDS$species)
##write.csv(BC.NMDS$species, "ASV NMDS coordinates for biplot.csv")
##write.csv(BC.NMDS$points, "UniqueID NMDS coordinates for biplot.csv")
BC.NMDS_sample_coords <- scores(BC.NMDS, display = "sites") %>%
  as.data.frame() %>%
  #make my row names into a column
  rownames_to_column("UniqueID") %>%
  #join data frame with metadata
  full_join(., samp.coords, by = "UniqueID")
View(BC.NMDS_sample_coords)
##write.csv(BC.NMDS_sample_coords, "ReVisualizing NMDS based on Bart Suggestions.csv")

Mature_Coords <- subset(BC.NMDS_sample_coords, Age=="Mature")
View(Mature_Coords)
#write.csv(Mature_Coords, "Mature Coords Revisualizing.csv")
ggplot(Mature_Coords, aes(x=NMDS1, y=NMDS2, color=Age, shape=Depth)) +
  geom_point(alpha=0.9, size=5, color='darkblue') +
  theme_bw() +
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 20)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  xlim(-0.92,0.72) + ylim(-0.85,0.8) +
  theme(legend.position="") +
  annotate("text", x=-0.9, y=0.8, label = "B", size=7) +
  annotate("text", x=-0.75, y=-0.85, label="Mature microbiomes", size=7)
ggsave("NMDS Mature.tiff")

Juvenile_Coords <- subset(BC.NMDS_sample_coords, Age=="Juvenile")
#write.csv(Juvenile_Coords, "Juvenile Coords Revisualizing.csv")
View(Juvenile_Coords)
ggplot(Juvenile_Coords, aes(x=NMDS1, y=NMDS2, color=Age, shape=Depth)) +
  geom_point(alpha=0.9, size=5, color='yellow3') +
  theme_bw() +
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 20)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  xlim(-0.92,0.72) + ylim(-0.85,0.8) +
  theme(legend.position="") +
  annotate("text", x=-0.9, y=0.8, label = "A", size=7) +
  annotate("text", x=-0.75, y=-0.85, label="Juvenile microbiomes", size=7) 
ggsave("NMDS juvenile.tiff")

ggplot(BC.NMDS_sample_coords, aes(x=NMDS1, y=NMDS2, color=Age, shape=Depth)) +
  geom_point(alpha=0.9, size=5) +
  theme_bw() +
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 20)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  xlim(-0.92,0.72) + ylim(-0.85,0.8) 
ggsave("NMDS Mature.tiff")



#ReRun Diversity Metrics by combining all juveniles vs all matures---------
DivMet <- read.csv("031923 - ReDone Diversity Metrics - For Barplots.csv")
View(DivMet)
ggplot(data=DivMet, aes(x=))

ASVRich <- subset(DivMet, DiversityMetrics=="ASV Richness")
View(ASVRich)

ggplot(data=ASVRich, aes(x=Age, y=Average, fill=Age)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() +
  scale_fill_manual(values=c('yellow3','darkblue')) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  ylab("ASV Richness") +
  geom_errorbar(aes(ymin=Average, ymax=Average+Stdev), width=.2,
                position=position_dodge(.9)) +
  ylim(0,1500) + 
  annotate("text", x=2, y=1450, label = "***", size=8) +
  annotate("text", x=0.5, y=1490, label = "A", size=8)

Evenness <- subset(DivMet, DiversityMetrics=="Pielou Evenness")

ggplot(data=Evenness, aes(x=Age, y=Average, fill=Age)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() +
  scale_fill_manual(values=c('yellow3','darkblue')) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  ylab("Pielou's Evenness") +
  geom_errorbar(aes(ymin=Average, ymax=Average+Stdev), width=.2,
                position=position_dodge(.9)) +
  ylim(0,0.85) + 
  annotate("text", x=2, y=0.85, label = "***", size=8) +
  annotate("text", x=0.5, y=0.85, label = "B", size=8)

Disper <- subset(DivMet, DiversityMetrics=="Beta Dispersion")
ggplot(data=Disper, aes(x=Age, y=Average, fill=Age)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() +
  scale_fill_manual(values=c('yellow3','darkblue')) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  ylab("Beta Dispersion") +
  geom_errorbar(aes(ymin=Average, ymax=Average+Stdev), width=.2,
                position=position_dodge(.9)) +
  ylim(0,0.45) + 
  annotate("text", x=1, y=0.45, label = "**", size=8) +
  annotate("text", x=0.5, y=0.45, label = "C", size=8) 

#Diversity Metrics violin plots ----------
DiversityMetrics <- read.csv("Diversity Metrics.csv")
View(DiversityMetrics)

ggplot(data=DiversityMetrics, aes(y=Depth, x=ASVRichness, fill=Age)) +
  geom_violin(position=position_dodge(0.9)) +
  geom_boxplot(binaxis='y', width=.05, 
               position=position_dodge(0.9), color="white", outlier.shape = NA) +
  theme_bw() +
  scale_fill_manual(values=c('yellow3','darkblue')) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) +
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor=element_blank()) +
  annotate("text", x=500, y=3.4, label = "(A)", size=12) +
  annotate("text", x=1635, y=3.2, label = "A", size=10) +
  annotate("text", x=1635, y=2.8, label = "B", size=10) +
  annotate("text", x=1635, y=2.2, label = "C", size=10) +
  annotate("text", x=1635, y=1.8, label = "B", size=10) +
  annotate("text", x=1635, y=1.2, label = "C", size=10) +
  annotate("text", x=1635, y=0.8, label = "B", size=10)

ggplot(data=DiversityMetrics, aes(y=Depth, x=SimpsonEvenness, fill=Age)) +
  geom_violin(position=position_dodge(0.9)) +
  geom_boxplot(binaxis='y', width=.05, 
               position=position_dodge(0.9), color="white", outlier.shape = NA) +
  theme_bw() +
  scale_fill_manual(values=c('yellow3','darkblue')) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor=element_blank()) +
  annotate("text", x=0.01, y=3.4, label = "(E)", size=12) +
  annotate("text", x=0.17, y=3.2, label = "AC", size=10) +
  annotate("text", x=0.17, y=2.8, label = "B", size=10) +
  annotate("text", x=0.17, y=2.2, label = "C", size=10) +
  annotate("text", x=0.17, y=1.8, label = "B", size=10) +
  annotate("text", x=0.17, y=1.2, label = "C", size=10) +
  annotate("text", x=0.17, y=0.8, label = "BC", size=10) +
  xlab ("Simpson Index")

ggplot(data=DiversityMetrics, aes(y=Depth, x=ShannonDiversity, fill=Age)) +
  geom_violin(position=position_dodge(0.9)) +
  geom_boxplot(binaxis='y', width=.05, 
               position=position_dodge(0.9), color="white", outlier.shape = NA) +
  theme_bw() +
  scale_fill_manual(values=c('yellow3','darkblue')) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) +
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor=element_blank()) +
  annotate("text", x=3, y=3.3, label = "(D)", size=12) +
  annotate("text", x=6.2, y=3.2, label = "A", size=10) +
  annotate("text", x=6.2, y=2.8, label = "B", size=10) +
  annotate("text", x=6.2, y=2.2, label = "A", size=10) +
  annotate("text", x=6.2, y=1.8, label = "B", size=10) +
  annotate("text", x=6.2, y=1.2, label = "A", size=10) +
  annotate("text", x=6.2, y=0.8, label = "B", size=10) +
  xlab("Shannon Diversity")

ggplot(data=DiversityMetrics, aes(y=Depth, x=PielouEvenness, fill=Age)) +
  geom_violin(position=position_dodge(0.9)) +
  geom_boxplot(binaxis='y', width=.05, 
               position=position_dodge(0.9), color="white", outlier.shape = NA) +
  theme_bw() +
  scale_fill_manual(values=c('yellow3','darkblue')) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) +
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor=element_blank()) +
  annotate("text", x=0.56 , y=3.4, label = "(B)", size=12) +
  annotate("text", x=0.84 , y=3.2, label = "A", size=10) +
  annotate("text", x=0.84 , y=2.8, label = "B", size=10) +
  annotate("text", x=0.84 , y=2.2, label = "A", size=10) +
  annotate("text", x=0.84 , y=1.8, label = "B", size=10) +
  annotate("text", x=0.84 , y=1.2, label = "A", size=10) +
  annotate("text", x=0.84 , y=0.8, label = "B", size=10) +
  xlab("Pielou's Evenness")

ggplot(data=DiversityMetrics, aes(y=Depth, x=BetaDisper, fill=Age)) +
  geom_violin(position=position_dodge(0.9)) +
  geom_boxplot(binaxis='y', width=.05, 
               position=position_dodge(0.9), color="white", outlier.shape = NA) +
  theme_bw() +
  scale_fill_manual(values=c('yellow3','darkblue')) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) +
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor=element_blank()) +
  xlab("Beta Dispersion") +
  annotate("text", x=0.1 , y=3.4, label = "(C)", size=12) +
  annotate("text", x=0.63 , y=3.2, label = "A", size=10) +
  annotate("text", x=0.63 , y=2.8, label = "B", size=10) +
  annotate("text", x=0.63 , y=2.2, label = "A", size=10) +
  annotate("text", x=0.63 , y=1.8, label = "AB", size=10) +
  annotate("text", x=0.63 , y=1.2, label = "A", size=10) +
  annotate("text", x=0.63 , y=0.8, label = "B", size=10)
#
##CCA Plotting ---------
Photophys <- read.csv("Giant Kelp Photophysiology Data.csv")
Photophys %>% remove_rownames %>% column_to_rownames(var="UniqueID") 
Photophys <- column_to_rownames(Photophys, "UniqueID")
View(Photophys)

RA.data <- read.csv("PercentRelativeAbundance - Edited and Ready For CCA.csv")
RA.data %>% remove_rownames %>% column_to_rownames(var="UniqueID") 
RA.data <- column_to_rownames(RA.data, "UniqueID")
View(RA.data)

vare.cca <- cca(RA.data ~ PC+PN+CN+mgChl_to_mgC+SurfaceArea+FvFm , data=Photophys)
vare.cca
plot(vare.cca)
summary(vare.cca)

ASVcoord <- read.csv("CCA ASV Coordinates.csv")
View(ASVcoord)
Physcoord <- read.csv("CCA Physiological Coordinates.csv")
View(Physcoord)

ggplot(ASVcoord) +
  geom_point(mapping = aes(x=CCA1, y=CCA2, color=Age, shape=Depth), size=4, alpha=0.6) +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  scale_color_manual(values = c('yellow3', 'darkblue')) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text=element_text(size=20),
        legend.title = element_text(size=20)) +
  geom_segment(data=Physcoord,
               aes(x=0, xend=CCA1, y=0, yend=CCA2),
               type="closed", color="gray", alpha=0.5) +
  geom_text(data = Physcoord, aes(x = CCA1, y = CCA2, label = Physiology),
            size = 5) +
  xlab("CCA1 (33.3%)") + ylab("CCA2 (21.6%)")
ggsave("CCA.tiff")
##Corncob retry (4/27/23) - with just top genera > 1%. 2 comparisons (mat vs juv, surf vs subsurf)-----
View(gathered_rawreads)
justgenera_rawreads = subset(gathered_rawreads, select = c(ASV, genus, UniqueID, RawReadCount))
View(justgenera_rawreads)
genus_gathered_rawreads_agg <- aggregate(justgenera_rawreads$RawReadCount, by=list(genus=justgenera_rawreads$genus, UniqueID=justgenera_rawreads$UniqueID), FUN=sum)
View(genus_gathered_rawreads_agg)
genus_gathered_rawreads_agg_sub <- subset(genus_gathered_rawreads_agg, genus %in%
                                   c('Afipia',
                                     'Agaribacter',
                                     'Alcanivorax',
                                     'Algitalea',
                                     'Aliivibrio',
                                     'Alteromonadaceae',
                                     'Amylibacter',
                                     'Aureispira',
                                     'Bacteriovoracaceae',
                                     'Bdellovibrionaceae',
                                     'Blastopirellula',
                                     'Brucella',
                                     'Cellvibrionaceae',
                                     'Clade_I',
                                     'Clade_Ia',
                                     'Cocleimonas',
                                     'Colwelliaceae',
                                     'Croceitalea',
                                     'Cryomorphaceae',
                                     'Cyanobiaceae',
                                     'DEV007',
                                     'Dokdonia',
                                     'Flavicella',
                                     'Flavobacteriaceae',
                                     'Fusobacteriaceae',
                                     'Granulosicoccaceae',
                                     'Granulosicoccus',
                                     'Halieaceae',
                                     'Haloferula',
                                     'Halomonadaceae',
                                     'Halomonas',
                                     'Hyphomonadaceae',
                                     'Kiritimatiellaceae',
                                     'Lentimonas',
                                     'Leucothrix',
                                     'Lewinella',
                                     'Litorimonas',
                                     'Lutimonas',
                                     'Marinomonas',
                                     'Microtrichaceae',
                                     'Moraxellaceae',
                                     'NS4_marine_group',
                                     'OM60(NOR5)_clade',
                                     'Persicirhabdus',
                                     'Phycisphaera',
                                     'Phycisphaeraceae',
                                     'Pirellulaceae',
                                     'Planktomarina',
                                     'Polaribacter',
                                     'Propionigenium',
                                     'Pseudoalteromonadaceae',
                                     'Pseudoalteromonas',
                                     'Psychrobacter',
                                     'Psychromonadaceae',
                                     'Psychromonas',
                                     'Puniceicoccaceae',
                                     'Rapidithrix',
                                     'Rhizobiaceae',
                                     'Rhodobacteraceae',
                                     'Rhodopirellula',
                                     'Robiginitomaculum',
                                     'Roseibacillus',
                                     'Rubidimonas',
                                     'Rubinisphaeraceae',
                                     'Rubripirellula',
                                     'Rubritalea',
                                     'Rubritaleaceae',
                                     'Saprospiraceae',
                                     'SAR92_clade',
                                     'Sedimentitalea',
                                     'Stenotrophomonas',
                                     'Sva0996_marine_group',
                                     'Synechococcus_CC9902',
                                     'Tenacibaculum',
                                     'Terasakiellaceae',
                                     'Thiotrichaceae',
                                     'Vibrio',
                                     'Vibrionaceae',
                                     'Wenyingzhuangia',
                                     'Xanthobacteraceae'))
View(genus_gathered_rawreads_agg_sub)

RA.0.5 <- read.csv("Percent Relative Abundance Reads - For NMDS Test.csv")


##write.csv(genus_gathered_rawreads_agg_sub, "Top 0.5% Genera for Barplots and Corncob.csv")
desired_genera_summedrawreads <- tidyr::spread(genus_gathered_rawreads_agg_sub, genus, x)
View(desired_genera_summedrawreads) ##THIS IS OUR DESIRED GENERA DATAFRAME FROM WHICH WE'RE WORKING
seqtab.reretry <- desired_genera_summedrawreads
View(seqtab.reretry)
metadata.reretry <- read.csv("metadata_juvenileandmature_justgenera.csv")
taxa.reretry <- read.csv("genus_metadata_top1%.csv")
seqtab.reretry.edit <- seqtab.reretry %>% tibble::column_to_rownames("UniqueID")
View(seqtab.reretry.edit)
seqtab.reretry.mat <- as.matrix(seqtab.reretry.edit)
seqtab.reretry.mat.tab <- otu_table(seqtab.reretry.mat, taxa_are_rows = FALSE)
seqtab.reretry.mat.tab

taxa.reretry.seq <- taxa.reretry %>% tibble::column_to_rownames("genus")
taxa.reretry.seq.mat <- as.matrix(taxa.reretry.seq)
tax.reretry.seq.mat.tab <- tax_table(taxa.reretry.seq.mat)
tax.reretry.seq.mat.tab


metadata.reretry <- metadata.reretry %>% tibble::column_to_rownames("UniqueID")
met.reretry <- sample_data(metadata.reretry)
met.reretry

attempt.phylo <- phyloseq(seqtab.reretry.mat.tab, met.reretry, tax.reretry.seq.mat.tab)
attempt.phylo

depth.ctrlage.genus <- differentialTest(formula = ~ depth_surf_subsurf + age,
                                        phi.formula = ~ depth_surf_subsurf,
                                        formula_null = ~ age,
                                        phi.formula_null = ~ 1,
                                        data = attempt.phylo, 
                                        test = "Wald", boot = FALSE, 
                                        fdr_cutoff = 0.05)   
depth.ctrlage.genus$all_models
plot(depth.ctrlage.genus) +
theme(legend.position = "none") +
  ylab(NULL)

age.ctrlage.genus <- differentialTest(formula = ~ age + depth_surf_subsurf,
                                        phi.formula = ~ age ,
                                        formula_null = ~ depth_surf_subsurf,
                                        phi.formula_null = ~ 1,
                                        data = attempt.phylo, 
                                        test = "Wald", boot = FALSE, 
                                        fdr_cutoff = 0.05)  

age.ctrlage.genus$all_models
plot(age.ctrlage.genus) +
  theme(legend.position = "none") +
  ylab(NULL)

#Plot corncob genera for depth and age ---------
genera.corncob <- read.csv("Corncob Genera For Manuscript 042723.csv")
View(genera.corncob)
genera.corncob.depth.abundance <- subset(genera.corncob, Comparison=="Depth" & Analysis=="Abundance")
View(genera.corncob.depth.abundance)
genera.corncob.depth.abundance$Genus <- factor(genera.corncob.depth.abundance$Genus, levels=c('Rubritalea',
                                                                                              'Roseibacillus',
                                                                                              'Persicirhabdus',
                                                                                              'Rubripirellula',
                                                                                              'Phycisphaera',
                                                                                              'Psychrobacter',
                                                                                              'OM60(NOR5)_clade',
                                                                                              'Leucothrix',
                                                                                              'Propionigenium',
                                                                                              'Lewinella',
                                                                                              'Dokdonia',
                                                                                              'Hellea',
                                                                                              'Amylibacter',
                                                                                              'Afipia' ))
ggplot(genera.corncob.depth.abundance, aes(x=Average, y=Genus, color=Association)) +
  geom_point() +
  geom_pointrange(aes(xmin=Average-Stdev, xmax=Average+Stdev)) + theme_bw() +
  geom_vline(xintercept=c(0,0), linetype="dotted") +
  xlim(-0.66, 1.66) +
  scale_color_manual(values=c("gray", "blue")) +
  xlab("Regression Coefficient") +
  theme(axis.text.y = element_text(face="italic")) +
  ggtitle("Differential abundance by depth (surface vs subsurface)") +
  theme(plot.title = element_text(hjust = 0.5))

genera.corncob.depth.variability <- subset(genera.corncob, Comparison=="Depth" & Analysis=="Variability")
View(genera.corncob.depth.variability)
genera.corncob.depth.variability$Genus <- factor(genera.corncob.depth.variability$Genus, levels=c('Rubritalea',
                                                                                              'Roseibacillus',
                                                                                              'Persicirhabdus',
                                                                                              'Rubripirellula',
                                                                                              'Phycisphaera',
                                                                                              'Psychrobacter',
                                                                                              'OM60(NOR5)_clade',
                                                                                              'Leucothrix',
                                                                                              'Propionigenium',
                                                                                              'Lewinella',
                                                                                              'Dokdonia',
                                                                                              'Hellea',
                                                                                              'Amylibacter',
                                                                                              'Afipia' ))
ggplot(genera.corncob.depth.variability, aes(x=Average, y=Genus, color=Association)) +
  geom_point() +
  geom_pointrange(aes(xmin=Average-Stdev, xmax=Average+Stdev)) + theme_bw() +
  geom_vline(xintercept=c(0,0), linetype="dotted") +
  xlim(-1.8, 1.9) +
  scale_color_manual(values=c("gray", "blue")) +
  xlab("Regression Coefficient") +
  theme(axis.text.y = element_text(face="italic")) +
  ggtitle("Differential variability by depth (surface vs subsurface)") +
  theme(plot.title = element_text(hjust = 0.5))

genera.corncob.age.abundance <- subset(genera.corncob, Comparison=="Age" & Analysis=="Abundance")
View(genera.corncob.age.abundance)
genera.corncob.age.abundance$Genus <- factor(genera.corncob.age.abundance$Genus, levels=c('Rubritalea',
                                                                                              'Rubripirellula',
                                                                                              'Blastopirellula',
                                                                                              'Phycisphaera',
                                                                                              'OM60(NOR5)_clade',
                                                                                              'Leucothrix',
                                                                                              'Psychrobacter',
                                                                                              'Granulosicoccus',
                                                                                              'Propionigenium',
                                                                                              'Synechococcus_CC902',
                                                                                              'Lewinella',
                                                                                              'Dokdonia',
                                                                                              'Croceitalea',
                                                                                              'Flavicella',
                                                                                              'Amylibacter',
                                                                                              'Afipia' ))
ggplot(genera.corncob.age.abundance, aes(x=Average, y=Genus, color=Association)) +
  geom_point() +
  geom_pointrange(aes(xmin=Average-Stdev, xmax=Average+Stdev)) + theme_bw() +
  geom_vline(xintercept=c(0,0), linetype="dotted") +
  xlim(-0.66, 1.66) +
  scale_color_manual(values=c("gray", "blue")) +
  xlab("Regression Coefficient") +
  theme(axis.text.y = element_text(face="italic"))+
  ggtitle("Differential abundance by age (juvenile vs mature)") +
  theme(plot.title = element_text(hjust = 0.5))

genera.corncob.age.variability <- subset(genera.corncob, Comparison=="Age" & Analysis=="Variability")
View(genera.corncob.age.variability)
genera.corncob.age.variability$Genus <- factor(genera.corncob.age.variability$Genus, levels=c('Rubritalea',
                                                                                          'Rubripirellula',
                                                                                          'Blastopirellula',
                                                                                          'Phycisphaera',
                                                                                          'OM60(NOR5)_clade',
                                                                                          'Leucothrix',
                                                                                          'Psychrobacter',
                                                                                          'Granulosicoccus',
                                                                                          'Propionigenium',
                                                                                          'Synechococcus_CC902',
                                                                                          'Lewinella',
                                                                                          'Dokdonia',
                                                                                          'Croceitalea',
                                                                                          'Flavicella',
                                                                                          'Amylibacter',
                                                                                          'Afipia' ))
ggplot(genera.corncob.age.variability, aes(x=Average, y=Genus, color=Association)) +
  geom_point() +
  geom_pointrange(aes(xmin=Average-Stdev, xmax=Average+Stdev)) + theme_bw() +
  geom_vline(xintercept=c(0,0), linetype="dotted") +
  xlim(-1.8, 1.9) +
  scale_color_manual(values=c("gray", "blue")) +
  xlab("Regression Coefficient") +
  theme(axis.text.y = element_text(face="italic"))+
  ggtitle("Differential variability by age (juvenile vs mature)") +
  theme(plot.title = element_text(hjust = 0.5))

##
##Corncob Surface vs Subsurface--------
library(corncob)
library(phyloseq)
seqtab.retry <- read.csv("seqtab_juvenileandmature_sequencereads.csv")
metadata.retry <- read.csv("metadata_juvenileandmature.csv")
taxa.retry <- read.csv("asv_taxa_sequencereplacedASV.csv")
View(metadata.retry)
View(taxa.retry)
seqtab.retry.edit <- seqtab.retry %>% tibble::column_to_rownames("SampleID")
View(seqtab.retry.edit)
seqtab.retry.mat <- as.matrix(seqtab.retry.edit)
seqtab.retry.mat.tab <- otu_table(seqtab.retry.mat, taxa_are_rows = FALSE)
seqtab.retry.mat.tab

taxa.retry.seq <- taxa.retry %>% tibble::column_to_rownames("Sequence")
taxa.retry.seq.mat <- as.matrix(taxa.retry.seq)
tax.retry.seq.mat.tab <- tax_table(taxa.retry.seq.mat)
tax.retry.seq.mat.tab


metadata.retry <- metadata.retry %>% tibble::column_to_rownames("SampleID")
View(metadata)
met.retry <- sample_data(metadata.retry)
met.retry

attempt.phylo.fam <- phyloseq(seqtab.retry.mat.tab, met.retry, tax.retry.seq.mat.tab)
attempt.phylo.fam

corncob.phylo.Family <- attempt.phylo.fam %>% tax_glom("Family")
corncob.phylo.Family

corncob.phylo.Genus <- attempt.phylo.fam %>% tax_glom("Genus")
corncob.phylo.Genus


Surf_SubSurf_DiffTest_Family_condition.dv <- differentialTest(formula = ~ depth_surf_subsurf ,
                                                 formula_null = ~ ,
                                                 phi.formula = ~ depth_surf_subsurf + age,
                                                 phi.formula_null = ~ depth_surf_subsurf + age,
                                                 data = corncob.phylo.Family, 
                                                 test = "LRT", boot = FALSE, 
                                                 fdr_cutoff = 0.05)   
Surf_SubSurf_DiffTest_Family_condition.dv$all_models
Surf_SubSurf_DiffTest_Genus_condition.dv
otu_to_taxonomy(OTU = Surf_SubSurf_DiffTest_Family_condition.dv$significant_taxa, data = corncob.phylo.Family)
plot(Surf_SubSurf_DiffTest_Family_condition.dv, level = "Family") + 
  theme(axis.text.y = element_text(face= "italic")) +
  theme(text = element_text(size = 6)) +
  geom_point(aes(colour = x >0), size = 3) +
  scale_colour_manual(values = setNames(c('forestgreen','darkslateblue'),c(T, F)))+
  theme(legend.position = "none") +
  ylab(NULL)

Surf_SubSurf_DiffTest_Genus_condition.dv <- differentialTest(formula = ~ depth_surf_subsurf,
                                                              formula_null = ~ 1,
                                                              phi.formula = ~ depth_surf_subsurf + age,
                                                              phi.formula_null = ~ depth_surf_subsurf,
                                                              data = corncob.phylo.Genus, 
                                                              test = "LRT", boot = FALSE, 
                                                              fdr_cutoff = 0.05)   
Surf_SubSurf_DiffTest_Genus_condition.dv$all_models
otu_to_taxonomy(OTU = Surf_SubSurf_DiffTest_Genus_condition.dv$significant_taxa, data = corncob.phylo.Genus)
plot(Surf_SubSurf_DiffTest_Genus_condition.dv, level = "Genus") + 
  theme(axis.text.y = element_text(face= "italic")) +
  theme(text = element_text(size = 6)) +
  geom_point(aes(colour = x >0), size = 3) +
  scale_colour_manual(values = setNames(c('forestgreen','darkslateblue'),c(T, F)))+
  theme(legend.position = "none") +
  ylab(NULL)

corncob.surf.subsurf <- read.csv("032123_Corncob_SurfacevsSubsurface.csv")
View(corncob.surf.subsurf)

corncob.surf.subsurf$Family <- factor(corncob.surf.subsurf$Family,
                                      levels=c( "Rubritaleaceae",
                                                "Puniceicoccaceae",
                                                "Pirellulaceae",
                                                "Phycisphaeraceae",
                                                "Fusobacteriaceae",
                                                "Cyanobiaceae",
                                                "Psychromodaceae",
                                                "Xanthobacteraceae" ))
ggplot(corncob.surf.subsurf, aes(x=Family, y=Average, 
                        ymin=Average-Stdev, ymax=Average+Stdev)) +
  geom_pointrange(aes(color=Depth), size=1) +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  geom_hline(yintercept=0, linetype='dashed', col='gray', size=1) +
  theme(axis.text.x = element_text(face= "italic")) +
  theme(axis.text.x = element_text(angle = 0),
        axis.text.y = element_text(angle=0)) +
  theme(axis.text = element_text(size = 5)) +
  theme(axis.text = element_text(size = 13)) +
  theme(axis.title = element_text(size = 13)) +
  scale_color_manual(values = c('darkgreen', 'purple3')) +
  labs(x="", y="Differential Variation") +
  theme(axis.text.y = element_text(face = "italic")) +
  ylim(-0.18, 0.31) +
  coord_flip() 


##Need to consult the hard drive and find the correct asv coords file

##Corncob by Age -------
View(corncob.phylo.Family)
Age_DiffTest_Family_condition.dv <- differentialTest(formula = ~ age + depth,
                                                              formula_null = ~ depth,
                                                              phi.formula = ~ age + depth,
                                                              phi.formula_null = ~ age + depth,
                                                              data = corncob.phylo.Family, 
                                                              test = "LRT", boot = FALSE, 
                                                              fdr_cutoff = 0.05)   
Age_DiffTest_Family_condition.dv$all_models
otu_to_taxonomy(OTU = Age_DiffTest_Family_condition.dv$significant_taxa, data = corncob.phylo.Family)
plot(Age_DiffTest_Family_condition.dv, level = "Family") + 
  theme(axis.text.y = element_text(face= "italic")) +
  theme(text = element_text(size = 6)) +
  geom_point(aes(colour = x >0), size = 3) +
  scale_colour_manual(values = setNames(c('forestgreen','darkslateblue'),c(T, F)))+
  theme(legend.position = "none") +
  ylab(NULL)

##ASVcoords for Barplots ----------
RA.ASV <- read.csv("032223_RA_MicrobiomeOnly_ForASVMerge.csv")
View(RA.ASV)
ASV.coord <- read.csv("032223_ASVcoords_ForMerge.csv")
View(ASV.coord)
ASV.merge <- merge(RA.ASV, ASV.coord, by="ASV")
View(ASV.merge)
write.csv(ASV.merge, "ASV merge attempt.csv")

#032223 - Emily Corncob Edits -------------
#Must test for each depth separately
seqtab.microbiome <- read.csv("seqtab_microbiomes_sequencereads.csv")
metadata.microbiome <- read.csv("metadata_microbiome_ForCorncob.csv")
taxa.microbiome <- read.csv("asv_taxa_sequencereplacedASV.csv")

seqtab.microbiome.edit <- seqtab.microbiome %>% tibble::column_to_rownames("SampleID")
seqtab.microbiome.mat <- as.matrix(seqtab.microbiome.edit)
seqtab.microbiome.mat.tab <- otu_table(seqtab.microbiome.mat, taxa_are_rows = FALSE)
seqtab.microbiome.mat.tab

taxa.microbiome.seq <- taxa.microbiome %>% tibble::column_to_rownames("Sequence")
taxa.microbiome.seq.mat <- as.matrix(taxa.microbiome.seq)
tax.microbiome.seq.mat.tab <- tax_table(taxa.microbiome.seq.mat)
tax.microbiome.seq.mat.tab


metadata.microbiome <- metadata.microbiome %>% tibble::column_to_rownames("SampleID")
met.microbiome <- sample_data(metadata.microbiome)
met.microbiome

microbiome.phylo.fam <- phyloseq(seqtab.microbiome.mat.tab, met.microbiome, tax.microbiome.seq.mat.tab)
microbiome.phylo.fam

corncob.microbiome.phylo.Family <- microbiome.phylo.fam %>% tax_glom("Family")
corncob.microbiome.phylo.Family

corncob.microbiome.phylo.Phylum <- microbiome.phylo.fam %>% tax_glom("Phylum")
corncob.microbiome.phylo.Phylum

?differentialTest
DepthwithAge <- differentialTest(formula = ~ Depth + Age,
                        phi.formula = ~ Depth,
                        formula_null = ~ Age, phi.formula_null = ~ 1,
                        data = corncob.microbiome.phylo.Family,
                        test = "Wald", boot = FALSE, 
                        fdr_cutoff = 0.05)

plot(DepthwithAge, level="Family")

AgewithDepth <-  differentialTest(formula = ~ Age + Depth,
                                  phi.formula = ~ Age,
                                  formula_null = ~ Depth, phi.formula_null = ~ 1,
                                  data = corncob.microbiome.phylo.Family,
                                  test = "Wald", boot = FALSE, 
                                  fdr_cutoff = 0.05, 
                                  full_output = TRUE)
plot(AgewithDepth, level="Family")
AgewithDepth$significant_models
summary(AgewithDepth)

#04202023 RawReads Stacked Barplots Build dataframe structure----------
#It would be best to make a table of relative abundance, raw reads, and stacked barplots per family
#First we need raw reads with family information per sample

rawreads <- read.csv("Raw Reads with UniqueID and ASV.csv")
rel.ab <- read.csv("Microbiome Only - Relative Abundance for Barplots and Corncob.csv")
View(rel.ab)
View(rawreads)
View(ASV.coord)

rawreads_asvcoords <- merge(ASV.coord, rawreads, by = "ASV")
rel.ab_asvcoords <- merge(ASV.coord, rel.ab, by = "ASV")
View(rel.ab_asvcoords)
View(rawreads_asvcoords) 

gathered_rawreads <- tidyr::gather(rawreads_asvcoords, "UniqueID", "RawReadCount", 8:120)
gathered_rel.ab <- tidyr::gather(rel.ab_asvcoords, "UniqueID", "RelativeAbundance", 8:93)
View(gathered_rel.ab)
View(gathered_rawreads)

gathered_rawreads_agg <- aggregate(gathered_rawreads$RawReadCount, by=list(family=gathered_rawreads$family, UniqueID=gathered_rawreads$UniqueID), FUN=sum)
gathered_rel.ab_agg <- aggregate(gathered_rel.ab$RelativeAbundance, by=list(genus=gathered_rel.ab$genus, UniqueID=gathered_rel.ab$UniqueID), FUN=sum)

View(gathered_rawreads_agg)
View(gathered_rel.ab_agg)
gathered_rawreads_agg$RawReadSum <- gathered_rawreads_agg$x
View(gathered_rawreads_agg)
gathered_rawreads_agg <- gathered_rawreads_agg[,-3]
View(gathered_rawreads_agg)

spread_raw_reads_sum <- tidyr::spread(gathered_rawreads_agg, UniqueID, RawReadSum)
View(spread_raw_reads_sum)

##write.csv(spread_raw_reads_sum, "RawReads Summaries per UI for stackedbarplot.csv")
write.csv(gathered_rel.ab_agg, "Gathered Genus RA%.csv")
rawreads.sbp <- read_csv("RawReads Summaries per UI for stackedbarplot.csv")
View(rawreads.sbp)
gather.rawreads.sbp <- tidyr::gather(rawreads.sbp, "UI","RawReadCounts", 2:114)
View(gather.rawreads.sbp)

samp.coords.sbp <- read_csv("Sample Coordinates - Juvenile and Mature Samples - For Stacked Barplot.csv")
View(samp.coords.sbp)
stackedbarplot.readyforsubset <- merge(samp.coords.sbp, gather.rawreads.sbp, by = "UI")
View(stackedbarplot.readyforsubset)
###write.csv(stackedbarplot.readyforsubset, "Stackedbarplot data - rawreads only.csv")
#Generated surface juvenile stackedbarplots of rawreads --------
#First, subset juvenile surface blades 
stackedbarplot.readyforsubset <- read_csv("Stackedbarplot data - rawreads only.csv")
surfacejuv.rawreads.family <- subset(stackedbarplot.readyforsubset, SampleCoords=="surfaceJuvenile")
View(surfacejuv.rawreads.family)
#now you need to subset your desired families + the "all other families" category


surfacejuv.domfam.rawreads <- subset(surfacejuv.rawreads.family, family %in% 
                                          c('Microtrichaceae',
                                            'Cryomorphaceae',
                                            'Cyclobacteriaceae',
                                            'Flammeovirgaceae',
                                            'Flavobacteriaceae',
                                            'Saprospiraceae',
                                            'Bacteriovoracaceae',
                                            'Bdellovibrioceae',
                                            'JG30-KF-CM45',
                                            'Cyanobiaceae',
                                            'Marinomodaceae',
                                            'Fusobacteriaceae',
                                            'Phycisphaeraceae',
                                            'Pirellulaceae',
                                            'Rubinisphaeraceae',
                                            'Alcanivoracaceae',
                                            'Alteromodaceae',
                                            'Cellvibrioceae',
                                            'Clade_I',
                                            'Colwelliaceae',
                                            'Granulosicoccaceae',
                                            'Halieaceae',
                                            'Halomodaceae',
                                            'Hyphomodaceae',
                                            'Moraxellaceae',
                                            'Porticoccaceae',
                                            'Pseudoalteromodaceae',
                                            'Psychromodaceae',
                                            'Rhizobiaceae',
                                            'Rhodobacteraceae',
                                            'Sphingomodaceae',
                                            'Terasakiellaceae',
                                            'Thiotrichaceae',
                                            'Vibrioceae',
                                            'Xanthobacteraceae',
                                            'Xanthomodaceae',
                                            'DEV007',
                                            'Kiritimatiellaceae',
                                            'Puniceicoccaceae',
                                            'Rubritaleaceae',
                                            'AllOtherFamilies'))

surfacejuv.domfam.rawreads$family <- factor(surfacejuv.domfam.rawreads$family, levels=c('Microtrichaceae',
                                                                                  'Cryomorphaceae',
                                                                                  'Cyclobacteriaceae',
                                                                                  'Flammeovirgaceae',
                                                                                  'Flavobacteriaceae',
                                                                                  'Saprospiraceae',
                                                                                  'Bacteriovoracaceae',
                                                                                  'Bdellovibrioceae',
                                                                                  'JG30-KF-CM45',
                                                                                  'Cyanobiaceae',
                                                                                  'Marinomodaceae',
                                                                                  'Fusobacteriaceae',
                                                                                  'Phycisphaeraceae',
                                                                                  'Pirellulaceae',
                                                                                  'Rubinisphaeraceae',
                                                                                  'Alcanivoracaceae',
                                                                                  'Alteromodaceae',
                                                                                  'Cellvibrioceae',
                                                                                  'Clade_I',
                                                                                  'Colwelliaceae',
                                                                                  'Granulosicoccaceae',
                                                                                  'Halieaceae',
                                                                                  'Halomodaceae',
                                                                                  'Hyphomodaceae',
                                                                                  'Moraxellaceae',
                                                                                  'Porticoccaceae',
                                                                                  'Pseudoalteromodaceae',
                                                                                  'Psychromodaceae',
                                                                                  'Rhizobiaceae',
                                                                                  'Rhodobacteraceae',
                                                                                  'Sphingomodaceae',
                                                                                  'Terasakiellaceae',
                                                                                  'Thiotrichaceae',
                                                                                  'Vibrioceae',
                                                                                  'Xanthobacteraceae',
                                                                                  'Xanthomodaceae',
                                                                                  'DEV007',
                                                                                  'Kiritimatiellaceae',
                                                                                  'Puniceicoccaceae',
                                                                                  'Rubritaleaceae',
                                                                                  'AllOtherFamilies'))

View(surfacejuv.domfam.rawreads)



surfacejuv.rawreads.family <- subset(stackedbarplot.readyforsubset, SampleCoords=="surfaceJuvenile")
View(surfacejuv.rawreads.family)
#now you need to subset your desired families + the "all other families" category

#Generated middle juvenile stackedbarplots of rawreads --------
#First, subset juvenile middle blades 
middlejuv.rawreads.family <- subset(stackedbarplot.readyforsubset, SampleCoords=="middleJuvenile")
View(middlejuv.rawreads.family)
#now you need to subset your desired families + the "all other families" category
middlejuv.domfam.rawreads <- subset(middlejuv.rawreads.family, family %in% 
                                       c('Microtrichaceae',
                                         'Cryomorphaceae',
                                         'Cyclobacteriaceae',
                                         'Flammeovirgaceae',
                                         'Flavobacteriaceae',
                                         'Saprospiraceae',
                                         'Bacteriovoracaceae',
                                         'Bdellovibrioceae',
                                         'JG30-KF-CM45',
                                         'Cyanobiaceae',
                                         'Marinomodaceae',
                                         'Fusobacteriaceae',
                                         'Phycisphaeraceae',
                                         'Pirellulaceae',
                                         'Rubinisphaeraceae',
                                         'Alcanivoracaceae',
                                         'Alteromodaceae',
                                         'Cellvibrioceae',
                                         'Clade_I',
                                         'Colwelliaceae',
                                         'Granulosicoccaceae',
                                         'Halieaceae',
                                         'Halomodaceae',
                                         'Hyphomodaceae',
                                         'Moraxellaceae',
                                         'Porticoccaceae',
                                         'Pseudoalteromodaceae',
                                         'Psychromodaceae',
                                         'Rhizobiaceae',
                                         'Rhodobacteraceae',
                                         'Sphingomodaceae',
                                         'Terasakiellaceae',
                                         'Thiotrichaceae',
                                         'Vibrioceae',
                                         'Xanthobacteraceae',
                                         'Xanthomodaceae',
                                         'DEV007',
                                         'Kiritimatiellaceae',
                                         'Puniceicoccaceae',
                                         'Rubritaleaceae',
                                         'AllOtherFamilies'))

middlejuv.domfam.rawreads$family <- factor(middlejuv.domfam.rawreads$family, levels=c('Microtrichaceae',
                                                                                        'Cryomorphaceae',
                                                                                        'Cyclobacteriaceae',
                                                                                        'Flammeovirgaceae',
                                                                                        'Flavobacteriaceae',
                                                                                        'Saprospiraceae',
                                                                                        'Bacteriovoracaceae',
                                                                                        'Bdellovibrioceae',
                                                                                        'JG30-KF-CM45',
                                                                                        'Cyanobiaceae',
                                                                                        'Marinomodaceae',
                                                                                        'Fusobacteriaceae',
                                                                                        'Phycisphaeraceae',
                                                                                        'Pirellulaceae',
                                                                                        'Rubinisphaeraceae',
                                                                                        'Alcanivoracaceae',
                                                                                        'Alteromodaceae',
                                                                                        'Cellvibrioceae',
                                                                                        'Clade_I',
                                                                                        'Colwelliaceae',
                                                                                        'Granulosicoccaceae',
                                                                                        'Halieaceae',
                                                                                        'Halomodaceae',
                                                                                        'Hyphomodaceae',
                                                                                        'Moraxellaceae',
                                                                                        'Porticoccaceae',
                                                                                        'Pseudoalteromodaceae',
                                                                                        'Psychromodaceae',
                                                                                        'Rhizobiaceae',
                                                                                        'Rhodobacteraceae',
                                                                                        'Sphingomodaceae',
                                                                                        'Terasakiellaceae',
                                                                                        'Thiotrichaceae',
                                                                                        'Vibrioceae',
                                                                                        'Xanthobacteraceae',
                                                                                        'Xanthomodaceae',
                                                                                        'DEV007',
                                                                                        'Kiritimatiellaceae',
                                                                                        'Puniceicoccaceae',
                                                                                        'Rubritaleaceae',
                                                                                        'AllOtherFamilies'))

View(middlejuv.domfam.rawreads)


#Generated bottom juvenile stackedbarplots of rawreads --------
#First, subset juvenile middle blades 
bottomjuv.rawreads.family <- subset(stackedbarplot.readyforsubset, SampleCoords=="bottomJuvenile")
View(bottomjuv.rawreads.family)
#now you need to subset your desired families + the "all other families" category
bottomjuv.domfam.rawreads <- subset(bottomjuv.rawreads.family, family %in% 
                                      c('Microtrichaceae',
                                        'Cryomorphaceae',
                                        'Cyclobacteriaceae',
                                        'Flammeovirgaceae',
                                        'Flavobacteriaceae',
                                        'Saprospiraceae',
                                        'Bacteriovoracaceae',
                                        'Bdellovibrioceae',
                                        'JG30-KF-CM45',
                                        'Cyanobiaceae',
                                        'Marinomodaceae',
                                        'Fusobacteriaceae',
                                        'Phycisphaeraceae',
                                        'Pirellulaceae',
                                        'Rubinisphaeraceae',
                                        'Alcanivoracaceae',
                                        'Alteromodaceae',
                                        'Cellvibrioceae',
                                        'Clade_I',
                                        'Colwelliaceae',
                                        'Granulosicoccaceae',
                                        'Halieaceae',
                                        'Halomodaceae',
                                        'Hyphomodaceae',
                                        'Moraxellaceae',
                                        'Porticoccaceae',
                                        'Pseudoalteromodaceae',
                                        'Psychromodaceae',
                                        'Rhizobiaceae',
                                        'Rhodobacteraceae',
                                        'Sphingomodaceae',
                                        'Terasakiellaceae',
                                        'Thiotrichaceae',
                                        'Vibrioceae',
                                        'Xanthobacteraceae',
                                        'Xanthomodaceae',
                                        'DEV007',
                                        'Kiritimatiellaceae',
                                        'Puniceicoccaceae',
                                        'Rubritaleaceae',
                                        'AllOtherFamilies'))

bottomjuv.domfam.rawreads$family <- factor(bottomjuv.domfam.rawreads$family, levels=c('Microtrichaceae',
                                                                                      'Cryomorphaceae',
                                                                                      'Cyclobacteriaceae',
                                                                                      'Flammeovirgaceae',
                                                                                      'Flavobacteriaceae',
                                                                                      'Saprospiraceae',
                                                                                      'Bacteriovoracaceae',
                                                                                      'Bdellovibrioceae',
                                                                                      'JG30-KF-CM45',
                                                                                      'Cyanobiaceae',
                                                                                      'Marinomodaceae',
                                                                                      'Fusobacteriaceae',
                                                                                      'Phycisphaeraceae',
                                                                                      'Pirellulaceae',
                                                                                      'Rubinisphaeraceae',
                                                                                      'Alcanivoracaceae',
                                                                                      'Alteromodaceae',
                                                                                      'Cellvibrioceae',
                                                                                      'Clade_I',
                                                                                      'Colwelliaceae',
                                                                                      'Granulosicoccaceae',
                                                                                      'Halieaceae',
                                                                                      'Halomodaceae',
                                                                                      'Hyphomodaceae',
                                                                                      'Moraxellaceae',
                                                                                      'Porticoccaceae',
                                                                                      'Pseudoalteromodaceae',
                                                                                      'Psychromodaceae',
                                                                                      'Rhizobiaceae',
                                                                                      'Rhodobacteraceae',
                                                                                      'Sphingomodaceae',
                                                                                      'Terasakiellaceae',
                                                                                      'Thiotrichaceae',
                                                                                      'Vibrioceae',
                                                                                      'Xanthobacteraceae',
                                                                                      'Xanthomodaceae',
                                                                                      'DEV007',
                                                                                      'Kiritimatiellaceae',
                                                                                      'Puniceicoccaceae',
                                                                                      'Rubritaleaceae',
                                                                                      'AllOtherFamilies'))

View(bottomjuv.domfam.rawreads)


#Generated surface mature stackedbarplots of rawreads --------
#First, subset mature surface blades 
surfacemat.rawreads.family <- subset(stackedbarplot.readyforsubset, SampleCoords=="surfaceMature")
View(surfacemat.rawreads.family)
#now you need to subset your desired families + the "all other families" category
surfacemat.domfam.rawreads <- subset(surfacemat.rawreads.family, family %in% 
                                      c('Microtrichaceae',
                                        'Cryomorphaceae',
                                        'Cyclobacteriaceae',
                                        'Flammeovirgaceae',
                                        'Flavobacteriaceae',
                                        'Saprospiraceae',
                                        'Bacteriovoracaceae',
                                        'Bdellovibrioceae',
                                        'JG30-KF-CM45',
                                        'Cyanobiaceae',
                                        'Marinomodaceae',
                                        'Fusobacteriaceae',
                                        'Phycisphaeraceae',
                                        'Pirellulaceae',
                                        'Rubinisphaeraceae',
                                        'Alcanivoracaceae',
                                        'Alteromodaceae',
                                        'Cellvibrioceae',
                                        'Clade_I',
                                        'Colwelliaceae',
                                        'Granulosicoccaceae',
                                        'Halieaceae',
                                        'Halomodaceae',
                                        'Hyphomodaceae',
                                        'Moraxellaceae',
                                        'Porticoccaceae',
                                        'Pseudoalteromodaceae',
                                        'Psychromodaceae',
                                        'Rhizobiaceae',
                                        'Rhodobacteraceae',
                                        'Sphingomodaceae',
                                        'Terasakiellaceae',
                                        'Thiotrichaceae',
                                        'Vibrioceae',
                                        'Xanthobacteraceae',
                                        'Xanthomodaceae',
                                        'DEV007',
                                        'Kiritimatiellaceae',
                                        'Puniceicoccaceae',
                                        'Rubritaleaceae',
                                        'AllOtherFamilies'))

surfacemat.domfam.rawreads$family <- factor(surfacemat.domfam.rawreads$family, levels=c('Microtrichaceae',
                                                                                      'Cryomorphaceae',
                                                                                      'Cyclobacteriaceae',
                                                                                      'Flammeovirgaceae',
                                                                                      'Flavobacteriaceae',
                                                                                      'Saprospiraceae',
                                                                                      'Bacteriovoracaceae',
                                                                                      'Bdellovibrioceae',
                                                                                      'JG30-KF-CM45',
                                                                                      'Cyanobiaceae',
                                                                                      'Marinomodaceae',
                                                                                      'Fusobacteriaceae',
                                                                                      'Phycisphaeraceae',
                                                                                      'Pirellulaceae',
                                                                                      'Rubinisphaeraceae',
                                                                                      'Alcanivoracaceae',
                                                                                      'Alteromodaceae',
                                                                                      'Cellvibrioceae',
                                                                                      'Clade_I',
                                                                                      'Colwelliaceae',
                                                                                      'Granulosicoccaceae',
                                                                                      'Halieaceae',
                                                                                      'Halomodaceae',
                                                                                      'Hyphomodaceae',
                                                                                      'Moraxellaceae',
                                                                                      'Porticoccaceae',
                                                                                      'Pseudoalteromodaceae',
                                                                                      'Psychromodaceae',
                                                                                      'Rhizobiaceae',
                                                                                      'Rhodobacteraceae',
                                                                                      'Sphingomodaceae',
                                                                                      'Terasakiellaceae',
                                                                                      'Thiotrichaceae',
                                                                                      'Vibrioceae',
                                                                                      'Xanthobacteraceae',
                                                                                      'Xanthomodaceae',
                                                                                      'DEV007',
                                                                                      'Kiritimatiellaceae',
                                                                                      'Puniceicoccaceae',
                                                                                      'Rubritaleaceae',
                                                                                      'AllOtherFamilies'))

View(surfacemat.domfam.rawreads)


#Generated middle mature stackedbarplots of rawreads --------
#First, subset mature surface blades 
middlemat.rawreads.family <- subset(stackedbarplot.readyforsubset, SampleCoords=="middleMature")
View(middlemat.rawreads.family)
#now you need to subset your desired families + the "all other families" category
middlemat.domfam.rawreads <- subset(middlemat.rawreads.family, family %in% 
                                       c('Microtrichaceae',
                                         'Cryomorphaceae',
                                         'Cyclobacteriaceae',
                                         'Flammeovirgaceae',
                                         'Flavobacteriaceae',
                                         'Saprospiraceae',
                                         'Bacteriovoracaceae',
                                         'Bdellovibrioceae',
                                         'JG30-KF-CM45',
                                         'Cyanobiaceae',
                                         'Marinomodaceae',
                                         'Fusobacteriaceae',
                                         'Phycisphaeraceae',
                                         'Pirellulaceae',
                                         'Rubinisphaeraceae',
                                         'Alcanivoracaceae',
                                         'Alteromodaceae',
                                         'Cellvibrioceae',
                                         'Clade_I',
                                         'Colwelliaceae',
                                         'Granulosicoccaceae',
                                         'Halieaceae',
                                         'Halomodaceae',
                                         'Hyphomodaceae',
                                         'Moraxellaceae',
                                         'Porticoccaceae',
                                         'Pseudoalteromodaceae',
                                         'Psychromodaceae',
                                         'Rhizobiaceae',
                                         'Rhodobacteraceae',
                                         'Sphingomodaceae',
                                         'Terasakiellaceae',
                                         'Thiotrichaceae',
                                         'Vibrioceae',
                                         'Xanthobacteraceae',
                                         'Xanthomodaceae',
                                         'DEV007',
                                         'Kiritimatiellaceae',
                                         'Puniceicoccaceae',
                                         'Rubritaleaceae',
                                         'AllOtherFamilies'))

middlemat.domfam.rawreads$family <- factor(middlemat.domfam.rawreads$family, levels=c('Microtrichaceae',
                                                                                        'Cryomorphaceae',
                                                                                        'Cyclobacteriaceae',
                                                                                        'Flammeovirgaceae',
                                                                                        'Flavobacteriaceae',
                                                                                        'Saprospiraceae',
                                                                                        'Bacteriovoracaceae',
                                                                                        'Bdellovibrioceae',
                                                                                        'JG30-KF-CM45',
                                                                                        'Cyanobiaceae',
                                                                                        'Marinomodaceae',
                                                                                        'Fusobacteriaceae',
                                                                                        'Phycisphaeraceae',
                                                                                        'Pirellulaceae',
                                                                                        'Rubinisphaeraceae',
                                                                                        'Alcanivoracaceae',
                                                                                        'Alteromodaceae',
                                                                                        'Cellvibrioceae',
                                                                                        'Clade_I',
                                                                                        'Colwelliaceae',
                                                                                        'Granulosicoccaceae',
                                                                                        'Halieaceae',
                                                                                        'Halomodaceae',
                                                                                        'Hyphomodaceae',
                                                                                        'Moraxellaceae',
                                                                                        'Porticoccaceae',
                                                                                        'Pseudoalteromodaceae',
                                                                                        'Psychromodaceae',
                                                                                        'Rhizobiaceae',
                                                                                        'Rhodobacteraceae',
                                                                                        'Sphingomodaceae',
                                                                                        'Terasakiellaceae',
                                                                                        'Thiotrichaceae',
                                                                                        'Vibrioceae',
                                                                                        'Xanthobacteraceae',
                                                                                        'Xanthomodaceae',
                                                                                        'DEV007',
                                                                                        'Kiritimatiellaceae',
                                                                                        'Puniceicoccaceae',
                                                                                        'Rubritaleaceae',
                                                                                        'AllOtherFamilies'))

View(middlemat.domfam.rawreads)


#Generated bottom mature stackedbarplots of rawreads --------
#First, subset mature surface blades 
bottommat.rawreads.family <- subset(stackedbarplot.readyforsubset, SampleCoords=="bottomMature")
View(bottommat.rawreads.family)
#now you need to subset your desired families + the "all other families" category
bottommat.domfam.rawreads <- subset(bottommat.rawreads.family, family %in% 
                                      c('Microtrichaceae',
                                        'Cryomorphaceae',
                                        'Cyclobacteriaceae',
                                        'Flammeovirgaceae',
                                        'Flavobacteriaceae',
                                        'Saprospiraceae',
                                        'Bacteriovoracaceae',
                                        'Bdellovibrioceae',
                                        'JG30-KF-CM45',
                                        'Cyanobiaceae',
                                        'Marinomodaceae',
                                        'Fusobacteriaceae',
                                        'Phycisphaeraceae',
                                        'Pirellulaceae',
                                        'Rubinisphaeraceae',
                                        'Alcanivoracaceae',
                                        'Alteromodaceae',
                                        'Cellvibrioceae',
                                        'Clade_I',
                                        'Colwelliaceae',
                                        'Granulosicoccaceae',
                                        'Halieaceae',
                                        'Halomodaceae',
                                        'Hyphomodaceae',
                                        'Moraxellaceae',
                                        'Porticoccaceae',
                                        'Pseudoalteromodaceae',
                                        'Psychromodaceae',
                                        'Rhizobiaceae',
                                        'Rhodobacteraceae',
                                        'Sphingomodaceae',
                                        'Terasakiellaceae',
                                        'Thiotrichaceae',
                                        'Vibrioceae',
                                        'Xanthobacteraceae',
                                        'Xanthomodaceae',
                                        'DEV007',
                                        'Kiritimatiellaceae',
                                        'Puniceicoccaceae',
                                        'Rubritaleaceae',
                                        'AllOtherFamilies'))

bottommat.domfam.rawreads$family <- factor(bottommat.domfam.rawreads$family, levels=c('Microtrichaceae',
                                                                                      'Cryomorphaceae',
                                                                                      'Cyclobacteriaceae',
                                                                                      'Flammeovirgaceae',
                                                                                      'Flavobacteriaceae',
                                                                                      'Saprospiraceae',
                                                                                      'Bacteriovoracaceae',
                                                                                      'Bdellovibrioceae',
                                                                                      'JG30-KF-CM45',
                                                                                      'Cyanobiaceae',
                                                                                      'Marinomodaceae',
                                                                                      'Fusobacteriaceae',
                                                                                      'Phycisphaeraceae',
                                                                                      'Pirellulaceae',
                                                                                      'Rubinisphaeraceae',
                                                                                      'Alcanivoracaceae',
                                                                                      'Alteromodaceae',
                                                                                      'Cellvibrioceae',
                                                                                      'Clade_I',
                                                                                      'Colwelliaceae',
                                                                                      'Granulosicoccaceae',
                                                                                      'Halieaceae',
                                                                                      'Halomodaceae',
                                                                                      'Hyphomodaceae',
                                                                                      'Moraxellaceae',
                                                                                      'Porticoccaceae',
                                                                                      'Pseudoalteromodaceae',
                                                                                      'Psychromodaceae',
                                                                                      'Rhizobiaceae',
                                                                                      'Rhodobacteraceae',
                                                                                      'Sphingomodaceae',
                                                                                      'Terasakiellaceae',
                                                                                      'Thiotrichaceae',
                                                                                      'Vibrioceae',
                                                                                      'Xanthobacteraceae',
                                                                                      'Xanthomodaceae',
                                                                                      'DEV007',
                                                                                      'Kiritimatiellaceae',
                                                                                      'Puniceicoccaceae',
                                                                                      'Rubritaleaceae',
                                                                                      'AllOtherFamilies'))

View(bottommat.domfam.rawreads)


#Stacked barplot ggplot -------
ggplot(bottommat.domfam.rawreads, aes(x=FrondandBlade, y=RawReadCounts, fill=family)) +
  theme_bw() +
  geom_bar(position="stack", stat="identity", color="black") +
  xlab("Bottom Mature Samples") + ylab("Raw Read Count") +
  scale_fill_manual(values = c("gray75", "lightblue", "lightblue2",
                               "lightblue3", "lightblue4", "dodgerblue",
                               "orange1", "orange3", "yellow1",
                               "green3", "pink1", "pink4",
                               "red1", "red3", "darkred",
                               "aquamarine", "aquamarine3", "cyan",
                               "cyan3", "darkseagreen", "darkseagreen2",
                               "royalblue1", "royalblue4", "slateblue1",
                               "slateblue3", "turquoise", "turquoise1",
                               "turquoise4", "tan", "tan1",
                               "tan3", "goldenrod", "goldenrod1",
                               "goldenrod4", "azure", "azure2",
                               "deeppink", "deeppink3", "bisque", 
                               "bisque3", "black", "white")) +
  ylim(0,140000)


##Barplot of Genera based on corncob results ----------
corn.genera.ra <- read.csv("Corncob genera for barplot and anova comparison.csv")
View(corn.genera.ra)
corn.genera.depth.ra <- subset(corn.genera.ra, Category=="Subsurface" | Category=="Surface")

ggplot(corn.genera.depth.ra, aes(x=Average, y=Genera, fill=Category)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(xmin=Average, xmax=Average+Stdev), width=.2,
                position=position_dodge(.9)) +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) 

Gen.aov <- read.csv("ANOVA Percent RA Genera Averages by depth and age >1%.csv")
View(Gen.aov)
Afipia.aov <- subset(Gen.aov, genus=="Afipia")
View(Afipia.aov)

# Compute the analysis of variance
Afipia.aov.sum <- aov(RelativeAbundance ~ Category, data = Afipia.aov)
# Summary of the analysis
summary(Afipia.aov.sum)

Amylibacter <- subset(Gen.aov, genus=="Amylibacter")
Amylibacter.aov <- aov(RelativeAbundance~Category, data=Amylibacter)
summary(Amylibacter.aov)

Dokdonia <- subset(Gen.aov, genus=="Dokdonia")
Dokdonia.aov <- aov(RelativeAbundance~Category, data=Dokdonia)
summary(Dokdonia.aov)

Hellea <- subset(Gen.aov, genus=="Hellea")
View(Hellea)
Hellea.aov <- aov(RelativeAbundance~Category, data=Hellea)
summary(Hellea.aov)

Leucothrix <- subset(Gen.aov, genus=="Leucothrix")
Leuco.aov <- aov(RelativeAbundance~Category, data=Leucothrix)
summary(Leuco.aov)

Lewinella <- subset(Gen.aov, genus=="Lewinella")
Lew.aov <- aov(RelativeAbundance~Category, data=Lewinella)
summary(Lew.aov)

OM60 <- subset(Gen.aov, genus=="OM60(NOR5)_clade")
OM60.aov <- aov(RelativeAbundance~Category, data=OM60)
summary(OM60.aov)

Persicirhabdus <- subset(Gen.aov, genus=="Persicirhabdus")
Persicirhabdus.aov <- aov(RelativeAbundance~Category, data=Persicirhabdus)
summary(Persicirhabdus.aov)

Phycisphaera <- subset(Gen.aov, genus=="Phycisphaera")
Phycis.aov <- aov(RelativeAbundance~Category, data=Phycisphaera)
summary(Phycis.aov)

Propionigenium <- subset(Gen.aov, genus=="Propionigenium")
Prop.aov <- aov(RelativeAbundance~Category, data=Propionigenium)
summary(Prop.aov)

Psychrobacter <- subset(Gen.aov, genus=="Psychrobacter")
Psych.aov <- aov(RelativeAbundance~Category, data=Psychrobacter)
summary(Psych.aov)

Roseibacillus <- subset(Gen.aov, genus=="Roseibacillus")
Rose.aov <- aov(RelativeAbundance~Category, data=Roseibacillus)
summary(Rose.aov)

Rubripirellula <- subset(Gen.aov, genus=="Rubripirellula")
Rubri.aov <- aov(RelativeAbundance~Category, data=Rubripirellula)
summary(Rubri.aov)

Rubritalea <- subset(Gen.aov, genus=="Rubritalea")
Rubri.aov <- aov(RelativeAbundance~Category, data=Rubritalea)
summary(Rubri.aov)

View(Gen.aov)
Genus.Depth <- subset(Gen.aov, genus=="Afipia" |
                        genus=="Amylibacter" |
                        genus=="Dokdonia" |
                        genus=="Hellea" |
                        genus=="Leucothrix" |
                        genus=="Lewinella" |
                        genus=="OM60(NOR5)_clade" |
                        genus=="Persicirhabdus" |
                        genus=="Phycisphaera" |
                        genus=="Propionigenium" |
                        genus=="Psychrobacter"|
                        genus=="Roseibacillus" |
                        genus=="Rubripirellula" |
                        genus=="Rubritalea")

#Age corncob barplots
View(corn.genera.ra)
Age.corn.genera.ra <- subset(corn.genera.ra, Category=="Juvenile" | Category=="Mature")
View(Age.corn.genera.ra)

ggplot(Age.corn.genera.ra, aes(x=Average, y=Genera, fill=Category)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(xmin=Average, xmax=Average+Stdev), width=.2,
                position=position_dodge(.9)) +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) 

View(Gen.aov)
Genus.Age <- subset(Gen.aov, genus=="Afipia" |
                        genus=="Amylibacter" |
                        genus=="Blastopirellula" |
                        genus=="Croceitalea" |
                        genus=="Dokdonia" |
                        genus=="Flavicella" |
                        genus=="Granulosicoccus" |
                        genus=="Leucothrix" |
                        genus=="Lewinella" |
                        genus=="OM60(NOR5)_clade" |
                        genus=="Phycisphaera"|
                        genus=="Propionigenium" |
                        genus=="Psychrobacter" |
                        genus=="Rubripirellula" |
                        genus=="Rubritalea" |
                        genus=="Synechococcus_CC9902")
View(Genus.Age)

Afipia.age <- subset(Genus.Age, genus=="Afipia")
Afipia.aov.age <- aov(RelativeAbundance~Age, data=Afipia.age)
summary(Afipia.aov.age)

Amyl.age <- subset(Genus.Age, genus=="Amylibacter")
Amyl.aov.age <- aov(RelativeAbundance~Age, data=Amyl.age)
summary(Amyl.aov.age)

Blast.age <- subset(Genus.Age, genus=="Blastopirellula")
Blast.aov.age <- aov(RelativeAbundance~Age, data=Blast.age)
summary(Blast.aov.age)

Croc.age <- subset(Genus.Age, genus=="Croceitalea")
Croc.aov.age <- aov(RelativeAbundance~Age, data=Croc.age)
summary(Croc.aov.age)

Dok.age <- subset(Genus.Age, genus=="Dokdonia")
Dok.aov.age <- aov(RelativeAbundance~Age, data=Dok.age)
summary(Dok.aov.age)

Flav.age <- subset(Genus.Age, genus=="Flavicella")
Flav.aov.age <- aov(RelativeAbundance~Age, data=Flav.age)
summary(Flav.aov.age)

Gran.age <- subset(Genus.Age, genus=="Granulosicoccus")
Gran.aov.age <- aov(RelativeAbundance~Age, data=Gran.age)
summary(Gran.aov.age)

Leuc.age <- subset(Genus.Age, genus=="Leucothrix")
Leuc.aov.age <- aov(RelativeAbundance~Age, data=Leuc.age)
summary(Leuc.aov.age)

Lew.age <- subset(Genus.Age, genus=="Lewinella")
Lew.aov.age <- aov(RelativeAbundance~Age, data=Lew.age)
summary(Lew.aov.age)

OM.age <- subset(Genus.Age, genus=="OM60(NOR5)_clade")
OM.aov.age <- aov(RelativeAbundance~Age, data=OM.age)
summary(OM.aov.age)

Phyc.age <- subset(Genus.Age, genus=="Phycisphaera")
Phyc.aov.age <- aov(RelativeAbundance~Age, data=Phyc.age)
summary(Phyc.aov.age)

Prop.age <- subset(Genus.Age, genus=="Propionigenium")
View(Prop.age)
Prop.aov.age <- aov(RelativeAbundance~Age, data=Prop.age)
summary(Prop.aov.age)

Psych.age <- subset(Genus.Age, genus=="Psychrobacter")
Psych.aov.age <- aov(RelativeAbundance~Age, data=Psych.age)
summary(Psych.aov.age)

Rub.age <- subset(Genus.Age, genus=="Rubripirellula")
Rub.aov.age <- aov(RelativeAbundance~Age, data=Rub.age)
summary(Rub.aov.age)

Rubri.age <- subset(Genus.Age, genus=="Rubritalea")
Rubri.aov.age <- aov(RelativeAbundance~Age, data=Rubri.age)
summary(Rubri.aov.age)

Syn.age <- subset(Genus.Age, genus=="Synechococcus_CC9902")
Syn.aov.age <- aov(RelativeAbundance~Age, data=Syn.age)
summary(Syn.aov.age)




















#Corncob of Genera > 0.5%--------
corncob.dat <- read.csv("Corncob output for plotting genera > 0.5% RA.csv")
View(corncob.dat)
corn.depth.ra.abundance <- subset(corncob.dat, Analysis=="Depth" & Differential=="Abundance")
View(corn.depth.ra.abundance)

corn.depth.ra.abundance$Genus <- factor(corn.depth.ra.abundance$Genus, levels=c('Rubritalea',
                                                                                'Roseibacillus',
                                                                                'Persicirhabdus',
                                                                                'Rubripirellula',
                                                                                'Phycisphaera',
                                                                                'Psychrobacter',
                                                                                'OM60(NOR5)_clade',
                                                                                'Leucothrix',
                                                                                'Agaribacter',
                                                                                'Propionigenium',
                                                                                'Wenyingzhuangia',
                                                                                'Lewinella',
                                                                                'Dokdonia',
                                                                                'Algitalea',
                                                                                'Clade_Ia',
                                                                                'Amylibacter',
                                                                                'Afipia'))

ggplot(corn.depth.ra.abundance, aes(x=Average, y=Genus, color=Association)) +
  geom_point(size=3) +
  geom_pointrange(aes(xmin=Average-Stdev, xmax=Average+Stdev)) + theme_bw() +
  geom_vline(xintercept=c(0,0), linetype="dotted") +
  scale_color_manual(values=c("gray", "blue")) +
  xlab("Regression Coefficient") +
  theme(axis.text.y = element_text(face="italic")) +
  ggtitle("Differential abundance: surface vs subsurface") +
  theme(plot.title = element_text(hjust = 0.5)) + xlim(-1,3) +
  theme(text = element_text(size=20)) +
  annotate("text", x=3, y=17, label = "A", size=8)
ggsave("Diff abundance depth.tiff")

corn.depth.ra.variability <- subset(corncob.dat, Analysis=="Depth" & Differential=="Variability")
View(corn.depth.ra.variability)
corn.depth.ra.variability$Genus <- factor(corn.depth.ra.variability$Genus, levels=c('Rubritalea',
                                                                                'Roseibacillus',
                                                                                'Persicirhabdus',
                                                                                'Rubripirellula',
                                                                                'Phycisphaera',
                                                                                'Psychrobacter',
                                                                                'OM60(NOR5)_clade',
                                                                                'Leucothrix',
                                                                                'Agaribacter',
                                                                                'Propionigenium',
                                                                                'Wenyingzhuangia',
                                                                                'Lewinella',
                                                                                'Dokdonia',
                                                                                'Algitalea',
                                                                                'Clade_Ia',
                                                                                'Amylibacter',
                                                                                'Afipia'))
ggplot(corn.depth.ra.variability, aes(x=Average, y=Genus, color=Association)) +
  geom_point(size=3) +
  geom_pointrange(aes(xmin=Average-Stdev, xmax=Average+Stdev)) + theme_bw() +
  geom_vline(xintercept=c(0,0), linetype="dotted") +
  scale_color_manual(values=c("gray", "blue")) +
  xlab("Regression Coefficient") +
  theme(axis.text.y = element_text(face="italic")) +
  ggtitle("Differential variability: surface vs subsurface") +
  theme(plot.title = element_text(hjust = 0.5)) + xlim(-2,2.2) +
  theme(text = element_text(size=20)) +
  annotate("text", x=2, y=17, label = "B", size=8)
ggsave("Diff variability depth.tiff")

corn.age.ra.abundance <- subset(corncob.dat, Analysis=="Age" & Differential=="Abundance")
View(corn.age.ra.abundance)
corn.age.ra.abundance$Genus <- factor(corn.age.ra.abundance$Genus, levels=c('Persicirhabdus',
                                                                            'Rubripirellula',
                                                                            'Rhodopirellula',
                                                                            'Phycisphaera',
                                                                            'Blastopirellula',
                                                                            'Psychrobacter',
                                                                            'OM60(NOR5)_clade',
                                                                            'Leucothrix',
                                                                            'Granulosicoccus',
                                                                            'Aliivibrio',
                                                                            'Agaribacter',
                                                                                    'Propionigenium',
                                                                                    'Synechococcus_CC902',
                                                                                    'NS4_marine_group',
                                                                                    'Lewinella',
                                                                                    'Flavicella',
                                                                                    'Dokdonia',
                                                                                    'Croceitalea',
                                                                                    'Brucella',
                                                                                    'Amylibacter',
                                                                                    'Afipia'))
ggplot(corn.age.ra.abundance, aes(x=Average, y=Genus, color=Association)) +
  geom_point(size=3) +
  geom_pointrange(aes(xmin=Average-Stdev, xmax=Average+Stdev)) + theme_bw() +
  geom_vline(xintercept=c(0,0), linetype="dotted") +
  scale_color_manual(values=c("gray", "blue")) +
  xlab("Regression Coefficient") +
  theme(axis.text.y = element_text(face="italic")) +
  ggtitle("Differential abundance: juvenile vs mature") +
  theme(plot.title = element_text(hjust = 0.5)) + xlim(-1,3) +
  theme(text = element_text(size=20)) +
  annotate("text", x=3, y=21, label = "C", size=8)
ggsave("Diff abundance age.tiff")

corn.age.ra.variability <- subset(corncob.dat, Analysis=="Age" & Differential=="Variability")
corn.age.ra.variability$Genus <- factor(corn.age.ra.variability$Genus, levels=c('Persicirhabdus',
                                                                            'Rubripirellula',
                                                                            'Rhodopirellula',
                                                                            'Phycisphaera',
                                                                            'Blastopirellula',
                                                                            'Psychrobacter',
                                                                            'OM60(NOR5)_clade',
                                                                            'Leucothrix',
                                                                            'Granulosicoccus',
                                                                            'Aliivibrio',
                                                                            'Agaribacter',
                                                                            'Propionigenium',
                                                                            'Synechococcus_CC902',
                                                                            'NS4_marine_group',
                                                                            'Lewinella',
                                                                            'Flavicella',
                                                                            'Dokdonia',
                                                                            'Croceitalea',
                                                                            'Brucella',
                                                                            'Amylibacter',
                                                                            'Afipia'))
ggplot(corn.age.ra.variability, aes(x=Average, y=Genus, color=Association)) +
  geom_point(size=3) +
  geom_pointrange(aes(xmin=Average-Stdev, xmax=Average+Stdev)) + theme_bw() +
  geom_vline(xintercept=c(0,0), linetype="dotted") +
  scale_color_manual(values=c("gray", "blue")) +
  xlab("Regression Coefficient") +
  theme(axis.text.y = element_text(face="italic")) +
  ggtitle("Differential variability: juvenile vs mature") +
  theme(plot.title = element_text(hjust = 0.5)) + xlim(-2,2.2) +
  theme(text = element_text(size=20)) +
  annotate("text", x=2, y=21, label = "D", size=8)
ggsave("Diff variability age.tiff")
#Barplots of Genera from Corncob >0.5% RA-------
Depth.bar <- read.csv("Depth Specific Barplot based on Corncob.csv")
plot.Depth.bar <- subset(Depth.bar, genus=="Amylibacter" |
                           genus=="Clade_Ia" |
                           genus=="Algitalea" |
                           genus=="Dokdonia" |
                           genus=="Lewinella" |
                           genus=="Propionigenium" |
                           genus=="OM60(NOR5)_clade" |
                           genus=="Psychrobacter" |
                           genus=="Rubripirellula" |
                           genus=="Roseibacillus")
plot.Depth.bar$genus <- factor(plot.Depth.bar$genus, levels=c('Roseibacillus',
                                                              'Rubripirellula',
                                                              'Psychrobacter',
                                                              'OM60(NOR5)_clade',
                                                              'Propionigenium',
                                                              'Lewinella',
                                                              'Dokdonia',
                                                              'Algitalea',
                                                              'Clade_Ia',
                                                              'Amylibacter'))
View(plot.Depth.bar)
ggplot(plot.Depth.bar, aes(x=average, y=genus, fill=depth)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values=c("purple2", "forestgreen"))+
  geom_errorbar(aes(xmin=average, xmax=average+stdev), width=.2,
                position=position_dodge(.9)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  xlab("% Relative Abundance") + ylab("Genus") +
  theme(axis.text.y = element_text(face="italic")) +
  annotate("text", x=1.7, y=10, label="*", size=12, color="forestgreen") +
  annotate("text", x=1.7, y=9, label="***", size=12, color="purple2") +
  annotate("text", x=1.7, y=8, label="***", size=12, color="forestgreen") +
  annotate("text", x=1.7, y=7, label="*", size=12, color="purple2") +
  annotate("text", x=1.7, y=6, label="**", size=12, color="purple2") +
  annotate("text", x=2.7, y=5, label="***", size=12, color="purple2") +
  annotate("text", x=1.7, y=4, label="***", size=12, color="purple2") +
  annotate("text", x=2.7, y=3, label="*", size=12, color="forestgreen") +
  annotate("text", x=3.8, y=2, label="*", size=12, color="forestgreen") +
  annotate("text", x=6, y=1, label="*", size=12, color="forestgreen") +
  annotate("text", x=15, y=10, label="A", size=12, color="black") +
  theme(text = element_text(size=20))
ggsave("Barplot post corncob depth.tiff")

Age.bar <- read.csv("Age specific barplot based on Corncob.csv")
View(Age.bar)
plot.Age.bar <- subset(Age.bar, genus=="Brucella" |
                         genus=="Dokdonia" |
                         genus=="Flavicella" |
                         genus=="Lewinella" |
                         genus=="NS4_marine_group" |
                         genus=="Synechococcus_CC902" |
                         genus=="Propionigenium" |
                         genus=="Leucothrix" |
                         genus=="OM60(NOR5)_clade" |
                         genus=="Psychrobacter" |
                         genus=="Blastopirellula" |
                         genus=="Persicirhabdus")
plot.Age.bar$genus <- factor(plot.Age.bar$genus, levels=c('Persicirhabdus',
                                                          'Blastopirellula',
                                                          'Psychrobacter',
                                                          'OM60(NOR5)_clade',
                                                          'Leucothrix',
                                                          'Propionigenium',
                                                          'Synechococcus_CC902',
                                                          'NS4_marine_group',
                                                          'Lewinella',
                                                          'Flavicella',
                                                          'Dokdonia',
                                                          'Brucella'))

View(plot.Age.bar)
ggplot(plot.Age.bar, aes(x=average, y=genus, fill=age)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(xmin=average, xmax=average+stdev), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values=c("goldenrod3", "turquoise4")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  xlab("% Relative Abundance") + ylab("Genus") +
  theme(axis.text.y = element_text(face="italic")) +
  annotate("text", x=3, y=11, label="***", size=12, color="turquoise4") +
  annotate("text", x=3, y=10, label="*", size=12, color="turquoise4") +
  annotate("text", x=6, y=9, label="**", size=12, color="goldenrod3") +
  annotate("text", x=3, y=8, label="**", size=12, color="turquoise4") +
  annotate("text", x=3, y=7, label="*", size=12, color="goldenrod3") +
  annotate("text", x=6, y=6, label="***", size=12, color="turquoise4") +
  annotate("text", x=6, y=5, label="*", size=12, color="turquoise4") +
  annotate("text", x=3, y=4, label="*", size=12, color="turquoise4") +
  annotate("text", x=3, y=3, label="*", size=12, color="turquoise4") +
  annotate("text", x=16, y=2, label="*", size=12, color="goldenrod3") +
  annotate("text", x=6, y=1, label="*", size=12, color="goldenrod3") +
  annotate("text", x=15, y=11, label="B", size=12, color="black") +
  theme(text = element_text(size=20))
ggsave("barplot post corncob age.tiff")

##ASV Specific ANOVA for biplot--------
ASV.biplot <- read.csv("ASV for biplot anova.csv")
View(ASV.biplot)
ASV148 <- ASV.biplot[ , c(1:4)]    
View(ASV148)
ASV148.aov <- aov(ASV148~Age, data=ASV148)
summary(ASV148.aov)

ASV61 <- ASV.biplot[ , c(1:3,5)]    
View(ASV61)
ASV61.aov <- aov(ASV61~Depth, data=ASV61)
summary(ASV61.aov)

ASV49 <- ASV.biplot[ , c(1:3,6)]    
ASV49.aov.depth <- aov(ASV49~Depth, data=ASV49)
ASV49.aov.age <- aov(ASV49~Age, data=ASV49)
summary(ASV49.aov.age)
summary(ASV49.aov.depth) #Dokdonia depth only story has meaning

ASV101 <- ASV.biplot[ , c(1:3,7)]    
View(ASV101)
ASV101.aov.depth <- aov(ASV101~Depth, data=ASV101)
ASV101.aov.age <- aov(ASV101~Age, data=ASV101)
summary(ASV101.aov.age)
summary(ASV101.aov.depth)

ASV98 <- ASV.biplot[,c(1:3,8)]
View(ASV98)
ASV98.aov.depth <- aov(Propionigenium_98~Depth, data=ASV98)
ASV98.aov.age <- aov(Propionigenium_98~Age, data=ASV98)
summary(ASV98.aov.age)
summary(ASV98.aov.depth)

ASV32 <- ASV.biplot[,c(1:3,9)]
View(ASV32)
ASV32.aov.depth <- aov(Leucothrix_32~Depth, data=ASV32)
ASV32.aov.age <- aov(Leucothrix_32~Age, data=ASV32)
summary(ASV32.aov.age)
summary(ASV32.aov.depth)

ASV33 <- ASV.biplot[,c(1:3,10)]
View(ASV33)
ASV33.aov.depth <- aov(Leucothrix_33~Depth, data=ASV33)
ASV33.aov.age <- aov(Leucothrix_33~Age, data=ASV33)
summary(ASV33.aov.age)
summary(ASV33.aov.depth)

ASV142 <- ASV.biplot[,c(1:3,11)]
View(ASV142)
ASV142.aov.depth <- aov(NOR5_142~Depth, data=ASV142)
summary(ASV142.aov.depth)

ASV525 <- ASV.biplot[,c(1:3,12)]
View(ASV525)
ASV525.aov.depth <- aov(NOR5_525~Depth, data=ASV525)
summary(ASV525.aov.depth)

ASV86 <- ASV.biplot[,c(1:3,13)]
View(ASV86)
ASV86.aov.depth <- aov(NOR5_86~Depth, data=ASV86)
summary(ASV86.aov.depth)

ASV110 <- ASV.biplot[,c(1:3,14)]
View(ASV110)
ASV110.aov.depth <- aov(Psychrobacter_110~Depth, data=ASV110)
ASV110.aov.age <- aov(Psychrobacter_110~Age, data=ASV110)
summary(ASV110.aov.age)
summary(ASV110.aov.depth)

ASV91 <- ASV.biplot[,c(1:3,15)]
View(ASV91)
ASV91.aov.depth <- aov(Psychrobacter_91~Depth, data=ASV91)
ASV91.aov.age <- aov(Psychrobacter_91~Age, data=ASV91)
summary(ASV91.aov.age)
summary(ASV91.aov.depth)

ASV95 <- ASV.biplot[,c(1:3,16)]
View(ASV95)
ASV95.aov.depth <- aov(Psychrobacter_95~Depth, data=ASV95)
ASV95.aov.age <- aov(Psychrobacter_95~Age, data=ASV95)
summary(ASV95.aov.age)
summary(ASV95.aov.depth)

ASV10 <- ASV.biplot[,c(1:3,17)]
View(ASV10)
ASV10.aov.age <- aov(Blastopirellula_10~Age, data=ASV10)
summary(ASV10.aov.age)

ASV11 <- ASV.biplot[,c(1:3,18)]
View(ASV11)
ASV11.aov.age <- aov(Blastopirellula_11~Age, data=ASV11)
summary(ASV11.aov.age)

ASV20 <- ASV.biplot[,c(1:3,19)]
View(ASV20)
ASV20.aov.age <- aov(Blastopirellula_20~Age, data=ASV20)
summary(ASV20.aov.age)

ASV52 <- ASV.biplot[,c(1:3,20)]
View(ASV52)
ASV52.aov.age <- aov(Blastopirellula_52~Age, data=ASV52)
summary(ASV52.aov.age)

ASV8 <- ASV.biplot[,c(1:3,21)]
View(ASV8)
ASV8.aov.age <- aov(Blastopirellula_8~Age, data=ASV8)
summary(ASV8.aov.age)

ASV109 <- ASV.biplot[,c(1:3,22)]
View(ASV109)
ASV109.aov.depth <- aov(Rubripirellula_109~Depth, data=ASV109)
summary(ASV109.aov.depth)

ASV120 <- ASV.biplot[,c(1:3,23)]
View(ASV120)
ASV120.aov.depth <- aov(Rubripirellula_120~Depth, data=ASV120)
summary(ASV120.aov.depth)

ASV28 <- ASV.biplot[,c(1:3,24)]
View(ASV28)
ASV28.aov.depth <- aov(Rubripirellula_28~Depth, data=ASV28)
summary(ASV28.aov.depth)

ASV24 <- ASV.biplot[,c(1:3,25)]
View(ASV24)
ASV24.aov.age <- aov(Persicirhabdus_24~Age, data=ASV24)
summary(ASV24.aov.age)

ASV26 <- ASV.biplot[,c(1:3,26)]
View(ASV26)
ASV26.aov.age <- aov(Persicirhabdus_26~Age, data=ASV26)
summary(ASV26.aov.age)

ASV59 <- ASV.biplot[,c(1:3,27)]
View(ASV59)
ASV59.aov.age <- aov(Persicirhabdus_59~Age, data=ASV59)
summary(ASV59.aov.age)

ASV82 <- ASV.biplot[,c(1:3,28)]
View(ASV82)
ASV82.aov.age <- aov(Persicirhabdus_82~Age, data=ASV82)
summary(ASV82.aov.age)

ASV19 <- ASV.biplot[,c(1:3,29)]
View(ASV19)
ASV19.aov.depth <- aov(Roseibacillus_19~Depth, data=ASV19)
summary(ASV19.aov.depth)

ASV43 <- ASV.biplot[,c(1:3,30)]
View(ASV43)
ASV43.aov.depth <- aov(Roseibacillus_43~Depth, data=ASV43)
summary(ASV43.aov.depth)

ASV67 <- ASV.biplot[,c(1:3,31)]
View(ASV67)
ASV67.aov.depth <- aov(Roseibacillus_67~Depth, data=ASV67)
summary(ASV67.aov.depth)

ASV73 <- ASV.biplot[,c(1:3,32)]
View(ASV73)
ASV73.aov.depth <- aov(Roseibacillus_73~Depth, data=ASV73)
summary(ASV73.aov.depth)
##
##ASV CCA-----
ASV.dat <- read.csv("ASV for biplot anova.csv")
View(ASV.dat)
ASV.dat %>% remove_rownames %>% column_to_rownames(var="UniqueID") 
ASV.dat <- column_to_rownames(ASV.dat, "UniqueID")
View(ASV.dat)

RA.dat <- read.csv("Relative Abundance for ASV CCA biplot.csv")
RA.dat %>% remove_rownames %>% column_to_rownames(var="UniqueID") 
RA.dat <- column_to_rownames(RA.dat, "UniqueID")
View(RA.dat)

ASV.cca <- cca(RA.dat ~ Brucella_148+
                 Amylibacter_61+
                 Dokdonia_49+
                 Lewinella_101+
                 Propionigenium_98+
                 Leucothrix_32+
                 Leucothrix_33+
                 NOR5_142+
                 NOR5_525+
                 NOR5_86+
                 Psychrobacter_110+
                 Psychrobacter_91+
                 Psychrobacter_95+
                 Blastopirellula_10+
                 Blastopirellula_11+
                 Blastopirellula_20+
                 Blastopirellula_20+
                 Blastopirellula_8+
                 Rubripirellula_109+
                 Rubripirellula_120+
                 Rubripirellula_28+
                 Persicirhabdus_24+
                 Persicirhabdus_26+
                 Persicirhabdus_59+
                 Persicirhabdus_82+
                 Roseibacillus_19+
                 Roseibacillus_43+
                 Roseibacillus_67+
                 Roseibacillus_73 , data=ASV.dat)
plot(ASV.cca)
summary(ASV.cca)
View(ASV.cca)

ASV.cca.coords <- read.csv("CCA coordinates ASV segments.csv")
View(ASV.cca.coords)
ASV.cca.coords.depth <- subset(ASV.cca.coords, ASV=="Amylibacter_61"|
                                 ASV=="Dokdonia_49" |
                                 ASV=="Lewinella_101"|
                                 ASV=="Propionigenium_98"|
                                 ASV=="NOR5_142"|
                                 ASV=="NOR5_525"|
                                 ASV=="NOR5_86" |
                                 ASV=="Psychrobacter_110" |
                                 ASV=="Psychrobacter_91" |
                                 ASV=="Psychrobacter_95" |
                                 ASV=="Rubripirellula_109" |
                                 ASV=="Rubripirellula_120" |
                                 ASV=="Rubripirellula_28" |
                                 ASV=="Roseibacillus_19" |
                                 ASV=="Roseibacillus_43" |
                                 ASV=="Roseibacillus_67" |
                                 ASV=="Roseibacillus_73")
View(ASV.cca.coords.depth)
UI.cca.coords <- read.csv("UniqueID Coordinate for ASV biplot CCA.csv")
View(UI.cca.coords)

ggplot(UI.cca.coords) +
  geom_point(mapping = aes(x=CCA1, y=CCA2, color=Depth, shape=Age), size=4, alpha=0.6) +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  scale_color_manual(values = c('purple2', 'forestgreen')) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18))+
  geom_segment(data=ASV.cca.coords.depth,
               aes(x=0, xend=CCA1, y=0, yend=CCA2),
               type="closed", color="black") +
  geom_text(data = ASV.cca.coords.depth, aes(x = CCA1, y = CCA2, label = ASV),
            size = 4) +
  xlab("CCA1 (33.9%)") + ylab("CCA2 (30.6%)")

##ASV NMDS Biplot
UI.NMDS.coords <- read.csv("UniqueID NMDS coordinates for biplot.csv")
View(UI.NMDS.coords)
ASV.NMDS.coords <- read.csv("NMDS biplot Age ASV coordinates.csv")
View(ASV.NMDS.coords)
ggplot(UI.NMDS.coords) +
  geom_point(mapping = aes(x=NMDS1, y=NMDS2, color=Age, shape=Depth), size=6, alpha=0.4) +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  scale_color_manual(values = c('goldenrod3', 'turquoise3')) +
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 20))+
  geom_segment(data=ASV.NMDS.coords,
               aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               type="closed", color="gray", alpha=0.75) +
  geom_text(data = ASV.NMDS.coords, aes(x = NMDS1, y = NMDS2, label = ASVLabel),
            size = 3) + xlim(-1.01, 1.55) +
  annotate("text", x=-1, y=1, label="D", size=12, color="black")
ggsave("NMDS.biplot.depth.tiff")

ASV.depth.NMDS.coords <- read.csv("Depth specific ASV NMDS biplot coordinates.csv")
ggplot(UI.NMDS.coords) +
  geom_point(mapping = aes(x=NMDS1, y=NMDS2, color=Depth, shape=Age), size=6, alpha=0.4) +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  scale_color_manual(values = c('purple2', 'forestgreen')) +
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 20))+
  geom_segment(data=ASV.depth.NMDS.coords,
               aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               type="closed", color="gray", alpha=0.75) +
  geom_text(data = ASV.depth.NMDS.coords, aes(x = NMDS1, y = NMDS2, label = ASVLabel),
            size = 3) + xlim(-1.1, 1.55)  +
  annotate("text", x=-1, y=1, label="C", size=12, color="black")
ggsave("NMDS.biplot.age.tiff")

#060723 Re-DoCalculations------------
#Creating photophsyiology plots---------
PhotoPhys.New <- read.csv("PhotophysiologyData 060723 Reworked.csv")
View(PhotoPhys.New)
ggplot(PhotoPhys.New, aes(x=SurfaceAreacm, y=Comparison, fill=Label)) +
  geom_boxplot(alpha=0.8, position=position_dodge()) +
  theme_bw() +
  coord_flip() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c('yellow3', 'darkblue', 'forestgreen', 'purple3')) +
  scale_color_manual(values = c('yellow3', 'darkblue', 'forestgreen', 'purple3')) +
  annotate("text", x=500, y=0.5, label = "A", size=12) +
  xlab(bquote('Surface Area '(cm^2))) + ylab("")

ggplot(PhotoPhys.New, aes(x=SurfaceAreacm, y=Comparison, fill=Label)) +
  geom_violin(position=position_dodge()) +
  geom_boxplot(binaxis='y', width=.05, 
               position=position_dodge(0.9), color="white", outlier.shape = NA) +
  theme_bw() +
  coord_flip() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c('yellow3', 'darkblue', 'forestgreen', 'purple3')) +
  scale_color_manual(values = c('yellow3', 'darkblue', 'forestgreen', 'purple3')) +
  annotate("text", x=500, y=0.5, label = "(A)", size=12) +
  xlab(bquote('Surface Area '(cm^2))) + ylab("")

SA.aov <- aov(SurfaceAreacm ~ Label, data = PhotoPhys.New)
summary(SA.aov)
TukeyHSD(SA.aov)

ggplot(PhotoPhys.New, aes(x=FvFm, y=Comparison, fill=Label)) +
  geom_violin(position=position_dodge()) +
  geom_boxplot(binaxis='y', width=.05, 
               position=position_dodge(0.9), color="white", outlier.shape = NA) +
  theme_bw() +
  coord_flip() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  scale_fill_manual(values = c('yellow3', 'darkblue', 'forestgreen', 'purple3')) +
  scale_color_manual(values = c('yellow3', 'darkblue', 'forestgreen', 'purple3')) +
  xlim(0.219, 0.72) +
  annotate("text", x=0.7, y=0.78, label = "A", size=10) +
  annotate("text", x=0.7, y=1.23, label = "B", size=10) +
  annotate("text", x=0.7, y=1.78, label = "A", size=10) +
  annotate("text", x=0.7, y=2.23, label = "B", size=10) +
  annotate("text", x=0.72, y=0.5, label = "(B)", size=12) +
  xlab("Maximum Potential Quantum Efficiency (Fv/Fm)") +
  ylab("")
  

FvFm.aov <- aov(FvFm ~ Label, data = PhotoPhys.New)
summary(FvFm.aov)
TukeyHSD(FvFm.aov)

ggplot(PhotoPhys.New, aes(x=ChlaC, y=Comparison, fill=Label)) +
  geom_violin(position=position_dodge(),) +
  geom_boxplot(binaxis='y', width=.05, 
               position=position_dodge(0.9), color="white", outlier.shape = NA) +
  theme_bw() +
  coord_flip() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  scale_fill_manual(values = c('yellow3', 'darkblue', 'forestgreen', 'purple3')) +
  scale_color_manual(values = c('yellow3', 'darkblue', 'forestgreen', 'purple3')) +
  xlim(0.001, 0.013) +
  annotate("text", x=0.012, y=0.79, label = "A", size=10) +
  annotate("text", x=0.012, y=1.22, label = "C", size=10) +
  annotate("text", x=0.012, y=1.78, label = "BC", size=10) +
  annotate("text", x=0.012, y=2.22, label = "AB", size=10) +
  annotate("text", x=0.013, y=0.5, label = "(C)", size=12) +
  xlab("Chla:C (mg)") + ylab("")
  
  ChlaC.aov <- aov(ChlaC ~ Label, data = PhotoPhys.New)
  summary(ChlaC.aov)
  TukeyHSD(ChlaC.aov)
  

  ggplot(PhotoPhys.New, aes(x=CN, y=Comparison, fill=Label)) +
    geom_violin(position=position_dodge()) +
    geom_boxplot(binaxis='y', width=.05, 
                 position=position_dodge(0.9), color="white", outlier.shape = NA) +
    theme_bw() +
    coord_flip() +
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text=element_text(size=20),
          legend.title = element_text(size=20),
          legend.position="none") +
    scale_fill_manual(values = c('yellow3', 'darkblue', 'forestgreen', 'purple3')) +
    scale_color_manual(values = c('yellow3', 'darkblue', 'forestgreen', 'purple3')) +
    xlim(6, 30) +
    annotate("text", x=29, y=0.77, label = "A", size=10) +
    annotate("text", x=29, y=1.23, label = "A", size=10) +
    annotate("text", x=29, y=1.77, label = "B", size=10) +
    annotate("text", x=29, y=2.23, label = "C", size=10) +
    annotate("text", x=30, y=0.53, label = "(D)", size=12) +
    xlab("C:N (%)") + ylab("")
  
  CN.aov <- aov(CN ~ Label, data = PhotoPhys.New)
  summary(CN.aov)
  TukeyHSD(CN.aov)
#Re-Do Diversity Metrics -----------
Rework.DivMet <- read.csv("Reworked-060723-DiversityMetrics-Barplot.csv")
View(Rework.DivMet)  
ggplot(Rework.DivMet, aes(x=AlphaDiversity, y=Comparison, fill=Label)) +
  geom_violin(position=position_dodge()) +
  geom_boxplot(binaxis='y', width=.05, 
               position=position_dodge(0.9), color="white", outlier.shape = NA) +
  theme_bw() +
  coord_flip() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=20),
        legend.title = element_text(size=20),
        legend.position="none") +
  scale_fill_manual(values = c('yellow3', 'darkblue', 'forestgreen', 'purple3')) +
  scale_color_manual(values = c('yellow3', 'darkblue', 'forestgreen', 'purple3')) +
  xlim(480, 2200) +
  xlab("ASV Richness") + ylab("") +
  annotate("text", x=2150, y=0.78, label = "A", size=10) +
  annotate("text", x=2150, y=1.23, label = "C", size=10) +
  annotate("text", x=2150, y=1.78, label = "B", size=10) +
  annotate("text", x=2150, y=2.23, label = "AB", size=10) +
  annotate("text", x=2200, y=0.5, label = "(A)", size=12) 
  

  ASV.aov <- aov(AlphaDiversity ~ Label, data = Rework.DivMet)
  summary(ASV.aov)
  TukeyHSD(ASV.aov)
  
  ggplot(Rework.DivMet, aes(x=PielouEvenness, y=Comparison, fill=Label)) +
    geom_violin(position=position_dodge()) +
    geom_boxplot(binaxis='y', width=.05, 
                 position=position_dodge(0.9), color="white", outlier.shape = NA) +
    theme_bw() +
    coord_flip() +
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text=element_text(size=20),
          legend.title = element_text(size=20),
          legend.position = "none") +
    scale_fill_manual(values = c('yellow3', 'darkblue', 'forestgreen', 'purple3')) +
    scale_color_manual(values = c('yellow3', 'darkblue', 'forestgreen', 'purple3')) +
    xlim(0.55, 0.85) +
    xlab("Pielou's Evenness") + ylab("") +
    annotate("text", x=0.84, y=0.79, label = "A", size=10) +
    annotate("text", x=0.84, y=1.23, label = "C", size=10) +
    annotate("text", x=0.84, y=1.79, label = "B", size=10) +
    annotate("text", x=0.84, y=2.23, label = "B", size=10) +
    annotate("text", x=0.85, y=0.5, label = "(B)", size=12) 
  
  PE.aov <- aov(PielouEvenness ~ Label, data = Rework.DivMet)
  summary(PE.aov)
  TukeyHSD(PE.aov)
  
 

  ggplot(Rework.DivMet, aes(x=BetaDisper, y=Comparison, fill=Label)) +
    geom_violin(position=position_dodge()) +
    geom_boxplot(binaxis='y', width=.05, 
                 position=position_dodge(0.9), color="white", outlier.shape = NA) +
    theme_bw() +
    coord_flip() +
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text=element_text(size=20),
          legend.title = element_text(size=20)) +
    scale_fill_manual(values = c('yellow3', 'darkblue', 'forestgreen', 'purple3')) +
    scale_color_manual(values = c('yellow3', 'darkblue', 'forestgreen', 'purple3')) +
    xlim(0.12, 0.6) +
    xlab("Beta Dispersion") + ylab("") +
    annotate("text", x=0.59, y=0.78, label = "A", size=10) +
    annotate("text", x=0.59, y=1.23, label = "B", size=10) +
    annotate("text", x=0.59, y=1.78, label = "C", size=10) +
    annotate("text", x=0.59, y=2.23, label = "BC", size=10) +
    annotate("text", x=0.6, y=0.52, label = "(C)", size=12) 
  
  BD.aov <- aov(BetaDisper ~ Label, data = Rework.DivMet)
  summary(BD.aov)
  TukeyHSD(BD.aov)
  
#050723 Re-Do NMDS ----------
  #Based on Bart's and Joey's suggestions

NMDS.coords <- read.csv("ReVisualizing NMDS based on Bart Suggestions.csv")
View(NMDS.coords)  
ggplot(NMDS.coords, aes(x=NMDS1, y=NMDS2, color=Age, shape=Age)) +
  geom_point(alpha=0.9, size=5) +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=20),
        legend.title = element_text(size=20)) +
  scale_color_manual(values=c("yellow3", "darkblue")) +
  stat_ellipse(geom = "polygon",
               aes(fill = Age),  
               alpha = 0.25) +
  scale_fill_manual(values=c("yellow3", "darkblue")) +
  xlim(-1, 1.1) + ylim(-1, 1.2) 
  

ggplot(NMDS.coords, aes(x=NMDS1, y=NMDS2, color=Depth, shape=Depth)) +
  geom_point(alpha=0.9, size=5) +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=20),
        legend.title = element_text(size=20)) +
  scale_color_manual(values=c("forestgreen", "purple2")) +
  stat_ellipse(geom = "polygon",
               aes(fill = Depth), 
               alpha = 0.25) +
  scale_fill_manual(values=c("forestgreen", "purple2")) +
  xlim(-1, 1.1) + ylim(-1, 1.2) +
  annotate("text", x=-1, y=1, label = "B", size=10)

#NMDS Just Mature ----------
mature.blade.rarreads <- read.csv("062323 - Mature Only NMDS.csv")
View(mature.blade.rarreads)
mature.blade.rarreads.tidy <- mature.blade.rarreads %>%
  tibble::column_to_rownames("UniqueID")
BC.NMDS.shell<- function(dataframe) {
  metaMDS(dataframe, distance="bray", k=2, trymax=1000)
} 
mature.BC.NMDS <- BC.NMDS.shell(mature.blade.rarreads.tidy)  
View(mature.BC.NMDS$points)
View(mature.BC.NMDS$species)
#write.csv(mature.BC.NMDS$points, "NMDS Only Mature points ordination.csv")
#write.csv(mature.BC.NMDS$species, "NMDS Only Mature ASV coordinates for biplot.csv")
Mat.Samp.coords <- read.csv("NMDS Only Mature points ordination.csv")
Mat.ASV.coords <- read.csv("NMDS ASV coords mature only for biplot.csv")
View(Mat.ASV.coords)
ggplot(Mat.Samp.coords) +
  geom_point(mapping = aes(x=NMDS1, y=NMDS2, shape=Depth), color="darkblue", size=6, alpha=0.4) +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text=element_text(size=20),
        legend.title = element_text(size=20)) +
  geom_segment(data=Mat.ASV.coords,
               aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               type="closed", color="gray", alpha=0.5) +
  geom_text(data = Mat.ASV.coords, aes(x = NMDS1, y = NMDS2, label = Label),
            size = 4) + xlim(-1.2, 1)
