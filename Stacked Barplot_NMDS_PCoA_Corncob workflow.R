##The purpose of this R Script is to setup all of the figures necessary for the meeting with Lizzy & Emily
##First - let's make some stacked barplots --------
PRA <- read.csv("Percent Relative Abundance - Microbiome Only.csv")
View(PRA)
FamilyASV <- read.csv("ASV & Family Coordinates.csv")
View(FamilyASV)

PRA.FamilyASV <- merge(FamilyASV, PRA, by = "AmpliconSequenceVariant") #this worked but now we need to get rid of the annoying NAs
View(PRA.FamilyASV) 
####write.csv(PRA.FamilyASV, "Percent Relative Abundance of Bacteria Families - Microbiome Only.csv")

##Now I want to sum all Families with the same name
library(tidyr)
library(dplyr)

PR.Family <- PRA.FamilyASV[, -c(1:2)] #get rid of ASV and phylum columns
View(PR.Family)

PR.Family.gather <- tidyr::gather(PR.Family, "UniqueID", "RelativeAbundance", 2:87) #perfect gather
View(PR.Family.gather)

#Now, we want to calculate averages and use those to determine which will stay in the stacked barplot and which will be other
#To do that, we need to merge a sample coords file
Sample.Coords <- read.csv("Sample Coordinates - Microbiome Only.csv")
View(Sample.Coords)
PRA.Family.gather.samps <- merge(Sample.Coords, PR.Family.gather, by = "UniqueID")
View(PRA.Family.gather.samps)


#Sum Families per age per depth
PRA.Family.Sum.Samp <- PRA.Family.gather.samps %>% group_by(Family, UniqueID, Depth, Age, FrondandBlade) %>% 
  summarise(across(c(RelativeAbundance),sum),
            .groups = 'drop') %>% as.data.frame()  #sum calculations
View(PRA.Family.Sum.Samp)

PRA.Family.Avg.Samp <- PRA.Family.Sum.Samp %>% group_by(Family, Depth, Age, FrondandBlade) %>%
  summarise(across(c(RelativeAbundance),mean),
            .groups = 'drop') %>% as.data.frame() #average calculations
View(PRA.Family.Avg.Samp)
#write.csv(PRA.Family.Avg.Samp, "Percent Relative Abundance Family Averages per category - for observing.csv")
#Now, I need to analyze this in excel to figure out which ones have % RA > 1
#Then, I'll set up the stacked barplots.
#The following bacterial families have %RA > 0.5%
#Microtrichaceae, Cryomorphaceae, Cyclobacteriaceae, Flammeovirgaceae, Flavobacteriaceae,
#Saprospiraceae, Bacteriovoracaceae, Bdellovibrioceae, JG30-KF-CM45, Cyanobiaceae, Marinomodaceae,
#Fusobacteriaceae, Phycisphaeraceae, Pirellulaceae, Rubinisphaeraceae, Alcanivoracaceae, Alteromodaceae
#Cellvibrioceae, Clade_I, Colwelliaceae, Granulosicoccaceae, Halieaceae, Halomodaceae, Hyphomodaceae,
#Moraxellaceae, Porticoccaceae, Pseudoalteromodaceae, Psychromodaceae, Rhizobiaceae, Rhodobacteraceae, Sphingomodaceae
#Terasakiellaceae, Thiotrichaceae, Vibrioceae, Xanthobacteraceae, Xanthomodaceae, DEV007, Kiritimatiellaceae, Puniceicoccaceae, Rubritaleaceae


Selected.PRA.Family.Sum.samps <- subset(PRA.Family.Sum.Samp, Family %in% 
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
                                               'Alcanivoracaceae1',
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
                                               'Families < 0.5%'))
View(Selected.PRA.Family.Sum.samps)
####write.csv(Selected.PRA.Family.Sum.samps, "Families RA > 0.5%.csv") 

RA.sbp <- read.csv("Families RA > 0.5%.csv")
View(RA.sbp)
RA.sbp.juvenile.surface <- subset(RA.sbp, Depth=="surface" & Age=="Juvenile")
View(RA.sbp.juvenile.surface)
RA.sbp.juvenile.surface$UniqueID<- as.character(RA.sbp.juvenile.surface$UniqueID)
RA.sbp.juvenile.surface$Family <- factor(RA.sbp.juvenile.surface$Family, levels=c('Microtrichaceae',
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
                                                                           'OtherFamilies'))

library(ggplot2)
ggplot(RA.sbp.juvenile.surface, aes(fill=Family, y=RelativeAbundance, x=FrondandBlade)) +
  geom_bar(position="stack", stat="identity", color="black") +
  theme_bw() + xlab("Juvenile Surface Samples") + ylab("% Relative Abundance") +
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
                               "bisque3", "black", "white"))

RA.sbp.juvenile.middle <- subset(RA.sbp, Depth=="middle" & Age=="Juvenile")
RA.sbp.juvenile.middle$UniqueID<- as.character(RA.sbp.juvenile.middle$UniqueID)
RA.sbp.juvenile.middle$Family <- factor(RA.sbp.juvenile.middle$Family, levels=c('Microtrichaceae',
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
                                                                                  'OtherFamilies'))

View(RA.sbp.juvenile.middle)
ggplot(RA.sbp.juvenile.middle, aes(fill=Family, y=RelativeAbundance, x=FrondandBlade)) +
  geom_bar(position="stack", stat="identity", color="black") +
  theme_bw() + xlab("Juvenile Middle Samples") + ylab("% Relative Abundance") +
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
                               "bisque3", "black", "white"))

RA.sbp.juvenile.bottom <- subset(RA.sbp, Depth=="bottom" & Age=="Juvenile")
RA.sbp.juvenile.bottom$UniqueID<- as.character(RA.sbp.juvenile.bottom$UniqueID)
RA.sbp.juvenile.bottom$Family <- factor(RA.sbp.juvenile.bottom$Family, levels=c('Microtrichaceae',
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
                                                                                'OtherFamilies'))

ggplot(RA.sbp.juvenile.bottom, aes(fill=Family, y=RelativeAbundance, x=FrondandBlade)) +
  geom_bar(position="stack", stat="identity", color="black") +
  theme_bw() + xlab("Juvenile Bottom Samples") + ylab("% Relative Abundance") +
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
                               "bisque3", "black", "white"))

RA.sbp.mature.bottom <- subset(RA.sbp, Depth=="bottom" & Age=="Mature")
RA.sbp.mature.bottom$UniqueID<- as.character(RA.sbp.mature.bottom$UniqueID)
RA.sbp.mature.bottom$Family <- factor(RA.sbp.mature.bottom$Family, levels=c('Microtrichaceae',
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
                                                                                'OtherFamilies'))

ggplot(RA.sbp.mature.bottom, aes(fill=Family, y=RelativeAbundance, x=FrondandBlade)) +
  geom_bar(position="stack", stat="identity", color="black") +
  theme_bw() + xlab("Mature Bottom Samples") + ylab("% Relative Abundance") +
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
                               "bisque3", "black", "white"))

RA.sbp.mature.middle <- subset(RA.sbp, Depth=="middle" & Age=="Mature")
RA.sbp.mature.middle$UniqueID<- as.character(RA.sbp.mature.middle$UniqueID)
RA.sbp.mature.middle$Family <- factor(RA.sbp.mature.middle$Family, levels=c('Microtrichaceae',
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
                                                                            'OtherFamilies'))

ggplot(RA.sbp.mature.middle, aes(fill=Family, y=RelativeAbundance, x=FrondandBlade)) +
  geom_bar(position="stack", stat="identity", color="black") +
  theme_bw() + xlab("Mature Middle Samples") + ylab("% Relative Abundance") +
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
                               "bisque3", "black", "white"))

RA.sbp.mature.surface <- subset(RA.sbp, Depth=="surface" & Age=="Mature")
View(RA.sbp.mature.surface)
RA.sbp.mature.surface$UniqueID<- as.character(RA.sbp.mature.surface$UniqueID)
RA.sbp.mature.surface$Family <- factor(RA.sbp.mature.surface$Family, levels=c('Microtrichaceae',
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
                                                                            'OtherFamilies'))

ggplot(RA.sbp.mature.surface, aes(fill=Family, y=RelativeAbundance, x=FrondandBlade)) +
  geom_bar(position="stack", stat="identity", color="black") +
  theme_bw() + xlab("Mature Surface Samples") + ylab("% Relative Abundance") +
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
                               "bisque3", "black", "white"))

#Genus level stack barplots -------
RA <- read.csv("Percent Relative Abundance - Microbiome Only.csv")
View(PRA)
GenusASV <- read.csv("Genus.Species.Coordinates.csv")
View(GenusASV)

PRA.Genus.ASV <- merge(GenusASV, PRA, by = "AmpliconSequenceVariant") #this worked but now we need to get rid of the annoying NAs
View(PRA.Genus.ASV) 
##write.csv(PRA.Genus.ASV, "Percent Relative Abundance of Bacteria Genera - Microbiome Only.csv")

PR.Genus <- PRA.Genus.ASV[, -c(1:2)] #get rid of ASV and phylum columns
View(PR.Genus)

PR.Genus.gather <- tidyr::gather(PR.Genus, "UniqueID", "RelativeAbundance", 2:87) #perfect gather
View(PR.Genus.gather)

#Now, we want to calculate averages and use those to determine which will stay in the stacked barplot and which will be other
#To do that, we need to merge a sample coords file
Sample.Coords <- read.csv("Sample Coordinates - Microbiome Only.csv")
View(Sample.Coords)
PRA.Genus.gather.samps <- merge(Sample.Coords, PR.Genus.gather, by = "UniqueID")
View(PRA.Genus.gather.samps)

#Sum Genera per age per depth
library(dplyr)
PRA.Genus.Sum.Samp <- PRA.Genus.gather.samps %>% group_by(genus, UniqueID, Depth, Age, FrondandBlade) %>% 
  summarise(across(c(RelativeAbundance),sum),
            .groups = 'drop') %>% as.data.frame()  #sum calculations
View(PRA.Genus.Sum.Samp)

PRA.Genus.Avg.Samp <- PRA.Genus.Sum.Samp %>% group_by(genus, Depth, Age, UniqueID, FrondandBlade) %>%
  summarise(across(c(RelativeAbundance),mean),
            .groups = 'drop') %>% as.data.frame() #average calculations
View(PRA.Genus.Avg.Samp)

##average by depth and age - no unique id
Genus.Avg <- PRA.Genus.Sum.Samp %>% group_by(genus, Depth, Age) %>%
  summarise(across(c(RelativeAbundance),mean),
            .groups = 'drop') %>% as.data.frame() #average calculations

##write.csv(Genus.Avg, "Percent RA Genera Averages by depth and age >1%.csv")
###write.csv(PRA.Genus.Avg.Samp, "Percent Relative Abundance Genera Averages per category - for observing.csv")
library(tidyr)
PRA.Genus.Sum.Samp.sub <- subset(PRA.Genus.Sum.Samp, genus %in%
                                   c('Afipia',
                                     'Amylibacter',
                                     'Blastopirellula',
                                     'Croceitalea',
                                     'Dokdonia',
                                     'Flavicella',
                                     'Granulosicoccus',
                                     'Halomos',
                                     'Hellea',
                                     'Lentimos',
                                     'Leucothrix',
                                     'Lewinella',
                                     'Litorimos',
                                     'OM60(NOR5)_clade',
                                     'Persicirhabdus',
                                     'Phycisphaera',
                                     'Propionigenium',
                                     'Pseudoalteromos',
                                     'Psychrobacter',
                                     'Psychromos',
                                     'Roseibacillus',
                                     'Rubidimos',
                                     'Rubripirellula',
                                     'Rubritalea',
                                     'Synechococcus_CC9902',
                                     'Wenyingzhuangia'))

View(PRA.Genus.Sum.Samp.sub)
##write.csv(PRA.Genus.Sum.Samp.sub, "Genus subset greater than 1% RA.csv")
Genus.SBP <- read.csv("Genus subset greater than 1% RA.csv")
View(Genus.SBP)
#genus surface juvenile blades
Genus.SBP.surface.juv <- subset(Genus.SBP, Depth=="surface" & Age=="Juvenile")
View(Genus.SBP.surface.juv)
Genus.SBP.surface.juv$FrondandBlade<- as.character(Genus.SBP.surface.juv$FrondandBlade)
Genus.SBP.surface.juv$genus <- factor(Genus.SBP.surface.juv$genus, levels=c('Lewinella',
                                                                            'Rubidimos',
                                                                            'Croceitalea',
                                                                            'Dokdonia',
                                                                            'Flavicella',
                                                                            'Wenyingzhuangia',
                                                                            'Synechococcus_CC9902',
                                                                            'Propionigenium',
                                                                            'Phycisphaera',
                                                                            'Blastopirellula',
                                                                            'Rubripirellula',
                                                                            'Hellea',
                                                                            'Litorimos',
                                                                            'Amylibacter',
                                                                            'Afipia',
                                                                            'Pseudoalteromos',
                                                                            'Psychromos',
                                                                            'OM60(NOR5)_clade',
                                                                            'Granulosicoccus',
                                                                            'Halomos',
                                                                            'Psychrobacter',
                                                                            'Leucothrix',
                                                                            'Lentimos',
                                                                            'Persicirhabdus',
                                                                            'Roseibacillus',
                                                                            'Rubritalea',
                                                                            'OtherGenera'))
ggplot(Genus.SBP.surface.juv, aes(fill=genus, y=RelativeAbundance, x=FrondandBlade)) +
  geom_bar(position="stack", stat="identity", color="black") +
  theme_bw() + xlab("Juvenile Surface Samples") + ylab("% Relative Abundance") +
  scale_fill_manual(values = c("cornflowerblue", "cyan", "darkcyan",
                               "blue", "blue4", "darkslateblue",
                               "gray75", "floralwhite", "indianred",
                               "lightcoral", "indianred4", "bisque4",
                               "bisque1", "burlywood", "burlywood1",
                               "palegreen", "palegreen4", "olivedrab",
                               "lawngreen", "green3", "yellowgreen",
                               "darkolivegreen1", "chocolate", "chocolate4",
                               "yellow1", "yellow3", "black")) 

#genus surface mature blades
Genus.SBP.surface.mat <- subset(Genus.SBP, Depth=="surface" & Age=="Mature")
Genus.SBP.surface.mat$FrondandBlade<- as.character(Genus.SBP.surface.mat$FrondandBlade)
Genus.SBP.surface.mat$genus <- factor(Genus.SBP.surface.mat$genus, levels=c('Lewinella',
                                                                            'Rubidimos',
                                                                            'Croceitalea',
                                                                            'Dokdonia',
                                                                            'Flavicella',
                                                                            'Wenyingzhuangia',
                                                                            'Synechococcus_CC9902',
                                                                            'Propionigenium',
                                                                            'Phycisphaera',
                                                                            'Blastopirellula',
                                                                            'Rubripirellula',
                                                                            'Hellea',
                                                                            'Litorimos',
                                                                            'Amylibacter',
                                                                            'Afipia',
                                                                            'Pseudoalteromos',
                                                                            'Psychromos',
                                                                            'OM60(NOR5)_clade',
                                                                            'Granulosicoccus',
                                                                            'Halomos',
                                                                            'Psychrobacter',
                                                                            'Leucothrix',
                                                                            'Lentimos',
                                                                            'Persicirhabdus',
                                                                            'Roseibacillus',
                                                                            'Rubritalea',
                                                                            'OtherGenera'))
ggplot(Genus.SBP.surface.mat, aes(fill=genus, y=RelativeAbundance, x=FrondandBlade)) +
  geom_bar(position="stack", stat="identity", color="black") +
  theme_bw() + xlab("Mature Surface Samples") + ylab("% Relative Abundance") +
  scale_fill_manual(values = c("cornflowerblue", "cyan", "darkcyan",
                               "blue", "blue4", "darkslateblue",
                               "gray75", "floralwhite", "indianred",
                               "lightcoral", "indianred4", "bisque4",
                               "bisque1", "burlywood", "burlywood1",
                               "palegreen", "palegreen4", "olivedrab",
                               "lawngreen", "green3", "yellowgreen",
                               "darkolivegreen1", "chocolate", "chocolate4",
                               "yellow1", "yellow3", "black")) 

#genus middle juvenile blades
Genus.SBP.middle.juv <- subset(Genus.SBP, Depth=="middle" & Age=="Juvenile")
Genus.SBP.middle.juv$FrondandBlade<- as.character(Genus.SBP.middle.juv$FrondandBlade)
Genus.SBP.middle.juv$genus <- factor(Genus.SBP.middle.juv$genus, levels=c('Lewinella',
                                                                            'Rubidimos',
                                                                            'Croceitalea',
                                                                            'Dokdonia',
                                                                            'Flavicella',
                                                                            'Wenyingzhuangia',
                                                                            'Synechococcus_CC9902',
                                                                            'Propionigenium',
                                                                            'Phycisphaera',
                                                                            'Blastopirellula',
                                                                            'Rubripirellula',
                                                                            'Hellea',
                                                                            'Litorimos',
                                                                            'Amylibacter',
                                                                            'Afipia',
                                                                            'Pseudoalteromos',
                                                                            'Psychromos',
                                                                            'OM60(NOR5)_clade',
                                                                            'Granulosicoccus',
                                                                            'Halomos',
                                                                            'Psychrobacter',
                                                                            'Leucothrix',
                                                                            'Lentimos',
                                                                            'Persicirhabdus',
                                                                            'Roseibacillus',
                                                                            'Rubritalea',
                                                                            'OtherGenera'))
ggplot(Genus.SBP.middle., aes(fill=genus, y=RelativeAbundance, x=FrondandBlade)) +
  geom_bar(position="stack", stat="identity", color="black") +
  theme_bw() + xlab("Juvenile Middle Samples") + ylab("% Relative Abundance") +
  scale_fill_manual(values = c("cornflowerblue", "cyan", "darkcyan",
                               "blue", "blue4", "darkslateblue",
                               "gray75", "floralwhite", "indianred",
                               "lightcoral", "indianred4", "bisque4",
                               "bisque1", "burlywood", "burlywood1",
                               "palegreen", "palegreen4", "olivedrab",
                               "lawngreen", "green3", "yellowgreen",
                               "darkolivegreen1", "chocolate", "chocolate4",
                               "yellow1", "yellow3", "black")) 

#genus bottom juvenile
Genus.SBP.bottom.juv <- subset(Genus.SBP, Depth=="bottom" & Age=="Mature")
Genus.SBP.middle.mat$FrondandBlade<- as.character(Genus.SBP.middle.mat$FrondandBlade)
Genus.SBP.middle.mat$genus <- factor(Genus.SBP.middle.mat$genus, levels=c('Lewinella',
                                                                          'Rubidimos',
                                                                          'Croceitalea',
                                                                          'Dokdonia',
                                                                          'Flavicella',
                                                                          'Wenyingzhuangia',
                                                                          'Synechococcus_CC9902',
                                                                          'Propionigenium',
                                                                          'Phycisphaera',
                                                                          'Blastopirellula',
                                                                          'Rubripirellula',
                                                                          'Hellea',
                                                                          'Litorimos',
                                                                          'Amylibacter',
                                                                          'Afipia',
                                                                          'Pseudoalteromos',
                                                                          'Psychromos',
                                                                          'OM60(NOR5)_clade',
                                                                          'Granulosicoccus',
                                                                          'Halomos',
                                                                          'Psychrobacter',
                                                                          'Leucothrix',
                                                                          'Lentimos',
                                                                          'Persicirhabdus',
                                                                          'Roseibacillus',
                                                                          'Rubritalea',
                                                                          'OtherGenera'))
ggplot(Genus.SBP.middle.mat, aes(fill=genus, y=RelativeAbundance, x=FrondandBlade)) +
  geom_bar(position="stack", stat="identity", color="black") +
  theme_bw() + xlab("Mature Middle Samples") + ylab("% Relative Abundance") +
  scale_fill_manual(values = c("cornflowerblue", "cyan", "darkcyan",
                               "blue", "blue4", "darkslateblue",
                               "gray75", "floralwhite", "indianred",
                               "lightcoral", "indianred4", "bisque4",
                               "bisque1", "burlywood", "burlywood1",
                               "palegreen", "palegreen4", "olivedrab",
                               "lawngreen", "green3", "yellowgreen",
                               "darkolivegreen1", "chocolate", "chocolate4",
                               "yellow1", "yellow3", "black")) 

#genus middle mature
Genus.SBP.middle.mat <- subset(Genus.SBP, Depth=="middle" & Age=="Mature")
Genus.SBP.middle.mat$FrondandBlade<- as.character(Genus.SBP.middle.mat$FrondandBlade)
Genus.SBP.middle.mat$genus <- factor(Genus.SBP.middle.mat$genus, levels=c('Lewinella',
                                                                          'Rubidimos',
                                                                          'Croceitalea',
                                                                          'Dokdonia',
                                                                          'Flavicella',
                                                                          'Wenyingzhuangia',
                                                                          'Synechococcus_CC9902',
                                                                          'Propionigenium',
                                                                          'Phycisphaera',
                                                                          'Blastopirellula',
                                                                          'Rubripirellula',
                                                                          'Hellea',
                                                                          'Litorimos',
                                                                          'Amylibacter',
                                                                          'Afipia',
                                                                          'Pseudoalteromos',
                                                                          'Psychromos',
                                                                          'OM60(NOR5)_clade',
                                                                          'Granulosicoccus',
                                                                          'Halomos',
                                                                          'Psychrobacter',
                                                                          'Leucothrix',
                                                                          'Lentimos',
                                                                          'Persicirhabdus',
                                                                          'Roseibacillus',
                                                                          'Rubritalea',
                                                                          'OtherGenera'))
ggplot(Genus.SBP.middle.mat, aes(fill=genus, y=RelativeAbundance, x=FrondandBlade)) +
  geom_bar(position="stack", stat="identity", color="black") +
  theme_bw() + xlab("Mature Middle Samples") + ylab("% Relative Abundance") +
  scale_fill_manual(values = c("cornflowerblue", "cyan", "darkcyan",
                               "blue", "blue4", "darkslateblue",
                               "gray75", "floralwhite", "indianred",
                               "lightcoral", "indianred4", "bisque4",
                               "bisque1", "burlywood", "burlywood1",
                               "palegreen", "palegreen4", "olivedrab",
                               "lawngreen", "green3", "yellowgreen",
                               "darkolivegreen1", "chocolate", "chocolate4",
                               "yellow1", "yellow3", "black")) 

#genus bottom mature
Genus.SBP.bottom.mat <- subset(Genus.SBP, Depth=="bottom" & Age=="Mature")
Genus.SBP.bottom.mat$FrondandBlade<- as.character(Genus.SBP.bottom.mat$FrondandBlade)
Genus.SBP.bottom.mat$genus <- factor(Genus.SBP.bottom.mat$genus, levels=c('Lewinella',
                                                                          'Rubidimos',
                                                                          'Croceitalea',
                                                                          'Dokdonia',
                                                                          'Flavicella',
                                                                          'Wenyingzhuangia',
                                                                          'Synechococcus_CC9902',
                                                                          'Propionigenium',
                                                                          'Phycisphaera',
                                                                          'Blastopirellula',
                                                                          'Rubripirellula',
                                                                          'Hellea',
                                                                          'Litorimos',
                                                                          'Amylibacter',
                                                                          'Afipia',
                                                                          'Pseudoalteromos',
                                                                          'Psychromos',
                                                                          'OM60(NOR5)_clade',
                                                                          'Granulosicoccus',
                                                                          'Halomos',
                                                                          'Psychrobacter',
                                                                          'Leucothrix',
                                                                          'Lentimos',
                                                                          'Persicirhabdus',
                                                                          'Roseibacillus',
                                                                          'Rubritalea',
                                                                          'OtherGenera'))
ggplot(Genus.SBP.bottom.mat, aes(fill=genus, y=RelativeAbundance, x=FrondandBlade)) +
  geom_bar(position="stack", stat="identity", color="black") +
  theme_bw() + xlab("Bottom Mature Samples") + ylab("% Relative Abundance") +
  scale_fill_manual(values = c("cornflowerblue", "cyan", "darkcyan",
                               "blue", "blue4", "darkslateblue",
                               "gray75", "floralwhite", "indianred",
                               "lightcoral", "indianred4", "bisque4",
                               "bisque1", "burlywood", "burlywood1",
                               "palegreen", "palegreen4", "olivedrab",
                               "lawngreen", "green3", "yellowgreen",
                               "darkolivegreen1", "chocolate", "chocolate4",
                               "yellow1", "yellow3", "black")) 

#genus bottom juvenile
Genus.SBP.bottom.juv <- subset(Genus.SBP, Depth=="bottom" & Age=="Juvenile")
Genus.SBP.bottom.juv$FrondandBlade<- as.character(Genus.SBP.bottom.juv$FrondandBlade)
Genus.SBP.bottom.juv$genus <- factor(Genus.SBP.bottom.juv$genus, levels=c('Lewinella',
                                                                          'Rubidimos',
                                                                          'Croceitalea',
                                                                          'Dokdonia',
                                                                          'Flavicella',
                                                                          'Wenyingzhuangia',
                                                                          'Synechococcus_CC9902',
                                                                          'Propionigenium',
                                                                          'Phycisphaera',
                                                                          'Blastopirellula',
                                                                          'Rubripirellula',
                                                                          'Hellea',
                                                                          'Litorimos',
                                                                          'Amylibacter',
                                                                          'Afipia',
                                                                          'Pseudoalteromos',
                                                                          'Psychromos',
                                                                          'OM60(NOR5)_clade',
                                                                          'Granulosicoccus',
                                                                          'Halomos',
                                                                          'Psychrobacter',
                                                                          'Leucothrix',
                                                                          'Lentimos',
                                                                          'Persicirhabdus',
                                                                          'Roseibacillus',
                                                                          'Rubritalea',
                                                                          'OtherGenera'))
ggplot(Genus.SBP.bottom.juv, aes(fill=genus, y=RelativeAbundance, x=FrondandBlade)) +
  geom_bar(position="stack", stat="identity", color="black") +
  theme_bw() + xlab("Bottom Juvenile Samples") + ylab("% Relative Abundance") +
  scale_fill_manual(values = c("cornflowerblue", "cyan", "darkcyan",
                               "blue", "blue4", "darkslateblue",
                               "gray75", "floralwhite", "indianred",
                               "lightcoral", "indianred4", "bisque4",
                               "bisque1", "burlywood", "burlywood1",
                               "palegreen", "palegreen4", "olivedrab",
                               "lawngreen", "green3", "yellowgreen",
                               "darkolivegreen1", "chocolate", "chocolate4",
                               "yellow1", "yellow3", "black")) 

##
#Second - Ordination Plots -------
Ordi <- read.csv("Relative Abundance - For Ordinations.csv")
View(Ordi)

library(ggplot2)
library(ggcorrplot)
library(corrr)
library(FactoMineR)
library(usethis)
library(devtools)
library("factoextra")
library(psych)

#PCA first
data_normalized <- scale(Ordi)
head(data_normalized)
View(data_normalized)
corr_matrix <- as.data.frame(cor(data_normalized))
View(corr_matrix)

data.pca <- princomp(~ ., data = corr_matrix, cor = TRUE)
summary(data.pca)
Ordi.coords <- data.pca$loadings[, 1:2]
View(Ordi.coords)
fviz_eig(data.pca, addlabels = TRUE) ##x axis explains (37.2%) & yaxis explains (24.1%)
#xlim is -1,1 & ylim is -1,1
# Graph of the variables
fviz_pca_var(data.pca, col.var = "black")


screeplot(data.pca, type = "l", npcs = 15)
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)

fviz_pca_var(data.pca, col.var = "black")ox.pri
View(data.pca$scores)
write.csv(data.pca$scores, "Principal Component - Coordinates.csv")

#ggplot PCA coords
PC.coords <- read.csv("Principal Component - Coordinates.csv")
View(PC.coords)
PC.coords.mature <- subset(PC.coords, Age=="Mature")
PC.coords.juvenile <- subset(PC.coords, Age=="Juvenile")
View(PC.coords.mature)
ggplot(PC.coords.juvenile, aes(x=Comp.1, y=Comp.2, color=Age, shape=Depth)) +
  geom_point(alpha=0.7, size=5, color='yellow3') +
  theme_bw() +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  xlim(-12.8, 7.6) + ylim(-7.7,11) +
  xlab("PC1 (37.2%)") + ylab("PC2 (24.1%)")

#Now I want to retry the principal component without scaling
corr_matrix_noscale <- as.data.frame(cor(Ordi))
data.pca_noscale <- princomp(~ ., data = corr_matrix_noscale, cor = TRUE)
summary(data.pca_noscale)
Ordi.coords_noscale <- data.pca_noscale$loadings[, 1:2]
View(Ordi.coords_noscale)
fviz_eig(data.pca_noscale, addlabels = TRUE) ##x axis explains (37.2%) & yaxis explains (24.1%)
#xlim is -1,1 & ylim is -1,1
# Graph of the variables
fviz_pca_var(data.pca_noscale, col.var = "black")
data.pca_noscale

screeplot(data.pca_noscale, type = "l", npcs = 15)
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)

View(data.pca_noscale$scores)
write.csv(data.pca_noscale$scores, "Principal Component - Not Scaled - Coordinates.csv")

PC.coords_notscaled <- read.csv("Principal Component - Not Scaled - Coordinates.csv")
View(PC.coords_notscaled)
PC.coords.mature_notscaled <- subset(PC.coords_notscaled, Age=="Mature")
PC.coords.juvenile_notscaled <- subset(PC.coords_notscaled, Age=="Juvenile")
View(PC.coords.mature)
ggplot(PC.coords.juvenile_notscaled, aes(x=Comp.1, y=Comp.2, color=Age, shape=Depth)) +
  geom_point(alpha=0.7, size=5, color='yellow3') +
  theme_bw() +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  xlim(-12.8, 7.6) + ylim(-7.7,11) +
  xlab("PC1 (37.2%)") + ylab("PC2 (24.1%)") ##same results

#PCoA --------
library(ecodist)
library(vegan)
library(ggplot2)
data("varespec")
View(varespec)
Ordi.PCoA <- read.csv("Percent Relative Abundance - For PCoA - Ordination.csv")
View(Ordi.PCoA)
?vegdist
Ordi.bray <- vegdist(Ordi.PCoA, method = "bray") # dissimilarity matrix using bray-curtis distance indices on the varespec dataset native to vegan
pcoaVS <- pco(Ordi.bray, negvals = "zero", dround = 0) # if negvals = 0 sets all negative eigenvalues to zero; if = "rm" corrects for negative eigenvalues using method 1 of Legendre and Anderson 1999
pcoaVS$vector
plot(pcoaVS$vectors[,1], pcoaVS$vectors[,2], type = "n", xlab = "PCoA1", ylab = "PCoA2",
     axes = TRUE, main = "PCoA (ecodist) on varespec data")

text(pcoaVS$vectors[,1], pcoaVS$vectors[,2], labels(Ordi.bray), 
     cex = 0.9, xpd = TRUE)
pcoaVS$values # eigenvalue for each component. This is a measure of the variance explained by each dimension
pcoaVS$vectors # eigenvectors. Each column contains the scores for that dimension.
View(pcoaVS$values)
write.csv(pcoaVS$vectors, "PCoA - Bray Curtis - Coordinates.csv")

PCoA.coords <- read.csv("PCoA - Bray Curtis - Coordinates.csv")
View(PCoA.coords)
PCoA.coords.mature <- subset(PCoA.coords, Age=="Mature")
PCoA.coords.juvenile <- subset(PCoA.coords, Age=="Juvenile")

ggplot(PCoA.coords.juvenile, aes(x=PCoA1, y=PCoA2, color=Age, shape=Depth)) +
  geom_point(alpha=0.7, size=5, color='yellow3') +
  theme_bw() +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  xlim(-0.1, 0.09) + ylim(-0.09,0.12) 

ggplot(PCoA.coords.mature, aes(x=PCoA1, y=PCoA2, color=Age, shape=Depth)) +
  geom_point(alpha=0.7, size=5, color='darkblue') +
  theme_bw() +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  xlim(-0.1, 0.09) + ylim(-0.09,0.12) 

#PCoA with unifrac and phyloseq--------
library(phyloseq)
#First, we need to create a phyloseq object
?ordinate
ASVMatrix <- read.csv("ASVMatrix - For Phyloseq Unifrac PCoA.csv")
View(ASVMatrix)
TaxaMatrix <- read.csv("Species Coordinates.csv")
View(TaxaMatrix)
SampleMatrix <- read.csv("Sample Coordinates - Microbiome Only.csv")
View(SampleMatrix)

ASVMatrix <- ASVMatrix %>%
  tibble::column_to_rownames("AmpliconSequenceVariant")
View(ASVMatrix)
TaxaMatrix <- TaxaMatrix %>% 
  tibble::column_to_rownames("AmpliconSequenceVariant")
SampleMatrix <- SampleMatrix %>% 
  tibble::column_to_rownames("UniqueID") 
View(SampleMatrix)

asv_mat <- as.matrix(ASVMatrix)
tax_mat <- as.matrix(TaxaMatrix)

OTU = otu_table(asv_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(SampleMatrix)

phylo.object <- phyloseq(OTU, TAX, samples)
phylo.object


library("ape")
random_tree = rtree(ntaxa(phylo.object), rooted=TRUE, tip.label=taxa_names(phylo.object))
plot(random_tree)
physeq1 = merge_phyloseq(phylo.object, samples, random_tree)
physeq1

ord_PCoA_unifrac = ordinate(physeq1, "PCoA", "unifrac", weighted=TRUE)
ord_PCoA_jaccard = ordinate(physeq1, "PCoA", "jaccard", weighted=TRUE)
ord_PCoA_bray = ordinate(physeq1, "PCoA", "bray", weighted=TRUE)

ord_PCoA_unifrac$values
ord_PCoA_jaccard
View(ord_PCoA_jaccard$values)
View(ord_PCoA_jaccard$vectors)
ord_PCoA_bray$values
###write.csv(ord_PCoA_unifrac$vectors, "Unifrac PCoA Weighted.csv")
##write.csv(ord_PCoA_jaccard$vectors, "Jaccard PCoA.csv")
###write.csv(ord_PCoA_bray$vectors, "Bray PCoA.csv")
Unifrac.weighted <- read.csv("Unifrac PCoA Weighted.csv")
View(Unifrac.weighted)

Unifrac.weighted.mature <- subset(Unifrac.weighted, Age=="Mature")

ggplot(Unifrac.weighted.mature, aes(x=Axis.1, y=Axis.2, color=Age, shape=Depth)) +
  geom_point(alpha=0.7, size=5, color='darkblue') +
  theme_bw() +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  xlim(-0.34, 0.2) + ylim(-0.21,0.22) +
  xlab("PCoA1 (24.39% of Variance)") +
  ylab("PCoA2 (15.14% of Variance)")

Unifrac.weighted.juvenile <- subset(Unifrac.weighted, Age=="Juvenile")

ggplot(Unifrac.weighted.juvenile, aes(x=Axis.1, y=Axis.2, color=Age, shape=Depth)) +
  geom_point(alpha=0.7, size=5, color='yellow3') +
  theme_bw() +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  xlim(-0.34, 0.2) + ylim(-0.21,0.22) +
  xlab("PCoA1 (24.39% of Variance)") +
  ylab("PCoA2 (15.14% of Variance)")

ggplot(Unifrac.weighted , aes(x=Axis.1, y=Axis.2, color=Depth, shape=Depth)) +
  geom_point(alpha=0.7, size=5) +
  theme_bw() +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  scale_color_manual(values=c('forestgreen', 'purple3')) +
  xlab("PC1 (24.39% of Variance)") +
  ylab("PC2 (15.14% of Variance)") +
  stat_ellipse(geom = "polygon", alpha=0.2,
               aes(color = Depth, fill=Depth)) +
  scale_fill_manual(values=c('forestgreen', 'purple3')) +
  theme(legend.text = element_text(size=14)) +
  theme(legend.title = element_text(size=14)) 

ggplot(Unifrac.weighted , aes(x=Axis.1, y=Axis.2, color=Age, shape=Age)) +
  geom_point(alpha=0.7, size=5) +
  theme_bw() +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  scale_color_manual(values=c('yellow3', 'darkblue')) +
  xlab("PC1 (24.39% of Variance)") +
  ylab("PC2 (15.14% of Variance)") +
  stat_ellipse(geom = "polygon", alpha=0.2,
               aes(color = Age, fill=Age)) +
  scale_fill_manual(values=c('yellow3', 'darkblue')) +
  theme(legend.text = element_text(size=14)) +
  theme(legend.title = element_text(size=14)) 



#PCoA with Jaccard and phyloseq --------
jacc.coords <- read.csv("Jaccard PCoA.csv")
jacc.coords.mat <- subset(jacc.coords, Age=="Mature")
jacc.coords.juv <- subset(jacc.coords, Agee="Juvenile")

View(jacc.coords)
ggplot(jacc.coords.mat, aes(x=Axis.1, y=Axis.2, color=Age, shape=Depth)) +
  geom_point(alpha=0.7, size=5, color='darkblue') +
  theme_bw() +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  xlim(-0.37, 0.38) + ylim(-0.44, 0.31) +
  xlab("PCoA1 (16.11% of Variation)") +
  ylab("PCoA2 (14.00 % of Variation")


ggplot(jacc.coords.juv, aes(x=Axis.1, y=Axis.2, color=Age, shape=Depth)) +
  geom_point(alpha=0.7, size=5, color='yellow3') +
  theme_bw() +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  xlim(-0.37, 0.38) + ylim(-0.44, 0.31) +
  xlab("PCoA1 (16.11% of Variation)") +
  ylab("PCoA2 (14.00 % of Variation")

bray.coords <- read.csv("Bray PCoA.csv")
View(bray.coords)
bray.coords.mat <- subset(bray.coords, Age=="Mature")
bray.coords.juv <- subset(bray.coords, Age=="Juvenile")
View(bray.coords.juv)

ggplot(bray.coords.juv, aes(x=Axis.1, y=Axis.2, color=Age, shape=Depth)) +
  geom_point(alpha=0.7, size=5, color='yellow3') +
  theme_bw() +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  xlim(-0.4, 0.36) + ylim(-0.32, 0.42) +
  xlab("PCoA1 (24.07% of Variation)") +
  ylab("PCoA2 (20.47 % of Variation)")

ggplot(bray.coords.mat, aes(x=Axis.1, y=Axis.2, color=Age, shape=Depth)) +
  geom_point(alpha=0.7, size=5, color='darkblue') +
  theme_bw() +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  xlim(-0.4, 0.36) + ylim(-0.32, 0.42) +
  xlab("PCoA1 (24.07% of Variation)") +
  ylab("PCoA2 (20.47 % of Variation)")

library(ggplot2)
