#Load Libraries --------------------------------
library(permute)
library(lattice)
library(tidyverse)
library(vegan)
library(plyr)
library(readr)
library(ggplot2)
library(gapminder)
library(lubridate)
library(gridExtra)
library(Rcpp)
library(dada2)
library(DESeq2)
library(phyloseq)
library(gdata)


#Histogram of ASV Reads-----
UniqueID_RawReadCount <- read_csv("UniqueID_RawReadCount.csv")
View(UniqueID_RawReadCount)
hist(UniqueID_RawReadCount$TotalRawReads, breaks = 20)

#Rarefaction --------------------------------
#Import the metadata of your Unique ID groupings 

sample_info_tab <- read.csv("UniqueID_grouping_NoTimeZero.csv")
View(sample_info_tab)

#Import your RAW ASV reads that have not yet been rarefied
#Read in already transposed asv raw reads table
asv_raw_good_reads_only <- read_csv("asv_raw_good_reads_only_NoTimeZero.csv") 
View(asv_raw_good_reads_only)

#Get the Unique ID's into the categorical far left column so that the dataframe is ONLY NUMBERS.
asv_raw_good_reads_only %>% remove_rownames %>% column_to_rownames(var="UniqueID")
View(asv_raw_good_reads_only)
asv2 <- column_to_rownames(asv_raw_good_reads_only, "UniqueID")
View(asv2)
#Rarefaction curve plotting
rarecurve(asv2, step=100, lwd=2 ,col=sample_info_tab$color2 ,ylab="ASVs", label=F) 
abline(v=(min(rowSums(asv2))))

#Actually rarefying the data using rrarefy
specnumber(asv2)
spAbund <- rowSums(asv2)
spAbund
raremin <- min(rowSums(asv2))
raremin
sRare <- rrarefy(asv2, raremin)
View(sRare) 
write.csv(sRare, "~/Kelp/Kelp_Data/16S_Library/16S Analysis No Time Zero\\rarefied_ASVs_per_sample_NoTimeZero.csv",row.names=TRUE)

#Gamma Diversity ----------------
Gamma_Diversity_NoTimeZero <- read_csv("Gamma_Diversity_NoTimeZero.csv")
View(Gamma_Diversity_NoTimeZero)

ggplot(Gamma_Diversity_NoTimeZero, aes(x=depth, y=gamma_diversity, fill=sample_type, color=sample_type)) +
  geom_boxplot(fatten=4) + theme_classic() + coord_flip()

ggplot(Gamma_Diversity_NoTimeZero, aes(x=depth, y=gamma_diversity, fill=sample_type, color=sample_type)) +
  geom_jitter(size=4) + theme_classic() + coord_flip()

ggplot(Gamma_Diversity_NoTimeZero, aes(x=depth, y=gamma_diversity, fill=sample_type, color=sample_type)) +
  geom_point(size=5, alpha=0.5) + theme_classic() + coord_flip()

# Alpha Diversity -------------
Alpha_Diversity_NoTimeZero <- read_csv("AlphaDiversity_NoTimeZero.csv")
View(Alpha_Diversity_NoTimeZero)
ggplot(Alpha_Diversity_NoTimeZero, aes(x=depth, y=Alpha_Diversity, fill=sampling_type)) +
  geom_boxplot() + theme_classic() + coord_flip()

#AOV and Tukey of Alpha Diversity Samples 
Alpha.lm <- lm(Alpha_Diversity_NoTimeZero$Alpha_Diversity~Alpha_Diversity_NoTimeZero$depth_and_sampling_type, data=Alpha_Diversity_NoTimeZero)
Alpha.av <- aov(Alpha.lm)
summary(Alpha.av)
Alpha.Tukey <- TukeyHSD(Alpha.av)
Alpha.Tukey

#BetaDisper with these rarefied data--------------
# Beta Disper
rar_counts<- read_csv("rarefied_ASVs_per_sample_NoTimeZero.csv")

rar_counts %>% remove_rownames %>% column_to_rownames(var="UniqueID")
rar_counts2 <- column_to_rownames(rar_counts, "UniqueID")
View(rar_counts2)
dis <- vegdist(rar_counts2)
dis
groups <- factor(c(rep(1,13),
                   rep(2,15),
                   rep(3,13),
                   rep(4,1),
                   rep(5,2),
                   rep(6,14),
                   rep(7,14),
                   rep(8,2),
                   rep(9,2),
                   rep(10,2),
                   rep(11,17),
                   rep(12,15),
                   rep(13,12),
                   rep(14,2),
                   rep(15,2)), labels = c("surface_CF 2 wks old",
                                          "surface_CF 4 wks old",
                                          "surface_SA 2 wks old",
                                          "surface_water 2 wk old",
                                          "surface_water 4 wk old",
                                          "middle_CF 2 wks old",
                                          "middle_CF 4 wks old",
                                          "middle_SA 2 wks old",
                                          "middle_water 2 wk old",
                                          "middle_water 4 wk old",
                                          "bottom_CF 2 wks old",
                                          "bottom_CF 4 wks old",
                                          "bottom_SA 2 wks old",
                                          "bottom_water 2 wk old",
                                          "bottom_water 4 wk old"))
mod <- betadisper(dis, groups)
mod

## Permutation test for F
permutest(mod, pairwise = TRUE, permutations = 99)

## Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)
summary(mod.HSD)
## Plot the groups and distances to centroids on the
## first two PCoA axes
plot(mod)
summary(mod)
boxplot(mod, lwd=1, rotate=45)
axis(side=2, labels=FALSE)

df <- data.frame(Distance_to_centroid=mod$distances, Group=mod$group)
df.PCoA <- data.frame(PCoA = mod$vectors, Group=mod$group)
View(df.PCoA)
summary(df.PCoA)
write.csv(df.PCoA, "~/Kelp/Kelp_Data/16S_Library/16S Analysis No Time Zero\\betadisper_PCoA.csv",row.names=TRUE)

groups <- mod$group
groups
View(df)
write.csv(df, "~/Kelp/Kelp_Data/16S_Library/16S Analysis No Time Zero\\betadisper.csv",row.names=TRUE)

betadisper <- read_csv("betadisper.csv")
View(betadisper)

ggplot(betadisper, aes(x=Distance_to_centroid, y=Group, fill=sampling_type)) +
  geom_boxplot()  + theme_classic()

betadisper_PCoA <- read_csv("betadisper_PCoA.csv")
View(betadisper_PCoA)


ggplot(betadisper_PCoA, aes(x=PCoA1, y=PCoA2)) +
  geom_point(aes(color=sampling_type, shape=depth, alpha=0.5, size=5)) + theme_classic() + theme(axis.text=element_text(size=18),
                                                                                                 axis.title=element_text(size=20))

#Load Alpha Diversity and Photophysiology data------
Alpha_Diversity_NoTimeZero <- read_csv("~/Kelp/Kelp_Data/16S_Library/16S Analysis No Time Zero/Alpha Diversity_NoTimeZero.csv")
View(Alpha_Diversity_NoTimeZero)

ggplot(Alpha_Diversity_NoTimeZero, aes(x=Alpha, y=CN)) + geom_point(aes(size=2, color=sampling_type, shape=depth)) + theme_classic() + geom_smooth(method=lm)
summary(lm(Alpha_Diversity_NoTimeZero$Alpha~Alpha_Diversity_NoTimeZero$CN))

ggplot(Alpha_Diversity_NoTimeZero, aes(x=Alpha, y=Surface_Area)) + geom_point(aes(size=2, color=sampling_type, shape=depth)) + theme_classic() + geom_smooth(method=lm)
summary(lm(Alpha_Diversity_NoTimeZero$Alpha~Alpha_Diversity_NoTimeZero$Surface_Area))

ggplot(Alpha_Diversity_NoTimeZero, aes(x=Alpha, y=ChlC)) + geom_point(aes(size=2, color=sampling_type, shape=depth)) + theme_classic() + geom_smooth(method=lm)
summary(lm(Alpha_Diversity_NoTimeZero$Alpha~Alpha_Diversity_NoTimeZero$ChlC))

ggplot(Alpha_Diversity_NoTimeZero, aes(x=Alpha, y=FvFmAvg)) + geom_point(aes(size=2, color=sampling_type, shape=depth)) + theme_classic() + geom_smooth(method=lm)
summary(lm(Alpha_Diversity_NoTimeZero$Alpha~Alpha_Diversity_NoTimeZero$FvFmAvg))

#NMDS with these Rarefied samples------------
BC.NMDS.shell<- function(dataframe) {
  metaMDS(dataframe, distance="bray", k=2, trymax=1000)
}
BC.NMDS <- BC.NMDS.shell(rar_counts2)
Blade_meta <- read_csv("~/Kelp/Kelp_Data/16S_Library/16S Analysis No Time Zero/Blade_meta_NoTimeZero.csv")
BC.NMDS_sample_coords <- scores(BC.NMDS, display = "sites") %>%
  as.data.frame() %>%
  #make my row names into a column
  rownames_to_column("UniqueID") %>%
  #join data frame with metadata
  full_join(., Blade_meta, by = "UniqueID")
View(BC.NMDS_sample_coords)

write.csv(BC.NMDS_sample_coords, "~/Kelp/Kelp_Data/16S_Library/16S Analysis No Time Zero\\NMDS_NoTimeZero.csv",row.names=TRUE)

NMDS_NoTimeZero <- read.csv("~/Kelp/Kelp_Data/16S_Library/16S Analysis No Time Zero/NMDS_NoTimeZero.csv")
View(NMDS_NoTimeZero)
ggplot(NMDS_NoTimeZero, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(size=4, alpha=0.5, shape=depth, color=sampling_type)) + theme_classic()

#NMDS with water Temperature
NMDS_NoTimeZero_WaterTemp <- read_csv("NMDS_NoTimeZero_WaterTemp.csv")
View(NMDS_NoTimeZero_WaterTemp)

ggplot(NMDS_NoTimeZero_WaterTemp, aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(shape=depth_and_sampling_type, color=WaterTemperature_ADCP), size=3, stroke=1.5) +
  scale_shape_manual(values=c(0, 7, 14, 1, 10, 13, 2, 6, 11)) +
  scale_color_gradient(low="red", high="blue") + theme_classic() +
  xlim(-1.8,0.75) + ylim(-1,1)

#Let's try visualizing this a different way

ggplot(NMDS_NoTimeZero_WaterTemp, aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(shape=sampling_type, size=depth, color=WaterTemperature_ADCP), alpha=0.65) +
  scale_color_gradient(low="blue1", high="red2") + theme_classic()+
  xlim(-1.8,0.75) + ylim(-1,1) +
  theme(text = element_text(size=20),
        axis.text.x = element_text(hjust=1))
  

#NMDS with CN
NMDS_NoTimeZero_CN <- read_csv("NMDS_NoTimeZero_CN.csv")
View(NMDS_NoTimeZero_CN)
ggplot(NMDS_NoTimeZero_CN, aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(shape=depth_and_sampling_type, color=CN), size=3, stroke=1.5) +
  scale_shape_manual(values=c(0, 7, 14, 1, 10, 13, 2, 6, 11)) +
  scale_color_gradient(low="green", high="orange") + theme_classic() +
  xlim(-1.8,0.75) + ylim(-1,1)

ggplot(NMDS_NoTimeZero_CN, aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(shape=sampling_type, size=depth, color=CN), alpha=0.65) +
  scale_color_gradient(low="green4", high="orange") + theme_classic()+
  xlim(-1.8,0.75) + ylim(-1,1) +
  theme(text = element_text(size=20),
        axis.text.x = element_text(hjust=1))

#NMDS with Surface Area
NMDS_NoTimeZero_SurfaceArea <- read_csv("NMDS_NoTimeZero_SurfaceArea.csv")
View(NMDS_NoTimeZero_SurfaceArea)
ggplot(NMDS_NoTimeZero_SurfaceArea, aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(shape=depth_and_sampling_type, color=Surface_Area), size=3, stroke=1.5) +
  scale_shape_manual(values=c(0, 7, 14, 1, 10, 13, 2, 6, 11)) +
  scale_color_gradient(low="yellow", high="purple") + theme_classic() +
  xlim(-1.8,0.75) + ylim(-1,1)

ggplot(NMDS_NoTimeZero_SurfaceArea, aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(shape=sampling_type, size=depth, color=Surface_Area), alpha=0.8) +
  scale_color_gradient(low="yellow", high="purple") + theme_classic()+
  xlim(-1.8,0.75) + ylim(-1,1) +
  theme(text = element_text(size=20),
        axis.text.x = element_text(hjust=1))

#NMDS with ChlC
NMDS_NoTimeZero_ChlC <- read_csv("NMDS_NoTimeZero_ChlC.csv")
View(NMDS_NoTimeZero_ChlC)
ggplot(NMDS_NoTimeZero_ChlC, aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(shape=depth_and_sampling_type, color=ChlC), size=3, stroke=1.5) +
  scale_shape_manual(values=c(0, 7, 14, 1, 10, 13, 2, 6, 11)) +
  scale_color_gradient(low="green3", high="pink") + theme_classic() +
  xlim(-1.8,0.75) + ylim(-1,1)

ggplot(NMDS_NoTimeZero_ChlC, aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(shape=sampling_type, size=depth, color=ChlC), alpha=0.8) +
  scale_color_gradient(low="pink", high="green4") + theme_classic()+
  xlim(-1.8,0.75) + ylim(-1,1) +
  theme(text = element_text(size=20),
        axis.text.x = element_text(hjust=1))

#NMDS with FvFm Avg
NMDS_NoTimeZero_FvFmAvg <- read_csv("NMDS_NoTimeZero_FvFmAvg.csv")
View(NMDS_NoTimeZero_FvFmAvg)
ggplot(NMDS_NoTimeZero_FvFmAvg, aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(shape=depth_and_sampling_type, color=FvFmAvg), size=4, stroke=1.5) +
  scale_shape_manual(values=c(0, 7, 14, 1, 10, 13, 2, 6, 11)) +
  scale_color_gradient(low="purple", high="orange") + theme_classic() +
  xlim(-1.8,0.75) + ylim(-1,1)

ggplot(NMDS_NoTimeZero_FvFmAvg, aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(shape=sampling_type, size=depth, color=FvFmAvg), alpha=0.65) +
  scale_color_gradient(low="purple", high="orange") + theme_classic()+
  xlim(-1.8,0.75) + ylim(-1,1) +
  theme(text = element_text(size=20),
        axis.text.x = element_text(hjust=1))

#Regression Analysis Re-Do Based On Holly's Comments------------
Alpha.Temp.Surf <- subset(Alpha_Diversity_NoTimeZero, depth=="surface")
View(Alpha.Temp.Surf)
Alpha.Temp.Mid <- subset(Alpha_Diversity_NoTimeZero, depth=="middle")
View(Alpha.Temp.Mid)
Alpha.Temp.Bott <- subset(Alpha_Diversity_NoTimeZero, depth=="bottom")
View(Alpha.Temp.Bott)

#LTER Temperature All Samples
ggplot(Alpha_Diversity_NoTimeZero, aes(x=WaterTemperature_ADCP, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type, shape=depth)) + theme_classic() + geom_smooth(method=lm)
summary(lm(Alpha_Diversity_NoTimeZero$Alpha~Alpha_Diversity_NoTimeZero$WaterTemperature_ADCP))
#LTER Temperature Surface Samples
ggplot(Alpha.Temp.Surf, aes(x=WaterTemperature_ADCP, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type), shape=0) + theme_classic() + geom_smooth(method=lm)
summary(lm(Alpha.Temp.Surf$Alpha~Alpha.Temp.Surf$WaterTemperature_ADCP))
#LTER Temperature Mid-Water Samples
ggplot(Alpha.Temp.Mid, aes(x=WaterTemperature_ADCP, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type), shape=24) + theme_classic() + geom_smooth(method=lm)
summary(lm(Alpha.Temp.Mid$Alpha~Alpha.Temp.Mid$WaterTemperature_ADCP))
#LTER Temperature Bottom Samples
ggplot(Alpha.Temp.Bott, aes(x=WaterTemperature_ADCP, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type), shape=1) + theme_classic() + geom_smooth(method=lm)
summary(lm(Alpha.Temp.Bott$Alpha~Alpha.Temp.Bott$WaterTemperature_ADCP))

#Further Subset to CF 2 wks old and SA 2 wks old, then re-run analysis
Alpha.Temp.Surf.2wkold <- subset(Alpha.Temp.Surf, sampling_type=="CF 2 wks old" | sampling_type== "SA 2 wks old")
View(Alpha.Temp.Surf.2wkold)
Alpha.Temp.Mid.2wkold <- subset(Alpha.Temp.Mid, sampling_type=="CF 2 wks old" | sampling_type== "SA 2 wks old")
View(Alpha.Temp.Mid.2wkold)
Alpha.Temp.Bott.2wkold <- subset(Alpha.Temp.Bott, sampling_type=="CF 2 wks old" | sampling_type== "SA 2 wks old")
View(Alpha.Temp.Bott.2wkold)
#All 2 wk old samples combined
Alpha.Temp.2wkold <- subset(Alpha_Diversity_NoTimeZero, sampling_type=="CF 2 wks old" | sampling_type== "SA 2 wks old")
View(Alpha.Temp.2wkold)

#LTER Temp 2 wk old Surface
ggplot(Alpha.Temp.Surf.2wkold, aes(x=WaterTemperature_ADCP, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type), shape=0) + theme_classic() + geom_smooth(method=lm)
summary(lm(Alpha.Temp.Surf.2wkold$Alpha~Alpha.Temp.Surf.2wkold$WaterTemperature_ADCP))
#LTER Temp 2 wk old Middle
ggplot(Alpha.Temp.Mid.2wkold, aes(x=WaterTemperature_ADCP, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type), shape=24) + theme_classic() + geom_smooth(method=lm)
summary(lm(Alpha.Temp.Mid.2wkold$Alpha~Alpha.Temp.Mid.2wkold$WaterTemperature_ADCP))
#LTER Temp 2 wk old Bottom
ggplot(Alpha.Temp.Bott.2wkold, aes(x=WaterTemperature_ADCP, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type), shape=1) + theme_classic() + geom_smooth(method=lm)
summary(lm(Alpha.Temp.Bott.2wkold$Alpha~Alpha.Temp.Bott.2wkold$WaterTemperature_ADCP))
#All 2 wk old samples combined
ggplot(Alpha.Temp.2wkold, aes(x=WaterTemperature_ADCP, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type, shape=depth)) + theme_classic() + geom_smooth(method=lm)
summary(lm(Alpha.Temp.2wkold$Alpha~Alpha.Temp.2wkold$WaterTemperature_ADCP))

#CN
ggplot(Alpha_Diversity_NoTimeZero, aes(x=CN, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type, shape=depth)) + theme_classic() + geom_smooth(method=lm)
summary(lm(Alpha_Diversity_NoTimeZero$Alpha~Alpha_Diversity_NoTimeZero$CN))
#CN Surface Samples
ggplot(Alpha.Temp.Surf, aes(x=CN, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type), shape=0) + theme_classic() + geom_smooth(method=lm)
summary(lm(Alpha.Temp.Surf$Alpha~Alpha.Temp.Surf$CN))
#CN Mid Samples
ggplot(Alpha.Temp.Mid, aes(x=CN, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type), shape=24) + theme_classic() + geom_smooth(method=lm)
summary(lm(Alpha.Temp.Mid$Alpha~Alpha.Temp.Mid$CN))
summary(lm(Alpha.Temp.Mid$CN~Alpha.Temp.Mid$Alpha))
#CN Bottom Samples
ggplot(Alpha.Temp.Bott, aes(x=CN, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type), shape=1) + theme_classic() + geom_smooth(method=lm)
summary(lm(Alpha.Temp.Bott$Alpha~Alpha.Temp.Bott$CN))

ggplot(Alpha_Diversity_NoTimeZero, aes(x=depth, y=CN, fill=sampling_type)) +
  geom_boxplot() + coord_flip() + theme_classic()

CN.lm <- lm(Alpha_Diversity_NoTimeZero$CN~Alpha_Diversity_NoTimeZero$category, data=Alpha_Diversity_NoTimeZero)
CN.av <- aov(CN.lm)
summary(CN.lm)
CN.Tukey <- TukeyHSD(CN.av)
CN.Tukey

#ChlC
#Surface ChlC
Surf.ChlC <- ggplot(Alpha.Temp.Surf, aes(x=ChlC, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type), shape=0) + theme_classic() + geom_smooth(method=lm) + theme(legend.position="none")
Surf.ChlC
summary(lm(Alpha.Temp.Surf$ChlC~Alpha.Temp.Surf$Alpha))
#Mid-Water ChlC
Mid.ChlC <- ggplot(Alpha.Temp.Mid, aes(x=ChlC, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type), shape=24) + theme_classic() + geom_smooth(method=lm) + theme(legend.position="none")
Mid.ChlC
summary(lm(Alpha.Temp.Mid$ChlC~Alpha.Temp.Mid$Alpha))
#Bottom-Water ChlC
Bottom.ChlC <- ggplot(Alpha.Temp.Bott, aes(x=ChlC, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type), shape=1) + theme_classic() + geom_smooth(method=lm) + theme(legend.position="none")
Bottom.ChlC
summary(lm(Alpha.Temp.Bott$ChlC~Alpha.Temp.Bott$Alpha))
#All Combined ChlC
All.ChlC <- ggplot(Alpha_Diversity_NoTimeZero, aes(x=ChlC, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type, shape=depth)) + theme_classic() + geom_smooth(method=lm)
summary(lm(Alpha_Diversity_NoTimeZero$ChlC~Alpha_Diversity_NoTimeZero$Alpha))

grid.arrange(Surf.ChlC, Mid.ChlC, Bottom.ChlC, All.ChlC, ncol=2)

#Boxplot of ChlC
ggplot(Alpha_Diversity_NoTimeZero, aes(x=depth, y=ChlC, fill=sampling_type)) +
  geom_boxplot() + coord_flip() + theme_classic() + theme(axis.text=element_text(size=18),
                                                          axis.title=element_text(size=20))

ChlC.lm <- lm(Alpha_Diversity_NoTimeZero$ChlC~Alpha_Diversity_NoTimeZero$category, data=Alpha_Diversity_NoTimeZero)
ChlC.av <- aov(ChlC.lm)
summary(ChlC.lm)
ChlC.Tukey <- TukeyHSD(ChlC.av)
ChlC.Tukey

#Surface Area
Surf.SA <- ggplot(Alpha.Temp.Surf, aes(x=Surface_Area, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type), shape=0) + theme_classic() + geom_smooth(method=lm) + theme(legend.position="none")
Surf.SA
summary(lm(Alpha.Temp.Surf$Surface_Area~Alpha.Temp.Surf$Alpha))
#Mid-Water Surface Area
Mid.SA <- ggplot(Alpha.Temp.Mid, aes(x=Surface_Area, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type), shape=24) + theme_classic() + geom_smooth(method=lm) + theme(legend.position="none")
Mid.SA
summary(lm(Alpha.Temp.Mid$Surface_Area~Alpha.Temp.Mid$Alpha))
#Bottom-Water Surface Area
Bottom.SA <- ggplot(Alpha.Temp.Bott, aes(x=Surface_Area, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type), shape=1) + theme_classic() + geom_smooth(method=lm) + theme(legend.position="none")
Bottom.SA
summary(lm(Alpha.Temp.Bott$Surface_Area~Alpha.Temp.Bott$Alpha))
#All Combined Surface Area
All.SA <- ggplot(Alpha_Diversity_NoTimeZero, aes(x=Surface_Area, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type, shape=depth)) + theme_classic() + geom_smooth(method=lm)
All.SA
summary(lm(Alpha_Diversity_NoTimeZero$Surface_Area~Alpha_Diversity_NoTimeZero$Alpha))

grid.arrange(Surf.SA, Mid.SA, Bottom.SA, All.SA, ncol=2)

#Boxplot of Surface Area
ggplot(Alpha_Diversity_NoTimeZero, aes(x=depth, y=Surface_Area, fill=sampling_type)) +
  geom_boxplot() + coord_flip() + theme_classic() + theme(axis.text=element_text(size=18),
                                                          axis.title=element_text(size=20))

SA.lm <- lm(Alpha_Diversity_NoTimeZero$Surface_Area~Alpha_Diversity_NoTimeZero$category, data=Alpha_Diversity_NoTimeZero)
SA.av <- aov(SA.lm)
summary(SA.lm)
SA.Tukey <- TukeyHSD(SA.av)
SA.Tukey

#Surface FvFm Averages
SU.FvFmAvg <- ggplot(Alpha.Temp.Surf, aes(x=FvFmAvg, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type), shape=0) + theme_classic() + geom_smooth(method=lm) + theme(legend.position="none")
SU.FvFmAvg
summary(lm(Alpha.Temp.Surf$FvFmAvg~Alpha.Temp.Surf$Alpha))
summary(lm(Alpha.Temp.Surf$Alpha~Alpha.Temp.Surf$FvFmAvg))
#Mid-Water FvFm Avg
Mid.FvFmAvg <- ggplot(Alpha.Temp.Mid, aes(x=FvFmAvg, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type), shape=24) + theme_classic() + geom_smooth(method=lm) + theme(legend.position="none")
Mid.FvFmAvg
summary(lm(Alpha.Temp.Mid$FvFmAvg~Alpha.Temp.Mid$Alpha))
#Bottom-Water Surface Area
Bottom.FvFmAvg <- ggplot(Alpha.Temp.Bott, aes(x=FvFmAvg, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type), shape=1) + theme_classic() + geom_smooth(method=lm) + theme(legend.position="none")
Bottom.FvFmAvg
summary(lm(Alpha.Temp.Bott$FvFmAvg~Alpha.Temp.Bott$Alpha))
#All Combined FvFm
All.FvFm <- ggplot(Alpha_Diversity_NoTimeZero, aes(x=FvFmAvg, y=Alpha)) +
  geom_point(aes(size=2, color=sampling_type, shape=depth)) + theme_classic() + geom_smooth(method=lm)
All.FvFm
summary(lm(Alpha_Diversity_NoTimeZero$FvFmAvg~Alpha_Diversity_NoTimeZero$Alpha))

#Boxplot of all FvFm
 ggplot(Alpha_Diversity_NoTimeZero, aes(x=depth, y=FvFmAvg, fill=sampling_type)) +
  geom_boxplot() + coord_flip() + theme_classic() + theme(axis.text=element_text(size=18),
                                                          axis.title=element_text(size=20))

grid.arrange(SU.FvFmAvg, Mid.FvFmAvg, Bottom.FvFmAvg, All.FvFm, ncol=2)

FvFmAvg.lm <- lm(Alpha_Diversity_NoTimeZero$FvFmAvg~Alpha_Diversity_NoTimeZero$category, data=Alpha_Diversity_NoTimeZero)
FvFmAvg.av <- aov(FvFmAvg.lm)
summary(FvFmAvg.lm)
FvFmAvg.Tukey <- TukeyHSD(FvFmAvg.av)
FvFmAvg.Tukey

#GLMs------------------------------
#First, let's check the histogram of each of these variables
AlphaDiversity_NoTimeZero <- read_csv("AlphaDiversity_NoTimeZero.csv")
View(AlphaDiversity_NoTimeZero)
hist(AlphaDiversity_NoTimeZero$Alpha) #Alpha Diversity is normally distributed
hist(AlphaDiversity_NoTimeZero$WaterTemperature_ADCP) #not normally distributed
hist(AlphaDiversity_NoTimeZero$CN) #not normally distributed
hist(AlphaDiversity_NoTimeZero$Surface_Area) #pretty normally distributed
hist(AlphaDiversity_NoTimeZero$ChlC) #not normally distributed although I'm supposed to try log transforming
hist(AlphaDiversity_NoTimeZero$FvFmAvg) #not normally distributed

fit.gaussian <- glm(data=AlphaDiversity_NoTimeZero, Alpha~WaterTemperature_ADCP+CN+Surface_Area+ChlC+FvFmAvg, family=gaussian())
summary(fit.gaussian)

fit.poisson <- glm(data=AlphaDiversity_NoTimeZero, Alpha~WaterTemperature_ADCP+CN+Surface_Area+ChlC+FvFmAvg, family=poisson())
summary(fit.poisson)
confint(fit.poisson)
exp(coef(fit.poisson))
exp(confint(fit.poisson))
predict(fit.poisson, type = "response")
residuals(fit.poisson, type="deviance")
plot(fit.poisson)



pairs(AlphaDiversity_NoTimeZero[,6:10], lower.panel=NULL, pch=19)

fit.poisson.All <- glm(data=AlphaDiversity_NoTimeZero, 
                        Alpha~FvFmAvg+
                          CN+
                          ChlC+
                          Surface_Area+
                          WaterTemperature_ADCP+
                          FvFmAvg*CN+ 
                          FvFmAvg*ChlC+
                          FvFmAvg*Surface_Area+
                          FvFmAvg*WaterTemperature_ADCP+
                          CN*ChlC+
                          CN*Surface_Area+
                          CN*WaterTemperature_ADCP+
                          CN*Surface_Area +
                          CN*WaterTemperature_ADCP+
                          ChlC*Surface_Area+
                          ChlC*WaterTemperature_ADCP+
                          Surface_Area*WaterTemperature_ADCP+
                          FvFmAvg*CN*ChlC+
                          FvFmAvg*CN*Surface_Area+
                          FvFmAvg*CN*WaterTemperature_ADCP+
                          FvFmAvg*ChlC*Surface_Area+
                          FvFmAvg*ChlC*WaterTemperature_ADCP+
                          FvFmAvg*Surface_Area*WaterTemperature_ADCP+
                          CN*ChlC*Surface_Area+
                          CN*ChlC*WaterTemperature_ADCP+
                          ChlC*Surface_Area*WaterTemperature_ADCP+
                          FvFmAvg*CN*ChlC*Surface_Area+
                          FvFmAvg*ChlC*Surface_Area*WaterTemperature_ADCP+
                          FvFmAvg*CN*Surface_Area*WaterTemperature_ADCP+
                          FvFmAvg*CN*ChlC*WaterTemperature_ADCP+
                          CN*ChlC*Surface_Area*WaterTemperature_ADCP+
                          FvFmAvg*CN*ChlC*Surface_Area*WaterTemperature_ADCP
                          
                          ,family=poisson())
summary(fit.poisson.All)

summary(lm(AlphaDiversity_NoTimeZero$CN~AlphaDiversity_NoTimeZero$WaterTemperature_ADCP))
summary(lm(AlphaDiversity_NoTimeZero$Surface_Area~AlphaDiversity_NoTimeZero$WaterTemperature_ADCP))
summary(lm(AlphaDiversity_NoTimeZero$ChlC~AlphaDiversity_NoTimeZero$WaterTemperature_ADCP))
summary(lm(AlphaDiversity_NoTimeZero$FvFmAvg~AlphaDiversity_NoTimeZero$WaterTemperature_ADCP))

summary(lm(AlphaDiversity_NoTimeZero$Surface_Area~AlphaDiversity_NoTimeZero$CN))
summary(lm(AlphaDiversity_NoTimeZero$ChlC~AlphaDiversity_NoTimeZero$CN))
summary(lm(AlphaDiversity_NoTimeZero$FvFmAvg~AlphaDiversity_NoTimeZero$CN))

summary(lm(AlphaDiversity_NoTimeZero$ChlC~AlphaDiversity_NoTimeZero$Surface_Area))
summary(lm(AlphaDiversity_NoTimeZero$FvFmAvg~AlphaDiversity_NoTimeZero$Surface_Area))


summary(lm(AlphaDiversity_NoTimeZero$FvFmAvg~AlphaDiversity_NoTimeZero$ChlC))

#GLM for betadisper -----------------
betadisper_withenvparameters <- read_csv("betadisper_withenvparameters.csv")
View(betadisper_withenvparameters)

fit.betadisper.poisson.All <- glm(data=betadisper_withenvparameters, 
                       Distance_to_centroid~FvFmAvg+
                         CN+
                         ChlC+
                         Surface_Area+
                         WaterTemperature_ADCP+
                         FvFmAvg*CN+ 
                         FvFmAvg*ChlC+
                         FvFmAvg*Surface_Area+
                         FvFmAvg*WaterTemperature_ADCP+
                         CN*ChlC+
                         CN*Surface_Area+
                         CN*WaterTemperature_ADCP+
                         CN*Surface_Area +
                         CN*WaterTemperature_ADCP+
                         ChlC*Surface_Area+
                         ChlC*WaterTemperature_ADCP+
                         Surface_Area*WaterTemperature_ADCP+
                         FvFmAvg*CN*ChlC+
                         FvFmAvg*CN*Surface_Area+
                         FvFmAvg*CN*WaterTemperature_ADCP+
                         FvFmAvg*ChlC*Surface_Area+
                         FvFmAvg*ChlC*WaterTemperature_ADCP+
                         FvFmAvg*Surface_Area*WaterTemperature_ADCP+
                         CN*ChlC*Surface_Area+
                         CN*ChlC*WaterTemperature_ADCP+
                         ChlC*Surface_Area*WaterTemperature_ADCP+
                         FvFmAvg*CN*ChlC*Surface_Area+
                         FvFmAvg*ChlC*Surface_Area*WaterTemperature_ADCP+
                         FvFmAvg*CN*Surface_Area*WaterTemperature_ADCP+
                         FvFmAvg*CN*ChlC*WaterTemperature_ADCP+
                         CN*ChlC*Surface_Area*WaterTemperature_ADCP+
                         FvFmAvg*CN*ChlC*Surface_Area*WaterTemperature_ADCP
                       
                       ,family=poisson())
summary(fit.betadisper.poisson.All)


#Indicator Species Analysis ---------------
help("indicspecies")
library(indicspecies)
pc <- read_csv("rarefied_ASVs_per_sample_NoTimeZero.csv")

abund <- pc[,3:ncol(pc)]
time <- pc$depth_and_sampling_type
memory.limit(size=1)

inv = multipatt(abund, time, func = "r.g", control = how(nperm=1))

summary(inv)