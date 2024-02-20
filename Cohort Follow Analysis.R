library(readr)
library(ggplot2)
library(ggpattern)
library(dplyr)
library(forcats)
#NMDS of kelp blade vs water column -----------------

NMDS <- read_csv("NMDS Water vs Kelp.csv")
View(NMDS)

ggplot(NMDS, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(shape=timepoint, color=samplecoords, shape=timepoint, size=5, alpha=0.75)) +
  theme_classic() +
  scale_color_manual(values=c("#0000FF","#0000FF",
                              "#339933","#339933",
                              "#FF6600","#FF6600")) +
  theme(axis.text = element_text(size=20)) +
  theme(axis.title = element_text(size=18))

ggplot(NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(shape=depth, colour=(factor(source))), size=7, alpha=0.7) +
  theme_classic() +
  theme(axis.text = element_text(size=20)) +
  theme(axis.title = element_text(size=20)) 
  theme(plot.title = element_text(size=20))
#Kelp has different community composition than water.
  
#NMDS of kelp blades only ------------------
NMDS.blades <- read_csv("NMDS_CohortFollow_BladesOnly.csv")
View(NMDS.blades)

NMDS.Juvenile.blades <- subset(NMDS.blades, sampletype=="T1 Juv")
View(NMDS.Juvenile.blades)
NMDS.Mature.blades <- subset(NMDS.blades, sampletype=="T2 Mat")
View(NMDS.Mature.blades)

ggplot(NMDS.Juvenile.blades, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(shape=depth, color=depth, size=5, alpha=0.75)) +
  theme_classic() +
  scale_color_manual(values=c("#3333FF","#66FF33",
                              "#FF6600")) +
  theme(axis.text = element_text(size=20)) +
  theme(axis.title = element_text(size=18)) +
  xlim(-1.2, 0.7) + ylim(-0.75, 0.62)
  
ggplot(NMDS.Mature.blades, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(shape=depth, color=depth, size=5, alpha=0.75)) +
  theme_classic() +
  scale_color_manual(values=c("#003366","#339900",
                              "#CC3300")) +
  theme(axis.text = element_text(size=20)) +
  theme(axis.title = element_text(size=18)) +
  xlim(-1.2, 0.7) + ylim(-0.75, 0.62)

ggplot(NMDS.blades, aes(NMDS1, NMDS2)) +
  geom_point(aes(shape=depth, colour=(factor(sampletype))), size=7, alpha=0.7) +
  theme_classic() +
  theme(axis.text = element_text(size=20)) +
  theme(axis.title = element_text(size=20)) +
  scale_color_manual(values=c("#CC6600","#3300CC")) 
#Average by Class----------------
ClassAgg <- read_csv("Class Level Aggregate Sums per sample for Averaging.csv")
View(ClassAgg)

library(dplyr)
library(tidyr)

ClassSetUp <- gather(ClassAgg, "Class", "Class_Relative_Abundance", 3:95)
View(ClassSetUp)
Class.Avg.RA <- aggregate(data=ClassSetUp, Class_Relative_Abundance ~ depth_and_sampling_type + Class, FUN=mean, na.rm=TRUE)
View(Class.Avg.RA)
write.csv(Class.Avg.RA,"~/Kelp/Kelp_Data/16S_Library/Cohort Follow Story Only/Cohort Follow\\Averaged Class Aggregate Data.csv", row.names = FALSE)


Class.Stdev.RA <- aggregate(data=ClassSetUp, Class_Relative_Abundance ~ depth_and_sampling_type + Class, FUN=sd, na.rm=TRUE)
View(Class.Stdev.RA)
write.csv(Class.Stdev.RA,"~/Kelp/Kelp_Data/16S_Library/Cohort Follow Story Only/Cohort Follow\\Standard Deviation Class Aggregate Data.csv", row.names = FALSE)

#Plotting Averages and Standard Deviation of Class > 5% RA
Kelp.vs.Water <- read_csv("Class Level Kelp vs Water Calculation SetUp.csv")
View(Kelp.vs.Water)

Kelp.vs.Water.surface <- subset(Kelp.vs.Water, depth=="surface")
View(Kelp.vs.Water.surface)

ggplot(Kelp.vs.Water.surface, aes(x=Class, y=Average, group=sample_type, color=source, shape=sample_type)) +
  geom_pointrange(aes(ymin=Average-Stdev, ymax=Average+Stdev), width=.01,
                               position=position_dodge(0.5), size=1.4) + 
  coord_flip() + theme_classic()

ggplot(Kelp.vs.Water.surface, aes(fill=sample_type, pattern=source, y=Average, x=Class)) +
  geom_bar_pattern(stat='identity', position=position_dodge(preserve = "single"),
                   color="black",
                   pattern_fill="black",
                   pattern_angle = 45,
                   pattern_density=0.1,
                   pattern_spacing=0.025,
                   pattern_key_scale_factor = 0.6) +
  geom_errorbar(aes(ymin=Average, ymax=Average+Stdev), width=.3,
                position=position_dodge(.9)) +
  scale_pattern_manual(values=c(blade = "none", water = "stripe")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) + theme_minimal() +
  ylim(0,40)

Kelp.vs.Water.middle <- subset(Kelp.vs.Water, depth=="middle")
View(Kelp.vs.Water.middle)

ggplot(Kelp.vs.Water.middle, aes(fill=sample_type, pattern=source, y=Average, x=Class)) +
  geom_bar_pattern(stat='identity', position=position_dodge(preserve = "single"),
                   color="black",
                   pattern_fill="black",
                   pattern_angle = 45,
                   pattern_density=0.1,
                   pattern_spacing=0.025,
                   pattern_key_scale_factor = 0.6) +
  geom_errorbar(aes(ymin=Average, ymax=Average+Stdev), width=.3,
                position=position_dodge(.9)) +
  scale_pattern_manual(values=c(blade = "none", water = "stripe")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) + theme_minimal() +
  ylim(0,40)

Kelp.vs.Water.bottom <- subset(Kelp.vs.Water, depth=="bottom")
View(Kelp.vs.Water.bottom)

ggplot(Kelp.vs.Water.bottom, aes(fill=sample_type, pattern=source, y=Average, x=Class)) +
  geom_bar_pattern(stat='identity', position=position_dodge(preserve = "single"),
                   color="black",
                   pattern_fill="black",
                   pattern_angle = 45,
                   pattern_density=0.1,
                   pattern_spacing=0.025,
                   pattern_key_scale_factor = 0.6) +
  geom_errorbar(aes(ymin=Average, ymax=Average+Stdev), width=.3,
                position=position_dodge(.9)) +
  scale_pattern_manual(values=c(blade = "none", water = "stripe")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) + theme_minimal() +
  ylim(0,40)

#Family Level NMDS and Indicator Arrows ------------------
Family.Samples <- read_csv("Family Sample Coords Cohort Follow Only.csv")
View(Family.Samples)
Family.Coords <- read_csv("Indicator Families Species Coords Cohort Follow Only.csv")
View(Family.Coords)
ggplot(Family.Samples, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(shape=depth, color=BladeAge, size=5, alpha=0.5))+
  theme_bw()+ 
  theme(axis.text = element_text(size=20)) +
  theme(axis.title = element_text(size=18)) +
  scale_color_manual(values=c("#CC6600", "#3300CC")) +
  geom_segment(data=Family.Coords, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), alpha=0.2, arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data=Family.Coords, aes(x=NMDS1, y=NMDS2, label=Family), size=3) 
#Family Level Relative Abundance BarPlots-----------

FamilyAgg <- read_csv("Family Level Aggregates for Averaging and Stdev.csv")
View(FamilyAgg)
FamilySetUp <- gather(FamilyAgg, "Family", "Family_Relative_Abundance", 3:317)
View(FamilySetUp)
Family.Avg.RA <- aggregate(data=FamilySetUp, Family_Relative_Abundance ~ depth_and_sampling_type + Family, FUN=mean, na.rm=TRUE)
View(Family.Avg.RA)
write.csv(Family.Avg.RA,"~/Kelp/Kelp_Data/16S_Library/Cohort Follow Story Only/Cohort Follow\\Averaged Family Aggregate Data.csv", row.names = FALSE)

Family.Stdev.RA <- aggregate(data=FamilySetUp, Family_Relative_Abundance ~ depth_and_sampling_type + Family, FUN=sd, na.rm=TRUE)
View(Family.Stdev.RA)
write.csv(Family.Stdev.RA,"~/Kelp/Kelp_Data/16S_Library/Cohort Follow Story Only/Cohort Follow\\Stdev Family Aggregate Data.csv", row.names = FALSE)

FamilySetUp <- read_csv("Family Level Blade Calculation SetUp For Barplots.csv")
View(FamilySetUp)

FamilySetUp.Test <- subset(FamilySetUp, Family=="Alteromodaceae" | Family == "Cellvibrioceae")
View(FamilySetUp.Test)
ggplot(FamilySetUp.Test, aes(fill=sample_type, pattern=depth, y=Average, x=Family)) +
  geom_bar_pattern(width=.1, stat='identity', position=position_dodge(preserve = "single"),
                   color="black",
                   pattern_fill="black",
                   pattern_angle = 45,
                   pattern_density=0.1,
                   pattern_spacing=0.025,
                   pattern_key_scale_factor = 0.6) +
  geom_errorbar(aes(ymin=Average, ymax=Average+Stdev), width=.1,
                position=position_dodge(.9)) +
  scale_pattern_manual(values=c(surface = "none", middle = "stripe", bottom = "weave")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) + theme_minimal()

ggplot(FamilySetUp, aes(x=Family, y=Average, fill=sample_type)) +
  geom_bar(width=0.85, color="black", stat="identity", position=position_dodge(0.85)) + 
  coord_flip() +
  scale_fill_manual(values = c("#FFCC00", "#CC9900", "#99FF66", "#339900", "#99CCFF", "#3399FF")) +
  theme_minimal() +
  geom_errorbar(aes(ymin=Average, ymax=Average+Stdev), width=.85,
                position=position_dodge(0.85))# +
  theme(legend.position = "none")

FamilySetUpV2 <- read_csv("Family Level Blade Calculation SetUp For Barplots version2.csv")

ggplot(FamilySetUpV2, aes(x=Family, y=Average, fill=sample_type)) +
  geom_bar(width=0.85, color="black", stat="identity", position=position_dodge(0.85)) + 
  coord_flip() +
  scale_fill_manual(values = c("#FFCC00", "#FF6600", "#cc3300", "#99CCFF", "#3399ff", "#003366")) +
  theme_minimal() +
  geom_errorbar(aes(ymin=Average, ymax=Average+Stdev), width=.85,
                position=position_dodge(0.85)) #+
  theme(legend.position = "none")

#Family Level Ternary Plots --------------
library(ggtern)
JuvenileTern <- read.csv("Juvenile Blade Ternary Plot Cohort Follow.csv")  
View(JuvenileTern)
head(JuvenileTern)

ggtern(JuvenileTern, aes(x=ï..surface, y=middle, z=bottom, color=family)) +
  geom_point(aes(size=average), alpha=0.75) +
  theme_bvbw() +
  scale_size(range=c(1,10)) +
  #scale_shape_manual(values=c(21, 22, 23, 24, 25, 3, 4 )) +
  scale_color_manual(values = c("#33CCFF",
                                "#0099FF",
                                "#CCFF00",
                                "#66FF00",
                                "#999999",
                                "#0066FF",
                                "#000000",
                                "#0033FF",
                                "#FF33FF",
                                "#CC33CC",
                                "#000099",
                                "#FF0000",
                                "#FFCC66",
                                "#FF9900",
                                "#FF99FF",
                                "#660000",
                                "#00CC00",
                                "#6633FF",
                                "#006666",
                                "#993366")) +
  theme(legend.position = c(-1,0))

MatureTern <- read.csv("Mature Blade Ternary Plot Cohort Follow.csv")  
View(MatureTern)

ggtern(MatureTern, aes(x=ï..surface, y=middle, z=bottom)) +
  geom_point(aes(color=family, size=average), alpha=0.75) +
  theme_bvbw() +
  scale_size(range=c(1,10)) +
  #scale_shape_manual(values=c(21, 22, 23, 24, 25, 3, 4 )) +
  scale_color_manual(values = c("#33CCFF",
                                "#0099FF",
                                "#CCFF00",
                                "#66FF00",
                                "#999999",
                                "#0066FF",
                                "#000000",
                                "#0033FF",
                                "#FF33FF",
                                "#CC33CC",
                                "#000099",
                                "#FF0000",
                                "#FFCC66",
                                "#FF9900",
                                "#FF99FF",
                                "#660000",
                                "#00CC00",
                                "#6633FF",
                                "#006666",
                                "#993366")) +
  theme(legend.position = c(0,0.8))

ImportantFamilies <- read_csv("Juvenile and Mature Ternary Important Families.csv")
View(ImportantFamilies)

ggtern(ImportantFamilies, aes(x=surface, y=middle, z=bottom)) +
  geom_point(aes(color=family, size=average, shape=age), alpha=0.75) +
  theme_bvbw() +
  scale_size(range=c(1,10)) +
  #scale_shape_manual(values=c(21, 22, 23, 24, 25, 3, 4 )) +
  scale_color_manual(values = c("#33CCFF",
                                "#0099FF",
                                "#CCFF00",
                                "#999999",
                                "#0066FF",
                                "#FF33FF",
                                "#FF0000",
                                "#FFCC66",
                                "#FF9900",
                                "#FF99FF",
                                "#6633FF",
                                "#006666",
                                "#993366"))+
  theme(legend.position = c(0,-1))

#HeatMap Family Level ------------------------
HeatMapSetUp <- read_csv("HeatMapSetUp.csv")
View(HeatMapSetUp)


ggplot(HeatMapSetUp, aes(x=Depth, y=Family, fill=Log2AbsoluteChange)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#FF0000",
                       mid = "#FFFFFF",
                       high = "#000099") +
  geom_text(aes(label=Mature), color="black", size=5) +
  coord_fixed() + theme_minimal() 

ggplot(HeatMapSetUp, aes(x=Depth, y=Family, fill=Juvenile)) +
  geom_tile(color="black") +
  scale_fill_continuous(limits=c(0,14)) +
  scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn")) +
  geom_text(aes(label = Juvenile), color = "black", size = 5) +
  coord_fixed() + theme_minimal()

ggplot(HeatMapSetUp, aes(x=Depth, y=Family, fill=Mature)) +
  geom_tile(color="black") +
  scale_fill_continuous(limits=c(0,14)) +
  scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn")) +
  geom_text(aes(label = Mature), color = "black", size = 5) +
  coord_fixed() + theme_minimal()
