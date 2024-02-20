#load packages -------------
library(ape)
library(permute)
library(lattice)
library(vegan)
library(readr)
library(tidyverse)
library(ggplot2)
library(GGally)
library(CCA)
library(easyGgplot2)
library(ggpubr)
#Rarefied Abundance PCoA Run of T1 and T2 blade microbiomes and water samples ---------------
Rar.Mic <- read_csv("Rarefied Microbiome Counts for PCoA Setup.csv")
View(Rar.Mic)

mic.dst <- dist(Rar.Mic[,3:12724])
mic.dst
mic.bd <- betadisper(mic.dst, Rar.Mic$SampleCoords)
mic.bd
anova(mic.bd)
permutest(mic.bd)
labs <- paste("Dimension", 3:12724, "(", 
              round(100*mic.bd$eig / sum(mic.bd$eig), 2), "%)")

plot(mic.bd, cex=2, pch=15:17,
     main="Rarefied Bacterial Counts: Giant Kelp and Water Samples", cex.lab=1.25,
     xlab=labs[1], ylab=labs[2],
     hull=FALSE, ellipse=TRUE, conf=0.68, lwd=2)

boxplot(mic.bd, xlab="Site",  col=c("gray", "red", "green"))

plot(mic.bd, main="Giant Kelp", 
     hull=FALSE, ellipse=TRUE, conf=0.5)

#Percent Relative Abundance PCoA Trial ---------------
RA.Mic <- read_csv("Percent Relative Abundance Counts for PCoA Setup.csv")
View(RA.Mic)

RA.dst <- dist(RA.Mic[,2:7824])
RA.dst
RA.bd <- betadisper(RA.dst, RA.Mic$SampleCoords)
RA.bd
anova(RA.bd)
permutest(RA.bd)
labs <- paste("Dimension", 2:7824, "(", 
              round(100*RA.bd$eig / sum(RA.bd$eig), 2), "%)")

plot(RA.bd, cex=2, pch=15:17,
     main="Percent Relative Abundance: Giant Kelp and Water Samples", cex.lab=1.25,
     xlab=labs[1], ylab=labs[2],
     hull=FALSE, ellipse=TRUE, conf=0.68, lwd=2)

boxplot(RA.bd, xlab="Site",  col=c("gray", "red", "green"))

plot(RA.bd, main="Giant Kelp", 
     hull=FALSE, ellipse=TRUE, conf=0.5)
RA.PCoA.df <- data.frame(Distance_to_centroid=RA.bd$vectors, Group=RA.bd$group)
View(RA.PCoA.df)
write.csv(RA.PCoA.df, "~/Kelp/Kelp_Data/16S_Library/July 22 2022 - Re Doing Some Figures for Prospectus and Paper\\PCoA Coordinates of Percent Relative Abundance.csv",row.names=TRUE)

Norm.df.PCoA <- data.frame(PCoA = Norm.bd$vectors, Group=Norm.bd$group)

#PCoA - Percent Relative Abundance - Visualization -----------------
PCoAData <- read_csv("PCoA Coordinates Water Excluded of Percent Relative Abundance.csv")
View(PCoAData)

ggplot(PCoAData, aes(x=PCoA1, y=PCoA2, shape=Age, color=Depth)) +
  geom_point(size=7, alpha=0.9) +
  stat_ellipse(type="norm", linetype=1) +
  theme_classic() +
  labs(x = "PCoA1  (32.0% Variance)", y = "PCoA2  (16.5% Variance)") +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  theme(text=element_text(size=18)) +
  labs(title = "Giant Kelp Microbiome % Relative Abundance") +
  theme(plot.title=element_text(hjust=0.5))

#NMDS Check to Make Sure My NMDS' all Make Sense--------------
Mic.NMDS<- read_csv("Percent Relative Abundance Counts for NMDS SetUp.csv")
View(Mic.NMDS)
Mic.NMDS.Run <- metaMDS(Mic.NMDS, distance="bray")
Mic.NMDS.Run
plot(Mic.NMDS.Run)
Mic.NMDS.data.scores = as.data.frame(scores(Mic.NMDS.Run))
View(Mic.NMDS.data.scores)
write.csv(Mic.NMDS.data.scores, "~/Kelp/Kelp_Data/16S_Library/July 22 2022 - Re Doing Some Figures for Prospectus and Paper\\NMDS Coordinates of Percent Relative Abundance.csv",row.names=TRUE)
#Plotting NMDS ---------------------------
NMDS.Values <- read_csv("NMDS Values no water from Percent Relative Abundance.csv")
View(NMDS.Values)
NMDS.CF <- subset(NMDS.Values, SampleCoords=="Bottom T1 Juvenile" |
                    SampleCoords=="Bottom T2 Mature" |
                    SampleCoords=="Middle T1 Juvenile" |
                    SampleCoords=="Middle T2 Mature" |
                    SampleCoords=="Surface T1 Juvenile" |
                    SampleCoords=="Surface T2 Mature")
View(NMDS.CF)
ggplot(NMDS.CF, aes(x=NMDS1, y=NMDS2, shape=Depth, color=Age)) +
  geom_point(size=5, alpha=0.5) +
  theme_classic() +
  scale_color_manual(values = c("#CC6600", "#3300CC")) +
  theme(text=element_text(size=18)) +
  labs(title = "Giant Kelp Microbiome % Relative Abundance") +
  theme(plot.title=element_text(hjust=0.5)) +
  xlim(-0.7, 1) + ylim(-1.3, 1.3) +
  theme_bw()

#Subsetting just mature blades
NMDS.Mature <- subset(NMDS.Values, Age=="Mature")
View(NMDS.Mature)
Mat.NMDS <- ggplot(NMDS.Mature, aes(x=NMDS1, y=NMDS2, shape=Depth, color=Age)) +
  geom_point(size=5, alpha=0.7) +
  theme_bw() +
  scale_color_manual(values = c("#3300CC")) +
  theme(text=element_text(size=14)) +
  theme(plot.title=element_text(hjust=0.5)) +
  xlim(-0.7, 1) + ylim(-1.3, 1.3) +
  theme(legend.position='none') +
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic", hjust=0.5),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  )+
  labs(title="Mature Microbiomes",
       x ="NMDS1", y = "NMDS2") 


NMDS.Juvenile <- subset(NMDS.CF, Age =="Juvenile")
View(NMDS.Juvenile)

Juv.NMDS <- ggplot(NMDS.Juvenile, aes(x=NMDS1, y=NMDS2, shape=Depth, color=Age)) +
  geom_point(size=5, alpha=0.7) +
  theme_bw() +
  scale_color_manual(values = c("#CC6600")) +
  theme(text=element_text(size=14)) +
  theme(plot.title=element_text(hjust=0.5)) +
  xlim(-0.7, 1) + ylim(-1.3, 1.3)+
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic", hjust=0.5),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  )+
  labs(title="Juvenile Microbiomes",
       x ="NMDS1", y = "NMDS2") +
  theme(legend.position='none')


ggarrange(Juv.NMDS, Mat.NMDS, ncol=2, nrow=1)

#Plotting NMDS with water
All.NMDS <- read_csv("NMDSCoordinates of Percent Relative Abundance.csv")
View(as.data.frame(All.NMDS))
ggplot(All.NMDS, aes(x=NMDS1, y=NMDS2, color=Source)) +
  geom_point(size=7, alpha=0.7) +
  theme_classic() +
  theme(text=element_text(size=20)) +
  theme(plot.title=element_text(hjust=0.5)) +
  theme_bw() +
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic", hjust=0.5),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  ) 

#NMDS Photophysiology Specific ------------------
Phot.NMDS <- read_csv("Use for NMDS for Photophysiology Data Instead of Microbiome.csv")
View(Phot.NMDS)

Phot.NMDS.Run <- metaMDS(Phot.NMDS, distance="bray")
Phot.NMDS.Run
plot(Phot.NMDS.Run)
Phot.NMDS.data.scores = as.data.frame(scores(Phot.NMDS.Run))
View(Phot.NMDS.data.scores)  
write.csv(Phot.NMDS.data.scores, "~/Kelp/Kelp_Data/16S_Library/July 22 2022 - Re Doing Some Figures for Prospectus and Paper\\NMDS Coordinates of Photophysiology Characteristics.csv",row.names=TRUE)

Phot.NMDS.Plot <- read_csv("Photophysiology NMDS Coordinates for Plotting.csv")
View(Phot.NMDS.Plot)

Phot.NMDS.Plot.CF <- subset(Phot.NMDS.Plot, SampleCoords=="Bottom T1 Juvenile" |
                              SampleCoords == "Bottom T2 Mature" |
                              SampleCoords == "Middle T1 Juvenile" |
                              SampleCoords == "Middle T2 Mature" |
                              SampleCoords == "Surface T1 Juvenile" |
                              SampleCoords == "Surface T2 Mature")
ggplot() +
  geom_point(data=Phot.NMDS.Plot.CF, aes(x=NMDS1, y=NMDS2, color=Age, shape=Depth, alpha=0.7, size=6)) +
  scale_color_manual(values=c("#CC6600", "#3300CC"))+ 
  theme_classic() +
  theme(text=element_text(size=18)) +
  theme(plot.title=element_text(hjust=0.5)) 

#ANOSIM analysis ----------------------------------------
library(permute)
library(lattice)
library(vegan)
?anosim()
Microbiome <- read.csv("Microbiome for ANOSIM.csv")
Microbiome.env <- read.csv("Microbiome Envs for ANOSIM.csv")
View(Microbiome.env)
Microbiome.dist <- vegdist(Microbiome)
Microbiome.ano <- with(Microbiome.env, anosim(Microbiome.dist, Source))
summary(Microbiome.ano)
plot(Microbiome.ano)
#Water vs Kelp (R = 0.942, P = 0.001)
#ANOSIM of just blades mature and juvenile-----------------
CF.Microbiome <- read.csv("Cohort Follow Blades Only Microbiome for ANOSIM.csv")
CF.Microbiome.env <- read.csv("Cohort Follow Blades Only Microbiome Envs for ANOSIM.csv")
CF.Microbiome.env
View(CF.Microbiome.env)
CF.Microbiome.dist <- vegdist(CF.Microbiome)
CF.Microbiome.ano <- with(CF.Microbiome.env, anosim(CF.Microbiome.dist, Age))
summary(CF.Microbiome.ano)
plot(CF.Microbiome.ano)        
#all distinct blades and ages is R = 0.66, p =0.001
#Mature Blades Only ANOSIM
Mat.Microbiome <- read.csv("Mature Blades Only Microbiome for ANOSIM.csv")
Mat.Microbiome.env <- read.csv("Mature Blades Only Microbiome Envs for ANOSIM.csv")
Mat.Microbiome.env
Mat.Microbiome.dist <- vegdist(Mat.Microbiome)
Mat.Microbiome.ano <- with(Mat.Microbiome.env, anosim(Mat.Microbiome.dist, DepthBin))
summary(Mat.Microbiome.ano)
plot(Mat.Microbiome.ano)

#ANOSIM Photophysiology ----------------
Phy.Mic <- read.csv("ANOSIM NMDS for Photophysiology Data Instead of Microbiome.csv")
Phy.Mic.env <- read.csv("ANOSIM Photophysiology NMDS Coordinates.csv")
Phy.Mic.env
Phy.Microbiome.dist <- vegdist(Phy.Mic)
Phy.Microbiome.ano <- with(Phy.Mic.env, anosim(Phy.Microbiome.dist, SampleCoords))
summary(Phy.Microbiome.ano)
plot(CF.Microbiome.ano)   
#CCA Analysis ReDO------
install.packages("Envfit")
library(Envfit)
cca.data <- read_csv("SetUp For CCA with Percent Relative Abundances.csv")
View(cca.data)
cca.data %>% remove_rownames %>% column_to_rownames(var="UniqueID")
cca.data2 <- column_to_rownames(cca.data, "UniqueID")
View(cca.data2)
CCA_Env_Data <- read_csv("CCA Env Data SetUp.csv")
View(CCA_Env_Data)
CCA_Env_Data %>% remove_rownames %>% column_to_rownames(var="UI")
CCA_Env_Data2 <- column_to_rownames(CCA_Env_Data, "UI")
View(CCA_Env_Data2)

RA.cca <- cca(cca.data2 ~ PC + PN + CN + SurfaceArea + mgChl_to_mgC + FvFm, data=CCA_Env_Data2)
str(RA.cca)
plot(RA.cca)
RA.cca
View(RA.cca)
library(mgcv)
library(nlme)

ccavectors <- as.matrix(scores(RA.cca, display = "bp", scaling = "species")*2.34713) %>%
  as.data.frame()
summary(ccavectors)
View(ccavectors)
write.csv(ccavectors, "~/Kelp/Kelp_Data/16S_Library/July 22 2022 - Re Doing Some Figures for Prospectus and Paper\\EnvironmentalVectorsForCCA.csv",row.names=TRUE)

NMDS.CCA.Run <- metaMDS(cca.data2, distance="bray")
NMDS.CCA.data.scores = as.data.frame(scores(NMDS.CCA.Run))
View(NMDS.CCA.data.scores)
write.csv(NMDS.CCA.data.scores, "~/Kelp/Kelp_Data/16S_Library/July 22 2022 - Re Doing Some Figures for Prospectus and Paper\\BladeNMDSCoordinates.csv",row.names=TRUE)


#Plotting CCA
BladeNMDSCoordinates_NoT2Juvenile <- read_csv("BladeNMDSCoordinates_NoT2Juvenile.csv")
View(BladeNMDSCoordinates_NoT2Juvenile)
EnvironmentalVectorsForCCA <- read_csv("EnvironmentalVectorsForCCA.csv")
View(EnvironmentalVectorsForCCA)
EnvironmentalVectorsForCCA

ggplot() +
  # arrows for species
  geom_segment(data = EnvironmentalVectorsForCCA, aes(x = 0, y = 0, xend = CCA1, yend = CCA2), alpha = .9, arrow = arrow(length = unit(0, "cm"))) +
  # samples
  geom_point(data = BladeNMDSCoordinates_NoT2Juvenile, aes(x = NMDS1, y = NMDS2, color = depth, shape=age), alpha=0.9, size=6) +
  # label arrows
  # labels
  labs(x = "CCA1",
       y = "CCA2") +
  theme_classic() +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))

##For the histogram - raw reads setup---------------------
Hist.data <- read.csv("ForHistogram - RawReadsSetup.csv")
View(Hist.data)
Hist.data.blades <- subset(Hist.data, sampletype=="Juvenile" | sampletype=="Mature")
View(Hist.data.blades)
ggplot2.histogram(data=Hist.data.blades, xName='RawReadsSum',
                  groupColors=c("#CC6600", "#3300CC"), 
                  groupName='sampletype', legendPosition="top",
                  alpha=0.1, binwidth=5000) +
  theme_bw() +
  labs(title="Giant Kelp Blades",
       x ="Reads per Sample", y = "count") +
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic", hjust=0.5),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  ) +
  labs(color = "Age") +
  labs(fill = "Age") +
  theme(axis.text = element_text(size = 14)) 
  
Hist.data.water <- subset(Hist.data, sampletype=="1" | sampletype=="2")
View(Hist.data.water)
Hist.data.water <- as.data.frame(Hist.data.water)
ggplot2.histogram(data=Hist.data.water, xName='RawReadsSum',
                  groupColors=c("#CC6600", "#3300CC"), 
                  groupName='sampletype', legendPosition="top",
                  alpha=0.1, binwidth=5000) +
  theme_bw() +
  labs(title="Seawater Microbial Community",
       x ="Reads per Sample", y = "count") +
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic", hjust=0.5),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  ) +
  labs(color = "Timepoint") +
  labs(fill = "Timepoint") +
  theme(axis.text = element_text(size = 14)) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE))


##Calculate Faith's phylogenetic diversity using picante ----------
remotes::install_github("ropensci/phylocomr")
library(picante)
library(phylocomr)
library(ape)

ph_ecovolve(speciation = 0.05, extinction = 0.005, time_units = 50)

taxa_file <- system.file("examples/taxa", package = "phylocomr")
View(as.data.frame(taxa_file))
phylo_file <- system.file("examples/phylo", package = "phylocomr")
(taxa_str <- readLines(taxa_file))
