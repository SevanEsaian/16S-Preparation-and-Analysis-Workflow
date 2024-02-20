#load required libraries -------
library(ggplot2)
library(vegan)
library(dplyr)
#Figure 1A visualizing photophysiology data----------
photophys <- read.csv("Giant Kelp Photophysiology Data.csv")
View(photophys)
ggplot(photophys, aes(x = SurfaceArea, y = Depth, fill = Age)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_point(size = 3, aes(color = Age), alpha = 0.7, position=position_dodge(width=0.75)) +  # Add geom_jitter for points
  scale_fill_manual(values = c("yellow3", "blue3")) +
  scale_color_manual(values=c("yellow4", "blue4")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    x = expression("Surface Area ( cm"^"2"~")"),
    fill = NULL) +
  annotate("text", x = 0.15, y = 3.5, label = "(A)", hjust = 0, vjust = 1, size = 10) +
  annotate("text", x = 80, y = 3.25, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 80, y = 2.9, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 80, y = 2.25, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 80, y = 1.9, label = "b", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 80, y = 1.25, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 80, y = 0.9, label = "a", hjust = 0, vjust = 1, size = 8) 

surfacearea_anova_result <- aov(SurfaceArea ~ AgeDepth, data = photophys)
summary(surfacearea_anova_result)
surfacearea_tukey_test <- TukeyHSD(surfacearea_anova_result)
surfacearea_tukey_test

#Figure 1B visualizing fvfm --------
ggplot(photophys, aes(x = FvFm, y = Depth, fill = Age)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_point(size = 3, aes(color = Age), alpha = 0.7, position=position_dodge(width=0.75)) +  # Add geom_jitter for points
  scale_fill_manual(values = c("yellow3", "blue3")) +
  scale_color_manual(values=c("yellow4", "blue4")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    x = expression("Maximum Quantum Yield of Photosystem II (Fv/Fm)"),
    fill = NULL) +
  annotate("text", x = 0.2, y = 3.5, label = "(B)", hjust = 0, vjust = 1, size = 10) +
  annotate("text", x = 0.7, y = 3.25, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.7, y = 2.85, label = "b", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.7, y = 2.25, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.7, y = 1.85, label = "b", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.7, y = 1.25, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.7, y = 0.85, label = "c", hjust = 0, vjust = 1, size = 8)


fvfm_anova_result <- aov(FvFm ~ AgeDepth, data = photophys)
summary(fvfm_anova_result)
fvfm_tukey_test <- TukeyHSD(fvfm_anova_result)
fvfm_tukey_test

# Figure 1C chl:c ---------
ggplot(photophys, aes(x = mgChl_to_mgC, y = Depth, fill = Age)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_point(size = 3, aes(color = Age), alpha = 0.7, position=position_dodge(width=0.75)) +  # Add geom_jitter for points
  scale_fill_manual(values = c("yellow3", "blue3")) +
  scale_color_manual(values=c("yellow4", "blue4")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    x = expression("Chla:C (mg)"),
    fill = NULL) +
  annotate("text", x = 0.0005, y = 3.5, label = "(C)", hjust = 0, vjust = 1, size = 10) +
  annotate("text", x = 0.011, y = 3.25, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.011, y = 2.85, label = "b", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.011, y = 2.25, label = "c", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.011, y = 1.85, label = "b", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.0105, y = 1.25, label = "a,c", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.011, y = 0.85, label = "b", hjust = 0, vjust = 1, size = 8)
  
  
chlc_anova_result <- aov(mgChl_to_mgC ~ AgeDepth, data = photophys)
summary(chlc_anova_result)
chlc_tukey_test <- TukeyHSD(chlc_anova_result)
chlc_tukey_test

#Figure 1D C:N -------
ggplot(photophys, aes(x = CN, y = Depth, fill = Age)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_point(size = 3, aes(color = Age), alpha = 0.7, position=position_dodge(width=0.75)) +  # Add geom_jitter for points
  scale_fill_manual(values = c("yellow3", "blue3")) +
  scale_color_manual(values=c("yellow4", "blue4")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    x = expression("C:N (%)"),
    fill = NULL) +
  annotate("text", x = 4.5, y = 3.5, label = "(D)", hjust = 0, vjust = 1, size = 10) +
  annotate("text", x = 28.5, y = 3.25, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 28.5, y = 2.85, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 28.5, y = 2.25, label = "b", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 28.5, y = 1.85, label = "b", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 28.5, y = 1.25, label = "b", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 28.5, y = 0.85, label = "b", hjust = 0, vjust = 1, size = 8)

cn_anova_result <- aov(CN ~ AgeDepth, data = photophys)
summary(cn_anova_result)
cn_tukey_test <- TukeyHSD(cn_anova_result)
cn_tukey_test

#NMDS ---------
NMDS.coords <- read.csv("ReVisualizing NMDS based on Bart Suggestions.csv")
View(NMDS.coords)
ggplot(NMDS.coords, aes(x=NMDS1, y=NMDS2)) +
  geom_point(size=5, aes(color=Age, shape=Depth), alpha=0.75) +
  theme_bw() +
  scale_color_manual(values=c("yellow4", "blue4")) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  stat_ellipse(geom = "polygon",
               aes(fill = Age),  
               alpha = 0.25) +
  scale_fill_manual(values=c("yellow3", "blue3"))

#Figure 3A ASV Richness---------
DiversityMetrics <- read.csv("Diversity Metrics.csv")
View(DiversityMetrics)
ggplot(DiversityMetrics, aes(x = ASVRichness, y = Depth, fill = Age)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_point(size = 3, aes(color = Age), alpha = 0.7, position=position_dodge(width=0.75)) +  # Add geom_jitter for points
  scale_fill_manual(values = c("yellow3", "blue3")) +
  scale_color_manual(values=c("yellow4", "blue4")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  annotate("text", x = 470, y = 3.5, label = "(A)", hjust = 0, vjust = 1, size = 10) +
  annotate("text", x = 1620, y = 3.25, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 1620, y = 2.85, label = "b", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 1620, y = 2.25, label = "c", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 1620, y = 1.85, label = "b", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 1620, y = 1.25, label = "c", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 1620, y = 0.85, label = "b", hjust = 0, vjust = 1, size = 8) +
  xlab("ASV Richness")
  

asv_anova_result <- aov(ASVRichness ~ Age*Depth, data = DiversityMetrics)
summary(asv_anova_result)
asv_tukey_test <- TukeyHSD(asv_anova_result)
asv_tukey_test

#Figure 3B Pielou Evenness-----
ggplot(DiversityMetrics, aes(x = PielouEvenness, y = Depth, fill = Age)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_point(size = 3, aes(color = Age), alpha = 0.7, position=position_dodge(width=0.75)) +  # Add geom_jitter for points
  scale_fill_manual(values = c("yellow3", "blue3")) +
  scale_color_manual(values=c("yellow4", "blue4")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  annotate("text", x = 0.55, y = 3.5, label = "(B)", hjust = 0, vjust = 1, size = 10) +
  annotate("text", x = 0.84, y = 3.25, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.84, y = 2.85, label = "b", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.84, y = 2.25, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.84, y = 1.85, label = "b", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.84, y = 1.25, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.84, y = 0.85, label = "c", hjust = 0, vjust = 1, size = 8) +
  xlab("Pielou's Evenness")

pe_anova_result <- aov(PielouEvenness ~ Age*Depth, data = DiversityMetrics)
summary(pe_anova_result)
pe_tukey_test <- TukeyHSD(pe_anova_result)
pe_tukey_test

#Figure 3C Betadisper --------
ggplot(DiversityMetrics, aes(x = BetaDisper, y = Depth, fill = Age)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_point(size = 3, aes(color = Age), alpha = 0.7, position=position_dodge(width=0.75)) +  # Add geom_jitter for points
  scale_fill_manual(values = c("yellow3", "blue3")) +
  scale_color_manual(values=c("yellow4", "blue4")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  xlab("Beta Dispersion") +
  annotate("text", x = 0.15, y = 3.5, label = "(C)", hjust = 0, vjust = 1, size = 10) +
  annotate("text", x = 0.62, y = 3.25, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.62, y = 2.85, label = "b", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.62, y = 2.25, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.61, y = 1.85, label = "a,b", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.62, y = 1.25, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.62, y = 0.85, label = "b", hjust = 0, vjust = 1, size = 8)



bd_anova_result <- aov(BetaDisper ~ Age*Depth, data = DiversityMetrics)
summary(bd_anova_result)
bd_tukey_test <- TukeyHSD(bd_anova_result)
bd_tukey_test

#Supp Figure Simpson Diversity ----------
ggplot(DiversityMetrics, aes(x = SimpsonEvenness, y = Depth, fill = Age)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_point(size = 3, aes(color = Age), alpha = 0.7, position=position_dodge(width=0.75)) +  # Add geom_jitter for points
  scale_fill_manual(values = c("yellow3", "blue3")) +
  scale_color_manual(values=c("yellow4", "blue4")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) + xlab("Simpson Diversity") +
  annotate("text", x = 0.013, y = 3.5, label = "(A)", hjust = 0, vjust = 1, size = 10) +
  annotate("text", x = 0.15, y = 3.25, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.15, y = 2.85, label = "b", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.145, y = 2.25, label = "a,c", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.15, y = 1.85, label = "b", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.15, y = 1.25, label = "c", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 0.145, y = 0.85, label = "b,c", hjust = 0, vjust = 1, size = 8)


simpson_anova_result <- aov(SimpsonEvenness ~ Age*Depth, data = DiversityMetrics)
summary(simpson_anova_result)
simpson_tukey_test <- TukeyHSD(simpson_anova_result)
simpson_tukey_test

#Supp Figure Shannon Diversity --------
ggplot(DiversityMetrics, aes(x = ShannonDiversity, y = Depth, fill = Age)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_point(size = 3, aes(color = Age), alpha = 0.7, position=position_dodge(width=0.75)) +  # Add geom_jitter for points
  scale_fill_manual(values = c("yellow3", "blue3")) +
  scale_color_manual(values=c("yellow4", "blue4")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) + xlab("Shannon Diversity") +
  annotate("text", x = 3.3, y = 3.5, label = "(B)", hjust = 0, vjust = 1, size = 10) +
  annotate("text", x = 6.1, y = 3.25, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 6.1, y = 2.85, label = "b", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 6.1, y = 2.25, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 6.1, y = 1.85, label = "b", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 6.1, y = 1.25, label = "a", hjust = 0, vjust = 1, size = 8) +
  annotate("text", x = 6.1, y = 0.85, label = "b", hjust = 0, vjust = 1, size = 8)


shannon_anova_result <- aov(ShannonDiversity ~ Age*Depth, data = DiversityMetrics)
summary(shannon_anova_result)
shannon_tukey_test <- TukeyHSD(shannon_anova_result)
shannon_tukey_test

#CCA Figure 4 -----
ASV.coord <- read.csv("CCA ASV Coordinates.csv")
View(ASV.coord)
Photo.coords <- read.csv("CCA Physiological Coordinates.csv")
View(Photo.coords)
ggplot(ASV.coord, aes(x=CCA1, y=CCA2)) +
  geom_point(aes(color=Age, shape=Depth), size=4) +
  scale_color_manual(values=c("yellow3", "blue4")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  geom_segment(data=Photo.coords,
               aes(x=0, xend=CCA1, y=0, yend=CCA2), color="gray", alpha=0.5) +
  geom_text(data=Photo.coords, aes(x=CCA1, y=CCA2, label=Physiology), size=5) +
  xlab("CCA1 (33.3%)") + ylab("CCA2 (21.6%)")


#PERMANOVA
data <- data.frame(
  Sample1 = c(10, 5, 0, 20),
  Sample2 = c(15, 7, 2, 18),
  # ... add more samples (rows) and taxa (columns)
  row.names = c("TaxonA", "TaxonB", "TaxonC", "TaxonD")
)
View(data)

group <- data.frame(
  Groups = c("Group1", "Group2", "Group1", "Group2")  # Assigning group names to samples
)
View(group)

#PERMANOVA ------
ASV_ra <- read.csv("ASV Relative Abundance for PERMANOVA.csv")
View(ASV_ra)
Samp_Info <- read.csv("Sample Information.csv")
View(Samp_Info)
merged_ASVra_SampInfo <- merge(Samp_Info, ASV_ra, by = "UniqueID")
View(merged_ASVra_SampInfo)

#AgeDepth
subset_AgeDepth_data <- merged_ASVra_SampInfo[, c(4, 6:7708), drop = FALSE]
View(subset_data)
subset_AgeDepth_data$AgeDepth <- as.factor(subset_data$AgeDepth)
data_for_permanova_AgeDepth <- subset_AgeDepth_data[, -1]
subset_AgeDepth_data$AgeDepth <- as.factor(subset_AgeDepth_data$AgeDepth)
permanova_result_AgeDepth <- adonis(data_for_permanova_AgeDepth ~ AgeDepth, data = subset_AgeDepth_data)
permanova_result_AgeDepth$aov.tab

#SurfaceSubsurface
subset_SurfSubSurf_data <- merged_ASVra_SampInfo[, c(5, 6:7708), drop = FALSE]
View(subset_SurfSubSurf_data)
subset_SurfSubSurf_data$Surface_Subsurface <- as.factor(subset_SurfSubSurf_data$Surface_Subsurface)
data_for_permanova_SurfSubSurf <- subset_SurfSubSurf_data[, -1]
subset_SurfSubSurf_data$Surface_Subsurface <- as.factor(subset_SurfSubSurf_data$Surface_Subsurface)
permanova_result_SurfSubSurf <- adonis(data_for_permanova_SurfSubSurf ~ Surface_Subsurface, data = subset_SurfSubSurf_data)
permanova_result_SurfSubSurf$aov.tab

#Depth
subset_Depth_data <- merged_ASVra_SampInfo[, c(2, 6:7708), drop = FALSE]
View(subset_Depth_data)
subset_Depth_data$Depth <- as.factor(subset_Depth_data$Depth)
data_for_permanova_Depth <- subset_Depth_data[, -1]
View(data_for_permanova_Depth)
subset_Depth_data$Depth <- as.factor(subset_Depth_data$Depth)
View(subset_Depth_data)
permanova_result_Depth <- adonis(data_for_permanova_Depth ~ Depth, data = subset_Depth_data)
permanova_result_Depth$aov.tab

#Age
subset_Age_data <- merged_ASVra_SampInfo[, c(3, 6:7708), drop = FALSE]
View(subset_Age_data)
subset_Age_data$Age <- as.factor(subset_Age_data$Age)
data_for_permanova_Age <- subset_Age_data[, -1]
View(data_for_permanova_Age)
subset_Age_data$Age <- as.factor(subset_Age_data$Age)
View(subset_Age_data)
permanova_result_Age <- adonis(data_for_permanova_Age ~ Age, data = subset_Age_data)
permanova_result_Age$aov.tab

#Water
KvW <- read.csv("Kelp vs Water Relative Abundance.csv")
View(KvW)
KvW$Source <- as.factor(KvW$Source)
data_for_permanova_KvW <- KvW[, -1]
permanova_result_KvW <- adonis(data_for_permanova_KvW ~ Source, data = KvW)
permanova_result_KvW$aov.tab

#Physiology
Photo.PERM <- read.csv("Photophysiology and Relative Abundance for PERMANOVA.csv")
View(Photo.PERM)


perm_data_for_permanova <- Photo.PERM[, c("Age", "SurfaceArea", colnames(Photo.PERM)[9:ncol(Photo.PERM)])]
perm_data_for_permanova$Age <- as.factor(perm_data_for_permanova$Age)
permanova_result <- adonis(
  as.dist(perm_data_for_permanova[, 3:ncol(perm_data_for_permanova)]),
  perm_data_for_permanova[, c("Age", "SurfaceArea")],
  permutations = 999
)
summary(permanova_result)
#
#ANOSIM--------
#AgeDepth
View(subset_AgeDepth_data)
dist_matrix_AgeDepth <- vegdist(subset_AgeDepth_data[, -1], method = "bray")
anosim_result_AgeDepth <- anosim(dist_matrix_AgeDepth, subset_AgeDepth_data$AgeDepth)
print(anosim_result_AgeDepth)

#Surface vs Subsurface
View(subset_SurfSubSurf_data)
dist_matrix_SurfSubSurf <- vegdist(subset_SurfSubSurf_data[, -1], method = "bray")
anosim_result_SurfSubSurf <- anosim(dist_matrix_SurfSubSurf, subset_SurfSubSurf_data$Surface_Subsurface)
print(anosim_result_SurfSubSurf)

#Depth
View(subset_Depth_data)
dist_matrix_Depth <- vegdist(subset_Depth_data[, -1], method = "bray")
anosim_result_Depth <- anosim(dist_matrix_Depth, subset_Depth_data$Depth)
print(anosim_result_Depth)

#Age
View(subset_Age_data)
dist_matrix_Age <- vegdist(subset_Age_data[, -1], method = "bray")
anosim_result_Age <- anosim(dist_matrix_Age, subset_Age_data$Age)
print(anosim_result_Age)

#Water
dist_matrix_KvW <- vegdist(KvW [, -1], method="bray")
anosim_result_KvW <- anosim(dist_matrix_KvW, KvW$Source)
print(anosim_result_KvW)

#Photophsyiology Surface Area
Photo.PERM <- read.csv("NMDS coords for photophysiology and permanova.csv")
View(Photo.PERM)
result_perm_SurfaceArea <- adonis(cbind(NMDS1, NMDS2) ~ Age * Depth * SurfaceArea, data = Photo.PERM)
result_perm_SurfaceArea$aov.tab

#FvFm
result_perm_FvFm <- adonis(cbind(NMDS1, NMDS2) ~ Age * Depth * FvFm, data = Photo.PERM)
result_perm_FvFm$aov.tab

#ChlC
result_perm_ChlC <- adonis(cbind(NMDS1, NMDS2) ~ Age * Depth * ChlaC, data = Photo.PERM)
result_perm_ChlC$aov.tab

#CN
result_perm_CN <- adonis(cbind(NMDS1, NMDS2) ~ Age * Depth * CN, data = Photo.PERM)
result_perm_CN$aov.tab
#
#PERMDISP ------
#AgeDepth
dist_matrix_AgeDepth <- vegdist(subset_AgeDepth_data[, -1], method = "bray")
permdisp_result_AgeDepth <- betadisper(dist_matrix_AgeDepth, subset_AgeDepth_data$AgeDepth)
permutation_test_result_AgeDepth <- permutest(permdisp_result_AgeDepth, pairwise = TRUE)
permutation_test_result_AgeDepth$tab

#Surface vs Subsurface
dist_matrix_SurfSubSurf <- vegdist(subset_SurfSubSurf_data[, -1], method = "bray")
permdisp_result_SurfSubSurf <- betadisper(dist_matrix_SurfSubSurf, subset_SurfSubSurf_data$Surface_Subsurface)
permutation_test_result_SurfSubSurf <- permutest(permdisp_result_SurfSubSurf, pairwise = TRUE)
permutation_test_result_SurfSubSurf$tab

#Depth
dist_matrix_Depth <- vegdist(subset_Depth_data[, -1], method = "bray")
permdisp_result_Depth <- betadisper(dist_matrix_Depth, subset_Depth_data$Depth)
permutation_test_result_Depth <- permutest(permdisp_result_Depth, pairwise = TRUE)
permutation_test_result_Depth$tab

#Age
dist_matrix_Age <- vegdist(subset_Age_data[, -1], method = "bray")
permdisp_result_Age <- betadisper(dist_matrix_Age, subset_Age_data$Age)
permutation_test_result_Age <- permutest(permdisp_result_Age, pairwise = TRUE)
permutation_test_result_Age$tab

#Water
dist_matrix_KvW <- vegdist(KvW[, -1], method="bray")
permdisp_result_KvW <- betadisper(dist_matrix_KvW, KvW$Source)
permutation_test_result_KvW <- permutest(permdisp_result_KvW, pairwise = TRUE)
permutation_test_result_KvW$tab

#Principal Component Analysis ------
#Unifrac Weighted Principal Component Analysis
PC.coord <- read.csv("Principal Component - Not Scaled - Coordinates.csv")
View(PC.coord)
ggplot(PC.coord, aes(x=Comp.1, y=Comp.2)) +
  geom_point(size=5, aes(color=Age, shape=Depth), alpha=0.75) +
  theme_bw() +
  scale_color_manual(values=c("yellow3", "blue4")) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  xlab("PC1") + ylab("PC2") +
  annotate("text", x = -13, y = 10, label = "(A)", hjust = 0, vjust = 1, size = 10)

PCoA.coord <- read.csv("PCoA Coordinates of Percent Relative Abundance.csv")
View(PCoA.coord)
ggplot(PCoA.coord, aes(x=PCoA1, y=PCoA2)) +
  geom_point(size=5, aes(color=Age, shape=Depth), alpha=0.75) +
  theme_bw() +
  scale_color_manual(values=c("yellow3", "blue4")) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  annotate("text", x = -1, y = 1.75, label = "(B)", hjust = 0, vjust = 1, size = 10)

