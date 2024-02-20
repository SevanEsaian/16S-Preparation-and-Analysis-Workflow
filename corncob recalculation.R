#Load packages --------
library(phyloseq)
library(corncob)
library(dplyr)
library(tidyr)
library(ggplot2)

#Upload necessary files and run corncob---------
metadata <- read.csv("metadata_juvenileandmature.csv")
View(metadata)
speciescoords <- read.csv("species_coords.csv")
View(speciescoords)
rawreads <- read.csv("ASVRawReads with Euks and NA removed.csv")
View(rawreads)
seqtab <- rawreads %>% tibble::column_to_rownames("UniqueID")
seqtab.matrix <- as.matrix(seqtab)
View(seqtab.matrix)
seqtab.matrix.table <- otu_table(seqtab.matrix, taxa_are_rows = FALSE)
seqtab.matrix.table

taxa <- speciescoords %>% tibble::column_to_rownames("ASV")
taxa.matrix <- as.matrix(taxa)
taxa.matrix.table <- tax_table(taxa.matrix)
View(taxa.matrix.table)
taxa.genus <- speciescoords %>% tibble::column_to_rownames("genus")

metadata <- metadata %>% tibble::column_to_rownames("UniqueID")
metadata.sample <- sample_data(metadata)
metadata.sample

phylo.object <- phyloseq(seqtab.matrix.table, metadata.sample, taxa.matrix.table)
phylo.object

phylo.object.genus <- phylo.object %>%
  tax_glom("genus")
phylo.object.genus
Taxa.sheet <- as.data.frame(tax_table(phylo.object.genus))
View(Taxa.sheet)
subset(Taxa.sheet, genus=="Psuedoalteromos")
taxa_table(phylo.object.genus)

diff.ab.age <- differentialTest(formula = ~ age,
                                phi.formula = ~ age,
                                formula_null = ~1, 
                                phi.formula_null = ~ age,
                                data=phylo.object,
                                test="Wald", fdr_cutoff=0.05)
diff.ab.age.genus <- differentialTest(formula = ~ age,
                                phi.formula = ~ age,
                                formula_null = ~1, 
                                phi.formula_null = ~ age,
                                data=phylo.object.genus,
                                test="Wald", fdr_cutoff=0.05)
plot(diff.ab.age)
plot(diff.ab.age.genus)
diff.ab.age.genus.plotting <- bbdml(formula = ASV8166 ~ age, phi.formula= ~ age,
                        data=phylo.object.genus)
summary(diff.ab.age.genus.plotting)
plot(diff.ab.age.genus.plotting, color="age")

plot(diff.ab.age.genus, level=c("genus")) 

plot(depth.ctrlage.genus) +
  theme(legend.position = "none") +
  ylab(NULL)
diff.ab.age.genus$significant_taxa
summary(diff.ab.age.genus)  

otu_to_taxonomy(OTU = diff.ab.age.genus$significant_taxa, data=phylo.object.genus)

diff.ab.depth <- differentialTest(formula = ~ depth,
                                  phi.formula = ~ depth,
                                  formula_null = ~1, 
                                  phi.formula_null = ~ depth,
                                  data=phylo.object,
                                  test="Wald", fdr_cutoff=0.05)
diff.ab.depth.genus <- differentialTest(formula = ~ depth,
                                  phi.formula = ~ depth,
                                  formula_null = ~1, 
                                  phi.formula_null = ~ depth,
                                  data=phylo.object.genus,
                                  test="Wald", fdr_cutoff=0.05)
plot(diff.ab.depth.genus, level=c("genus")) 
diff.ab.depth.genus.plotting <- bbdml(formula = ASV11861~ depth, phi.formula= ~ depth,
                                    data=phylo.object.genus)

summary(diff.ab.depth.genus.plotting)

otu_to_taxonomy(OTU = diff.ab.depth.genus$significant_taxa, data=phylo.object.genus)
diff.ab.depth.genus$p_fdr
diff.ab.depth.genus$significant_models
diff.ab.depth.genus$significant_taxa



Pseudoalteromos.age <- bbdml(formula=Pseudoalteromos ~ age, phi.formula = ~age, 
                             data=phylo.object)


#calculate average value per genera and cross reference to corncob ---------
RA <- read.csv("Relative_Abundance.csv")
View(RA)
taxa <- read.csv("species_coords.csv")
View(taxa)
RA.taxa <- merge(taxa, RA, by="ASV")
View(RA.taxa)
#write.csv(RA.taxa, "ASV relative abundance species coordinates.csv")
RA.taxa.gather <- tidyr::gather(RA.taxa, "UniqueID", "RelativeAbundance", 8:93)
View(RA.taxa.gather)
RA.taxa.gather.genus <- RA.taxa.gather[c(-1, -2, -3, -4, -5, -6)]
View(RA.taxa.gather.genus)

RA.taxa.genus.sum <- aggregate(RA.taxa.gather.genus$RelativeAbundance, list(RA.taxa.gather.genus$genus), FUN=sum)
View(RA.taxa.genus.sum)

RA.taxa.sum.spread <- tidyr::spread(RA.taxa.gather.genus, "UniqueID", "RelativeAbundance")
genus.test.pivotwider <- RA.taxa.gather.genus %>%
  tidyr::pivot_wider(names_from = "genus",
                     values_from="RelativeAbundance",
                     values_fn=sum)
#write.csv(genus.test.pivotwider, "Genus summed RA per sample for prevlance.csv")
?spread
View(RA.taxa.sum.spread)

RA.taxa.genus.average <- aggregate(RA.taxa.gather.genus$RelativeAbundance, list(RA.taxa.gather.genus$genus), FUN=mean)
View(RA.taxa.genus.average)
##Subset desired genera separately --------------
View(RA.taxa.gather.sub)
RA.taxa.gather.sub <- subset(RA.taxa.gather, genus=="Pseudoalteromos" |
                        genus=="Psychromos" |
                        genus=="Blastopirellula" |
                        genus=="OM43_clade" |
                        genus=="Rubripirellula" |
                        genus=="Lutimos" |
                        genus=="Synechococcus_CC9902" |
                        genus=="Roseibacillus" |
                        genus=="Lentimos" |
                        genus=="Flavicella" |
                        genus=="Alcanivorax" |
                        genus=="Stenotrophomos" |
                        genus=="Halomos" |
                        genus=="Flavicella"| 
                        genus=="Lewinella" )
View(RA.taxa.gather.sub)
RA.taxa.gather.sub.genus <- RA.taxa.gather.sub[c(-1, -2, -3, -4, -5, -6)]
View(RA.taxa.gather.sub.genus)

RA.taxa.gather.sub.genus.sum <- RA.taxa.gather.sub.genus %>% group_by(genus,UniqueID) %>% 
  summarise(RelativeAbundance=sum(RelativeAbundance),
            .groups = 'drop') %>%
  as.data.frame()
View(RA.taxa.gather.sub.genus)
View(RA.taxa.gather.sub.genus.sum)
RA.taxa.gather.genus.spread.sum <- tidyr::spread(RA.taxa.gather.sub.genus.sum, "genus", "RelativeAbundance")
View(RA.taxa.gather.genus.spread.sum)

#spread the genera and merge it with metadata
RA.taxa.genus.spread <- tidyr::spread(RA.taxa.gather.sub.genus.sum, "genus", "RelativeAbundance")
View(RA.taxa.genus.spread)
#write.csv(RA.taxa.genus.spread, "summed spread subset genera relative abundances.csv")


#check all asvs for the hits you have found find the ones that are significant--------

#Amylibacter just for genera plot ASVs aren't significant because %RA is so low--------
Amylibacter.depth <- bbdml(formula=ASV5404 ~ depth, phi.formula = ~depth,  data=phylo.object.genus)
summary(Amylibacter.depth)
plot(Amylibacter.depth, color="depth")
#ASV5404
ASV5404.depth <- bbdml(formula=ASV5404 ~ depth, phi.formula = ~depth, data=phylo.object)
#ASV61
ASV61.depth <- bbdml(formula=ASV61 ~ depth, phi.formula = ~depth, data=phylo.object)
#
#Lewinella -------
Lewinella.depth <- bbdml(formula=ASV4540 ~ depth, phi.formula = ~depth,  data=phylo.object.genus)
summary(Lewinella.depth)
plot(Lewinella.depth, color="depth")
diff.ab.age.genus.lewinella <- bbdml(formula = ASV4540 ~ age, phi.formula= ~ age, data=phylo.object.genus)
summary(diff.ab.age.genus.lewinella)#significant abundance
plot(diff.ab.age.genus.lewinella, color="age")

#ASV101
diff.ab.age.101 <- bbdml(formula=ASV101 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.101 <- bbdml(formula=ASV101 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV1010
diff.ab.age.1010 <- bbdml(formula=ASV1010 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.1010 <- bbdml(formula=ASV1010 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV1031 #highly favorable towards mature samples and surface samples
diff.ab.age.1031 <- bbdml(formula=ASV1031 ~ age, phi.formula = ~ age, data=phylo.object)
summary(diff.ab.age.1031)
plot(diff.ab.age.1031, color="age") 
diff.ab.depth.1031 <- bbdml(formula=ASV1031 ~ depth, phi.formula = ~ depth, data=phylo.object)
summary(diff.ab.depth.1031)
plot(diff.ab.depth.1031, color="depth")
#ASV1190
diff.ab.age.1190 <- bbdml(formula=ASV1190 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.1190 <- bbdml(formula=ASV1190 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV1321
diff.ab.age.1321 <- bbdml(formula=ASV1321 ~ age, phi.formula = ~ age, data=phylo.object)
diff.depth.age.1321 <- bbdml(formula=ASV1321 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV1627
diff.ab.age.1627 <- bbdml(formula=ASV1627 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.1627 <- bbdml(formula=ASV1627 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV1878
diff.ab.age.1878 <- bbdml(formula=ASV1878 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.1878 <- bbdml(formula=ASV1878 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV2207
diff.ab.age.2207 <- bbdml(formula=ASV2207 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.2207 <- bbdml(formula=ASV2207 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV2224
diff.ab.age.2224 <- bbdml(formula=ASV2224 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.2224 <- bbdml(formula=ASV2224 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV2352
diff.ab.age.2352 <- bbdml(formula=ASV2352 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.2352 <- bbdml(formula=ASV2352 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV3018
diff.ab.age.3018 <- bbdml(formula=ASV3018 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.3018 <- bbdml(formula=ASV3018 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV3175
diff.ab.age.3175 <- bbdml(formula=ASV3175 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.3175 <- bbdml(formula=ASV3175 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV3262 # highly favorable towards mature samples and subsurface samples
diff.ab.age.3262 <- bbdml(formula=ASV3262 ~ age, phi.formula = ~ age, data=phylo.object)
summary(diff.ab.age.3262)
plot(diff.ab.age.3262, color="age")
diff.ab.depth.3262 <- bbdml(formula=ASV3262 ~ depth, phi.formula = ~ depth, data=phylo.object)
summary(diff.ab.depth.3262)
plot(diff.ab.depth.3262, color="depth")
#ASV3318
diff.ab.age.3318 <- bbdml(formula=ASV3318 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.3318 <- bbdml(formula=ASV3318 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV3411
diff.ab.age.3411 <- bbdml(formula=ASV3411 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.3411 <- bbdml(formula=ASV3411 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV3423
diff.ab.age.3423 <- bbdml(formula=ASV3423 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.3423 <- bbdml(formula=ASV3423 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV361
diff.ab.age.361 <- bbdml(formula=ASV361 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.361 <- bbdml(formula=ASV361 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV3615
diff.ab.age.3615 <- bbdml(formula=ASV3615 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.3615 <- bbdml(formula=ASV3615 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV3688
diff.ab.age.3688 <- bbdml(formula=ASV3688 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.3688 <- bbdml(formula=ASV3688 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV3780
diff.ab.age.3780 <- bbdml(formula=ASV3780 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.3780 <- bbdml(formula=ASV3780 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV4529
diff.ab.age.4529 <- bbdml(formula=ASV4529 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.4529 <- bbdml(formula=ASV4529 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV4540 # mature favorable and surface samples
diff.ab.age.4540 <- bbdml(formula=ASV4540 ~ age, phi.formula = ~ age, data=phylo.object)
summary(diff.ab.age.4540)
plot(diff.ab.age.4540, color="age")
diff.ab.depth.4540 <- bbdml(formula=ASV4540 ~ depth, phi.formula = ~ depth, data=phylo.object)
summary(diff.ab.depth.4540)
plot(diff.ab.depth.4540, color="depth")
#ASV4889
diff.ab.age.4899 <- bbdml(formula=ASV4899 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.4899 <- bbdml(formula=ASV4899 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV544
diff.ab.age.544 <- bbdml(formula=ASV544 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.544 <- bbdml(formula=ASV544 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV566
diff.ab.age.566 <- bbdml(formula=ASV566 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.566 <- bbdml(formula=ASV566 ~ depth, phi.formula = ~ depth, data=phylo.object)
summary(diff.ab.depth.566)
#ASV574
diff.ab.age.574 <- bbdml(formula=ASV574 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.574 <- bbdml(formula=ASV574 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV691
diff.ab.age.691 <- bbdml(formula=ASV691 ~ age, phi.formula = ~ age, data=phylo.object)
#ASV7029 #highly favorable towards mature and subsurface samples
diff.ab.age.7029 <- bbdml(formula=ASV7029 ~ age, phi.formula = ~ age, data=phylo.object)
summary(diff.ab.age.7029)
plot(diff.ab.age.7029, color="age")
diff.ab.depth.7029 <- bbdml(formula=ASV7029 ~ depth, phi.formula = ~ depth, data=phylo.object)
summary(diff.ab.depth.7029)
plot(diff.ab.depth.7029, color="depth")
#ASV7154 #favorable towards mature - very low %RA - and subsurface (very low %RA)
diff.ab.age.7154 <- bbdml(formula=ASV7154 ~ age, phi.formula = ~ age, data=phylo.object)
summary(diff.ab.age.7154)
plot(diff.ab.age.7154, color="age")
diff.ab.depth.7154 <- bbdml(formula=ASV7154 ~ depth, phi.formula = ~ depth, data=phylo.object)
summary(diff.ab.depth.7154)
plot(diff.ab.depth.7154, color="depth")
#ASV7206
diff.ab.age.7206 <- bbdml(formula=ASV7206 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.7206 <- bbdml(formula=ASV7206 ~ depth, phi.formula = ~ depth, data=phylo.object)
#ASV7213
diff.ab.age.7213 <- bbdml(formula=ASV7213 ~ age, phi.formula = ~ age, data=phylo.object)
#ASV7310 # significant for subsurface blades
diff.ab.age.7310 <- bbdml(formula=ASV7310 ~ age, phi.formula = ~ age, data=phylo.object)
diff.ab.depth.7310 <- bbdml(formula=ASV7310 ~ depth, phi.formula = ~ depth, data=phylo.object)
summary(diff.ab.depth.7310)
plot(diff.ab.depth.7310, color="depth")
#ASV736
diff.ab.age.736 <- bbdml(formula=ASV736 ~ age, phi.formula = ~ age, data=phylo.object)
#ASV7636
diff.ab.age.7636 <- bbdml(formula=ASV7636 ~ age, phi.formula = ~ age, data=phylo.object)
#ASV803
diff.ab.age.803 <- bbdml(formula=ASV803 ~ age, phi.formula = ~ age, data=phylo.object)
#ASV901
diff.ab.age.901 <- bbdml(formula=ASV901 ~ age, phi.formula = ~ age, data=phylo.object)

#Synechococcus_CC9902 genus level juvenile significant:------- 
Synech.age <- bbdml(formula=ASV7846 ~ age, phi.formula = ~age, data=phylo.object.genus)
summary(Synech.age)
plot(Synech.age, color="age") # genus juvenile significant
#ASV7846
diff.ab.age.7846 <- bbdml(formula=ASV7846~age, phi.formula=~age, data=phylo.object)
summary(diff.ab.age.7846)
#ASV1002 - model failed
diff.ab.age.1002 <- bbdml(formula = ASV1002 ~ age, phi.formula= ~ age, data=phylo.object)
#ASV1225 - model failed
diff.ab.age.1225 <- bbdml(formula = ASV1225 ~ age, phi.formula= ~ age, data=phylo.object)
#ASV16 - no significant
diff.ab.age.16 <- bbdml(formula = ASV16 ~ age, phi.formula= ~ age, data=phylo.object)
summary(diff.ab.age.16)
#ASV459 - model failed
diff.ab.age.459 <- bbdml(formula = ASV459 ~ age, phi.formula= ~ age, data=phylo.object)
#ASV5 - model failed
diff.ab.age.5 <- bbdml(formula = ASV5 ~ age, phi.formula= ~ age, data=phylo.object)
#ASV6547 - model failed
diff.ab.age.6547 <- bbdml(formula = ASV6547 ~ age, phi.formula= ~ age, data=phylo.object)
#ASV7672 - model failed
diff.ab.age.7672 <- bbdml(formula = ASV7672 ~ age, phi.formula= ~ age, data=phylo.object)
#ASV9 - model failed
diff.ab.age.9 <- bbdml(formula = ASV9 ~ age, phi.formula= ~ age, data=phylo.object)
#Pseudoalteromos genus juvenile significant and subsurface significant-------
Pseudo.depth <- bbdml(formula=ASV1065~age, phi.formula=~age, data=phylo.object.genus)
summary(Pseudo.depth)
plot(Pseudo.depth, color="depth")
Pseudo.age <- bbdml(formula=ASV1065 ~ age, phi.formula = ~age, data=phylo.object.genus)
summary(Pseudo.age)
plot(Synech.age, color="age")
Pseudo.depth <- bbdml(formula=ASV1065 ~ depth, phi.formula = ~depth, data=phylo.object.genus)
summary(Pseudo.depth)
plot(Pseudo.depth, color="depth")
#ASV1065 - model worked favors mature and subsurface
diff.ab.age.1065 <- bbdml(formula = ASV1065 ~ age, phi.formula= ~ age, data=phylo.object.genus)
plot(diff.ab.age.1065, color = "age", B=50)
summary(diff.ab.age.1065)
diff.ab.depth.1065 <- bbdml(formula = ASV1065 ~ depth, phi.formula= ~ depth, data=phylo.object.genus)
plot(diff.ab.depth.1065, color = "depth", B=50)
summary(diff.ab.depth.1065)
#ASV112 - model failed - both
diff.ab.age.112 <- bbdml(formula = ASV112 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.112, color = "age", B=50)
summary(diff.ab.age.112)
diff.ab.depth.112 <- bbdml(formula = ASV112 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.112, color = "depth", B=50)
summary(diff.ab.depth.112)
#ASV1257 - model failed - both
diff.ab.age.1257 <- bbdml(formula = ASV1257 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.112, color = "age", B=50)
summary(diff.ab.age.112)
diff.ab.depth.1257 <- bbdml(formula = ASV1257 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.112, color = "depth", B=50)
summary(diff.ab.depth.112)
#ASV1687 - model failed - both
diff.ab.age.1687 <- bbdml(formula = ASV1687 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.112, color = "age", B=50)
summary(diff.ab.age.112)
diff.ab.depth.1687 <- bbdml(formula = ASV1687 ~ depth, phi.formula= ~ depth,  data=phylo.object)
plot(diff.ab.depth.112, color = "depth", B=50)
summary(diff.ab.depth.112)
#ASV2165 - subsurface significant
diff.ab.age.2165 <- bbdml(formula = ASV2165 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.2165, color = "age", B=50)
summary(diff.ab.age.2165)
diff.ab.depth.2165 <- bbdml(formula = ASV2165 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.2165, color = "depth", B=50)
summary(diff.ab.depth.2165)
#ASV2217 - model fails for both
diff.ab.age.2217 <- bbdml(formula = ASV2217 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.2217, color = "age", B=50)
summary(diff.ab.age.2217)
diff.ab.depth.2217 <- bbdml(formula = ASV2217 ~ depth, phi.formula= ~ depth,  data=phylo.object)
plot(diff.ab.depth.2217, color = "depth", B=50)
summary(diff.ab.depth.2217)
#ASV23 - model fails for both
diff.ab.age.23 <- bbdml(formula = ASV23 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.23, color = "age", B=50)
summary(diff.ab.age.23)
diff.ab.depth.23 <- bbdml(formula = ASV23 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.23, color = "depth", B=50)
summary(diff.ab.depth.23)
#ASV2538 - model fails for both
diff.ab.age.2538 <- bbdml(formula = ASV2538 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.2538, color = "age", B=50)
summary(diff.ab.age.2538)
diff.ab.depth.2538 <- bbdml(formula = ASV2538 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.2538, color = "depth", B=50)
summary(diff.ab.depth.2538)
#ASV2831 - not significant significant
diff.ab.age.2831 <- bbdml(formula = ASV2831 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.2831, color = "age", B=50)
summary(diff.ab.age.2831)
diff.ab.depth.2831 <- bbdml(formula = ASV2831 ~ depth, phi.formula= ~ depth,  data=phylo.object)
plot(diff.ab.depth.2831, color = "depth", B=50)
summary(diff.ab.depth.2831)
#ASV3586 - not significant
diff.ab.age.3586 <- bbdml(formula = ASV3586 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.3586, color = "age", B=50)
summary(diff.ab.age.3586)
diff.ab.depth.3586 <- bbdml(formula = ASV3586 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.3586, color = "depth", B=50)
summary(diff.ab.depth.3586)
#ASV5564 - model fails
diff.ab.age.5564 <- bbdml(formula = ASV5564 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.5564, color = "age", B=50)
summary(diff.ab.age.5564)
diff.ab.depth.5564 <- bbdml(formula = ASV5564 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.5564, color = "depth", B=50)
summary(diff.ab.depth.5564)
#ASV634 - not significant
diff.ab.age.634 <- bbdml(formula = ASV634 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.634, color = "age", B=50)
summary(diff.ab.age.634)
diff.ab.depth.634 <- bbdml(formula = ASV634 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.634, color = "depth", B=50)
summary(diff.ab.depth.634)
#ASV656 - not significant
diff.ab.age.656 <- bbdml(formula = ASV656 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.656, color = "age", B=50)
summary(diff.ab.age.656)
diff.ab.depth.656 <- bbdml(formula = ASV656 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.656, color = "depth", B=50)
summary(diff.ab.depth.656)
#ASV88 - model fails  both
diff.ab.age.88 <- bbdml(formula = ASV88 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.88, color = "age", B=50)
summary(diff.ab.age.88)
diff.ab.depth.88 <- bbdml(formula = ASV88 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.88, color = "depth", B=50)
summary(diff.ab.depth.88)



#Psychromos --------
#ASV12 - juvenile significant and depth subsurface
Psycro.genus <- bbdml(formula=ASV1711 ~ age, phi.formula = ~age, data=phylo.object.genus)
summary(Psycro.genus)
plot(Psycro.genus, color="age")
Psychro.depth <- bbdml(formula=ASV1711~depth, phi.formula=~depth, data=phylo.object.genus)
summary(Psychro.depth)
plot(Psychro.depth, color="depth")
#ASV1711 mature significant and surface significant
diff.ab.age.1711 <- bbdml(formula=ASV1711 ~ age, phi.formula=~age, data=phylo.object)
summary(diff.ab.age.1711)
plot(diff.ab.age.1711, color="age")
diff.ab.depth.1711 <- bbdml(formula=ASV1711~depth, phi.formula=~depth, data=phylo.object)
summary(diff.ab.depth.1711)
plot(diff.ab.depth.1711, color="depth")
#ASV12 - not significant 
diff.ab.age.12 <- bbdml(formula = ASV12 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.12, color = "age", B=50)
summary(diff.ab.age.12)
diff.ab.depth.12 <- bbdml(formula = ASV12 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.12, color = "depth", B=50)
summary(diff.ab.depth.12)
#ASV1425 - no significance
diff.ab.age.1425 <- bbdml(formula = ASV1425 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.1425, color = "age", B=50)
summary(diff.ab.age.1425)
diff.ab.depth.1425 <- bbdml(formula = ASV1425 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.1425, color = "depth", B=50)
summary(diff.ab.depth.1425)
#ASV1651 - no significance
diff.ab.age.1651 <- bbdml(formula = ASV1651 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.1651, color = "age", B=50)
summary(diff.ab.age.1651)
diff.ab.depth.1651 <- bbdml(formula = ASV1651 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.1651, color = "depth", B=50)
summary(diff.ab.depth.1651)
#ASV1652 - no significance
diff.ab.age.1652 <- bbdml(formula = ASV1652 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.1652, color = "age", B=50)
summary(diff.ab.age.1652)
diff.ab.depth.1652 <- bbdml(formula = ASV1652 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.1652, color = "depth", B=50)
summary(diff.ab.depth.1652)
#ASV1711 - significant mature and surface 
diff.ab.age.1711 <- bbdml(formula = ASV1711 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.1711, color = "age", B=50)
summary(diff.ab.age.1711)
diff.ab.depth.1711 <- bbdml(formula = ASV1711 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.1711, color = "depth", B=50)
summary(diff.ab.depth.1711)
#ASV2402 - no significance
diff.ab.age.2402 <- bbdml(formula = ASV2402 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.2402, color = "age", B=50)
summary(diff.ab.age.2402)
diff.ab.depth.2402 <- bbdml(formula = ASV2402 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.2402, color = "depth", B=50)
summary(diff.ab.depth.2402)
#ASV2525 - no significance
diff.ab.age.2525 <- bbdml(formula = ASV2525 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.2525, color = "age", B=50)
summary(diff.ab.age.2525)
diff.ab.depth.2525 <- bbdml(formula = ASV2525 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.2525, color = "depth", B=50)
summary(diff.ab.depth.2525)
#ASV296 - no significance
diff.ab.age.296 <- bbdml(formula = ASV296 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.296, color = "age", B=50)
summary(diff.ab.age.296)
diff.ab.depth.296 <- bbdml(formula = ASV296 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.296, color = "depth", B=50)
summary(diff.ab.depth.296)
#ASV311 - no significance
diff.ab.age.311 <- bbdml(formula = ASV311 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.311, color = "age", B=50)
summary(diff.ab.age.311)
diff.ab.depth.311 <- bbdml(formula = ASV311 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.311, color = "depth", B=50)
summary(diff.ab.depth.311)
#ASV3335 - no significance
diff.ab.age.3335 <- bbdml(formula = ASV3335 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.3335, color = "age", B=50)
summary(diff.ab.age.3335)
diff.ab.depth.3335 <- bbdml(formula = ASV296 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.3335, color = "depth", B=50)
summary(diff.ab.depth.3335)
#ASV3468 - no significance
diff.ab.age.3468 <- bbdml(formula = ASV3468 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.3468, color = "age", B=50)
summary(diff.ab.age.3468)
diff.ab.depth.3468 <- bbdml(formula = ASV3468 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.3468, color = "depth", B=50)
summary(diff.ab.depth.3468)
#ASV367 - no significance
diff.ab.age.367 <- bbdml(formula = ASV367 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.367, color = "age", B=50)
summary(diff.ab.age.367)
diff.ab.depth.367 <- bbdml(formula = ASV367 ~ depth, phi.formula= ~ depth,  data=phylo.object)
plot(diff.ab.depth.367, color = "depth", B=50)
summary(diff.ab.depth.367)
#ASV3704 - not significant
diff.ab.age.296 <- bbdml(formula = ASV296 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.296, color = "age", B=50)
summary(diff.ab.age.296)
diff.ab.depth.296 <- bbdml(formula = ASV296 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.296, color = "depth", B=50)
summary(diff.ab.depth.296)
#ASV3773 - not significant
diff.ab.age.3773 <- bbdml(formula = ASV3773 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.3773, color = "age", B=50)
summary(diff.ab.age.3773)
diff.ab.depth.3373 <- bbdml(formula = ASV3373 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.3373, color = "depth", B=50)
summary(diff.ab.depth.3373)
#ASV3785 - not significant
diff.ab.age.3785 <- bbdml(formula = ASV3785 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.3785, color = "age", B=50)
summary(diff.ab.age.3785)
diff.ab.depth.3785 <- bbdml(formula = ASV3785 ~ depth, phi.formula= ~ depth,  data=phylo.object)
plot(diff.ab.depth.3785, color = "depth", B=50)
summary(diff.ab.depth.3785)
#ASV4500 - not significant
diff.ab.age.4500 <- bbdml(formula = ASV4500 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.4500, color = "age", B=50)
summary(diff.ab.age.4500)
diff.ab.depth.4500 <- bbdml(formula = ASV4500 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.4500, color = "depth", B=50)
summary(diff.ab.depth.4500)
#ASV4636 - subsurface significant
diff.ab.age.4636 <- bbdml(formula = ASV4636 ~ age, phi.formula= ~ age,  data=phylo.object)
plot(diff.ab.age.4636, color = "age", B=50)
summary(diff.ab.age.4636)
diff.ab.depth.4636 <- bbdml(formula = ASV4636 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.4636, color = "depth", B=50)
summary(diff.ab.depth.4636)
#ASV4775 - no significant
diff.ab.age.4775 <- bbdml(formula = ASV4775 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.4775, color = "age", B=50)
summary(diff.ab.age.4775)
diff.ab.depth.4775 <- bbdml(formula = ASV4775 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.4775, color = "depth", B=50)
summary(diff.ab.depth.4775)
#ASV5716 - no significance
diff.ab.age.5716 <- bbdml(formula = ASV5716 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.5716, color = "age", B=50)
summary(diff.ab.age.5716)
diff.ab.depth.5716 <- bbdml(formula = ASV5716 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.5716, color = "depth", B=50)
summary(diff.ab.depth.5716)
#ASV6 - no significance
diff.ab.age.6 <- bbdml(formula = ASV6 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.6, color = "age", B=50)
summary(diff.ab.age.6)
diff.ab.depth.6 <- bbdml(formula = ASV6 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.6, color = "depth", B=50)
summary(diff.ab.depth.6)
#ASV6148 - significant mature and subsurface
diff.ab.age.6148 <- bbdml(formula = ASV6148 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.6148, color = "age", B=50)
summary(diff.ab.age.6148)
diff.ab.depth.6148 <- bbdml(formula = ASV6148 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.6148, color = "depth", B=50)
summary(diff.ab.depth.6148)
#ASV617 - no significance
diff.ab.age.617 <- bbdml(formula = ASV617 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.617, color = "age", B=50)
summary(diff.ab.age.617)
diff.ab.depth.617 <- bbdml(formula = ASV617 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.617, color = "depth", B=50)
summary(diff.ab.depth.617)
#ASV7040 - mature signifiance and subsurface significance
diff.ab.age.7040 <- bbdml(formula = ASV7040 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.7040, color = "age", B=50)
summary(diff.ab.age.7040)
diff.ab.depth.7040 <- bbdml(formula = ASV7040 ~ depth, phi.formula= ~ depth,  data=phylo.object)
plot(diff.ab.depth.7040, color = "depth", B=50)
summary(diff.ab.depth.7040)
#ASV758 - no significance
diff.ab.age.758 <- bbdml(formula = ASV758 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.758, color = "age", B=50)
summary(diff.ab.age.758)
diff.ab.depth.758 <- bbdml(formula = ASV758 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.758, color = "depth", B=50)
summary(diff.ab.depth.758)
#ASV801 - no significance
diff.ab.age.801 <- bbdml(formula = ASV801 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.801, color = "age", B=50)
summary(diff.ab.age.801)
diff.ab.depth.801 <- bbdml(formula = ASV801 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.801, color = "depth", B=50)
summary(diff.ab.depth.801)
#ASV970 - not significant
diff.ab.age.970 <- bbdml(formula = ASV970 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.970, color = "age", B=50)
summary(diff.ab.age.970)
diff.ab.depth.970 <- bbdml(formula = ASV970 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.970, color = "depth", B=50)
summary(diff.ab.depth.970)

#
#Blastopirellula -----------
Blastopirellula.age <- bbdml(formula=ASV7999 ~ age, phi.formula = ~age,  data=phylo.object.genus)
summary(Blastopirellula.age)
plot(Blastopirellula.age, color="age")
Blastopirellula.depth <- bbdml(formula=ASV7999 ~ age, phi.formula = ~depth,  data=phylo.object.genus)
summary(Blastopirellula.depth)
plot(Blastopirellula.depth, color="depth")
#ASV7999 - mature and subsurface
diff.ab.age.7999 <- bbdml(formula = ASV7999 ~ age, phi.formula= ~ age, data=phylo.object.genus)
plot(diff.ab.age.7999, color = "age", B=50)
summary(diff.ab.age.7999)
diff.ab.depth.1351 <- bbdml(formula = ASV1351 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.1351, color = "depth", B=50)
summary(diff.ab.depth.1351)
#ASV1351 - subsurface
diff.ab.age.1351 <- bbdml(formula = ASV1351 ~ age, phi.formula= ~ age, data=phylo.object)
summary(diff.ab.age.1351)
diff.ab.depth.1351 <- bbdml(formula = ASV1351 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.1351, color = "depth", B=50)
summary(diff.ab.depth.1351)
#ASV796 - no significance
diff.ab.age.796 <- bbdml(formula = ASV796 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.796, color = "age", B=50)
summary(diff.ab.age.796)
diff.ab.depth.796 <- bbdml(formula = ASV796 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.796, color = "depth", B=50)
summary(diff.ab.depth.796)
#ASV1077 - mature significant
diff.ab.age.1077 <- bbdml(formula = ASV1077 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.1077, color = "age", B=50)
summary(diff.ab.age.1077)
diff.ab.depth.1077 <- bbdml(formula = ASV12 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.1077, color = "depth", B=50)
summary(diff.ab.depth.1077)
#ASV1409 - no significance
diff.ab.age.1409 <- bbdml(formula = ASV1409 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.1409, color = "age", B=50)
summary(diff.ab.age.1409)
diff.ab.depth.1409 <- bbdml(formula = ASV12 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.1409, color = "depth", B=50)
summary(diff.ab.depth.1409)
#ASV762 - no significance
diff.ab.age.762 <- bbdml(formula = ASV762 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.762, color = "age", B=50)
summary(diff.ab.age.762)
diff.ab.depth.762 <- bbdml(formula = ASV762 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.762, color = "depth", B=50)
summary(diff.ab.depth.762)
#ASV959 - subsurface significant
diff.ab.age.959 <- bbdml(formula = ASV959 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.959, color = "age", B=50)
summary(diff.ab.age.959)
diff.ab.depth.959 <- bbdml(formula = ASV959 ~ depth, phi.formula= ~ depth,  data=phylo.object)
plot(diff.ab.depth.959, color = "depth", B=50)
summary(diff.ab.depth.959)
#ASV871 - no significance
diff.ab.age.871 <- bbdml(formula = ASV871 ~ age, phi.formula= ~ age,  data=phylo.object)
plot(diff.ab.age.871, color = "age", B=50)
summary(diff.ab.age.871)
diff.ab.depth.871 <- bbdml(formula = ASV871 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.871, color = "depth", B=50)
summary(diff.ab.depth.871)
#ASV996 - no significance
diff.ab.age.996 <- bbdml(formula = ASV996 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.996, color = "age", B=50)
summary(diff.ab.age.996)
diff.ab.depth.996 <- bbdml(formula = ASV996 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.996, color = "depth", B=50)
summary(diff.ab.depth.996)
#ASV707 - no significance
diff.ab.age.707 <- bbdml(formula = ASV707 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.707, color = "age", B=50)
summary(diff.ab.age.707)
diff.ab.depth.707 <- bbdml(formula = ASV707 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.707, color = "depth", B=50)
summary(diff.ab.depth.707)
#ASV338 - no significance
diff.ab.age.338 <- bbdml(formula = ASV338 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.338, color = "age", B=50)
summary(diff.ab.age.338)
diff.ab.depth.338 <- bbdml(formula = ASV338 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.338, color = "depth", B=50)
summary(diff.ab.depth.338)
#ASV754 - no significance
diff.ab.age.754 <- bbdml(formula = ASV754 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.754, color = "age", B=50)
summary(diff.ab.age.754)
diff.ab.depth.754 <- bbdml(formula = ASV754 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.754, color = "depth", B=50)
summary(diff.ab.depth.754)
#ASV679 - no significance
diff.ab.age.679 <- bbdml(formula = ASV679 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.679, color = "age", B=50)
summary(diff.ab.age.679)
diff.ab.depth.679 <- bbdml(formula = ASV679 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.679, color = "depth", B=50)
summary(diff.ab.depth.679)
#ASV703 - no significance
diff.ab.age.703 <- bbdml(formula = ASV703 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.703, color = "age", B=50)
summary(diff.ab.age.703)
diff.ab.depth.703 <- bbdml(formula = ASV703 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.703, color = "depth", B=50)
summary(diff.ab.depth.703)
#ASV185 - mature significant and surface significant
diff.ab.age.185 <- bbdml(formula = ASV185 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.185, color = "age", B=50)
summary(diff.ab.age.185)
diff.ab.depth.185 <- bbdml(formula = ASV185 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.185, color = "depth", B=50)
summary(diff.ab.depth.185)
#ASV427 - no significance
diff.ab.age.427 <- bbdml(formula = ASV427 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.427, color = "age", B=50)
summary(diff.ab.age.427)
diff.ab.depth.427 <- bbdml(formula = ASV427 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.427, color = "depth", B=50)
summary(diff.ab.depth.427)
#ASV332 - no significance
diff.ab.age.332 <- bbdml(formula = ASV332 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.332, color = "age", B=50)
summary(diff.ab.age.332)
diff.ab.depth.332 <- bbdml(formula = ASV332 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.332, color = "depth", B=50)
summary(diff.ab.depth.332)
#ASV377 - no significance
diff.ab.age.377 <- bbdml(formula = ASV337 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.337, color = "age", B=50)
summary(diff.ab.age.337)
diff.ab.depth.337 <- bbdml(formula = ASV337 ~ depth, phi.formula= ~ depth,  data=phylo.object)
plot(diff.ab.depth.337, color = "depth", B=50)
summary(diff.ab.depth.337)
#ASV225 - not significant
diff.ab.age.225 <- bbdml(formula = ASV225 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.225, color = "age", B=50)
summary(diff.ab.age.225)
diff.ab.depth.225 <- bbdml(formula = ASV225 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.225, color = "depth", B=50)
summary(diff.ab.depth.225)
#ASV99 - not significant
diff.ab.age.99 <- bbdml(formula = ASV99 ~ age, phi.formula= ~ age,  data=phylo.object)
plot(diff.ab.age.99, color = "age", B=50)
summary(diff.ab.age.99)
diff.ab.depth.99 <- bbdml(formula = ASV99 ~ depth, phi.formula= ~ depth,  data=phylo.object)
plot(diff.ab.depth.99, color = "depth", B=50)
summary(diff.ab.depth.99)
#ASV34 - not significant
diff.ab.age.34 <- bbdml(formula = ASV34 ~ age, phi.formula= ~ age,  data=phylo.object)
plot(diff.ab.age.34, color = "age", B=50)
summary(diff.ab.age.34)
diff.ab.depth.34 <- bbdml(formula = ASV34 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.34, color = "depth", B=50)
summary(diff.ab.depth.34)
#ASV244 - not significant
diff.ab.age.244 <- bbdml(formula = ASV244 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.244, color = "age", B=50)
summary(diff.ab.age.244)
diff.ab.depth.244 <- bbdml(formula = ASV244 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.244, color = "depth", B=50)
summary(diff.ab.depth.244)
#ASV20 - not significant
diff.ab.age.20 <- bbdml(formula = ASV20 ~ age, phi.formula= ~ age,  data=phylo.object)
plot(diff.ab.age.20, color = "age", B=50)
summary(diff.ab.age.20)
diff.ab.depth.20 <- bbdml(formula = ASV20 ~ depth, phi.formula= ~ depth,  data=phylo.object)
plot(diff.ab.depth.20, color = "depth", B=50)
summary(diff.ab.depth.20)
#ASV83 - no significance
diff.ab.age.83 <- bbdml(formula = ASV83 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.83, color = "age", B=50)
summary(diff.ab.age.83)
diff.ab.depth.83 <- bbdml(formula = ASV83 ~ depth, phi.formula= ~ depth,  data=phylo.object)
plot(diff.ab.depth.83, color = "depth", B=50)
summary(diff.ab.depth.83)
#ASV52 - no significance
diff.ab.age.52 <- bbdml(formula = ASV52 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.52, color = "age", B=50)
summary(diff.ab.age.52)
diff.ab.depth.52 <- bbdml(formula = ASV52 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.52, color = "depth", B=50)
summary(diff.ab.depth.52)
#ASV8 - no significance
diff.ab.age.8 <- bbdml(formula = ASV8 ~ age, phi.formula= ~ age,  data=phylo.object)
plot(diff.ab.age.8, color = "age", B=50)
summary(diff.ab.age.8)
diff.ab.depth.8 <- bbdml(formula = ASV8 ~ depth, phi.formula= ~ depth,  data=phylo.object)
plot(diff.ab.depth.8, color = "depth", B=50)
summary(diff.ab.depth.8)
#ASV11 - subsurface significant
diff.ab.age.11 <- bbdml(formula = ASV11 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.11, color = "age", B=50)
summary(diff.ab.age.11)
diff.ab.depth.11 <- bbdml(formula = ASV11 ~ depth, phi.formula= ~ depth,  data=phylo.object)
plot(diff.ab.depth.11, color = "depth", B=50)
summary(diff.ab.depth.11)
#ASV10 - no significance
diff.ab.age.10 <- bbdml(formula = ASV10 ~ age, phi.formula= ~ age,   data=phylo.object)
plot(diff.ab.age.10, color = "age", B=50)
summary(diff.ab.age.10)
diff.ab.depth.10 <- bbdml(formula = ASV10 ~ depth, phi.formula= ~ depth,  data=phylo.object)
plot(diff.ab.depth.10, color = "depth", B=50)
summary(diff.ab.depth.10)
#OM43 --------
OM43.age <- bbdml(formula=ASV294~age, phi.formula=~age, data=phylo.object.genus)
summary(OM43.age)
plot(OM43.age, color="age")
OM43.depth <- bbdml(formula=ASV294~depth, phi.formula=~depth, data=phylo.object.genus)
summary(OM43.depth)
plot(OM43.age, color="depth")
#ASV294 Mature significant
diff.ab.age.294 <- bbdml(formula = ASV294 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.294, color = "age", B=50)
summary(diff.ab.age.294)
#ASV1404 - subsurface significant 
diff.ab.age.1404 <- bbdml(formula = ASV1404 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.1404, color = "age", B=50)
summary(diff.ab.age.1404)
diff.ab.depth.1404 <- bbdml(formula = ASV1404 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.1404, color = "depth", B=50)
summary(diff.ab.depth.1404)
#Lutimos #Genus level mature and subsurface------
Lutimos.age <- bbdml(formula=ASV11861~age, phi.formula=~age, data=phylo.object.genus)
summary(Lutimos.age)
plot(Lutimos.age, color="age")
Lutimos.depth <- bbdml(formula=ASV11861~age, phi.formula=~age, data=phylo.object.genus)
summary(Lutimos.depth)
plot(Lutimos.depth, color="depth")
#ASV11861 mature significant and subsurface significant - 
ASV.11861.age <- bbdml(formula=ASV11861~age, phi.formula=~age, data=phylo.object)
summary(ASV.11861.age)
plot(ASV.11861.age, color="age")
ASV11861.depth <- bbdml(formula=ASV11861~depth, phi.formula=~depth, data=phylo.object)
summary(ASV11861.depth)
plot(ASV11861.depth, color="depth")
#ASV116 - no significance
diff.ab.age.116 <- bbdml(formula = ASV116 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.116, color = "age", B=50)
summary(diff.ab.age.116)
diff.ab.depth.116 <- bbdml(formula = ASV116 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.116, color = "depth", B=50)
summary(diff.ab.depth.116)
#ASV1619 - no significance
diff.ab.age.1619 <- bbdml(formula = ASV1619 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.1619, color = "age", B=50)
summary(diff.ab.age.1619)
diff.ab.depth.1619 <- bbdml(formula = ASV1619 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.1619, color = "depth", B=50)
summary(diff.ab.depth.1619)
#ASV2501 - no significance
diff.ab.age.2501 <- bbdml(formula = ASV2501 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.2501, color = "age", B=50)
summary(diff.ab.age.2501)
diff.ab.depth.2501 <- bbdml(formula = ASV2501 ~ depth, phi.formula= ~ depth,  data=phylo.object)
plot(diff.ab.depth.2501, color = "depth", B=50)
summary(diff.ab.depth.2501)
#ASV3136 - no significance
diff.ab.age.3136 <- bbdml(formula = ASV3136 ~ age, phi.formula= ~ age,  data=phylo.object)
plot(diff.ab.age.3136, color = "age", B=50)
summary(diff.ab.age.3136)
diff.ab.depth.3136 <- bbdml(formula = ASV3136 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.3136, color = "depth", B=50)
summary(diff.ab.depth.3136)
#ASV4077 - no significance
diff.ab.age.4077 <- bbdml(formula = ASV4077 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.4077, color = "age", B=50)
summary(diff.ab.age.4077)
diff.ab.depth.4077 <- bbdml(formula = ASV4077 ~ depth, phi.formula= ~ depth,  data=phylo.object)
plot(diff.ab.depth.4077, color = "depth", B=50)
summary(diff.ab.depth.4077)
#ASV54 - no significance
diff.ab.age.54 <- bbdml(formula = ASV54 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.54, color = "age", B=50)
summary(diff.ab.age.54)
diff.ab.depth.54 <- bbdml(formula = ASV54 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.54, color = "depth", B=50)
summary(diff.ab.depth.54)
#ASV6102 - no significance
diff.ab.age.6102 <- bbdml(formula = ASV6102 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.6102, color = "age", B=50)
summary(diff.ab.age.6102)
diff.ab.depth.6102 <- bbdml(formula = ASV6102 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.6102, color = "depth", B=50)
summary(diff.ab.depth.6102)
#ASV663 - no significance
diff.ab.age.663 <- bbdml(formula = ASV663 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.663, color = "age", B=50)
summary(diff.ab.age.663)
diff.ab.depth.663 <- bbdml(formula = ASV663 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.663, color = "depth", B=50)
summary(diff.ab.depth.663)
#Rubripirellula (genus level mature significant and subsurface--------
Rubri.10363.age <- bbdml(formula=ASV10363~age, phi.formula=~age, data=phylo.object.genus)
summary(Rubri.10363.age)
plot(Rubri.10363.age, color="age")
Rubri.10363.depth <- bbdml(formula=ASV10363~age, phi.formula=~age, data=phylo.object.genus)
summary(Rubri.10363.depth)
plot(Rubri.10363.depth, color="depth")
#ASV10363
diff.ab.age.10363 <- bbdml(formula=ASV10363~age, phi.formula=~age, data=phylo.object)
summary(diff.ab.age.10363)
diff.ab.depth.10363 <- bbdml(formula=ASV10363~depth, phi.formula=~depth, data=phylo.object)
summary(diff.ab.depth.10363)
#ASV2459 - no significance
diff.ab.age.2459 <- bbdml(formula = ASV2459 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.2459, color = "age", B=50)
summary(diff.ab.age.2459)
diff.ab.depth.2459 <- bbdml(formula = ASV2459 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.2459, color = "depth", B=50)
summary(diff.ab.depth.2459)
#ASV1272 - no significance
diff.ab.age.1272 <- bbdml(formula = ASV1272 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.1272, color = "age", B=50)
summary(diff.ab.age.1272)
diff.ab.depth.1272 <- bbdml(formula = ASV1272 ~ depth, phi.formula= ~ depth,  data=phylo.object)
plot(diff.ab.depth.1272, color = "depth", B=50)
summary(diff.ab.depth.1272)
#ASV749 - no significance
diff.ab.age.749 <- bbdml(formula = ASV749 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.749, color = "age", B=50)
summary(diff.ab.age.749)
diff.ab.depth.749 <- bbdml(formula = ASV749 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.749, color = "depth", B=50)
summary(diff.ab.depth.749)
#ASV817 - no significance
diff.ab.age.817 <- bbdml(formula = ASV817 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.817, color = "age", B=50)
summary(diff.ab.age.817)
diff.ab.depth.817 <- bbdml(formula = ASV817 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.817, color = "depth", B=50)
summary(diff.ab.depth.817)
#ASV793 - no significance
diff.ab.age.793 <- bbdml(formula = ASV793 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.793, color = "age", B=50)
summary(diff.ab.age.793)
diff.ab.depth.793 <- bbdml(formula = ASV793 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.793, color = "depth", B=50)
summary(diff.ab.depth.793)
#ASV109 - no significance
diff.ab.age.109 <- bbdml(formula = ASV109 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.109, color = "age", B=50)
summary(diff.ab.age.109)
diff.ab.depth.109 <- bbdml(formula = ASV109 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.109, color = "depth", B=50)
summary(diff.ab.depth.109)
#ASV120 - no significance
diff.ab.age.120 <- bbdml(formula = ASV120 ~ age, phi.formula= ~ age,  data=phylo.object)
plot(diff.ab.age.120, color = "age", B=50)
summary(diff.ab.age.120)
diff.ab.depth.120 <- bbdml(formula = ASV120 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.120, color = "depth", B=50)
summary(diff.ab.depth.120)
#ASV28 - no significance
diff.ab.age.28 <- bbdml(formula = ASV28 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.28, color = "age", B=50)
summary(diff.ab.age.28)
diff.ab.depth.28 <- bbdml(formula = ASV28 ~ depth, phi.formula= ~ depth, data=phylo.object)
plot(diff.ab.depth.28, color = "depth", B=50)
summary(diff.ab.depth.28)
#Roseibacillus Mature significant & subsurface -------
Rose.age <- bbdml(formula=ASV8113~age, phi.formula=~age, data=phylo.object.genus)
summary(Rose.age)
plot(Rose.age, color="age")
Rose.depth <- bbdml(formula=ASV8113~depth, phi.formula=~depth, data=phylo.object.genus)
summary(Rose.depth)
plot(Rose.age, color="depth")
#ASV8113 only subsurface significant
diff.ab.age.8113 <- bbdml(formula=ASV8113~age, phi.formula=~age, data=phylo.object)
diff.ab.depth.8113 <- bbdml(formula=ASV8113~depth, phi.formula=~depth, data=phylo.object)
summary(diff.ab.depth.8113)
plot(diff.ab.depth.8113, color="depth")
#ASV1812 - no significance
diff.ab.age.1812 <- bbdml(formula = ASV1812 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.1812, color = "age", B=50)
diff.ab.depth.1812 <- bbdml(formula = ASV1812 ~ depth, phi.formula= ~ depth, data=phylo.object)
#ASV1497 - no significance
diff.ab.age.1497 <- bbdml(formula = ASV1497 ~ age, phi.formula= ~ age,data=phylo.object)
plot(diff.ab.age.1497, color = "age", B=50)
summary(diff.ab.age.1497)
diff.ab.depth.1497 <- bbdml(formula = ASV1497 ~ depth, phi.formula= ~ depth,data=phylo.object)
#ASV1228 - no significance
diff.ab.age.1228 <- bbdml(formula = ASV1228 ~ age, phi.formula= ~ age,data=phylo.object)
plot(diff.ab.age.1228, color = "age", B=50)
summary(diff.ab.age.1228)
diff.ab.depth.1228 <- bbdml(formula = ASV1228 ~ depth, phi.formula= ~ depth,data=phylo.object)
#ASV847 - no significance
diff.ab.age.847 <- bbdml(formula = ASV847 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.847, color = "age", B=50)
summary(diff.ab.age.847)
diff.ab.depth.847 <- bbdml(formula = ASV847 ~ depth, phi.formula= ~ depth, data=phylo.object)
#ASV27 - no significance
diff.ab.age.27 <- bbdml(formula = ASV27 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.27, color = "age", B=50)
summary(diff.ab.age.27)
diff.ab.depth.27 <- bbdml(formula = ASV27 ~ depth, phi.formula= ~ depth, data=phylo.object)
#ASV932 - no significance
diff.ab.age.932 <- bbdml(formula = ASV932 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.932, color = "age", B=50)
summary(diff.ab.age.932)
diff.ab.depth.932 <- bbdml(formula = ASV932 ~ depth, phi.formula= ~ depth, data=phylo.object)
#ASV406 - no significance
diff.ab.age.406 <- bbdml(formula = ASV406 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.406, color = "age", B=50)
summary(diff.ab.age.406)
diff.ab.depth.406 <- bbdml(formula = ASV406 ~ depth, phi.formula= ~ depth, data=phylo.object)
#ASV585 - no significance
diff.ab.age.585 <- bbdml(formula = ASV585 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.585, color = "age", B=50)
summary(diff.ab.age.585)
diff.ab.depth.585 <- bbdml(formula = ASV585 ~ depth, phi.formula= ~ depth, data=phylo.object)
#ASV651 - no significance
diff.ab.age.651 <- bbdml(formula = ASV651 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.651, color = "age", B=50)
summary(diff.ab.age.651)
diff.ab.depth.651 <- bbdml(formula = ASV651 ~ depth, phi.formula= ~ depth, data=phylo.object)
#ASV486 - no significance
diff.ab.age.486 <- bbdml(formula = ASV486 ~ age, phi.formula= ~ age,  data=phylo.object)
plot(diff.ab.age.486, color = "age", B=50)
summary(diff.ab.age.486)
diff.ab.depth.486 <- bbdml(formula = ASV486 ~ depth, phi.formula= ~ depth,  data=phylo.object)
#ASV298 - no significance
diff.ab.age.298 <- bbdml(formula = ASV298 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.298, color = "age", B=50)
summary(diff.ab.age.298)
diff.ab.depth.298 <- bbdml(formula = ASV298 ~ depth, phi.formula= ~ depth, data=phylo.object)
#ASV297 - subsurface significant
diff.ab.age.297 <- bbdml(formula = ASV297 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.297, color = "age", B=50)
summary(diff.ab.age.297)
diff.ab.depth.297 <- bbdml(formula = ASV297 ~ depth, phi.formula= ~ depth, data=phylo.object)
summary(diff.ab.depth.297)
plot(diff.ab.depth.297, color="depth")
#ASV186 - no significance
diff.ab.age.186 <- bbdml(formula = ASV186 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.186, color = "age", B=50)
summary(diff.ab.age.186)
diff.ab.depth.186 <- bbdml(formula = ASV186 ~ depth, phi.formula= ~ depth, data=phylo.object)
#ASV64 - no significance
diff.ab.age.64 <- bbdml(formula = ASV64 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.64, color = "age", B=50)
summary(diff.ab.age.64)
diff.ab.depth.64 <- bbdml(formula = ASV64 ~ depth, phi.formula= ~ depth, data=phylo.object)
summary(diff.ab.depth.64)
#ASV184 - no significance
diff.ab.age.184 <- bbdml(formula = ASV184 ~ age, phi.formula= ~ age,  data=phylo.object)
plot(diff.ab.age.184, color = "age", B=50)
summary(diff.ab.age.184)
diff.ab.depth.184 <- bbdml(formula = ASV184 ~ depth, phi.formula= ~ depth,  data=phylo.object)
summary(diff.ab.depth.184)
#ASV122 - no significance
diff.ab.age.122 <- bbdml(formula = ASV122 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.122, color = "age", B=50)
summary(diff.ab.age.122)
diff.ab.depth.122 <- bbdml(formula = ASV122 ~ depth, phi.formula= ~ depth, data=phylo.object)
summary(diff.ab.depth.122)
#ASV90 - no significance
diff.ab.age.90 <- bbdml(formula = ASV90 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.90, color = "age", B=50)
summary(diff.ab.age.90)
diff.ab.depth.90 <- bbdml(formula = ASV90 ~ depth, phi.formula= ~ depth, data=phylo.object)
#ASV37 - no significance
diff.ab.age.37 <- bbdml(formula = ASV37 ~ age, phi.formula= ~ age,  data=phylo.object)
plot(diff.ab.age.37, color = "age", B=50)
summary(diff.ab.age.37)
diff.ab.depth.37 <- bbdml(formula = ASV37 ~ depth, phi.formula= ~ depth,  data=phylo.object)
#ASV73 - no significance 
diff.ab.age.73 <- bbdml(formula = ASV73 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.73, color = "age", B=50)
summary(diff.ab.age.73)
diff.ab.depth.73 <- bbdml(formula = ASV73 ~ depth, phi.formula= ~ depth, data=phylo.object)
#ASV67 - mature significant
diff.ab.age.67 <- bbdml(formula = ASV67 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.67, color = "age", B=50)
summary(diff.ab.age.67)
diff.ab.depth.67 <- bbdml(formula = ASV67 ~ depth, phi.formula= ~ depth, data=phylo.object)
summary(diff.ab.depth.67)
#ASV43 - no significance
diff.ab.age.43 <- bbdml(formula = ASV43 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.43, color = "age", B=50)
summary(diff.ab.age.43)
diff.ab.depth.43 <- bbdml(formula = ASV43 ~ depth, phi.formula= ~ depth, data=phylo.object)
#ASV19 - no significance
diff.ab.age.19 <- bbdml(formula = ASV19 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.19, color = "age", B=50)
summary(diff.ab.age.19)
diff.ab.depth.19 <- bbdml(formula = ASV19 ~ depth, phi.formula= ~ depth, data=phylo.object)
#Lentimos Juvenile significant-------
Lentimos.age <- bbdml(formula=ASV10160~age, phi.formula=~age, data=phylo.object.genus)
summary(Lentimos.age)
plot(Lentimos.age, color="age")
#ASV10160 - juvenile significant
diff.ab.age.10160 <- bbdml(formula=ASV10160~age, phi.formula=~age, data=phylo.object)
summary(diff.ab.age.10160)
plot(diff.ab.age.10160, color="age")
#ASV2184 - no significance
diff.ab.age.2184 <- bbdml(formula = ASV2184 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.2184, color = "age", B=50)
summary(diff.ab.age.2184)
#ASV769 - no significance
diff.ab.age.769 <- bbdml(formula = ASV769 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.769, color = "age", B=50)
summary(diff.ab.age.769)
#ASV781 - no significance
diff.ab.age.781 <- bbdml(formula = ASV781 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.781, color = "age", B=50)
summary(diff.ab.age.781)
#ASV154 - no significance
diff.ab.age.154 <- bbdml(formula = ASV154 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.154, color = "age", B=50)
summary(diff.ab.age.154)
#ASV157 - no significance
diff.ab.age.157 <- bbdml(formula = ASV157 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.157, color = "age", B=50)
summary(diff.ab.age.157)
#ASV145 - no significance
diff.ab.age.145 <- bbdml(formula = ASV145 ~ age, phi.formula= ~ age,  data=phylo.object)
plot(diff.ab.age.145, color = "age", B=50)
summary(diff.ab.age.145)
#ASV115 - no significance
diff.ab.age.115 <- bbdml(formula = ASV115 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.115, color = "age", B=50)
summary(diff.ab.age.115)
#Flavicella -----------
Flavi.age <- bbdml(formula=ASV8421~age, phi.formula=~age, data=phylo.object.genus)
summary(Flavi.age)
plot(Flavi.age, color="age")
Flavi.depth <- bbdml(formula=ASV8421~depth, phi.formula=~depth, data=phylo.object.genus)
summary(Flavi.depth)
plot(Flavi.age, color="depth")
#ASV8421 - mature significant
diff.ab.age.8421 <- bbdml(formula=ASV8421~age, phi.formula=~age, data=phylo.object)
summary(diff.ab.age.8421)
plot(diff.ab.age.8421, color="age")
#ASV2208 - no significance
diff.ab.age.2208 <- bbdml(formula = ASV2208 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.2208, color = "age", B=50)
summary(diff.ab.age.2208)
#ASV2852 - no significance
diff.ab.age.2852 <- bbdml(formula = ASV2852 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.2852, color = "age", B=50)
summary(diff.ab.age.2852)
#ASV1267 - no significance
diff.ab.age.1267 <- bbdml(formula = ASV1267 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.1267, color = "age", B=50)
summary(diff.ab.age.1267)
#ASV863 - no significance
diff.ab.age.863 <- bbdml(formula = ASV863 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.863, color = "age", B=50)
summary(diff.ab.age.863)
#ASV161 - no significance 
diff.ab.age.161 <- bbdml(formula = ASV161 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.161, color = "age", B=50)
summary(diff.ab.age.161)
#ASV42 - no significance
diff.ab.age.42 <- bbdml(formula = ASV42 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.42, color = "age", B=50)
summary(diff.ab.age.42)
#Alcanivorax Mature significant-------
Alcani.age <- bbdml(formula=ASV689~age, phi.formula=~age, data=phylo.object.genus)
summary(Alcani.age)
plot(Alcani.age, color="age")
#ASV689 - Mature significant
diff.ab.age.689 <- bbdml(formula=ASV689~age, phi.formula=~age, data=phylo.object)
summary(diff.ab.age.689)
plot(diff.ab.age.689, color="age")
#ASV1286 - no significance
diff.ab.age.1286 <- bbdml(formula = ASV1286 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.1286, color = "age", B=50)
summary(diff.ab.age.1286)
#ASV1596 - no significance
diff.ab.age.1596 <- bbdml(formula = ASV1596 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.1596, color = "age", B=50)
summary(diff.ab.age.1596)
#ASV2952 - no significance
diff.ab.age.2952 <- bbdml(formula = ASV2952 ~ age, phi.formula= ~ age,  data=phylo.object)
plot(diff.ab.age.2952, color = "age", B=50)
summary(diff.ab.age.2952)
#ASV4548 - no significance
diff.ab.age.4548 <- bbdml(formula = ASV4548 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.4548, color = "age", B=50)
summary(diff.ab.age.4548)
#ASV670 - no significance
diff.ab.age.670<- bbdml(formula = ASV670 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.670, color = "age", B=50)
summary(diff.ab.age.670)
#ASV689 - mature significant
diff.ab.age.689 <- bbdml(formula = ASV689 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.689, color = "age", B=50)
summary(diff.ab.age.689)
#Stenotrophomos -------
Stenot.age <- bbdml(formula=ASV5269~age, phi.formula=~age, data=phylo.object.genus)
summary(Stenot.age)
plot(Stenot.age, color="age")
#ASV5269 - 
diff.ab.age.5269 <- bbdml(formula=ASV5269~age, phi.formula = ~age, data=phylo.object)
summary(diff.ab.age.5269)
plot(diff.ab.age.5269, color="depth")
#ASV222 - no significance
diff.ab.age.222<- bbdml(formula = ASV222 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.222, color = "age", B=50)
summary(diff.ab.age.222)
#ASV236 - no significance
diff.ab.age.236 <- bbdml(formula = ASV236 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.236, color = "age", B=50)
summary(diff.ab.age.236)
#ASV237 - no significance
diff.ab.age.237 <- bbdml(formula = ASV237 ~ age, phi.formula= ~ age,  data=phylo.object)
plot(diff.ab.age.237, color = "age", B=50)
summary(diff.ab.age.237)
#ASV2570 - no significance
diff.ab.age.2570 <- bbdml(formula = ASV2570 ~ age, phi.formula= ~ age, data=phylo.object)
plot(diff.ab.age.2570, color = "age", B=50)
summary(diff.ab.age.2570)
#Halomos ------
#ASV13 - no significance
Halo.depth <- bbdml(formula=ASV2816~depth, phi.formula=~depth, data=phylo.object.genus)
summary(Halo.depth)
plot(Halo.depth, color="depth")
diff.ab.depth.13 <- bbdml(formula = ASV13 ~ depth, phi.formula= ~ depth,
                          data=phylo.object)
plot(diff.ab.depth.13, color = "depth", B=50)
summary(diff.ab.depth.13)

#ASV2816 - subsurface significant
diff.ab.depth.2816 <- bbdml(formula = ASV2816 ~ depth, phi.formula= ~ depth,
                          data=phylo.object)
plot(diff.ab.depth.2816, color = "depth", B=50)
summary(diff.ab.depth.2816)

#ASV4130 - no significance
diff.ab.depth.4130 <- bbdml(formula = ASV4130 ~ depth, phi.formula= ~ depth,
                          data=phylo.object)
plot(diff.ab.depth.4130, color = "depth", B=50)
summary(diff.ab.depth.4130)


#Flavicella --------
#ASV6543 - no significance
diff.ab.depth.6543 <- bbdml(formula = ASV6543 ~ depth, phi.formula= ~ depth,
                          data=phylo.object)
plot(diff.ab.depth.6543, color = "depth", B=50)
summary(diff.ab.depth.6543)

#ASV2577 - no significance
diff.ab.depth.2577 <- bbdml(formula = ASV2577 ~ depth, phi.formula= ~ depth,
                          data=phylo.object)
plot(diff.ab.depth.2577, color = "depth", B=50)
summary(diff.ab.depth.2577)

#ASV2208 - no significance
diff.ab.depth.2208 <- bbdml(formula = ASV2208 ~ depth, phi.formula= ~ depth,
                          data=phylo.object)
plot(diff.ab.depth.2208, color = "depth", B=50)
summary(diff.ab.depth.2208)

#ASV2852 - no significance
diff.ab.depth.2852 <- bbdml(formula = ASV2852 ~ depth, phi.formula= ~ depth,
                          data=phylo.object)
plot(diff.ab.depth.2852, color = "depth", B=50)
summary(diff.ab.depth.2852)

#ASV1267 - subsurface significant
diff.ab.depth.1267 <- bbdml(formula = ASV1267 ~ depth, phi.formula= ~ depth,
                          data=phylo.object)
plot(diff.ab.depth.1267, color = "depth", B=50)
summary(diff.ab.depth.1267)

#ASV863 - subsurface significant
diff.ab.depth.863 <- bbdml(formula = ASV863 ~ depth, phi.formula= ~ depth,
                          data=phylo.object)
plot(diff.ab.depth.863, color = "depth", B=50)
summary(diff.ab.depth.863)

#ASV161 - no significance
diff.ab.depth.161 <- bbdml(formula = ASV161 ~ depth, phi.formula= ~ depth,
                          data=phylo.object)
plot(diff.ab.depth.161, color = "depth", B=50)
summary(diff.ab.depth.161)

#ASV42 - no significance
diff.ab.depth.42 <- bbdml(formula = ASV42 ~ depth, phi.formula= ~ depth,
                          data=phylo.object)
plot(diff.ab.depth.42, color = "depth", B=50)
summary(diff.ab.depth.42)

##Now I'm going to try and calculate log fold change using DeSeq2 package ---------
library(remotes)
library(DESeq2)
library(ggplot2)

DE_ra <- read.csv("DESeq2_relativeabundance.csv")
View(DE_ra)
DE_ra <- DE_ra %>% tibble::column_to_rownames("Genera")
DE_meta <- read.csv("DESeq2_metadata.csv")
View(DE_meta)

dds.depth <- DESeqDataSetFromMatrix(countData=round(DE_ra), 
                              colData=DE_meta, 
                              design=~depth)
dds.depth
dds.depth.DESeq <- DESeq(dds.depth)
res <- results(dds.depth.DESeq)
View(results(dds.depth.DESeq, tidy=TRUE)) 

dds.age <- DESeqDataSetFromMatrix(countData=round(DE_ra), 
                                    colData=DE_meta, 
                                    design=~age)
dds.age
dds.age.DESeq <- DESeq(dds.age)
res.age <- results(dds.age.DESeq)
View(results(dds.age.DESeq, tidy=TRUE)) 

View(res.age)

##ggplot DESeq2 outputs
depth.log <- read.csv("DESeq2_depth_log2foldchange.csv")
View(depth.log)
depth.log.sub <- subset(depth.log, Genera =="Halomos" |
                          Genera =="Pseudoalteromos" |
                          Genera =="Psychromos" |
                          Genera =="Blastopirellula" |
                          Genera =="Flavicella" |
                          Genera == "Rubripirellula" |
                          Genera =="OM43_clade" |
                          Genera =="Lutimos" |
                          Genera=="Lewinella")
View(depth.log.sub)
ggplot(depth.log.sub, aes(x=log2FoldChange, y=factor(Genera, level=c('Psychromos',
                                                                      'Pseudoalteromos',
                                                                      'OM43_clade',
                                                                      'Halomos',
                                                                      'Rubripirellula',
                                                                      'Blastopirellula',
                                                                      'Lutimos',
                                                                      'Flavicella')),fill=DepthAssociation)) +
  geom_bar(stat="identity") +
  theme_bw() + 
  xlim(-1.7, 3) + 
  ylab("Genera") + xlab("Log2 Fold Change") +
  scale_fill_manual(values=c("forestgreen", "purple2"))

age.log <- read.csv("DESeq2_age_log2foldchange.csv")
View(age.log)
age.log.sub <- subset(age.log, Genera =="Synechococcus_CC9902" |
                          Genera =="Pseudoalteromos" |
                          Genera =="Psychromos" |
                          Genera =="Blastopirellula" |
                          Genera=="OM43_clade" |
                        Genera=="Rubripirellula" |
                        Genera=="Lutimos" |
                        Genera =="Roseibacillus" |
                        Genera =="Lentimos" |
                        Genera=="Flavicella" |
                        Genera=="Alcanivorax" |
                        Genera=="Stenotrophomos")
View(age.log.sub)
ggplot(age.log.sub, aes(x=log2FoldChange, y=Genera)) +
  geom_bar(stat="identity") +
  theme_bw() + xlim(-1.7, 3)

View(age.log.sub)
ggplot(age.log.sub, aes(x=log2FoldChange, y=factor(Genera, level=c('Roseibacillus',
                                                                   'Lentimos',
                                                                   'Stenotrophomos',
                                                                   'Psychromos',
                                                                   'Pseudoalteromos',
                                                                   'OM43_clade',
                                                                   'Alcanivorax',
                                                                   'Rubripirellula',
                                                                   'Blastopirellula',
                                                                   'Synechococcus_CC9902',
                                                                   'Lutimos',
                                                                   'Flavicella')), fill=AgeAssociation)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  ylab("Genera") + xlab("Log2 Fold Change") +
  scale_fill_manual(values=c("yellow3", "darkblue"))
  
       
#Biplot for ASVs --------------
ASV.coords.biplot <- read.csv("ASVcoordsforbiplot.csv")
UI.coords.biplot <- read.csv("UniqueID NMDS coordinates for biplot.csv")
View(UI.coords.biplot)
View(ASV.coords.biplot)
ASV.coords.biplot.depth <- subset(ASV.coords.biplot, Depth =="Surface" |
                                    Depth == "Subsurface")
View(ASV.coords.biplot.depth)
ASV.coords.biplot.age <- subset(ASV.coords.biplot, Age=="Mature" |
                                  Age=="Juvenile")
View(ASV.coords.biplot.age)

ggplot(UI.coords.biplot) + 
  geom_point(aes(x=NMDS1, y=NMDS2, color=Age, shape=Age), alpha=0.33, size=6) +
  scale_color_manual(values=c("yellow3", "darkblue")) +
  geom_segment(data=ASV.coords.biplot,
               aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               type="closed", color="gray", alpha=0.7) +
  geom_text(data = ASV.coords.biplot, aes(x = NMDS1, y = NMDS2, label = ASV),
            size = 3) +
  theme_bw() +
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 20),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  xlim(-1.3, 1.4)

ggplot(UI.coords.biplot) + 
  geom_point(aes(x=NMDS1, y=NMDS2, color=Depth, shape=Age), alpha=0.3, size=5) +
  geom_segment(data=ASV.coords.biplot.depth,
               aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               type="closed", color="gray", alpha=0.5) +
  geom_text(data = ASV.coords.biplot.depth, aes(x = NMDS1, y = NMDS2, label = ASV),
            size = 2.5) +
  theme_bw() +
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 20)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.7) + ylim(-1.5, 1.5)
#
##Corncob plotting data-----
corn.ready <- read.csv("Corncob calculated values for plotting genus level.csv")
View(corn.ready)
corn.ready.age <- subset(corn.ready, AnalysisType=="Age")
View(corn.ready.age)
corn.ready.depth <- subset(corn.ready, AnalysisType=="Depth")
View(corn.ready.depth)

ggplot(corn.ready.age, aes(x=average, y=factor(genus, level=c('Roseibacillus',
                                                              'Lentimos',
                                                              'Rubripirellula',
                                                              'Blastopirellula',
                                                              'Stenotrophomos',
                                                              'Alcanivorax',
                                                              'OM43_clade',
                                                              'Psychromos',
                                                              'Pseudoalteromos',
                                                              'Synechococcus_CC9902',
                                                              'Lutimos',
                                                              'Flavicella',
                                                              'Lewinella')), color=Association)) +
  geom_point() +
  geom_pointrange(size=1 ,aes(xmin=average-stdev, xmax=average+stdev)) + theme_bw() +
  geom_vline(xintercept=c(0,0), linetype="dotted") +
  scale_color_manual(values=c("yellow3", "darkblue")) +
  xlab("Regression Coefficient") + ylab("Genus") +
  theme(axis.text.y = element_text(face="italic"),
        axis.text = element_text(size=15),
        axis.title = element_text(size=18),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12)) +
  xlim(-4, 3.5) +
  annotate("text", x = -3, y = 13, label = "A", size=12)



ggplot(corn.ready.depth, aes(x=average, y=factor(genus, level=c('Rubripirellula',
                                                                'Blastopirellula',
                                                                'Halomos',
                                                                'OM43_clade',
                                                                'Psychromos',
                                                                'Pseudoalteromos',
                                                                'Lutimos',
                                                                'Flavicella',
                                                                'Lewinella')), color=Association)) +
  geom_point() +
  geom_pointrange(size=1, aes(xmin=average-stdev, xmax=average+stdev)) + theme_bw() +
  geom_vline(xintercept=c(0,0), linetype="dotted") +
  scale_color_manual(values=c("forestgreen", "purple2")) +
  xlab("Regression Coefficient") + ylab("Genus") +
  theme(axis.text.y = element_text(face="italic"),
        axis.text = element_text(size=15),
        axis.title = element_text(size=18),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12)) +
  xlim(-4, 3.5) +
  annotate("text", x = -3, y = 9, label = "B", size=12)


#ASV1065 - model worked favors mature and subsurface
diff.ab.age.5187 <- bbdml(formula = ASV5187 ~ age, phi.formula= ~ age,
                          data=phylo.object.genus)
plot(diff.ab.age.5187, color = "age", B=50)
summary(diff.ab.age.5187)
diff.ab.depth.1065 <- bbdml(formula = ASV1065 ~ depth, phi.formula= ~ depth,
                            data=phylo.object.genus)
plot(diff.ab.depth.1065, color = "depth", B=50)
summary(diff.ab.depth.1065)

##After meeting with Lizzy, it's been decided to calculate the following:
#1) Prevalence calculation (e.g., 15 samples)
#2) Subset prevalent genera and look at corncob
#3) Determine which data are extremely different and piece those out as well


##Lizzy has suggested subsetting only mature and doing a depth-wise comparison. --------
mat.metadata <- read.csv("Mature_metadata_juvenileandmature.csv")
View(mat.metadata)
speciescoords <- read.csv("species_coords.csv")
View(speciescoords)
mat.rawreads <- read.csv("Mature_ASVRawReads with Euks and NA removed.csv")
View(mat.rawreads)
mat.seqtab <- mat.rawreads %>% tibble::column_to_rownames("UniqueID")
mat.seqtab.matrix <- as.matrix(mat.seqtab)
View(mat.seqtab.matrix)
mat.seqtab.matrix.table <- otu_table(mat.seqtab.matrix, taxa_are_rows = FALSE)
mat.seqtab.matrix.table

mat.taxa <- speciescoords %>% tibble::column_to_rownames("ASV")
mat.taxa.matrix <- as.matrix(mat.taxa)
mat.taxa.matrix.table <- tax_table(mat.taxa.matrix)
View(mat.taxa.matrix.table)

mat.metadata <- mat.metadata %>% tibble::column_to_rownames("UniqueID")
mat.metadata.sample <- sample_data(mat.metadata)
mat.metadata.sample

mat.phylo.object <- phyloseq(mat.seqtab.matrix.table, mat.metadata.sample, mat.taxa.matrix.table)
mat.phylo.object

mat.phylo.object.genus <- mat.phylo.object %>%
  tax_glom("genus")
mat.phylo.object.genus
mat.Taxa.sheet <- as.data.frame(tax_table(mat.phylo.object.genus))
View(mat.Taxa.sheet)

mat.diff.ab.depth.genus <- differentialTest(formula = ~ depth,
                                  phi.formula = ~ depth,
                                  formula_null = ~1, 
                                  phi.formula_null = ~ depth,
                                  data=mat.phylo.object,
                                  test="Wald", fdr_cutoff=0.05)
plot(mat.diff.ab.depth.genus)
plot(mat.diff.ab.depth.genus, level=c("genus")) 
summary(mat.diff.ab.depth.genus)
mat.diff.ab.depth.genus
otu_to_taxonomy(OTU = mat.diff.ab.depth.genus$significant_taxa, data=mat.phylo.object.genus)
mat.diff.ab.depth$p_fdr
diff.ab.depth.genus$significant_models
diff.ab.depth.genus$significant_taxa

##Replotting corncob coordinates------------
corncob.new.new <- read.csv("062123_Updated_CorncobCoords_ForPlotting.csv")
View(corncob.new.new)
corncob.new.age <- subset(corncob.new, AnalysisType=="Age")
View(corncob.new.age)
ggplot(Corncob.age.redone, aes(x=Average, y=factor(ASV, level=c('Roseibacillus ASV6762',
                                                              'Lentimonas ASV7032',
                                                              'Rubripirellula ASV6855',
                                                              'Blastopirellula ASV7364',
                                                              'Blastopirellula ASV4953',
                                                              'Blastopirellula ASV4017',
                                                              'Blastopirellula ASV3712',
                                                              'Blastopirellula ASV2621',
                                                              'Blastopirellula ASV1641',
                                                              'Blastopirellula ASV185',
                                                              'Psychromonas ASV3704',
                                                              'Psychromonas ASV1711',
                                                              'Lewinella ASV3262',
                                                              'Lewinella ASV1031')), color=Association)) +
  geom_point() +
  geom_pointrange(size=1 ,aes(xmin=Average-Stdev, xmax=Average+Stdev)) + theme_bw() +
  geom_vline(xintercept=c(0,0), linetype="dotted") +
  scale_color_manual(values=c("yellow3", "darkblue")) +
  xlab("Regression Coefficient") + ylab("Genus") +
  theme(axis.text.y = element_text(face="italic"),
        axis.text = element_text(size=15),
        axis.title = element_text(size=18),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))

Corncob.asv.redone <- read.csv("ASV_ReDone_ForBiplot_LizzyMeeting.csv")
View(Corncob.asv.redone)
Corncob.age.redone <- subset(Corncob.asv.redone, Comparison=="Age")
View(Corncob.age.redone)
ggplot(Corncob.age.redone, aes(x=Average, y=factor(ASV, level=c('Blastopirellula ASV4017',
                                                                  'Blastopirellula ASV3712',
                                                                  'Blastopirellula ASV2621',
                                                                  'Blastopirellula ASV1641',
                                                                  'Blastopirellula ASV185',
                                                                  'Psychromonas ASV3704',
                                                                  'Psychromonas ASV1711',
                                                                  'Lewinella ASV3262',
                                                                'Lewinella ASV1031')), color=Association)) +
  geom_point() +
  geom_pointrange(size=1 ,aes(xmin=Average-Stdev, xmax=Average+Stdev)) + theme_bw() +
  geom_vline(xintercept=c(0,0), linetype="dotted") +
  scale_color_manual(values=c("yellow3", "darkblue")) +
  xlab("Regression Coefficient") + ylab("ASV") +
  theme(axis.text.y = element_text(face="italic"),
        axis.text = element_text(size=15),
        axis.title = element_text(size=18),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))


corncob.new.depth <- subset(corncob.new, AnalysisType=="Depth")
View(corncob.new.depth)
ggplot(corncob.new.depth, aes(x=average, y=factor(genus, level=c('Blastopirellula',
                                                                 'Granulosicoccus',
                                                                 'Pseudoalteromonas',
                                                                 'Leucothrix',
                                                                 'Halomonas',
                                                                 'Propionigenium',
                                                                 'Synechococcus_CC9902',
                                                                 'Flavicella',
                                                                 'Afipia')), color=Association)) +
  geom_point() +
  geom_pointrange(size=1, aes(xmin=average-stdev, xmax=average+stdev)) + theme_bw() +
  geom_vline(xintercept=c(0,0), linetype="dotted") +
  scale_color_manual(values=c("forestgreen", "purple2")) +
  xlab("Regression Coefficient") + ylab("Genus") +
  theme(axis.text.y = element_text(face="italic"),
        axis.text = element_text(size=15),
        axis.title = element_text(size=18),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))

Corncob.depth.redone <- subset(Corncob.asv.redone, Comparison=="Depth")
View(Corncob.depth.redone)
ggplot(Corncob.depth.redone, aes(x=Average, y=factor(ASV, level=c('Blastopirellula ASV5806',
                                                               'Blastopirellula ASV1641',
                                                               'Blastopirellula ASV185',
                                                               'Pseudoalteromonas ASV1065',
                                                               'Pseudoalteromonas ASV656',
                                                               'Propionigenium ASV631')), color=Association)) +
  geom_point() +
  geom_pointrange(size=1 ,aes(xmin=Average-Stdev, xmax=Average+Stdev)) + theme_bw() +
  geom_vline(xintercept=c(0,0), linetype="dotted") +
  scale_color_manual(values=c("forestgreen", "purple2")) +
  xlab("Regression Coefficient") + ylab("ASV") +
  theme(axis.text.y = element_text(face="italic"),
        axis.text = element_text(size=15),
        axis.title = element_text(size=18),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))


corncob.new.depth <- subset(corncob.new, AnalysisType=="Depth")
View(corncob.new.depth)
ggplot(corncob.new.depth, aes(x=average, y=factor(genus, level=c('Blastopirellula',
                                                                 'Granulosicoccus',
                                                                 'Pseudoalteromonas',
                                                                 'Leucothrix',
                                                                 'Halomonas',
                                                                 'Propionigenium',
                                                                 'Synechococcus_CC9902',
                                                                 'Flavicella',
                                                                 'Afipia')), color=Association)) +
  geom_point() +
  geom_pointrange(size=1, aes(xmin=average-stdev, xmax=average+stdev)) + theme_bw() +
  geom_vline(xintercept=c(0,0), linetype="dotted") +
  scale_color_manual(values=c("forestgreen", "purple2")) +
  xlab("Regression Coefficient") + ylab("Genus") +
  theme(axis.text.y = element_text(face="italic"),
        axis.text = element_text(size=15),
        axis.title = element_text(size=18),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))


#Individual ASV testing--------
#Blastopirellula------
ASV2 <- bbdml(formula=ASV8~depth, phi.formula=~depth, data=mat.phylo.object)
summary(OM43.age)
plot(OM43.age, color="age")
ASV2
ASV8 <- bbdml(formula=ASV8~depth, phi.formula=~depth, data=mat.phylo.object)

ASV10 <- bbdml(formula=ASV10~depth, phi.formula=~depth, data=mat.phylo.object)
ASV11 <- bbdml(formula=ASV11~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV11)
plot(ASV11, color="depth")
ASV20 <- bbdml(formula=ASV20~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV20)
plot(ASV20, color="depth")
ASV34 <- bbdml(formula=ASV34~depth, phi.formula=~depth, data=mat.phylo.object)
ASV52 <- bbdml(formula=ASV52~depth, phi.formula=~depth, data=mat.phylo.object)
ASV83 <- bbdml(formula=ASV83~depth, phi.formula=~depth, data=mat.phylo.object)
ASV99 <- bbdml(formula=ASV99~depth, phi.formula=~depth, data=mat.phylo.object)
ASV185 <- bbdml(formula=ASV185~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV185)
plot(ASV185, color="depth")
ASV225 <- bbdml(formula=ASV225~depth, phi.formula=~depth, data=mat.phylo.object)
ASV244 <- bbdml(formula=ASV244~depth, phi.formula=~depth, data=mat.phylo.object)
ASV274 <- bbdml(formula=ASV274~depth, phi.formula=~depth, data=mat.phylo.object)
ASV285 <- bbdml(formula=ASV285~depth, phi.formula=~depth, data=mat.phylo.object)
ASV332 <- bbdml(formula=ASV332~depth, phi.formula=~depth, data=mat.phylo.object)
ASV338 <- bbdml(formula=ASV338~depth, phi.formula=~depth, data=mat.phylo.object)
ASV377 <- bbdml(formula=ASV377~depth, phi.formula=~depth, data=mat.phylo.object)
ASV427 <- bbdml(formula=ASV427~depth, phi.formula=~depth, data=mat.phylo.object)
ASV565 <- bbdml(formula=ASV565~depth, phi.formula=~depth, data=mat.phylo.object)
ASV679 <- bbdml(formula=ASV679~depth, phi.formula=~depth, data=mat.phylo.object)
ASV696 <- bbdml(formula=ASV696~depth, phi.formula=~depth, data=mat.phylo.object)
ASV703 <- bbdml(formula=ASV703~depth, phi.formula=~depth, data=mat.phylo.object)
ASV707 <- bbdml(formula=ASV707~depth, phi.formula=~depth, data=mat.phylo.object)
ASV725 <- bbdml(formula=ASV725~depth, phi.formula=~depth, data=mat.phylo.object)
ASV754 <- bbdml(formula=ASV754~depth, phi.formula=~depth, data=mat.phylo.object)
ASV762 <- bbdml(formula=ASV762~depth, phi.formula=~depth, data=mat.phylo.object)
ASV796 <- bbdml(formula=ASV796~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV796)
plot(ASV796, color="depth")
ASV797 <- bbdml(formula=ASV797~depth, phi.formula=~depth, data=mat.phylo.object)
ASV871 <- bbdml(formula=ASV871~depth, phi.formula=~depth, data=mat.phylo.object)
ASV959 <- bbdml(formula=ASV959~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV959)
plot(ASV959, color="depth")
ASV996 <- bbdml(formula=ASV996~depth, phi.formula=~depth, data=mat.phylo.object)
ASV1077 <- bbdml(formula=ASV1077~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV1077)
ASV1319 <- bbdml(formula=ASV1319~depth, phi.formula=~depth, data=mat.phylo.object)
ASV1351 <- bbdml(formula=ASV1351~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV1351)
plot(ASV1351, color="depth")
ASV1409 <- bbdml(formula=ASV1409~depth, phi.formula=~depth, data=mat.phylo.object)
ASV1641 <- bbdml(formula=ASV1641~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV1641)
plot(ASV1641, color="depth")
ASV1838 <- bbdml(formula=ASV1838~depth, phi.formula=~depth, data=mat.phylo.object)
ASV1856 <- bbdml(formula=ASV1856~depth, phi.formula=~depth, data=mat.phylo.object)
ASV1858 <- bbdml(formula=ASV1858~depth, phi.formula=~depth, data=mat.phylo.object)
ASV2095 <- bbdml(formula=ASV2095~depth, phi.formula=~depth, data=mat.phylo.object)
ASV2114 <- bbdml(formula=ASV2114~depth, phi.formula=~depth, data=mat.phylo.object)
ASV2238 <- bbdml(formula=ASV2238~depth, phi.formula=~depth, data=mat.phylo.object)
ASV2254 <- bbdml(formula=ASV2254~depth, phi.formula=~depth, data=mat.phylo.object)
ASV2283 <- bbdml(formula=ASV2283~depth, phi.formula=~depth, data=mat.phylo.object)
ASV2366 <- bbdml(formula=ASV2366~depth, phi.formula=~depth, data=mat.phylo.object)
ASV2621 <- bbdml(formula=ASV2621~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV2621)
ASV2748 <- bbdml(formula=ASV2748~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV2748)
ASV2760 <- bbdml(formula=ASV2760~depth, phi.formula=~depth, data=mat.phylo.object)
ASV2793 <- bbdml(formula=ASV2793~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV2793)
plot(ASV2793, color="depth")
ASV2794 <- bbdml(formula=ASV2794~depth, phi.formula=~depth, data=mat.phylo.object)
ASV3187 <- bbdml(formula=ASV3187~depth, phi.formula=~depth, data=mat.phylo.object)
ASV3231 <- bbdml(formula=ASV3231~depth, phi.formula=~depth, data=mat.phylo.object)
ASV3255 <- bbdml(formula=ASV3255~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV3255)
plot(ASV3255, color="depth")
ASV3485 <- bbdml(formula=ASV3485~depth, phi.formula=~depth, data=mat.phylo.object)
ASV3674 <- bbdml(formula=ASV3674~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV3674)
ASV3712 <- bbdml(formula=ASV3712~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV3712)
plot(ASV3712, color="depth")
ASV3886 <- bbdml(formula=ASV3886~depth, phi.formula=~depth, data=mat.phylo.object)
ASV4017 <- bbdml(formula=ASV4017~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV4017)
ASV4044 <- bbdml(formula=ASV4044~depth, phi.formula=~depth, data=mat.phylo.object)
ASV4164 <- bbdml(formula=ASV4164~depth, phi.formula=~depth, data=mat.phylo.object)
ASV4272 <- bbdml(formula=ASV4272~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV4272)
ASV4530 <- bbdml(formula=ASV4530~depth, phi.formula=~depth, data=mat.phylo.object)
ASV4651 <- bbdml(formula=ASV4651~depth, phi.formula=~depth, data=mat.phylo.object)
ASV4731 <- bbdml(formula=ASV4731~depth, phi.formula=~depth, data=mat.phylo.object)
ASV4771 <- bbdml(formula=ASV4771~depth, phi.formula=~depth, data=mat.phylo.object)
ASV4772 <- bbdml(formula=ASV4772~depth, phi.formula=~depth, data=mat.phylo.object)
ASV4774 <- bbdml(formula=ASV4774~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV4774)
ASV4805 <- bbdml(formula=ASV4805~depth, phi.formula=~depth, data=mat.phylo.object)
ASV4830 <- bbdml(formula=ASV4830~depth, phi.formula=~depth, data=mat.phylo.object)
ASV4953 <- bbdml(formula=ASV4953~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV4953)
ASV5002 <- bbdml(formula=ASV5002~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV5002)
ASV5075 <- bbdml(formula=ASV5075~depth, phi.formula=~depth, data=mat.phylo.object)
ASV5090 <- bbdml(formula=ASV5090~depth, phi.formula=~depth, data=mat.phylo.object)
ASV5113 <- bbdml(formula=ASV5113~depth, phi.formula=~depth, data=mat.phylo.object)
ASV5234 <- bbdml(formula=ASV5234~depth, phi.formula=~depth, data=mat.phylo.object)
ASV5279 <- bbdml(formula=ASV5279~depth, phi.formula=~depth, data=mat.phylo.object)
ASV5518 <- bbdml(formula=ASV5518~depth, phi.formula=~depth, data=mat.phylo.object)
ASV5652 <- bbdml(formula=ASV5652~depth, phi.formula=~depth, data=mat.phylo.object)
ASV5657 <- bbdml(formula=ASV5657~depth, phi.formula=~depth, data=mat.phylo.object)
ASV5703 <- bbdml(formula=ASV5703~depth, phi.formula=~depth, data=mat.phylo.object)
ASV5806 <- bbdml(formula=ASV5806~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV5806)
plot(ASV5806, color="depth")
ASV5810 <- bbdml(formula=ASV5810~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV5810)
ASV5879 <- bbdml(formula=ASV5879~depth, phi.formula=~depth, data=mat.phylo.object)
ASV6366 <- bbdml(formula=ASV6366~depth, phi.formula=~depth, data=mat.phylo.object)
ASV6416 <- bbdml(formula=ASV6416~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV6416)
ASV6457 <- bbdml(formula=ASV6457~depth, phi.formula=~depth, data=mat.phylo.object)
ASV6679 <- bbdml(formula=ASV6679~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV6679)
ASV6778 <- bbdml(formula=ASV6778~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV6778)
ASV7007 <- bbdml(formula=ASV7007~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV7007)
plot(ASV7007, color="depth")
ASV7131 <- bbdml(formula=ASV7131~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV7131)
ASV7351 <- bbdml(formula=ASV7351~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV7351)
ASV7364 <- bbdml(formula=ASV7364~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV7364)
ASV7432 <- bbdml(formula=ASV7432~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV7432)
plot(ASV7432, color="depth")
ASV7814 <- bbdml(formula=ASV7814~depth, phi.formula=~depth, data=mat.phylo.object)
ASV7933 <- bbdml(formula=ASV7933~depth, phi.formula=~depth, data=mat.phylo.object)
ASV7999 <- bbdml(formula=ASV7999~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV7999)
plot(ASV7999, color="depth")
ASV8150 <- bbdml(formula=ASV8150~depth, phi.formula=~depth, data=mat.phylo.object)
ASV8189 <- bbdml(formula=ASV8189~depth, phi.formula=~depth, data=mat.phylo.object)
ASV8312 <- bbdml(formula=ASV8312~depth, phi.formula=~depth, data=mat.phylo.object)
ASV8458 <- bbdml(formula=ASV8458~depth, phi.formula=~depth, data=mat.phylo.object)
ASV9166 <- bbdml(formula=ASV9166~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV9166)
ASV9364 <- bbdml(formula=ASV9364~depth, phi.formula=~depth, data=mat.phylo.object)
ASV9513 <- bbdml(formula=ASV9513~depth, phi.formula=~depth, data=mat.phylo.object)
ASV9601 <- bbdml(formula=ASV9601~depth, phi.formula=~depth, data=mat.phylo.object)
ASV9603 <- bbdml(formula=ASV9603~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV9603)
plot(ASV9603, color="depth")
ASV10006 <- bbdml(formula=ASV10006~depth, phi.formula=~depth, data=mat.phylo.object)
ASV10070 <- bbdml(formula=ASV10070~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV10070)
ASV10073 <- bbdml(formula=ASV10073~depth, phi.formula=~depth, data=mat.phylo.object)
ASV10159 <- bbdml(formula=ASV10159~depth, phi.formula=~depth, data=mat.phylo.object)
ASV10345 <- bbdml(formula=ASV10345~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV10345)
ASV10471 <- bbdml(formula=ASV10471~depth, phi.formula=~depth, data=mat.phylo.object)
ASV10672 <- bbdml(formula=ASV10672~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV10672)
ASV10757 <- bbdml(formula=ASV10757~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV10757)
plot(ASV10757, color="depth")
ASV11294 <- bbdml(formula=ASV11294~depth, phi.formula=~depth, data=mat.phylo.object)
ASV11867 <- bbdml(formula=ASV11867~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV11867)
ASV12241 <- bbdml(formula=ASV12241~depth, phi.formula=~depth, data=mat.phylo.object)

#Granulosicoccus------
ASV31 <-bbdml(formula=ASV31~depth, phi.formula=~depth, data=mat.phylo.object)
ASV36 <- bbdml(formula=ASV36~depth, phi.formula=~depth, data=mat.phylo.object)
ASV51 <- bbdml(formula=ASV51~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV51)
ASV100 <- bbdml(formula=ASV100~depth, phi.formula=~depth, data=mat.phylo.object)
ASV123 <- bbdml(formula=ASV123~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV123)
ASV124 <- bbdml(formula=ASV124~depth, phi.formula=~depth, data=mat.phylo.object)
ASV129 <- bbdml(formula=ASV129~depth, phi.formula=~depth, data=mat.phylo.object)
ASV178 <- bbdml(formula=ASV178~depth, phi.formula=~depth, data=mat.phylo.object)
ASV277 <- bbdml(formula=ASV277~depth, phi.formula=~depth, data=mat.phylo.object)
ASV305 <- bbdml(formula=ASV305~depth, phi.formula=~depth, data=mat.phylo.object)
ASV334 <- bbdml(formula=ASV334~depth, phi.formula=~depth, data=mat.phylo.object)
ASV360 <- bbdml(formula=ASV360~depth, phi.formula=~depth, data=mat.phylo.object)
ASV382 <- bbdml(formula=ASV382~depth, phi.formula=~depth, data=mat.phylo.object)
ASV1281 <- bbdml(formula=ASV1281~depth, phi.formula=~depth, data=mat.phylo.object)
ASV1531 <- bbdml(formula=ASV1531~depth, phi.formula=~depth, data=mat.phylo.object)
ASV1737 <- bbdml(formula=ASV1737~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV1737)
ASV1792 <- bbdml(formula=ASV1792~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV1792)
ASV1957 <- bbdml(formula=ASV1957~depth, phi.formula=~depth, data=mat.phylo.object)
ASV2422 <- bbdml(formula=ASV2422~depth, phi.formula=~depth, data=mat.phylo.object)
ASV2425 <- bbdml(formula=ASV2425~depth, phi.formula=~depth, data=mat.phylo.object)
ASV3343 <- bbdml(formula=ASV3343~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV3343)
plot(ASV3343, color="depth")
ASV3889 <- bbdml(formula=ASV3889~depth, phi.formula=~depth, data=mat.phylo.object)
ASV7912 <- bbdml(formula=ASV7912~depth, phi.formula=~depth, data=mat.phylo.object)
ASV8069 <- bbdml(formula=ASV8069~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV8069)
plot(ASV8069, color="depth")
ASV10359 <- bbdml(formula=ASV10359~depth, phi.formula=~depth, data=mat.phylo.object)
ASV12636 <- bbdml(formula=ASV12636~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV12636)

#Afipia -------
ASV30 <- bbdml(formula=ASV30~depth, phi.formula=~depth, data=mat.phylo.object)
ASV436 <- bbdml(formula=ASV436~depth, phi.formula=~depth, data=mat.phylo.object)
ASV1155 <- bbdml(formula=ASV1155~depth, phi.formula=~depth, data=mat.phylo.object)
ASV2359 <- bbdml(formula=ASV2359~depth, phi.formula=~depth, data=mat.phylo.object)
ASV4064 <- bbdml(formula=ASV4064~depth, phi.formula=~depth, data=mat.phylo.object)
ASV10367 <- bbdml(formula=ASV10367~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV10367)
ASV12050 <- bbdml(formula=ASV12050~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV12050)
ASV436 <- bbdml(formula=ASV436~depth, phi.formula=~depth, data=mat.phylo.object)
ASV1155 <- bbdml(formula=ASV1155~depth, phi.formula=~depth, data=mat.phylo.object)
ASV2359 <- bbdml(formula=ASV2359~depth, phi.formula=~depth, data=mat.phylo.object)
ASV4064 <- bbdml(formula=ASV4064~depth, phi.formula=~depth, data=mat.phylo.object)
ASV10367 <- bbdml(formula=ASV10367~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV10367)
ASV12050 <- bbdml(formula=ASV12050~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV12050)
#Flavicella --------
ASV42 <- bbdml(formula=ASV42~depth, phi.formula=~depth, data=mat.phylo.object)
ASV161 <- bbdml(formula=ASV161~depth, phi.formula=~depth, data=mat.phylo.object)
ASV275 <- bbdml(formula=ASV275~depth, phi.formula=~depth, data=mat.phylo.object)
ASV863 <- bbdml(formula=ASV863~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV863)
ASV1267 <- bbdml(formula=ASV1267~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV1267)
ASV2208 <- bbdml(formula=ASV2208~depth, phi.formula=~depth, data=mat.phylo.object)
ASV2577 <- bbdml(formula=ASV2577~depth, phi.formula=~depth, data=mat.phylo.object)
ASV2852 <- bbdml(formula=ASV2852~depth, phi.formula=~depth, data=mat.phylo.object)
ASV5115 <- bbdml(formula=ASV5115~depth, phi.formula=~depth, data=mat.phylo.object)
ASV5155 <- bbdml(formula=ASV5155~depth, phi.formula=~depth, data=mat.phylo.object)
ASV5312 <- bbdml(formula=ASV5312~depth, phi.formula=~depth, data=mat.phylo.object)
ASV6356 <- bbdml(formula=ASV6356~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV6356)
plot(ASV6356, color="depth")
ASV6543 <- bbdml(formula=ASV6543~depth, phi.formula=~depth, data=mat.phylo.object)
ASV6610 <- bbdml(formula=ASV6610~depth, phi.formula=~depth, data=mat.phylo.object)
ASV6813 <- bbdml(formula=ASV6813~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV6813)
ASV7669 <- bbdml(formula=ASV7669~depth, phi.formula=~depth, data=mat.phylo.object)
ASV7740 <- bbdml(formula=ASV7740~depth, phi.formula=~depth, data=mat.phylo.object)
ASV8421 <- bbdml(formula=ASV8421~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV8421)
plot(ASV8421, color="depth")
ASV12779 <- bbdml(formula=ASV12779~depth, phi.formula=~depth, data=mat.phylo.object)

#Synechococcus_CC9902---------
ASV5 <- bbdml(formula=ASV5~depth, phi.formula=~depth, data=mat.phylo.object)
ASV9 <- bbdml(formula=ASV9~depth, phi.formula=~depth, data=mat.phylo.object)
ASV16 <- bbdml(formula=ASV16~depth, phi.formula=~depth, data=mat.phylo.object)
ASV262 <- bbdml(formula=ASV262~depth, phi.formula=~depth, data=mat.phylo.object)
ASV459 <- bbdml(formula=ASV459~depth, phi.formula=~depth, data=mat.phylo.object)
ASV1002 <- bbdml(formula=ASV1002~depth, phi.formula=~depth, data=mat.phylo.object)
ASV1225 <- bbdml(formula=ASV1225~depth, phi.formula=~depth, data=mat.phylo.object)
ASV6547 <- bbdml(formula=ASV6547~depth, phi.formula=~depth, data=mat.phylo.object)
ASV7672 <- bbdml(formula=ASV7672~depth, phi.formula=~depth, data=mat.phylo.object)
ASV7846 <- bbdml(formula=ASV7846~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV7846)
plot(ASV7846, color="depth")

#Propionigenium---------
ASV60 <- bbdml(formula=ASV60~depth, phi.formula=~depth, data=mat.phylo.object)
ASV98 <- bbdml(formula=ASV98~depth, phi.formula=~depth, data=mat.phylo.object)
ASV396 <- bbdml(formula=ASV396~depth, phi.formula=~depth, data=mat.phylo.object)
ASV596 <- bbdml(formula=ASV596~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV596)
ASV631 <- bbdml(formula=ASV631~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV631)
plot(ASV631, color="depth")

#Halomonas -------
ASV3 <- bbdml(formula=ASV3~depth, phi.formula=~depth, data=mat.phylo.object)
ASV13 <- bbdml(formula=ASV13~depth, phi.formula=~depth, data=mat.phylo.object)
ASV2816 <- bbdml(formula=ASV2816~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV2816)
plot(ASV2816, color="depth")
ASV4130 <- bbdml(formula=ASV4130~depth, phi.formula=~depth, data=mat.phylo.object)
ASV4878 <- bbdml(formula=ASV4878~depth, phi.formula=~depth, data=mat.phylo.object)

#Leucothrix --------
ASV32 <- bbdml(formula=ASV32~depth, phi.formula=~depth, data=mat.phylo.object)
ASV33 <- bbdml(formula=ASV33~depth, phi.formula=~depth, data=mat.phylo.object)
ASV94 <- bbdml(formula=ASV94~depth, phi.formula=~depth, data=mat.phylo.object)
ASV136 <- bbdml(formula=ASV136~depth, phi.formula=~depth, data=mat.phylo.object)
ASV447 <- bbdml(formula=ASV447~depth, phi.formula=~depth, data=mat.phylo.object)
ASV490 <- bbdml(formula=ASV490~depth, phi.formula=~depth, data=mat.phylo.object)
ASV738 <- bbdml(formula=ASV738~depth, phi.formula=~depth, data=mat.phylo.object)
ASV859 <- bbdml(formula=ASV859~depth, phi.formula=~depth, data=mat.phylo.object)
ASV1058 <- bbdml(formula=ASV1058~depth, phi.formula=~depth, data=mat.phylo.object)
ASV2751 <- bbdml(formula=ASV2751~depth, phi.formula=~depth, data=mat.phylo.object)
ASV2893 <- bbdml(formula=ASV2893~depth, phi.formula=~depth, data=mat.phylo.object)
ASV3008 <- bbdml(formula=ASV3008~depth, phi.formula=~depth, data=mat.phylo.object)
ASV3399 <- bbdml(formula=ASV3399~depth, phi.formula=~depth, data=mat.phylo.object)
ASV5913 <- bbdml(formula=ASV5913~depth, phi.formula=~depth, data=mat.phylo.object)
ASV7300 <- bbdml(formula=ASV7300~depth, phi.formula=~depth, data=mat.phylo.object)
ASV10805 <- bbdml(formula=ASV10805~depth, phi.formula=~depth, data=mat.phylo.object)
ASV12067 <- bbdml(formula=ASV12067~depth, phi.formula=~depth, data=mat.phylo.object)
ASV12844 <- bbdml(formula=ASV12844~depth, phi.formula=~depth, data=mat.phylo.object)
summary(ASV12844)
plot(ASV12844, color="depth")

#NMDS setup-------------
All.NMDS <- read.csv("UniqueID NMDS coordinates for biplot.csv")
View(All.NMDS)
ASV.Age <- read.csv("AgeNMDS Biplot After Meeting with Lizzy 071023.csv")
View(ASV.Age)



ggplot(All.NMDS) +
  geom_point(mapping = aes(x=NMDS1, y=NMDS2, shape=Age, color=Age), size=6, alpha=0.4) +
  theme_bw() +
  scale_color_manual(values=c("yellow3", "darkblue")) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text=element_text(size=20),
        legend.title = element_text(size=20)) +
  geom_segment(data=ASV.Age,
               aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               type="closed", color="gray", alpha=0.5) +
  geom_text(data = ASV.Age, aes(x = NMDS1, y = NMDS2, label = ASV),
            size = 4) + xlim(-1.3, 1.5)
