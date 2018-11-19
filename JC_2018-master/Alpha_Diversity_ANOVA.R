#Edited July5 2018
#Alpha Diversity Analysis for Ecology MS

#Import spreadsheets with alpha diversity data

library(readr)
ROOTSalpha_diversity_analyses_nomito_nochloro_forANOVA <- read_delim("~/Dropbox/MS_Jatropha_Ecology/Code/AlphaDiversityANOVA/ROOTSalpha_diversity_analyses_nomito_nochloro_forANOVA.txt", 
                                                                     "\t", escape_double = FALSE, trim_ws = TRUE)
View(ROOTSalpha_diversity_analyses_nomito_nochloro_forANOVA)

LEAVESalpha_divesity_analyses_nomito_nochloro_forANOVA <- read_delim("~/Dropbox/MS_Jatropha_Ecology/Code/AlphaDiversityANOVA/LEAVESalpha_divesity_analyses_nomito_nochloro_forANOVA.txt", 
                                                                     "\t", escape_double = FALSE, trim_ws = TRUE)
View(LEAVESalpha_divesity_analyses_nomito_nochloro_forANOVA)

SOILalpha_diversity_analyses_nomito_nochloro_forANOVA <- read_delim("~/Dropbox/MS_Jatropha_Ecology/Final code/AlphaDiversityANOVA/SOILalpha_diversity_analyses_nomito_nochloro_forANOVA.txt", 
                                                                    "\t", escape_double = FALSE, trim_ws = TRUE)
View(SOILalpha_diversity_analyses_nomito_nochloro_forANOVA)

all_live_alpha_forANOVA <- read_delim("~/Dropbox/MS_Jatropha_Ecology/Final code/AlphaDiversityANOVA/all_live_alpha_forANOVA.txt", 
                                      "\t", escape_double = FALSE, col_types = cols(Description = col_factor(levels = c("Leaf", 
                                                                                                                        "Root", "Soil")), PD_whole_tree = col_number(), 
                                                                                    SampleID = col_character(), Site = col_factor(levels = c("Asientos", 
                                                                                                                                             "Ejido", "Divisa")), chao1 = col_number()), 
                                      trim_ws = TRUE)
View(all_live_alpha_forANOVA)

all_live_alpha_forANOVA_transformed <- read_delim("~/Dropbox/MS_Jatropha_Ecology/Final code/AlphaDiversityANOVA/all_live_alpha_forANOVA_transformed.txt", 
                                                  "\t", escape_double = FALSE, col_types = cols(LNChao1 = col_number(), 
                                                                                                PD_whole_tree = col_number(), SampleID = col_character(), 
                                                                                                Site = col_factor(levels = c("Asientos", 
                                                                                                                             "Ejido", "Divisa")), chao1 = col_number()), 
                                                  trim_ws = TRUE)
View(all_live_alpha_forANOVA_transformed)

#All soil,leaves, and roots (from live soil)
AllANOVAdata<-(all_live_alpha_forANOVA)

bartlett.test(chao1 ~ interaction(Description,Site), data = AllANOVAdata)
#	Bartlett test of homogeneity of variances
#data:  chao1 by interaction(Description, Site)
#Bartlett's K-squared = 57.359, df = 8, p-value = 1.532e-09
#Unequal variance

bartlett.test(chao1 ~ Description, data = AllANOVAdata)
#Bartlett's K-squared = 36.378, df = 2, p-value = 1.261e-08

bartlett.test(chao1 ~ Site, data = AllANOVAdata)
#Bartlett's K-squared = 4.817, df = 2, p-value = 0.08995

kruskal.test(chao1 ~ Site, data=AllANOVAdata)
#Kruskal-Wallis chi-squared = 6.1062, df = 2, p-value = 0.04721

kruskal.test(chao1 ~ Description, data=AllANOVAdata)
#Kruskal-Wallis chi-squared = 117.49, df = 2, p-value < 2.2e-16

kruskal.test(chao1 ~ Description*Site, data=AllANOVAdata)

aov.out.all = aov(chao1 ~ Description * Site, data = AllANOVAdata)
options(show.signif.stars = F)
aov.out.all
summary(aov.out.all)
#                  Df    Sum Sq  Mean Sq F value   Pr(>F)
#Description        2 124315037 62157518 300.015  < 2e-16
#Site               2   5711779  2855889  13.784 3.22e-06
#Description:Site   4   2013274   503318   2.429   0.0502
#Residuals        149  30870033   207181  


#All soil, roots, and leaves transformed
AllANOVAdata_transformed<-(all_live_alpha_forANOVA_transformed)

bartlett.test(LNChao1 ~ interaction(Description,Site), data = AllANOVAdata_transformed)
#Not normal

library(car)
leveneTest(LNChao1 ~ Description*Site, data = AllANOVAdata_transformed)
#Not normal

#Leaves between greenhouse and field
LeafANOVAdata<-(LEAVESalpha_divesity_analyses_nomito_nochloro_forANOVA)
LeafANOVAdata
#Roughly even repliate numbers? Yes
replications(chao1 ~ Location * Soil, data = LeafANOVAdata)

#difference of variances? Not Significant. p=0.18
bartlett.test(chao1 ~ interaction(Location,Soil), data = LeafANOVAdata)
#Bartlett's K-squared = 7.5863, df = 5, p-value = 0.1806

aov.out.leaf = aov(chao1 ~ Location * Soil, data = LeafANOVAdata)
options(show.signif.stars = F)
aov.out.leaf
summary(aov.out.leaf)
#Df  Sum Sq Mean Sq F value  Pr(>F)
#Location       1 2237687 2237687  41.919 2.1e-07
#Soil           2  408068  204034   3.822  0.0318
#Location:Soil  2  413186  206593   3.870  0.0306
#Residuals     34 1814967   53381 

TukeyHSD(aov.out.leaf)
#$Location
#                   diff      lwr      upr p adj
#Greenhouse-Field 478.4553 328.2751 628.6354 2e-07

#$Soil
#                   diff        lwr       upr     p adj
#Divisa-Asientos -113.9704 -332.0346 104.09384 0.4155548
#Ejido-Asientos  -245.8914 -463.9557 -27.82721 0.0242689
#Ejido-Divisa    -131.9210 -353.9868  90.14469 0.3245852

#Graphical Representation
boxplot(chao1 ~ Location * Soil, data = LeafANOVAdata)
boxplot(chao1 ~ Soil, data = LeafANOVAdata) #least variance in Divisa
boxplot(chao1 ~ Location, data = LeafANOVAdata)
#intractions boxplot showing diversity higher in greenhouse except in Divisa
with(LeafANOVAdata, interaction.plot(x.factor=Soil, trace.factor = Location, response = chao1, fun=mean, type = "b", legend = T, ylab = "Chao1", main= "Leaf Interaction Plot", pch = c(1,19)))

#Leaf OTUS
#difference of variances? Significant. p=0.002
bartlett.test(observed_otus ~ interaction(Location,Soil), data = LeafANOVAdata)




#Roots
RootANOVAdata<-(ROOTSalpha_diversity_analyses_nomito_nochloro_forANOVA)
RootANOVAdata

#Roughly even repliate numbers? Too many field samples thows things off a bit
replications(chao1 ~ Location * Soil, data = RootANOVAdata) # 2x # of Field samples than GH in Divisa and Ejido
replications(chao1 ~ Soil, data = RootANOVAdata) #Similar numbers in each soil
replications(chao1 ~ Location, data = RootANOVAdata) # 2x as many field as greenhouse samples

#difference of variances? Not Significant. p=0.129
bartlett.test(chao1 ~ interaction(Location,Soil), data = RootANOVAdata)

aov.out.root = aov(chao1 ~ Location * Soil, data = RootANOVAdata)
options(show.signif.stars = F)
aov.out.root
summary(aov.out.root)
#              Df   Sum Sq Mean Sq F value  Pr(>F)
#Location       1  2521639 2521639   9.115 0.00451
#Soil           2  1639140  819570   2.962 0.06374
#Location:Soil  2  1132524  566262   2.047 0.14315
#Residuals     38 10513027  276659                

TukeyHSD(aov.out.root)
#$Location
#                   diff       lwr       upr     p adj
#Greenhouse-Field -513.9764 -858.6189 -169.3339 0.0045128

#$Soil
#                   diff       lwr       upr     p adj
#Divisa-Asientos   87.69567 -398.3920 573.78337 0.8990747
#Ejido-Asientos  -347.62512 -826.6083 131.35805 0.1933073
#Ejido-Divisa    -435.32079 -896.3497  25.70811 0.0675296


#Graphical Representation
boxplot(chao1 ~ Location * Soil, data = RootANOVAdata)
boxplot(chao1 ~ Soil, data = RootANOVAdata) #least variance in Divisa
boxplot(chao1 ~ Location, data = RootANOVAdata)
#intractions boxplot showing diversity much higher in Field except in El Ejido
with(RootANOVAdata, interaction.plot(x.factor=Soil, trace.factor = Location, response = chao1, fun=mean, type = "b", legend = T, ylab = "Chao1", main= " Root Interaction Plot", pch = c(1,19)))

#Root observed OTUs
#difference of variances? Not Significant. p=0.353
bartlett.test(observed_otus ~ interaction(Location,Soil), data = RootANOVAdata)
#Bartlett's K-squared = 5.5413, df = 5, p-value = 0.3534

aov.out.root.otus = aov(observed_otus ~ Location * Soil, data = RootANOVAdata)
options(show.signif.stars = F)
aov.out.root.otus
summary(aov.out.root.otus)
#> summary(aov.out.root.otus)
#               Df  Sum Sq Mean Sq F value Pr(>F)
#Location       1  693632  693632   4.084 0.0504
#Soil           2 1385331  692665   4.079 0.0249
#Location:Soil  2  868723  434362   2.558 0.0908

TukeyHSD(aov.out.root.otus)

#Figure
boxplot(observed_otus ~ Location * Soil, data = RootANOVAdata)
boxplot(observed_otus ~ Soil, data = RootANOVAdata) #
boxplot(observed_otus1 ~ Location, data = RootANOVAdata)
#
with(RootANOVAdata, interaction.plot(x.factor=Soil, trace.factor = Location, response = observed_otus, fun=mean, type = "b", legend = T, ylab = "observed_otus", main= " Root Interaction Plot", pch = c(1,19)))

#Soil
SoilANOVAdata<-(SOILalpha_diversity_analyses_nomito_nochloro_forANOVA)
SoilANOVAdata

bartlett.test(chao1 ~ interaction(Site,Treatment), data = SoilANOVAdata)
#Bartlett's K-squared = 1.4107, df = 5, p-value = 0.9231

aov.out.soil.chao1 = aov(chao1 ~ Site * Treatment, data = SoilANOVAdata)
options(show.signif.stars = F)
aov.out.soil.otus
summary(aov.out.soil.chao1)

boxplot(chao1 ~ Site * Treatment, data = SoilANOVAdata)
boxplot(chao1 ~ Site, data = SoilANOVAdata) 
boxplot(chao1 ~ Treatment, data = SoilANOVAdata)

with(SoilANOVAdata, interaction.plot(x.factor=Site, trace.factor = Treatment, response = chao1, fun=mean, type = "b", legend = T, ylab = "Chao1", main= "Soil Alpha", pch = c(1,19)))


################Figures for manuscript##############
library(ggplot2)

AllAlphaPlot <- ggplot(all_live_alpha_forANOVA, aes(x=Description, y=chao1, fill=Site)) + 
  geom_boxplot() + theme_linedraw() + labs(y= "Chao1") + scale_fill_grey(start=0.5, end=1)  + theme(text = element_text(size = 15))
AllAlphaPlot

LeavesAlphaPlot <- ggplot(LeafANOVAdata, aes(x=Soil, y=chao1, fill=Location)) + 
  geom_boxplot() + theme_linedraw() + labs(y= "Chao1") + scale_fill_grey(start=0.5, end=1) + theme(text = element_text(size = 15))
LeavesAlphaPlot

RootsAlphaPlot <- ggplot(RootANOVAdata, aes(x=Soil, y=chao1, fill=Location)) + 
  geom_boxplot() + theme_linedraw() + labs(y= "Chao1") + scale_fill_grey(start=0.5, end=1) + theme(text = element_text(size = 15))
RootsAlphaPlot


