#Bringing bioms from QIIME into phyloseq for making figures

# https://web.stanford.edu/class/bios221/labs/phyloseq/lab_phyloseq.html
#http://joey711.github.io/phyloseq/import-data#import_biom

setwd("/Users/kimberlymighell/Dropbox/JCALL_NOV16")

library(phyloseq)
library(ggplot2)

# Load color palettes
library(RColorBrewer)

#Define source files
otu_table_field_roots_biom<-"field_endosymbionts_soil/roots/otu_table_field_roots_nomito_nochloro.biom"

otu_table_field_leaves_biom<-"field_endosymbionts_soil/leaves/otu_table_field_leaves_nomito_nochloro.biom"

otu_table_GHF_leaves_biom <- "Leaves/GHFieldComp/otu_table_leaves_comp_nomito_nochloro.biom"

otu_table_GHF_roots_biom <- "Roots/GHFieldComp/otu_table_roots_comp_nomito_nochloro.biom"

otu_table_GH_leaves_biom <- "GH/all/leaves/live/otu_table_GH_all_leaves_nomito_nochloro.biom"

otu_table_GH_roots_biom <-"GH/all/roots/live/otu_table_GH_all_roots_nomito_nochloro.biom"

tree_file<-"rep_set.tre"


#Import Biom
#RF = all root DNA from field sites

RFbiom <- import_biom(otu_table_field_roots_biom)

LFbiom <- import_biom(otu_table_field_leaves_biom)

GHFLbiom <- import_biom(otu_table_GHF_leaves_biom)

GHFRbiom <- import_biom(otu_table_GHF_roots_biom)

GHRbiom <- import_biom(otu_table_GH_roots_biom)

GHLbiom <- import_biom(otu_table_GH_leaves_biom)

#import map file
jcall_map <-import_qiime_sample_data("merged_jcall_final.txt")

#Merge biom, map, tree into one phyloseq object

RF_merge<- merge_phyloseq(RFbiom,jcall_map,tree_file)
LF_merge<- merge_phyloseq(LFbiom,jcall_map,tree_file)

GHFL_merge<- merge_phyloseq(GHFLbiom,jcall_map,tree_file)

GHFR_merge<- merge_phyloseq(GHFRbiom,jcall_map,tree_file)

GHL_merge<- merge_phyloseq(GHLbiom,jcall_map,tree_file)

GHR_merge<- merge_phyloseq(GHRbiom,jcall_map,tree_file)

#Rename Ranks to Taxonomy Categorie

colnames(tax_table(RF_merge)) <- c(Rank1 = 'Kingdom', Rank2 = "Phylum", Rank3 = "Class", 
                                 Rank4 = "Order", Rank5 = "Family", 
                                 Rank6 = "Genus", Rank7 = "Species")

colnames(tax_table(LF_merge)) <- c(Rank1 = 'Kingdom', Rank2 = "Phylum", Rank3 = "Class", 
                                   Rank4 = "Order", Rank5 = "Family", 
                                   Rank6 = "Genus", Rank7 = "Species")

colnames(tax_table(GHFL_merge)) <- c(Rank1 = 'Kingdom', Rank2 = "Phylum", Rank3 = "Class", 
                                   Rank4 = "Order", Rank5 = "Family", 
                                   Rank6 = "Genus", Rank7 = "Species")

colnames(tax_table(GHFR_merge)) <- c(Rank1 = 'Kingdom', Rank2 = "Phylum", Rank3 = "Class", 
                                     Rank4 = "Order", Rank5 = "Family", 
                                     Rank6 = "Genus", Rank7 = "Species")


colnames(tax_table(GHL_merge)) <- c(Rank1 = 'Kingdom', Rank2 = "Phylum", Rank3 = "Class", 
                                     Rank4 = "Order", Rank5 = "Family", 
                                     Rank6 = "Genus", Rank7 = "Species")

colnames(tax_table(GHR_merge)) <- c(Rank1 = 'Kingdom', Rank2 = "Phylum", Rank3 = "Class", 
                                     Rank4 = "Order", Rank5 = "Family", 
                                     Rank6 = "Genus", Rank7 = "Species")

#See Rank names
#rank_names(AA_rare_top10)

#Rarefy to 9522 reads.  rngseed = TRUE --> is set to default of 711. Note: this removed 15145 OTUs because they weren't included in random subsampling
RF_rare <- rarefy_even_depth(RF_merge, sample.size = 9522, 
                             rngseed = TRUE, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)

#Note:19776OTUs were removed because they are no longer present in any sample after random subsampling
LF_rare <- rarefy_even_depth(LF_merge, sample.size = 7067, 
                             rngseed = TRUE, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)

#16047 OTUs were removed because they are no longer present in any sample after random subsampling
GHFR_rare <- rarefy_even_depth(GHFR_merge, sample.size = 11124, 
                             rngseed = TRUE, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)

#19114 OTUs were removed because they are no longer present in any sample after random subsampling
GHFL_rare <- rarefy_even_depth(GHFL_merge, sample.size = 7866, 
                               rngseed = TRUE, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
#20313 OTUs were removed because they are no longer present in any sample after random subsampling
GHR_rare <- rarefy_even_depth(GHR_merge, sample.size = 23688, 
                               rngseed = TRUE, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)

#20446 OTUs were removed because they are no longer present in any sample after random subsampling
GHL_rare <- rarefy_even_depth(GHL_merge, sample.size = 17636, 
                               rngseed = TRUE, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)


#Keep only the 10 most abundant phyla
#Roots
phylum.sum = tapply(taxa_sums(RF_rare), tax_table(RF_rare)[, "Phylum"], sum, na.rm = TRUE)
top10phyla = names(sort(phylum.sum, TRUE))[1:10]
RF_rare_top10 = prune_taxa((tax_table(RF_rare)[, "Phylum"] %in% top10phyla), RF_rare)

#Leaves
phylum.sum.leaves = tapply(taxa_sums(LF_rare), tax_table(LF_rare)[, "Phylum"], sum, na.rm = TRUE)
top10phylaleaves = names(sort(phylum.sum.leaves, TRUE))[1:10]
LF_rare_top10 = prune_taxa((tax_table(LF_rare)[, "Phylum"] %in% top10phyla), LF_rare)

#GHFR
phylum.sum.GHFR = tapply(taxa_sums(GHFR_rare), tax_table(GHFR_rare)[, "Phylum"], sum, na.rm = TRUE)
top10phylaGHFR = names(sort(phylum.sum.GHFR, TRUE))[1:10]
GHFR_rare_top10 = prune_taxa((tax_table(GHFR_rare)[, "Phylum"] %in% top10phyla), GHFR_rare)

#GHFL
phylum.sum.GHFL = tapply(taxa_sums(GHFL_rare), tax_table(GHFL_rare)[, "Phylum"], sum, na.rm = TRUE)
top10phylaGHFL = names(sort(phylum.sum.GHFL, TRUE))[1:10]
GHFL_rare_top10 = prune_taxa((tax_table(GHFL_rare)[, "Phylum"] %in% top10phyla), GHFL_rare)

#GHL
phylum.sum.GHL = tapply(taxa_sums(GHL_rare), tax_table(GHL_rare)[, "Phylum"], sum, na.rm = TRUE)
top10phylaGHL = names(sort(phylum.sum.GHL, TRUE))[1:10]
GHL_rare_top10 = prune_taxa((tax_table(GHL_rare)[, "Phylum"] %in% top10phyla), GHL_rare)

#GHR
phylum.sum.GHR = tapply(taxa_sums(GHR_rare), tax_table(GHR_rare)[, "Phylum"], sum, na.rm = TRUE)
top10phylaGHR = names(sort(phylum.sum.GHR, TRUE))[1:10]
GHR_rare_top10 = prune_taxa((tax_table(GHR_rare)[, "Phylum"] %in% top10phyla), GHR_rare)


#See variable names
sample_variables(RF_rare)

sample_variables(LF_rare)

sample_variables(GHFL_rare)

sample_variables(GHFR_rare)

sample_variables(GHL_rare)

sample_variables(GHR_rare)

#Collapse top10 phyla Samples by Site and Transform to percentages of total

CollapsedbySite<-merge_samples(RF_rare, "Site", fun=mean)
CollapsedbySite_percent = transform_sample_counts(CollapsedbySite, function(x) 100 * x/sum(x))

CollapsedbySite_top10<-merge_samples(RF_rare_top10, "Site", fun=mean)
CollapsedbySite_top10_percent = transform_sample_counts(CollapsedbySite_top10, function(x) 100 * x/sum(x))

CollapsedbySite_top10_leaves<-merge_samples(LF_rare_top10, "Site", fun=mean)
CollapsedbySite_top10_leaves_percent = transform_sample_counts(CollapsedbySite_top10_leaves, function(x) 100 * x/sum(x))

CollapsedbyExpType_top10_GHFL<-merge_samples(GHFL_rare_top10, "ExpType", fun=mean)
CollapsedbyExpType_top10_GHFL_percent = transform_sample_counts(CollapsedbyExpType_top10_GHFL, function(x) 100 * x/sum(x))

CollapsedbyExpType_top10_GHFR<-merge_samples(GHFR_rare_top10, "ExpType", fun=mean)
CollapsedbyExpType_top10_GHFR_percent = transform_sample_counts(CollapsedbyExpType_top10_GHFR, function(x) 100 * x/sum(x))

CollapsedbySite_top10_GHL<-merge_samples(GHL_rare_top10, "Site", fun=mean)
CollapsedbySite_top10_GHL_percent = transform_sample_counts(CollapsedbySite_top10_GHL, function(x) 100 * x/sum(x))

CollapsedbySite_top10_GHR<-merge_samples(GHR_rare_top10, "Site", fun=mean)
CollapsedbySite_top10_GHR_percent = transform_sample_counts(CollapsedbySite_top10_GHR, function(x) 100 * x/sum(x))


#Merge OTUs by phyla
CollapsedbySite_percent_glom<-tax_glom(CollapsedbySite_percent, taxrank = "Phylum")

CollapsedbySite_percent_glom_top10<-tax_glom(CollapsedbySite_top10_percent, taxrank = "Phylum")

CollapsedbySite_percent_glom_top10_leaves<-tax_glom(CollapsedbySite_top10_leaves_percent, taxrank = "Phylum")

CollapsedbyExpType_percent_glom_top10_GHFL<-tax_glom(CollapsedbyExpType_top10_GHFL_percent, taxrank = "Phylum")

CollapsedbyExpType_percent_glom_top10_GHFR<-tax_glom(CollapsedbyExpType_top10_GHFR_percent, taxrank = "Phylum")

CollapsedbySite_percent_glom_top10_GHL<-tax_glom(CollapsedbySite_top10_GHL_percent, taxrank = "Phylum")

CollapsedbySite_percent_glom_top10_GHR<-tax_glom(CollapsedbySite_top10_GHR_percent, taxrank = "Phylum")


# Repair the merged values associated with each site after merge (see: http://joey711.github.io/phyloseq-demo/Restroom-Biogeography)
sample_data(CollapsedbySite)$Site <- levels(sample_data(CollapsedbySite)$Site)

sample_data(CollapsedbySite_top10)$Site <- levels(sample_data(CollapsedbySite_top10)$Site)

sample_data(CollapsedbySite_top10_leaves)$Site <- levels(sample_data(CollapsedbySite_top10_leaves)$Site)

sample_data(CollapsedbyExpType_top10_GHFL)$Site <- levels(sample_data(CollapsedbyExpType_top10_GHFL)$ExpType)

sample_data(CollapsedbyExpType_top10_GHFR)$Site <- levels(sample_data(CollapsedbyExpType_top10_GHFR)$ExpType)

sample_data(CollapsedbySite_top10_GHL)$Site <- levels(sample_data(CollapsedbySite_top10_GHL)$Site)

sample_data(CollapsedbySite_top10_GHR)$Site <- levels(sample_data(CollapsedbySite_top10_GHR)$Site)

#bar chart relative abundance



#PLot top 10 Phyla
RF_10_plot<-plot_bar(CollapsedbySite_percent_glom_top10, fill = "Phylum")

RF_10_plot+theme_classic()+labs(x="Site", y="Relative Abundance")+ggtitle("Relative abundance within field roots of 10 most abundant phyla")


LF_10_plot<-plot_bar(CollapsedbySite_percent_glom_top10_leaves, fill = "Phylum")

LF_10_plot+theme_classic()+labs(x="Site", y="Relative Abundance")+ggtitle("Relative Abundance within Field Leaves of 10 Most Abundant Phyla")


GHFL_10_plot<-plot_bar(CollapsedbyExpType_percent_glom_top10_GHFL, fill = "Phylum")

GHFL_10_plot+theme_classic()+labs(x="Experimental Location", y="Relative Abundance")+ggtitle("Relative Abundance within Greenhouse and Field Leaves of 10 Most Abundant Phyla")


GHFR_10_plot<-plot_bar(CollapsedbyExpType_percent_glom_top10_GHFR, fill = "Phylum")

GHFR_10_plot+theme_classic()+labs(x="Experimental Location", y="Relative Abundance")+ggtitle("Relative Abundance  within Greenhouse and Field Roots of 10 Most Abundant Phyla")


GHL_10_plot<-plot_bar(CollapsedbySite_percent_glom_top10_GHL, fill = "Phylum")

GHL_10_plot+theme_classic()+labs(x="Soil Origin", y="Relative Abundance")+ggtitle("Relative Abundance within Greenhouse Leaves of 10 Most Abundant Phyla")


GHR_10_plot<-plot_bar(CollapsedbySite_percent_glom_top10_GHR, fill = "Phylum")

GHR_10_plot+theme_classic()+labs(x="Soil Origin", y="Relative Abundance")+ggtitle("Relative Abundance within Greenhouse roots of 10 Most Abundant Phyla")



p<-plot_bar(CollapsedbySite_percent_glom, fill = "Phylum")
+scale_fill_manual(values= phylumPalette)
p
p+theme_classic()+labs(x="Site", y="Relative Abundance")+ggtitle("Relative Abundance\n within 10 Most Abundant Phyla")
