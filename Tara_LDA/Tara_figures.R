# See lines 69-77 to specify paths to figure, code and data folders.

library(grid)
library(gridExtra)
library(ggplot2)
library(ggforce)
library(ape)
library(car)
library(kriging)
library(devtools)
library(sirt)
# To plot the kriged maps:
library(scales)
library(gstat)
library(colorspace)
# To plot shapes using readOGR (not available on cluster):
library(rgdal)
library(sp)
library(rgeos)
library(maptools)
# To use raster manipulating functions:
library(raster)
library(ggtree)
library(marmap)
# To use as.grob:
library(ggplotify)
# To use varpart and rda:
library(vegan)
# For linear mixed models:
library(nlme)
# For sourcing from url:
library(devtools)

noArcticNoBiomark = 0
noCoastalArcticNoBiomark = 0
ArcticOnly = 0
noLagoon = 1

if (noLagoon)
{
  noLagoon_insert = "_noLagoon"
  lagoon_stations = c("TARA_043","TARA_046","TARA_049","TARA_054","TARA_057","TARA_113","TARA_114","TARA_115","TARA_116","TARA_117","TARA_118","TARA_119","TARA_120","TARA_121")
} else
  noLagoon_insert = ""

if (!noArcticNoBiomark && !noCoastalArcticNoBiomark && !ArcticOnly)
{
  noArcticNoBiomark_insert = ""
  Arctic_insert = "/All_stations"
} else if (ArcticOnly)
{
  noArcticNoBiomark_insert = "_ArcticOnly"
  Arctic_insert = "/Arctic"
} else if (noArcticNoBiomark)
{
  noArcticNoBiomark_insert = "_noArcticNoBiomark"
  Arctic_insert = "/All_stations_but_Arctic"
} else if (noCoastalArcticNoBiomark)
{
  noArcticNoBiomark_insert = "_noCoastalArcticNoBiomark"
  Arctic_insert = "/Arctic2stations"
} else if (NorthAtlantic)
{
  noArcticNoBiomark_insert = "_NorthAtlantic"
  Arctic_insert = "/NorthAtlantic"
}

########################
# Path to figure folder:
figure_folder = paste0("/Users/guilhemsommeria-klein/Desktop/Manuscrits/Tara/Figures/Gibbs",Arctic_insert)
# Path to saved results:
code_folder = "/Users/guilhemsommeria-klein/Desktop/Code/Projects/eDNA_LDA"
results_folder = paste0(code_folder,"/Tara_LDA/Saved_results")
# Path to local data folder:
data_folder = "/Users/guilhemsommeria-klein/Desktop/Post-doc_ENS/Donnees_Tara"
# Path to cluster data folder - workspace 3:
data_folder_workspace3 = "/Users/guilhemsommeria-klein/Desktop/Serveur_Bioclust.kingdoms.workspace3/Donnees_Tara"
########################

groups_to_remove = c("undetermined_eukaryote","other_filosan_*","other_myzozoan_*","undetermined_alveolate","environmental_samples",
                       "Unknown","other_eukaryotes_*","other_prasinophytes","other_core_chlorophytes","other_endomyxan_*",
                       "other_MAST_*","undetermined_archaeplastida",
                       "other_Lobosa","undetermined_incertae_sedis","other_heterotrophic_stramenopiles_*","undetermined_stramenopiles","other_lower_Fungi_*",
                       "other_Ochrophyta_*","undetermined_MAST","undetermined_opisthokonts","undetermined_filosan","undetermined_holozoan",
                       "undetermined_ochrophyta","undetermined_thecofilosan","undetermined_fungi","other_Discoba")
# groups_to_remove_psbO = c("Trebouxiophyceae","Chrysophyceae","Ciliophora","Euglenida")

taxo_groups = readRDS(paste0(results_folder,"/taxo_groups_modif_2plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
taxo_groups_psbO = readRDS(paste0(results_folder,"/psbO_taxo_groups_modif_2plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
taxo_groups_V4 = readRDS(paste0(results_folder,"/V4_taxo_groups_modif_2plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
####
diversity = readRDS(paste0(results_folder,"/diversity_2plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
# names(diversity) = taxo_groups # Recently modified! Check for compatiblity.
diversity_psbO = readRDS(paste0(results_folder,"/psbO_diversity_2plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
diversity_V4 = readRDS(paste0(results_folder,"/V4_diversity_2plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
names(diversity_V4) = taxo_groups_V4 # Recently modified! Check for compatiblity.
####
relativeAbund = readRDS(paste0(results_folder,"/groups_relativeAbund",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
####
load(paste0(results_folder,"/group_sizes_byStationByDepth_2plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"0.Rdata"))
div_threshold = 100
######
selected_groups = !(taxo_groups %in% groups_to_remove) & diversity > div_threshold
selected_groups_psbO = selected_groups_psbO.5nd = selected_groups_psbO.5 = selected_groups_psbO.4nd = selected_groups_psbO.4 = rep(F,length(taxo_groups_psbO))
selected_groups_psbO.5[1:5] = T
selected_groups_psbO.5nd[c(1,3:6)] = T
selected_groups_psbO.4[1:4] = T
selected_groups_psbO.4nd[c(1,3:5)] = T
selected_groups_psbO[1:6] = T
selected_groups_V4 = rep(F,length(taxo_groups_V4))
selected_groups_V4[1:41] = T
######
tot_reads = readRDS(paste0(results_folder,"/Total.read.numbers.rds"))

##########
dominant_function = readRDS(paste0(results_folder,"/dominant_function_2plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
dominant_function0 = dominant_function[selected_groups & diversity>div_threshold]
# dominant_function0[dominant_function0 %in% c("other metazoa","gelatineous_carnivores_filterers","copepoda","pteropoda")] = "metazoa"
# dominant_function0[dominant_function0 %in% c("unknown","photohost")] = NA
dominant_function0[dominant_function0 == "unknown"] = NA
dominant_function0[dominant_function0 == "photohost"] = "Collodaria"
dominant_function0[dominant_function0 == "pteropoda"] = "Pteropoda"
dominant_function0[dominant_function0 == "copepoda"] = "Copepoda"
dominant_function0[dominant_function0 == "gelatineous_carnivores_filterers"] = "Gel. carn. filterers"
dominant_function0[taxo_groups[selected_groups & diversity>div_threshold] == "Dinophyceae"] = "Dinophyceae"
dominant_function0[taxo_groups[selected_groups & diversity>div_threshold] == "Bacillariophyta"] = "Bacillariophyta"
dominant_function0[taxo_groups[selected_groups & diversity>div_threshold] == "Chrysophyceae"] = "Other phototrophs"
dominant_function0[dominant_function0 == "phototroph"] = "Other phototrophs"
dominant_function0[dominant_function0 == "phagotroph"] = "Phagotrophs"
dominant_function0[dominant_function0 == "other metazoa"] = "Other metazoa"
dominant_function0[dominant_function0 == "parasite"] = "Parasites"

dominant_function1 = dominant_function0[!dominant_function0 %in% c("Dinophyceae","Collodaria")]
dominant_function1[dominant_function1 %in% c("Bacillariophyta","Other phototrophs")] = "Phototrophs"
dominant_function1[dominant_function1 %in% c("Gel. carn. filterers","Other metazoa","Pteropoda","Copepoda")] = "Metazoans"
# alpha = rep(1,length(taxo_groups[selected_groups]))
# alpha[dominant_function0 %in% c("copepoda","pteropoda")] = 0 
# Changes how the factors are stored (order of levels()) so that geom_boxplot plots them by deceasing number of groups:
# dominant_function0 = factor(dominant_function0,names(sort(table(as.factor(dominant_function0)),decreasing = T)))
# Changes how the factors are stored (order of levels()) so that geom_boxplot plots them in the specified order:
if (div_threshold == 100)
{
  dominant_function0 = factor(dominant_function0,c("Bacillariophyta","Other phototrophs","Collodaria","Phagotrophs","Dinophyceae","Gel. carn. filterers","Pteropoda","Copepoda","Other metazoa","Parasites"))
  point_groups = c("Dinophyceae","Copepoda","Collodaria","Bacillariophyta","Pteropoda")
  
  dominant_function1 = factor(dominant_function1,c("Phagotrophs","Metazoans","Parasites","Phototrophs"))
} else if (div_threshold == 1000)
{
  dominant_function0 = factor(dominant_function0,c("Bacillariophyta","Other phototrophs","Collodaria","Phagotrophs","Dinophyceae","Gel. carn. filterers","Copepoda","Parasites"))
  point_groups = c("Dinophyceae","Copepoda","Collodaria","Bacillariophyta")
} 
#######
data.folder_name = paste0(data_folder,"/18S_V9_TARA_CompleteSizeRange_byStationByDepth_AllTAxa",noArcticNoBiomark_insert,noLagoon_insert)
load(paste0(data.folder_name,"/taxo_ref.Rdata"))
taxo_groups_OTUs = taxo_ref$taxogroup2
taxo_groups_unmodified = levels(as.factor(taxo_groups_OTUs))[sort.int(table(as.factor(taxo_groups_OTUs)),index.return = T,decreasing = T)$ix]
#######

#########
optimalK_prevalence.min.crossValid.complete = readRDS(paste0(results_folder,"/optimalK_prevalence.min.crossValid_Gibbs10sampleFolds2-35t_iter1000thin25burnin500_2plusOTUs_noLagoon.rds"))
optimalK_prevalence.min.crossValid.allTaxa = optimalK_prevalence.min.crossValid.complete[[1]]  
optimalK_prevalence.min.crossValid = optimalK_prevalence.min.crossValid.complete[[2]]   
####
optimalK_prevalence.min.crossValid_psbO = readRDS(paste0(results_folder,"/psbO_optimalK_prevalence.min.crossValid_Gibbs10sampleFolds2-30t_iter1000thin25burnin2000_2plusOTUs_noLagoon.rds"))
####
optimalK_prevalence.min.crossValid_V4 = readRDS(paste0(results_folder,"/V4_optimalK_prevalence.min.crossValid_Gibbs10sampleFolds2-30t_iter1000thin25burnin2000_2plusOTUs_noLagoon.rds"))
#######

######
mean_sim = readRDS(paste0(results_folder,"/mean_sim_optimalK_Gibbs.prevalence.min.crossValid10sampleFolds_2plusOTUs_noLagoon.rds"))
mean_sim.allTaxa = 0.894643970320553 
normalized.VI.within.group = readRDS(paste0(results_folder,"/Normalized_VI_within.groups_100reals_different.stations.means.rds"))
mean.normalized.VI.within.group = readRDS(paste0(results_folder,"/Normalized_VI_mean.within.groups_100reals_different.stations.means.rds"))
mean.normalized.VI.allTaxa = 0.8196249
R.hat.500.result = readRDS(paste0(results_folder,"/Rhat.burnin.500_nb_iter3000_nb_real100_occurrence.rds"))
R.hat.500.allTaxa = R.hat.500.result[[1]]
R.hat.500 = R.hat.500.result[[2]]
R.hat.1000.result = readRDS(paste0(results_folder,"/Rhat.burnin.1000_nb_iter3000_nb_real100_occurrence.rds"))
R.hat.1000.allTaxa = R.hat.1000.result[[1]]
R.hat.1000 = R.hat.1000.result[[2]]
R.hat.2000.result = readRDS(paste0(results_folder,"/Rhat.burnin.2000_nb_iter3000_nb_real100_occurrence.rds"))
R.hat.2000.allTaxa = R.hat.2000.result[[1]]
R.hat.2000 = R.hat.2000.result[[2]]
# SUR.DCM_VI_over_K = readRDS(results_folder,"/VI.over.K_Gibbs.prevalence.min.crossValid10sampleFolds.rds")
SUR.DCM_Normalized.VI = readRDS(paste0(results_folder,"/Normalized.VI_SUR.DCM_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
######

#######
Moran.I = readRDS(paste0(results_folder,"/Moran.I_inverse.squares_weighted.mean.over.topics_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
I_square.observed_w.mean_allTaxa = Moran.I[[1]]
I_square.observed_w.mean = Moran.I[[2]]
I_square.p.value_w.mean_allTaxa = Moran.I[[3]]
I_square.p.value_w.mean = Moran.I[[4]]
autocorr.scale = readRDS(paste0(results_folder,"/MoranI_charach.autocorr.scale_20increment_2plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
charac_scale = autocorr.scale[[1]]
charac_scale.allTaxa = autocorr.scale[[2]]
slope.intercept = readRDS(paste0(results_folder,"/MoranI_slope.intercept_20increment_2plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
slope = slope.intercept[[1]]
intercept = slope.intercept[[2]]
basin_I.result = readRDS(paste0(results_folder,"/Basin_I_noLagoon.rds"))
basin_I = basin_I.result[[1]]
lat_I.result = readRDS(paste0(results_folder,"/Lat_I_sigma2.25",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
lat_I = lat_I.result[[1]]
#####
Moran.I.104 = readRDS(paste0(results_folder,"/104_Moran.I_inverse.squares_weighted.mean.over.topics_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
I_square.observed_w.mean.104 = Moran.I.104[[1]]
I_square.p.value_w.mean.104 = Moran.I.104[[2]]
diversity.104 = readRDS(paste0(results_folder,"/Rarefied_diversity.104.105.106.rds"))[[1]]
Moran.I.105 = readRDS(paste0(results_folder,"/105_Moran.I_inverse.squares_weighted.mean.over.topics_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
I_square.observed_w.mean.105 = Moran.I.105[[1]]
I_square.p.value_w.mean.105 = Moran.I.105[[2]]
diversity.105 = readRDS(paste0(results_folder,"/Rarefied_diversity.104.105.106.rds"))[[2]]
######
mean_prevalence_results = readRDS(paste0(results_folder,"/mean_OTU.prevalence.200.250.300.350.500.1000.rds"))
mean_prevalence = mean_prevalence_results[[1]]
mean_prevalence.200 = mean_prevalence_results[[2]]
mean_prevalence.250 = mean_prevalence_results[[3]]
mean_prevalence.300 = mean_prevalence_results[[4]]
mean_prevalence.350 = mean_prevalence_results[[5]]
mean_prevalence.500 = mean_prevalence_results[[6]]
mean_prevalence.1000 = mean_prevalence_results[[7]]
median_prevalence_results = readRDS(paste0(results_folder,"/median_OTU.prevalence.200.250.300.350.500.1000.rds"))
median_prevalence = median_prevalence_results[[1]]
median_prevalence.200 = median_prevalence_results[[2]]
median_prevalence.250 = median_prevalence_results[[3]]
median_prevalence.300 = median_prevalence_results[[4]]
median_prevalence.350 = median_prevalence_results[[5]]
median_prevalence.500 = median_prevalence_results[[6]]
median_prevalence.1000 = median_prevalence_results[[7]]
Moran.I.200 = readRDS(paste0(results_folder,"/200_Moran.I_inverse.squares_weighted.mean.over.topics_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
I_square.observed_w.mean.200 = Moran.I.200[[1]]
I_square.p.value_w.mean.200 = Moran.I.200[[2]]
charac_scale.200 = readRDS(paste0(results_folder,"/200_MoranI_charach.autocorr.scale_20increment_2plusOTUs_noLagoon.rds"))
diversity.200 = readRDS(paste0(results_folder,"/Total.read.numbers.200_diversity.200.rds"))[[2]]
tot_reads.200 = readRDS(paste0(results_folder,"/Total.read.numbers.200_diversity.200.rds"))[[1]]
Moran.I.250 = readRDS(paste0(results_folder,"/250_Moran.I_inverse.squares_weighted.mean.over.topics_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
I_square.observed_w.mean.250 = Moran.I.250[[1]]
I_square.p.value_w.mean.250 = Moran.I.250[[2]]
diversity.250 = readRDS(paste0(results_folder,"/Total.read.numbers.250_diversity.250.rds"))[[2]]
tot_reads.250 = readRDS(paste0(results_folder,"/Total.read.numbers.250_diversity.250.rds"))[[1]]
Moran.I.300 = readRDS(paste0(results_folder,"/300_Moran.I_inverse.squares_weighted.mean.over.topics_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
I_square.observed_w.mean.300 = Moran.I.300[[1]]
I_square.p.value_w.mean.300 = Moran.I.300[[2]]
diversity.300 = readRDS(paste0(results_folder,"/Total.read.numbers.300_diversity.300.rds"))[[2]]
tot_reads.300 = readRDS(paste0(results_folder,"/Total.read.numbers.300_diversity.300.rds"))[[1]]
Moran.I.350 = readRDS(paste0(results_folder,"/350_Moran.I_inverse.squares_weighted.mean.over.topics_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
I_square.observed_w.mean.350 = Moran.I.350[[1]]
I_square.p.value_w.mean.350 = Moran.I.350[[2]]
diversity.350 = readRDS(paste0(results_folder,"/Total.read.numbers.350_diversity.350.rds"))[[2]]
tot_reads.350 = readRDS(paste0(results_folder,"/Total.read.numbers.350_diversity.350.rds"))[[1]]
Moran.I.500 = readRDS(paste0(results_folder,"/500_Moran.I_inverse.squares_weighted.mean.over.topics_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
I_square.observed_w.mean.500 = Moran.I.500[[1]]
I_square.p.value_w.mean.500 = Moran.I.500[[2]]
diversity.500 = readRDS(paste0(results_folder,"/Total.read.numbers.500_diversity.500.rds"))[[2]]
tot_reads.500 = readRDS(paste0(results_folder,"/Total.read.numbers.500_diversity.500.rds"))[[1]]
#######
# Moran.I.200.random = readRDS(results_folder,"/200.random.1-5_Moran.I_inverse.squares_weighted.mean.over.topics_Gibbs.prevalence.min.crossValid10sampleFolds.rds")
# I_square.observed_w.mean.200.random.5reals = Moran.I.200.random[[1]]
# I_square.p.value_w.mean.200.random.5reals = Moran.I.200.random[[2]]
Moran.I.250.random = readRDS(paste0(results_folder,"/250.random_Moran.I_inverse.squares_weighted.mean.over.topics_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
I_square.observed_w.mean.250.random = Moran.I.250.random[[1]]
I_square.p.value_w.mean.250.random = Moran.I.250.random[[2]]
Moran.I.350.random = readRDS(paste0(results_folder,"/350.random_Moran.I_inverse.squares_weighted.mean.over.topics_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
I_square.observed_w.mean.350.random = Moran.I.350.random[[1]]
I_square.p.value_w.mean.350.random = Moran.I.350.random[[2]]
Moran.I.200.random = readRDS(paste0(results_folder,"/200.random.1-10_Moran.I_inverse.squares_weighted.mean.over.topics_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
I_square.observed_w.mean.200.random.10reals = Moran.I.200.random[[1]]
I_square.p.value_w.mean.200.random.10reals = Moran.I.200.random[[2]]
Moran.I.400.random = readRDS(paste0(results_folder,"/400.random.1-10_Moran.I_inverse.squares_weighted.mean.over.topics_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
I_square.observed_w.mean.400.random.10reals = Moran.I.400.random[[1]]
I_square.p.value_w.mean.400.random.10reals = Moran.I.400.random[[2]]
Moran.I.800.random = readRDS(paste0(results_folder,"/800.random.1-10_Moran.I_inverse.squares_weighted.mean.over.topics_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
I_square.observed_w.mean.800.random.10reals = Moran.I.800.random[[1]]
I_square.p.value_w.mean.800.random.10reals = Moran.I.800.random[[2]]
mean_prevalence_results = readRDS(paste0(results_folder,"/mean_OTU.prevalence.250.350.random.rds"))
mean_prevalence.250.random = mean_prevalence_results[[1]]
mean_prevalence.350.random = mean_prevalence_results[[1]]
median_prevalence_results = readRDS(paste0(results_folder,"/median_OTU.prevalence.250.350.random.rds"))
median_prevalence.250.random = median_prevalence_results[[1]]
median_prevalence.350.random = median_prevalence_results[[1]]
#####
Moran.I.random.groups = readRDS(paste0(results_folder,"/random.group.div_Moran.I_inverse.squares_weighted.mean.over.topics_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
I_square.observed_w.mean.random.groups = Moran.I.random.groups[[1]]
I_square.p.value_w.mean.random.groups = Moran.I.random.groups[[2]]
lat_I_random.groups = readRDS(paste0(results_folder,"/random.group.div_lat_I_sigma2.25_noLagoon.rds"))[[1]]
basin_I_random.groups = readRDS(paste0(results_folder,"/random.group.div_basin_I_noLagoon.rds"))[[1]]
#####
shannon_results = readRDS(paste0(results_folder,"/Shannon.per.station_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
shannon_SUR = unlist(lapply(shannon_results[[2]],mean))
shannon_DCM = unlist(lapply(shannon_results[[3]],mean))
simpson_results = readRDS(paste0(results_folder,"/Simpson.per.station_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
simpson_SUR = 1-unlist(lapply(simpson_results[[2]],mean))
invsimpson_SUR = 1/unlist(lapply(simpson_results[[2]],mean))
shannon_total_results = readRDS(paste0(results_folder,"/Shannon.per.group_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
shannon_total_SUR = shannon_total_results[[2]]
shannon_total_DCM = shannon_total_results[[3]]
simpson_total_results = readRDS(paste0(results_folder,"/Simpson.per.group_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
simpson_total_SUR = 1-simpson_total_results[[2]]
invsimpson_total_SUR = 1/simpson_total_results[[2]]
nb_dominants_results = readRDS(paste0(results_folder,"/Nb.dominant.assemblages.per.group_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
nb_dominants = nb_dominants_results[[1]]
nb_absolute_dominants = nb_dominants_results[[2]]
######
Moran.I_V9.psbO = readRDS(paste0(results_folder,"/Moran.I_psbO.stations_inverse.squares_weighted.mean.over.topics_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
I_square.observed_w.mean_V9.psbO = Moran.I_V9.psbO[[1]]
I_square.p.value_w.mean_V9.psbO = Moran.I_V9.psbO[[2]]
charac_scale_V9.psbO = readRDS(paste0(results_folder,"/MoranI_psbO.stations_charach.autocorr.scale_20increment_2plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
lat_I_V9.psbO = readRDS(paste0(results_folder,"/Lat_I_psbO.stations_sigma2.25",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
basin_I_V9.psbO = readRDS(paste0(results_folder,"/Basin_I_psbO.stations",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
######
Moran.I_V9.V4 = readRDS(paste0(results_folder,"/Moran.I_V4.stations_inverse.squares_weighted.mean.over.topics_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
I_square.observed_w.mean_V9.V4 = Moran.I_V9.V4[[1]]
I_square.p.value_w.mean_V9.V4 = Moran.I_V9.V4[[2]]
charac_scale_V9.V4 = readRDS(paste0(results_folder,"/MoranI_V4.stations_charach.autocorr.scale_20increment_2plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
lat_I_V9.V4 = readRDS(paste0(results_folder,"/Lat_I_V4.stations_sigma2.25",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
basin_I_V9.V4 = readRDS(paste0(results_folder,"/Basin_I_V4.stations",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
######
Moran.I_psbO = readRDS(paste0(results_folder,"/psbO_Moran.I_inverse.squares_weighted.mean.over.topics_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
I_square.observed_w.mean_psbO = Moran.I_psbO[[1]]
I_square.p.value_w.mean_psbO = Moran.I_psbO[[2]]
charac_scale_psbO = readRDS(paste0(results_folder,"/psbO_MoranI_charach.autocorr.scale_20increment_2plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
lat_I_psbO = readRDS(paste0(results_folder,"/psbO_lat_I_sigma2.25",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
basin_I_psbO = readRDS(paste0(results_folder,"/psbO_basin_I",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
####
Moran.I_V4 = readRDS(paste0(results_folder,"/V4_Moran.I_inverse.squares_weighted.mean.over.topics_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
I_square.observed_w.mean_V4 = Moran.I_V4[[1]]
I_square.p.value_w.mean_V4 = Moran.I_V4[[2]]
charac_scale_V4 = readRDS(paste0(results_folder,"/V4_MoranI_charach.autocorr.scale_20increment_2plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
lat_I_V4 = readRDS(paste0(results_folder,"/V4_lat_I_sigma2.25",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
basin_I_V4 = readRDS(paste0(results_folder,"/V4_basin_I",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
######
Normalized_VI_V9.V4.psbO = readRDS(paste0(results_folder,"/Normalized.VI_marker.comparison.SUR_V9-V4_V9-psbO_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
Normalized_VI_V9.V4 = Normalized_VI_V9.V4.psbO[[1]]
Normalized_VI_V9.psbO = Normalized_VI_V9.V4.psbO[[2]]
Normalized_VI_V4.psbO = readRDS(paste0(results_folder,"/Normalized.VI_group.comparison.SUR_V4.psbO_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
Normalized_VI_V4 = Normalized_VI_V4.psbO[[1]]
Normalized_VI_psbO = Normalized_VI_V4.psbO[[2]]

# Old:
# varpart.indSelec = readRDS(paste0(results_folder,"/varpart_lda_environment-currents_bothDirectionsIndependentSelection_PCA0",abiotic_pca_insert,biotic_pca_insert,stdzation_insert,noArcticNoBiomark_insert,noLagoon_insert,".rds"))
#
# Variable selection + individually selected axes:
# varpart.indSelec = readRDS(paste0(results_folder,"/varpart_lda_Gibbs.prevalence.min.crossValid10sampleFolds",
#                                   "_separate.SUR.DCM_both.directions.independent.selection_PCA0_abioticPCA_bioticPCA_eigenvalueThres0.8",noLagoon_insert,"_(0ter).rds"))
# Variable selection without preliminary tests:
# varpart.indSelec = readRDS(paste0(results_folder,"/varpart_lda_Gibbs.prevalence.min.crossValid10sampleFolds",
#                                   "_separate.SUR.DCM_both.directions.independent.selection_PCA0_no.indiv.signif.axes_no.prelim.global.test_abioticPCA_bioticPCA_eigenvalueThres0.8.rds"))
# Only variable selection:
# varpart.indSelec = readRDS(paste0(results_folder,"/varpart_lda_Gibbs.prevalence.min.crossValid10sampleFolds",
                                  # "_separate.SUR.DCM_both.directions.independent.selection_PCA0_no.indiv.signif.axes_abioticPCA_bioticPCA_eigenvalueThres0.8.rds"))
# Only individually selected axes:
# varpart.indSelec = readRDS(paste0(results_folder,"/varpart_lda_Gibbs.prevalence.min.crossValid10sampleFolds",
#                                   "_separate.SUR.DCM_PCA0_abioticPCA_bioticPCA_eigenvalueThres0.8.rds"))
varpart.indSelec = readRDS(paste0(results_folder,"/varpart_lda_Gibbs.prevalence.min.crossValid10sampleFolds",
                                  "_separate.SUR.DCM_PCA0_BH.correction_abioticPCA_bioticPCA_eigenvalueThres0.8.rds"))
# No selection whatsoever:
# varpart.indSelec = readRDS(paste0(results_folder,"/varpart_lda_Gibbs.prevalence.min.crossValid10sampleFolds",
#                                   "_separate.SUR.DCM_PCA0_no.indiv.signif.axes_abioticPCA_bioticPCA_eigenvalueThres0.8.rds"))
varpart.env.spatial = varpart.indSelec[[1]]
varpart.env.spatial[[1]][varpart.env.spatial[[1]]<0] = 0
varpart.env.spatial[[2]][varpart.env.spatial[[2]]<0] = 0
varpart.env.spatial.pval = varpart.indSelec[[2]]

sub.varpart.biotic.abiotic = varpart.indSelec[[3]]
sub.varpart.biotic.abiotic[[1]][sub.varpart.biotic.abiotic[[1]]<0] = 0
sub.varpart.biotic.abiotic[[2]][sub.varpart.biotic.abiotic[[2]]<0] = 0
sub.varpart.biotic.abiotic.pval = varpart.indSelec[[4]]

# sub.varpart.SUR.DCM = varpart.indSelec[[5]]
# sub.varpart.SUR.DCM[sub.varpart.SUR.DCM<0] = 0
# sub.varpart.SUR.DCM.pval = varpart.indSelec[[6]]

varpart.biotic.abiotic = varpart.indSelec[[5]]
varpart.biotic.abiotic[[1]][varpart.biotic.abiotic[[1]]<0] = 0
varpart.biotic.abiotic[[2]][varpart.biotic.abiotic[[2]]<0] = 0
varpart.biotic.abiotic.pval = varpart.indSelec[[6]]

#########
Normalized_VI = readRDS(paste0(results_folder,"/Normalized.VI_group.comparison_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
div_threshold = 100
NormalizedVI_pcoa = vector(length = 3, mode = "list")
names(NormalizedVI_pcoa) = c("SUR","DCM","All")
for (i_case in 1:3)
{
  VI.pcoa = pcoa(as.dist(Normalized_VI[[i_case]][selected_groups & diversity>div_threshold,selected_groups & diversity>div_threshold]))
  NormalizedVI_pcoa[[i_case]] = VI.pcoa$vectors#[-69,]
  # VI_pcoa = pcoa.all(as.dist(Normalized_VI[[1]][selected_groups,selected_groups]))
}

########
normalized.VI.reals = readRDS(paste0(results_folder,"/Normalized_VI_10reals_different.stations.means.rds"))
normalized.VI.best.real = normalized.VI.reals[[1]] 
normalized.VI.10.reals = normalized.VI.reals[[2]]

# Normalized_VI_consistenbasin_I_withint = readRDS(paste0(data_folder,"/Normalized.VI.consistent_group.comparison_Gibbs.prevalence.min.crossValid10sampleFolds.rds"))
# NormalizedVI.consistent_pcoa = vector(length = 3, mode = "list")
# names(NormalizedVI.consistent_pcoa) = c("SUR","DCM","All")
# for (i_case in 3)
# {
#   VI.pcoa = pcoa(as.dist(Normalized_VI_consistent[[i_case]][selected_groups & diversity>div_threshold,selected_groups & diversity>div_threshold]))
#   NormalizedVI.consistent_pcoa[[i_case]] = VI.pcoa$vectors#[-69,]
#   # VI_pcoa = pcoa.all(as.dist(Normalized_VI[[1]][selected_groups,selected_groups]))
# }

devtools::source_url("https://github.com/guilhemSK/Useful_functions/raw/main/Plotting_functions/cor_plot.R")
devtools::source_url("https://github.com/guilhemSK/Useful_functions/raw/main/Plotting_functions/box_plot.R")
devtools::source_url("https://github.com/guilhemSK/Useful_functions/raw/main/Plotting_functions/bar_plot.R")
devtools::source_url("https://github.com/guilhemSK/Useful_functions/raw/main/Plotting_functions/curv_plot.R")

# PCoA Shepard
{
  plot.PCoA.Shepard = ggplot(data = data.frame(x = as.vector(as.dist(Normalized_VI[[3]][selected_groups & diversity>div_threshold,selected_groups & diversity>div_threshold][-69,-69])),
                                                      y = as.vector(dist(NormalizedVI_pcoa[[3]])))) + 
    geom_point(aes(x,y), na.rm=T) +
    theme_bw() +
    # ggtitle(LETTERS[1]) +
    # geom_hline(yintercept = 0, linetype = "dashed") +
    # geom_vline(xintercept = 0, linetype = "dashed") +
    # ggtitle("Norm. V.I. between clades' biogeographies") +
    theme(axis.title=element_text(size=22),
          axis.text=element_text(size=22),
          text = element_text(size=22),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm")) +
    labs(x="Full distances", y="Reduced-space distances")
    # geom_vline(xintercept = prop_within_OTU_null, linetype = "dashed") +
  
  pdf(paste0(figure_folder,"/PCoA_Shepard_allAxes.pdf"))
  grid::grid.draw(plot.PCoA.Shepard)
  # grid::grid.draw(plot.VI.NMDS.biplot[[3]])
  dev.off()
  
  pdf(paste0(figure_folder,"/PCoA_variance.pdf"))
  plot(1:length(VI.pcoa$values$Eigenvalues),VI.pcoa$values$Eigenvalues/sum(VI.pcoa$values$Eigenvalues),ann=F)
  title(x = "", y="Variance")
  dev.off()
  
  ##########
}

# t-SNE:
# library(tsne)
# test = tsne(as.dist(Normalized_VI[[1]][selected_groups & diversity>div_threshold,selected_groups & diversity>div_threshold]),
#             initial_config = NULL, k = 2, initial_dims = 10, perplexity = 10,max_iter = 1000, min_cost = 0, epoch_callback = NULL, whiten = TRUE,epoch=100)

# library(DistatisR)
# dis.matrix = Normalized_VI[[3]][selected_groups,selected_groups]
# dis.matrix[upper.tri(dis.matrix)] = t(dis.matrix)[upper.tri(dis.matrix)]
# diag(dis.matrix) = 0
# test = mmds(dis.matrix)

# names(diversity) = taxo_groups
# taxo_groups_V4[selected_groups_V4][diversity[taxo_groups_V4[selected_groups_V4]] < diversity_V4[selected_groups_V4]]

# NMDS
{
  library(vegan)
  i_case = 3
  VI.NMDS = metaMDS(as.dist(Normalized_VI[[i_case]][selected_groups & diversity>div_threshold,selected_groups & diversity>div_threshold]), 
                      k=3, 
                      maxit = 100000, 
                      sratmax = 0.999999999,
                      sfgrmin = 10^-9)
  
  pdf(paste0(figure_folder,"/NMDS_stress.vs.k.pdf"))
  plot(c(6,5,4,3,2),c(0.106015,0.1216638,0.143212,0.1762817,0.2371242),ann=F)
  title(x = "Nb of axes", y="Stress value")
  dev.off()
  
  plot.NMDS.Shepard = ggplot(data = data.frame(x = as.vector(as.dist(Normalized_VI[[3]][selected_groups & diversity>div_threshold,selected_groups & diversity>div_threshold][-69,-69])),
                                               y = as.vector(dist(VI.NMDS2$points[-69,])))) + 
    geom_point(aes(x,y), na.rm=T) +
    theme_bw() +
    # ggtitle(LETTERS[1]) +
    # geom_hline(yintercept = 0, linetype = "dashed") +
    # geom_vline(xintercept = 0, linetype = "dashed") +
    # ggtitle("Norm. V.I. between clades' biogeographies") +
    theme(axis.title=element_text(size=22),
          axis.text=element_text(size=22),
          text = element_text(size=22),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm"),
          legend.position="bottom",
          legend.text=element_text(size=20),
          legend.title=element_text(size=22)) +
    labs(x="Full distances", y="Reduced-space distances")
  # geom_vline(xintercept = prop_within_OTU_null, linetype = "dashed") +
  
  pdf(paste0(figure_folder,"/NMDS_Shepard_3axes.pdf"))
  grid::grid.draw(plot.NMDS.Shepard)
  # grid::grid.draw(plot.VI.NMDS.biplot[[3]])
  dev.off()
  
  plot.VI.NMDS.biplot = list()
  for (k in 1:2)
  {
    if (k==1)
    {
      axis1 = 1
      axis2 = 2
    } else if (k==2)
    {
      axis1 = 2
      axis2 = 3
    } else if (k==3)
    {
      axis1 = 3
      axis2 = 4
    }
    plot.VI.NMDS.biplot[[k]] = ggplot(data = data.frame(x = VI.NMDS2$points[,axis1][-69],
                                                   y = VI.NMDS2$points[,axis2][-69], 
                                                   col = log10(size_relativeAbund)[selected_groups & diversity>div_threshold][-69],
                                                   size = log10(as.vector(diversity))[selected_groups & diversity>div_threshold][-69])) + 
      geom_point(aes(x,y, colour = col, size = size), na.rm=T) +
      theme_bw() +
      # ggtitle(LETTERS[1]) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      # ggtitle("Norm. V.I. between clades' biogeographies") +
      theme(axis.title=element_text(size=22),
            axis.text=element_text(size=22),
            text = element_text(size=22),
            plot.title=element_text(hjust=0, size=24),
            plot.margin=unit(c(0.5,12,0.5,3),"mm"),
            legend.position="bottom",
            legend.text=element_text(size=20),
            legend.title=element_text(size=22)) +
      labs(x=paste("NMDS axis",axis1), y=paste("NMDS axis",axis2)) +
      # geom_vline(xintercept = prop_within_OTU_null, linetype = "dashed") +
      scale_colour_gradient2(guide = "colourbar", name = expression("Body size ("*mu*"m)"), low = "blue", high = "red", mid = "grey",
                             midpoint = median(log10(size_relativeAbund)[selected_groups & diversity>div_threshold][-69]),
                             breaks = log10(c(30,100,700)), labels = c(30,100,700)) +
      scale_size(range = c(2, 5), name = "Diversity (#OTUs)", 
                 breaks = log10(c(101,1000,70000)), 
                 labels = c("100","1,000","70,000")) +
      guides(colour = guide_colorbar(barwidth = 12, barheight = 0.8, title.position="top"),
             size = guide_legend(title.position="top"))
  }
  # pdf(paste0(figure_folder,"/NMDS_axis1vs2_colSizeDiversity_k=4.pdf"))
  pdf(paste0(figure_folder,"/NMDS_axis1vs2_colSizeDiversity_k=3_2.pdf"))
  grid::grid.draw(plot.VI.NMDS.biplot[[1]])
  grid::grid.draw(plot.VI.NMDS.biplot[[2]])
  # grid::grid.draw(plot.VI.NMDS.biplot[[3]])
  dev.off()
}

############################# Main text figures:

# Old fig. 2 autocorr. scale:
{
  Moran.step = readRDS(paste0(results_folder,"/MoranI_spatialAutocorr_20increment_2plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
  increment.allTaxa = Moran.step[[1]][[1]]
  I_step_observed_w.mean.allTaxa = Moran.step[[1]][[2]]
  # I_step_relative_w.mean.allTaxa = Moran.step[[1]][[3]]
  # I_step_p.value_w.mean.allTaxa = Moran.step[[1]][[4]]
  
  # col.vect = rep("black",nb_increment)
  # col.vect[I_step_p.value_w.mean.allTaxa[,1] > 0.05 | is.na(I_step_p.value_w.mean.allTaxa[,1])] = "red"
  
  y.smooth = smooth.spline(increment.allTaxa[,1],I_step_observed_w.mean.allTaxa[,1],df=7)$y
  lin.regr = lm(y.smooth[2:7] ~ increment.allTaxa[2:7,1])
  b.SUR = summary(lin.regr)$coefficients[1,1]
  a.SUR = summary(lin.regr)$coefficients[2,1]
  charac_scale.allTaxa.SUR = -b.SUR/a.SUR
  
  y.smooth = smooth.spline(increment.allTaxa[,2],I_step_observed_w.mean.allTaxa[,2],df=7)$y
  lin.regr = lm(y.smooth[2:7] ~ increment.allTaxa[2:7,2])
  b.DCM = summary(lin.regr)$coefficients[1,1]
  a.DCM = summary(lin.regr)$coefficients[2,1]
  charac_scale.allTaxa.DCM = -b.DCM/a.DCM
  
  Moran.step.observed.SUR.DCM.plot.w.spline = ggplot(data = data.frame(x = c(increment.allTaxa[,1],increment.allTaxa[,2]), 
                                                                       y = c(I_step_observed_w.mean.allTaxa[,1],I_step_observed_w.mean.allTaxa[,2]), 
                                                                       y.smooth = c(smooth.spline(increment.allTaxa[,1],I_step_observed_w.mean.allTaxa[,1],df=7)$y,
                                                                                    smooth.spline(increment.allTaxa[,2],I_step_observed_w.mean.allTaxa[,2],df=7)$y),
                                                                       col = c(rep("Surface",length(increment.allTaxa[,1])),rep("DCM",length(increment.allTaxa[,2]))))) +
    labs(x="Spatial scale (km)", y="\nSpatial autocorrelation") +
    geom_abline(slope = a.SUR,intercept = b.SUR, linetype = "dashed", size = 1) +
    geom_abline(slope = a.DCM, intercept = b.DCM, linetype = "dashed", size = 1, col = "red") +
    geom_point(aes(x,y,colour=factor(col)), size = 1.8) +
    geom_line(aes(x,y.smooth, colour=factor(col)), size = 1) +
    # geom_point(aes(x,y), size = 0.6, col = col.vect) +
    # coord_trans(x = "log10") +
    # scale_x_log10() +
    geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
    theme_bw() +
    # ggtitle(paste(taxon,"-",ifelse(taxon=="AllTaxa",sum(diversity),as.vector(diversity)[i_taxon]),"OTUs")) +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          plot.margin=unit(c(1,1,1,0.5),"mm"),
          legend.position=c(0.7,0.7),
          legend.text=element_text(size=22),
          legend.title=element_blank()) +
    # breaks: reordering legend labels to have Surface first:
    scale_colour_manual(values=c("red","black"),breaks=c("Surface","DCM"))
  # guides(fill = guide_legend(title = NULL, direction = "vertical"))
  
  pdf(paste0(figure_folder,"/Fig2_autocorr.scale.panel_SUR.DCM.pdf"),width=6,height=7)
  print(Moran.step.observed.SUR.DCM.plot.w.spline)
  dev.off()
}

# Old fig. 2 JSD Surf. DCM plot
{
  JSD_object = readRDS(paste0(results_folder,"/JSD_2plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
  JSD.allTaxa = JSD_object[[1]][[1]]
  coord_JSD.allTaxa = JSD_object[[1]][[2]]
  
  plot.x = coord_JSD.allTaxa[,1]
  plot.x.order = sort.int(coord_JSD.allTaxa[,1], decreasing = F, index.return = T)
  plot.y = JSD.allTaxa[,2]
  
  JSD.obs.over.exp.latitude.plot = ggplot(data = data.frame(x = plot.x, y = plot.y, x.smooth = plot.x.order$x, y.smooth = smooth.spline(plot.x.order$x,plot.y[plot.x.order$ix],df=5)$y)) +
    # geom_smooth(aes(x,y),method='lm') +
    geom_point(aes(x,y), size = 1.8) +
    geom_line(aes(x.smooth,y.smooth), size = 1) +
    theme_bw() +
    # ggtitle(paste(taxon,"-",ifelse(taxon=="AllTaxa",sum(diversity),as.vector(diversity)[i_taxon]),"OTUs")) +
    theme(axis.title=element_text(size=22),
          axis.text = element_text(size=22),
          # plot.title=element_text(hjust=0, size=15),
          plot.margin=unit(c(1,1,1,0.5),"mm")) +
    labs(x="Latitude", y="Observed over expected\n J-S divergence between Surface and DCM")
  
  pdf(paste0(figure_folder,"/Fig2_JSD.SUR.DCM.pdf"),width=6,height=7)
  print(JSD.obs.over.exp.latitude.plot)
  dev.off()
}

# Old fig. 2 barplot:
{
  nb_topics = 16
  nb_iter = 1000
  nb_real = 100
  thin = 25
  burnin = 2000
  
  data.folder.allTaxa = paste0(data_folder.workspace3,"/18S_V9_TARA_CompleteSizeRange_byStationByDepth_AllTaxa",noArcticNoBiomark_insert,noLagoon_insert)
  load(paste0(data.folder.allTaxa,"/coord.Rdata"))
  spatial_topicmix_kriged = readRDS(paste0(data.folder.allTaxa,"/Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics",
                                           nb_topics,"_nb_iter",nb_iter,"_nb_real",nb_real,"_meanPosteriorDistributedLlh_thin",thin,"_burnin",burnin,"_occurrence/1st_closestToMean_realization/Spatial_topicmix_kriged.rds"))
  # selecting the z.pred columns in all topics:
  documents = unlist(lapply(spatial_topicmix_kriged,function(g) g$z.pred))
  # setting one topic per column
  documents = matrix(documents,ncol=nb_topics,dimnames=list(rownames(spatial_topicmix_kriged[[1]]),paste0("assemblage",1:nb_topics)))
  
  color.pal1 = colorRampPalette(c("#7F0000","red","darkorange1","darkgoldenrod1","yellow"),space = "Lab")
  color.pal2 = colorRampPalette(c("#7FFF7F","darkgreen"),space = "Lab")
  color.pal3 = colorRampPalette(c("azure4","antiquewhite3","aliceblue"),space = "Lab")
  color.pal4 = colorRampPalette(c("#00007F","#007FFF","cyan"),space = "Lab")
  col = c(color.pal1(6),color.pal2(3),color.pal3(3),color.pal4(4))
  topic_reordering = c(1,3,9,7,
                       10,8,4,6,
                       11,16,15,12,
                       14,5,13,2)
  
  stations_depths = as.data.frame(matrix(nrow=nrow(documents),ncol=2))
  colnames(stations_depths) = c("Station","Depth")
  for (station_depth_index in 1:nrow(documents))
  {
    stations_depths[station_depth_index,] = c(strsplit(rownames(documents)[station_depth_index],split=" ")[[1]][1],
                                              strsplit(rownames(documents)[station_depth_index],split=" ")[[1]][2])
  }
  
  # pdf(paste0(figure_folder,"/AllTaxa_16t_Gibbs_SurDCM_barplot.pdf"),height=4,width = 8)
  pdf(paste0(figure_folder,"/AllTaxa_16t_Gibbs_SurDCM_barplot_latitudeOrder.pdf"),height=4,width = 8)
  par(lwd = 0.3)
  barplot(t(documents[stations_depths[,2] == "SUR",topic_reordering][sort.int(coord$y[stations_depths[,2] == "SUR"],index.return = T, decreasing = T)$ix,]),
          col=col, axisnames=F, legend.text = F, args.legend = list(x="topright",bty = "n"),ann=F, space = 0)
  barplot(t(documents[stations_depths[,2] == "DCM",topic_reordering][sort.int(coord$y[stations_depths[,2] == "DCM"],index.return = T, decreasing = T)$ix,]),
          col=col,  axisnames=F, legend.text = F, args.legend = list(x="topright",bty = "n"),ann=F, space = 0)
  barplot(t(documents[,topic_reordering][sort.int(coord$y,index.return = T, decreasing = T)$ix,]),
          col=col,  axisnames=F, legend.text = F, args.legend = list(x="topright",bty = "n"),ann=F, space = 0)
  dev.off()
}

# Old fig. 3 - v1:
{
  dominant_function0 = dominant_function0[-69]
  # alpha = rep(1,length(taxo_groups[selected_groups]))
  # alpha[dominant_function0 %in% c("copepoda","pteropoda")] = 0 
  # Changes how the factors are stored (order of levels()) so that geom_boxplot plots them by deceasing number of groups:
  # dominant_function0 = factor(dominant_function0,names(sort(table(as.factor(dominant_function0)),decreasing = T)))
  # Changes how the factors are stored (order of levels()) so that geom_boxplot plots them in the specified order:
  if (div_threshold == 100)
  {
    dominant_function0 = factor(dominant_function0,c("Bacillariophyta","Other phototrophs","Collodaria","Phagotrophs","Dinophyceae","Gel. carn. filterers","Pteropoda","Copepoda","Other metazoa","Parasites"))
    point_groups = c("Dinophyceae","Copepoda","Collodaria","Bacillariophyta","Pteropoda")
  } else if (div_threshold == 1000)
  {
    dominant_function0 = factor(dominant_function0,c("Bacillariophyta","Other phototrophs","Collodaria","Phagotrophs","Dinophyceae","Gel. carn. filterers","Copepoda","Parasites"))
    point_groups = c("Dinophyceae","Copepoda","Collodaria","Bacillariophyta")
  } 
  
  functions0 = c("Collodaria","Pteropoda","Copepoda","Gel. carn. filterers","Dinophyceae","Bacillariophyta","Other phototrophs","Phagotrophs","Other metazoa","Parasites")
  # "Bacillariophyta"   "Collodaria"   "Copepoda"   "Dinophyceae"   "Gel. carn. filterers" "Other phototrophs"    "Parasites"    "Phagotrophs" 
  ten_colors = c("darkorange","cadetblue","darkblue","darkturquoise","deeppink1","darkgreen","chartreuse2","firebrick2","dodgerblue1","darkgoldenrod1")
  
  library(gridExtra)
  pdf(paste0(figure_folder,"/Normalized_PCoA.axes1.2.vs.struct.div.size_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2*3/2,height=12/3*4/1.2*3/2)
  for (i_case in 1:2)
  {
    # Axis 1 vs 2 functions:
    ########################
    plot.VI.PCoA.biplot = ggplot(data = data.frame(x = NormalizedVI_pcoa[[3]][,1], y = NormalizedVI_pcoa[[3]][,2])) + 
      geom_point(aes(x,y, colour = dominant_function0), size = 3, na.rm=T) +
      theme_bw() +
      ggtitle(LETTERS[1]) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      # ggtitle("Norm. V.I. between clades' biogeographies") +
      theme(axis.title=element_text(size=22),
            axis.text=element_text(size=22),
            text = element_text(size=22),
            plot.title=element_text(hjust=0, size=24),
            plot.margin=unit(c(0.5,15,0.5,0.5),"mm"),
            legend.position="bottom",
            legend.text=element_text(size=22),
            legend.title=element_text(size=16, hjust = 10)) +
      labs(x="PCoA axis 1", y="PCoA axis 2") +
      # geom_vline(xintercept = prop_within_OTU_null, linetype = "dashed") +
      scale_colour_manual(values = setNames(ten_colors,functions0), 
                          na.value = "grey50", na.translate = FALSE, guide = "colourbar", name = "Functional group") +
      # geom_smooth(method='lm') +
      guides(colour = guide_legend(title=NULL, nrow = 5, ncol = 2))
    g1 = ggplotGrob(plot.VI.PCoA.biplot)
    
    # Plots axes vs. structure:
    ###########################
    plot.norm.VI.PCoA.axis1.vs.spatial.autocorr = ggplot(data = data.frame(x = I_square.observed_w.mean[selected_groups & diversity>div_threshold,i_case][-69],
                                                                           y = NormalizedVI_pcoa[[3]][,1])) +
      geom_point(aes(x,y)) +
      # scale_x_log10() +
      # scale_y_log10() +
      theme_bw() +
      ggtitle(LETTERS[2]) +
      # geom_hline(yintercept = 1, linetype="dashed") +
      # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
      # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
      theme(axis.text = element_text(size=22),
            axis.title=element_text(size=22),
            # axis.title.x=element_text(vjust = 45),
            plot.title=element_text(hjust=0, size=24),
            plot.margin=unit(c(0.5,7,0.5,0.5),"mm")) +
      #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
      labs(x=paste("Short-distance spatial autocorr.",c("(Surface)","(DCM)")[i_case]), y="PCoA axis 1") +
      geom_smooth(aes(x,y),method='lm',col="black")
    g2 = ggplotGrob(plot.norm.VI.PCoA.axis1.vs.spatial.autocorr)
    
    plot.norm.VI.PCoA.axis1.vs.SUR.DCM.sim = ggplot(data = data.frame(x = 1 - SUR.DCM_Normalized.VI[selected_groups & diversity>div_threshold][-69],
                                                                      y = NormalizedVI_pcoa[[3]][,1])) +
      geom_point(aes(x,y)) +
      # scale_x_log10() +
      # scale_y_log10() +
      theme_bw() +
      ggtitle(LETTERS[3]) +
      # geom_hline(yintercept = 1, linetype="dashed") +
      # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
      # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
      theme(axis.text = element_text(size=22),
            axis.title=element_text(size=22),
            # axis.title.x=element_text(vjust = 45),
            plot.title=element_text(hjust=0, size=24),
            plot.margin=unit(c(0.5,7,0.5,0.5),"mm")) +
      #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
      labs(x="Similarity between surface and DCM", y="PCoA axis 1") +
      geom_smooth(aes(x,y),method='lm',col="black")
    g3 = ggplotGrob(plot.norm.VI.PCoA.axis1.vs.SUR.DCM.sim)
    
    plot.norm.VI.PCoA.axis2.vs.spatial.scale = ggplot(data = data.frame(x = charac_scale[selected_groups & diversity>div_threshold,i_case][-69],
                                                                        y = NormalizedVI_pcoa[[3]][,2])) +
      geom_point(aes(x,y)) +
      # scale_x_log10() +
      # scale_y_log10() +
      theme_bw() +
      ggtitle(LETTERS[4]) +
      # geom_hline(yintercept = 1, linetype="dashed") +
      # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
      # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
      theme(axis.text = element_text(size=22),
            axis.title=element_text(size=22),
            # axis.title.x=element_text(vjust = 45),
            plot.title=element_text(hjust=0, size=24),
            plot.margin=unit(c(0.5,7,0.5,0.5),"mm")) +
      #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
      labs(x=paste("Charac. scale of spatial autocorr.",c("(Surface)","(DCM)")[i_case]), y="PCoA axis 2") +
      geom_smooth(aes(x,y),method='lm',col="black")
    g4 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.spatial.scale)
    
    # Plots str. vs. diversity:
    ##########################
    x = as.vector(diversity)[selected_groups & diversity>div_threshold][-69]
    y = I_square.observed_w.mean[selected_groups & diversity>div_threshold,i_case][-69]
    # lin.regr = lm(y ~ x)
    # b = summary(lin.regr)$coefficients[1,1]
    # a = summary(lin.regr)$coefficients[2,1]
    plot.spatial.autocorr.vs.diversity = ggplot(data = data.frame(x = x, y = y)) +
      geom_point(aes(x,y)) +
      # geom_point(aes(x,y2), col = "red") +
      scale_x_log10() +
      # scale_y_log10() +
      # geom_abline(intercept = b, slope = a) +
      theme_bw() +
      ggtitle(LETTERS[5]) +
      # geom_hline(yintercept = 1, linetype="dashed") +
      # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
      # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
      theme(axis.text = element_text(size=22),
            axis.title=element_text(size=22),
            # axis.title.x=element_text(vjust = 45),
            plot.title=element_text(hjust=0, size=24),
            plot.margin=unit(c(0.5,7,0.5,0.5),"mm")) +
      #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
      labs(x="Clade diversity", y=paste("Short-distance spatial autocorr.",c("(Surface)","(DCM)")[i_case])) +
      geom_smooth(data = data.frame(x = x[x<2000], y = y[x<2000]), aes(x,y),method='lm',col="black",se=T) +
      geom_segment(aes(x = 2000, xend = max(x), y = 0.705, yend = 0.705), linetype="dashed")
    # geom_smooth(aes(x,y2),method='lm',col="red")
    g5 = ggplotGrob(plot.spatial.autocorr.vs.diversity)
    
    x = as.vector(diversity)[selected_groups & diversity>div_threshold][-69]
    y = 1 - SUR.DCM_Normalized.VI[selected_groups & diversity>div_threshold][-69]
    # lin.regr = lm(y ~ log10(x))
    # b = summary(lin.regr)$coefficients[1,1]
    # a = summary(lin.regr)$coefficients[2,1]
    plot.SUR.DCM.sim.vs.diversity = ggplot(data = data.frame(x = x, y = y)) +
      # geom_point(aes(x,y1)) +
      geom_point(aes(x,y)) +
      scale_x_log10() +
      # scale_y_log10() +
      # geom_abline(intercept = b, slope = a) +
      theme_bw() +
      ggtitle(LETTERS[6]) +
      # geom_hline(yintercept = 1, linetype="dashed") +
      # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
      # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
      theme(axis.text = element_text(size=22),
            axis.title=element_text(size=22),
            # axis.title.x=element_text(vjust = 45),
            plot.title=element_text(hjust=0, size=24),
            plot.margin=unit(c(0.5,7,0.5,0.5),"mm")) +
      #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
      labs(x="Clade diversity", y="Similarity between surface and DCM") +
      geom_smooth(data = data.frame(x = x[x<2000], y = y[x<2000]), aes(x,y),method='lm',col="black",se=T) +
      geom_segment(aes(x = 2000, xend = max(x), y = 0.3175, yend = 0.3175), linetype="dashed")
    # geom_smooth(aes(x,y),method='lm')
    g6 = ggplotGrob(plot.SUR.DCM.sim.vs.diversity)
    
    # x = as.vector(diversity)[selected_groups & diversity>div_threshold][-69]
    # y = charac_scale[selected_groups & diversity>div_threshold,i_case][-69]
    # plot.spatial.scale.vs.diversity = ggplot(data = data.frame(x = x, y = y)) +
    #   # geom_point(aes(x,y1)) +
    #   geom_point(aes(x,y)) +
    #   scale_x_log10() +
    #   # scale_y_log10() +
    #   # geom_abline(intercept = b, slope = a) +
    #   theme_bw() +
    #   ggtitle(LETTERS[6]) +
    #   # geom_hline(yintercept = 1, linetype="dashed") +
    #   # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    #   # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    #   theme(axis.text = element_text(size=22),
    #         axis.title=element_text(size=22),
    #         # axis.title.x=element_text(vjust = 45),
    #         plot.title=element_text(hjust=0, size=24),
    #         plot.margin=unit(c(0.5,7,0.5,0.5),"mm")) +
    #   #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    #   labs(x="Clade diversity", y=paste("Charac. scale of spatial autocorr. -",c("Surf.","DCM")[i_case])) +
    #   geom_smooth(data = data.frame(x = x[x<2000], y = y[x<2000]), aes(x,y), method='lm',col="black",se=T)
    #   # geom_smooth(aes(x,y),method='lm')
    # g7 = ggplotGrob(plot.spatial.scale.vs.diversity)
    
    # Body size/boxplots:
    ##########################
    
    plot.norm.VI.PCoA.axis2.vs.size = ggplot(data = data.frame(x = size_relativeAbund[selected_groups & diversity>div_threshold][-69],
                                                               y = NormalizedVI_pcoa[[3]][,2])) +
      # geom_point(aes(x,y1)) +
      geom_point(aes(x,y)) +
      scale_x_log10() +
      # scale_y_log10() +
      theme_bw() +
      ggtitle(LETTERS[7]) +
      # geom_hline(yintercept = 1, linetype="dashed") +
      # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
      # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
      theme(axis.text = element_text(size=22),
            axis.title=element_text(size=22),
            # axis.title.x=element_text(vjust = 45),
            plot.title=element_text(hjust=0, size=24),
            plot.margin=unit(c(0.5,7,0.5,0.5),"mm")) +
      #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
      labs(x="Clade mean body size", y="PCoA axis 2") +
      geom_smooth(aes(x,y),method='lm',col="black")
    g7 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.size)
    
    boxplot.norm.VI.PCoA.axis2.functions = ggplot(data = data.frame(x = dominant_function0[!is.na(dominant_function0) & !dominant_function0 %in% point_groups],
                                                                    y = NormalizedVI_pcoa[[3]][,2][!is.na(dominant_function0) & !dominant_function0 %in% point_groups])) +
      scale_x_discrete(limits=levels(dominant_function0)) +
      geom_boxplot(aes(x,y)) +
      geom_point(data = data.frame(x = factor(point_groups),
                                   y = NormalizedVI_pcoa[[3]][,2][dominant_function0 %in% point_groups]),
                 aes(x,y)) +
      theme_bw() +
      ggtitle(LETTERS[8]) +
      # scale_y_log10() +
      # geom_hline(yintercept = 1, linetype = "dashed") +
      theme(axis.title=element_text(size=22),
            axis.text=element_text(size=22),
            plot.title=element_text(hjust=0, size=24),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            plot.margin=unit(c(0.5,7,0.5,0.5),"mm")) +
      labs(x="", y="PCoA axis 2")
    g8 = ggplotGrob(boxplot.norm.VI.PCoA.axis2.functions)
    
    boxplot.spatial.scale.functions = ggplot(data = data.frame(x = dominant_function0[!is.na(dominant_function0) & !dominant_function0 %in% point_groups],
                                                               y = charac_scale[selected_groups & diversity>div_threshold,i_case][-69][!is.na(dominant_function0) & !dominant_function0 %in% point_groups])) +
      scale_x_discrete(limits=levels(dominant_function0)) +
      geom_boxplot(aes(x,y)) +
      geom_point(data = data.frame(x = factor(point_groups),
                                   y = charac_scale[selected_groups & diversity>div_threshold,i_case][-69][dominant_function0 %in% point_groups]),
                 aes(x,y)) +
      theme_bw() +
      ggtitle(LETTERS[9]) +
      # scale_y_log10() +
      # geom_hline(yintercept = 1, linetype = "dashed") +
      theme(axis.title=element_text(size=22),
            axis.text=element_text(size=22),
            plot.title=element_text(hjust=0, size=24),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            plot.margin=unit(c(0.5,7,0.5,0.5),"mm")) +
      labs(x="", y=paste("Charac. scale of spatial autocorr.",c("(Surface)","(DCM)")[i_case]))
    g9 = ggplotGrob(boxplot.spatial.scale.functions)
    
    grid.arrange(grobs = list(g1,g2,g3,g4,g5,g6,g7,g8,g9),
                 # widths = rep(1,7),
                 layout_matrix = rbind(c(1, 2, 3),c(4, 5, 6),c(7, 8, 9)))
    # grid::grid.draw(gtable_rbind(gtable_cbind(g1, g2, g3), gtable_cbind(g4, g5), gtable_cbind(g6, g7)),recording=F)    
  }
  #grid::grid.newpage()
  dev.off()
}

# Old fig. 3 - v2:
{
  div_threshold = 100
  
  dominant_function0 = dominant_function0[-69]
  # alpha = rep(1,length(taxo_groups[selected_groups]))
  # alpha[dominant_function0 %in% c("copepoda","pteropoda")] = 0 
  # Changes how the factors are stored (order of levels()) so that geom_boxplot plots them by deceasing number of groups:
  # dominant_function0 = factor(dominant_function0,names(sort(table(as.factor(dominant_function0)),decreasing = T)))
  # Changes how the factors are stored (order of levels()) so that geom_boxplot plots them in the specified order:
  if (div_threshold == 100)
  {
    dominant_function0 = factor(dominant_function0,c("Bacillariophyta","Other phototrophs","Collodaria","Phagotrophs","Dinophyceae","Gel. carn. filterers","Pteropoda","Copepoda","Other metazoa","Parasites"))
    point_groups = c("Dinophyceae","Copepoda","Collodaria","Bacillariophyta","Pteropoda")
  } else if (div_threshold == 1000)
  {
    dominant_function0 = factor(dominant_function0,c("Bacillariophyta","Other phototrophs","Collodaria","Phagotrophs","Dinophyceae","Gel. carn. filterers","Copepoda","Parasites"))
    point_groups = c("Dinophyceae","Copepoda","Collodaria","Bacillariophyta")
  } 
  
  # file_name = paste0(figure_folder,"/Fig3.1_ggplot",abiotic_pca_insert,biotic_pca_insert,stdzation_insert,
  #                    "_independentVariableSelection_selected",if (div_threshold == 100) "100+1" else if (div_threshold == 1000) "1000",surf_DCM_insert,"_noNegativeAdjR2.pdf")
  # 
  
  functions0 = c("Collodaria","Pteropoda","Copepoda","Gel. carn. filterers","Dinophyceae","Bacillariophyta","Other phototrophs","Phagotrophs","Other metazoa","Parasites")
  # "Bacillariophyta"   "Collodaria"   "Copepoda"   "Dinophyceae"   "Gel. carn. filterers" "Other phototrophs"    "Parasites"    "Phagotrophs" 
  ten_colors = c("darkorange","cadetblue","darkblue","darkturquoise","deeppink1","darkgreen","chartreuse2","firebrick2","dodgerblue1","darkgoldenrod1")
  
  #############
  # SUR       #
  #############
  i_case = 1
  surf_DCM_insert = "_SUR"
  library(gridExtra)
  
  # Axis 1 vs 2 functions:
  ########################
  plot.VI.PCoA.biplot = ggplot(data = data.frame(x = NormalizedVI_pcoa[[3]][,1], y = NormalizedVI_pcoa[[3]][,2], 
                                                 col = log10(size_relativeAbund)[selected_groups & diversity>div_threshold][-69],
                                                 size = log10(as.vector(diversity))[selected_groups & diversity>div_threshold][-69])) + 
    geom_point(aes(x,y, colour = col, size = size), na.rm=T) +
    theme_bw() +
    ggtitle(LETTERS[1]) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    # ggtitle("Norm. V.I. between clades' biogeographies") +
    theme(axis.title=element_text(size=22),
          axis.text=element_text(size=22),
          text = element_text(size=22),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm"),
          legend.position="bottom",
          legend.text=element_text(size=20),
          legend.title=element_text(size=22)) +
    labs(x="PCoA axis 1", y="PCoA axis 2") +
    # geom_vline(xintercept = prop_within_OTU_null, linetype = "dashed") +
    scale_colour_gradient2(guide = "colourbar", name = expression("Body size ("*mu*"m)"), low = "blue", high = "red", mid = "grey",
                           midpoint = median(log10(size_relativeAbund)[selected_groups & diversity>div_threshold][-69]),
                           breaks = log10(c(30,100,700)), labels = c(30,100,700)) +
    scale_size(range = c(2, 5), name = "Diversity (#OTUs)", 
               breaks = log10(c(101,1000,70000)), 
               labels = c("100","1,000","70,000")) +
    guides(colour = guide_colorbar(barwidth = 12, barheight = 0.8, title.position="top"),
           size = guide_legend(title.position="top"))
  # guides(colour = guide_legend(title = NULL, nrow = 3))
  
  # Plots axes vs. structure:
  ###########################
  x = I_square.observed_w.mean[selected_groups & diversity>div_threshold,i_case][-69]
  y = NormalizedVI_pcoa[[3]][,1]
  cor.test = cor.test(x,y)
  plot.norm.VI.PCoA.axis1.vs.spatial.autocorr = ggplot(data = data.frame(x=x, y=y)) +
    geom_point(aes(x,y)) +
    # scale_x_log10() +
    # scale_y_log10() +
    theme_bw() +
    ggtitle(LETTERS[if (i_case == 1) 2 else 1]) +
    # geom_hline(yintercept = 1, linetype="dashed") +
    # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          # axis.title.x=element_text(vjust = 45),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm")) +
    #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    labs(x=paste("Short-distance spatial autocorr.",c("(Surface)","(DCM)")[i_case]), y="PCoA axis 1") +
    geom_smooth(aes(x,y),method='lm',col="black") +
    annotate(geom="text", 
             x=0.8*(max(x)-min(x))+min(x),
             y=0.2*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=7)
  
  x = 1 - SUR.DCM_Normalized.VI[selected_groups & diversity>div_threshold][-69]
  y = NormalizedVI_pcoa[[3]][,1]
  cor.test = cor.test(x,y)
  plot.norm.VI.PCoA.axis1.vs.SUR.DCM.sim = ggplot(data = data.frame(x=x, y=y)) +
    geom_point(aes(x,y)) +
    # scale_x_log10() +
    # scale_y_log10() +
    theme_bw() +
    ggtitle(LETTERS[3]) +
    # geom_hline(yintercept = 1, linetype="dashed") +
    # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          # axis.title.x=element_text(vjust = 45),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,10,0.5,3),"mm")) +
    #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    labs(x="Similarity between surface and DCM", y="PCoA axis 1") +
    geom_smooth(aes(x,y),method='lm',col="black") +
    annotate(geom="text", 
             x=0.8*(max(x)-min(x))+min(x),
             y=0.2*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=7)
  
  x = charac_scale[selected_groups & diversity>div_threshold,i_case][-69]
  y = NormalizedVI_pcoa[[3]][,2]
  cor.test = cor.test(x,y)
  plot.norm.VI.PCoA.axis2.vs.spatial.scale = ggplot(data = data.frame(x=x, y=y)) +
    geom_point(aes(x,y)) +
    # scale_x_log10() +
    # scale_y_log10() +
    theme_bw() +
    ggtitle(LETTERS[if (i_case == 1) 4 else 2]) +
    # geom_hline(yintercept = 1, linetype="dashed") +
    # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          # axis.title.x=element_text(vjust = 45),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm")) +
    #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    labs(x=paste("Spatial autocorr. scale",c("(km, Surface)","(km, DCM)")[i_case]), y="PCoA axis 2") +
    geom_smooth(aes(x,y),method='lm',col="black") +
    annotate(geom="text", 
             x=0.8*(max(x)-min(x))+min(x),
             y=0.2*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=7)
  
  # Plot axis 1 vs. diversity:
  ##########################
  
  x = as.vector(diversity)[selected_groups & diversity>div_threshold][-69]
  y = NormalizedVI_pcoa[[3]][,1]
  cor.test = cor.test(log10(x),y)
  plot.VI.PCoA.axis1.vs.diversity = ggplot(data = data.frame(x = x,y = y)) +
    geom_point(aes(x,y)) +
    scale_x_log10() +
    # scale_y_log10() +
    # geom_abline(intercept = b, slope = a) +
    theme_bw() +
    ggtitle(LETTERS[5]) +
    # geom_hline(yintercept = 1, linetype="dashed") +
    # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          # axis.title.x=element_text(vjust = 45),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm")) +
    #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    labs(x="Clade diversity", y="PCoA axis 1") +
    geom_smooth(data = data.frame(x = x[x<2000], y = y[x<2000]), aes(x,y),method='lm',col="black",se=T, inherit.aes = T) +
    geom_segment(aes(x = 2000, xend = max(x), y = 0.167, yend = 0.167), linetype="dashed") +
    annotate(geom="text", 
             x=(max(x)/min(x))^0.45*min(x),
             y=0.2*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=7)
  
  # x = as.vector(diversity)[selected_groups & diversity>div_threshold][-69]
  # y = charac_scale[selected_groups & diversity>div_threshold,i_case][-69]
  # plot.spatial.scale.vs.diversity = ggplot(data = data.frame(x = x, y = y)) +
  #   # geom_point(aes(x,y1)) +
  #   geom_point(aes(x,y)) +
  #   scale_x_log10() +
  #   # scale_y_log10() +
  #   # geom_abline(intercept = b, slope = a) +
  #   theme_bw() +
  #   ggtitle(LETTERS[6]) +
  #   # geom_hline(yintercept = 1, linetype="dashed") +
  #   # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
  #   # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
  #   theme(axis.text = element_text(size=22),
  #         axis.title=element_text(size=22),
  #         # axis.title.x=element_text(vjust = 45),
  #         plot.title=element_text(hjust=0, size=24),
  #         plot.margin=unit(c(0.5,7,0.5,0.5),"mm")) +
  #   #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
  #   labs(x="Clade diversity", y=paste("Charac. scale of spatial autocorr. -",c("Surf.","DCM")[i_case])) +
  #   geom_smooth(data = data.frame(x = x[x<2000], y = y[x<2000]), aes(x,y), method='lm',col="black",se=T)
  #   # geom_smooth(aes(x,y),method='lm')
  # g7 = ggplotGrob(plot.spatial.scale.vs.diversity)
  
  # Body size/boxplots:
  ##########################
  
  x = size_relativeAbund[selected_groups & diversity>div_threshold][-69]
  y = NormalizedVI_pcoa[[3]][,2]
  cor.test = cor.test(log10(x),y)
  plot.norm.VI.PCoA.axis2.vs.size = ggplot(data = data.frame(x = x, y = y)) +
    # geom_point(aes(x,y1)) +
    geom_point(aes(x,y)) +
    scale_x_log10() +
    # scale_y_log10() +
    theme_bw() +
    ggtitle(LETTERS[6]) +
    # geom_hline(yintercept = 1, linetype="dashed") +
    # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          # axis.title.x=element_text(vjust = 45),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm")) +
    #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    labs(x=expression("Clade mean body size ("*mu*"m)"), y="PCoA axis 2") +
    geom_smooth(aes(x,y),method='lm',col="black") +
    annotate(geom="text", 
             x=(max(x)/min(x))^0.45*min(x),
             y=0.2*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=7)
  
  boxplot.norm.VI.PCoA.axis2.functions = ggplot(data = data.frame(x = dominant_function0[!is.na(dominant_function0) & !dominant_function0 %in% point_groups],
                                                                  y = NormalizedVI_pcoa[[3]][,2][!is.na(dominant_function0) & !dominant_function0 %in% point_groups])) +
    scale_x_discrete(limits=levels(dominant_function0)) +
    geom_boxplot(aes(x,y)) +
    geom_point(data = data.frame(x = factor(point_groups),
                                 y = NormalizedVI_pcoa[[3]][,2][dominant_function0 %in% point_groups]),
               aes(x,y)) +
    theme_bw() +
    ggtitle(LETTERS[7]) +
    # scale_y_log10() +
    # geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.title=element_text(size=22),
          axis.text=element_text(size=22),
          plot.title=element_text(hjust=0, size=24),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin=unit(c(0.5,12,0.5,3),"mm")) +
    labs(x="", y="PCoA axis 2")

  # Plotting axis 1 vs. 2, plus 6 axes correlation panels:
  pdf(paste0(figure_folder,"/Normalized_PCoA.axes1.2.vs.struct.div.size_noFunction_",div_threshold,"plusOTUs",surf_DCM_insert,noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2*3/2,height=12/3*4/1.2*3/2)
  g1 = ggplotGrob(plot.VI.PCoA.biplot)
  g2 = ggplotGrob(plot.norm.VI.PCoA.axis1.vs.spatial.autocorr)
  g3 = ggplotGrob(plot.norm.VI.PCoA.axis1.vs.SUR.DCM.sim)
  g4 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.spatial.scale)
  g5 = ggplotGrob(plot.VI.PCoA.axis1.vs.diversity)
  g6 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.size)
  g7 = ggplotGrob(boxplot.norm.VI.PCoA.axis2.functions)
  grid.arrange(grobs = list(g1,g2,g3,g4,g5,g6,g7),
                   # widths = rep(1,7),
               layout_matrix = rbind(c(1, 2, 3),c(4, 5, 6),c(7, NA, NA)))
    # grid::grid.draw(gtable_rbind(gtable_cbind(g1, g2, g3), gtable_cbind(g4, g5), gtable_cbind(g6, g7)),recording=F)    
  #grid::grid.newpage()
  dev.off()
  
  # Plotting 6 axes correlation panels:
  pdf(paste0(figure_folder,"/Normalized_PCoA.axes1.2.vs.struct.div.size_noAxes1.2.Scatterplot_",div_threshold,"plusOTUs",surf_DCM_insert,noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2*3/2,height=12/3*4/1.2)
  g2 = ggplotGrob(plot.norm.VI.PCoA.axis1.vs.spatial.autocorr + ggtitle(LETTERS[1]))
  g3 = ggplotGrob(plot.norm.VI.PCoA.axis1.vs.SUR.DCM.sim + ggtitle(LETTERS[2]))
  g4 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.spatial.scale + ggtitle(LETTERS[3]))
  g5 = ggplotGrob(plot.VI.PCoA.axis1.vs.diversity + ggtitle(LETTERS[4]))
  g6 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.size + ggtitle(LETTERS[5]))
  g7 = ggplotGrob(boxplot.norm.VI.PCoA.axis2.functions + ggtitle(LETTERS[6]))
  grid.arrange(grobs = list(g2,g3,g4,g5,g6,g7),
               # widths = rep(1,7),
               layout_matrix = rbind(c(1, 2, 3),c(4, 5, 6)))
  # grid::grid.draw(gtable_rbind(gtable_cbind(g1, g2, g3), gtable_cbind(g4, g5), gtable_cbind(g6, g7)),recording=F)    
  #grid::grid.newpage()
  dev.off()
  
  #############
  # DCM       #
  #############
  i_case = 2
  surf_DCM_insert = "_DCM"
  pdf(paste0(figure_folder,"/Normalized_PCoA.axes1.2.vs.struct.div.size_noFunction_",div_threshold,"plusOTUs",surf_DCM_insert,noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2,height=12/3*4/1.2/2*2)
  
  x = I_square.observed_w.mean[selected_groups & diversity>div_threshold,i_case][-69]
  y = NormalizedVI_pcoa[[3]][,1]
  cor.test = cor.test(x,y)
  plot.norm.VI.PCoA.axis1.vs.spatial.autocorr = ggplot(data = data.frame(x = x, y = y)) +
    geom_point(aes(x,y)) +
    # scale_x_log10() +
    # scale_y_log10() +
    theme_bw() +
    ggtitle(LETTERS[if (i_case == 1) 2 else 1]) +
    # geom_hline(yintercept = 1, linetype="dashed") +
    # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          # axis.title.x=element_text(vjust = 45),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm")) +
    #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    labs(x=paste("Short-distance spatial autocorr.",c("(Surface)","(DCM)")[i_case]), y="PCoA axis 1") +
    geom_smooth(aes(x,y),method='lm',col="black") +
    annotate(geom="text", 
             x=0.8*(max(x)-min(x))+min(x),
              y=0.2*(max(y)-min(y))+min(y),
              label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
              size=7)
  g2 = ggplotGrob(plot.norm.VI.PCoA.axis1.vs.spatial.autocorr)
  
  x = charac_scale[selected_groups & diversity>div_threshold,i_case][-69]
  y = NormalizedVI_pcoa[[3]][,2]
  cor.test = cor.test(x,y)
  plot.norm.VI.PCoA.axis2.vs.spatial.scale = ggplot(data = data.frame(x = x, y = y)) +
    geom_point(aes(x,y)) +
    # scale_x_log10() +
    # scale_y_log10() +
    theme_bw() +
    ggtitle(LETTERS[if (i_case == 1) 4 else 2]) +
    # geom_hline(yintercept = 1, linetype="dashed") +
    # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          # axis.title.x=element_text(vjust = 45),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm")) +
    #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    labs(x=paste("Spatial autocorr. scale",c("(km, Surface)","(km, DCM)")[i_case]), y="PCoA axis 2") +
    geom_smooth(aes(x,y),method='lm',col="black") +
    annotate(geom="text", 
             x=0.8*(max(x)-min(x))+min(x),
             y=0.2*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=7)
  g4 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.spatial.scale)
  
  x = size_relativeAbund[selected_groups & diversity>div_threshold][-69]
  y = NormalizedVI_pcoa[[3]][,1]
  cor.test = cor.test(log10(x),y)
  plot.norm.VI.PCoA.axis1.vs.body.size = ggplot(data = data.frame(x = x, y = y)) +
    geom_point(aes(x,y)) +
    scale_x_log10() +
    # scale_y_log10() +
    theme_bw() +
    ggtitle(LETTERS[3]) +
    # geom_hline(yintercept = 1, linetype="dashed") +
    # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          # axis.title.x=element_text(vjust = 45),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm")) +
    #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    labs(x=expression("Clade mean body size ("*mu*"m)"), y="PCoA axis 1") +
    geom_smooth(aes(x,y),method='lm',col="black") +
    annotate(geom="text", 
             x=(max(x)/min(x))^0.8*min(x),
             y=0.9*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=7)
  g8 = ggplotGrob(plot.norm.VI.PCoA.axis1.vs.body.size)
  
  x = as.vector(diversity)[selected_groups & diversity>div_threshold][-69]
  y = NormalizedVI_pcoa[[3]][,2]
  cor.test = cor.test(log10(x),y)
  plot.norm.VI.PCoA.axis2.vs.diversity = ggplot(data = data.frame(x = x, y = y)) +
    geom_point(aes(x,y)) +
    scale_x_log10() +
    # scale_y_log10() +
    theme_bw() +
    ggtitle(LETTERS[4]) +
    # geom_hline(yintercept = 1, linetype="dashed") +
    # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          # axis.title.x=element_text(vjust = 45),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm")) +
    #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    labs(x="Clade diversity", y="PCoA axis 2") +
    geom_smooth(aes(x,y),method='lm',col="black") +
    annotate(geom="text", 
             x=(max(x)/min(x))^0.8*min(x),
             y=0.2*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=7)
  g9 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.diversity)
  
  grid.arrange(grobs = list(g2,g4,g8,g9),
               # widths = rep(1,7),
               # layout_matrix = matrix(nrow=1,data=c(1,2)))
               layout_matrix = rbind(c(1, 2),c(3, 4)))
  dev.off()
  
  ##########################
  # PCoA axis 1 vs. 2 only #
  ##########################
  pdf(paste0(figure_folder,"/Normalized_PCoA.axes1.vs.2_",div_threshold,"plusOTUs",surf_DCM_insert,noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=6.9,height=6.7)
  plot(ggplotGrob(plot.VI.PCoA.biplot + theme(plot.title = element_blank())))
  dev.off()
}

# Old fig. 3 - lat. sym.
{
  lat_JSD.result = readRDS(paste0(results_folder,"/Lat_JSD_sigma2.10",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
  lat_JSD = lat_JSD.result[[1]]
  lat_I.result = readRDS(paste0(results_folder,"/Lat_I_sigma2.25",noArcticNoBiomark_insert,noLagoon_insert,".rds"))
  lat_I = lat_I.result[[1]]
  
  i_case = 1
  surf_DCM_insert = "_SUR"
  library(gridExtra)
  
  ##########
  plot.lat.JSD.range = ggplot(data = data.frame(x = sigma2.range, y = 1 - lat_JSD.allTaxa.range[,i_case])) +
    # geom_point(aes(x,y1)) +
    geom_point(aes(x,y)) +
    # scale_x_log10() +
    # scale_y_log10() +
    theme_bw() +
    # ggtitle(LETTERS[6]) +
    # geom_hline(yintercept = 1, linetype="dashed") +
    # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          # axis.title.x=element_text(vjust = 45),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm")) +
    #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    labs(x="Sigma2", y="Latitudinal symmetry - 1-JSD")
  
  pdf(paste0(figure_folder,"/Lat.sym.JSD.allTaxa_sigma2.range_",div_threshold,"plusOTUs",surf_DCM_insert,noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=6.9,height=6.7)
  plot(ggplotGrob(plot.lat.JSD.range))
  dev.off()
  
  plot.lat.I.range = ggplot(data = data.frame(x = sigma2.range, y = lat_I.allTaxa.range[,i_case])) +
    # geom_point(aes(x,y1)) +
    geom_point(aes(x,y)) +
    # scale_x_log10() +
    # scale_y_log10() +
    theme_bw() +
    # ggtitle(LETTERS[6]) +
    # geom_hline(yintercept = 1, linetype="dashed") +
    # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          # axis.title.x=element_text(vjust = 45),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm")) +
    #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    labs(x="Sigma2", y="Latitudinal symmetry - I")
  
  pdf(paste0(figure_folder,"/Lat.sym.I.allTaxa_sigma2.range_",div_threshold,"plusOTUs",surf_DCM_insert,noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=6.9,height=6.7)
  plot(ggplotGrob(plot.lat.I.range))
  dev.off()
  
  ###########
  x = lat_I[selected_groups & diversity>div_threshold,i_case][-69]
  y = 1 - lat_JSD[selected_groups & diversity>div_threshold,i_case][-69]
  cor.test = cor.test(x,y)
  plot.lat.JSD.vs.I = ggplot(data = data.frame(x = x, y = y)) +
    # geom_point(aes(x,y1)) +
    geom_point(aes(x,y)) +
    # scale_x_log10() +
    # scale_y_log10() +
    theme_bw() +
    # ggtitle(LETTERS[6]) +
    # geom_hline(yintercept = 1, linetype="dashed") +
    # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          # axis.title.x=element_text(vjust = 45),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm")) +
    #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    labs(x="Latitudinal symmetry - I", y="Latitudinal symmetry - 1-JSD") +
    geom_smooth(aes(x,y),method='lm',col="black") +
    annotate(geom="text", 
             x=0.8*(max(x)-min(x))+min(x),
             y=0.2*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=7)
  pdf(paste0(figure_folder,"/Lat.sym_JSD.vs.I_",div_threshold,"plusOTUs",surf_DCM_insert,noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=6.9,height=6.7)
  plot(ggplotGrob(plot.lat.JSD.vs.I))
  dev.off()
  
  #########
  x = 1 - lat_JSD[selected_groups & diversity>div_threshold,i_case][-69]
  y = NormalizedVI_pcoa[[3]][,2]
  cor.test = cor.test(x,y)
  plot.norm.VI.PCoA.axis2.vs.lat.JSD = ggplot(data = data.frame(x = x, y = y)) +
    # geom_point(aes(x,y1)) +
    geom_point(aes(x,y)) +
    # scale_x_log10() +
    # scale_y_log10() +
    theme_bw() +
    # ggtitle(LETTERS[6]) +
    # geom_hline(yintercept = 1, linetype="dashed") +
    # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          # axis.title.x=element_text(vjust = 45),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm")) +
    #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    labs(x=expression("Latitudinal symmetry"), y="PCoA axis 2") +
    geom_smooth(aes(x,y),method='lm',col="black") +
    annotate(geom="text", 
             x=0.8*(max(x)-min(x))+min(x),
             y=0.2*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=7)
  
  pdf(paste0(figure_folder,"/Normalized_PCoA.axes2.vs.lat.sym.JSD_sigma2.",sigma2,"_",div_threshold,"plusOTUs",surf_DCM_insert,noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=6.9,height=6.7)
  plot(ggplotGrob(plot.norm.VI.PCoA.axis2.vs.lat.JSD))
  dev.off()
  
  #########
  x = lat_I[selected_groups & diversity>div_threshold,i_case][-69]
  y = NormalizedVI_pcoa[[3]][,2]
  cor.test = cor.test(x,y)
  plot.norm.VI.PCoA.axis2.vs.lat.I = ggplot(data = data.frame(x = x, y = y)) +
    # geom_point(aes(x,y1)) +
    geom_point(aes(x,y)) +
    # scale_x_log10() +
    # scale_y_log10() +
    theme_bw() +
    # ggtitle(LETTERS[6]) +
    # geom_hline(yintercept = 1, linetype="dashed") +
    # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          # axis.title.x=element_text(vjust = 45),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm")) +
    #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    labs(x=expression("Latitudinal symmetry"), y="PCoA axis 2") +
    geom_smooth(aes(x,y),method='lm',col="black") +
    annotate(geom="text", 
             x=0.8*(max(x)-min(x))+min(x),
             y=0.1*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=7)
  
  pdf(paste0(figure_folder,"/Normalized_PCoA.axis2.vs.lat.sym.I_sigma2.",sigma2,"_",div_threshold,"plusOTUs",surf_DCM_insert,noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=6.9,height=6.7)
  plot(ggplotGrob(plot.norm.VI.PCoA.axis2.vs.lat.I))
  dev.off()
  
  #########
  x = lat_I[selected_groups & diversity>div_threshold,i_case][-69]/I_square.observed_w.mean[selected_groups & diversity>div_threshold,i_case][-69]
  y = NormalizedVI_pcoa[[3]][,1]
  cor.test = cor.test(x,y)
  plot.norm.VI.PCoA.axis1.vs.norm.lat.I = ggplot(data = data.frame(x = x, y = y)) +
    # geom_point(aes(x,y1)) +
    geom_point(aes(x,y)) +
    # scale_x_log10() +
    # scale_y_log10() +
    theme_bw() +
    # ggtitle(LETTERS[6]) +
    # geom_hline(yintercept = 1, linetype="dashed") +
    # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          # axis.title.x=element_text(vjust = 45),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm")) +
    #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    labs(x="Latitudinal symmetry over spatial str.", y="PCoA axis 1") +
    geom_smooth(aes(x,y),method='lm',col="black") +
    annotate(geom="text", 
             x=0.8*(max(x)-min(x))+min(x),
             y=0.1*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=7)
  
  pdf(paste0(figure_folder,"/Normalized_PCoA.axis1.vs.lat.sym.I.over.autocorr_sigma2.",sigma2,"_",div_threshold,"plusOTUs",surf_DCM_insert,noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=6.9,height=6.7)
  plot(ggplotGrob(plot.norm.VI.PCoA.axis1.vs.norm.lat.I))
  dev.off()
}

# Old fig. 3 - basin
{
  i_case = 1
  surf_DCM_insert = "_SUR"
  library(gridExtra)
  
  x = basin_I[selected_groups & diversity>div_threshold,i_case][-69]
  y = NormalizedVI_pcoa[[3]][,1]
  cor.test = cor.test(x,y)
  plot.norm.VI.PCoA.axis1.vs.basin.I = ggplot(data = data.frame(x = x, y = y)) +
    # geom_point(aes(x,y1)) +
    geom_point(aes(x,y)) +
    # scale_x_log10() +
    # scale_y_log10() +
    theme_bw() +
    # ggtitle(LETTERS[6]) +
    # geom_hline(yintercept = 1, linetype="dashed") +
    # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          # axis.title.x=element_text(vjust = 45),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm")) +
    #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    labs(x=expression("Basin-scale similarity"), y="PCoA axis 1") +
    geom_smooth(aes(x,y),method='lm',col="black") +
    annotate(geom="text", 
             x=0.8*(max(x)-min(x))+min(x),
             y=0.1*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=7)
  
  pdf(paste0(figure_folder,"/Normalized_PCoA.axis1.vs.basin.I_",div_threshold,"plusOTUs",surf_DCM_insert,noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=6.9,height=6.7)
  plot(ggplotGrob(plot.norm.VI.PCoA.axis1.vs.basin.I))
  dev.off()
  
  #########
  div_threshold = 100
  
  x = basin_I[selected_groups & diversity>div_threshold,i_case][-69]/I_square.observed_w.mean[selected_groups & diversity>div_threshold,i_case][-69]
  y = NormalizedVI_pcoa[[3]][,1]
  cor.test = cor.test(x,y)
  plot.norm.VI.PCoA.axis1.vs.basin.I = ggplot(data = data.frame(x = x, y = y)) +
    # geom_point(aes(x,y1)) +
    geom_point(aes(x,y)) +
    # scale_x_log10() +
    # scale_y_log10() +
    theme_bw() +
    # ggtitle(LETTERS[6]) +
    # geom_hline(yintercept = 1, linetype="dashed") +
    # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          # axis.title.x=element_text(vjust = 45),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm")) +
    #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    labs(x="Basin-scale similarity over spatial str.", y="PCoA axis 1") +
    geom_smooth(aes(x,y),method='lm',col="black") +
    annotate(geom="text", 
             x=0.8*(max(x)-min(x))+min(x),
             y=0.1*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=7)
  
  pdf(paste0(figure_folder,"/Normalized_PCoA.axis1.vs.basin.I.over.autocorr_",div_threshold,"plusOTUs",surf_DCM_insert,noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=6.9,height=6.7)
  plot(ggplotGrob(plot.norm.VI.PCoA.axis1.vs.basin.I))
  dev.off()
  
  #########
  div_threshold = 1000
  
  x = basin_I[selected_groups & diversity>div_threshold,i_case][-69]
  y = lat_I[selected_groups & diversity>div_threshold,i_case][-69]
  cor.test = cor.test(x,y)
  plot.lat.I.vs.basin.I = ggplot(data = data.frame(x = x, y = y)) +
    # geom_point(aes(x,y1)) +
    geom_point(aes(x,y)) +
    # scale_x_log10() +
    # scale_y_log10() +
    theme_bw() +
    # ggtitle(LETTERS[6]) +
    # geom_hline(yintercept = 1, linetype="dashed") +
    # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          # axis.title.x=element_text(vjust = 45),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm")) +
    #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    labs(x="Basin-scale similarity", y="Latitudinal symmetry") +
    geom_smooth(aes(x,y),method='lm',col="black") +
    annotate(geom="text", 
             x=0.8*(max(x)-min(x))+min(x),
             y=0.1*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=7)
  
  pdf(paste0(figure_folder,"/Lat.I.vs.basin.I_",div_threshold,"plusOTUs",surf_DCM_insert,noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=6.9,height=6.7)
  plot(ggplotGrob(plot.lat.I.vs.basin.I))
  dev.off()
}
 
# Fig. 2 - PCoA axes biplot  
{
  plot.VI.PCoA.biplot = ggplot(data = data.frame(x = NormalizedVI_pcoa[[3]][,1], y = NormalizedVI_pcoa[[3]][,2], 
                                                 col = log10(size_relativeAbund)[selected_groups & diversity>div_threshold],
                                                 size = log10(as.vector(diversity))[selected_groups & diversity>div_threshold])) + 
    geom_point(aes(x,y, colour = col, size = size), na.rm=T) +
    theme_bw() +
    # ggtitle(LETTERS[1]) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    # ggtitle("Norm. V.I. between clades' biogeographies") +
    theme(axis.title=element_text(size=22),
          axis.text=element_text(size=22),
          # axis.text.x=element_text(hjust=seq(from=0.1,to=0.9,length.out=6)),
          # axis.text.y=element_text(margin=unit(c(-1,-1,-1,-1), "cm")),
          # axis.ticks.length=unit(-0.25, "cm"),
          text = element_text(size=22),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm"),
          # plot.margin=unit(c(5,5,5,5),"mm"),
          # legend.position=c(0.8,0.2),
          # legend.position = "none",
          legend.position = "bottom",
          legend.direction = "vertical",
          # legend.spacing.x = unit(0.1,"npc"),
          # legend.spacing.y = unit(-0.1,"npc"),
          legend.text=element_text(size=20),
          legend.title=element_text(size=22)) +
    labs(x="PCoA axis 1", y="PCoA axis 2") +
    # geom_vline(xintercept = prop_within_OTU_null, linetype = "dashed") +
    scale_colour_gradient2(guide = "colourbar", name = expression("Body size ("*mu*"m)"), low = "blue", high = "red", mid = "grey",
                           midpoint = median(log10(size_relativeAbund)[selected_groups & diversity>div_threshold]),
                           breaks = log10(c(30,100,700)), labels = c(30,100,700)) +
    scale_size(range = c(2, 5), name = "Diversity (#OTUs)", 
               breaks = log10(c(101,1000,70000)), 
               labels = c("100","1,000","70,000")) +
    guides(size = guide_legend(title.position="top"),
           colour = guide_colorbar(barwidth = 12, barheight = 0.8, title.position="top"))
  
  pdf(paste0(figure_folder,"/Fig2_PCoA.axes.only.with.vert.legend_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.VI.PCoA.biplot)
  dev.off()
  
  ########
  plot.VI.consistent.PCoA.biplot = ggplot(data = data.frame(x = NormalizedVI.consistent_pcoa[[3]][,1], y = NormalizedVI.consistent_pcoa[[3]][,2], 
                                                 col = log10(size_relativeAbund)[selected_groups & diversity>div_threshold],
                                                 size = log10(as.vector(diversity))[selected_groups & diversity>div_threshold])) + 
    geom_point(aes(x,y, colour = col, size = size), na.rm=T) +
    theme_bw() +
    # ggtitle(LETTERS[1]) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    # ggtitle("Norm. V.I. between clades' biogeographies") +
    theme(axis.title=element_text(size=22),
          axis.text=element_text(size=22),
          # axis.text.x=element_text(hjust=seq(from=0.1,to=0.9,length.out=6)),
          # axis.text.y=element_text(margin=unit(c(-1,-1,-1,-1), "cm")),
          # axis.ticks.length=unit(-0.25, "cm"),
          text = element_text(size=22),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,12,0.5,3),"mm"),
          # plot.margin=unit(c(5,5,5,5),"mm"),
          # legend.position=c(0.8,0.2),
          # legend.position = "none",
          legend.position = "bottom",
          legend.direction = "vertical",
          # legend.spacing.x = unit(0.1,"npc"),
          # legend.spacing.y = unit(-0.1,"npc"),
          legend.text=element_text(size=20),
          legend.title=element_text(size=22)) +
    labs(x="PCoA axis 1", y="PCoA axis 2") +
    # geom_vline(xintercept = prop_within_OTU_null, linetype = "dashed") +
    scale_colour_gradient2(guide = "colourbar", name = expression("Body size ("*mu*"m)"), low = "blue", high = "red", mid = "grey",
                           midpoint = median(log10(size_relativeAbund)[selected_groups & diversity>div_threshold]),
                           breaks = log10(c(30,100,700)), labels = c(30,100,700)) +
    scale_size(range = c(2, 5), name = "Diversity (#OTUs)", 
               breaks = log10(c(101,1000,70000)), 
               labels = c("100","1,000","70,000")) +
    guides(size = guide_legend(title.position="top"),
           colour = guide_colorbar(barwidth = 12, barheight = 0.8, title.position="top"))
  
  pdf(paste0(figure_folder,"/Fig2_VI.consistent_PCoA.axes.only.with.vert.legend_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.VI.consistent.PCoA.biplot)
  dev.off()
}  
 
# Fig. 2 - Clades map plots:   
{
  data_Federico = 1
  data_V4 = 0
  data_psbO = 0
  
  variable_groups = T
  variable_reals = F
  variable_Ks = F
  variable_alpha.delta = F
  
  if (data_Federico)
  {
    marker = "18S_V9"
    short_marker = ""
  } else if (data_V4)
  {
    marker = "18S_V4"
    short_marker = "V4_"
  } else if (data_psbO)
  {
    marker = "psbO"
    short_marker = "psbO_"
  }
  
  figure_folder.marker = paste0(figure_folder,"/",substr(short_marker,1,nchar(short_marker)-1))
  
  # group_vect = c("Bacillariophyta","Collodaria","Acantharea","Arthropoda","Diplonemida","Dinophyceae")
  # number of topics to be represented (needs to be < nb_topics-1, or equal to nb_topics):
  # nb_topics1_vect = c(11,13,8,10,9,12)
  
  # group_vect = c("Chordata","Vannellida","MAST-3,_12","RAD-C")
  # V9 groups:
  if (data_Federico)
  {
    if (variable_groups)
    {
    # group_vect = c("Chordata","Diplonemida","Arthropoda","Bacillariophyta","Dinophyceae","Collodaria","Haptophyta","MAST-3,_12")
      group_vect = c("Diplonemida","Arthropoda","Chordata","Mamiellophyceae","Collodaria","Bacillariophyta","MAST-3,_12","Dinophyceae","MALV-II","Haptophyta")
    # group_vect = taxo_groups[selected_groups][sort.int(as.vector(I_square.observed_w.mean[selected_groups,1]),decreasing = T,index.return = T)$ix]
    # group_vect = taxo_groups[selected_groups][sort.int(as.vector(NormalizedVI_pcoa[[3]][,1]),decreasing = T,index.return = T)$ix]
    # group_vect = taxo_groups[selected_groups][sort.int(as.vector(NormalizedVI_pcoa[[3]][,2]),decreasing = T,index.return = T)$ix]
    # group_vect = taxo_groups[selected_groups & diversity > 1000][sort.int(as.vector(NormalizedVI_pcoa[[3]][diversity[selected_groups] > 1000, 2]), decreasing = T,index.return = T)$ix]
    # group_vect = taxo_groups_V4[Normalized_VI_V9.V4 != 0][sort.int(Normalized_VI_V9.V4[Normalized_VI_V9.V4 != 0], decreasing = F,index.return = T)$ix]
    # group_vect = taxo_groups_psbO[Normalized_VI_V9.psbO != 0][sort.int(Normalized_VI_V9.psbO[Normalized_VI_V9.psbO != 0], decreasing = F,index.return = T)$ix]
      # group_vect = c("Diplonemida","Mamiellophyceae","Collodaria","MALV-II")
    } else if (variable_reals)
    {
    # real.number.vect = c(1,37,64)
    } else if (variable_Ks)
    {
    # K.vect = c(3,5,8,12)
    } else if (variable_alpha.delta)
    {
      alpha.delta.vect = rbind(c(50/nb_topics,0.1),
                             c(0.1,0.1),
                             c(0.05/nb_topics,0.1), # 1st_best_realization
                             c(1,1),
                             c(0.1,1), # Not available
                             c(0.05/nb_topics,9.4-nb_topics/5)) #1st_best_realization
    }
    # } else if (k_by_k)
    # {
    #   group = "AllTaxa"
    #   nb_topics = 16
    # }
  # psbO groups:
  } else if (data_psbO)
  {
    # group_vect = c("Bacillariophyta","Haptophyta")
    group_vect = taxo_groups_psbO[Normalized_VI_V9.psbO != 0][sort.int(Normalized_VI_V9.psbO[Normalized_VI_V9.psbO != 0], decreasing = F,index.return = T)$ix]
  # V4 groups:
  } else if (data_V4)
    # group_vect = c("Bacillariophyta","Haptophyta","Dinophyceae","Collodaria","MAST-3,_12")
    group_vect = taxo_groups_V4[Normalized_VI_V9.V4 != 0][sort.int(Normalized_VI_V9.V4[Normalized_VI_V9.V4 != 0], decreasing = F,index.return = T)$ix]
  # group_vect = c("Phaeodaria","Annelida","Mollusca","Bacillariophyta","Dinophyceae","Chlorophyceae")
  # group_vect = "AllTaxa"
  
  monochrom = F # For groups other than AllTaxa, shows assemblages in variations of a single color
  manual = F # Manual setting of colors; if both monochrom = F and manual = F, attempts automatic setting of colors based on latitude
  dominant = F # Show one assemablage per station only (the dominant one)
  k_by_k = T # Plot each assemblage individually on a color gradient
  
  # Only relevant for monochrom = T:
  fixed_nb_colors = T
  new_colors = F
  V9_Ks = F
  SUR_sorted_assemblages = T
  
  if (V9_Ks)
  {
    V9_Ks_insert = "V9.Ks_"
  } else
    V9_Ks_insert = ""
  
  if (SUR_sorted_assemblages)
  {
    SUR.sorted_insert = "SUR.sorted_"
  } else
    SUR.sorted_insert = ""
  
  if (dominant)
  {
    dominant_insert = "_dominant.only"
  } else
  {
    dominant_insert = ""
  }
  
  if (k_by_k)
  {
    k.by.k_insert = "_k.by.k"
  } else
    k.by.k_insert = ""
  
  # Colours of main taxonomic groups:
  # colours = c("burlywood1",#"darkgoldenrod1",
  #             colorRampPalette(c("darkgreen","chartreuse2"),space = "Lab")(4),
  #             #"cornsilk1",
  #             "mintcream",
  #             colorRampPalette(c("darkblue","darkturquoise"),space = "Lab")(4),
  #             colorRampPalette(c("firebrick2","darkgoldenrod1"),space = "Lab")(5),
  #             colorRampPalette(c("darkviolet","deeppink"),space = "Lab")(3),
  #             # "burlywood1",
  #             "darkseagreen",
  #             "grey")
  
  if (monochrom)
  {
    ground_colors = c("burlywood1",
                      colorRampPalette(c("darkgreen","chartreuse2"),space = "Lab")(4)[1:2],
                      colorRampPalette(c("darkblue","darkturquoise"),space = "Lab")(4)[1],
                      colorRampPalette(c("firebrick2","darkgoldenrod1"),space = "Lab")(5)[1],
                      colorRampPalette(c("darkviolet","deeppink"),space = "Lab")(3)[1:2],
                      "darkseagreen",
                      "mintcream",
                      colorRampPalette(c("firebrick2","darkgoldenrod1"),space = "Lab")(5)[3])
    ground_colors = desaturate(ground_colors, amount = 0.2)
    names(ground_colors) = c("Diplonemida","Arthropoda","Chordata","Collodaria","Dinophyceae","Bacillariophyta","MAST-3,_12","Haptophyta","Mamiellophyceae","MALV-II")
    
    if (new_colors)
    {
      groups_col.pal = list(c(ground_colors[1],"#794506","#8f6024","#ee902b","#957350","#dda46e"),
                            c(ground_colors[2],"#039606","#26ee2d","#62b861","#456444","#a8f5a4","#0b3308"),
                            c(ground_colors[3],"#aaf896","#33552a","#30cb10","#89af83"),
                            c(ground_colors[4],"#776de8","#464c94","#0040ff","#08399c","#98c3ff","#96b1cc","#1a76b8"),
                            # colorRampPalette(c(ground_colors[5],"firebrick1","coral","chocolate2","coral4","darkred"),space = "Lab"),
                            # colorRampPalette(c("tan1",ground_colors[5],"darkred"),space = "Lab"),
                            # colorRampPalette(c(ground_colors[5],"sienna1","sandybrown","peru","coral4","darkred"),space = "Lab"),
                            c(ground_colors[5],"#8c3429","#f55424","#eb7d42","#9e6a3c","#a66100","#cf991d"),
                            c(ground_colors[6],"#5c0980","#c118f0","#d670e6","#a33ea8","#820a7e","#e517cd","#e874d1"),
                            c(ground_colors[7],"#800962","#e815a2","#d16da7","#a83e73","#820a3c","#e5175c","#e85d80"),
                            c(ground_colors[8],"#22d48d","#85c7af","#3d997a","#09573f","#3d9940","#105714"),
                            c(ground_colors[9],"#a6fcd1","#36f796","#109d57","#5f9178","#41d2ae","#149474","#7ba"),
                            c(ground_colors[10],"#d28c67","#f2af8c","#e6631e","#bb480c","#eb974e","#d29967","#f2bc8c"))
    } else
      groups_col.pal = list(colorRampPalette(c(ground_colors[1],"darkgoldenrod2","#dda46e","#ee902b","chocolate","chocolate4","#bfa382","#eba557"),space = "Lab"),
                            ########
                            # colorRampPalette(c(ground_colors[2],"olivedrab3","darkolivegreen1","chartreuse2","aquamarine","lightseagreen"),space = "Lab"),
                            # colorRampPalette(c(ground_colors[2],"darkolivegreen1","olivedrab3","#62b861","#456444","#a8f5a4","#039606"),space = "Lab"),
                            colorRampPalette(c(ground_colors[2],"#bced68","#8fbf2e","#61b860","#5d945b","#9be399","#69d934"),space = "Lab"),
                            ########
                            # colorRampPalette(c(ground_colors[3],"springgreen3","darkolivegreen1","olivedrab3","aquamarine"),space = "Lab"),
                            # colorRampPalette(c(ground_colors[3],"darkolivegreen1","springgreen3","olivedrab3","#30cb10"),space = "Lab"),
                            colorRampPalette(c(ground_colors[3],"#c7f47b","#28b870","#92c133","#6ddb57"),space = "Lab"),
                            ########
                            # colorRampPalette(c(ground_colors[4],"#0941b3","#7ce6f3","slateblue","#5079f1","deepskyblue","#98c3ff","darkcyan","paleturquoise3"),space = "Lab"),
                            colorRampPalette(c(ground_colors[4],"#1348b4","#7ba5f9","#7e7ed8","#5079f1","#4985c1","#b0caff","#001ba3","#90a5d0",
                                               "#41417c","#616eb3","#4242b3","#b8b8ff","#2b2bf7"),space = "Lab"),
                            # colorRampPalette(c(ground_colors[4],"#08399c","cadetblue1","slateblue","royalblue","deepskyblue","#98c3ff","darkcyan","paleturquoise3"),space = "Lab"),
                            ########
                            # colorRampPalette(c(ground_colors[5],"darkorange1","goldenrod1","peru","coral4","darkred"),space = "Lab"),
                            # colorRampPalette(c(ground_colors[5],"goldenrod1","darkorange1","#9e2108","#b05b07","#cf991d","coral4","#9e6a3c"),space = "Lab"),
                            colorRampPalette(c(ground_colors[5],"goldenrod1","darkorange1","#9e2108","#b05b07","#cf991d","coral4","#9e6a3c",
                                               "#d58739","#e0b96c","#f7a654","#f48662"),space = "Lab"),
                            ########
                            # colorRampPalette(c(ground_colors[6],"#a33ea8","darkorchid1","plum1","palevioletred","hotpink","deeppink3"),space = "Lab"),
                            colorRampPalette(c(ground_colors[6],"#EE95C7","#B53CF2","#a33ea8","#e517cd","#CC1076","#EB94FF","#F65DA9"),space = "Lab"),
                            ########
                            # colorRampPalette(c(ground_colors[7],"darkorchid1","plum1","pink3","hotpink","deeppink3"),space = "Lab"),
                            # colorRampPalette(c(ground_colors[7],"#BF3DF1","#EE9AFF","#E58ECB","#E84E89","#E4149B","#CC1076","#DB6696"),space = "Lab"),
                            # colorRampPalette(c(ground_colors[7],"#BF3DF1","#EE9AFF","#E58ECB","#d75d8c","#d72d9b","#a11762","#e889b0"),space = "Lab"),
                            colorRampPalette(c(ground_colors[7],"#BF3DF1","#EE9AFF","#f57dd2","#d75d8c","#eb36ab","#b81a6f","#e889b0"),space = "Lab"),
                            ########
                            # colorRampPalette(c(ground_colors[8],"palegreen2","aquamarine3","cyan4","paleturquoise4","cadetblue2"),space = "Lab"))
                            colorRampPalette(c(ground_colors[8],"palegreen2","#3d997a","#85c7af","#147a8f","paleturquoise4","#22d48d"),space = "Lab"),
                            colorRampPalette(c(ground_colors[9],"#83d8ad","#37fa98","#c0fcde","#5f9178","#2bb894","#128266","#6cb2a0")[c(3,7,1,2,5,6,4,8)],space = "Lab"),
                            colorRampPalette(c(ground_colors[10],"#d28c67","#f2af8c","#e6631e","#bb480c","#eb974e","#d29967","#f2bc8c"),space = "Lab"))
    names(groups_col.pal) = c("Diplonemida","Arthropoda","Chordata","Collodaria","Dinophyceae","Bacillariophyta","MAST-3,_12","Haptophyta","Mamiellophyceae","MALV-II")                       
    
    if (fixed_nb_colors)
    {  
      # groups_nb_colors =  c(5,6,5,7,6,7,6,6)      
      groups_nb_colors =    c(8,7,5,14,12,8,8,7,8,8)
      names(groups_nb_colors) = c("Diplonemida","Arthropoda","Chordata","Collodaria","Dinophyceae","Bacillariophyta","MAST-3,_12","Haptophyta","Mamiellophyceae","MALV-II")
    }
  } else if (manual)
  {
    color.pal1 = colorRampPalette(c("#7F0000","red","darkorange1","darkgoldenrod1","yellow"),space = "Lab")
    color.pal2 = colorRampPalette(c("#7FFF7F","darkgreen"),space = "Lab")
    color.pal3 = colorRampPalette(c("azure4","antiquewhite3","aliceblue"),space = "Lab")
    color.pal4 = colorRampPalette(c("#00007F","#007FFF","cyan"),space = "Lab")
  } else if (k_by_k)
  {
    # col.pal = colorRampPalette(c("darkblue","firebrick2"),space = "Lab")
    library(paletteer)
  } else
  {
    # color.pal = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),space = "Lab")
    color.pal = colorRampPalette(c("#00007F", "#007FFF", "cyan", "darkgreen", "yellow", "darkorange1", "red", "#7F0000"),space = "Lab")
    # color.pal = colorRampPalette(c("#00007F", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),space = "Lab")
    color.pal.warm = colorRampPalette(c("#7F0000","red","darkorange1","darkgoldenrod1","yellow"),space = "Lab")
    color.pal.green = colorRampPalette(c("#7FFF7F","darkgreen"),space = "Lab")
    color.pal.blue = colorRampPalette(c("cyan","#007FFF","blue","#00007F"),space = "Lab")
    # color.pal = colorRampPalette(c("#00007F", "#007FFF", "cyan", "#7FFF7F", "darkgreen", "yellow", "darkgoldenrod1", "darkorange1", "red", "#7F0000"),space = "Lab")
    # color.pal = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),space = "Lab")
  }
  
  # wheel <- function(col, radius = 1, ...)
  #   pie(rep(1, length(col)), col = col, radius = radius, ...)
  # wheel(col)
  
  # Loading batymetric data for the specified range of longitudes and latitudes in degrees, whith resolution in minutes
  bat = getNOAA.bathy(-180, 180, -90, 90, res = 50, keep=T, 
                      # path = paste0(data_folder,"/marmap_ocean_04.21"))
                      # path = paste0(data_folder,"/marmap_ocean_06.21"))
                      path = paste0(data_folder,"/marmap_ocean_08.19"))
  
  tropic_bound = 23
  arctic_bound = 50
  
  # group = "AllTaxa"
  # group = "Bacillariophyta"
  # group = "Collodaria"
  # group = "Arthropoda"
  # group = "Diplonemida"
  # group = "Chordata"
  # group = "Mollusca"
  # group = "Acantharea"
  # group = "Dinophyceae"
  # group = "Bolidophyceae"
  
  #######
  data.folder_name = paste0(data_folder,"/18S_V9_TARA_CompleteSizeRange_byStationByDepth_AllTAxa",noArcticNoBiomark_insert,noLagoon_insert)
  load(paste0(data.folder_name,"/coord.Rdata"))
  stations_depths = as.data.frame(matrix(nrow=nrow(coord),ncol=2,dimnames = list(rownames(coord),c("Station","Depth"))))
  for (station_depth_index in 1:nrow(coord))
  {
    stations_depths[station_depth_index,1:2] = c(strsplit(rownames(coord)[station_depth_index],split=" ",fixed=T)[[1]][1],
                                                 strsplit(rownames(coord)[station_depth_index],split=" ",fixed=T)[[1]][2])
  }
  pie_positions.SUR = read.table(paste0(results_folder,"/Surf_pies_modified.txt"),sep=",",header=F)
  # pie_positions.SUR[c(55,57) - 1,] = pie_positions.SUR[c(57,55) - 1,]
  pie_positions.SUR[55:57 - 1,] = pie_positions.SUR[c(56,57,55) - 1,]
  pie_positions.SUR[56 - 1,] = pie_positions.SUR[56 - 1,] + c(1.5,0)
  pie_positions.SUR[45:46 - 1,] = pie_positions.SUR[46:45 - 1,]
  pie_positions.SUR[51 - 1,] = pie_positions.SUR[51 - 1,] + c(-2,3)
  rownames(pie_positions.SUR) = rownames(coord[stations_depths[,2] == "SUR",])[-44]
  #########
  
  if (!data_Federico)
  {
    data.folder_name = paste0(data_folder,"/",marker,"_TARA_CompleteSizeRange_byStationByDepth_AllTAxa",noArcticNoBiomark_insert,noLagoon_insert)
    load(paste0(data.folder_name,"/coord.Rdata"))
    stations_depths = as.data.frame(matrix(nrow=nrow(coord),ncol=2,dimnames = list(rownames(coord),c("Station","Depth"))))
    for (station_depth_index in 1:nrow(coord))
    {
      stations_depths[station_depth_index,1:2] = c(strsplit(rownames(coord)[station_depth_index],split=" ",fixed=T)[[1]][1],
                                                   strsplit(rownames(coord)[station_depth_index],split=" ",fixed=T)[[1]][2])
    }
  }
  
  # library(gridBase)
  
  # plot.map.group = list()
  if (variable_groups)
  {
    if (!k_by_k)
      page_length = ifelse(length(group_vect) %/% 4 == 0,length(group_vect),4)
    else
      # page_length = ifelse(nb_topics %/% 6 == 0,length(nb_topics),6)
      page_length = 6
  } else if (variable_reals)
  {
    page_length = ifelse(length(real.number.vect) %/% 4 == 0,length(real.number.vect),4)
  } else if (variable_Ks)
  {
    page_length = ifelse(length(K.vect) %/% 4 == 0,length(K.vect),4)
  } else if (variable_alpha.delta)
  {
    page_length = ifelse(length(alpha.delta.vect) %/% 4 == 0,length(alpha.delta.vect),4)
  }
  # pdf(paste0(figure_folder.marker,"/Fig2_",short_marker,V9_Ks_insert,if (monochrom) SUR.sorted_insert else "","maps.only.unif.bg",dominant_insert,"_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1_alt7.pdf"),
  # pdf(paste0(figure_folder.marker,"/FigS.all.selected.groups.decreasing.V9.similarity_",short_marker,V9_Ks_insert,if (monochrom) SUR.sorted_insert else "","maps.only.unif.bg",dominant_insert,"_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
  # pdf(paste0(figure_folder.marker,"/FigS.allTaxa_K16_real1-37-64_alpha0.1_delta0.1.pdf"),
  # pdf(paste0(figure_folder.marker,"/FigS.allTaxa_K3-20_real1_alpha0.1_delta0.1.pdf"),
  # pdf(paste0(figure_folder.marker,"/FigS.allTaxa_K3-5-8-12_real1_alpha0.1_delta0.1.pdf"),
  # pdf(paste0(figure_folder.marker,"/FigS.allTaxa_K16_real1_variable.alpha-delta.pdf"),
  # pdf(paste0(figure_folder.marker,"/FigS.allTaxa_K16_k.by.k.pdf"),
  pdf(paste0(figure_folder.marker,"/FigS_k_by_k_10groups.pdf"),
  # pdf(paste0(figure_folder.marker,"/FigS_Diplo-Mamiel-Collo-MALV.II.pdf"),
  # pdf(paste0(figure_folder.marker,"/FigS.all.selected.groups.decreasing.V4.similarity_",short_marker,V9_Ks_insert,if (monochrom) SUR.sorted_insert else "","maps.only.unif.bg",dominant_insert,"_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
  # pdf(paste0(figure_folder.marker,"/FigS.all.selected.groups.decreasing.Moran.I_",short_marker,V9_Ks_insert,if (monochrom) SUR.sorted_insert else "","maps.only.unif.bg",dominant_insert,"_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"), 
  # pdf(paste0(figure_folder.marker,"/FigS.all.selected.groups.decreasing.PCoA1_",short_marker,V9_Ks_insert,if (monochrom) SUR.sorted_insert else "","maps.only.unif.bg",dominant_insert,"_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"), 
  # pdf(paste0(figure_folder.marker,"/FigS.>1000.groups.decreasing.PCoA2_",short_marker,V9_Ks_insert,if (monochrom) SUR.sorted_insert else "","maps.only.unif.bg",dominant_insert,"_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      # width=7.5*2.2/1.2/2*1.9*3,height=12/3*4/1.2/2*page_length,onefile=T)
  width = if (k_by_k) 7.5*2.2/1.2/2*1.9*3 else 7.5*2.2/1.2/2*1.9,
  height = 12/3*4/1.2/2*page_length,
  onefile=T)
  # width=7.5*2.2/1.2/2*3,height=12/3*4/1.2/2
      # colormodel = "cmyk")
  if (k_by_k)
  {
    # layout_mat = matrix(rep(1:(page_length*3),each=2),nrow=page_length,byrow = T)
    # for (i in 1:nrow(layout_mat))
    #   rep(1:(page_length*3),each=2)
    layout(mat=matrix(rep(1:(page_length*3),each=2),nrow=page_length,byrow = T))
  } else
    layout(mat=matrix(rep(1:page_length,each=2),nrow=page_length,byrow = T))
  # i_group = 0
  for (group in group_vect)
  # for (real.number in real.number.vect)
  # for (nb_topics in K.vect)
  # for (i_param in 1:nrow(alpha.delta.vect))
  # for (k in 1:nb_topics)
  {
    if (variable_groups)
    {
      cat(group,"\n")
    } else if (variable_reals)
    {
      cat(real.number,"\n")
    } else if (variable_Ks)
    {
      cat(nb_topics,"\n")
    } else if (variable_alpha.delta)
    {
      alpha = alpha.delta.vect[i_param,1]
      delta = alpha.delta.vect[i_param,2]
      cat("alpha =",alpha,"- delta =",delta,"\n")
    }
      
    # i_group = i_group+1
      
    if (!variable_Ks && !variable_alpha.delta)
    {
      if (data_Federico || V9_Ks)
      {
        if (group != "AllTaxa")
          nb_topics = optimalK_prevalence.min.crossValid[taxo_groups == group]
        else
          nb_topics = optimalK_prevalence.min.crossValid.allTaxa
      } else if (data_psbO)
      {
        nb_topics = optimalK_prevalence.min.crossValid_psbO[taxo_groups_psbO == group]
      } else if (data_V4)
        nb_topics = optimalK_prevalence.min.crossValid_V4[taxo_groups_V4 == group]
    }
    
    # data.folder_name = paste0(data_folder_workspace3,"/",marker,"_TARA_CompleteSizeRange_byStationByDepth_",group,"_noLagoon")
    data.folder_name = paste0(data_folder,"/",marker,"_TARA_CompleteSizeRange_byStationByDepth_",group,"_noLagoon")
    if (variable_alpha.delta)
    {
      if (alpha == 0 && delta == 0)
        param.folder_name = paste0("Rtopicmodels_LDA_VEM_nb_topics",nb_topics,"_nb_real100_em_tol1e-06_var_tol1e-08_best_keep_occurrence")
      else if (alpha == 50/nb_topics && delta == 0.1)
        param.folder_name = paste0("Rtopicmodels_LDA_Gibbs_alpha50:K_delta0.1_nb_topics",nb_topics,"_nb_iter2000_nb_real100_occurrence")
      else if (alpha == 0.1 && delta == 0.1)
        param.folder_name = paste0("Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics",nb_topics,"_nb_iter2000_nb_real100_occurrence")
      else if (alpha == 0.05/nb_topics && delta == 0.1)
        param.folder_name = paste0("Rtopicmodels_LDA_Gibbs_alpha0.05:K_delta0.1_nb_topics",nb_topics,"_nb_iter2000_nb_real100_occurrence")
      else if (alpha == 1 && delta == 1)
        param.folder_name = paste0("Rtopicmodels_LDA_Gibbs_alpha1_delta1_nb_topics",nb_topics,"_nb_iter2000_nb_real100_occurrence")
      else if (alpha == 0.1 && delta == 1)
        param.folder_name = paste0("Rtopicmodels_LDA_Gibbs_alpha0.1_delta1_nb_topics",nb_topics,"_nb_iter2000_nb_real100_occurrence")
      else if (alpha == 0.05/nb_topics && delta == 9.4-nb_topics/5)
        param.folder_name = paste0("Rtopicmodels_LDA_Gibbs_alpha0.05:K_delta",9.4-nb_topics/5,"_nb_topics",nb_topics,"_nb_iter2000_nb_real100_occurrence")
    } else
      param.folder_name = paste0("Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics",nb_topics,"_nb_iter1000_nb_real100_meanPosteriorDistributedLlh_thin25_burnin2000_occurrence")
    
    if (variable_alpha.delta)
      save.file_name = paste0(data.folder_name,"/",param.folder_name,"/1st_best_realization/Spatial_topicmix_kriged.rds")
    else
      save.file_name = paste0(data.folder_name,"/",param.folder_name,"/1st_closestToMean_realization/Spatial_topicmix_kriged.rds")
    
    if (variable_reals)
    {
      if (real.number == real.number.vect[1])
      {
        # spatial_topicmix_kriged = readRDS(paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics",nb_topics,"_nb_iter250_nb_real100_best_thin2_burnin0_occurrence/1st_best_realization/Spatial_topicmix_kriged.rds"))
        spatial_topicmix_kriged = readRDS(save.file_name)
        # selecting the z.pred columns in all topics:
        documents = unlist(lapply(spatial_topicmix_kriged,function(g) g$z.pred))
        # setting one topic per column
        documents = matrix(documents,ncol=nb_topics,dimnames=list(rownames(spatial_topicmix_kriged[[1]]),paste0("assemblage",1:nb_topics)))
        rownames.documents = rownames(documents)
        documents_real1 = documents
        
        save.file.allReals_name = paste0(data.folder_name,"/",param.folder_name,"/",param.folder_name,".Rdata")
        load(save.file.allReals_name)
        Ordered_realizations = readRDS(paste0(data.folder_name,"/",param.folder_name,"/Ordered_realizations.rds"))
      }
      if (real.number != 1)
      {
        documents = Result[[Ordered_realizations$ix[real.number]]]@gamma
        rownames(documents) = rownames.documents
      } else
        documents = documents_real1
    } else if (file.exists(save.file_name))
    {
      spatial_topicmix_kriged = readRDS(save.file_name)
      # selecting the z.pred columns in all topics:
      documents = unlist(lapply(spatial_topicmix_kriged,function(g) g$z.pred))
      # setting one topic per column
      documents = matrix(documents,ncol=nb_topics,dimnames=list(rownames(spatial_topicmix_kriged[[1]]),paste0("assemblage",1:nb_topics)))
    } else if (variable_alpha.delta)
    {
      save.file.allReals_name = paste0(data.folder_name,"/",param.folder_name,"/",param.folder_name,".Rdata")
      
      local.data.folder_name = paste0(data_folder,"/",marker,"_TARA_CompleteSizeRange_byStationByDepth_",group,"_noLagoon")
      local.save.file_name = paste0(local.data.folder_name,"/",param.folder_name,"/documents_1st_real.rds")
      if (file.exists(local.save.file_name))
      {
        documents = readRDS(local.save.file_name)
        Ordered_realizations = readRDS(paste0(local.data.folder_name,"/",param.folder_name,"/Ordered_realizations.rds"))
      } else if (file.exists(save.file.allReals_name))
      {
        load(save.file.allReals_name)
        Ordered_realizations = readRDS(paste0(data.folder_name,"/",param.folder_name,"/Ordered_realizations.rds"))
        documents = Result[[Ordered_realizations$ix[1]]]@gamma
        rownames(documents) = rownames(coord)
        
        if (!dir.exists(paste0(local.data.folder_name,"/",param.folder_name,"/")))
          dir.create(paste0(local.data.folder_name,"/",param.folder_name,"/"))
        saveRDS(documents,file = local.save.file_name)
        saveRDS(Ordered_realizations, file = paste0(local.data.folder_name,"/",param.folder_name,"/Ordered_realizations.rds"))
      } else
      {
        cat("Missing file for alpha =",alpha,"- delta =",delta,"\n")
        next
      }
    } else
    {
      if (variable_Ks)
        cat("Missing file for",nb_topics,"assemblages\n")
      else if (variable_groups)
        cat("Missing file for",group,"\n")
      next
    }
    
    # Concatenating the spatial coord. with documents into spatial_topicmix_kriged_all_topics:
    # spatial_topicmix_kriged_all_topics = cbind(data.frame(x = spatial_topicmix_kriged[[1]]$x, y = spatial_topicmix_kriged[[1]]$y),documents)
    spatial_topicmix_kriged_all_topics = cbind(data.frame(x = coord$x[rownames(coord) %in% rownames(documents)],
                                                          y = coord$y[rownames(coord) %in% rownames(documents)]),documents)
    
    # # Storing all assemblages into a single dataframe:
    # spatial_topicmix_kriged_all_topics = data.frame(x = spatial_topicmix_kriged[[1]]$x, y = spatial_topicmix_kriged[[1]]$y)
    # rownames(spatial_topicmix_kriged_all_topics) = rownames(spatial_topicmix_kriged[[1]])
    # for (k in 1:nb_topics)
    # {
    #   spatial_topicmix_kriged[[k]]$z.pred[spatial_topicmix_kriged[[k]]$z.pred < 0.001] = 0
    #   spatial_topicmix_kriged_all_topics = cbind(spatial_topicmix_kriged_all_topics,spatial_topicmix_kriged[[k]]$z.pred)
    #   colnames(spatial_topicmix_kriged_all_topics)[2+k] = paste0("z.pred",k)
    # }
    
    # cat("\nPlotting ...")
    
    # plot.map.group[[i_group]] = plot.map(pie_positions.SUR,spatial_topicmix_kriged_all_topics,nb_topics,bat) 
    
    # Sorting topics by their abundance in surface samples:
    documents_SUR = documents[rownames(documents) %in% rownames(pie_positions.SUR),]
    sorted_SUR_assemblages = sort.int(colSums(documents_SUR),decreasing=T,index.return=T)
    
    blues = c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
    greys = c(grey(0.6), grey(0.93), grey(0.99))
    
    if (!k_by_k)
    {
      # plot(bat, image = T, land = T, lty = 0, bty = "n", xaxt = "n", yaxt = "n", ann = F, bpal = list(c(0, max(bat), greys[1]), c(min(bat), 0, blues))) #plot map without isobaths
      plot(bat, image = T, land = T, lty = 0, bty = "n", xaxt = "n", yaxt = "n", ann = F, bpal = list(c(0, max(bat), greys[1]), c(min(bat), 0, blues[3])))
      plot(bat, lwd = 0.2, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline
      # par(mar = c(0,0,0,0))
      # present_topics = apply(spatial_topicmix_kriged_all_topics[,3:(2+(nb_topics0))],1,function(g) which(g>0))
      # pure_stations = unlist(lapply(present_topics,function(g) length(g))) == 1
      # points(cbind(spatial_topicmix_kriged_all_topics$x[stations_depths[,2] == "SUR" & pure_stations],spatial_topicmix_kriged_all_topics$y[stations_depths[,2] == "SUR" & pure_stations]),
      #        pch = 21, cex=1.2, bg = col[unlist(present_topics[stations_depths[,2] == "SUR" & pure_stations])])
      # pie.colors.matrix = data.frame(matrix(nrow = length(which(stations_depths[rownames(coord) %in% rownames(documents),2] == "SUR")), ncol=nb_topics0, data = NA))
      # pie.slices.matrix = data.frame(matrix(nrow = length(which(stations_depths[rownames(coord) %in% rownames(documents),2] == "SUR")), ncol=nb_topics0, data = NA))
      
      if (variable_groups)
        title(group, cex.main = 2)
      else if (variable_Ks)
        title(paste(nb_topics,"assemblages"), cex.main = 2)
      else if (variable_alpha.delta)
        title(paste("alpha =",alpha,"- delta =",delta), cex.main = 2)
      # else if (k_by_k)
      #   title(paste("Assemblage",k), cex.main = 5)
      # title(paste(group,"-", nb_absolute_dominants[taxo_groups == group,1], "dominant assemblages"), cex.main = 2)
      # title(paste(group,"- V9 dissimilarity =", format(Normalized_VI_V9.V4[taxo_groups_V4 == group],digits=3)), cex.main = 2)
      # title(paste(group,"- V9 dissimilarity =", format(Normalized_VI_V9.psbO[taxo_groups_psbO == group],digits=3)), cex.main = 2)
      
      if (!dominant)
      {  
        if (monochrom)
        {
          col = vector(length = nb_topics, mode = "character")
          col.pal = groups_col.pal[[which(names(groups_col.pal) == group)]]
          if (fixed_nb_colors)
          {
            nb_colors = groups_nb_colors[names(groups_nb_colors) == group]
            if (SUR_sorted_assemblages)
            {
              col[sorted_SUR_assemblages$ix[1:nb_colors]] = if (new_colors) col.pal else col.pal(nb_colors)
              col[sorted_SUR_assemblages$ix[(nb_colors+1):nb_topics]] = rep("grey",nb_topics - nb_colors)
            } else
            {
              col[1:nb_colors] = if (new_colors) col.pal else col.pal(nb_colors)
              col[(nb_colors+1):nb_topics] = rep("grey",nb_topics - nb_colors)
            } 
            
            if (data_psbO)
            {
              if (V9_Ks)
              {
                if (group == "Bacillariophyta")
                {
                  col[1:8] = c(col[c(2,3,1)],
                               "grey",
                               col[c(6,5,7)],
                               col[4]) # 7 colors
                  # col[8:10] = c("red","blue","green")
                } else if (group == "Haptophyta")
                {
                  col[1:6] = c(col[c(2,1,3,5,4)],
                               "grey")
                } else if (group == "Dinophyceae")   
                  col = col
              } else
              {
                if (group == "Bacillariophyta")
                {
                  col[1:8] = col[c(3,2,1,6,8,5,7,4)]
                } else if (group == "Haptophyta")
                {
                  col[1:7] = c(col[c(2,1,4)],"grey",col[c(3,5,7)]) # 6
                } else if (group == "Dinophyceae")  
                  col = col
              }
            } else if (data_V4)
            {
              if (V9_Ks)
              {
                if (group == "Bacillariophyta")
                {
                  # col[1:7] = col[c(1,2,3,6,7,4,5)]
                  col[1:7] = c(col[c(1,2,3,6)],"grey",col[c(4,5)])
                } else if (group == "Haptophyta")
                {
                  col[1:7] = c(col[c(1,2)],"grey",col[c(3,4,5)],
                               col[6])
                } else if (group == "Dinophyceae")
                {
                  col[1:6] = c(col[c(1,2,5,3)],rep("grey",2))
                }
              } else
              {
                if (group == "Bacillariophyta")
                {
                  col[1:9] = col[c(1,3,2,6,4,7,5,8,9)]
                } else if (group == "Haptophyta")
                {
                  col[1:7] = c(col[c(1,2)],"grey",col[c(6,5,3,4)])
                } else if (group == "Dinophyceae")
                {
                  col[1:14] = c(col[c(1,2,3,6,10,4,5,8,12)],"grey",col[11],rep("grey",3))
                } else if (group == "Collodaria")
                {
                  col[1:13] = c(col[c(1,2,10,6,5)],rep("grey",3),col[c(3,9)],rep("grey",2),col[13])
                  # col[1:9] = c(col[c(1,6,7,2,5,4)],"grey","grey","grey") #(3),9,8
                } else if (group == "MAST-3,_12")
                  col[1:8] = col[1:8]
                # col = col
              }
            }
          } else
          {
            col[sorted_SUR_assemblages$ix[sorted_SUR_assemblages$x > 1]] = col.pal(length(which(sorted_SUR_assemblages$x > 1)))#[sample(length(which(sorted_SUR_assemblages$x > 1)))]
            col[sorted_SUR_assemblages$ix[sorted_SUR_assemblages$x < 1]] = rep("grey",length(which(sorted_SUR_assemblages$x < 1)))
            
            cat(length(which(sorted_SUR_assemblages$x > 1)),"colors\n")
          }
        } else if (manual)
        {
          # col = c(color.pal1(6),color.pal2(3),color.pal3(3),rev(color.pal4(4)))
          # col = c("#7e170f","#e30a14","#eb611d","#f4961b",
          #         "#fdc611","#f3e600","#94c56a","#47ad42",
          #         "#06652e","#838b8b","#ccbfb0","#f0f8fe",
          #         "#6cc6d9","#509fd8","#2e58a5","#282e67")
          
          if (nb_topics == 16)
          {
            col = c("#7e170f","#e30a14","#eb611d","#f4961b",
                    muted("#fdc611",l=83,c=200),"#f3e600",muted("#94c56a",l=83,c=60),"#47ad42",
                    "#06652e","#838b8b","#ccbfb0","#f0f8fe",
                    "#6cc6d9","#509fd8","#2e58a5","#282e67")
          } else if (nb_topics == 14)
            col = c("#7e170f","#eb611d","#f4961b",
                    muted("#fdc611",l=83,c=200),"#f3e600",muted("#94c56a",l=83,c=60),"#47ad42",
                    "#06652e","#838b8b","#ccbfb0","#f0f8fe",
                    "#6cc6d9","#509fd8","#2e58a5")
          
          
          # col[sorted_SUR_assemblages$ix] = col
          if (variable_alpha.delta)
          {
            if (alpha == 50/nb_topics && delta == 0.1)
            {
              col = col[c(14,7,9,12,
                          5,3,6,10,
                          4,11,13,8,
                          1,2)]
            } else if (alpha == 0.05/nb_topics && delta == 0.1)
            {
              col = col[c(1,12,2,6,
                          14,4,3,5,
                          7,13,11,8,
                          9,10)]
            } else if (alpha == 1 && delta == 1)
            {
              col = col[c(1,12,11,10,
                          8,3,7,9,
                          5,2,4,6,
                          13,14)]
            }
          } else if (variable_reals)
          {
            if (real.number == 1)
            {
              col = col[c(1,13,2,7,
                          15,8,4,6,
                          3,5,9,12,
                          14,16,11,10)]
            } else if (real.number == 37)
            {
              col = col[c(14,7,10,13,
                          12,11,5,16,
                          6,8,3,9,
                          2,15,1,4)]
            } else if (real.number == 64)
            {
              col = col[c(10,8,5,15,
                          14,11,3,1,
                          4,12,9,6,
                          7,13,2,16)]
            }
          } else if (variable_Ks)
          {
            if (nb_topics == 3)
              col = col[c(1,13,3)]
            else if (nb_topics == 5)
              col = col[c(1,13,8,3,
                          15)]
            else if (nb_topics == 8)
              col = col[c(1,13,15,3,
                          8,4,7,12)]
            else if (nb_topics == 12)
              col = col[c(1,13,3,15,
                          7,8,4,6,
                          14,12,11,10)]
          }
          # } else
          #   topic_reordering = 1:9
        } else
        {  
          # spatial_topicmix_kriged_all_topics[,2+(1:nb_topics)] = spatial_topicmix_kriged_all_topics[,2 + sorted_SUR_assemblages$ix[1:nb_topics]]
          
          # Computing the latitude of topics' barycentres in surface:
          mean_abs_latitude = vector(length = ncol(documents), mode = "numeric")
          for (k in 1:ncol(documents))
          {
            # mean_latitude[k] = weighted.mean(coord$y[rownames(coord) %in% rownames(documents) & stations_depths[,2] == "SUR"],
            #                                  w=documents[stations_depths[rownames(coord) %in% rownames(documents),2] == "SUR",k]/sum(documents[stations_depths[rownames(coord) %in% rownames(documents),2] == "SUR",k]))
            mean_abs_latitude[k] = weighted.mean(abs(coord$y)[rownames(coord) %in% rownames(documents) & stations_depths[,2] == "SUR"],
                                                 w=documents[stations_depths[rownames(coord) %in% rownames(documents),2] == "SUR",k]/sum(documents[stations_depths[rownames(coord) %in% rownames(documents),2] == "SUR",k]))  
          }
          
          col = vector(length = nb_topics, mode = "character")
          # Dark red color for the most abundant topic:
          col[1] = "#7F0000"
          # sorted_latitudes = sort.int(mean_latitude[sorted_SUR_assemblages$ix[1:nb_topics1]][-1],decreasing=T,index.return=T)
          # Sorting the nb_topics1 most abundant topics by the absolute latitude of their barycentre, excluding the most abundant one:
          sorted_abs_latitudes = sort.int(mean_abs_latitude[sorted_SUR_assemblages$ix[1:nb_topics]][-1],decreasing=T,index.return=T)
          # Blue palette for the topics above 50°:
          if (length(which(sorted_abs_latitudes$x>arctic_bound)) == 1)
          {
            col[sorted_abs_latitudes$ix[sorted_abs_latitudes$x>arctic_bound]+1] = "cyan"
          } else if (length(which(sorted_abs_latitudes$x>arctic_bound)) == 2)
          {
            col[sorted_abs_latitudes$ix[sorted_abs_latitudes$x>arctic_bound]+1] = c("cyan","#007FFF")
          } else if (length(which(sorted_abs_latitudes$x>arctic_bound)) > 2)
            col[sorted_abs_latitudes$ix[sorted_abs_latitudes$x>arctic_bound]+1] = color.pal.blue(length(which(sorted_abs_latitudes$x>arctic_bound)))
          # Red-yellow palette for topics between 50° and 20°, exclusing the dark red color already used for the most abundant topic:
          if (length(which(sorted_abs_latitudes$x<arctic_bound & sorted_abs_latitudes$x>tropic_bound)) > 0)
            col[sorted_abs_latitudes$ix[sorted_abs_latitudes$x<arctic_bound & sorted_abs_latitudes$x>tropic_bound]+1] = color.pal.warm(length(which(sorted_abs_latitudes$x<arctic_bound & sorted_abs_latitudes$x>tropic_bound))+1)[-1]
          # Green palette for topics between 0° and 20°:
          if (length(which(sorted_abs_latitudes$x<tropic_bound)) > 0)
            col[sorted_abs_latitudes$ix[sorted_abs_latitudes$x<tropic_bound]+1] = color.pal.green(length(which(sorted_abs_latitudes$x<tropic_bound)))
        }
        
        pie.colors.matrix = matrix(nrow = length(which(rownames(pie_positions.SUR) %in% rownames(documents))), 
                                   ncol=nb_topics, 
                                   data = NA,
                                   dimnames = list(rownames(pie_positions.SUR)[rownames(pie_positions.SUR) %in% rownames(documents)],
                                                   NULL))
        pie.slices.matrix = matrix(nrow = length(which(rownames(pie_positions.SUR) %in% rownames(documents))), 
                                   ncol=nb_topics, 
                                   data = NA,
                                   dimnames = list(rownames(pie_positions.SUR)[rownames(pie_positions.SUR) %in% rownames(documents)],
                                                   NULL))
        # pie.colors.matrix = matrix(nrow = 2, 
        #                            ncol=nb_topics, 
        #                            data = NA,
        #                            dimnames = list(c("TARA_148 SUR","TARA_149 SUR"),
        #                                            NULL))
        # pie.slices.matrix = matrix(nrow = 2, 
        #                            ncol=nb_topics, 
        #                            data = NA,
        #                            dimnames = list(c("TARA_148 SUR","TARA_149 SUR"),
        #                                            NULL))
        
        ii_pie = 0
        for (i_pie in which(rownames(documents) %in% rownames(pie_positions.SUR)))
          # for (i_pie in which(row.names %in% c("TARA_148 SUR","TARA_149 SUR")))
          # for (i_pie in which(stations_depths[rownames(coord) %in% rownames(documents),2] == "SUR"))
        {
          ii_pie = ii_pie+1
          station_topics = which(documents[i_pie,1:nb_topics] > 0)
          pie.colors.matrix[ii_pie,] = c(col[station_topics],rep(NA,(nb_topics)-length(station_topics)))
          pie.slices.matrix[ii_pie,] = c(as.numeric(documents[i_pie,station_topics]),rep(NA,(nb_topics)-length(station_topics)))
          if (length(which(pie.colors.matrix[ii_pie,] == "grey")) > 1)
          {
            pie.slices.matrix[ii_pie,which(pie.colors.matrix[ii_pie,] == "grey")[1]] = sum(pie.slices.matrix[ii_pie,which(pie.colors.matrix[ii_pie,] == "grey")])
            pie.slices.matrix[ii_pie,which(pie.colors.matrix[ii_pie,] == "grey")[-1]] = NA
            pie.colors.matrix[ii_pie,which(pie.colors.matrix[ii_pie,] == "grey")[-1]] = NA
          }
        }
        #####
        # space.pies(spatial_topicmix_kriged_all_topics$x[stations_depths[rownames(coord) %in% row.names,2] == "SUR"],
        #            spatial_topicmix_kriged_all_topics$y[stations_depths[rownames(coord) %in% row.names,2] == "SUR"],
        space.pies(spatial_topicmix_kriged_all_topics$x[rownames(documents) %in% rownames(pie_positions.SUR)],
                   spatial_topicmix_kriged_all_topics$y[rownames(documents) %in% rownames(pie_positions.SUR)],
                   pie.slices = pie.slices.matrix, pie.colors = pie.colors.matrix,
                   # pie.radius=3, pie.space=0.01,
                   pie.radius=6, pie.space=0.4,
                   link=TRUE, seg.lwd=1, seg.col=1, seg.lty=1, 
                   coord= pie_positions.SUR[rownames(pie_positions.SUR) %in% rownames(documents),])
        
        # space.pies(spatial_topicmix_kriged_all_topics$x[row.names %in% c("TARA_148 SUR","TARA_149 SUR")],
        #            spatial_topicmix_kriged_all_topics$y[row.names %in% c("TARA_148 SUR","TARA_149 SUR")],
        #            pie.slices = data.frame(pie.slices.matrix), pie.colors = data.frame(pie.colors.matrix),
        #            # pie.radius=3, pie.space=0.01,
        #            pie.radius=6, pie.space=0.4,
        #            link=TRUE, seg.lwd=1, seg.col=1, seg.lty=1, 
        #            coord= pie_positions.SUR[rownames(pie_positions.SUR) %in% c("TARA_148 SUR","TARA_149 SUR"),])
      } else if (dominant)
      {
        dominant_k = vector(length = length(which(rownames(pie_positions.SUR) %in% rownames(documents))), mode = "numeric")
        ii = 1
        for (i in which(rownames(documents) %in% rownames(pie_positions.SUR)))
        {
          dominant_k[ii] = sort.int(documents[i,],decreasing=T,index.return = T)$ix[1]
          names(dominant_k)[ii] = rownames(documents)[i]
          ii = ii+1
        }
        
        if (monochrom)
        {
          if (fixed_nb_colors && !nb_colors > length(levels(as.factor(dominant_k))))
          { 
            col.array = c(col.pal(nb_colors),rep("grey",nb_topics-nb_colors))
          } else
          {
            col.array = col.pal(length(levels(as.factor(dominant_k))))#[sample(length(levels(as.factor(dominant_k))))]
            cat(length(levels(as.factor(dominant_k))),"colors\n")
          }
          
          col = vector(length = nb_topics, mode = "character")
          i = 1
          # Assigning colors in the order of the most common dominant assemblage to the least common
          for (k in names(table(as.factor(dominant_k))))
          {
            k = as.numeric(k)
            if (k != 0)
            {
              col[k] = col.array[i]
              i = i+1
            }
          }
        }
        
        # pvpick_Hellinger = pvpick(pvclust_Hellinger,alpha = 0.96)
        # nb_clust = length(pvpick_Hellinger$clusters)
        # 
        # grp_Hellinger = rep(1,nrow(documents))
        # names(grp_Hellinger) = rownames(documents)
        # for (i_clust in 1:nb_clust)
        # {
        #   grp_Hellinger[names(grp_Hellinger) %in% pvpick_Hellinger$clusters[[i_clust]]] = i_clust+1
        # }
        
        # Plotting stations associated with each dominant assemblage  
        for (k in levels(as.factor(dominant_k)))
        {
          k = as.numeric(k)
          if (k !=0)
          {
            pies.coord = pie_positions.SUR[rownames(pie_positions.SUR) %in% rownames(documents),][dominant_k == k,]
            stations.coord = coord[rownames(coord) %in% names(dominant_k),][dominant_k == k,]
            segments(x0=stations.coord[,2],y0=stations.coord[,1],x1=pies.coord[,1],y1=pies.coord[,2])
            points(pie_positions.SUR[rownames(pie_positions.SUR) %in% rownames(documents),][dominant_k == k,],
                   pch = 21, cex=7, bg = col[k])
          }
        }
      } 
    } else if (k_by_k)
    {
      par(mar = c(1,1,1,1), oma=c(1,1,10,1))
      for (k in 1:nb_topics)
      {
        cat(k,"\n")
        
        # plot(bat, image = T, land = T, lty = 0, bty = "n", xaxt = "n", yaxt = "n", ann = F, bpal = list(c(0, max(bat), greys[1]), c(min(bat), 0, blues))) #plot map without isobaths
        plot(bat, image = T, land = T, lty = 0, bty = "n", xaxt = "n", yaxt = "n", ann = F, bpal = list(c(0, max(bat), greys[1]), c(min(bat), 0, blues[3])))
        plot(bat, lwd = 0.2, deep = 0, shallow = 0, step = 0, add = TRUE)
        
        # col = vector(length = length(which(rownames(pie_positions.SUR) %in% rownames(documents))), mode = "numeric")
        # for (i_station in 1:length(which(rownames(pie_positions.SUR) %in% rownames(documents))))
        
        # plot = ggplot(data = data.frame(x=pies.coord$x,y=pies.coord$y,z.pred=documents_SUR[,k])) +
        #                 geom_point(aes(x,y,colour=z.pred)) +
        #                 scale_colour_gradientn(colours = c("darkblue","firebrick2"), space = "Lab")
        
        nColor = 20
        colors = paletteer_c(palette = "viridis::viridis", n = nColor)
        # colors = colorRampPalette(c("darkblue","white","red"),space = "Lab")(nColor)
        # colors = topo.colors(nColor)
        
        if (group == "AllTaxa")
        {
          assemblage_order = c(1,3,9,7,
                               10,8,4,6,
                               11,16,15,12,
                               2,13,5,14)
        } else
          assemblage_order = 1:nb_topics
        
        # Cut the continuous range into bins, including dummy 0 and 1 at the beginning to scale the min and max of the color range: 
        # cut.doc.0.1 = as.numeric(cut(c(0,1,documents_SUR[,assemblage_order[k]]),nColor))[-c(1,2)]
        cut.doc.0.1 = as.numeric(cut(documents_SUR[,assemblage_order[k]],nColor))
        
        pies.coord = pie_positions.SUR[rownames(pie_positions.SUR) %in% rownames(documents),]
        stations.coord = coord[rownames(coord) %in% rownames(documents) & rownames(coord) %in% rownames(pie_positions.SUR) & stations_depths[,2] == "SUR",]
        segments(x0=stations.coord[,2],y0=stations.coord[,1],x1=pies.coord[,1],y1=pies.coord[,2])
        points(pie_positions.SUR[rownames(pie_positions.SUR) %in% rownames(documents),],
               pch = 21, cex=7, 
               # bg = col.pal(nrow(documents_SUR))[findInterval(documents_SUR[,k], sort(documents_SUR[,k]))])
               bg = colors[cut.doc.0.1])
      }
      mtext(group, outer=TRUE,  cex=9, line=-0.5, adj = 0)
      grid.newpage()
    }
  }
  # plot.new()
  # vps = baseViewports()
  # pushViewport(vps$figure)
  # print(g2)
  # vp1 = plotViewport(c(0,0,0,0))
  # print(g2,vp = vp1) 
  dev.off()

  ################
  
  
  # pdf(paste0(figure_folder,"/Normalized_PCoA.axis1.vs.basin.I.over.autocorr_",div_threshold,"plusOTUs",surf_DCM_insert,noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
  #     width=6.9,height=6.7)
  # plot(ggplotGrob(plot.norm.VI.PCoA.axis1.vs.basin.I))
  # dev.off()
}

# Fig. 3
{
  size = 1.2
  
  x = as.vector(diversity)[selected_groups]
  y = NormalizedVI_pcoa[[3]][,1]
  x.2000 = x[x<2000]
  y.2000 = y[x<2000]
  x.cor.pos = 0.8
  y.cor.pos = 0.05
  cor.test = cor.test(log10(x.2000),y.2000,na.rm=T,method="spearman")
  plot.norm.VI.PCoA.axis1.vs.diversity = cor.plot(x = x,
                                                  y = y,
                                                  x.lab ="Diversity (#OTUs), log scale",
                                                  y.lab = "Biogeographic axis 1",
                                                  y.lab.hjust = 0.5,
                                                  size.factor = size,
                                                  x.log = T,
                                                  y.log = F,
                                                  fit = F,
                                                  mar.vect = c(5,15,5.5,5)) +
    geom_smooth(data = data.frame(x = x.2000, y = y.2000),aes(x,y),method='lm',col="black",se=T) +
    # geom_smooth(aes(x[x<2000], y[x<2000]),method='lm',col="black",se=T, inherit.aes = T) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    annotate(geom="text", 
             x=(max(x)/min(x))^x.cor.pos*min(x),
             y=y.cor.pos*(max(y)-min(y))+min(y),
             label=bquote(atop(rho[S]==.(format(cor.test$estimate,digits=2,nsmall=2)),P==.(format(cor.test$p.value,digits=1)))),
             size=7.5*size*0.9)
  
  # excluded_groups = "Porifera"
  x = size_relativeAbund[selected_groups]
  y = NormalizedVI_pcoa[[3]][,2]
  plot.norm.VI.PCoA.axis2.vs.body.size = cor.plot(x = x,
                                                  y = y,
                                                  # excluded.points = which(taxo_groups[selected_groups] %in% excluded_groups),
                                                  x.lab = expression("Body size ("*mu*"m), log scale"),
                                                  y.lab = "Biogeographic axis 2",
                                                  y.lab.hjust = 0.5,
                                                  size.factor = size,
                                                  x.log = T,
                                                  y.log = F,
                                                  fit = T,
                                                  fit.display = "spearman.p",
                                                  fit.size.factor = 0.9,
                                                  x.cor.pos = 0.8,
                                                  y.cor.pos = 0.05,
                                                  mar.vect = c(5,5,3,5))
  
  excluded_groups = "Porifera"
  detached_groups = c("RAD-C","Ascomycota")
  boxplot.norm.VI.PCoA.axis2.functions = box.plot(classes = factor(dominant_function1, c("Phototrophs","Phagotrophs","Metazoans","Parasites")),
                                                            values = NormalizedVI_pcoa[[3]][!names(NormalizedVI_pcoa[[3]][,2]) %in% c("Dinophyceae","Collodaria"),2],
                                                            detached.points = which(names(dominant_function1) %in% detached_groups),
                                                            excluded.points = which(names(dominant_function1) %in% excluded_groups),
                                                            y.lab = "Biogeographic axis 2",
                                                            y.lab.hjust = 0.39,
                                                            size.factor = size,
                                                            x.angle = 16,
                                                            fit = T,
                                                            fit.size.factor = 0.85,
                                                            x.cor.pos=0.85,y.cor.pos=0.1,
                                                            mar.vect = c(5,5,1,5)) +
    theme(axis.title.y = element_blank(),
          axis.text.x = element_text(vjust = 0.8, hjust = 0.7))
  
  pdf(paste0(figure_folder,"/Fig3_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1_excluded.groups.pdf"),
      width=7.5*2.2/1.2/2*3,height=12/3*4/1.2/2,onefile=F)
  g1 = ggplotGrob(plot.norm.VI.PCoA.axis1.vs.diversity)
  g2 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.body.size)
  g3 = ggplotGrob(boxplot.norm.VI.PCoA.axis2.functions)
  grid.arrange(grobs = list(g1, g2, g3), nrow = 1, ncol = 3)
  dev.off()
  
  ##################
  cor.test(as.vector(diversity)[selected_groups & as.vector(diversity) < 2000],
           NormalizedVI_pcoa[[3]][as.vector(diversity)[selected_groups] < 2000,1],
           # method="kendall")
           method="spearman")
           # method="pearson")
  cor.test(size_relativeAbund[selected_groups & !taxo_groups %in% "Porifera"],
           NormalizedVI_pcoa[[3]][!taxo_groups[selected_groups] %in% "Porifera",2],
           # method="kendall")
           method="spearman")
           # method="pearson")
}

# Fig. 3bis
{
  div_threshold = 100
  
  excluded_groups = c("RAD-C","Porifera")
  detached_groups = NULL
  x = size_relativeAbund[selected_groups]
  y = I_square.observed_w.mean[selected_groups,1]
  plot.spatial.autocorr.vs.body.size = cor.plot(x = x,
                                                y = y,
                                                excluded.points = which(taxo_groups[selected_groups] %in% excluded_groups),
                                                x.lab = expression("Body size ("*mu*"m)"),
                                                y.lab = "Short-distance spatial autocorrelation",
                                                y.lab.hjust = 0.5,
                                                x.log = T,
                                                y.log = F,
                                                fit = T,
                                                fit.display = "pearson.spearman",
                                                x.cor.pos = 0.6,
                                                y.cor.pos = 0.1,
                                                mar.vect = c(5,5,1,5))
  
  boxplot.spatial.autocorr.functions = box.plot(classes = factor(dominant_function1, c("Phototrophs","Phagotrophs","Metazoans","Parasites")),
                                                values = I_square.observed_w.mean[selected_groups & !rownames(I_square.observed_w.mean) %in% c("Dinophyceae","Collodaria"),1],
                                                detached.points = if (!is.null(detached_groups)) which(names(dominant_function1) %in% detached_groups) else NULL,
                                                excluded.points = if (!is.null(excluded_groups)) which(names(dominant_function1) %in% excluded_groups) else NULL,
                                                y.lab = "Short-distance spatial autocorrelation",
                                                y.lab.hjust = 0.39,
                                                mar.vect = c(5,7,-7,5))
  
  excluded_groups = c("Porifera","RAD-C")
  detached_groups = NULL
  x = size_relativeAbund[selected_groups]
  y = charac_scale[selected_groups,1]
  plot.autocorr.scale.vs.body.size = cor.plot(x = x,
                                              y = y,
                                              excluded.points = which(taxo_groups[selected_groups] %in% excluded_groups),
                                              x.lab = expression("Body size ("*mu*"m)"),
                                              y.lab = "Scale of biogeographic organization (km)",
                                              y.lab.hjust = 0.5,
                                              x.log = T,
                                              y.log = F,
                                              fit = T,
                                              fit.display = "pearson.spearman",
                                              x.cor.pos = 0.7,
                                              y.cor.pos = 0.1,
                                              mar.vect = c(5,5,1,5))
  
  boxplot.autocorr.scale.functions = box.plot(classes = factor(dominant_function1, c("Phototrophs","Phagotrophs","Metazoans","Parasites")),
                                              values = charac_scale[selected_groups & !rownames(charac_scale) %in% c("Dinophyceae","Collodaria"),1],
                                              detached.points = if (!is.null(detached_groups)) which(names(dominant_function1) %in% detached_groups) else NULL,
                                              excluded.points = if (!is.null(excluded_groups)) which(names(dominant_function1) %in% excluded_groups) else NULL,
                                              y.lab = "Scale of biogeographic organization (km)",
                                              y.lab.hjust = 0.39,
                                              mar.vect = c(5,7,-7,5))
  
  excluded_groups = NULL
  detached_groups = NULL
  x = size_relativeAbund[selected_groups]
  y = basin_I[selected_groups,1]/I_square.observed_w.mean[selected_groups,1]
  plot.basin.autocorr.vs.body.size = cor.plot(x = x,
                                              y = y,
                                              excluded.points = which(taxo_groups[selected_groups] %in% excluded_groups),
                                              x.lab = expression("Body size ("*mu*"m)"),
                                              y.lab = "Whithin-basin homogeneity",
                                              y.lab.hjust = 0.5,
                                              x.log = T,
                                              y.log = F,
                                              fit = T,
                                              fit.display = "pearson.spearman",
                                              x.cor.pos = 0.6,
                                              y.cor.pos = 0.2,
                                              mar.vect = c(5,5,1,5))
  
  boxplot.basin.autocorr.functions = box.plot(classes = factor(dominant_function1, c("Phototrophs","Phagotrophs","Metazoans","Parasites")),
                                              values = (basin_I[,1]/I_square.observed_w.mean[,1])[selected_groups & !rownames(basin_I) %in% c("Dinophyceae","Collodaria")],
                                              detached.points = if (!is.null(detached_groups)) which(names(dominant_function1) %in% detached_groups) else NULL,
                                              excluded.points = if (!is.null(excluded_groups)) which(names(dominant_function1) %in% excluded_groups) else NULL,
                                              y.lab = "Whithin-basin homogeneity",
                                              y.lab.hjust = 0.39,
                                              mar.vect = c(5,7,-7,5))
  
  excluded_groups = NULL
  detached_groups = NULL
  x = size_relativeAbund[selected_groups]
  y = lat_I[selected_groups,1]/I_square.observed_w.mean[selected_groups,1]
  plot.lat.autocorr.vs.body.size = cor.plot(x = x,
                                            y = y,
                                            excluded.points = which(taxo_groups[selected_groups] %in% excluded_groups),
                                            x.lab = expression("Body size ("*mu*"m)"),
                                            y.lab = "Latitudinal symmetry",
                                            y.lab.hjust = 0.5,
                                            x.log = T,
                                            y.log = F,
                                            fit = T,
                                            fit.display = "pearson.spearman",
                                            x.cor.pos = 0.8,
                                            y.cor.pos = 0.85,
                                            mar.vect = c(5,5,1,5))
  
  boxplot.lat.autocorr.functions = box.plot(classes = factor(dominant_function1, c("Phototrophs","Phagotrophs","Metazoans","Parasites")),
                                            values = (lat_I[,1]/I_square.observed_w.mean[,1])[selected_groups & !rownames(lat_I) %in% c("Dinophyceae","Collodaria")],
                                            detached.points = if (!is.null(detached_groups)) which(names(dominant_function1) %in% detached_groups) else NULL,
                                            excluded.points = if (!is.null(excluded_groups)) which(names(dominant_function1) %in% excluded_groups) else NULL,
                                            y.lab = "Latitudinal symmetry",
                                            y.lab.hjust = 0.39,
                                            mar.vect = c(5,7,-7,5))
  
  pdf(paste0(figure_folder,"/Fig3.bis_no.axes_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1_excluded.groups_stars.pdf"),
      width=7.5*2.2/1.2/2*2,height=12/3*4/1.2/2*4,onefile=F)
  g1 = ggplotGrob(plot.spatial.autocorr.vs.body.size)
  g2 = ggplotGrob(boxplot.spatial.autocorr.functions)
  g3 = ggplotGrob(plot.autocorr.scale.vs.body.size)
  g4 = ggplotGrob(boxplot.autocorr.scale.functions)
  g5 = ggplotGrob(plot.basin.autocorr.vs.body.size)
  g6 = ggplotGrob(boxplot.basin.autocorr.functions)
  g7 = ggplotGrob(plot.lat.autocorr.vs.body.size)
  g8 = ggplotGrob(boxplot.lat.autocorr.functions)
  grid.arrange(grobs = list(g1, g2, g3, g4, g5, g6, g7, g8), nrow = 4, ncol = 2)
  dev.off()
  
  ##################
  cor.test(size_relativeAbund[selected_groups & !taxo_groups %in% c("RAD-C","Porifera")],
           I_square.observed_w.mean[selected_groups & !taxo_groups %in% c("RAD-C","Porifera"),1],
           method="kendall")
           # method="spearman")
  cor.test(size_relativeAbund[selected_groups & !taxo_groups %in% c("RAD-C","Porifera")],
           charac_scale[selected_groups & !taxo_groups %in% c("RAD-C","Porifera"),1],
           method="kendall")
           # method="spearman")
  cor.test(size_relativeAbund[selected_groups],
           basin_I[selected_groups,1],
           method="kendall")
           # method="spearman")
  cor.test(size_relativeAbund[selected_groups],
           lat_I[selected_groups,1],
           method="kendall")
           # method="spearman")
  ####################
  # x = log(size_relativeAbund[selected_groups & !taxo_groups %in% c("Porifera","RAD-C","RAD-A","Collodaria")])
  # y = charac_scale[selected_groups & !taxo_groups %in% c("Porifera","RAD-C","RAD-A","Collodaria"),1]
  x1 = charac_scale[selected_groups,1]
  x2 = basin_I[selected_groups,1]
  x3 = lat_I[selected_groups,1]
  x4 = I_square.observed_w.mean[selected_groups,1]
  y = log(size_relativeAbund[selected_groups])
  model0 = lm(y ~ x1 + x2 + x3)
  model1 = lm(y ~ x2 + x3)
  anova(model0)
  anova(model1,model0)
  summary(model0)
  ######
  g = factor(dominant_function1[!names(dominant_function1) %in% c("Porifera")],c("Phagotrophs","Phototrophs","Metazoans","Parasites"))
  Meta.vs.Photo = c(0,1,-1,0)
  Phago.vs.Meta = c(1,0,-1,0)
  Phago.vs.Photo = c(1,-1,0,0)
  contrasts(g) = cbind(Meta.vs.Photo,Phago.vs.Meta,Phago.vs.Photo)
  y = charac_scale[selected_groups & !taxo_groups %in% c("Dinophyceae","Collodaria") & !taxo_groups %in% c("Porifera"),1]
  model = lm(y ~ g)
  anova(model)
  summary(model)
  #
  model = aov(y ~ g)
  CList = list("Meta.vs.Photo" = 1,
               "Phago.vs.Meta" = 2,
               "Phago.vs.Photo" = 3)
  summary(model, split=list(g=CList))
  summary.lm(model)
  # 
  Meta.vs.rest = c(1,1,-3,1)
  contrasts(g) = Meta.vs.rest
  model = aov(y ~ g)
  CList = list("Meta.vs.rest" = 1)
  summary(model, split=list(g=CList))
  #
  Meta.vs.Photo = c(0,1,-1,0)
  contrasts(g) = Meta.vs.Photo
  model = aov(y ~ g)
  CList = list("Meta.vs.Photo" = 1)
  summary(model, split=list(g=CList))
  
  ######
  x.control = NormalizedVI_pcoa[[3]][,1]
  x.control1 = I_square.observed_w.mean[selected_groups,1]
  y = NormalizedVI_pcoa[[3]][,2]
  y1 = charac_scale[selected_groups,1]
  y1.norm = charac_scale[selected_groups,1]/I_square.observed_w.mean[selected_groups,1]
  y2 = lat_I[selected_groups,1]
  y2.norm = lat_I[selected_groups,1]/I_square.observed_w.mean[selected_groups,1]
  y3 = basin_I[selected_groups,1]
  y3.norm = basin_I[selected_groups,1]/I_square.observed_w.mean[selected_groups,1]
  x2 = log10(size_relativeAbund)[selected_groups]
  pcor.test(y,x2,x.control)
  pcor.test(y1,x2,x.control)
  pcor.test(y2,x2,x.control)
  pcor.test(y3,x2,x.control)
  #
  y = NormalizedVI_pcoa[[3]][!taxo_groups[selected_groups] %in% c("Porifera","Dinophyceae","Collodaria"),2]
  y1 = charac_scale[selected_groups & !taxo_groups %in% c("Porifera","Dinophyceae","Collodaria"),1]
  y2 = lat_I[selected_groups & !taxo_groups %in% c("Porifera","Dinophyceae","Collodaria"),1]
  y3 = basin_I[selected_groups & !taxo_groups %in% c("Porifera","Dinophyceae","Collodaria"),1]
  x2 = log10(size_relativeAbund)[selected_groups & !taxo_groups %in% c("Porifera","Dinophyceae","Collodaria")]
  x.control = NormalizedVI_pcoa[[3]][!taxo_groups[selected_groups] %in% c("Porifera","Dinophyceae","Collodaria"),1]
  g = factor(dominant_function1[!names(dominant_function1) %in% c("Porifera")],c("Phagotrophs","Phototrophs","Metazoans","Parasites"))
  Meta.vs.Photo = c(0,1,-1,0)
  Phago.vs.Meta = c(1,0,-1,0)
  Phago.vs.Photo = c(1,-1,0,0)
  contrasts(g) = cbind(Meta.vs.Photo,Phago.vs.Meta,Phago.vs.Photo)
  # anova(lm(y ~ g + x.control))
  anova(lm(y1 ~ g + x.control))
  anova(lm(y2 ~ g + x.control))
  anova(lm(y3 ~ g + x.control))
  anova(lm(y1 ~ g + x.control + x2))
  #######
  y = (colSums(varpart.env.spatial[[1]][2:3,])/(colSums(varpart.env.spatial[[1]][1:2,])))[selected_groups & taxo_groups != "Porifera"]
  x = log10(size_relativeAbund)[selected_groups & taxo_groups != "Porifera"]
  x.control = NormalizedVI_pcoa[[3]][taxo_groups[selected_groups] != "Porifera",1]
  pcor.test(y,x,x.control)
  anova(lm(y ~ x + x.control))
  summary(lm(y ~ x + x.control))
  varpart(y,x,x.control)
  anova(rda(y,x,x.control))
  #
  y = (colSums(varpart.env.spatial[[1]][2:3,])/(colSums(varpart.env.spatial[[1]][1:2,])))[selected_groups & !taxo_groups %in% c("Porifera","Dinophyceae","Collodaria")]
  x = log10(size_relativeAbund)[selected_groups & !taxo_groups %in% c("Porifera","Dinophyceae","Collodaria")]
  x.control = NormalizedVI_pcoa[[3]][!taxo_groups[selected_groups] %in% c("Porifera","Dinophyceae","Collodaria"),1]
  g = factor(dominant_function1[!names(dominant_function1) %in% c("Porifera")],c("Phagotrophs","Phototrophs","Metazoans","Parasites"))
  Meta.vs.Photo = c(0,1,-1,0)
  Phago.vs.Meta = c(1,0,-1,0)
  Phago.vs.Photo = c(1,-1,0,0)
  contrasts(g) = cbind(Meta.vs.Photo,Phago.vs.Meta,Phago.vs.Photo)
  anova(lm(y ~ g + x.control))
  anova(lm(y ~ g + x.control + x))
}

# Old fig. 4 v1:
{
  div_threshold = 100
  
  ##########
  dominant_function0 = dominant_function0[-69]
  # alpha = rep(1,length(taxo_groups[selected_groups]))
  # alpha[dominant_function0 %in% c("copepoda","pteropoda")] = 0 
  # Changes how the factors are stored (order of levels()) so that geom_boxplot plots them by deceasing number of groups:
  # dominant_function0 = factor(dominant_function0,names(sort(table(as.factor(dominant_function0)),decreasing = T)))
  # Changes how the factors are stored (order of levels()) so that geom_boxplot plots them in the specified order:
  if (div_threshold == 100)
  {
    dominant_function0 = factor(dominant_function0,c("Bacillariophyta","Other phototrophs","Collodaria","Phagotrophs","Dinophyceae","Gel. carn. filterers","Pteropoda","Copepoda","Other metazoa","Parasites"))
    point_groups = c("Dinophyceae","Copepoda","Collodaria","Bacillariophyta","Pteropoda")
  } else if (div_threshold == 1000)
  {
    dominant_function0 = factor(dominant_function0,c("Bacillariophyta","Other phototrophs","Collodaria","Phagotrophs","Dinophyceae","Gel. carn. filterers","Copepoda","Parasites"))
    point_groups = c("Dinophyceae","Copepoda","Collodaria","Bacillariophyta")
  } 
  
  # file_name = paste0(figure_folder,"/Fig3.1_ggplot",abiotic_pca_insert,biotic_pca_insert,stdzation_insert,
  #                    "_independentVariableSelection_selected",if (div_threshold == 100) "100+1" else if (div_threshold == 1000) "1000",surf_DCM_insert,"_noNegativeAdjR2.pdf")
  # 
  
  functions0 = c("Collodaria","Pteropoda","Copepoda","Gel. carn. filterers","Dinophyceae","Bacillariophyta","Other phototrophs","Phagotrophs","Other metazoa","Parasites")
  # "Bacillariophyta"   "Collodaria"   "Copepoda"   "Dinophyceae"   "Gel. carn. filterers" "Other phototrophs"    "Parasites"    "Phagotrophs" 
  ten_colors = c("darkorange","cadetblue","darkblue","darkturquoise","deeppink1","darkgreen","chartreuse2","firebrick2","dodgerblue1","darkgoldenrod1")
  
  library(gridExtra)
  pdf(paste0(figure_folder,"/Fig4_",div_threshold,"plusOTUs_logFractionRatio_selected100+1.pdf"),
      width=7.5*2.2/1.2*3/2,height=12/3*4/1.2*3/2)
  for (i_case in 1:2)
  {
    # i_case = 1  
    # if (i_case == 1)
    # {
    #   surf_DCM_insert = "_SUR" 
    # } else if (i_case == 2)
    #   surf_DCM_insert = "_DCM"
    
    # Tot. explained variance: 
    #########################
    plot.tot.var.vs.autocorr = ggplot(data=data.frame(x=I_square.observed_w.mean[selected_groups & diversity>div_threshold,i_case][-69],
                                                      y=colSums(varpart.env.spatial[[i_case]])[selected_groups & diversity>div_threshold][-69])) +
      geom_point(aes(x,y)) +
      # scale_y_log10() +
      theme_bw() +
      ggtitle(LETTERS[1]) +
      # geom_hline(yintercept = 1, linetype="dashed") +
      # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
      # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
      theme(axis.text = element_text(size=22),
            axis.title=element_text(size=22),
            # axis.title.x=element_text(vjust = 45),
            plot.title=element_text(hjust=0, size=24),
            plot.margin=unit(c(0.5,7,0.5,0.5),"mm")) +
      labs(x=paste("Short-distance spatial autocorr.",c("(Surface)","(DCM)")[i_case]),
           y=paste("Total variance explained\n by envir. and connectivity",c("(Surface)","(DCM)")[i_case])) +
      # labs(x=paste("Number of",if (i_case == 1) "surface" else "DCM","community types"),y="") +
      geom_smooth(aes(x,y),method='lm',col="black") 
    g1 = ggplotGrob(plot.tot.var.vs.autocorr)
    
    plot.tot.var.vs.SUR.DCM.sim = ggplot(data=data.frame(x=1 - SUR.DCM_Normalized.VI[selected_groups & diversity>div_threshold][-69],
                                                      y=colSums(varpart.env.spatial[[i_case]])[selected_groups & diversity>div_threshold][-69])) +
      geom_point(aes(x,y)) +
      # scale_y_log10() +
      theme_bw() +
      ggtitle(LETTERS[2]) +
      # geom_hline(yintercept = 1, linetype="dashed") +
      # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
      # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
      theme(axis.text = element_text(size=22),
            axis.title=element_text(size=22),
            # axis.title.x=element_text(vjust = 45),
            plot.title=element_text(hjust=0, size=24),
            plot.margin=unit(c(0.5,7,0.5,0.5),"mm")) +
      labs(x="Ssimilarity between surface and DCM",
           y=paste("Total variance explained\n by envir. and connectivity",c("(Surface)","(DCM)")[i_case])) +
      # labs(x=paste("Number of",if (i_case == 1) "surface" else "DCM","community types"),y="") +
      geom_smooth(aes(x,y),method='lm',col="black") 
    g2 = ggplotGrob(plot.tot.var.vs.SUR.DCM.sim)
    
    plot.tot.var.vs.PCoA.axis1 = ggplot(data=data.frame(x=NormalizedVI_pcoa[[3]][,1],
                                                            y=colSums(varpart.env.spatial[[i_case]])[selected_groups & diversity>div_threshold][-69])) +
      geom_point(aes(x,y)) +
      # scale_y_log10() +
      theme_bw() +
      ggtitle(LETTERS[3]) +
      # geom_hline(yintercept = 1, linetype="dashed") +
      # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
      # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
      theme(axis.text = element_text(size=22),
            axis.title=element_text(size=22),
            # axis.title.x=element_text(vjust = 45),
            plot.title=element_text(hjust=0, size=24),
            plot.margin=unit(c(0.5,7,0.5,0.5),"mm")) +
      labs(x="PCoA axis 1",
           y=paste("Total variance explained\n by envir. and connectivity",c("(Surface)","(DCM)")[i_case])) +
      # labs(x=paste("Number of",if (i_case == 1) "surface" else "DCM","community types"),y="") +
      geom_smooth(aes(x,y),method='lm',col="black") 
    g3 = ggplotGrob(plot.tot.var.vs.PCoA.axis1)
    
    # Env. connectivity fractions 1
    ###############################

    # tot_sorting_indSelec = sort.int(colSums(varpart.env.spatial[,selected_groups]), index.return = T, decreasing = T, na.last = T, method = "radix")$ix
    tot_sorting = sort.int(colSums(varpart.env.spatial[[i_case]])[selected_groups & diversity>div_threshold][-69], index.return = T, decreasing = T, na.last = T, method = "radix")$ix
    
    # "Pure biotic","Mixed biotic-abiotic","Mixed biotic-currents","Mixed biotic-abiotic-currents","Pure abiotic","Mixed abiotic-currents","Pure currents"
    barplot_data.frame = rbind(varpart.env.spatial[[i_case]][1,],varpart.env.spatial[[i_case]][2,],varpart.env.spatial[[i_case]][3,])
    dimnames(barplot_data.frame) = list(c("Purely by environment","Jointly by connectivity and envir.","Purely by connectivity"), taxo_groups_unmodified)
    barplot_data.frame = barplot_data.frame[,selected_groups & diversity > div_threshold][,-69][,tot_sorting]
    barplot.env.connectivity = ggplot(data = data.table::melt(barplot_data.frame)) +
      geom_bar(aes(x = Var2, y = value, fill = Var1), stat="identity") +
      theme_bw() +
      ggtitle(LETTERS[4]) +
      theme(axis.text = element_text(size=22),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
          axis.title=element_text(size=22),
          # axis.title.x=element_text(vjust = 45),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,7,20,0.5),"mm"),
          legend.position=c(0.5,-0.1),
          legend.text=element_text(size=22),
          legend.title=element_text(size=16, hjust = 10)) +
      labs(x="", y = paste("\nExplained variance",c("(Surface)","(DCM)")[i_case])) +
      scale_fill_manual(values=terrain.colors(3)) +
      guides(fill = guide_legend(title = NULL, direction = "vertical"))
    g4 = ggplotGrob(barplot.env.connectivity) 
    
    # plot.connectivity.vs.env = ggplot(data=data.frame(x=varpart.env.spatial[[i_case]][1,][selected_groups & diversity>div_threshold][-69],
    #                                                                       y=varpart.env.spatial[[i_case]][3,][selected_groups & diversity>div_threshold][-69])) +
    #   geom_point(aes(x,y)) +
    #   # scale_x_log10() +
    #   # scale_y_log10() +
    #   theme_bw() +
    #   ggtitle(LETTERS[5]) +
    #   # geom_hline(yintercept = 1, linetype="dashed") +
    #   # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    #   # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    #   theme(axis.text = element_text(size=22),
    #         axis.title=element_text(size=22),
    #         # axis.title.x=element_text(vjust = 45),
    #         plot.title=element_text(hjust=0, size=24),
    #         plot.margin=unit(c(0.5,0.5,0.5,0.5),"mm")) +
    #   #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    #   labs(x=paste("Variance purely explained by envir.",c("(Surface)","(DCM)")[i_case]), 
    #        y=paste("Variance purely explained\n by connectivity",c("(Surface)","(DCM)")[i_case])) +
    #   geom_smooth(aes(x,y),method='lm',col="black")
    # g5 = ggplotGrob(plot.connectivity.vs.env)
    
    plot.connectivity.env.ratio.vs.spatial.scale = ggplot(data=data.frame(x=charac_scale[selected_groups & diversity>div_threshold, i_case][-69],
                                                                 y=(varpart.env.spatial[[i_case]][3,]/varpart.env.spatial[[i_case]][1,])[selected_groups & diversity>div_threshold][-69])) +
      geom_point(aes(x,y)) +
      # scale_x_log10() +
      scale_y_log10() +
      theme_bw() +
      ggtitle(LETTERS[5]) +
      geom_hline(yintercept = 1, linetype="dashed") +
      # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
      # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
      theme(axis.text = element_text(size=22),
            axis.title=element_text(size=22),
            # axis.title.x=element_text(vjust = 45),
            plot.title=element_text(hjust=0, size=24),
            plot.margin=unit(c(0.5,7,0.5,0.5),"mm")) +
      #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
      labs(x=paste("Charac. scale of spatial autocorr.",c("(Surface)","(DCM)")[i_case]),
           y=paste("Ratio of purely connectivity-explained\n over envir.-explained variance",c("(Surface)","(DCM)")[i_case])) +
      geom_smooth(aes(x,y),method='lm',col="black")
    g5 = ggplotGrob(plot.connectivity.env.ratio.vs.spatial.scale)
    
    plot.connectivity.env.ratio.vs.PCoA.axis2 = ggplot(data=data.frame(x=NormalizedVI_pcoa[[3]][,2],
                                                                       y=(varpart.env.spatial[[i_case]][3,]/varpart.env.spatial[[i_case]][1,])[selected_groups & diversity>div_threshold][-69])) +
      geom_point(aes(x,y)) +
      # scale_x_log10() +
      scale_y_log10() +
      theme_bw() +
      ggtitle(LETTERS[6]) +
      geom_hline(yintercept = 1, linetype="dashed") +
      # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
      # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
      theme(axis.text = element_text(size=22),
            axis.title=element_text(size=22),
            # axis.title.x=element_text(vjust = 45),
            plot.title=element_text(hjust=0, size=24),
            plot.margin=unit(c(0.5,7,0.5,0.5),"mm")) +
      #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
      labs(x=paste("PCoA axis 2"), 
           y=paste("Ratio of purely connectivity-explained\n over envir.-explained variance",c("(Surface)","(DCM)")[i_case])) +
      geom_smooth(aes(x,y),method='lm',col="black")
    g6 = ggplotGrob(plot.connectivity.env.ratio.vs.PCoA.axis2)
    
    # Env. connectivity fractions 2
    ###############################
    
    plot.connectivity.env.ratio.vs.size = ggplot(data=data.frame(x=size_relativeAbund[selected_groups & diversity>div_threshold][-69],
                                                                 y=(varpart.env.spatial[[i_case]][3,]/varpart.env.spatial[[i_case]][1,])[selected_groups & diversity>div_threshold][-69])) +
      geom_point(aes(x,y)) +
      scale_x_log10() +
      scale_y_log10() +
      theme_bw() +
      ggtitle(LETTERS[7]) +
      geom_hline(yintercept = 1, linetype="dashed") +
      # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
      # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
      theme(axis.text = element_text(size=22),
            axis.title=element_text(size=22),
            # axis.title.x=element_text(vjust = 45),
            plot.title=element_text(hjust=0, size=24),
            plot.margin=unit(c(0.5,7,0.5,0.5),"mm")) +
      #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
      labs(x="Mean size (micron)", 
           y=paste("Ratio of purely connectivity-explained\n over envir.-explained variance",c("(Surface)","(DCM)")[i_case])) +
      geom_smooth(aes(x,y),method='lm',col="black")
    g7 = ggplotGrob(plot.connectivity.env.ratio.vs.size)
    
    boxplot.connectivity.env.ratio.functions = ggplot(data = data.frame(x = dominant_function0[!is.na(dominant_function0) & !dominant_function0 %in% point_groups],
                                                                    y = (varpart.env.spatial[[i_case]][3,]/varpart.env.spatial[[i_case]][1,])[selected_groups & diversity>div_threshold][-69]
                                                                    [!is.na(dominant_function0) & !dominant_function0 %in% point_groups])) +
      scale_x_discrete(limits=levels(dominant_function0)) +
      geom_boxplot(aes(x,y)) +
      scale_y_log10() +
      geom_point(data = data.frame(x = factor(point_groups),
                                   y = (varpart.env.spatial[[i_case]][3,]/varpart.env.spatial[[i_case]][1,])[selected_groups & diversity>div_threshold][-69]
                                   [dominant_function0 %in% point_groups]),
                 aes(x,y)) +
      theme_bw() +
      ggtitle(LETTERS[8]) +
      # scale_y_log10() +
      geom_hline(yintercept = 1, linetype = "dashed") +
      theme(axis.title=element_text(size=22),
            axis.text=element_text(size=22),
            plot.title=element_text(hjust=0, size=24),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            plot.margin=unit(c(0.5,7,-2,0.5),"mm")) +
      labs(x="", y=paste("Ratio of purely connectivity-explained\n over envir.-explained variance",c("(Surface)","(DCM)")[i_case]))
    g8 = ggplotGrob(boxplot.connectivity.env.ratio.functions)
    
    tot_sorting = sort.int(colSums(varpart.biotic.abiotic[[i_case]])[selected_groups & diversity>div_threshold][-69], index.return = T, decreasing = T, na.last = T, method = "radix")$ix
    
    # "Pure biotic","Mixed biotic-abiotic","Mixed biotic-currents","Mixed biotic-abiotic-currents","Pure abiotic","Mixed abiotic-currents","Pure currents"
    barplot_data.frame = rbind(varpart.biotic.abiotic[[i_case]][1,],varpart.biotic.abiotic[[i_case]][2,],varpart.biotic.abiotic[[i_case]][3,])
    dimnames(barplot_data.frame) = list(c("Purely by abiotic conditions","Jointly by biotic and abiotic conditions","Purely by biotic conditions"), taxo_groups_unmodified)
    barplot_data.frame = barplot_data.frame[,selected_groups & diversity > div_threshold][,-69][,tot_sorting]
    barplot.biotic.abiotic = ggplot(data = data.table::melt(barplot_data.frame)) +
      geom_bar(aes(x = Var2, y = value, fill = Var1), stat="identity") +
      theme_bw() +
      ggtitle(LETTERS[9]) +
      theme(axis.text = element_text(size=22),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title=element_text(size=22),
            # axis.title.x=element_text(vjust = 45),
            plot.title=element_text(hjust=0, size=24),
            plot.margin=unit(c(0.5,7,20,0.5),"mm"),
            legend.position=c(0.5,-0.1),
            legend.text=element_text(size=22),
            legend.title=element_text(size=16, hjust = 10)) +
      labs(x="", y = paste("\nExplained variance",c("(Surface)","(DCM)")[i_case])) +
      scale_fill_manual(values= rev(c("deepskyblue4","cadetblue","chartreuse3"))) +
      guides(fill = guide_legend(title = NULL, nrow = 3))
    g9 = ggplotGrob(barplot.biotic.abiotic) 
    
    ##############
    # boxplot.ratio.biotic.abiotic.functions = ggplot(data = data.frame(x = dominant_function0[!is.na(dominant_function0) & !dominant_function0 %in% point_groups],
    #                                                                   y = (varpart.biotic.abiotic[[i_case]][3,]/varpart.biotic.abiotic[[i_case]][1,])
    #                                                                   [selected_groups & diversity>div_threshold & (varpart.biotic.abiotic[[i_case]][3,]/varpart.biotic.abiotic[[i_case]][1,]) < 30]
    #                                                                   [!is.na(dominant_function0) & !dominant_function0 %in% point_groups])) +
    #   scale_x_discrete(limits=levels(dominant_function0)) +
    #   geom_boxplot(aes(x,y)) +
    #   geom_point(data = data.frame(x = factor(point_groups),
    #                                y = (varpart.biotic.abiotic[[i_case]][3,]/varpart.biotic.abiotic[[i_case]][1,])
    #                                [selected_groups & diversity>div_threshold & (varpart.biotic.abiotic[[i_case]][3,]/varpart.biotic.abiotic[[i_case]][1,]) < 30]
    #                                [dominant_function0 %in% point_groups]),
    #              aes(x,y)) +
    #   theme_bw() +
    #   ggtitle(LETTERS[4]) +
    #   scale_y_log10() +
    #   geom_hline(yintercept = 1, linetype = "dashed") +
    #   theme(axis.title=element_text(size=22),
    #         axis.text=element_text(size=22),
    #         axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    #         plot.title=element_text(hjust=0, size=24),
    #         plot.margin=unit(c(10,1,1,10),"mm")) +
    #   labs(x="", y=paste("Ratio of",if (i_case == 1) "surface" else "DCM","variance\n purely explained by biotic vs. abiotic envir."))
    # labs(x="",y="")
    
    ##############
    # file_name = paste0(figure_folder,"/Fig3.1_ggplot",abiotic_pca_insert,biotic_pca_insert,stdzation_insert,
    #                    "_independentVariableSelection_selected",if (div_threshold == 100) "100+1" else if (div_threshold == 1000) "1000",surf_DCM_insert,"_noNegativeAdjR2.pdf")
    # 
    ############  
    # g1 = ggplotGrob(boxplot.ratio.currents.env.functions)
    # g2 = ggplotGrob(plot.currents.env.ratio.vs.K)
    # g3 = ggplotGrob(plot.currents.env.ratio.vs.size)
    # #grid::grid.newpage()
    # pdf(file_name, width=7.5*2.5,height=12/3*2)
    # grid::grid.draw(gtable_cbind(g1, g2, g3),recording=F)
    # dev.off()
    # ############
    # g1 = ggplotGrob(plot.currents.env.ratio.vs.K)
    # g2 = ggplotGrob(plot.currents.env.ratio.vs.size)
    # g3 = ggplotGrob(boxplot.ratio.currents.env.functions)
    # g4 = ggplotGrob(boxplot.ratio.biotic.abiotic.functions)
    #grid::grid.newpage()
    # pdf(file_name, width=7.5*2.2/1.2,height=12/3*4/1.2)
    # grid::grid.draw(gtable_rbind(gtable_cbind(g1, g2), gtable_cbind(g3, g4)),recording=F)
    grid.arrange(grobs = list(g1,g2,g3,g4,g5,g6,g7,g8,g9),
                 # widths = rep(1,7),
                 layout_matrix = rbind(c(1, 2, 3),c(4, 5, 6),c(7, 8, 9)))
  }
  dev.off()
  ##############
  
  # plots = list(boxplot.ratio.currents.env.functions,plot.currents.env.ratio.vs.K,plot.currents.env.ratio.vs.size)
  # grobs = list()  
  # heights = list()
  # for (i in 1:length(plots))
  # {
  #   grobs[[i]] = ggplotGrob(plots[[i]])
  #   heights[[i]] = grobs[[i]]$heights[2:5]
  # }
  # maxheight = do.call(grid::unit.pmax, heights)
  # for (i in 1:length(grobs))
  # {
  #   grobs[[i]]$heights[2:5] = as.list(maxheight)
  # }
  # plot = do.call("grid.arrange", c(grobs, nrow = 1))
  
  #plot = grid.arrange(grobs = list(boxplot.ratio.currents.env.functions,plot.currents.env.ratio.vs.K,plot.currents.env.ratio.vs.size), ncol = 3, layout_matrix = matrix(1:3,nrow=1,byrow=T))
  #ggsave(file_name, plot=plot,width=7.5*2.5,height=12/3*2)
}

# Fig. 4:
{
  div_threshold = 100
  
  excluded_groups = NULL
  ##### Individually selected axes only:
  # Two first outliers for pure+mixed fractions:
  excluded_groups = c("Porifera", "Spumellaria")
  # Two first outliers for pure fractions:
  excluded_groups = c("Spumellaria", "Nemertea")
  ##### Individually selected axes only, with Benjamini-Hochberg correction:
  # First outliers for pure+mixed fractions:
  excluded_groups = c("RAD-C", "Collodaria", "Vannellida") # "Porifera" is removed. 
  # Also "Collodaria" for varpart_ratio vs. PCoA2.
  ##### Variable selection only:
  excluded_groups = "MALV-V"
  ##### No variable selection whatsoever (pure+mixed):
  excluded_groups = c("RAD-C", "MALV-IV","Prasinophyceae_Clade_7","Rhodophyta")
  # excluded_groups = c("RAD-C","Prasinophyceae_Clade_7","Rhodophyta")
  
  g1 = g1_bis = g2 = g3 = g4 = list()
  for (i_case in 1:2)
  {
    excluded_groups = NULL
    if (i_case == 2)
      excluded_groups = "Ctenophora"
    x = NormalizedVI_pcoa[[3]][,1]
    y = colSums(varpart.env.spatial[[i_case]])[selected_groups]
    plot.tot.var.vs.PCoA.axis1 = cor.plot(x = x,
                                          y = y,
                                          excluded.points = if (!is.null(excluded_groups)) which(taxo_groups[selected_groups] %in% excluded_groups) else NULL, 
                                          y.lab = paste("Total variance explained\n by Moran maps and environment"),#c("","(DCM)")[i_case]),
                                          x.lab = "PCoA axis 1",
                                          y.lab.hjust = 0.5,
                                          x.log = F,
                                          y.log = F,
                                          fit = T,
                                          fit.display = "pearson.p",
                                          # fit.display = "pearson.spearman",
                                          x.cor.pos = 0.8,
                                          y.cor.pos = 0.1,
                                          mar.vect = c(5,5,1,5))
                                          # letter = LETTERS[1])
    g1[[i_case]] = ggplotGrob(plot.tot.var.vs.PCoA.axis1)
    
    # x = I_square.observed_w.mean[selected_groups & diversity>div_threshold,i_case]
    # plot.tot.var.vs.spat.autocorr = cor.plot(x = x,
    #                                          y = y,
    #                                          y.lab = paste("Total variance explained\n by connectivity and local environment"),#c("","(DCM)")[i_case]),
    #                                          x.lab = "Amount of biogeographic structure",
    #                                          y.lab.hjust = 0.5,
    #                                          x.log = F,
    #                                          y.log = F,
    #                                          fit = T,
    #                                          x.cor.pos = 0.8,
    #                                          y.cor.pos = 0.1,
    #                                          mar.vect = c(5,5,5,9))
    # g1_bis[[i_case]] = ggplotGrob(plot.tot.var.vs.spat.autocorr)
    
    # Individually selected axes only, with Benjamini-Hochberg correction - pure+mixed:
    # excluded_groups = c("RAD-C", "Collodaria")
    excluded_groups = NULL
    if (i_case == 1)
      excluded_groups = "RAD-C"
    x = NormalizedVI_pcoa[[3]][,2]
    # y = (varpart.env.spatial[[i_case]][3,]/varpart.env.spatial[[i_case]][1,])[selected_groups & diversity>div_threshold]
    y = (colSums(varpart.env.spatial[[i_case]][2:3,])/(colSums(varpart.env.spatial[[i_case]][1:2,])))[selected_groups]
    plot.connectivity.env.ratio.vs.PCoA.axis2 = cor.plot(x = x,
                                                         y = y,
                                                         excluded.points = if (!is.null(excluded_groups)) which(taxo_groups[selected_groups] %in% excluded_groups) else NULL,
                                                         # y.lab = paste("Ratio of variance purely explained\n by connectivity over local environment"),#c("","(DCM)")[i_case]),
                                                         y.lab = paste("Ratio of variance explained\n by Moran maps over environment (log scale)"),#c("","(DCM)")[i_case]),
                                                         x.lab = "PCoA axis 2",
                                                         y.lab.hjust = 0.5,
                                                         x.log = F,
                                                         y.log = T,
                                                         fit = T,
                                                         # fit.display = "pearson.spearman",
                                                         fit.display = "spearman.p",
                                                         x.cor.pos = 0.8,
                                                         y.cor.pos = 0.1,
                                                         mar.vect = c(5,5,1,5))
                                                         # letter = LETTERS[2])
    g2[[i_case]] = ggplotGrob(plot.connectivity.env.ratio.vs.PCoA.axis2 + 
                                # geom_vline(xintercept = 1, linetype="dashed"))
                                geom_hline(yintercept = 1, linetype="dashed"))
    
    # Individually selected axes only, with Benjamini-Hochberg correction - pure+mixed:
    # excluded_groups = c("RAD-C", "Collodaria")
    excluded_groups = NULL
    if (i_case == 1)
      excluded_groups = "RAD-C"
    x = size_relativeAbund[selected_groups]
    # y = (varpart.env.spatial[[i_case]][3,]/varpart.env.spatial[[i_case]][1,])[selected_groups & diversity>div_threshold]
    y = (colSums(varpart.env.spatial[[i_case]][2:3,])/(colSums(varpart.env.spatial[[i_case]][1:2,])))[selected_groups]
    fraction.ratio.vs.size = cor.plot(x = x,
                                      y = y,
                                      excluded.points = if (!is.null(excluded_groups)) which(taxo_groups[selected_groups] %in% excluded_groups) else NULL,
                                      # y.lab = paste("Ratio of variance purely explained\n by connectivity over local environment"),#c("","(DCM)")[i_case]),
                                      y.lab = paste("Ratio of variance explained\n by Moran maps over environment (log scale)"),#c("","(DCM)")[i_case]),
                                      x.lab = expression("Body size ("*mu*"m), log scale"),
                                      y.lab.hjust = 0.5,
                                      x.log = T,
                                      y.log = T,
                                      fit = T,
                                      fit.display = "spearman.p",
                                      # fit.display = "pearson.spearman",
                                      x.cor.pos = 0.8,
                                      y.cor.pos = 0.1,
                                      mar.vect = c(5,5,1,5))
                                      # letter = LETTERS[3])
    g3[[i_case]] = ggplotGrob(fraction.ratio.vs.size + geom_hline(yintercept = 1, linetype="dashed"))
    
    # Individually selected axes only, with Benjamini-Hochberg correction - pure+mixed:
    excluded_groups = NULL
    if (i_case == 1)
      excluded_groups = "RAD-C"
    # values = (varpart.env.spatial[[i_case]][3,]/varpart.env.spatial[[i_case]][1,])[selected_groups & diversity>div_threshold
    values = (colSums(varpart.env.spatial[[i_case]][2:3,])/colSums(varpart.env.spatial[[i_case]][1:2,]))[selected_groups
                                                                                                         & !colnames(varpart.env.spatial[[i_case]]) %in% c("Dinophyceae","Collodaria")]
                                                                                                         # & !colnames(varpart.env.spatial[[i_case]]) %in% excluded_groups]
    fraction.ratio.vs.function = box.plot(classes = factor(dominant_function1, #[!names(dominant_function1) %in% excluded_groups],
                                                           c("Phototrophs","Phagotrophs","Metazoans","Parasites")),
                                          values = values,
                                          excluded.points = if (!is.null(excluded_groups)) which(names(dominant_function1) %in% excluded_groups) else NULL, 
                                          # y.lab = paste("Ratio of variance purely explained\n by connectivity over local environment"),#c("","(DCM)")[i_case]),
                                          y.lab = paste("Ratio of variance explained\n by Moran maps over environment (log scale)"),#c("","(DCM)")[i_case]),
                                          y.log = T,
                                          # y.lim = if (i_case == 1) c(min(values),4) else NULL,
                                          x.angle = 14,
                                          fit = T,
                                          x.cor.pos=0.7,
                                          y.cor.pos = 0.1,
                                          # y.lab.hjust = 0.01,
                                          y.lab.hjust = NULL,
                                          mar.vect = c(5,5,1,5))
    
                                                    # letter = LETTERS[4])
      # annotate(geom="text",
      #          x = c(1,2,3) - 0.5,
      #          y = values[c("Pelagophyceae","Spumellaria","Nemertea")],
      #          label = c("Pelagophyceae","Spumellaria","Nemertea"),
      #          # label = bquote("Median"==.(format(median.adjR2.per.MEM[k],digits=2))),
      #          hjust = 0,
      #          size=9)
    
    g4[[i_case]] = ggplotGrob(fraction.ratio.vs.function + 
                                geom_hline(yintercept = 1, linetype="dashed") +
                                theme(axis.title.y = element_blank(),
                                      axis.text.x = element_text(vjust = 0.8, hjust = 0.7)))
  }
  
  # pdf(paste0(figure_folder,"/Fig4_pure+mixed_",div_threshold,"plusOTUs_PCoA.axis.1.x_selected100+1_apusozoans_(0).pdf"),
  #     # pdf(paste0(figure_folder,"/Fig4_dis.MEM_",div_threshold,"plusOTUs_PCoA.axes.x_selected100+1.pdf"),
  #     width=7.5*2.2/1.2/2*3,height=12/3*4/1.2/2)
  # for (i_case in 1:2)
  # {
  #   # grid.arrange(grobs = list(g1[[i_case]],g3[[i_case]],g4[[i_case]]),
  #   grid.arrange(grobs = list(g1_bis[[i_case]],g3[[i_case]],g4[[i_case]]),
  #                nrow = 1, ncol = 3)
  # }
  # dev.off()
  
  # pdf(paste0(figure_folder,"/Fig4_",div_threshold,"plusOTUs_PCoA.axes.x_selected100+1_no.selection_indiv.signif.axes.pdf"),
  # pdf(paste0(figure_folder,"/Fig4_pure+mixed_",div_threshold,"plusOTUs_PCoA.axes.x_selected100+1_indiv.signif.axes_no.var.selec_noPorifSpumell.pdf"),
  # pdf(paste0(figure_folder,"/Fig4_",div_threshold,"plusOTUs_PCoA.axes.x_selected100+1_indiv.signif.axes_no.var.selec_noSpumellNemert.pdf"),
  
  # pdf(paste0(figure_folder,"/Fig4_pure+mixed_",div_threshold,"plusOTUs_PCoA.axes.x_selected100+1_no.indiv.signif.axes_no.var.select_noRADC-MALVIV-Prasino-Rodo_(1).pdf"),
  pdf(paste0(figure_folder,"/Fig4_pure+mixed_",div_threshold,"plusOTUs_PCoA.axes.x_selected100+1_indiv.signif.axes_BH.correction_no.var.select_excluded.RAD-C.pdf"),
  # pdf(paste0(figure_folder,"/Fig4_pure_",div_threshold,"plusOTUs_PCoA.axes.x_selected100+1_indiv.signif.axes_BH.correction_no.var.select.pdf"),
  # pdf(paste0(figure_folder,"/Fig4_",div_threshold,"plusOTUs_PCoA.axes.x_selected100+1_no.indiv.signif.axes_no.var.select.pdf"),
  # pdf(paste0(figure_folder,"/Fig4_pure+mixed_",div_threshold,"plusOTUs_PCoA.axes.x_selected100+1_no.indiv.signif.axes_no.prelim.global.test.pdf"),
  # pdf(paste0(figure_folder,"/Fig4_",div_threshold,"plusOTUs_PCoA.axes.x_selected100+1_no.indiv.signif.axes_no.prelim.global.test.pdf"),
  # pdf(paste0(figure_folder,"/Fig4_pure+mixed_",div_threshold,"plusOTUs_PCoA.axes.x_selected100+1_var.selection.only.pdf"),
  # pdf(paste0(figure_folder,"/Fig4_",div_threshold,"plusOTUs_PCoA.axes.x_selected100+1_apusozans_(0).pdf"),
  # pdf(paste0(figure_folder,"/Fig4_dis.MEM_",div_threshold,"plusOTUs_PCoA.axes.x_selected100+1.pdf"),
      width=7.5*2.2/1.2,height=12/3*4/1.2)
  for (i_case in 1:2)
  {
    grid.arrange(grobs = list(g1[[i_case]],g2[[i_case]],g3[[i_case]],g4[[i_case]]),
                 # widths = rep(1,7),
                 # heights = c(0.8,0.8,1),
                 layout_matrix = rbind(c(1, 2),c(3, 4)))
  }
  dev.off()
  
  ##########################
  cor.test(NormalizedVI_pcoa[[3]][,1],
           colSums(varpart.env.spatial[[1]])[selected_groups],
           # method="kendall")
           method="spearman")
           # method="pearson")
  cor.test((colSums(varpart.env.spatial[[1]][2:3,])/(colSums(varpart.env.spatial[[1]][1:2,])))[selected_groups & !taxo_groups %in% "RAD-C"],
           NormalizedVI_pcoa[[3]][!taxo_groups[selected_groups] %in% "RAD-C",2],
           # method="kendall")
           method="spearman")
           # method="pearson")
  cor.test((colSums(varpart.env.spatial[[1]][2:3,])/(colSums(varpart.env.spatial[[1]][1:2,])))[selected_groups & !taxo_groups %in% "RAD-C"],
           size_relativeAbund[selected_groups & !taxo_groups %in% "RAD-C"],
           # method="kendall")
           method="spearman")
           # method="pearson")
  
  ############################
  excluded_groups = NULL
  ##### Individually selected axes only:
  # First outliers for pure+mixed fractions:
  excluded_groups = c("Porifera", "Spumellaria")
  # First outliers for pure fractions:
  excluded_groups = c("Spumellaria", "Nemertea")#,"MALV-III","Bicoecea","Ascomycota")
  ##### Individually selected axes only, with Benjamini-Hochberg correction:
  # "Porifera" and "Oomycota" are removed. 
  # First outliers for pure+mixed fractions - varfrac vs. log_size:
  excluded_groups = c("RAD-C", "Collodaria", "Vannellida") 
  # First outliers for pure+mixed fractions - varfrac vs. PCoA2:
  excluded_groups = c("RAD-C","Porifera")#, "Collodaria") 
  # First outliers for pure fractions:
  excluded_groups = c("Streptophyta")#,"Picomonadida")
  ##### Variable selection only:
  excluded_groups = "MALV-V"
  ##### No variable selection whatsoever (pure+mixed):
  excluded_groups = c("RAD-C", "MALV-IV")
  
  i_case = 1
  x  = log10(as.vector(diversity))[selected_groups & diversity>div_threshold
                                   # & !taxo_groups %in% c("Dinophyceae","Collodaria")
                                   & !taxo_groups %in% excluded_groups]
  x = log10(size_relativeAbund)[selected_groups & diversity>div_threshold
                                & !taxo_groups %in% c("Dinophyceae","Collodaria")
                                & !taxo_groups %in% excluded_groups]
  y = NormalizedVI_pcoa[[3]][,2][!names(NormalizedVI_pcoa[[3]][,2]) %in% excluded_groups]
                                 # & !names(NormalizedVI_pcoa[[3]][,2]) %in% c("Dinophyceae","Collodaria")] 
  y = log10((varpart.env.spatial[[i_case]][3,]/varpart.env.spatial[[i_case]][1,])[selected_groups & diversity>div_threshold
                                                                            & !colnames(varpart.env.spatial[[i_case]]) %in% c("Dinophyceae","Collodaria")
                                                                            & !colnames(varpart.env.spatial[[i_case]]) %in% excluded_groups])
  y = log10((colSums(varpart.env.spatial[[i_case]][2:3,])/colSums(varpart.env.spatial[[i_case]][1:2,]))[selected_groups & diversity>div_threshold
                                                                                                  & !colnames(varpart.env.spatial[[i_case]]) %in% c("Dinophyceae","Collodaria")
                                                                                                  & !colnames(varpart.env.spatial[[i_case]]) %in% excluded_groups])
  f = factor(dominant_function1[!names(dominant_function1) %in% excluded_groups],c("Phototrophs","Phagotrophs","Metazoans","Parasites"))
  f = factor(dominant_function1,c("Phagotrophs","Phototrophs","Metazoans","Parasites"))
  # f = factor(dominant_function1,c("Metazoans","Phagotrophs","Phototrophs","Parasites"))
  f1 = f2 = as.vector(f)
  # f1[f1 %in% c("Phagotrophs","Parasites")] = "Phago-Para"
  f2[f2 %in% c("Metazoans","Parasites")] = "Meta-Para"
  anova = lm(y[!is.infinite(y)] ~ f2[!is.infinite(y)])
  plot(anova)
  anova(anova)
  summary(anova)
  # TukeyHSD(anova)
  
  ############# Fig. 4 with abiotic only:
  g1 = g2 = g3 = g4 = list()
  for (i_case in 1)
  {
    x = NormalizedVI_pcoa[[3]][,1]
    y = colSums(varpart.env.spatial[[i_case]])[selected_groups & diversity>div_threshold]
    plot.tot.var.vs.PCoA.axis1 = cor.plot(x = x,
                                          y = y,
                                          y.lab = paste("Total variance explained\n by travel times and environment",c("","(DCM)")[i_case]),
                                          x.lab = "PCoA axis 1",
                                          y.lab.hjust = 0.5,
                                          x.log = F,
                                          y.log = F,
                                          fit = T,
                                          x.cor.pos = 0.8,
                                          y.cor.pos = 0.1,
                                          mar.vect = c(5,5,5,9))
    g1[[i_case]] = ggplotGrob(plot.tot.var.vs.PCoA.axis1)
    
    x = NormalizedVI_pcoa[[3]][,2]
    # y = (varpart.env.spatial[[i_case]][3,]/varpart.env.spatial[[i_case]][1,])[selected_groups & diversity>div_threshold]
    y = (colSums(varpart.env.spatial[[i_case]][2:3,])/(colSums(varpart.biotic.abiotic[[i_case]][1:2,])))[selected_groups & diversity>div_threshold]
    plot.connectivity.env.ratio.vs.PCoA.axis2 = cor.plot(x = x,
                                                         y = y,
                                                         # y.lab = paste("Ratio of variance purely explained\n by travel times over environment",c("","(DCM)")[i_case]),
                                                         y.lab = paste("Ratio of variance explained\n by travel times over abiotic environment",c("","(DCM)")[i_case]),
                                                         x.lab = "PCoA axis 2",
                                                         y.lab.hjust = 0.5,
                                                         x.log = F,
                                                         y.log = T,
                                                         fit = T,
                                                         x.cor.pos = 0.8,
                                                         y.cor.pos = 0.1,
                                                         mar.vect = c(5,5,5,5))
    g2[[i_case]] = ggplotGrob(plot.connectivity.env.ratio.vs.PCoA.axis2 + 
                                # geom_vline(xintercept = 1, linetype="dashed"))
                                geom_hline(yintercept = 1, linetype="dashed"))
    
    x = size_relativeAbund[selected_groups & diversity>div_threshold]
    # y = (varpart.env.spatial[[i_case]][3,]/varpart.env.spatial[[i_case]][1,])[selected_groups & diversity>div_threshold]
    y = (colSums(varpart.env.spatial[[i_case]][2:3,])/(colSums(varpart.biotic.abiotic[[i_case]][1:2,])))[selected_groups & diversity>div_threshold]
    fraction.ratio.vs.size = cor.plot(x = x,
                                      y = y,
                                      # y.lab = paste("Ratio of variance purely explained\n by travel times over environment",c("","(DCM)")[i_case]),
                                      y.lab = paste("Ratio of variance explained\n by travel times over abiotic environment",c("","(DCM)")[i_case]),
                                      x.lab = expression("Body size ("*mu*"m)"),
                                      y.lab.hjust = 0.5,
                                      x.log = T,
                                      y.log = T,
                                      fit = T,
                                      x.cor.pos = 0.8,
                                      y.cor.pos = 0.1,
                                      mar.vect = c(5,5,5,5))
    g3[[i_case]] = ggplotGrob(fraction.ratio.vs.size + geom_hline(yintercept = 1, linetype="dashed"))
    
    # values = (varpart.env.spatial[[i_case]][3,]/varpart.env.spatial[[i_case]][1,])[selected_groups & diversity>div_threshold & !colnames(varpart.env.spatial[[i_case]]) %in% c("Dinophyceae","Collodaria")]
    values = (colSums(varpart.env.spatial[[i_case]][2:3,])/(colSums(varpart.biotic.abiotic[[i_case]][1:2,])))[selected_groups & diversity>div_threshold & !colnames(varpart.env.spatial[[i_case]]) %in% c("Dinophyceae","Collodaria")]
    fraction.ratio.vs.function = functional.boxplot(classes = factor(dominant_function1, c("Phototrophs","Phagotrophs","Metazoans","Parasites")),
                                                    values = values,
                                                    # y.lab = paste("Ratio of variance purely explained\n by travel times over environment",c("","(DCM)")[i_case]),
                                                    y.lab = paste("Ratio of variance explained\n by travel times over abiotic environment",c("","(DCM)")[i_case]),
                                                    y.log = T,
                                                    y.lab.hjust = 0.01,
                                                    mar.vect = c(5,5,-7,5))
    g4[[i_case]] = ggplotGrob(fraction.ratio.vs.function + geom_hline(yintercept = 1, linetype="dashed"))
  }
  
  pdf(paste0(figure_folder,"/Fig4_pure+mixed_abiotic.env.only_",div_threshold,"plusOTUs_PCoA.axes.x_selected100+1.pdf"),
      # pdf(paste0(figure_folder,"/Fig4_",div_threshold,"plusOTUs_PCoA.axes.x_selected100+1.pdf"),
      # pdf(paste0(figure_folder,"/Fig4_dis.MEM_",div_threshold,"plusOTUs_PCoA.axes.x_selected100+1.pdf"),
      width=7.5*2.2/1.2,height=12/3*4/1.2)
  for (i_case in 1)
  {
    grid.arrange(grobs = list(g1[[i_case]],g2[[i_case]],g3[[i_case]],g4[[i_case]]),
                 # widths = rep(1,7),
                 # heights = c(0.8,0.8,1),
                 layout_matrix = rbind(c(1, 2),c(3, 4)))
  }
  dev.off()
}

###################################### Main SI figures:

# Mean sim vs. log diversity
{
  plot.mean.sim.vs.diversity =  cor.plot(y = mean_sim[selected_groups],
                                         x = as.vector(diversity)[selected_groups],
                                         # y.lab = paste0("Amount of biogeographic structure",ifelse(i_case == 1,""," (DCM)")),
                                         y.lab = "Mean similarity between posterior samples",
                                         x.lab = "Diversity (#OTUs), log scale",
                                         x.log = T,
                                         y.log = F,
                                         fit = F,
                                         # fit.display="pearson.spearman",
                                         x.cor.pos = 0.8,
                                         y.cor.pos = 0.1) +
    geom_vline(xintercept = 2000, linetype = "dashed")
  
  pdf(paste0(figure_folder,"/Fig.SI_mean.sim.vs.diversity.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*2/1.2,onefile=F)
  print(plot.mean.sim.vs.diversity)
  dev.off()
}

# Old Fig. S axes interpretation:
{
  i_case = 1
  
  y = I_square.observed_w.mean[selected_groups,i_case]
  x = NormalizedVI_pcoa[[3]][,1]
  plot.norm.VI.PCoA.axis1.vs.autocorr = cor.plot(x = x,
                                                 y = y,
                                                 # y.lab = paste0("Amount of biogeographic structure",ifelse(i_case == 1,""," (DCM)")),
                                                 y.lab = "Short-distance spatial autocorrelation",
                                                 x.lab = "PCoA axis 1",
                                                 x.log = F,
                                                 y.log = F,
                                                 fit = T,
                                                 fit.display="pearson.spearman",
                                                 x.cor.pos = 0.8,
                                                 y.cor.pos = 0.1)
  
  y = charac_scale[selected_groups,i_case]
  x = NormalizedVI_pcoa[[3]][,2]
  plot.norm.VI.PCoA.axis2.vs.autocorr.scale = cor.plot(x = x,
                                                       y = y,
                                                       # y.lab = paste(ifelse(i_case == 1,"Scale of biogeographic organization (km)",
                                                       #                      "Scale of biogeo. organization(km; DCM)")),
                                                       y.lab = "Scale of biogeographic organization (km)",
                                                       x.lab = "PCoA axis 2",
                                                       x.log = F,
                                                       y.log = F,
                                                       fit = T,
                                                       fit.display="pearson.spearman",
                                                       x.cor.pos = 0.2,
                                                       y.cor.pos = 0.85)
  
  y = basin_I[selected_groups & diversity>div_threshold,i_case]/I_square.observed_w.mean[selected_groups & diversity>div_threshold,i_case]
  # y = basin_I[selected_groups,i_case]
  x = NormalizedVI_pcoa[[3]][,2]
  plot.norm.VI.PCoA.axis2.vs.basin.I.within = cor.plot(x = x,
                                                       y = y,
                                                       # y.lab = paste0("Within-basin homogeneity",ifelse(i_case == 1,""," (DCM)")),
                                                       y.lab = "Within-basin homogeneity",
                                                       x.lab = "PCoA axis 2",
                                                       x.log = F,
                                                       y.log = F,
                                                       fit = T,
                                                       fit.display="pearson.spearman",
                                                       x.cor.pos = 0.8,
                                                       y.cor.pos = 0.1)
  
  # y = basin_I_between[selected_groups & diversity>div_threshold,i_case]/I_square.observed_w.mean[selected_groups & diversity>div_threshold,i_case]
  # # x = basin_I[selected_groups & diversity>div_threshold,1]
  # x = NormalizedVI_pcoa[[3]][,2]
  # plot.norm.VI.PCoA.axis2.vs.basin.I.between = cor.plot(x = x,
  #                                                       y = y,
  #                                                       y.lab = "Between-basin similarity",
  #                                                       x.lab = "PCoA axis 2",
  #                                                       x.log = F,
  #                                                       y.log = F,
  #                                                       fit = T,
  #                                                       x.cor.pos = 0.8,
  #                                                       y.cor.pos = 0.1)
  # 
  # y = basin_I_contrast[selected_groups & diversity>div_threshold,i_case]/I_square.observed_w.mean[selected_groups & diversity>div_threshold,i_case]
  # # x = basin_I_contrast[selected_groups & diversity>div_threshold,1]
  # x = NormalizedVI_pcoa[[3]][,2]
  # plot.norm.VI.PCoA.axis2.vs.basin.I.contrast = cor.plot(x = x,
  #                                                        y = y,
  #                                                        y.lab = "Basin structure",
  #                                                        x.lab = "PCoA axis 2",
  #                                                        x.log = F,
  #                                                        y.log = F,
  #                                                        fit = T,
  #                                                        x.cor.pos = 0.2,
  #                                                        y.cor.pos = 0.1)
  
  y = lat_I[selected_groups & diversity>div_threshold,i_case]/I_square.observed_w.mean[selected_groups & diversity>div_threshold,i_case]
  # y = lat_I[selected_groups,i_case]
  x = NormalizedVI_pcoa[[3]][,2]
  plot.norm.VI.PCoA.axis2.vs.lat.I.sym = cor.plot(x = x,
                                                  y = y,
                                                  # y.lab = paste0("Latitudinal symmetry",ifelse(i_case == 1,""," (DCM)")),
                                                  y.lab = "Latitudinal symmetry",
                                                  x.lab = "PCoA axis 2",
                                                  x.log = F,
                                                  y.log = F,
                                                  fit = T,
                                                  fit.display="pearson.spearman",
                                                  x.cor.pos = 0.85,
                                                  y.cor.pos = 0.85)
  
  # y = lat_I_nonsym[selected_groups & diversity>div_threshold,i_case]/I_square.observed_w.mean[selected_groups & diversity>div_threshold,i_case]
  # # x = lat_I_nonsym[selected_groups & diversity>div_threshold,1]
  # x = NormalizedVI_pcoa[[3]][,2]
  # plot.norm.VI.PCoA.axis2.vs.lat.I.nonsym = cor.plot(x = x,
  #                                                    y = y,
  #                                                    y.lab = "Latitudinal homogeneity",
  #                                                    x.lab = "PCoA axis 2",
  #                                                    x.log = F,
  #                                                    y.log = F,
  #                                                    fit = T,
  #                                                    x.cor.pos = 0.8,
  #                                                    y.cor.pos = 0.1)
  # 
  # y = lat_I_all[selected_groups & diversity>div_threshold,i_case]/I_square.observed_w.mean[selected_groups & diversity>div_threshold,i_case]
  # # x = lat_I_all[selected_groups & diversity>div_threshold,1]
  # x = NormalizedVI_pcoa[[3]][,2]
  # plot.norm.VI.PCoA.axis2.vs.lat.I.all = cor.plot(x = x,
  #                                                 y = y,
  #                                                 y.lab = "Latitudinal structure",
  #                                                 x.lab = "PCoA axis 2",
  #                                                 x.log = F,
  #                                                 y.log = F,
  #                                                 fit = T,
  #                                                 x.cor.pos = 0.8,
  #                                                 y.cor.pos = 0.1)
  
  pdf(paste0(figure_folder,"/Fig.SI_axes.interpretation_normalized_",ifelse(i_case==1,"","DCM_"),div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1_apusozoans.pdf"),
      width=7.5*2.2/1.2*2/2,height=12/3*4/1.2,onefile=F)
  #par(mfrow=c(5,3))
  # layout(mat=rbind(c(1, 1, 2),c(3, 3, 4),c(5, 6, 7),c(8, 8, 9),c(10, 10, 11)))
  # plot(plot.map.group[[1]])
  # print(plot.norm.VI.PCoA.axis2.vs.body.size)
  # plot(plot.map.group[[2]],newpage=T)
  # plot(plot.norm.VI.PCoA.axis2.vs.autocorr.scale,newpage=T)
  # plot(plot.norm.VI.PCoA.axis1.vs.diversity,newpage=T)
  # plot(plot.norm.VI.PCoA.axis1.vs.autocorr,newpage=T)
  # plot(plot.VI.PCoA.biplot,newpage=T)
  # plot(plot.map.group[[3]],newpage=T)
  # plot(plot.norm.VI.PCoA.axis2.vs.basin.I,newpage=T)
  # plot(plot.map.group[[4]],newpage=T)
  # plot(plot.norm.VI.PCoA.axis2.vs.lat.I,newpage=T)
  # g1 = ggplotGrob(plot.map.group[[1]])
  g1 = ggplotGrob(plot.norm.VI.PCoA.axis1.vs.autocorr)
  g2 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.autocorr.scale)
  g3 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.basin.I.within)
  # g4 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.basin.I.between)
  # g5 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.basin.I.contrast)
  g6 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.lat.I.sym)
  # g7 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.lat.I.nonsym)
  # g8 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.lat.I.all)
  # g8 = ggplotGrob(plot.map.group[[1]])
  # g9 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.basin.I)
  # g10 = ggplotGrob(plot.map.group[[1]])
  # g11 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.lat.I)
  # grid.draw(arrangeGrob(grobs = list(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11),
  grid.arrange(grobs = list(g1,g2,g3,g6),
               # widths = rep(1,7),
               # heights = c(0.8,0.8,1),
               layout_matrix = rbind(c(1, 2),c(3, 4)),
               respect=T)
  # grid.draw(arrangeGrob(grobs = list(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11),
  #                       layout_matrix = rbind(c(1, 1, 2),c(3, 3, 4),c(5, 6, 7),c(8, 8, 9),c(10, 10, 11)),
  #                       respect=T))
  dev.off()
}

# Fig. S axis 1 interpretation:
{
  plot.pcoa1.vs.autocorr = list()
  plot.pcoa1.vs.shannon = list()
  plot.pcoa1.vs.shannon.total = list()
  plot.pcoa1.vs.nb.dominants = list()
  for (i_case in 1:2)
  {
    plot.pcoa1.vs.autocorr[[i_case]] = cor.plot(y = NormalizedVI_pcoa[[3]][,1],
                                                x = I_square.observed_w.mean[selected_groups,i_case],
                                                # y.lab = paste0("Amount of biogeographic structure",ifelse(i_case == 1,""," (DCM)")),
                                                y.lab = "Biogeographic axis 1",
                                                x.lab = "Short-distance spatial autocorrelation",
                                                x.log = F,
                                                y.log = F,
                                                fit = T,
                                                fit.display="pearson.p",
                                                x.cor.pos = 0.8,
                                                y.cor.pos = 0.1) + 
      ylim(c(-0.3,0.3)) 
      # xlim(min(I_square.observed_w.mean[selected_groups,1],I_square.observed_w.mean[selected_groups,2]),
      #      max(I_square.observed_w.mean[selected_groups,1],I_square.observed_w.mean[selected_groups,2]))
    
    plot.pcoa1.vs.shannon[[i_case]] = cor.plot(y = NormalizedVI_pcoa[[3]][,1],
                                               x = if (i_case == 1) shannon_SUR[selected_groups] 
                                               else if (i_case == 2) shannon_DCM[selected_groups],
                                               # excluded.points = which(taxo_groups[selected_groups] == "Porifera"),
                                               y.lab = "Biogeographic axis 1",
                                               x.lab = "Mean Shannon entropy of\nassemblage proportions across stations",
                                               y.lab.hjust = 0.5,
                                               x.log = F,
                                               y.log = F,
                                               fit = T,
                                               fit.display="pearson.p",
                                               x.cor.pos = 0.8,
                                               y.cor.pos = 0.8,
                                               mar.vect = c(5,5,1,5)) + 
      ylim(c(-0.3,0.3)) 
      # xlim(min(shannon_SUR[selected_groups],shannon_DCM[selected_groups]),
      #      max(shannon_SUR[selected_groups],shannon_DCM[selected_groups]))
    
    plot.pcoa1.vs.shannon.total[[i_case]] = cor.plot(y = NormalizedVI_pcoa[[3]][,1],
                                                     x = if (i_case == 1) shannon_total_SUR[selected_groups]
                                                     else if  (i_case == 2) shannon_total_DCM[selected_groups],
                                                     y.lab = "Biogeographic axis 1",
                                                     x.lab = "Overall Shannon entropy of\nassemblage proportions in the biogeography     ",
                                                     y.lab.hjust = 0.5,
                                                     x.log = F,
                                                     y.log = F,
                                                     fit = T,
                                                     fit.display="pearson.p",
                                                     x.cor.pos = 0.8,
                                                     y.cor.pos = 1.1,
                                                     mar.vect = c(5,5,1,5)) +
      ylim(c(-0.3,0.3)) 
      # xlim(min(shannon_total_SUR[selected_groups],shannon_total_DCM[selected_groups]),
      #      max(shannon_total_SUR[selected_groups],shannon_total_DCM[selected_groups]))
    
    # plot.pcoa1.vs.nb.dominants[[i_case]] = cor.plot(y = NormalizedVI_pcoa[[3]][,1],
    #                                                 x = nb_dominants[selected_groups,i_case],
    #                                                 y.lab = "PCoA axis 1",
    #                                                 x.lab = "Number of dominant assemblages\n in the biogeography",
    #                                                 y.lab.hjust = 0.5,
    #                                                 x.log = F,
    #                                                 y.log = F,
    #                                                 fit = T,
    #                                                 fit.display="pearson.p",
    #                                                 x.cor.pos = 0.8,
    #                                                 y.cor.pos = 1.1,
    #                                                 mar.vect = c(5,5,1,5)) + 
    #   ylim(c(-0.3,0.3)) +
    #   xlim(min(nb_dominants[selected_groups,1],nb_dominants[selected_groups,2]),
    #        max(nb_dominants[selected_groups,1],nb_dominants[selected_groups,2]))
  }
  
  # plot.pcoa1.vs.nb.absolute.dominants = cor.plot(x = nb_absolute_dominants[selected_groups,1],
  #                                                y = NormalizedVI_pcoa[[3]][,1],
  #                                                x.lab = "Number of dominant (>50%) assemblages",
  #                                                y.lab = "PCoA axis 1",
  #                                                y.lab.hjust = 0.5,
  #                                                x.log = F,
  #                                                y.log = F,
  #                                                fit = T,
  #                                                x.cor.pos = 0.8,
  #                                                y.cor.pos = 0.8,
  #                                                mar.vect = c(5,5,1,5))
  
  pdf(paste0(figure_folder,"/Fig.SI_axis.1.interpretation_",ifelse(i_case==1,"","DCM_"),div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1_apusozoans.pdf"),
      width=7.5*2.2/1.2*2/2,height=12/3*2/1.2*2,onefile=F)
  g1 = ggplotGrob(plot.pcoa1.vs.autocorr)
  g2 = ggplotGrob(plot.pcoa1.vs.shannon)
  g3 = ggplotGrob(plot.pcoa1.vs.shannon.total)
  g4 = ggplotGrob(plot.pcoa1.vs.nb.dominants)
  grid.arrange(grobs = list(g1,g2,g3,g4),
               layout_matrix = rbind(c(1, 2),c(3, 4)),
               respect=T)
  dev.off()
  
  pdf(paste0(figure_folder,"/Fig.SI_axis.1.interpretation_SUR.DCM",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2*2/2,height=12/3*2/1.2*3,onefile=F)
  g1 = ggplotGrob(plot.pcoa1.vs.autocorr[[1]])
  g2 = ggplotGrob(plot.pcoa1.vs.autocorr[[2]] + theme(axis.title.y=element_blank()))
  g3 = ggplotGrob(plot.pcoa1.vs.shannon[[1]])
  g4 = ggplotGrob(plot.pcoa1.vs.shannon[[2]] + theme(axis.title.y=element_blank()))
  g5 = ggplotGrob(plot.pcoa1.vs.shannon.total[[1]])
  g6 = ggplotGrob(plot.pcoa1.vs.shannon.total[[2]] + theme(axis.title.y=element_blank()))
  grid.arrange(grobs = list(g1,g2,g3,g4,g5,g6),
               layout_matrix = rbind(c(1, 2),c(3, 4),c(5, 6)),
               respect=T)
  dev.off()
}

# Fig. S axis 2 interpretation:
{
  # i_case = 1
  plot.norm.VI.PCoA.axis2.vs.autocorr.scale = list()
  plot.norm.VI.PCoA.axis2.vs.basin.I.within = list()
  plot.norm.VI.PCoA.axis2.vs.lat.I.sym = list()
  for (i_case in 1:2)
  {
    plot.norm.VI.PCoA.axis2.vs.autocorr.scale[[i_case]] = cor.plot(x = charac_scale[selected_groups,i_case],
                                                                   y = NormalizedVI_pcoa[[3]][,2],
                                                                   # y.lab = paste(ifelse(i_case == 1,"Scale of biogeographic organization (km)",
                                                                   #                      "Scale of biogeo. organization(km; DCM)")),
                                                                   x.lab = "Scale of biogeographic organization (km)",
                                                                   y.lab = "Biogeographic axis 2",
                                                                   x.log = F,
                                                                   y.log = F,
                                                                   fit = T,
                                                                   fit.display="pearson.p",
                                                                   x.cor.pos = 0.8,
                                                                   y.cor.pos = 0.15,
                                                                   mar.vect=c(5,5,1,5)) + 
      ylim(c(-0.25,0.25)) 
      # xlim(min(charac_scale[selected_groups,1],charac_scale[selected_groups,2]),
      #      max(charac_scale[selected_groups,1],charac_scale[selected_groups,2]))
    
    plot.norm.VI.PCoA.axis2.vs.basin.I.within[[i_case]] = cor.plot(x = basin_I[selected_groups & diversity>div_threshold,i_case]/I_square.observed_w.mean[selected_groups & diversity>div_threshold,i_case],
                                                                   y = NormalizedVI_pcoa[[3]][,2],
                                                                   # y.lab = paste0("Within-basin homogeneity",ifelse(i_case == 1,""," (DCM)")),
                                                                   x.lab = "Within-basin autocorrelation",
                                                                   y.lab = "Biogeographic axis 2",
                                                                   x.log = F,
                                                                   y.log = F,
                                                                   fit = T,
                                                                   fit.display="pearson.p",
                                                                   x.cor.pos = 0.2,
                                                                   y.cor.pos = 0.8,
                                                                   mar.vect=c(5,5,1,5)) + 
      ylim(c(-0.25,0.25)) 
      # xlim(min(basin_I[selected_groups & diversity>div_threshold,1]/I_square.observed_w.mean[selected_groups & diversity>div_threshold,1],
      #          basin_I[selected_groups & diversity>div_threshold,2]/I_square.observed_w.mean[selected_groups & diversity>div_threshold,2]),
      #      max(basin_I[selected_groups & diversity>div_threshold,1]/I_square.observed_w.mean[selected_groups & diversity>div_threshold,1],
      #          basin_I[selected_groups & diversity>div_threshold,2]/I_square.observed_w.mean[selected_groups & diversity>div_threshold,2]))
    
    plot.norm.VI.PCoA.axis2.vs.lat.I.sym[[i_case]] = cor.plot(x = lat_I[selected_groups & diversity>div_threshold,i_case]/I_square.observed_w.mean[selected_groups & diversity>div_threshold,i_case],
                                                              y = NormalizedVI_pcoa[[3]][,2],
                                                              # y.lab = paste0("Latitudinal symmetry",ifelse(i_case == 1,""," (DCM)")),
                                                              x.lab = "Latitudinal autocorrelation",
                                                              y.lab = "Biogeographic axis 2",
                                                              x.log = F,
                                                              y.log = F,
                                                              fit = T,
                                                              fit.display="pearson.p",
                                                              x.cor.pos = 0.85,
                                                              y.cor.pos = 1,
                                                              mar.vect=c(5,5,1,5)) + 
      ylim(c(-0.25,0.25))
      # xlim(min(lat_I[selected_groups & diversity>div_threshold,1]/I_square.observed_w.mean[selected_groups & diversity>div_threshold,1],
      #          lat_I[selected_groups & diversity>div_threshold,2]/I_square.observed_w.mean[selected_groups & diversity>div_threshold,2]),
      #      max(lat_I[selected_groups & diversity>div_threshold,1]/I_square.observed_w.mean[selected_groups & diversity>div_threshold,1],
      #          lat_I[selected_groups & diversity>div_threshold,2]/I_square.observed_w.mean[selected_groups & diversity>div_threshold,2]))
  }
  
  pdf(paste0(figure_folder,"/Fig.SI_axis.2.interpretation_",ifelse(i_case==1,"","DCM_"),div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2*2/2,height=12/3*2/1.2*2,onefile=F)
  g1 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.autocorr.scale[[i_case]])
  g2 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.basin.I.within[[i_case]])
  g3 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.lat.I.sym[[i_case]])
  grid.arrange(grobs = list(g1,g2,g3),
               layout_matrix = rbind(c(1, 2),c(3, NA)),
               respect=T)
  dev.off()
  
  pdf(paste0(figure_folder,"/Fig.SI_axis.2.interpretation_SUR.DCM",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2*2/2,height=12/3*2/1.2*3,onefile=F)
  g1 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.autocorr.scale[[1]])
  g2 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.autocorr.scale[[2]] + theme(axis.title.y=element_blank()))
  g3 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.basin.I.within[[1]])
  g4 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.basin.I.within[[2]] + theme(axis.title.y=element_blank()))
  g5 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.lat.I.sym[[1]])
  g6 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.lat.I.sym[[2]] + theme(axis.title.y=element_blank()))
  grid.arrange(grobs = list(g1,g2,g3,g4,g5,g6),
               layout_matrix = rbind(c(1, 2),c(3, 4),c(5, 6)),
               respect=T)
  dev.off()
}

# Fig. S axis 2 interpretation - residuals:
{
  # i_case = 1
  plot.norm.VI.PCoA.axis2.vs.res.autocorr.scale = list()
  plot.norm.VI.PCoA.axis2.vs.res.basin.I.within = list()
  plot.norm.VI.PCoA.axis2.vs.res.lat.I.sym = list()
  for (i_case in 1:2)
  {
    plot.norm.VI.PCoA.axis2.vs.res.autocorr.scale[[i_case]] = cor.plot(x = lm(charac_scale[selected_groups,i_case] ~ NormalizedVI_pcoa[[3]][,1])$residuals,
                                                                       y = NormalizedVI_pcoa[[3]][,2],
                                                                       # y.lab = paste(ifelse(i_case == 1,"Scale of biogeographic organization (km)",
                                                                       #                      "Scale of biogeo. organization(km; DCM)")),
                                                                       x.lab = "Scale of biogeographic organization\n- residuals after regression on PCoA1",
                                                                       y.lab = "PCoA axis 2",
                                                                       x.log = F,
                                                                       y.log = F,
                                                                       fit = T,
                                                                       fit.display="pearson.p",
                                                                       x.cor.pos = 0.8,
                                                                       y.cor.pos = 0.15,
                                                                       mar.vect=c(5,5,1,5)) + 
      ylim(c(-0.25,0.25)) 
      # xlim(min(lm(charac_scale[selected_groups,1] ~ NormalizedVI_pcoa[[3]][,1])$residuals,
      #          lm(charac_scale[selected_groups,2] ~ NormalizedVI_pcoa[[3]][,1])$residuals),
      #      max(lm(charac_scale[selected_groups,1] ~ NormalizedVI_pcoa[[3]][,1])$residuals,
      #          lm(charac_scale[selected_groups,2] ~ NormalizedVI_pcoa[[3]][,1])$residuals))
    
    plot.norm.VI.PCoA.axis2.vs.res.basin.I.within[[i_case]] = cor.plot(x = lm(basin_I[selected_groups,i_case]/I_square.observed_w.mean[selected_groups,i_case] ~ NormalizedVI_pcoa[[3]][,1])$residuals,
                                                                       y = NormalizedVI_pcoa[[3]][,2],
                                                                       # y.lab = paste0("Within-basin homogeneity",ifelse(i_case == 1,""," (DCM)")),
                                                                       x.lab = "Within-basin homogeneity\n- residuals after regression on PCoA1",
                                                                       y.lab = "PCoA axis 2",
                                                                       x.log = F,
                                                                       y.log = F,
                                                                       fit = T,
                                                                       fit.display="pearson.p",
                                                                       x.cor.pos = 0.2,
                                                                       y.cor.pos = 0.8,
                                                                       mar.vect=c(5,5,1,5)) + 
      ylim(c(-0.25,0.25)) 
      # xlim(min(lm(basin_I[selected_groups,1]/I_square.observed_w.mean[selected_groups,1] ~ NormalizedVI_pcoa[[3]][,1])$residuals,
      #          lm(basin_I[selected_groups,2]/I_square.observed_w.mean[selected_groups,2] ~ NormalizedVI_pcoa[[3]][,1])$residuals),
      #      max(lm(basin_I[selected_groups,1]/I_square.observed_w.mean[selected_groups,1] ~ NormalizedVI_pcoa[[3]][,1])$residuals,
      #          lm(basin_I[selected_groups,2]/I_square.observed_w.mean[selected_groups,2] ~ NormalizedVI_pcoa[[3]][,1])$residuals))
    
    plot.norm.VI.PCoA.axis2.vs.res.lat.I.sym[[i_case]] = cor.plot(x = lm(lat_I[selected_groups,i_case]/I_square.observed_w.mean[selected_groups,i_case] ~ NormalizedVI_pcoa[[3]][,1])$residuals,
                                                                  y = NormalizedVI_pcoa[[3]][,2],
                                                                  # y.lab = paste0("Latitudinal symmetry",ifelse(i_case == 1,""," (DCM)")),
                                                                  x.lab = "Latitudinal symmetry\n- residuals after regression on PCoA1",
                                                                  y.lab = "PCoA axis 2",
                                                                  x.log = F,
                                                                  y.log = F,
                                                                  fit = T,
                                                                  fit.display="pearson.p",
                                                                  x.cor.pos = 0.85,
                                                                  y.cor.pos = 1,
                                                                  mar.vect=c(5,5,1,5)) + 
      ylim(c(-0.25,0.25)) 
      # xlim(min(lm(lat_I[selected_groups,1]/I_square.observed_w.mean[selected_groups,1] ~ NormalizedVI_pcoa[[3]][,1])$residuals,
      #          lm(lat_I[selected_groups,2]/I_square.observed_w.mean[selected_groups,2] ~ NormalizedVI_pcoa[[3]][,1])$residuals),
      #      max(lm(lat_I[selected_groups,1]/I_square.observed_w.mean[selected_groups,1] ~ NormalizedVI_pcoa[[3]][,1])$residuals,
      #          lm(lat_I[selected_groups,2]/I_square.observed_w.mean[selected_groups,2] ~ NormalizedVI_pcoa[[3]][,1])$residuals))
  }
  
  # pdf(paste0(figure_folder,"/Fig.SI_axis.2.interpretation_",ifelse(i_case==1,"","DCM_"),div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
  #     width=7.5*2.2/1.2*2/2,height=12/3*2/1.2*2,onefile=F)
  # g1 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.autocorr.scale[[i_case]])
  # g2 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.basin.I.within[[i_case]])
  # g3 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.lat.I.sym[[i_case]])
  # grid.arrange(grobs = list(g1,g2,g3),
  #              layout_matrix = rbind(c(1, 2),c(3, NA)),
  #              respect=T)
  # dev.off()
  
  pdf(paste0(figure_folder,"/Fig.SI_axis.2.residuals.interpretation_SUR.DCM",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2*2/2,height=12/3*2/1.2*3,onefile=F)
  g1 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.res.autocorr.scale[[1]])
  g2 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.res.autocorr.scale[[2]] + theme(axis.title.y=element_blank()))
  g3 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.res.basin.I.within[[1]])
  g4 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.res.basin.I.within[[2]] + theme(axis.title.y=element_blank()))
  g5 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.res.lat.I.sym[[1]])
  g6 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.res.lat.I.sym[[2]] + theme(axis.title.y=element_blank()))
  grid.arrange(grobs = list(g1,g2,g3,g4,g5,g6),
               layout_matrix = rbind(c(1, 2),c(3, 4),c(5, 6)),
               respect=T)
  dev.off()
}

# Fig. S biogeographic descriptors vs diversity
{
  i_case = 1
  plot.moran.I.vs.log.diversity = list()
  plot.shannon.vs.log.diversity = list()
  plot.shannon.total.vs.log.diversity = list()
  # for (i_case in 1)
  # {
  y = I_square.observed_w.mean[selected_groups,i_case]
  x = as.vector(diversity)[selected_groups]
  cor.test = cor.test(log10(x[x<2000]),y[x<2000],na.rm=T,method="spearman")
  x.cor.pos = 0.8
  y.cor.pos = 0.15
  plot.moran.I.vs.log.diversity[[i_case]] = cor.plot(y = y,
                                                     x = x,
                                                     # excluded.points = which(taxo_groups[selected_groups] %in% "RAD-A"),
                                                     y.lab = "Surface short-distance spatial autocorr.",
                                                     x.lab = "Diversity (#OTUs), log scale",
                                                     x.log = T,
                                                     y.log = F,
                                                     fit = F,
                                                     mar.vect=c(5,10,5,5)) +
    geom_smooth(data = data.frame(x = x[x<2000], y = y[x<2000]), aes(x,y),method='lm',col="black",se=T) +
    # geom_smooth(aes(x[x<2000], y[x<2000]),method='lm',col="black",se=T, inherit.aes = T) +
    geom_vline(xintercept = 2000, linetype="dashed") +
    annotate(geom="text", 
             x=(max(x)/min(x))^x.cor.pos*min(x),
             y=y.cor.pos*(max(y)-min(y))+min(y),
             label=bquote(atop(rho[S]==.(format(cor.test$estimate,digits=2,nsmall=2)),P==.(format(cor.test$p.value,digits=1)))),
             size=7.5) 
  
  y = if (i_case == 1) shannon_SUR[selected_groups] else if (i_case == 2) shannon_DCM[selected_groups]
  x = as.vector(diversity)[selected_groups]
  cor.test = cor.test(log10(x[x<2000]),y[x<2000],na.rm=T,method="spearman")
  x.cor.pos = 0.8
  y.cor.pos = 0.7
  plot.shannon.vs.log.diversity[[i_case]] = cor.plot(y = y,
                                                     x = x,
                                                     excluded.points = which(taxo_groups[selected_groups] %in% "Porifera"),
                                                     y.lab = "Mean assemblage Shannon entropy\n across surface stations",
                                                     x.lab = "Diversity (#OTUs), log scale",
                                                     x.log = T,
                                                     y.log = F,
                                                     fit = F) +
    geom_smooth(data = data.frame(x = x[x<2000], y = y[x<2000]), aes(x,y),method='lm',col="black",se=T) +
    # geom_smooth(aes(x[x<2000], y[x<2000]),method='lm',col="black",se=T, inherit.aes = T) +
    geom_vline(xintercept = 2000, linetype="dashed") +
    annotate(geom="text", 
             x=(max(x)/min(x))^x.cor.pos*min(x),
             y=y.cor.pos*(max(y)-min(y))+min(y),
             label=bquote(atop(rho[S]==.(format(cor.test$estimate,digits=2,nsmall=2)),P==.(format(cor.test$p.value,digits=1)))),
             size=7.5) 
  
  y = if (i_case == 1) shannon_total_SUR[selected_groups] else if  (i_case == 2) shannon_total_DCM[selected_groups]
  x = as.vector(diversity)[selected_groups]
  cor.test = cor.test(log10(x[x<2000]),y[x<2000],na.rm=T,method="spearman")
  x.cor.pos = 0.8
  y.cor.pos = 0.8
  plot.shannon.total.vs.log.diversity[[i_case]] = cor.plot(y = y,
                                                           x = x,
                                                           # excluded.points = which(taxo_groups[selected_groups] %in% "RAD-A"),
                                                           y.lab = "Overall assemblage Shannon entropy\n in the surface biogeography",
                                                           x.lab = "Diversity (#OTUs), log scale",
                                                           x.log = T,
                                                           y.log = F,
                                                           fit = F,
                                                           mar.vect=c(5,10,5,5)) +
    geom_smooth(data = data.frame(x = x[x<2000], y = y[x<2000]), aes(x,y),method='lm',col="black",se=T) +
    # geom_smooth(aes(x[x<2000], y[x<2000]),method='lm',col="black",se=T, inherit.aes = T) +
    geom_vline(xintercept = 2000, linetype="dashed") +
    annotate(geom="text", 
             x=(max(x)/min(x))^x.cor.pos*min(x),
             y=y.cor.pos*(max(y)-min(y))+min(y),
             label=bquote(atop(rho[S]==.(format(cor.test$estimate,digits=2,nsmall=2)),P==.(format(cor.test$p.value,digits=1)))),
             size=7.5) 
  # }
  pdf(paste0(figure_folder,"/Fig.SI_biographic.descriptors.vs.log.diversity_",ifelse(i_case==1,"","DCM_"),".pdf"),
      width=7.5*2.2/1.2*2/2,height=12/3*2/1.2*2,onefile=F)
  g1 = ggplotGrob(plot.moran.I.vs.log.diversity[[i_case]])
  g2 = ggplotGrob(plot.shannon.vs.log.diversity[[i_case]])
  g3 = ggplotGrob(plot.shannon.total.vs.log.diversity[[i_case]])
  grid.arrange(grobs = list(g1,g2,g3),
               layout_matrix = rbind(c(1, 2),c(3, NA)),
               respect=T)
  dev.off()
}

# Fig. S biogeographic descriptors vs body size - ecological function:
{
  library(ppcor)
  
  ##############
  
  i_case = 1
  plot.res.autocorr.scale.vs.res.log.body.size = list()
  plot.res.basin.I.within.vs.res.log.body.size = list()
  plot.res.lat.I.sym.vs.res.log.body.size = list()
  # for (i_case in 1)
  # {
  plot.res.autocorr.scale.vs.res.log.body.size[[i_case]] = cor.plot(y = lm(charac_scale[selected_groups,i_case] ~ NormalizedVI_pcoa[[3]][,1])$residuals,
                                                                    x = lm(log(size_relativeAbund)[selected_groups] ~ NormalizedVI_pcoa[[3]][,1])$residuals,
                                                                    # y.lab = paste(ifelse(i_case == 1,"Scale of biogeographic organization (km)",
                                                                    #                      "Scale of biogeo. organization(km; DCM)")),
                                                                    y.lab = "Scale of biogeographic organization\n- residuals after regression on PCoA1",
                                                                    x.lab = "Log body size\n- residuals after regression on PCoA1",
                                                                    x.log = F,
                                                                    y.log = F,
                                                                    fit = T,
                                                                    fit.display="pearson.spearman",
                                                                    x.cor.pos = 0.8,
                                                                    y.cor.pos = 0.15)
  
  plot.res.basin.I.within.vs.res.log.body.size[[i_case]] = cor.plot(y = lm(basin_I[selected_groups,i_case]/I_square.observed_w.mean[selected_groups,i_case] ~ NormalizedVI_pcoa[[3]][,1])$residuals,
                                                                    x = lm(log(size_relativeAbund[selected_groups]) ~ NormalizedVI_pcoa[[3]][,1])$residuals,
                                                                    # y.lab = paste0("Within-basin homogeneity",ifelse(i_case == 1,""," (DCM)")),
                                                                    excluded_groups = "RAD-A",
                                                                    y.lab = "Within-basin homogeneity\n- residuals after regression on PCoA1",
                                                                    x.lab = "Log body size\n- residuals after regression on PCoA1",
                                                                    x.log = F,
                                                                    y.log = F,
                                                                    fit = T,
                                                                    fit.display="pearson.spearman",
                                                                    x.cor.pos = 0.8,
                                                                    y.cor.pos = 0.1)
  
  plot.res.lat.I.sym.vs.res.log.body.size[[i_case]] = cor.plot(y = lm(lat_I[selected_groups,i_case]/I_square.observed_w.mean[selected_groups,i_case] ~ NormalizedVI_pcoa[[3]][,1])$residuals,
                                                               x = lm(log(size_relativeAbund[selected_groups]) ~ NormalizedVI_pcoa[[3]][,1])$residuals,
                                                               # y.lab = paste0("Latitudinal symmetry",ifelse(i_case == 1,""," (DCM)")),
                                                               y.lab = "Latitudinal symmetry\n- residuals after regression on PCoA1",
                                                               x.lab = "Log body size\n- residuals after regression on PCoA1",
                                                               x.log = F,
                                                               y.log = F,
                                                               fit = T,
                                                               fit.display="pearson.spearman",
                                                               x.cor.pos = 0.8,
                                                               y.cor.pos = 0.8)
  # }
  pdf(paste0(figure_folder,"/Fig.SI_residuals.biographic.descriptors.vs.residuals.log.body.size_",ifelse(i_case==1,"","DCM_"),div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2*2/2,height=12/3*2/1.2*2,onefile=F)
  g1 = ggplotGrob(plot.res.autocorr.scale.vs.res.log.body.size[[i_case]])
  g2 = ggplotGrob(plot.res.basin.I.within.vs.res.log.body.size[[i_case]])
  g3 = ggplotGrob(plot.res.lat.I.sym.vs.res.log.body.size[[i_case]])
  grid.arrange(grobs = list(g1,g2,g3),
               layout_matrix = rbind(c(1, 2),c(3, NA)),
               respect=T)
  dev.off()
  
  ##############
  i_case = 1
  plot.res.autocorr.scale.vs.log.body.size = list()
  plot.res.basin.I.within.vs.log.body.size = list()
  plot.res.lat.I.sym.vs.log.body.size = list()
  # for (i_case in 1)
  # {
  plot.res.autocorr.scale.vs.log.body.size[[i_case]] = cor.plot(y = lm(charac_scale[selected_groups,i_case] ~ NormalizedVI_pcoa[[3]][,1])$residuals,
                                                                    x = log(size_relativeAbund[selected_groups]),
                                                                    # y.lab = paste(ifelse(i_case == 1,"Scale of biogeographic organization (km)",
                                                                    #                      "Scale of biogeo. organization(km; DCM)")),
                                                                    excluded.points = which(taxo_groups[selected_groups] %in% "RAD-A"),
                                                                    y.lab = "Scale of surface biogeographic organization\n- residuals after regression on axis 1",
                                                                    x.lab = expression("Body size ("*mu*"m), log scale"),
                                                                    x.log = T,
                                                                    y.log = F,
                                                                    fit = T,
                                                                    fit.display="spearman.p",
                                                                    # fit.display="pearson.spearman",
                                                                    x.cor.pos = 0.8,
                                                                    y.cor.pos = 0.95)
  
  plot.res.basin.I.within.vs.log.body.size[[i_case]] = cor.plot(y = lm(basin_I[selected_groups,i_case]/I_square.observed_w.mean[selected_groups,i_case] ~ NormalizedVI_pcoa[[3]][,1])$residuals,
                                                                    x = log(size_relativeAbund[selected_groups]),
                                                                    excluded.points = which(taxo_groups[selected_groups] %in% "RAD-C"),
                                                                    # y.lab = paste0("Within-basin homogeneity",ifelse(i_case == 1,""," (DCM)")),
                                                                    y.lab = "Surface within-basin autocorrelation\n- residuals after regression on axis 1",
                                                                    x.lab = expression("Body size ("*mu*"m), log scale"),
                                                                    x.log = T,
                                                                    y.log = F,
                                                                    fit = T,
                                                                    fit.display="spearman.p",
                                                                    # fit.display="pearson.spearman",
                                                                    x.cor.pos = 0.8,
                                                                    y.cor.pos = 0.1)
  
  plot.res.lat.I.sym.vs.log.body.size[[i_case]] = cor.plot(y = lm(lat_I[selected_groups,i_case]/I_square.observed_w.mean[selected_groups,i_case] ~ NormalizedVI_pcoa[[3]][,1])$residuals,
                                                               x = log(size_relativeAbund[selected_groups]),
                                                               excluded.points = which(taxo_groups[selected_groups] %in% c("RAD-A","RAD-C")),
                                                               # y.lab = paste0("Latitudinal symmetry",ifelse(i_case == 1,""," (DCM)")),
                                                               y.lab = "Surface latitudinal autocorrelation\n- residuals after regression on axis 1",
                                                               x.lab = expression("Body size ("*mu*"m), log scale"),
                                                               x.log = T,
                                                               y.log = F,
                                                               fit = T,
                                                               fit.display="spearman.p",
                                                               # fit.display="pearson.spearman",
                                                               x.cor.pos = 0.8,
                                                               y.cor.pos = 0.8)
  # }
  # pdf(paste0(figure_folder,"/Fig.SI_residuals.biographic.descriptors.vs.log.body.size_",ifelse(i_case==1,"","DCM_"),"excluded.points.pdf"),
  #     width=7.5*2.2/1.2*2/2,height=12/3*2/1.2*2,onefile=F)
  # g1 = ggplotGrob(plot.res.autocorr.scale.vs.log.body.size[[i_case]])
  # g2 = ggplotGrob(plot.res.basin.I.within.vs.log.body.size[[i_case]])
  # g3 = ggplotGrob(plot.res.lat.I.sym.vs.log.body.size[[i_case]])
  # grid.arrange(grobs = list(g1,g2,g3),
  #              layout_matrix = rbind(c(1, 2),c(3, NA)),
  #              respect=T)
  # dev.off()
  
  ###############
  i_case = 1
  boxplot.res.autocorr.scale = list()
  boxplot.res.basin.I.within = list()
  boxplot.res.lat.I.sym = list()
  
  # g = factor(dominant_function1,c("Phagotrophs","Phototrophs","Metazoans","Parasites"))
  # Meta.vs.Photo = c(0,1,-1,0)
  # Phago.vs.Meta = c(1,0,-1,0)
  # Phago.vs.Photo = c(1,-1,0,0)
  # contrasts(g) = cbind(Meta.vs.Photo,Phago.vs.Meta,Phago.vs.Photo)
  # y = charac_scale[selected_groups & !taxo_groups %in% c("Dinophyceae","Collodaria"),1]
  # model = lm(y ~ g)
  # anova(model)
  # summary(model)
  
  # excluded_groups = "Porifera"
  # detached_groups = c("RAD-C","Ascomycota")
  boxplot.res.autocorr.scale[[i_case]] = box.plot(classes = factor(dominant_function1, c("Phototrophs","Phagotrophs","Metazoans","Parasites")),
                                                  values = lm(charac_scale[selected_groups & !taxo_groups %in% c("Dinophyceae","Collodaria"),i_case] ~ 
                                                                NormalizedVI_pcoa[[3]][!names(NormalizedVI_pcoa[[3]][,1]) %in% c("Dinophyceae","Collodaria"),1])$residuals,
                                                  # detached.points = which(names(dominant_function1) %in% detached_groups),
                                                  # excluded.points = which(names(dominant_function1) %in% excluded_groups),
                                                  excluded.points = which(names(dominant_function1) %in% "RAD-A"),
                                                  y.lab = "Scale of biogeographic organization\n- residuals after regression on axis 1",
                                                  # y.lab.hjust = 0.39,
                                                  x.angle = 14,
                                                  fit = T,
                                                  y.cor.pos = 0.9,
                                                  mar.vect = c(5,5,5,5)) +
    theme(axis.text.x = element_text(vjust = 0.8, hjust = 0.7))
  
  boxplot.res.basin.I.within[[i_case]] = box.plot(classes = factor(dominant_function1, c("Phototrophs","Phagotrophs","Metazoans","Parasites")),
                                                  values = lm(basin_I[selected_groups & !taxo_groups %in% c("Dinophyceae","Collodaria"),i_case]/
                                                                I_square.observed_w.mean[selected_groups & !taxo_groups %in% c("Dinophyceae","Collodaria"),i_case] ~ 
                                                                NormalizedVI_pcoa[[3]][!names(NormalizedVI_pcoa[[3]][,1]) %in% c("Dinophyceae","Collodaria"),1])$residuals,
                                                  # detached.points = which(names(dominant_function1) %in% detached_groups),
                                                  # excluded.points = which(names(dominant_function1) %in% excluded_groups),
                                                  excluded.points = which(names(dominant_function1) %in% "RAD-C"),
                                                  y.lab = "Within-basin autocorrelation\n- residuals after regression on axis 1",
                                                  # y.lab.hjust = 0.39,
                                                  x.angle = 14,
                                                  fit = T,
                                                  y.cor.pos = 0.95,
                                                  mar.vect = c(5,5,5,5)) +
    theme(axis.text.x = element_text(vjust = 0.8, hjust = 0.7))
  
  boxplot.res.lat.I.sym[[i_case]] = box.plot(classes = factor(dominant_function1, c("Phototrophs","Phagotrophs","Metazoans","Parasites")),
                                             values = lm(lat_I[selected_groups & !taxo_groups %in% c("Dinophyceae","Collodaria"),i_case]/
                                                           I_square.observed_w.mean[selected_groups & !taxo_groups %in% c("Dinophyceae","Collodaria"),i_case] ~ 
                                                           NormalizedVI_pcoa[[3]][!names(NormalizedVI_pcoa[[3]][,1]) %in% c("Dinophyceae","Collodaria"),1])$residuals,
                                             # detached.points = which(names(dominant_function1) %in% detached_groups),
                                             # excluded.points = which(names(dominant_function1) %in% excluded_groups),
                                             excluded.points = which(names(dominant_function1) %in% c("RAD-A","RAD-C")),
                                             y.lab = "Latitudinal autocorrelation\n- residuals after regression on axis 1",
                                             # y.lab.hjust = 0.39,
                                             x.angle = 14,
                                             fit = T,
                                             y.cor.pos = 0.9,
                                             mar.vect = c(5,5,5,5)) +
    theme(axis.text.x = element_text(vjust = 0.8, hjust = 0.7))
  
  # pdf(paste0(figure_folder,"/Fig.SI_residuals.biographic.descriptors.vs.ecological.function_",ifelse(i_case==1,"","DCM_"),"excluded.points.pdf"),
  #     width=7.5*2.2/1.2*2/2,height=12/3*2/1.2*2,onefile=F)
  # g1 = ggplotGrob(boxplot.res.autocorr.scale[[i_case]])
  # g2 = ggplotGrob(boxplot.res.basin.I.within[[i_case]])
  # g3 = ggplotGrob(boxplot.res.lat.I.sym[[i_case]])
  # grid.arrange(grobs = list(g1,g2,g3),
  #              layout_matrix = rbind(c(1, 2),c(3, NA)),
  #              respect=T)
  # dev.off()
  
  ##############
  i_case = 1
  
  pdf(paste0(figure_folder,"/Fig.SI_residuals.biographic.descriptors.SUR_vs_log.body.size.ecological.function.pdf"),
      width=7.5*2.2/1.2/2*2,height=12/3*2/1.2*3,onefile=F)
  g1 = ggplotGrob(plot.res.autocorr.scale.vs.log.body.size[[i_case]])
  g2 = ggplotGrob(boxplot.res.autocorr.scale[[i_case]] + theme(axis.title.y=element_blank()))
  g3 = ggplotGrob(plot.res.basin.I.within.vs.log.body.size[[i_case]])
  g4 = ggplotGrob(boxplot.res.basin.I.within[[i_case]] + theme(axis.title.y=element_blank()))
  g5 = ggplotGrob(plot.res.lat.I.sym.vs.log.body.size[[i_case]])
  g6 = ggplotGrob(boxplot.res.lat.I.sym[[i_case]] + theme(axis.title.y=element_blank()))
  grid.arrange(grobs = list(g1,g2,g3,g4,g5,g6),
               layout_matrix = rbind(c(1, 2),c(3, 4),c(5, 6)),
               respect=T)
  dev.off()
  
  #############
  # Without taking the residuals:
  
  i_case = 1
  plot.autocorr.scale.vs.log.body.size = list()
  plot.basin.I.within.vs.log.body.size = list()
  plot.lat.I.sym.vs.log.body.size = list()
  boxplot.autocorr.scale = list()
  boxplot.basin.I.within = list()
  boxplot.lat.I.sym = list()
  # for (i_case in 1)
  # {
  plot.autocorr.scale.vs.log.body.size[[i_case]] = cor.plot(y = charac_scale[selected_groups,i_case],
                                                            x = log(size_relativeAbund[selected_groups]),
                                                            # y.lab = paste(ifelse(i_case == 1,"Scale of biogeographic organization (km)",
                                                            #                      "Scale of biogeo. organization(km; DCM)")),
                                                            # excluded.points = which(taxo_groups[selected_groups] %in% "RAD-A"),
                                                            y.lab = "Scale of biogeographic organization (km)",
                                                            x.lab = expression("Body size ("*mu*"m), log scale"),
                                                            x.log = T,
                                                            y.log = F,
                                                            fit = T,
                                                            fit.display="spearman.p",
                                                            # fit.display="pearson.spearman",
                                                            x.cor.pos = 0.8,
                                                            y.cor.pos = 0.97)
  
  plot.basin.I.within.vs.log.body.size[[i_case]] = cor.plot(y = basin_I[selected_groups,i_case]/I_square.observed_w.mean[selected_groups,i_case],
                                                            x = log(size_relativeAbund[selected_groups]),
                                                            excluded.points = which(taxo_groups[selected_groups] %in% "RAD-C"),
                                                            # y.lab = paste0("Within-basin homogeneity",ifelse(i_case == 1,""," (DCM)")),
                                                            y.lab = "Within-basin homogeneity",
                                                            x.lab = expression("Body size ("*mu*"m), log scale"),
                                                            x.log = T,
                                                            y.log = F,
                                                            fit = T,
                                                            fit.display="spearman.p",
                                                            # fit.display="pearson.spearman",
                                                            x.cor.pos = 0.8,
                                                            y.cor.pos = 0.1)
  
  plot.lat.I.sym.vs.log.body.size[[i_case]] = cor.plot(y = lat_I[selected_groups,i_case]/I_square.observed_w.mean[selected_groups,i_case],
                                                       x = log(size_relativeAbund[selected_groups]),
                                                       # excluded.points = which(taxo_groups[selected_groups] %in% c("RAD-A","RAD-C")),
                                                       # y.lab = paste0("Latitudinal symmetry",ifelse(i_case == 1,""," (DCM)")),
                                                       y.lab = "Latitudinal symmetry",
                                                       x.lab = expression("Body size ("*mu*"m), log scale"),
                                                       x.log = T,
                                                       y.log = F,
                                                       fit = T,
                                                       fit.display="spearman.p",
                                                       # fit.display="pearson.spearman",
                                                       x.cor.pos = 0.8,
                                                       y.cor.pos = 0.9)
  
  boxplot.autocorr.scale[[i_case]] = box.plot(classes = factor(dominant_function1, c("Phototrophs","Phagotrophs","Metazoans","Parasites")),
                                              values = charac_scale[selected_groups & !taxo_groups %in% c("Dinophyceae","Collodaria"),i_case],
                                              # detached.points = which(names(dominant_function1) %in% detached_groups),
                                              # excluded.points = which(names(dominant_function1) %in% excluded_groups),
                                              # excluded.points = which(names(dominant_function1) %in% "RAD-A"),
                                              y.lab = "Scale of biogeographic organization",
                                              # y.lab.hjust = 0.39,
                                              fit = T,
                                              y.cor.pos = 0.05,
                                              mar.vect = c(5,5,5,5))
  
  boxplot.basin.I.within[[i_case]] = box.plot(classes = factor(dominant_function1, c("Phototrophs","Phagotrophs","Metazoans","Parasites")),
                                              values = basin_I[selected_groups & !taxo_groups %in% c("Dinophyceae","Collodaria"),i_case]/
                                                I_square.observed_w.mean[selected_groups & !taxo_groups %in% c("Dinophyceae","Collodaria"),i_case],
                                              # detached.points = which(names(dominant_function1) %in% detached_groups),
                                              # excluded.points = which(names(dominant_function1) %in% excluded_groups),
                                              excluded.points = which(names(dominant_function1) %in% "RAD-C"),
                                              y.lab = "Within-basin homogeneity",
                                              # y.lab.hjust = 0.39,
                                              fit = T,
                                              y.cor.pos = 0.1,
                                              mar.vect = c(5,5,5,5))
  
  boxplot.lat.I.sym[[i_case]] = box.plot(classes = factor(dominant_function1, c("Phototrophs","Phagotrophs","Metazoans","Parasites")),
                                         values = lat_I[selected_groups & !taxo_groups %in% c("Dinophyceae","Collodaria"),i_case]/
                                           I_square.observed_w.mean[selected_groups & !taxo_groups %in% c("Dinophyceae","Collodaria"),i_case],
                                         # detached.points = which(names(dominant_function1) %in% detached_groups),
                                         # excluded.points = which(names(dominant_function1) %in% excluded_groups),
                                         # excluded.points = which(names(dominant_function1) %in% c("RAD-A","RAD-C")),
                                         y.lab = "Latitudinal symmetry",
                                         # y.lab.hjust = 0.39,
                                         fit = T,
                                         y.cor.pos = 0.9,
                                         mar.vect = c(5,5,5,5)) 
  pdf(paste0(figure_folder,"/Fig.SI_biographic.descriptors.SUR_vs_log.body.size.ecological.function.pdf"),
      width=7.5*2.2/1.2/2*2,height=12/3*2/1.2*3,onefile=F)
  g1 = ggplotGrob(plot.autocorr.scale.vs.log.body.size[[i_case]])
  g2 = ggplotGrob(boxplot.autocorr.scale[[i_case]] + theme(axis.title.y=element_blank()))
  g3 = ggplotGrob(plot.basin.I.within.vs.log.body.size[[i_case]])
  g4 = ggplotGrob(boxplot.basin.I.within[[i_case]] + theme(axis.title.y=element_blank()))
  g5 = ggplotGrob(plot.lat.I.sym.vs.log.body.size[[i_case]])
  g6 = ggplotGrob(boxplot.lat.I.sym[[i_case]] + theme(axis.title.y=element_blank()))
  grid.arrange(grobs = list(g1,g2,g3,g4,g5,g6),
               layout_matrix = rbind(c(1, 2),c(3, 4),c(5, 6)),
               respect=T)
  dev.off()
  # }
}

# Fig. S variation partitioning
{
  g1 = g2 = list()
  for (i_case in 1)
  {
    # tot_sorting_indSelec = sort.int(colSums(varpart.env.spatial[,selected_groups]), index.return = T, decreasing = T, na.last = T, method = "radix")$ix
    tot_sorting = sort.int(colSums(varpart.env.spatial[[i_case]])[selected_groups & diversity>div_threshold], index.return = T, decreasing = T, na.last = T, method = "radix")$ix
    clades_annotations = c("Chordata","Diplonemida","Arthropoda","Bacillariophyta","Dinophyceae","Collodaria","Haptophyta","MAST-3,_12")
    clades_annotations_labels0 = rep("",length(which(selected_groups & diversity>div_threshold)))
    clades_annotations_labels = clades_annotations_labels0
    clades_annotations_labels[which(taxo_groups[selected_groups & diversity>div_threshold][tot_sorting] %in% clades_annotations)] = c("Chordata","Diplonemida","Arthropoda","Bacillariophyta","Dinophyceae","Collodaria","Haptophyta","MAST 3,12")
    
    # "Pure biotic","Mixed biotic-abiotic","Mixed biotic-currents","Mixed biotic-abiotic-currents","Pure abiotic","Mixed abiotic-currents","Pure currents"
    barplot_data.frame = rbind(varpart.env.spatial[[i_case]][3,],varpart.env.spatial[[i_case]][2,],varpart.env.spatial[[i_case]][1,])
    dimnames(barplot_data.frame) = list(c("Purely by Moran maps","Jointly by Moran maps and local environment","Purely by local environment"), taxo_groups_unmodified)
    barplot_data.frame = barplot_data.frame[,selected_groups & diversity > div_threshold][,tot_sorting]
    #
    barplot.env.connectivity = bar.plot(barplot_data.frame, 
                                       y.lab = paste("\nFraction of variance explained",c("","(DCM)")[i_case]),
                                       colors = rev(terrain.colors(3)),
                                       bar.width = 0.7,
                                       legend.position=c(0.7,0.7),
                                       mar.vect=c(5,5,-5,5),
                                       x.text = clades_annotations_labels,
                                       size.x.text = 16, angle.x.text = 60, vjust.x.text=1.2, hjust.x.text = 1.2)
    # Default parameter values: angle.x.text = 90, vjust.x.text=0.5, hjust.x.text = 0
    g1[[i_case]] = ggplotGrob(barplot.env.connectivity) #+ 
                                # theme_void() +
                                # geom_text(aes(label = c(rep(clades_annotations_labels0,2),clades_annotations_labels),
                                #               x = rep(1:length(which(selected_groups & diversity>div_threshold)),3),
                                #               y = -0.01),
                                #           angle = 60, size = 6, hjust = "right"))
                                # scale_x_discrete(limits=clades_annotations_labels,breaks=which(taxo_groups[selected_groups & diversity>div_threshold][tot_sorting] == clades_annotations))) 
    
    # tot_sorting = sort.int(colSums(varpart.biotic.abiotic[[i_case]])[selected_groups & diversity>div_threshold], index.return = T, decreasing = T, na.last = T, method = "radix")$ix
    # "Pure biotic","Mixed biotic-abiotic","Mixed biotic-currents","Mixed biotic-abiotic-currents","Pure abiotic","Mixed abiotic-currents","Pure currents"
    barplot_data.frame = rbind(varpart.biotic.abiotic[[i_case]][1,],varpart.biotic.abiotic[[i_case]][2,],varpart.biotic.abiotic[[i_case]][3,])
    dimnames(barplot_data.frame) = list(c("Purely by local abiotic conditions","Jointly by local biotic and abiotic conditions","Purely by local biotic conditions"), taxo_groups_unmodified)
    barplot_data.frame = barplot_data.frame[,selected_groups & diversity > div_threshold][,tot_sorting]
    #
    barplot.biotic.abiotic = bar.plot(barplot_data.frame, 
                                     y.lab = paste("\nFraction of variance explained",c("","(DCM)")[i_case]),
                                     colors = rev(c("deepskyblue4","cadetblue","chartreuse3")),
                                     bar.width = 0.7,
                                     legend.position=c(0.7,0.7),
                                     mar.vect=c(1,5,-5,5),
                                     x.text = clades_annotations_labels,
                                     size.x.text = 16, angle.x.text = 60, vjust.x.text=1.2, hjust.x.text = 1.2)
    g2[[i_case]] = ggplotGrob(barplot.biotic.abiotic) #+
                                # theme_void() +
                                # geom_text(aes(label = c(rep(clades_annotations_labels0,2),clades_annotations_labels),
                                #               x = rep(1:length(which(selected_groups & diversity>div_threshold)),3),
                                #               y = -0.01),
                                #           angle = 60, size = 6, hjust = "right")) 
  }
  
  # pdf(paste0(figure_folder,"/FigS_varpart_dis.MEM_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
  pdf(paste0(figure_folder,"/FigS_varpart_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1_BH.correction.pdf"),
      width=7.5*2.2/1.2,height=12/3*4/1.2)
  # for (i_case in 1:2)
  # {
    grid.arrange(grobs = list(g1[[i_case]], g2[[i_case]]), 
                 layout = rbind(c(1, 1),c(2, 2)))
  # }
  dev.off()
  
}

# Fig. S. sup. correlations diversity & body size:
{
  size = 1.2
  
  x = size_relativeAbund[selected_groups]
  y = NormalizedVI_pcoa[[3]][,1]
  plot.norm.VI.PCoA.axis1.vs.body.size = cor.plot(x = x,
                                                  y = y,
                                                  x.lab = expression("Body size ("*mu*"m), log scale"),
                                                  y.lab = "PCoA axis 1",
                                                  y.lab.hjust = 0.5,
                                                  size.factor = size,
                                                  x.log = T,
                                                  y.log = F,
                                                  fit = T,
                                                  fit.display = "spearman.p", 
                                                  x.cor.pos = 0.8,
                                                  y.cor.pos = 0.1,
                                                  mar.vect = c(5,5,1,5))
  
  boxplot.norm.VI.PCoA.axis1.functions = box.plot(classes = factor(dominant_function1, c("Phototrophs","Phagotrophs","Metazoans","Parasites")),
                                                            values = NormalizedVI_pcoa[[3]][!names(NormalizedVI_pcoa[[3]][,1]) %in% c("Dinophyceae","Collodaria"),1],
                                                            y.lab = "PCoA axis 1",
                                                            y.lab.hjust = 0.39,
                                                            size.factor = size,
                                                            fit = T,
                                                            fit.size.factor = 0.85,
                                                            mar.vect = c(5,5,1,5))
  
  x = as.vector(diversity)[selected_groups]
  y = NormalizedVI_pcoa[[3]][,2]
  # x.2000 = x[x<2000]
  # y.2000 = y[x<2000]
  # x.cor.pos = 0.8
  # y.cor.pos = 0.1
  # excluded_groups = "Porifera"
  plot.norm.VI.PCoA.axis2.vs.diversity = cor.plot(x = x,
                                                  y = y,
                                                  # excluded.points = which(names(NormalizedVI_pcoa[[3]][,2]) == excluded_groups),
                                                  x.lab ="Diversity (#OTUs), log scale",
                                                  y.lab = "PCoA axis 2",
                                                  y.lab.hjust = 0.5,
                                                  size.factor = size,
                                                  x.log = T,
                                                  y.log = F,
                                                  fit = T,
                                                  fit.display = "spearman.p", 
                                                  x.cor.pos = 0.8,
                                                  y.cor.pos = 0.1,
                                                  mar.vect = c(5,12,1,5))
    # geom_smooth(data = data.frame(x = x.2000, y = y.2000),aes(x,y),method='lm',col="black",se=T) +
    # # geom_smooth(aes(x[x<2000], y[x<2000]),method='lm',col="black",se=T, inherit.aes = T) +
    # geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    # annotate(geom="text", 
    #          x=(max(x)/min(x))^x.cor.pos*min(x),
    #          y=y.cor.pos*(max(y)-min(y))+min(y),
    #          label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
    #          size=8)
  
 
  
  pdf(paste0(figure_folder,"/FigS_PCoA.axes.div.size.ecology_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1_excluded.Porifera.pdf"),
      # width=7.5*2.2/1.2/2*3,height=12/3*4/1.2/2,onefile=F)
      width=7.5*2.2/1.2/2*3,height=12/3*4/1.2/2,onefile=F)
  g1 = ggplotGrob(plot.norm.VI.PCoA.axis1.vs.body.size)
  g2 = ggplotGrob(boxplot.norm.VI.PCoA.axis1.functions)
  g3 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.diversity)
  grid.arrange(grobs = list(g1, g2, g3), nrow = 1, ncol = 3)
  dev.off()
  
  # pdf(paste0(figure_folder,"/FigSup_correlations.no.axis_",div_threshold,"plusOTUs_logFractionRatio_selected100+1.pdf"),
  #     width=7.5*2.2/1.2*3/2,height=12/3*4/1.2*3/2)
  # for (i_case in 1:2)
  # {
  #   x = as.vector(diversity)[selected_groups & diversity>div_threshold][-69]
  #   y = I_square.observed_w.mean[selected_groups & diversity>div_threshold,i_case][-69]
  #   # lin.regr = lm(y ~ x)
  #   # b = summary(lin.regr)$coefficients[1,1]
  #   # a = summary(lin.regr)$coefficients[2,1]
  #   plot.spatial.autocorr.vs.diversity = ggplot(data = data.frame(x = x, y = y)) +
  #     geom_point(aes(x,y)) +
  #     # geom_point(aes(x,y2), col = "red") +
  #     scale_x_log10() +
  #     # scale_y_log10() +
  #     # geom_abline(intercept = b, slope = a) +
  #     theme_bw() +
  #     ggtitle(LETTERS[1]) +
  #     # geom_hline(yintercept = 1, linetype="dashed") +
  #     # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
  #     # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
  #     theme(axis.text = element_text(size=22),
  #           axis.title=element_text(size=22),
  #           # axis.title.x=element_text(vjust = 45),
  #           plot.title=element_text(hjust=0, size=24),
  #           plot.margin=unit(c(5,7,0.5,0.5),"mm")) +
  #     #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
  #     labs(x="Clade diversity", y=paste("Short-distance spatial autocorr.",c("(Surface)","(DCM)")[i_case])) +
  #     geom_smooth(data = data.frame(x = x[x<2000], y = y[x<2000]), aes(x,y),method='lm',col="black",se=T) +
  #     geom_segment(aes(x = 2000, xend = max(x), y = if (i_case == 1) 0.705 else 0.63, yend = if (i_case == 1) 0.705 else 0.63), linetype="dashed")
  #   # geom_smooth(aes(x,y2),method='lm',col="red")
  #   g1 = ggplotGrob(plot.spatial.autocorr.vs.diversity)
  #   
  #   x = as.vector(diversity)[selected_groups & diversity>div_threshold][-69]
  #   y = 1 - SUR.DCM_Normalized.VI[selected_groups & diversity>div_threshold][-69]
  #   # lin.regr = lm(y ~ log10(x))
  #   # b = summary(lin.regr)$coefficients[1,1]
  #   # a = summary(lin.regr)$coefficients[2,1]
  #   plot.SUR.DCM.sim.vs.diversity = ggplot(data = data.frame(x = x, y = y)) +
  #     # geom_point(aes(x,y1)) +
  #     geom_point(aes(x,y)) +
  #     scale_x_log10() +
  #     # scale_y_log10() +
  #     # geom_abline(intercept = b, slope = a) +
  #     theme_bw() +
  #     ggtitle(LETTERS[2]) +
  #     # geom_hline(yintercept = 1, linetype="dashed") +
  #     # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
  #     # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
  #     theme(axis.text = element_text(size=22),
  #           axis.title=element_text(size=22),
  #           # axis.title.x=element_text(vjust = 45),
  #           plot.title=element_text(hjust=0, size=24),
  #           plot.margin=unit(c(5,7,0.5,0.5),"mm")) +
  #     #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
  #     labs(x="Clade diversity", y="Similarity between surface and DCM") +
  #     geom_smooth(data = data.frame(x = x[x<2000], y = y[x<2000]), aes(x,y),method='lm',col="black",se=T) +
  #     geom_segment(aes(x = 2000, xend = max(x), y = 0.3175, yend = 0.3175), linetype="dashed")
  #   # geom_smooth(aes(x,y),method='lm')
  #   g2 = ggplotGrob(plot.SUR.DCM.sim.vs.diversity)
  #   
  #   x = size_relativeAbund[selected_groups & diversity>div_threshold][-69]
  #   y = charac_scale[selected_groups & diversity>div_threshold,i_case][-69]
  #   # lin.regr = lm(y ~ log10(x))
  #   # b = summary(lin.regr)$coefficients[1,1]
  #   # a = summary(lin.regr)$coefficients[2,1]
  #   plot.spatial.scale.vs.body.size = ggplot(data = data.frame(x = x, y = y)) +
  #     # geom_point(aes(x,y1)) +
  #     geom_point(aes(x,y)) +
  #     scale_x_log10() +
  #     # scale_y_log10() +
  #     # geom_abline(intercept = b, slope = a) +
  #     theme_bw() +
  #     ggtitle(LETTERS[if (i_case == 1) 3 else 2]) +
  #     # geom_hline(yintercept = 1, linetype="dashed") +
  #     # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
  #     # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
  #     theme(axis.text = element_text(size=22),
  #           axis.title=element_text(size=22),
  #           # axis.title.x=element_text(vjust = 45),
  #           plot.title=element_text(hjust=0, size=24),
  #           plot.margin=unit(c(5,7,0.5,0.5),"mm")) +
  #     #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
  #     labs(x=expression("Mean size ("*mu*"m)"), y=paste("Spatial autocorr. scale",c("(km, Surface)","(km, DCM)")[i_case])) +
  #     geom_smooth(aes(x,y),method='lm',col="black") 
  #   # geom_smooth(aes(x,y),method='lm')
  #   g3 = ggplotGrob(plot.spatial.scale.vs.body.size)
  #   
  #   boxplot.spatial.scale.functions = ggplot(data = data.frame(x = dominant_function0[!is.na(dominant_function0) & !dominant_function0 %in% point_groups],
  #                                                              y = charac_scale[selected_groups & diversity>div_threshold,i_case][-69][!is.na(dominant_function0) & !dominant_function0 %in% point_groups])) +
  #     scale_x_discrete(limits=levels(dominant_function0)) +
  #     geom_boxplot(aes(x,y)) +
  #     geom_point(data = data.frame(x = factor(point_groups),
  #                                  y = charac_scale[selected_groups & diversity>div_threshold,i_case][-69][dominant_function0 %in% point_groups]),
  #                aes(x,y)) +
  #     theme_bw() +
  #     ggtitle(LETTERS[if (i_case == 1) 4 else 3]) +
  #     # scale_y_log10() +
  #     # geom_hline(yintercept = 1, linetype = "dashed") +
  #     theme(axis.title=element_text(size=22),
  #           # axis.title.y = element_text(hjust = 0.6),
  #           axis.text=element_text(size=22),
  #           axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  #           plot.title = element_text(hjust=0, size=24),
  #           plot.margin=unit(c(5,7,0.5,0.5),"mm")) +
  #     labs(x="", y=paste("Spatial autocorr. scale",c("(km, Surface)","(km, DCM)")[i_case]))
  #   g4 = ggplotGrob(boxplot.spatial.scale.functions)
  #   
  #   plot.tot.var.vs.autocorr = ggplot(data=data.frame(x=I_square.observed_w.mean[selected_groups & diversity>div_threshold,i_case][-69],
  #                                                     y=colSums(varpart.env.spatial[[i_case]])[selected_groups & diversity>div_threshold][-69])) +
  #     geom_point(aes(x,y)) +
  #     # scale_y_log10() +
  #     theme_bw() +
  #     ggtitle(LETTERS[if (i_case == 1) 5 else 4]) +
  #     # geom_hline(yintercept = 1, linetype="dashed") +
  #     # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
  #     # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
  #     theme(axis.text = element_text(size=22),
  #           axis.title=element_text(size=22),
  #           # axis.title.x=element_text(vjust = 45),
  #           plot.title=element_text(hjust=0, size=24),
  #           plot.margin=unit(c(5,7,0.5,0.5),"mm")) +
  #     labs(x=paste("Short-distance spatial autocorr.",c("(Surface)","(DCM)")[i_case]),
  #          y=paste("Total variance explained\n by envir. and connectivity",c("(Surface)","(DCM)")[i_case])) +
  #     # labs(x=paste("Number of",if (i_case == 1) "surface" else "DCM","community types"),y="") +
  #     geom_smooth(aes(x,y),method='lm',col="black") 
  #   g5 = ggplotGrob(plot.tot.var.vs.autocorr)
  #   
  #   plot.tot.var.vs.SUR.DCM.sim = ggplot(data=data.frame(x=1 - SUR.DCM_Normalized.VI[selected_groups & diversity>div_threshold][-69],
  #                                                        y=colSums(varpart.env.spatial[[i_case]])[selected_groups & diversity>div_threshold][-69])) +
  #     geom_point(aes(x,y)) +
  #     # scale_y_log10() +
  #     theme_bw() +
  #     ggtitle(LETTERS[if (i_case == 1) 6 else 5]) +
  #     # geom_hline(yintercept = 1, linetype="dashed") +
  #     # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
  #     # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
  #     theme(axis.text = element_text(size=22),
  #           axis.title=element_text(size=22),
  #           # axis.title.x=element_text(vjust = 45),
  #           plot.title=element_text(hjust=0, size=24),
  #           plot.margin=unit(c(5,7,0.5,0.5),"mm")) +
  #     labs(x="Similarity between surface and DCM",
  #          y=paste("Total variance explained\n by envir. and connectivity",c("(Surface)","(DCM)")[i_case])) +
  #     # labs(x=paste("Number of",if (i_case == 1) "surface" else "DCM","community types"),y="") +
  #     geom_smooth(aes(x,y),method='lm',col="black") 
  #   g6 = ggplotGrob(plot.tot.var.vs.SUR.DCM.sim)
  #   
  #   plot.connectivity.env.ratio.vs.spatial.scale = ggplot(data=data.frame(x=charac_scale[selected_groups & diversity>div_threshold, i_case][-69],
  #                                                                         y=(varpart.env.spatial[[i_case]][3,]/varpart.env.spatial[[i_case]][1,])[selected_groups & diversity>div_threshold][-69])) +
  #     geom_point(aes(x,y)) +
  #     # scale_x_log10() +
  #     scale_y_log10() +
  #     theme_bw() +
  #     ggtitle(LETTERS[if (i_case == 1) 7 else 6]) +
  #     geom_hline(yintercept = 1, linetype="dashed") +
  #     # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
  #     # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
  #     theme(axis.text = element_text(size=22),
  #           axis.title=element_text(size=22),
  #           # axis.title.x=element_text(vjust = 45),
  #           plot.title=element_text(hjust=0, size=24),
  #           plot.margin=unit(c(5,7,0.5,0.5),"mm")) +
  #     #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
  #     labs(x=paste("Spatial autocorr. scale",c("(km, Surface)","(km, DCM)")[i_case]),
  #          y=paste("Ratio of purely connectivity-explained\n over envir.-explained variance",c("(Surface)","(DCM)")[i_case])) +
  #     geom_smooth(aes(x,y),method='lm',col="black")
  #   g7 = ggplotGrob(plot.connectivity.env.ratio.vs.spatial.scale)
  #   
  #   plot.connectivity.env.ratio.vs.size = ggplot(data=data.frame(x=size_relativeAbund[selected_groups & diversity>div_threshold][-69],
  #                                                                y=(varpart.env.spatial[[i_case]][3,]/varpart.env.spatial[[i_case]][1,])[selected_groups & diversity>div_threshold][-69])) +
  #     geom_point(aes(x,y)) +
  #     scale_x_log10() +
  #     scale_y_log10() +
  #     theme_bw() +
  #     ggtitle(LETTERS[if (i_case == 1) 8 else 7]) +
  #     geom_hline(yintercept = 1, linetype="dashed") +
  #     # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
  #     # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
  #     theme(axis.text = element_text(size=22),
  #           axis.title=element_text(size=22),
  #           # axis.title.x=element_text(vjust = 45),
  #           plot.title=element_text(hjust=0, size=24),
  #           plot.margin=unit(c(5,7,0.5,0.5),"mm")) +
  #     #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
  #     labs(x=expression("Mean size ("*mu*"m)"), 
  #          y=paste("Ratio of purely connectivity-explained\n over envir.-explained variance",c("(Surface)","(DCM)")[i_case])) +
  #     geom_smooth(aes(x,y),method='lm',col="black")
  #   g8 = ggplotGrob(plot.connectivity.env.ratio.vs.size)
  #   
  #   boxplot.connectivity.env.ratio.functions = ggplot(data = data.frame(x = dominant_function0[!is.na(dominant_function0) & !dominant_function0 %in% point_groups],
  #                                                                       y = (varpart.env.spatial[[i_case]][3,]/varpart.env.spatial[[i_case]][1,])[selected_groups & diversity>div_threshold][-69]
  #                                                                       [!is.na(dominant_function0) & !dominant_function0 %in% point_groups])) +
  #     scale_x_discrete(limits=levels(dominant_function0)) +
  #     geom_boxplot(aes(x,y)) +
  #     scale_y_log10() +
  #     geom_point(data = data.frame(x = factor(point_groups),
  #                                  y = (varpart.env.spatial[[i_case]][3,]/varpart.env.spatial[[i_case]][1,])[selected_groups & diversity>div_threshold][-69]
  #                                  [dominant_function0 %in% point_groups]),
  #                aes(x,y)) +
  #     theme_bw() +
  #     ggtitle(LETTERS[if (i_case == 1) 9 else 8]) +
  #     # scale_y_log10() +
  #     geom_hline(yintercept = 1, linetype = "dashed") +
  #     theme(axis.title=element_text(size=22),
  #           axis.text=element_text(size=22),
  #           plot.title=element_text(hjust=0, size=24),
  #           axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  #           plot.margin=unit(c(5,7,-2,0.5),"mm")) +
  #     labs(x="", y=paste("Ratio of purely connectivity-explained\n over envir.-explained variance",c("(Surface)","(DCM)")[i_case]))
  #   g9 = ggplotGrob(boxplot.connectivity.env.ratio.functions)
  #   
  #   if (i_case == 1)
  #   {
  #     grid.arrange(grobs = list(g1,g2,g3,g4,g5,g6,g7,g8,g9),
  #                  # widths = rep(1,7),
  #                  layout_matrix = rbind(c(1, 2, 3),c(4, 5, 6),c(7, 8, 9)))
  #   } else
  #   {
  #     grid.arrange(grobs = list(g1,g3,g4,g5,g6,g7,g8,g9),
  #                  # widths = rep(1,7),
  #                  layout_matrix = rbind(c(1, 2, 3),c(4, 5, 6),c(7, 8, NA)))
  #   }
  # }
  # dev.off()
}

# Fig. S. sup. correlations diversity & body size - residuals on log diversity and log body size:
{
  size = 1.2
  
  x = size_relativeAbund[selected_groups]
  y = NormalizedVI_pcoa[[3]][,1]
  plot.norm.VI.PCoA.axis1.vs.body.size = cor.plot(x = x,
                                                  y = lm(y ~ log(as.vector(diversity[selected_groups])))$residuals,
                                                  x.lab = expression("Body size ("*mu*"m), log scale"),
                                                  y.lab = "Biogeographic axis 1\n- residuals after regr. on log diversity",
                                                  y.lab.hjust = 0.5,
                                                  size.factor = size,
                                                  x.log = T,
                                                  y.log = F,
                                                  fit = T,
                                                  fit.display = "spearman.p", 
                                                  x.cor.pos = 0.8,
                                                  y.cor.pos = 0.9,
                                                  mar.vect = c(5,5,1,5))
  
  values = NormalizedVI_pcoa[[3]][!names(NormalizedVI_pcoa[[3]][,1]) %in% c("Dinophyceae","Collodaria"),1]
  values.res = lm(values ~ log(as.vector(diversity[selected_groups & !taxo_groups %in% c("Dinophyceae","Collodaria")])))$residuals
  boxplot.norm.VI.PCoA.axis1.functions = box.plot(classes = factor(dominant_function1, c("Phototrophs","Phagotrophs","Metazoans","Parasites")),
                                                  values = values.res,
                                                  y.lab = "Biogeographic axis 1\n- residuals after regr. on log diversity",
                                                  # y.lab.hjust = 0.39,
                                                  size.factor = size,
                                                  fit = T,
                                                  fit.size.factor = 0.85,
                                                  y.cor.pos=0.2,
                                                  mar.vect = c(5,5,1,5)) + theme(axis.title.y = element_blank())
  
  x = as.vector(diversity)[selected_groups]
  y = NormalizedVI_pcoa[[3]][,2]
  # x.2000 = x[x<2000]
  # y.2000 = y[x<2000]
  # x.cor.pos = 0.8
  # y.cor.pos = 0.1
  # excluded_groups = "Collodaria"
  plot.norm.VI.PCoA.axis2.vs.diversity = cor.plot(x = x,
                                                  y = lm(y ~ log(size_relativeAbund[selected_groups]))$residuals,
                                                  # excluded.points = which(names(NormalizedVI_pcoa[[3]][,2]) == excluded_groups),
                                                  x.lab ="Diversity (#OTUs), log scale",
                                                  y.lab = "Biogeographic axis 2\n- residuals after regr. on log body size  ",
                                                  y.lab.hjust = 0.5,
                                                  size.factor = size,
                                                  x.log = T,
                                                  y.log = F,
                                                  fit = T,
                                                  fit.display = "spearman.p", 
                                                  x.cor.pos = 0.8,
                                                  y.cor.pos = 0.35,
                                                  mar.vect = c(5,12,1,5))
  # geom_smooth(data = data.frame(x = x.2000, y = y.2000),aes(x,y),method='lm',col="black",se=T) +
  # # geom_smooth(aes(x[x<2000], y[x<2000]),method='lm',col="black",se=T, inherit.aes = T) +
  # geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
  # annotate(geom="text", 
  #          x=(max(x)/min(x))^x.cor.pos*min(x),
  #          y=y.cor.pos*(max(y)-min(y))+min(y),
  #          label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
  #          size=8)
  
  
  
  pdf(paste0(figure_folder,"/FigS_PCoA.axes.div.size.ecology_residuals.axis1.axis2_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1_excluded.Porifera.pdf"),
      # width=7.5*2.2/1.2/2*3,height=12/3*4/1.2/2,onefile=F)
      width=7.5*2.2/1.2/2*3,height=12/3*4/1.2/2,onefile=F)
  g1 = ggplotGrob(plot.norm.VI.PCoA.axis1.vs.body.size)
  g2 = ggplotGrob(boxplot.norm.VI.PCoA.axis1.functions)
  g3 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.diversity)
  grid.arrange(grobs = list(g1, g2, g3), nrow = 1, ncol = 3)
  dev.off()
}

# Fig. S. sup. correlations diversity & body size -residuals and non-residuals
{
  size = 1
  
  x = size_relativeAbund[selected_groups]
  y = NormalizedVI_pcoa[[3]][,1]
  plot.norm.VI.PCoA.axis1.vs.body.size = cor.plot(x = x,
                                                  y = y,
                                                  x.lab = expression("Body size ("*mu*"m), log scale"),
                                                  y.lab = "Biogeographic axis 1",
                                                  y.lab.hjust = 0.5,
                                                  size.factor = size,
                                                  x.log = T,
                                                  y.log = F,
                                                  fit = T,
                                                  fit.display = "spearman.p", 
                                                  x.cor.pos = 0.8,
                                                  y.cor.pos = 0.9,
                                                  mar.vect = c(5,5,1,5))
  
  x = size_relativeAbund[selected_groups]
  y = NormalizedVI_pcoa[[3]][,1]
  plot.norm.VI.PCoA.axis1.res.vs.body.size = cor.plot(x = x,
                                                  y = lm(y ~ log(as.vector(diversity[selected_groups])))$residuals,
                                                  x.lab = expression("Body size ("*mu*"m), log scale"),
                                                  y.lab = "Biogeographic axis 1\n- residuals after regres. on log diversity",
                                                  y.lab.hjust = 0.5,
                                                  size.factor = size,
                                                  x.log = T,
                                                  y.log = F,
                                                  fit = T,
                                                  fit.display = "spearman.p", 
                                                  x.cor.pos = 0.8,
                                                  y.cor.pos = 0.9,
                                                  mar.vect = c(5,5,1,5))
  
  boxplot.norm.VI.PCoA.axis1.functions = box.plot(classes = factor(dominant_function1, c("Phototrophs","Phagotrophs","Metazoans","Parasites")),
                                                  values = NormalizedVI_pcoa[[3]][!names(NormalizedVI_pcoa[[3]][,1]) %in% c("Dinophyceae","Collodaria"),1],
                                                  y.lab = "Biogeographic axis 1",
                                                  # y.lab.hjust = 0.39,
                                                  size.factor = size,
                                                  x.angle = 14,
                                                  fit = T,
                                                  # fit.size.factor = 0.85,
                                                  y.cor.pos=0.9,
                                                  mar.vect = c(5,5,1,5)) +
    theme(axis.text.x = element_text(vjust = 0.8, hjust = 0.7))
  
  values = NormalizedVI_pcoa[[3]][!names(NormalizedVI_pcoa[[3]][,1]) %in% c("Dinophyceae","Collodaria"),1]
  values.res = lm(values ~ log(as.vector(diversity[selected_groups & !taxo_groups %in% c("Dinophyceae","Collodaria")])))$residuals
  boxplot.norm.VI.PCoA.axis1.res.functions = box.plot(classes = factor(dominant_function1, c("Phototrophs","Phagotrophs","Metazoans","Parasites")),
                                                      values = values.res,
                                                      y.lab = "Biogeographic axis 1\n- residuals after regres. on log diversity",
                                                      # y.lab.hjust = 0.39,
                                                      size.factor = size,
                                                      x.angle = 14,
                                                      fit = T,
                                                      # fit.size.factor = 0.85,
                                                      y.cor.pos=0.2,
                                                      mar.vect = c(5,5,1,5)) +
    theme(axis.text.x = element_text(vjust = 0.8, hjust = 0.7))
  
  x = as.vector(diversity)[selected_groups]
  y = NormalizedVI_pcoa[[3]][,2]
  # x.2000 = x[x<2000]
  # y.2000 = y[x<2000]
  # x.cor.pos = 0.8
  # y.cor.pos = 0.1
  # excluded_groups = "Porifera"
  plot.norm.VI.PCoA.axis2.vs.diversity = cor.plot(x = x,
                                                  y = y,
                                                  # excluded.points = which(names(NormalizedVI_pcoa[[3]][,2]) == excluded_groups),
                                                  x.lab ="Diversity (#OTUs), log scale",
                                                  y.lab = "Biogeographic axis 2",
                                                  y.lab.hjust = 0.5,
                                                  size.factor = size,
                                                  x.log = T,
                                                  y.log = F,
                                                  fit = T,
                                                  fit.display = "spearman.p", 
                                                  x.cor.pos = 0.8,
                                                  y.cor.pos = 0.1,
                                                  mar.vect = c(5,12,1,5))
  
  x = as.vector(diversity)[selected_groups]
  y = NormalizedVI_pcoa[[3]][,2]
  # x.2000 = x[x<2000]
  # y.2000 = y[x<2000]
  # x.cor.pos = 0.8
  # y.cor.pos = 0.1
  # excluded_groups = "Collodaria"
  plot.norm.VI.PCoA.axis2.res.vs.diversity = cor.plot(x = x,
                                                  y = lm(y ~ log(size_relativeAbund[selected_groups]))$residuals,
                                                  # excluded.points = which(names(NormalizedVI_pcoa[[3]][,2]) == excluded_groups),
                                                  x.lab ="Diversity (#OTUs), log scale",
                                                  y.lab = "Biogeographic axis 2\n- residuals after regression on log body size  ",
                                                  y.lab.hjust = 0.5,
                                                  size.factor = size,
                                                  x.log = T,
                                                  y.log = F,
                                                  fit = T,
                                                  fit.display = "spearman.p", 
                                                  x.cor.pos = 0.8,
                                                  y.cor.pos = 0.35,
                                                  mar.vect = c(5,12,1,5))
  
  pdf(paste0(figure_folder,"/FigS_PCoA.axes.div.size.ecology_res.and.non-res_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1_excluded.Porifera.pdf"),
      # width=7.5*2.2/1.2/2*3,height=12/3*4/1.2/2,onefile=F)
      width=7.5*2.2/1.2/2*2,height=12/3*4/1.2/2*3,onefile=F)
  g1 = ggplotGrob(plot.norm.VI.PCoA.axis1.vs.body.size)
  g2 = ggplotGrob(plot.norm.VI.PCoA.axis1.res.vs.body.size)
  g3 = ggplotGrob(boxplot.norm.VI.PCoA.axis1.functions)
  g4 = ggplotGrob(boxplot.norm.VI.PCoA.axis1.res.functions)
  g5 = ggplotGrob(plot.norm.VI.PCoA.axis2.vs.diversity)
  g6 = ggplotGrob(plot.norm.VI.PCoA.axis2.res.vs.diversity)
  grid.arrange(grobs = list(g1, g2, g3, g4, g5, g6), nrow = 3, ncol = 2)
  dev.off()
}

# Fig. S. body size functional plot:
{
  values = size_relativeAbund[selected_groups & diversity>div_threshold & !taxo_groups %in% c("Dinophyceae","Collodaria")]
  body.size.vs.function = box.plot(classes = factor(dominant_function1, c("Phototrophs","Phagotrophs","Metazoans","Parasites")),
                                                  values = values,
                                                  y.lab = expression("Body size ("*mu*"m), log scale"),
                                                  size.factor = 0.8,
                                                  y.log = T,
                                                  x.angle = 14,
                                                  y.lab.hjust = 0.5,
                                                  mar.vect = c(5,5,1,5)) +
    theme(axis.text.x = element_text(vjust = 0.8, hjust = 0.7))
  
  pdf(paste0(figure_folder,"/FigS_body.size.diff.across.ecological.categories_",div_threshold,"plusOTUs_PCoA.axes.x_selected100+1_apusozoans.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2)
  print(body.size.vs.function)
  dev.off()
  
  # Older version:
  #################
  
  pdf(paste0(figure_folder,"/FigS8_body.size.functional.boxplot_",div_threshold,"plusOTUs_selected100+1.pdf"),
      width=7,height=7)
  boxplot.body.size.functions = ggplot(data = data.frame(x = dominant_function0[!is.na(dominant_function0) & !dominant_function0 %in% point_groups],
                                                         y = size_relativeAbund[selected_groups & diversity>div_threshold][-69]
                                                         [!is.na(dominant_function0) & !dominant_function0 %in% point_groups])) +
    scale_x_discrete(limits=levels(dominant_function0)) +
    geom_boxplot(aes(x,y)) +
    scale_y_log10() +
    geom_point(data = data.frame(x = factor(point_groups),
                                 y = size_relativeAbund[selected_groups & diversity>div_threshold][-69]
                                 [dominant_function0 %in% point_groups]),
               aes(x,y)) +
    theme_bw() +
    # ggtitle(LETTERS[8]) +
    # scale_y_log10() +
    # geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.title=element_text(size=22),
          axis.text=element_text(size=22),
          plot.title=element_text(hjust=0, size=24),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin=unit(c(0.5,1,-5,5),"mm")) +
    labs(x="", y="Mean size (micron)")
  grid::grid.draw(boxplot.body.size.functions)
  dev.off()
  # ggplotGrob(
  
  t.test(log10(size_relativeAbund)[selected_groups & diversity>div_threshold][dominant_function0 %in% c("Bacillariophyta","Other phototrophs")],
         log10(size_relativeAbund)[selected_groups & diversity>div_threshold][dominant_function0 %in% c("Gel. carn. filterers","Other metazoa","Pteropoda","Copepoda")])
  t.test(log10(size_relativeAbund)[selected_groups & diversity>div_threshold][dominant_function0 %in% c("Bacillariophyta","Other phototrophs")],
         log10(size_relativeAbund)[selected_groups & diversity>div_threshold][dominant_function0 %in% "Phagotrophs"])
  t.test(log10(size_relativeAbund)[selected_groups & diversity>div_threshold][dominant_function0 %in% c("Bacillariophyta","Other phototrophs")],
         log10(size_relativeAbund)[selected_groups & diversity>div_threshold][dominant_function0 %in% "Parasites"])
  t.test(log10(size_relativeAbund)[selected_groups & diversity>div_threshold][dominant_function0 %in% "Phagotrophs"],
         log10(size_relativeAbund)[selected_groups & diversity>div_threshold][dominant_function0 %in% c("Gel. carn. filterers","Other metazoa","Pteropoda","Copepoda")])
  t.test(log10(size_relativeAbund)[selected_groups & diversity>div_threshold][dominant_function0 %in% "Parasites"],
         log10(size_relativeAbund)[selected_groups & diversity>div_threshold][dominant_function0 %in% c("Gel. carn. filterers","Other metazoa","Pteropoda","Copepoda")])
  t.test(log10(size_relativeAbund)[selected_groups & diversity>div_threshold][dominant_function0 %in% "Parasites"],
         log10(size_relativeAbund)[selected_groups & diversity>div_threshold][dominant_function0 %in% "Phagotrophs"])
  
}

# Fig. S. - effect of ecological categories on axis 2 accounting for body size
{
  x = log10(size_relativeAbund)[selected_groups & !taxo_groups %in% c("Dinophyceae","Collodaria","Porifera")]
  # g1 = factor(dominant_function1,c("Phagotrophs","Metazoans","Parasites","Phototrophs"))
  g = factor(dominant_function1[names(dominant_function1) != "Porifera"],c("Phototrophs","Phagotrophs","Metazoans","Parasites"))
  # g = factor(dominant_function1,c("Phagotrophs","Phototrophs","Metazoans","Parasites"))
  # g = factor(dominant_function1,c("Metazoans","Phagotrophs","Phototrophs","Parasites"))
  # g = factor(dominant_function1,c("Parasites","Metazoans","Phagotrophs","Phototrophs"))
  # g = dominant_function1
  # To switch between default contrasts (i.e., no contrast, or "treatment contrasts") and Hebert contrasts:
  # options(contrasts=c("contr.helmert", "contr.poly"))
  # options(contrasts=c("contr.treatment", "contr.poly"))
  Phago.vs.rest = c(-3,1,1,1)
  Photo.vs.rest = c(1,-3,1,1)
  Meta.vs.rest = c(1,1,-3,1)
  Para.vs.rest = c(1,1,1,-3)
  # contrasts(g) = cbind(Para.vs.rest,Photo.vs.rest,Meta.vs.rest,Phago.vs.rest)
  # contrasts(g) = cbind(Para.vs.rest,Phago.vs.rest)
  # contrasts(g) = cbind(Meta.vs.rest,Photo.vs.rest,Phago.vs.rest)
  # contrasts(g) = cbind(Photo.vs.rest,Meta.vs.rest)
  #
  # The contrast c(0,1,-1,0) actually compares the Meta mean with the mean of Meta and Photo
  Meta.vs.Photo = c(0,1,-1,0)
  Phago.vs.Meta = c(1,0,-1,0)
  Phago.vs.Photo = c(1,-1,0,0)
  #
  Meta.vs.Photo = c(0,1,-1,0)
  Para.vs.Phago = c(1,0,0,-1)
  Meta.Para.vs.Photo.Phago = c(1,1,-1,-1)
  #
  Para.vs.Meta = c(0,0,1,-1)
  Para.vs.Photo = c(0,1,0,-1)
  Para.vs.Phago = c(1,0,0,-1)
  #
  contrasts(g) = cbind(Meta.vs.Photo,Phago.vs.Meta,Phago.vs.Photo)
  contrasts(g) = cbind(Meta.vs.Photo,Para.vs.Phago,Meta.Para.vs.Photo.Phago)
  contrasts(g) = Meta.vs.Photo
  contrasts(g) = Para.vs.Phago
  # contrasts(g) = cbind(Para.vs.Meta,Para.vs.Photo,Para.vs.Phago)
  # contrasts(g) = cbind(Meta.vs.Photo,Phago.vs.Photo,Para.vs.Photo)
  y = NormalizedVI_pcoa[[3]][!names(NormalizedVI_pcoa[[3]][,2]) %in% c("Dinophyceae","Collodaria","Porifera"),2]
  ####################
  # aov(): type I sum of squares, sensitive to ordering:
  ancova = aov(y ~ g + x)
  # CList = list("Para.vs.rest" = 1,
  #              "Photo.vs.rest" = 2,
  #              "Meta.vs.rest" = 3,
  #              "Phago.vs.rest" = 4)
  # CList = list("Para.vs.rest" = 1,
  #              "Phago.vs.rest" = 2)
  # CList = list("Meta.vs.rest" = 1,
  #              "Photo.vs.rest" = 2,
  #              "Phago.vs.rest" = 3)
  CList = list("Meta.vs.Photo" = 1,
               "Phago.vs.Meta" = 2,
               "Phago.vs.Photo" = 3)
  CList = list("Meta.vs.Photo" = 1,
               "Para.vs.Phago" = 2,
               "Meta.Para.vs.Photo.Phago" = 3)
  CList = list("Meta.vs.Photo" = 1)
  CList = list("Para.vs.Phago" = 1)
  # CList = list("Para.vs.Meta" = 1,
  #              "Para.vs.Photo" = 2,
  #              "Para.vs.Phago" = 3)
  # CList = list("Meta.vs.Photo" = 1,
  #              "Phago.vs.Photo" = 2,
  #              "Para.vs.Photo" = 3)
  summary(ancova,
          split=list(g=CList))
  # plot(ancova) allows to check the assumptions of the ANOVA 
  # (such as the variance not being significantly different between categories,
  # which can be tested by Fligner-Killeen test: fligner.test(y ~ g))
  #####################
  # car::Anova( , type = "II"): type II sum of squares, not sensitive to ordering:
  # it actually yields the second line of summary(aov(x + g)) and of summary(aov(g + x)).
  # Contrasts are needed to perform F-tests between specific groups (interpretation not easy though).
  car::Anova(ancova, type = "II")
  # Unlike aov and car::Anova, lm() shows coefficients, effect sizes and t-tests for scoefficients, 
  # but only a single global F-test:
  anova = lm(y ~ g)
  summary(anova)
  summary.lm(aov(y ~ x + g))
  # summary(lm()) gives the same result as summary.lm(aov());
  # anova(lm()) or summary.aov(lm()) gives the same result as sumary(aov()).
  # anova-aov focus on variance (F-tests), lm on coefficients (effect sizes, t-tests), but the model is the same.
  # anova(model1,model2) can be used to compare the variance explained by two models with an F-test.
  # -> if the F-test in non-significant, the variance explained by the two models can be considered the same, 
  # and so the simpler model is to be preferred. AIC can also be computed for the two models using AIC(model).
  # t-tests (lm) should be favored for detailed interpretation of levels once the F-test is significant.
  ancova = lm(y ~ x + g)
  summary(ancova)
  
  decreasing_sizes = sort.int(size_relativeAbund[selected_groups & !taxo_groups %in% c("Dinophyceae","Collodaria","apusozoan_*") & diversity>div_threshold],decreasing=F,index.return = T)
  size_class = vector(length = length(decreasing_sizes$ix), mode = "character")
  nb_classes = 2
  
  if (nb_classes == 4)
  {
    # 4 size classes:
    size_class[decreasing_sizes$ix[1:17]] = "24-35"
    size_class[decreasing_sizes$ix[18:(17*2)]] = "35-61"
    size_class[decreasing_sizes$ix[(17*2+1):(17*3)]] = "61-208"
    size_class[decreasing_sizes$ix[(17*3+1):67]] ="208-731"
    size_class = factor(size_class,c("24-35","35-61","61-208","208-731"))
  } else if (nb_classes == 3)
  {
    # 3 size classes:
    size_class[decreasing_sizes$ix[1:23]] = "24-40"
    size_class[decreasing_sizes$ix[24:(22+23)]] = "40-94"
    size_class[decreasing_sizes$ix[(22+23+1):67]] = "94-731"
    size_class = factor(size_class,c("24-40","40-94","94-731"))
  } else if (nb_classes == 2)
  {
    # 2 size classes:
    size_class[decreasing_sizes$ix[1:34]] = "24-61"
    size_class[decreasing_sizes$ix[35:67]] = "61-731"
    size_class = factor(size_class,c("24-61","61-731"))
  }
  
  axis2.functional.boxplot <- function(indices)
  {
    boxplot.2nd.axis.functions = ggplot(data = data.frame(x = dominant_function1[indices],
                                                          y = NormalizedVI_pcoa[[3]][!names(NormalizedVI_pcoa[[3]][,2]) %in% c("Dinophyceae","Collodaria"),2][indices])) +
      scale_x_discrete(limits=levels(dominant_function1)) +
      geom_boxplot(aes(x,y)) +
      # scale_y_log10() +
      # geom_point(data = data.frame(x = factor(point_groups),
      #                              y = NormalizedVI_pcoa[[3]][,1]
      #                              [dominant_function0 %in% point_groups]),
      # aes(x,y)) +
      theme_bw() +
      # ggtitle(LETTERS[8]) +
      # scale_y_log10() +
      # geom_hline(yintercept = 1, linetype = "dashed") +
      theme(axis.title=element_text(size=22),
            axis.text=element_text(size=22),
            plot.title=element_text(hjust=0, size=24),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            plot.margin=unit(c(0.5,1,-5,5),"mm")) +
      labs(x="", y="PCoA axis 2")
  }
  
  i = 0
  functional.boxplot = list()
  for (size in levels(size_class))
  {
    i = i+1
    functional.boxplot[[i]] = axis2.functional.boxplot(which(size_class == size))
  }
  
  g = list()
  for (i in 1:nb_classes)
  {
    g[[i]] = ggplotGrob(functional.boxplot[[i]] + ggtitle(paste(LETTERS[i],"-",levels(size_class)[i])))
  }
  if (nb_classes == 4)
  {
    pdf(paste0(figure_folder,"/Normalized.PCoA.axis2.vs.function.within.4.size.classes_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
        width=7.5*2.2/1.2,height=12/3*4/1.2)
    grid.arrange(grobs = list(g[[1]],g[[2]],g[[3]],g[[4]]),
                 # widths = rep(1,7),
                 layout_matrix = rbind(c(1, 2),c(3, 4)))
    dev.off()
  } else if (nb_classes == 3)
  {
    pdf(paste0(figure_folder,"/Normalized.PCoA.axis2.vs.function.within.3.size.classes_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
        width=7.5*2.2/1.2*3/2,height=12/3*2/1.2)
    grid.arrange(grobs = list(g[[1]],g[[2]],g[[3]]),
                 # widths = rep(1,7),
                 layout_matrix = t(matrix(c(1,2,3))))
    dev.off()
  } else if (nb_classes == 2)
  {
    pdf(paste0(figure_folder,"/Normalized.PCoA.axis2.vs.function.within.2.size.classes_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
        width=7.5*2.2/1.2,height=12/3*2/1.2)
    grid.arrange(grobs = list(g[[1]],g[[2]]),
                 # widths = rep(1,7),
                 layout_matrix = t(matrix(c(1,2))))
    dev.off()
  }
}

# Fig. S - axis 2 vs. body size within ecological categories
{
  x = size_relativeAbund[selected_groups & !taxo_groups %in% c("Dinophyceae","Collodaria")][dominant_function1 == "Phototrophs"]
  y = NormalizedVI_pcoa[[3]][,2][!rownames(NormalizedVI_pcoa[[3]]) %in% c("Dinophyceae","Collodaria")][dominant_function1 == "Phototrophs"]
  plot.phototrophs.norm.VI.PCoA.axis2.vs.size = cor.plot(x = x,
                                                       y = y,
                                                       y.lab = "Biogeographic axis 2",
                                                       x.lab = expression("Body size ("*mu*"m), log scale"),
                                                       y.lab.hjust = 0.5,
                                                       x.log = T,
                                                       y.log = F,
                                                       fit = T,
                                                       # fit.display = "spearman.p",
                                                       fit.display = "pearson.p",
                                                       x.cor.pos = 0.8,
                                                       y.cor.pos = 0.1,
                                                       mar.vect = c(5,5,5,5))
  g1 = ggplotGrob(plot.phototrophs.norm.VI.PCoA.axis2.vs.size + 
                    # ggtitle(paste(LETTERS[1],"-","Phototrophs")))
                    # ylim(-0.2,0.2) +
                    ggtitle("Phototrophs"))
  
  x = size_relativeAbund[selected_groups & diversity>div_threshold & !taxo_groups %in% c("Dinophyceae","Collodaria")][dominant_function1 == "Phagotrophs"]
  y = NormalizedVI_pcoa[[3]][,2][!rownames(NormalizedVI_pcoa[[3]]) %in% c("Dinophyceae","Collodaria")][dominant_function1 == "Phagotrophs"]
  plot.phagotrophs.norm.VI.PCoA.axis2.vs.size = cor.plot(x = x,
                                                         y = y,
                                                         # excluded.points = which(names(dominant_function1[dominant_function1 == "Phagotrophs"]) == "Phaeodaria"),
                                                         y.lab = "Biogeographic axis 2",
                                                         x.lab = expression("Body size ("*mu*"m), log scale"),
                                                         y.lab.hjust = 0.5,
                                                         x.log = T,
                                                         y.log = F,
                                                         fit = T,
                                                         # fit.display = "spearman.p",
                                                         fit.display = "pearson.p",
                                                         x.cor.pos = 0.8,
                                                         y.cor.pos = 0.1,
                                                         mar.vect = c(5,5,1,5))
  g2 = ggplotGrob(plot.phagotrophs.norm.VI.PCoA.axis2.vs.size + 
                    # ggtitle(paste(LETTERS[2],"-","Phagotrophs")))
                    ggtitle("Phagotrophs") +
                    # ylim(-0.2,0.2) +
                    theme(axis.title.y = element_blank()))
  
  x = size_relativeAbund[selected_groups & diversity>div_threshold & !taxo_groups %in% c("Dinophyceae","Collodaria")][dominant_function1 == "Metazoans"]
  y = NormalizedVI_pcoa[[3]][,2][!rownames(NormalizedVI_pcoa[[3]]) %in% c("Dinophyceae","Collodaria")][dominant_function1 == "Metazoans"]
  plot.metazoans.norm.VI.PCoA.axis2.vs.size = cor.plot(x = x,
                                                         y = y,
                                                         excluded.points = which(names(dominant_function1[dominant_function1 == "Metazoans"]) == "Porifera"),# c("Porifera","Chaetognatha")),
                                                         y.lab = "Biogeographic axis 2",
                                                         x.lab = expression("Body size ("*mu*"m), log scale"),
                                                         y.lab.hjust = 0.5,
                                                         x.log = T,
                                                         y.log = F,
                                                         fit = T,
                                                         # fit.display = "spearman.p",
                                                         fit.display = "pearson.p",
                                                         x.cor.pos = 0.8,
                                                         y.cor.pos = 0.2,
                                                         mar.vect = c(5,5,1,5))
  g3 = ggplotGrob(plot.metazoans.norm.VI.PCoA.axis2.vs.size + 
                    # ylim(-0.2,0.2) +
                    # ggtitle(paste(LETTERS[3],"-","Metazoans")) +
                    ggtitle("Metazoans"))
  
  x = size_relativeAbund[selected_groups & !taxo_groups %in% c("Dinophyceae","Collodaria")][dominant_function1 == "Parasites"]
  y = NormalizedVI_pcoa[[3]][,2][!rownames(NormalizedVI_pcoa[[3]]) %in% c("Dinophyceae","Collodaria")][dominant_function1 == "Parasites"]
  plot.parasites.norm.VI.PCoA.axis2.vs.size = cor.plot(x = x,
                                                       y = y,
                                                       # excluded.points = which(names(dominant_function1[dominant_function1 == "Parasites"]) == "MALV-IV"),
                                                       y.lab = "Biogeographic axis 2",
                                                       x.lab = expression("Body size ("*mu*"m), log scale"),
                                                       y.lab.hjust = 0.5,
                                                       x.log = T,
                                                       y.log = F,
                                                       fit = T,
                                                       # fit.display = "spearman.p",
                                                       fit.display = "pearson.p",
                                                       x.cor.pos = 0.8,
                                                       y.cor.pos = 0.2,
                                                       mar.vect = c(5,5,1,5))
  g4 = ggplotGrob(plot.parasites.norm.VI.PCoA.axis2.vs.size + 
                    # ggtitle(paste(LETTERS[4],"-","Parasites")))
                    ggtitle("Parasites") +
                    # ylim(-0.2,0.2) +
                    theme(axis.title.y = element_blank()))
  
  pdf(paste0(figure_folder,"/FigS_PCoA2.vs.size.within.categories_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1_excluded.group.pdf"),
      width=7.5*2.2/1.2,height=12/3*4/1.2)
  grid.arrange(grobs = list(g1,g2,g3,g4),
               # widths = rep(1,7),
               layout_matrix = rbind(c(1, 2),c(3, 4)))
  dev.off()
  
  cor.test(NormalizedVI_pcoa[[3]][diversity[selected_groups & diversity>div_threshold][-69]<2000,1],log10(diversity)[selected_groups & diversity>div_threshold][-69][diversity[selected_groups & diversity>div_threshold][-69]<2000])
  cor.test(NormalizedVI_pcoa[[3]][,1],log10(size_relativeAbund)[selected_groups & diversity>div_threshold][-69])
  cor.test(I_square.observed_w.mean[selected_groups & diversity>div_threshold,1][-69],log10(size_relativeAbund)[selected_groups & diversity>div_threshold][-69])
  
  cor.test(NormalizedVI_pcoa[[3]][,2],log10(size_relativeAbund)[selected_groups & diversity>div_threshold][-69])
  cor.test(NormalizedVI_pcoa[[3]][,2],log10(diversity)[selected_groups & diversity>div_threshold][-69])
  
  # Multiple comparisons should be perfomred with pairwise.t.test(x,y) or TukeyHSD(aov(y ~ x))
  t.test(NormalizedVI_pcoa[[3]][,2][dominant_function0 %in% c("Bacillariophyta","Other phototrophs")],
         NormalizedVI_pcoa[[3]][,2][dominant_function0 %in% c("Gel. carn. filterers","Other metazoa","Pteropoda","Copepoda")])
  t.test(NormalizedVI_pcoa[[3]][,2][dominant_function0 %in% c("Bacillariophyta","Other phototrophs")],
         NormalizedVI_pcoa[[3]][,2][dominant_function0 %in% "Phagotrophs"])
  t.test(NormalizedVI_pcoa[[3]][,2][dominant_function0 %in% c("Bacillariophyta","Other phototrophs")],
         NormalizedVI_pcoa[[3]][,2][dominant_function0 %in% "Parasites"])
  t.test(NormalizedVI_pcoa[[3]][,2][dominant_function0 %in% c("Gel. carn. filterers","Other metazoa","Pteropoda","Copepoda")],
         NormalizedVI_pcoa[[3]][,2][dominant_function0 %in% "Parasites"])
  t.test(NormalizedVI_pcoa[[3]][,2][dominant_function0 %in% c("Gel. carn. filterers","Other metazoa","Pteropoda","Copepoda")],
         NormalizedVI_pcoa[[3]][,2][dominant_function0 %in% "Phagotrophs"])
  t.test(NormalizedVI_pcoa[[3]][,2][dominant_function0 %in% "Phagotrophs"],
         NormalizedVI_pcoa[[3]][,2][dominant_function0 %in% "Parasites"])
}

# Fig. S - fractions ratio within ecological categories
{
  # plots = list(boxplot.ratio.currents.env.functions,plot.currents.env.ratio.vs.K,plot.currents.env.ratio.vs.size)
  # grobs = list()  
  # heights = list()
  # for (i in 1:length(plots))
  # {
  #   grobs[[i]] = ggplotGrob(plots[[i]])
  #   heights[[i]] = grobs[[i]]$heights[2:5]
  # }
  # maxheight = do.call(grid::unit.pmax, heights)
  # for (i in 1:length(grobs))
  # {
  #   grobs[[i]]$heights[2:5] = as.list(maxheight)
  # }
  # plot = do.call("grid.arrange", c(grobs, nrow = 1))
  
  #plot = grid.arrange(grobs = list(boxplot.ratio.currents.env.functions,plot.currents.env.ratio.vs.K,plot.currents.env.ratio.vs.size), ncol = 3, layout_matrix = matrix(1:3,nrow=1,byrow=T))
  #ggsave(file_name, plot=plot,width=7.5*2.5,height=12/3*2)
  
  fraction.ratio.vs.log.size <- function(x,y)
  {
    cor.test = cor.test(log10(x),log10(y))
    ggplot(data = data.frame(x = x, y = y)) +
      # geom_point(aes(x,y1)) +
      geom_point(aes(x,y)) +
      scale_x_log10() +
      scale_y_log10() +
      theme_bw() +
      # geom_hline(yintercept = 1, linetype="dashed") +
      # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
      # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
      theme(axis.text = element_text(size=22),
            axis.title=element_text(size=22),
            # axis.title.x=element_text(vjust = 45),
            plot.title=element_text(hjust=0, size=24),
            plot.margin=unit(c(0.5,12,0.5,3),"mm")) +
      #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
      labs(y="Ratio of purely connectivity-explained\n over envir.-explained variance (Surface)", x=expression("Clade mean body size ("*mu*"m)")) +
      geom_smooth(aes(x,y),method='lm',col="black") +
      annotate(geom="text", 
               x=(max(x)/min(x))^0.8*min(x),
               y=(max(y)/min(y))^0.2*min(y),
               label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
               size=7)
  }
  
  i_case = 1
  x = size_relativeAbund[selected_groups & diversity>div_threshold][-69][!is.na(dominant_function0) & dominant_function0 %in% c("Bacillariophyta","Other phototrophs")]
  y = (varpart.env.spatial[[i_case]][3,]/varpart.env.spatial[[i_case]][1,])[selected_groups & diversity>div_threshold][-69][!is.na(dominant_function0) & as.vector(dominant_function0) %in% c("Bacillariophyta","Other phototrophs")]
  plot.phototrophs.fraction.ratio.axis2.vs.size = fraction.ratio.vs.log.size(x,y)
  
  x = size_relativeAbund[selected_groups & diversity>div_threshold][-69][!is.na(dominant_function0) & as.vector(dominant_function0) == "Phagotrophs"]
  y = (varpart.env.spatial[[i_case]][3,]/varpart.env.spatial[[i_case]][1,])[selected_groups & diversity>div_threshold][-69][!is.na(dominant_function0) & as.vector(dominant_function0) == "Phagotrophs"]
  plot.phagotrophs.fraction.ratio.axis2.vs.size = fraction.ratio.vs.log.size(x,y)
  
  x = size_relativeAbund[selected_groups & diversity>div_threshold][-69][!is.na(dominant_function0) & dominant_function0 %in% c("Gel. carn. filterers","Pteropoda","Copepoda","Other metazoa")]
  y = (varpart.env.spatial[[i_case]][3,]/varpart.env.spatial[[i_case]][1,])[selected_groups & diversity>div_threshold][-69][!is.na(dominant_function0) & dominant_function0 %in% c("Gel. carn. filterers","Pteropoda","Copepoda","Other metazoa")]
  plot.metazoans.fraction.ratio.axis2.vs.size = fraction.ratio.vs.log.size(x,y)
  
  x = size_relativeAbund[selected_groups & diversity>div_threshold][-69][!is.na(dominant_function0) & dominant_function0 == "Parasites"]
  y = (varpart.env.spatial[[i_case]][3,]/varpart.env.spatial[[i_case]][1,])[selected_groups & diversity>div_threshold][-69][!is.na(dominant_function0) & dominant_function0 == "Parasites"]
  plot.parasites.fraction.ratio.axis2.vs.size = fraction.ratio.vs.log.size(x,y)
  
  pdf(paste0(figure_folder,"/Varpart.env.currents.fraction.ratio.vs.size.within.functional.group_",div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2,height=12/3*4/1.2)
  g1 = ggplotGrob(plot.phototrophs.fraction.ratio.axis2.vs.size + ggtitle(paste(LETTERS[1],"-","Phototrophs")))
  g2 = ggplotGrob(plot.phagotrophs.fraction.ratio.axis2.vs.size + ggtitle(paste(LETTERS[2],"-","Phagotrophs")))
  g3 = ggplotGrob(plot.metazoans.fraction.ratio.axis2.vs.size + ggtitle(paste(LETTERS[3],"-","Metazoans")))
  g4 = ggplotGrob(plot.parasites.fraction.ratio.axis2.vs.size + ggtitle(paste(LETTERS[4],"-","Parasites")))
  grid.arrange(grobs = list(g1,g2,g3,g4),
               # widths = rep(1,7),
               layout_matrix = rbind(c(1, 2),c(3, 4)))
  dev.off()
  
  t.test(log10(varpart.env.spatial[[1]][3,]/varpart.env.spatial[[1]][1,])[selected_groups & diversity>div_threshold][dominant_function0 %in% c("Bacillariophyta","Other phototrophs")],
         log10(varpart.env.spatial[[1]][3,]/varpart.env.spatial[[1]][1,])[selected_groups & diversity>div_threshold][dominant_function0 %in% c("Gel. carn. filterers","Other metazoa","Pteropoda","Copepoda")])
  t.test(log10(varpart.env.spatial[[1]][3,]/varpart.env.spatial[[1]][1,])[selected_groups & diversity>div_threshold][dominant_function0 %in% c("Bacillariophyta","Other phototrophs")],
         log10(varpart.env.spatial[[1]][3,]/varpart.env.spatial[[1]][1,])[selected_groups & diversity>div_threshold][dominant_function0 %in% "Phagotrophs"])
  t.test(log10(varpart.env.spatial[[1]][3,]/varpart.env.spatial[[1]][1,])[selected_groups & diversity>div_threshold][dominant_function0 %in% c("Bacillariophyta","Other phototrophs")],
         log10(varpart.env.spatial[[1]][3,]/varpart.env.spatial[[1]][1,])[selected_groups & diversity>div_threshold][dominant_function0 %in% "Parasites"])
  t.test(log10(varpart.env.spatial[[1]][3,]/varpart.env.spatial[[1]][1,])[selected_groups & diversity>div_threshold][dominant_function0 %in% c("Gel. carn. filterers","Other metazoa","Pteropoda","Copepoda")],
         log10(varpart.env.spatial[[1]][3,]/varpart.env.spatial[[1]][1,])[selected_groups & diversity>div_threshold][dominant_function0 %in% "Parasites"])
  t.test(log10(varpart.env.spatial[[1]][3,]/varpart.env.spatial[[1]][1,])[selected_groups & diversity>div_threshold][dominant_function0 %in% c("Gel. carn. filterers","Other metazoa","Pteropoda","Copepoda")],
         log10(varpart.env.spatial[[1]][3,]/varpart.env.spatial[[1]][1,])[selected_groups & diversity>div_threshold][dominant_function0 %in% "Phagotrophs"])
  t.test(log10(varpart.env.spatial[[1]][3,]/varpart.env.spatial[[1]][1,])[selected_groups & diversity>div_threshold][dominant_function0 %in% "Phagotrophs"],
         log10(varpart.env.spatial[[1]][3,]/varpart.env.spatial[[1]][1,])[selected_groups & diversity>div_threshold][dominant_function0 %in% "Parasites"])
}

# Fig. S - biotic-abiotic fraction ratio 
{
  # barplot_data.frame = rbind(varpart.biotic.abiotic[[i_case]][1,],varpart.biotic.abiotic[[i_case]][2,],varpart.biotic.abiotic[[i_case]][3,])
  # dimnames(barplot_data.frame) = list(c("Purely by abiotic conditions","Jointly by biotic and abiotic conditions","Purely by biotic conditions"), taxo_groups_unmodified)
  
  g1 = g2 = g3 = g4 = g5 = list()
  for (i_case in 1:2)
  {
    x = NormalizedVI_pcoa[[3]][,1]
    # y = (varpart.biotic.abiotic[[i_case]][3,]/varpart.biotic.abiotic[[i_case]][1,])[selected_groups & diversity>div_threshold]
    y = (colSums(varpart.biotic.abiotic[[i_case]][2:3,])/colSums(varpart.biotic.abiotic[[i_case]][1:2,]))[selected_groups & diversity>div_threshold]
    plot.biotic.abiotic.ratio.vs.PCoA.axis1 = cor.plot(x = x,
                                                       y = y,
                                                       # y.lab = paste("Ratio of variance purely explained\n by local biotic over abiotic conditions"),#c("","(DCM)")[i_case]),
                                                       y.lab = paste("Ratio of variance explained\nby biotic over abiotic conditions (log scale)"),#c("","(DCM)")[i_case]),
                                                       x.lab = "PCoA axis 1",
                                                       y.lab.hjust = 0.5,
                                                       x.log = F,
                                                       y.log = T,
                                                       fit = T,
                                                       fit.display = "spearman.p",
                                                       x.cor.pos = ifelse(i_case == 1,0.8,0.1),
                                                       y.cor.pos = ifelse(i_case == 1,0.1,0.1),
                                                       mar.vect = c(5,5,1,5))
    g1[[i_case]] = ggplotGrob(plot.biotic.abiotic.ratio.vs.PCoA.axis1 +
                                geom_hline(yintercept = 1, linetype="dashed"))
    
    x = NormalizedVI_pcoa[[3]][,2]
    # y = (varpart.biotic.abiotic[[i_case]][3,]/varpart.biotic.abiotic[[i_case]][1,])[selected_groups & diversity>div_threshold]
    y = (colSums(varpart.biotic.abiotic[[i_case]][2:3,])/colSums(varpart.biotic.abiotic[[i_case]][1:2,]))[selected_groups & diversity>div_threshold]
    plot.biotic.abiotic.ratio.vs.PCoA.axis2 = cor.plot(x = x,
                                                         y = y,
                                                         # y.lab = paste("Ratio of variance purely explained\n by local biotic over abiotic conditions"),#c("","(DCM)")[i_case]),
                                                         y.lab = paste("Ratio of variance explained\n by biotic over abiotic conditions (log scale)"),#c("","(DCM)")[i_case]),
                                                         x.lab = "PCoA axis 2",
                                                         y.lab.hjust = 0.5,
                                                         x.log = F,
                                                         y.log = T,
                                                         fit = T,
                                                         fit.display = "spearman.p",
                                                         x.cor.pos = ifelse(i_case == 1,0.8,0.1),
                                                         y.cor.pos = ifelse(i_case == 1,0.1,0.1),
                                                         mar.vect = c(5,5,1,5))
    g2[[i_case]] = ggplotGrob(plot.biotic.abiotic.ratio.vs.PCoA.axis2 + 
                                # geom_vline(xintercept = 1, linetype="dashed"))
                                geom_hline(yintercept = 1, linetype="dashed") +
                                theme(axis.title.y = element_blank()))
    
    x = diversity[selected_groups & diversity>div_threshold]
    # y = (varpart.biotic.abiotic[[i_case]][3,]/varpart.biotic.abiotic[[i_case]][1,])[selected_groups & diversity>div_threshold]
    y = (colSums(varpart.biotic.abiotic[[i_case]][2:3,])/colSums(varpart.biotic.abiotic[[i_case]][1:2,]))[selected_groups & diversity>div_threshold]
    fraction.biotic.abiotic.ratio.vs.diversity = cor.plot(x = x,
                                                          y = y,
                                                          # y.lab = paste("Ratio of variance purely explained\n by local biotic over abiotic conditions"),#c("","(DCM)")[i_case]),
                                                          y.lab = paste("Ratio of variance explained\n by biotic over abiotic conditions (log scale)"),#c("","(DCM)")[i_case]),
                                                          x.lab = "Diversity (#OTUs), log scale",
                                                          y.lab.hjust = 0.5,
                                                          x.log = T,
                                                          y.log = T,
                                                          fit = T,
                                                          fit.display = "spearman.p",
                                                          x.cor.pos = 0.8,
                                                          y.cor.pos = 0.85,
                                                          mar.vect = c(5,12,1,5))
    g3[[i_case]] = ggplotGrob(fraction.biotic.abiotic.ratio.vs.diversity + geom_hline(yintercept = 1, linetype="dashed"))
    
    x = size_relativeAbund[selected_groups & diversity>div_threshold]
    # y = (varpart.biotic.abiotic[[i_case]][3,]/varpart.biotic.abiotic[[i_case]][1,])[selected_groups & diversity>div_threshold]
    y = (colSums(varpart.biotic.abiotic[[i_case]][2:3,])/colSums(varpart.biotic.abiotic[[i_case]][1:2,]))[selected_groups & diversity>div_threshold]
    fraction.biotic.abiotic.ratio.vs.size = cor.plot(x = x,
                                                     y = y,
                                                     # y.lab = paste("Ratio of variance purely explained\n by local biotic over abiotic conditions"),#c("","(DCM)")[i_case]),
                                                     y.lab = paste("Ratio of variance explained\n by biotic over abiotic conditions (log scale)"),#c("","(DCM)")[i_case]),
                                                     x.lab = expression("Body size ("*mu*"m), log scale"),
                                                     y.lab.hjust = 0.5,
                                                     x.log = T,
                                                     y.log = T,
                                                     fit = T,
                                                     fit.display = "spearman.p",
                                                     x.cor.pos = 0.8,
                                                     y.cor.pos = 0.85,
                                                     mar.vect = c(5,5,1,5))
    g4[[i_case]] = ggplotGrob(fraction.biotic.abiotic.ratio.vs.size + 
                                geom_hline(yintercept = 1, linetype="dashed") +
                                theme(axis.title.y = element_blank()))
    
    # values = (varpart.biotic.abiotic[[i_case]][3,]/varpart.biotic.abiotic[[i_case]][1,])[selected_groups & diversity>div_threshold & !colnames(varpart.biotic.abiotic[[i_case]]) %in% c("Dinophyceae","Collodaria")]
    values = (colSums(varpart.biotic.abiotic[[i_case]][2:3,])/colSums(varpart.biotic.abiotic[[i_case]][1:2,]))[selected_groups & diversity>div_threshold & !colnames(varpart.biotic.abiotic[[i_case]]) %in% c("Dinophyceae","Collodaria")]
    fraction.biotic.abiotic.ratio.vs.function = box.plot(classes = factor(dominant_function1, c("Phototrophs","Phagotrophs","Metazoans","Parasites")),
                                                    values = values,
                                                    # y.lab = paste("Ratio of variance purely explained\n by local biotic over abiotic conditions"),#c("","(DCM)")[i_case]),
                                                    y.lab = paste("Ratio of variance explained\n by biotic over abiotic conditions (log scale)"),#c("","(DCM)")[i_case]),
                                                    y.log = T,
                                                    y.lab.hjust = NULL,
                                                    fit = T,
                                                    mar.vect = c(5,5,1,5))
    g5[[i_case]] = ggplotGrob(fraction.biotic.abiotic.ratio.vs.function + geom_hline(yintercept = 1, linetype="dashed"))
  }
  
  # pdf(paste0(figure_folder,"/FigS_biotic.abiotic.ratio_",div_threshold,"plusOTUs_selected100+1_indiv.signif.axes_no.var.select.pdf"),
  # pdf(paste0(figure_folder,"/FigS_pure+mixed_biotic.abiotic.ratio_",div_threshold,"plusOTUs_selected100+1_no.indiv.signif.axes_no.var.select.pdf"),
  pdf(paste0(figure_folder,"/FigS_pure+mixed_biotic.abiotic.ratio_",div_threshold,"plusOTUs_selected100+1_indiv.signif.axes.BH.correction.pdf"),
      # pdf(paste0(figure_folder,"/FigS_pure_biotic.abiotic.ratio_",div_threshold,"plusOTUs_selected100+1_indiv.signif.axes.BH.correction.pdf"),
      # pdf(paste0(figure_folder,"/FigS_biotic.abiotic.ratio_",div_threshold,"plusOTUs_selected100+1_apusozoans.pdf"),
      # pdf(paste0(figure_folder,"/FigS_biotic.abiotic.ratio_pure+mixed_",div_threshold,"plusOTUs_selected100+1_apusozoans.pdf"),
      # width=7.5*2.2/1.2*3/2,height=12/3*4/1.2)
      width=7.5*2.2/1.2/2*2,height=12/3*4/1.2/2*3)
  for (i_case in 1:2)
  {
    grid.arrange(grobs = list(g1[[i_case]],g2[[i_case]],g3[[i_case]],g4[[i_case]],g5[[i_case]]),
                 # widths = rep(1,7),
                 # heights = c(0.8,0.8,1),
                 layout_matrix = rbind(c(1, 2),c(3, 4),c(5, 6)))
  }
  dev.off()
  
  y = log10(varpart.biotic.abiotic[[i_case]][3,]/varpart.biotic.abiotic[[i_case]][1,])[selected_groups & diversity>div_threshold]
  y[is.infinite(y)] = NA
  
  t.test(y[!is.na(y) & dominant_function0 %in% c("Bacillariophyta","Other phototrophs")],
         y[!is.na(y) & dominant_function0 %in% c("Gel. carn. filterers","Other metazoa","Pteropoda","Copepoda")])
  t.test(y[!is.na(y) & dominant_function0 %in% c("Bacillariophyta","Other phototrophs")],
         y[!is.na(y) & dominant_function0 %in% "Phagotrophs"])
  t.test(y[!is.na(y) & dominant_function0 %in% c("Bacillariophyta","Other phototrophs")],
         y[!is.na(y) & dominant_function0 %in% "Parasites"])
  t.test(y[!is.na(y) & dominant_function0 %in% "Phagotrophs"],
         y[!is.na(y) & dominant_function0 %in% c("Gel. carn. filterers","Other metazoa","Pteropoda","Copepoda")])
  t.test(y[!is.na(y) & dominant_function0 %in% "Parasites"],
         y[!is.na(y) & dominant_function0 %in% c("Gel. carn. filterers","Other metazoa","Pteropoda","Copepoda")])
  t.test(y[!is.na(y) & dominant_function0 %in% "Parasites"],
         y[!is.na(y) & dominant_function0 %in% "Phagotrophs"])
}

# Fig. A1.2 - assemblage composition
{
  data.folder = paste0(data_folder,"/18S_V9_TARA_CompleteSizeRange_byStationByDepth_AllTaxa",noArcticNoBiomark_insert,noLagoon_insert)
  assemblage_composition = readRDS(paste0(data.folder,"/Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics16_nb_iter1000_nb_real100_meanPosteriorDistributedLlh_thin25_burnin2000_occurrence/1st_closestToMean_realization/assemblage_composition.rds"))
  
  functional.barplot = T
  taxo.barplot = F
  two.assemblage.richness = F
  
  # Setting the threshold for OTU presence-absence:
  threshold = vector(mode = "numeric", length = K)
  for (k in 1:K)
    threshold[k] = sort(assemblage_composition[,k],decreasing=T)[100000]*2
  OTU_presence = t(apply(assemblage_composition[,1:K],1, function(g) g >= threshold))
  # Setting OTUs with abundance below the threshold to 0:
  assemblage_composition[,1:K][!OTU_presence] = 0
  colnames(assemblage_composition) = c(assemblage_labels,colnames(assemblage_composition)[(K+1):ncol(assemblage_composition)])
  # Renormalizing assemblages:
  assemblage_composition[,1:K] = sweep(assemblage_composition[,1:K],2,colSums(assemblage_composition[,1:K]),"/")
  # Computing endemics and ubiquists:
  endemic_OTUs = unlist(lapply(apply(OTU_presence,1,which),length)) == 1
  # endemicity_assemblage = unlist(lapply(apply(OTU_presence,1,which),function(g) g[1]))
  # ubiquist_OTUs = apply(OTU_presence,1,all)
  assemblage_richness = unlist(lapply(apply(OTU_presence,2,which),length))
  assemblage_endemic_richness = unlist(lapply(apply(OTU_presence & endemic_OTUs,2,which),length))
  # OTU richness shared between 2 assemblages, for each assemblage:
  if (two.assemblage.richness)
  {
    twoAssemblage_OTUs = unlist(lapply(apply(OTU_presence,1,which),length)) == 2
    assemblage_shared_richness = matrix(nrow = K, ncol = K, data = 0)
    for (k in 1:K)
    {
      for (k1 in 1:K)
      {
        if (k1 != k)
          assemblage_shared_richness[k1,k] = length(which(twoAssemblage_OTUs & OTU_presence[,k] & OTU_presence[,k1]))
        else 
          assemblage_shared_richness[k1,k] = 0
      }
    }
  }
  # Computing the composition of assemblages in terms of taxonomic groups:
  if (taxo.barplot)
  {
    assemblage_group_composition = matrix(nrow = length(taxo_groups), ncol = K, dimnames = list(taxo_groups_unmodified,assemblage_labels), data = 0)
    # assemblage_group_endemic_proportion = matrix(nrow = length(taxo_groups), ncol = K, dimnames = list(taxo_groups_unmodified,assemblage_labels), data = 0)
    # assemblage_group_ubiquist_proportion = matrix(nrow = length(taxo_groups), ncol = K, dimnames = list(taxo_groups_unmodified,assemblage_labels), data = 0)
    for (group in taxo_groups_unmodified)
    {
      i_group = which(taxo_groups_unmodified == group)
      assemblage_group_composition[i_group,] =  colSums(assemblage_composition[assemblage_composition$taxogroup2 %in% group,1:K])
      # for (k in 1:K)
      # {
      #   assemblage_group_endemic_proportion[i_group,k] = sum(assemblage_composition[assemblage_composition$taxogroup2 %in% group & endemic_OTUs & endemicity_assemblage == k,k])
      #   assemblage_group_ubiquist_proportion[i_group,k] = sum(assemblage_composition[assemblage_composition$taxogroup2 %in% group & ubiquist_OTUs,k])
      # }
    }
    group_abundances = setNames(rowSums(assemblage_group_composition),rownames(assemblage_group_composition))
    selected_OTU_proportion_all = sort.int(group_abundances[selected_groups],decreasing = T,index.return = T)$ix
  }
  # Computing the composition of assemblages in terms of functions:
  if (functional.barplot)
  {
    # functions = c("phototroph","photohost","endophotosymbiont","phagotroph","copepoda","pteropoda","gelatineous_carnivores_filterers","other metazoa","parasite","unknown")
    # function_names = c("Phototroph","Photohost","Endophotosymbiont","Phagotroph","Copepoda","Pteropoda","Gel. carn. filterers","Other metazoa","Parasite","Unknown")
    assemblage_composition$Function[assemblage_composition$Function %in% "endophotosymbiont"] = "unknown"
    assemblage_composition$Function[assemblage_composition$Function %in% "pteropoda"] = "other metazoa"
    assemblage_composition$Function[assemblage_composition$taxogroup2 %in% "Dinophyceae"] = "dino"
    functions = c("phototroph","photohost","phagotroph","dino","copepoda","gelatineous_carnivores_filterers","other metazoa","parasite","unknown")
    function_names = c("Phototrophs","Collodaria","Phagotrophs","Dinophyceae","Copepoda","Gel. carn. filterers","Other metazoa","Parasites","Unknown")
    assemblage_function_composition = matrix(nrow = length(functions), ncol = K, dimnames = list(function_names,assemblage_labels), data = 0)
    # assemblage_function_endemic_proportion = matrix(nrow = length(functions), ncol = K, dimnames = list(function_names,assemblage_labels), data = 0)
    # assemblage_function_ubiquist_proportion = matrix(nrow = length(functions), ncol = K, dimnames = list(function_names,assemblage_labels), data = 0)
    for (fun in functions)
    {
      i_group = which(functions == fun)
      assemblage_function_composition[i_group,] =  colSums(assemblage_composition[assemblage_composition$Function %in% fun,1:K])
      # for (k in 1:K)
      # {
      #   assemblage_function_endemic_proportion[i_group,k] = sum(assemblage_composition[assemblage_composition$Function %in% fun & endemic_OTUs & endemicity_assemblage == k,k])
      #   assemblage_function_ubiquist_proportion[i_group,k] = sum(assemblage_composition[assemblage_composition$Function %in% fun & ubiquist_OTUs,k])
      # }
    }
    # function_colors = c("darkgreen","chartreuse2","firebrick2","firebrick","cadetblue","darkblue","darkturquoise","dodgerblue1","darkgoldenrod1","grey")
    function_colors = c("darkviolet","darkblue","firebrick","firebrick2",
                        colorRampPalette(c("darkgreen","chartreuse2"),space = "Lab")(4)[c(1,3,4)],
                        "darkgoldenrod1","grey")
  }
  
  # Assemblage colors and order:
  assemblage_reordering = c(1,3,9,7,
                            10,8,4,6,
                            11,16,15,12,
                            2,13,5,14)
  color.pal1 = colorRampPalette(c("#7F0000","red","darkorange1","darkgoldenrod1","yellow"),space = "Lab")
  color.pal2 = colorRampPalette(c("#7FFF7F","darkgreen"),space = "Lab")
  color.pal3 = colorRampPalette(c("azure4","antiquewhite3","aliceblue"),space = "Lab")
  color.pal4 = colorRampPalette(c("#00007F","#007FFF","cyan"),space = "Lab")
  col_assemblages = c(color.pal1(6),color.pal2(3),color.pal3(3),rev(color.pal4(4)))
  
  # Colours of main taxonomic groups:
  colours = c("burlywood1",#"darkgoldenrod1",
              colorRampPalette(c("darkgreen","chartreuse2"),space = "Lab")(4),
              #"cornsilk1",
              "mintcream",
              colorRampPalette(c("darkblue","darkturquoise"),space = "Lab")(4),
              colorRampPalette(c("firebrick2","darkgoldenrod1"),space = "Lab")(5),
              colorRampPalette(c("darkviolet","deeppink"),space = "Lab")(3),
              # "burlywood1",
              "darkseagreen",
              "grey")
  
  pdf(paste0(figure_folder,"/FigA1.2_selected100+1.pdf"), height = 14, width = 9.5)
  par(mfrow=c(2,1))
  # pdf(paste0(figure_folder,"/Assemblage_richness_AllTaxa_K",K,"_withEndemics.pdf"))
  par(mar=c(0.1,4.5,1.1,13.6))
  if (two.assemblage.richness)
  {
    barplot = rbind(assemblage_endemic_richness,assemblage_shared_richness[assemblage_reordering,],assemblage_richness-assemblage_endemic_richness-colSums(assemblage_shared_richness))[,assemblage_reordering]
  } else
  {
    # par(mar=c(0.1,4.1,1.1,11.2))
    barplot = rbind(assemblage_endemic_richness,assemblage_richness-assemblage_endemic_richness)[,assemblage_reordering]
  }
  barplot(barplot,
          legend.text=F,
          xaxt="n", cex.axis = 1.5, space = 0.5)
  for (k in 1:K)
  {
    barplot.k = barplot
    barplot.k[,-k] = 0
    barplot(barplot.k,
            legend.text=F,
            col = if (two.assemblage.richness) c(col_assemblages[k],col_assemblages[1:K],"black")
            else rep(col_assemblages[k],2),
            xaxt="n", yaxt="n", space = 0.5, add=T)
    barplot(barplot.k,
            legend.text=F,
            col = "black",
            density = if (two.assemblage.richness) c(0,rep(20,K),0) else c(20,0),
            xaxt="n", yaxt="n", space = 0.5, add=T)
  }
  title(ylab="#OTUs in assemblage",cex.lab=1.5)
  #
  par(mar=c(1.5,4.5,2.1,13.6))
  if (functional.barplot)
  {
    x = barplot(assemblage_function_composition[nrow(assemblage_function_composition):1,assemblage_reordering],
                col=rev(function_colors),
                legend.text = T, xaxt="n", space = 0.5, cex.axis = 1.5,
                args.legend = list(bty = "n", x = 37, y = 0.6, cex = 1.5))
    title(ylab="Relative contribution to assemblage",cex.lab=1.5)
    # labs = assemblage_labels[assemblage_reordering]
    labs = 1:16
    text(cex=1.5, x=x+0.9, y=1.035, labels = labs, xpd=TRUE, srt=0, pos=2)
  }
  if (taxo.barplot)
  {
    # Reduced to the 11 first groups:
    selected_OTU_proportion = rev(c(1,2,3,7,4,5,6,8,9,10,11))
    x = barplot(rbind(assemblage_group_composition[selected_groups,assemblage_reordering][selected_OTU_proportion,], 
                      Others = rep(1,K)-colSums(assemblage_group_composition[selected_groups,assemblage_reordering][selected_OTU_proportion,])),
                col=colours,
                legend.text = rev(c(rownames(assemblage_group_composition[selected_groups,assemblage_reordering][selected_OTU_proportion,]),"Others")),
                xaxt="n", space = 0.5, 
                # args.legend = list(bty = "n", x= 32.5, y=0.5))
                args.legend = list(bty = "n", x= 31.5, y=0.4, cex = 1.3, fill = colours),
                ylim=c(1,0))
    #colours = terrain.colors(nrow(silva_phylum_mat1))[sample(1:nrow(silva_phylum_mat1),nrow(silva_phylum_mat1))]
    #x = barplot(sweep(silva_phylum_mat1,2,colSums(silva_phylum_mat1),"/"),col=colours,
    #xaxt="n", space = 0.5)
    #legend("topright",rownames(silva_phylum_mat1),col=colours)
    title(ylab="Prevalence in biome",cex.lab=1.3)
    # labs = assemblage_labels[topic_reordering]
    labs = 1:16
    text(cex=1.5, x=x+0.7, y=-0.035, labels = labs, xpd=TRUE, srt=0, pos=2)
  }
  dev.off()
  ##########
  
  ##############
  # pdf(paste0(figure_folder,"/Assemblage_composition_functions_AllTaxa_K",K,"_stackedBarPlot.pdf"), width = 9.5)
  # par(mar=c(1.5,4.1,2.1,11.2))
  # # x = barplot(rbind(assemblage_function_composition[selected_functions,assemblage_reordering], Others = rep(1,K)-colSums(assemblage_function_composition[selected_functions,assemblage_reordering])),
  # #             col=c(rainbow_hcl(length(which(selected_functions)))[sample(1:length(which(selected_functions)),length(which(selected_functions)))],"grey"),
  # #             legend.text = T, xaxt="n", space = 0.5, 
  # #             args.legend = list(bty = "n", x= 25 + 1, y=1))
  # x = barplot(assemblage_function_composition[nrow(assemblage_function_composition):1,assemblage_reordering],
  #             col=rev(function_colors),
  #             legend.text = T, xaxt="n", space = 0.5, cex.axis = 1.3,
  #             args.legend = list(bty = "n", x = 34, y = 0.6, cex = 1.3))
  # # title(ylab="Prevalence in biome",cex.lab=1.3)
  # # labs = assemblage_labels[assemblage_reordering]
  # labs = 1:16
  # text(cex=1.5, x=x+0.8, y=1.035, labels = labs, xpd=TRUE, srt=0, pos=2)
  # dev.off()
  ###############
}

# Fig. A1.3 - bioregions
{
######## Complete:
#col = c(color.pal1(6),color.pal2(3),color.pal3(3),color.pal4(4))[c(4,11,10,7,3,2,6,5,14,8,9,12,13,16,15)]

# pdf(paste0(figure_folder,"/AllTaxa_Gibbs100r",nb_topics,"t_stationClusteringUPGMAHellinger_SURonly_",nb_clust,"clustersPVclust0.96.pdf"))
# # present_topics = apply(spatial_topicmix_kriged_all_topics[,3:(2+(nb_topics0))],1,function(g) which(g>0))
# #
# plot(bat, image = T, land = T, lty = 0, bty = "n", xaxt = "n", yaxt = "n", ann = F, bpal = list(c(0, max(bat), greys[1]), c(min(bat), 0, blues))) #plot map without isobaths
# plot(bat, lwd = 0.8, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline
# points(cbind(spatial_topicmix_kriged_all_topics$x[stations_depths[,2] == "SUR"],spatial_topicmix_kriged_all_topics$y[stations_depths[,2] == "SUR"]),
#        pch = 21, cex=1.2, bg = col[grp_Hellinger])
# dev.off()

##################

bat = getNOAA.bathy(-180, 180, -90, 90, res = 50, keep=T)

data.folder_name = paste0(data_folder,"/18S_V9_TARA_CompleteSizeRange_byStationByDepth_AllTAxa",noArcticNoBiomark_insert,noLagoon_insert)

load(paste0(data.folder_name,"/coord.Rdata"))
stations_depths = as.data.frame(matrix(nrow=nrow(coord),ncol=2,dimnames = list(rownames(coord),c("Station","Depth"))))
for (station_depth_index in 1:nrow(coord))
{
  stations_depths[station_depth_index,1:2] = c(strsplit(rownames(coord)[station_depth_index],split=" ",fixed=T)[[1]][1],
                                               strsplit(rownames(coord)[station_depth_index],split=" ",fixed=T)[[1]][2])
}

pie_positions.SUR = read.table(paste0(results_folder,"/Surf_pies_modified.txt"),sep=",",header=F)
# pie_positions.SUR[c(55,57) - 1,] = pie_positions.SUR[c(57,55) - 1,]
pie_positions.SUR[55:57 - 1,] = pie_positions.SUR[c(56,57,55) - 1,]
pie_positions.SUR[56 - 1,] = pie_positions.SUR[56 - 1,] + c(1.5,0)
pie_positions.SUR[45:46 - 1,] = pie_positions.SUR[46:45 - 1,]
pie_positions.SUR[51 - 1,] = pie_positions.SUR[51 - 1,] + c(-2,3)
rownames(pie_positions.SUR) = rownames(coord[stations_depths[,2] == "SUR",])[-44]

pie_positions.DCM = read.table(paste0(figure_folder,"/DCM_pies_modified.txt"),sep=",",header=F)
rownames(pie_positions.DCM) = rownames(coord[stations_depths[,2] == "DCM",])

# library(gridBase)

nb_topics = optimalK_prevalence.min.crossValid.allTaxa
spatial_topicmix_kriged = readRDS(paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics",nb_topics,"_nb_iter1000_nb_real100_meanPosteriorDistributedLlh_thin25_burnin2000_occurrence/1st_closestToMean_realization/Spatial_topicmix_kriged.rds"))

# selecting the z.pred columns in all topics:
documents = unlist(lapply(spatial_topicmix_kriged,function(g) g$z.pred))
# setting one topic per column
documents = matrix(documents,ncol=nb_topics,dimnames=list(rownames(spatial_topicmix_kriged[[1]]),paste0("assemblage",1:nb_topics)))
# Concatenating the spatial coord. with documents into spatial_topicmix_kriged_all_topics:
spatial_topicmix_kriged_all_topics = cbind(data.frame(x = spatial_topicmix_kriged[[1]]$x, y = spatial_topicmix_kriged[[1]]$y),documents)

Hellinger_stations = 1/sqrt(2)*dist(sqrt(documents), method = "euclidean", diag = FALSE, upper = FALSE)
# Hellinger_stations_SUR = 1/sqrt(2)*dist(sqrt(documents[stations_depths[,2] == "SUR",]), method = "euclidean", diag = FALSE, upper = FALSE)
upgma_Hellinger = agnes(Hellinger_stations, diss =T, method = "average", keep.diss =F, keep.data =F)
library(pvclust)

# Tree cut:
nb_clust = 13
grp_Hellinger = cutree(upgma_Hellinger, k = nb_clust)

# Bootstrap p-values:
# pvclust_Hellinger = pvclust(t(sqrt(documents)),
#                             method.dist = "euclidean",
#                             #method.dist = "cor",
#                             # method.dist = function(x) 1/sqrt(2)*dist(sqrt(x), method = "euclidean", diag = FALSE, upper = FALSE),
#                             # method.dist = "cor", 
#                             method.hclust = "average", 
#                             nboot=1000)
# pvpick_Hellinger = pvpick(pvclust_Hellinger,alpha = 0.96)
# nb_clust = length(pvpick_Hellinger$clusters)
# grp_Hellinger = rep(1,nrow(documents))
# names(grp_Hellinger) = rownames(documents)
# for (i_clust in 1:nb_clust)
# {
#   grp_Hellinger[names(grp_Hellinger) %in% pvpick_Hellinger$clusters[[i_clust]]] = i_clust+1
# }

color.pal1 = colorRampPalette(c("#7F0000","red","darkorange1","darkgoldenrod1","yellow"),space = "Lab")
color.pal2 = colorRampPalette(c("#7FFF7F","darkgreen"),space = "Lab")
color.pal3 = colorRampPalette(c("azure4","antiquewhite3","aliceblue"),space = "Lab")
color.pal4 = colorRampPalette(c("#00007F","#007FFF","cyan"),space = "Lab")
if (nb_clust == 13)
{
  col = c(color.pal1(6),color.pal2(3),color.pal3(3),color.pal4(4))[c(4,2,11,7,6,3,5,14,12,8,9,13,16)]
} else if (nb_clust == 14)
{
  col = c(color.pal1(6),color.pal2(3),color.pal3(3),color.pal4(4))[c(4,2,11,7,6,3,5,14,15,8,9,12,13,16)]
} else if (nb_clust == 15)
{
  col = c(color.pal1(6),color.pal2(3),color.pal3(3),color.pal4(4))[c(4,2,11,10,7,6,3,5,14,15,8,9,12,13,16)]
} else if (nb_clust == 16)
{
  col = c(color.pal1(6),color.pal2(3),color.pal3(3),color.pal4(4))[c(4,2,11,10,7,6,3,5,1,14,15,8,9,12,13,16)]
# col = c(color.pal1(6),color.pal2(3),color.pal3(3),color.pal4(7))
} else if (nb_clust == 17)
{
  col = c(color.pal1(6),color.pal2(3),color.pal3(3),color.pal4(5))[c(4,2,11,10,7,6,3,5,1,14,15,8,9,12,13,16,17)]
} else if (nb_clust == 18)
  col = c(color.pal1(6),color.pal2(3),color.pal3(3),color.pal4(5),"violet")[c(4,2,11,10,7,6,18,3,5,1,14,15,8,9,12,13,16,17)]
# if (nb_clust == 14) #cut tree
# {
#   col = c(color.pal1(6),color.pal2(3),color.pal3(3),color.pal4(4))[c(4,7,11,10,6,3,5,14,15,8,9,12,13,16)]
# } else if (nb_clust == 12) # pvclust
#   col = c("pink",color.pal1(6),color.pal2(3),color.pal3(3),color.pal4(4))[c(1,c(3,10,1,5,2,4,6,12,14,7,8,16)+1)]

# plot.map.group = list()
pdf(paste0(figure_folder,"/FigA1.3_bioregion.map.pdf"),
    width=7.5*2.2/1.2,height=12/3*4/1.2/2,onefile=T)

blues = c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys = c(grey(0.6), grey(0.93), grey(0.99))
# Plotting stations associated with each dominant assemblage  
# SUR:
# plot(bat, image = T, land = T, lty = 0, bty = "n", xaxt = "n", yaxt = "n", ann = F, bpal = list(c(0, max(bat), greys[1]), c(min(bat), 0, blues))) #plot map without isobaths
plot(bat, image = T, land = T, lty = 0, bty = "n", xaxt = "n", yaxt = "n", ann = F, bpal = list(c(0, max(bat), greys[1]), c(min(bat), 0, blues[3])))
plot(bat, lwd = 0.2, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline
pies.coord = pie_positions.SUR[rownames(pie_positions.SUR) %in% rownames(documents),]
stations.coord = coord[rownames(coord) %in% rownames(pie_positions.SUR) & rownames(coord) %in% rownames(documents),]
segments(x0=stations.coord[,2],y0=stations.coord[,1],x1=pies.coord[,1],y1=pies.coord[,2])
points(pie_positions.SUR[rownames(pie_positions.SUR) %in% rownames(documents),],
       pch = 21, cex=4.4, bg = col[grp_Hellinger[stations_depths[,2] == "SUR"]])
# DCM
plot(bat, image = T, land = T, lty = 0, bty = "n", xaxt = "n", yaxt = "n", ann = F, bpal = list(c(0, max(bat), greys[1]), c(min(bat), 0, blues[3])))
plot(bat, lwd = 0.2, deep = 0, shallow = 0, step = 0, add = TRUE) 
pies.coord = pie_positions.DCM[rownames(pie_positions.DCM) %in% rownames(documents),]
stations.coord = coord[rownames(coord) %in% rownames(pie_positions.DCM) & rownames(coord) %in% rownames(documents),]
segments(x0=stations.coord[,2],y0=stations.coord[,1],x1=pies.coord[,1],y1=pies.coord[,2])
points(pie_positions.DCM[rownames(pie_positions.DCM) %in% rownames(documents),],
       pch = 21, cex=4.4, bg = col[grp_Hellinger[stations_depths[,2] == "DCM"]])
# points(cbind(spatial_topicmix_kriged_all_topics$x[stations_depths[,2] == "SUR"],spatial_topicmix_kriged_all_topics$y[stations_depths[,2] == "SUR"]),
#        pch = 21, cex=1.2, bg = col[grp_Hellinger])
dev.off()

}

# First PCoA axis functional boxplot:
{
  pdf(paste0(figure_folder,"/PCoA.1st.axis.functional.boxplot_",div_threshold,"plusOTUs_selected100+1.pdf"),
      width=7,height=7)
  boxplot.1st.axis.functions = ggplot(data = data.frame(x = dominant_function0[!is.na(dominant_function0) & !dominant_function0 %in% point_groups],
                                                         y = NormalizedVI_pcoa[[3]][,1]
                                                         [!is.na(dominant_function0) & !dominant_function0 %in% point_groups])) +
    scale_x_discrete(limits=levels(dominant_function0)) +
    geom_boxplot(aes(x,y)) +
    # scale_y_log10() +
    geom_point(data = data.frame(x = factor(point_groups),
                                 y = NormalizedVI_pcoa[[3]][,1]
                                 [dominant_function0 %in% point_groups]),
               aes(x,y)) +
    theme_bw() +
    # ggtitle(LETTERS[8]) +
    # scale_y_log10() +
    # geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.title=element_text(size=22),
          axis.text=element_text(size=22),
          plot.title=element_text(hjust=0, size=24),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin=unit(c(0.5,1,-5,5),"mm")) +
    labs(x="", y="First PCoA axis")
  grid::grid.draw(boxplot.1st.axis.functions)
  dev.off()
  
  x = log10(diversity)[selected_groups & !taxo_groups %in% c("Dinophyceae","Collodaria") & diversity>div_threshold]
  # x2 = log10(size_relativeAbund)[selected_groups & !taxo_groups %in% c("Dinophyceae","Collodaria","apusozoan_*") & diversity>div_threshold]
  # g1 = factor(dominant_function1,c("Phagotrophs","Metazoans","Parasites","Phototrophs"))
  # g = factor(dominant_function1,c("Phototrophs","Phagotrophs","Metazoans","Parasites"))
  g = factor(dominant_function1,c("Phagotrophs","Phototrophs","Metazoans","Parasites"))
  # g = dominant_function1
  # contrasts(g) = contr.poly(4)
  Phago.vs.rest = c(-3,1,1,1)
  Photo.vs.rest = c(1,-3,1,1)
  Meta.vs.rest = c(1,1,-3,1)
  Para.vs.rest = c(1,1,1,-3)
  # To switch between default contrasts (i.e., no contrast, or "treatment contrasts") and Hebert contrasts:
  # options(contrasts=c("contr.helmert", "contr.poly"))
  # options(contrasts=c("contr.treatment", "contr.poly"))
  # contrasts(g) = cbind(Para.vs.rest,Photo.vs.rest,Meta.vs.rest,Phago.vs.rest)
  # contrasts(g) = cbind(Para.vs.rest,Phago.vs.rest)
  # contrasts(g) = cbind(Meta.vs.rest,Photo.vs.rest,Phago.vs.rest)
  # contrasts(g) = cbind(Photo.vs.rest,Meta.vs.rest)
  Meta.vs.Photo = c(0,1,-1,0)
  Phago.vs.Meta = c(1,0,-1,0)
  Phago.vs.Photo = c(1,-1,0,0)
  Para.vs.Meta = c(0,0,1,-1)
  Para.vs.Photo = c(0,1,0,-1)
  Para.vs.Phago = c(1,0,0,-1)
  contrasts(g) = cbind(Meta.vs.Photo,Phago.vs.Photo,Phago.vs.Meta)
  # contrasts(g) = cbind(Para.vs.Meta,Para.vs.Photo,Para.vs.Phago)
  # contrasts(g) = cbind(Meta.vs.Photo,Phago.vs.Photo,Para.vs.Photo)
  y = NormalizedVI_pcoa[[3]][!names(NormalizedVI_pcoa[[3]][,1]) %in% c("Dinophyceae","Collodaria"),1]
  ####################
  # aov (type I, sensitive to ordering) test:
  ancova = aov(y ~ g)
  # CList = list("Para.vs.rest" = 1,
  #              "Photo.vs.rest" = 2,
  #              "Meta.vs.rest" = 3,
  #              "Phago.vs.rest" = 4)
  # CList = list("Para.vs.rest" = 1,
  #              "Phago.vs.rest" = 2)
  # CList = list("Meta.vs.rest" = 1,
  #              "Photo.vs.rest" = 2,
  #              "Phago.vs.rest" = 3)
  CList = list("Meta.vs.Photo" = 1,
               "Phago.vs.Photo" = 2,
               "Phago.vs.Meta" = 3)
  # CList = list("Para.vs.Meta" = 1,
  #              "Para.vs.Photo" = 2,
  #              "Para.vs.Phago" = 3)
  # CList = list("Meta.vs.Photo" = 1,
  #              "Phago.vs.Photo" = 2,
  #              "Para.vs.Photo" = 3)
  summary(ancova,
          split=list(g=CList))
  #####################
  # lm (type II, not sensitive to ordering) test:
  ancova = lm(y ~ g + x)
  summary(ancova)
  ####################
  x = log10(as.vector(diversity))[selected_groups & diversity<2000]
  x2 = log10(size_relativeAbund)[selected_groups & diversity<2000]
  y = NormalizedVI_pcoa[[3]][diversity[selected_groups]<2000,1]
  model = lm(y ~ x + x2)
  varpart(y,x,x2)
  library(ppcor)
  pcor.test(y,x,x2)
  ####################
  y = NormalizedVI_pcoa[[3]][,2]
  x = log10(as.vector(diversity))[selected_groups]
  x2 = log10(size_relativeAbund)[selected_groups]
  pcor.test(y,x,x2)
  ####################
  
  t.test(NormalizedVI_pcoa[[3]][,1][dominant_function0 %in% c("Bacillariophyta","Other phototrophs")],
         NormalizedVI_pcoa[[3]][,1][dominant_function0 %in% c("Gel. carn. filterers","Other metazoa","Pteropoda","Copepoda")])
  t.test(NormalizedVI_pcoa[[3]][,1][dominant_function0 %in% c("Bacillariophyta","Other phototrophs")],
         NormalizedVI_pcoa[[3]][,1][dominant_function0 %in% "Phagotrophs"])
  t.test(NormalizedVI_pcoa[[3]][,1][dominant_function0 %in% c("Bacillariophyta","Other phototrophs")],
         NormalizedVI_pcoa[[3]][,1][dominant_function0 %in% "Parasites"])
  t.test(NormalizedVI_pcoa[[3]][,1][dominant_function0 %in% "Phagotrophs"],
         NormalizedVI_pcoa[[3]][,1][dominant_function0 %in% c("Gel. carn. filterers","Other metazoa","Pteropoda","Copepoda")])
  t.test(NormalizedVI_pcoa[[3]][,1][dominant_function0 %in% "Parasites"],
         NormalizedVI_pcoa[[3]][,1][dominant_function0 %in% c("Gel. carn. filterers","Other metazoa","Pteropoda","Copepoda")])
  t.test(NormalizedVI_pcoa[[3]][,1][dominant_function0 %in% "Parasites"],
         NormalizedVI_pcoa[[3]][,1][dominant_function0 %in% "Phagotrophs"])
}

# Fig. pres. summary sketch:
{
  pdf(paste0(figure_folder,"/Fig5_conclusion.sketch.pdf"),
      width=7,height=7)
  plot.VI.PCoA.biplot = ggplot(data = data.frame(x = NormalizedVI_pcoa[[3]][,1], y = NormalizedVI_pcoa[[3]][,2])) + 
    geom_point(aes(x,y), size = 3, na.rm=T) +
    theme_bw() +
    # ggtitle(LETTERS[1]) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    # ggtitle("Norm. V.I. between clades' biogeographies") +
    theme(axis.title=element_text(size=22),
          axis.text=element_text(size=22),
          text = element_text(size=22),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,0.5,0.5,0.5),"mm")) +
          # legend.position="bottom",
          # legend.text=element_text(size=22),
          # legend.title=element_text(size=16, hjust = 10)) +
    labs(x="PCoA axis 1", y="PCoA axis 2")
    # scale_colour_manual(values = setNames(ten_colors,functions0), 
    #                     na.value = "grey50", na.translate = FALSE, guide = "colourbar", name = "Functional group") +
    # geom_smooth(method='lm') +
    # guides(colour = guide_legend(title=NULL, nrow = 5, ncol = 2))
  grid::grid.draw(plot.VI.PCoA.biplot)
  dev.off()
}

# Fig. pres. axis 1 vs. diversity:
{
  pdf(paste0(figure_folder,"/Fig_axis1_div_pres.pdf"),
      width=7,height=7)
  x = as.vector(diversity)[selected_groups & diversity>div_threshold][-69]
  y = NormalizedVI_pcoa[[3]][,1]
  # lin.regr = lm(y ~ log10(x))
  # b = summary(lin.regr)$coefficients[1,1]
  # a = summary(lin.regr)$coefficients[2,1]
  plot.VI.PCoA.axis1.vs.diversity = ggplot(data = data.frame(x = x, y = y)) +
    # geom_point(aes(x,y1)) +
    geom_point(aes(x,y)) +
    scale_x_log10() +
    # scale_y_log10() +
    # geom_abline(intercept = b, slope = a) +
    theme_bw() +
    # ggtitle(LETTERS[6]) +
    # geom_hline(yintercept = 1, linetype="dashed") +
    # ggtitle("RDA adjR2 3 vs 2 abiotic PCA axes with stdzation") +
    # ggtitle("Pure abiotic vs. total abiotic adjR2 3 axes") +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          # axis.title.x=element_text(vjust = 45),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(0.5,10,0.5,0.5),"mm")) +
    #labs(x="Mean size (micron)", y=paste("Ratio of",c("surface","DCM")[i_case],"variance\n purely explained by currents vs. envir.")) +
    labs(x="Clade diversity", y="PCoA axis 1") +
    geom_smooth(data = data.frame(x = x[x<2000], y = y[x<2000]), aes(x,y),method='lm',col="black",se=T) +
    geom_segment(aes(x = 2000, xend = max(x), y = 0.167, yend = 0.167), linetype="dashed")
  # geom_smooth(aes(x,y),method='lm')
  grid::grid.draw(plot.VI.PCoA.axis1.vs.diversity)
  dev.off()
}

# Abiotic - biotic - MEM maps:
{
  data.folder_name = paste0(data_folder,"/18S_V9_TARA_CompleteSizeRange_byStationByDepth_AllTAxa",noArcticNoBiomark_insert,noLagoon_insert)
  load(paste0(data.folder_name,"/coord.Rdata"))
  load(paste0(data.folder_name,"/sample_ref.Rdata"))
  
  # Loading the WOA13 averaged-interpolated abiotic data (data sent by Paul Frémont and Olivier Jaillon):
  station_ref_all = read.table(paste0(data_folder,"/Abiotic_data/woa13_env_tara_all.csv"),sep=";",header=T,row.names=1)
  # Reordering the variables as in watteaux_selected_variables
  station_ref_all = station_ref_all[,c("T","no3","O2","o2s","aou","po4","si")]
  
  # selected_stations contains the stations effectively used for the decomposition, in the same order as in spatial_topicmix_kriged, with the names used in Iudicone's travel_time_matrix
  # i.e., converts rownames(coord) into the names in travel_time_matrix
  selected_stations.l = lapply(strsplit(rownames(coord),split=" ",fixed=T),strsplit,split="_",fixed=T)
  selected_stations = vector(length = nrow(coord), mode = "character")
  for (i in 1:nrow(coord))
  {
    selected_stations[i] = paste(selected_stations.l[[i]][[1]][2],selected_stations.l[[i]][[2]][1],sep = "_")
    station_split = strsplit(selected_stations[i],split="",fixed=T)[[1]]
    if (station_split[1] == "0")
    {
      j = 1
      while (station_split[j] == "0")
        j = j+1
      selected_stations[i] = paste(station_split[-(1:(j-1))],collapse = "")
    }
  }
  
  stations_depths = as.data.frame(matrix(nrow=nrow(coord),ncol=4))
  colnames(stations_depths) = c("Station","Depth","Station_depth","Station_number")
  for (station_depth_index in 1:nrow(coord))
  {
    stations_depths[station_depth_index,1:2] = c(strsplit(rownames(coord)[station_depth_index],split=" ",fixed=T)[[1]][1],
                                                 strsplit(rownames(coord)[station_depth_index],split=" ",fixed=T)[[1]][2])
    station_number = strsplit(stations_depths[station_depth_index,1],split="_",fixed=T)[[1]][2]
    # Removing the "0"s in front of the station number to conform to the rownames in sample_ref_indices  
    while (substr(station_number,1,1) == "0")
      station_number = substr(station_number,2,nchar(station_number))
    stations_depths[station_depth_index,3] = paste(c(station_number,stations_depths[station_depth_index,2]),collapse = "_")
    stations_depths[station_depth_index,4] = station_number
  }
  stations_names = levels(as.factor(stations_depths$Station))
  
  library(marmap)
  bat = getNOAA.bathy(-180, 180, -90, 90, res = 50, keep=T)
  blues = c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
  greys = c(grey(0.6), grey(0.93), grey(0.99))
  
  pie_positions.SUR = read.table(paste0(results_folder,"/Surf_pies_modified.txt"),sep=",",header=F)
  # pie_positions.SUR[c(55,57) - 1,] = pie_positions.SUR[c(57,55) - 1,]
  pie_positions.SUR[55:57 - 1,] = pie_positions.SUR[c(56,57,55) - 1,]
  pie_positions.SUR[56 - 1,] = pie_positions.SUR[56 - 1,] + c(1.5,0)
  pie_positions.SUR[45:46 - 1,] = pie_positions.SUR[46:45 - 1,]
  pie_positions.SUR[51 - 1,] = pie_positions.SUR[51 - 1,] + c(-2,3)
  rownames(pie_positions.SUR) = rownames(coord[stations_depths[,2] == "SUR",])[-44]
  
  # abiotic
  ###########
  federico_selected_variables = c("Iron.5m","Depth.nominal")
  # federico_selected_variables = "Iron.5m"
  # watteaux_selected_variables = c("Temperature","Latitude_woa13","Longitude_woa13","Nitrate","OxygenDissolved","OxygenSaturation","OxygenUtilization","Phosphate","Silicate","Distance_to_coast")
  watteaux_selected_variables = c("Temperature","Nitrate","OxygenDissolved","OxygenSaturation","OxygenUtilization","Phosphate","Silicate")
  # watteaux_selected_variables = c("Temperature","OxygenDissolved","OxygenSaturation","OxygenUtilization","Phosphate","Silicate","Distance_to_coast")
  colnames(station_ref_all) = watteaux_selected_variables
  
  abiotic_data_all = matrix(nrow = nrow(coord), ncol = length(federico_selected_variables),
                        dimnames = list(rownames(coord), federico_selected_variables), data = 0)
  for (i in 1:nrow(coord))
  {
    sample_ref_indices = sample_ref$Station.label == stations_depths[i,1] & sample_ref$Depth == stations_depths[i,2]
    if (length(which(sample_ref_indices)) > 1)
    {
      abiotic_data_all[i,] = colMeans(apply(sample_ref[sample_ref_indices,federico_selected_variables],2,as.numeric),na.rm=T)
    } else
      abiotic_data_all[i,] = as.numeric(as.vector(sample_ref[sample_ref_indices,federico_selected_variables]))
  }
  # colnames(station_ref_Arctic) = watteaux_selected_variables
  # abiotic_data[,(length(federico_selected_variables)+1):(length(federico_selected_variables)+length(watteaux_selected_variables))] = rbind(station_ref[stations_depths[1:(which(stations_depths[,3] == "155_SUR")-1),3],watteaux_selected_variables],
  #                                                                                                                                          station_ref_Arctic[stations_depths[which(stations_depths[,4] == "155")[1]:nrow(stations_depths),4],])
  abiotic_data_all = cbind(abiotic_data_all,station_ref_all[stations_depths[,3],watteaux_selected_variables])
  abiotic_data_all = abiotic_data_all[,colnames(abiotic_data_all) != "Depth.nominal"]
   
  library(ade4)
  
  abiotic_PCA = list()
  for (case in c("SUR","DCM"))
  {
    i_case = which(c("SUR","DCM") == case)
    abiotic_data = abiotic_data_all[stations_depths[,2] == case,]
    
    if (any(apply(abiotic_data,2,is.na)))
    {
      abiotic_data_trans = abiotic_data[,!(as.vector(unlist(lapply(apply(apply(abiotic_data,2,is.na),2,which),length))) > 3)]
    } else 
      abiotic_data_trans = abiotic_data
    
    # Centering and standardizing the abiotic data:
    # abiotic_data_trans = sweep(sweep(abiotic_data_trans,2,colMeans(abiotic_data_trans, na.rm = T),"-"),2,apply(abiotic_data_trans,2,sd,na.rm=T),"/")
    abiotic_data_trans = scale(abiotic_data_trans)
    
    # PCA on standardized ("normed") data, ie centered and divided by std. dev. for each column 
    # Removing stations where one abiotic value is missing:
    stations_to_remove = apply(apply(abiotic_data_trans,2,is.na),1,any)
    # Removing Latitude from the abiotic variables on which PCA is applied
    abiotic_PCA[[i_case]] = dudi.pca(as.data.frame(abiotic_data_trans[!stations_to_remove,colnames(abiotic_data_trans)!="Distance_to_coast"]), row.w = rep(1, nrow(abiotic_data_trans[!stations_to_remove,]))/nrow(abiotic_data_trans[!stations_to_remove,]), 
                           # col.w = rep(1, ncol(abiotic_data_trans) - 1),
                           # col.w = c(1,1/3,1/3,1/3,1,1/3,1/3,1/3,1,1),
                           col.w = rep(1,ncol(abiotic_data_trans)),
                           # col.w = c(1,1,1,1/3,1/3,1/3,1,1),
                           center = TRUE, scale = TRUE, scannf = F, nf = ncol(abiotic_data_trans))
  }
  
  pdf(paste0(figure_folder,"/FigS_3.surf.abiotic.pca.axes_marmap.pdf"),
      width=7.5*2.2/1.2,height=12/3*4/1.2/2*3,onefile=T)
  par(mfrow = c(3,1))
  for (k in 1:3)
  {
    # Plotting Surface:
    plot(bat, image = T, land = T, lty = 0, bty = "n", xaxt = "n", yaxt = "n", ann = F, bpal = list(c(0, max(bat), greys[1]), c(min(bat), 0, blues[3]))) #plot map without isobaths
    plot(bat, lwd = 0.2, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline
    z = as.data.frame(abiotic_PCA[[1]]$li)[rownames(abiotic_PCA[[1]]$li) %in% rownames(pie_positions.SUR),k]
    col = rgb(colorRamp(c("blue","grey","red"),space = "Lab")(rescale(z,to=c(0,1),from=c(-max(range(abs(z))),max(range(abs(z)))))),maxColorValue = 255)
    pies.coord = pie_positions.SUR[rownames(pie_positions.SUR) %in% rownames(abiotic_PCA[[1]]$li),]
    stations.coord = coord[rownames(coord) %in% rownames(pie_positions.SUR) & rownames(coord) %in% rownames(abiotic_PCA[[1]]$li),]
    segments(x0=stations.coord[,2],y0=stations.coord[,1],x1=pies.coord[,1],y1=pies.coord[,2])
    points(pies.coord,
           pch = 21, cex=7, bg = col)
    # title(paste("Abiotic PCA axis",k), cex.main = 2)
  }
  dev.off()
  
  temperature = list()
  for (case in c("SUR","DCM"))
  {
    i_case = which(c("SUR","DCM") == case)
    abiotic_data = abiotic_data_all[stations_depths[,2] == case,]
    
    if (any(apply(abiotic_data,2,is.na)))
    {
      abiotic_data_trans = abiotic_data[,!(as.vector(unlist(lapply(apply(apply(abiotic_data,2,is.na),2,which),length))) > 3)]
    } else 
      abiotic_data_trans = abiotic_data
    
    abiotic_data_trans = scale(abiotic_data_trans)
    
    # PCA on standardized ("normed") data, ie centered and divided by std. dev. for each column 
    # Removing stations where one abiotic value is missing:
    stations_to_remove = apply(apply(abiotic_data_trans,2,is.na),1,any)
    temperature[[i_case]] = abiotic_data_trans[!stations_to_remove,colnames(abiotic_data_trans)=="Temperature"]
    names(temperature[[i_case]]) = rownames(abiotic_data_trans[!stations_to_remove,])
  }
  
  pdf(paste0(figure_folder,"/FigS_temperature_marmap.pdf"),
      width=7.5*2.2/1.2,height=12/3*4/1.2/2*3,onefile=T)
  for (i_case in 1)
  {
  # Plotting Surface:
  plot(bat, image = T, land = T, lty = 0, bty = "n", xaxt = "n", yaxt = "n", ann = F, bpal = list(c(0, max(bat), greys[1]), c(min(bat), 0, blues[3]))) #plot map without isobaths
  plot(bat, lwd = 0.2, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline
  z = temperature[[i_case]][names(temperature[[i_case]]) %in% rownames(pie_positions.SUR)]
  col = rgb(colorRamp(c("blue","grey","red"),space = "Lab")(rescale(z,to=c(0,1),from=c(-max(range(abs(z))),max(range(abs(z)))))),maxColorValue = 255)
  pies.coord = pie_positions.SUR[rownames(pie_positions.SUR) %in% names(temperature[[i_case]]),]
  stations.coord = coord[rownames(coord) %in% rownames(pie_positions.SUR) & rownames(coord) %in% names(temperature[[i_case]]),]
  segments(x0=stations.coord[,2],y0=stations.coord[,1],x1=pies.coord[,1],y1=pies.coord[,2])
  points(pies.coord,
         pch = 21, cex=5, bg = col)
  # title(paste("Abiotic PCA axis",k), cex.main = 2)
  }
  dev.off()
  
  # biotic
  ##########
  relativeAbund_file = paste0(results_folder,"/groups_relativeAbund",noArcticNoBiomark_insert,noLagoon_insert,".rds")
  relativeAbund_all = readRDS(relativeAbund_file)
  
  relativeAbund_PCA0 = list()
  for (case in c("SUR","DCM"))
  {
    i_case = which(c("SUR","DCM") == case)
    relativeAbund = relativeAbund_all[stations_depths[,2] == case,]
    
    # PCA on standardized ("normed") data, ie centered and divided by std. dev. for each column 
    stations_to_remove = apply(apply(relativeAbund,2,is.na),1,any)
    relativeAbund_PCA0[[i_case]] = dudi.pca(as.data.frame(relativeAbund[!stations_to_remove,selected_groups]), row.w = rep(1, nrow(relativeAbund[!stations_to_remove,]))/nrow(relativeAbund[!stations_to_remove,]), 
                                  col.w = rep(1, ncol(relativeAbund[,selected_groups])),
                                  center = T, scale = T, scannf = F, nf = ncol(relativeAbund[,selected_groups]))
  }
  
  pdf(paste0(figure_folder,"/FigS_8.surf.biotic.pca.axes_marmap.pdf"),
      width=7.5*2.2/1.2*2,height=12/3*4/1.2/2*4,onefile=T)
  layout(rbind(c(1,2),c(3,4),c(5,6),c(7,8)))
  par(mar = 5*c(-1,-1,-1,-1))
  for (k in 1:8)
  {
    # Plotting Surface:
    plot(bat, image = T, land = T, lty = 0, bty = "n", xaxt = "n", yaxt = "n", ann = F, bpal = list(c(0, max(bat), greys[1]), c(min(bat), 0, blues[3]))) #plot map without isobaths
    plot(bat, lwd = 0.2, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline
    z = as.data.frame(relativeAbund_PCA0[[1]]$li)[rownames(relativeAbund_PCA0[[1]]$li) %in% rownames(pie_positions.SUR),k]
    col = rgb(colorRamp(c("blue","grey","red"),space = "Lab")(rescale(z,to=c(0,1),from=c(-max(range(abs(z))),max(range(abs(z)))))),maxColorValue = 255)
    pies.coord = pie_positions.SUR[rownames(pie_positions.SUR) %in% rownames(relativeAbund_PCA0[[1]]$li),]
    stations.coord = coord[rownames(coord) %in% rownames(pie_positions.SUR) & rownames(coord) %in% rownames(relativeAbund_PCA0[[1]]$li),]
    segments(x0=stations.coord[,2],y0=stations.coord[,1],x1=pies.coord[,1],y1=pies.coord[,2])
    points(pies.coord,
           pch = 21, cex=7, bg = col)
    # title(paste("Abiotic PCA axis",k), cex.main = 2)
  }
  dev.off()
  
  # MEM
  ##########
  # Need to compute the MEMs first (code from Multigroup_Tara_analyses.R)
  
  travel.folder_name = paste0(data_folder,"/Abiotic_data")
  travel_time_matrix_SUR = read.table(paste0(travel.folder_name,"/minAijji_tarrive_min_surface_10.csv"),sep="\t",header=T,row.names=1)
  travel_time_matrix_SUR = travel_time_matrix_SUR[,-ncol(travel_time_matrix_SUR)]
  travel_time_matrix_DCM = read.table(paste0(travel.folder_name,"/tarrive_min_75m_10.csv"),sep=";",header=T,row.names=1)
  for (i in 2:nrow(travel_time_matrix_DCM))
  {
    for (j in 1:(i-1))
    {
      travel_time_matrix_DCM[i,j] = travel_time_matrix_DCM[j,i] = min(travel_time_matrix_DCM[i,j],travel_time_matrix_DCM[j,i])
    }
  }
  
  selected_travel_time_matrix_SUR = travel_time_matrix_SUR[selected_stations,paste0("X",selected_stations)]
  selected_travel_time_matrix_SUR = selected_travel_time_matrix_SUR[stations_depths[,2] == "SUR",stations_depths[,2] == "SUR"]
  selected_travel_time_matrix_DCM = travel_time_matrix_DCM[selected_stations,paste0("X",selected_stations)]
  selected_travel_time_matrix_DCM = selected_travel_time_matrix_DCM[stations_depths[,2] == "DCM",stations_depths[,2] == "DCM"]
  
  coord_SUR = coord[stations_depths[,2] == "SUR",]
  coord_DCM = coord[stations_depths[,2] == "DCM",]
  
  #######################################
  # Building t_min-based MEMs
  ###########################
  tmin_thres_SUR = 2.1
  tmin_thres_DCM = 3.15
  
  truncated_Tmin_SUR = as.matrix(selected_travel_time_matrix_SUR)
  truncated_Tmin_SUR[is.nan(truncated_Tmin_SUR)] = 4*tmin_thres_SUR
  truncated_Tmin_SUR[truncated_Tmin_SUR > tmin_thres_SUR] = 4*tmin_thres_SUR
  diag(truncated_Tmin_SUR) = 4*tmin_thres_SUR
  
  truncated_Tmin_DCM = as.matrix(selected_travel_time_matrix_DCM)
  truncated_Tmin_DCM[is.nan(truncated_Tmin_DCM)] = 4*tmin_thres_DCM
  truncated_Tmin_DCM[truncated_Tmin_DCM > tmin_thres_DCM] = 4*tmin_thres_DCM
  diag(truncated_Tmin_DCM) = 4*tmin_thres_DCM
  
  devtools::source_url("https://github.com/guilhemSK/Useful_functions/raw/main/pcoa.all_fun.R")
  truncated_Tmin_SUR_pcoa = pcoa.all(as.dist(truncated_Tmin_SUR), rn = rownames(coord_SUR))
  truncated_Tmin_DCM_pcoa = pcoa.all(as.dist(truncated_Tmin_DCM), rn = rownames(coord_DCM))
  
  SUR_MEM = matrix(nrow = nrow(coord), ncol = ncol(truncated_Tmin_SUR_pcoa$vectors), dimnames = list(rownames(coord),1:ncol(truncated_Tmin_SUR_pcoa$vectors)), data = 0)
  SUR_MEM[rownames(truncated_Tmin_SUR_pcoa$vectors),] = truncated_Tmin_SUR_pcoa$vectors
  SUR_MEM_restricted = truncated_Tmin_SUR_pcoa$vectors
  colnames(SUR_MEM_restricted) = paste("SUR.MEM",1:ncol(SUR_MEM_restricted))
  
  DCM_MEM = matrix(nrow = nrow(coord), ncol = ncol(truncated_Tmin_DCM_pcoa$vectors), dimnames = list(rownames(coord),1:ncol(truncated_Tmin_DCM_pcoa$vectors)), data = 0)
  DCM_MEM[rownames(truncated_Tmin_DCM_pcoa$vectors),] = truncated_Tmin_DCM_pcoa$vectors
  DCM_MEM_restricted = truncated_Tmin_DCM_pcoa$vectors
  colnames(DCM_MEM_restricted) = paste("DCM.MEM",1:ncol(DCM_MEM_restricted))
  
  ##############
  library(scales)
  
  charac_dist.SUR_scale = readRDS(paste0(data_folder,"/MEM.charac.dist.scale_SUR.rds"))
  charac_tmin.SUR_scale = readRDS(paste0(data_folder,"/MEM.charac.tmin.scale_SUR.rds"))
  
  pdf(paste0(figure_folder,"/FigS_14first.tmin.surf.MEM_marmap.pdf"),
      width=7.5*2.2/1.2*2,height=12/3*4/1.2/2*4,onefile=T)
  layout(rbind(c(1,2),c(3,4),c(5,6),c(7,8)))
  par(mar = 5*c(-1,-1,-1,-1))
  for (k in 1:14)
  {
    # Plotting Surface:
    plot(bat, image = T, land = T, lty = 0, bty = "n", xaxt = "n", yaxt = "n", ann = F, bpal = list(c(0, max(bat), greys[1]), c(min(bat), 0, blues[3]))) #plot map without isobaths
    plot(bat, lwd = 0.2, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline
    z = as.data.frame(SUR_MEM_restricted)[rownames(SUR_MEM_restricted) %in% rownames(pie_positions.SUR),k]
    col = rgb(colorRamp(c("blue","grey","red"),space = "Lab")(rescale(z,to=c(0,1),from=c(-max(range(abs(z))),max(range(abs(z)))))),maxColorValue = 255)
    pies.coord = pie_positions.SUR[rownames(pie_positions.SUR) %in% rownames(SUR_MEM_restricted),]
    stations.coord = coord[rownames(coord) %in% rownames(pie_positions.SUR) & rownames(coord) %in% rownames(SUR_MEM_restricted),]
    segments(x0=stations.coord[,2],y0=stations.coord[,1],x1=pies.coord[,1],y1=pies.coord[,2])
    points(pies.coord,
           pch = 21, cex=7, bg = col)
    title(paste("MEM",k), cex.main = 4)
  }
  dev.off()

  # Need to compute median.adjR2.per.MEM first (see below # Fig. S. adj. R2 hist per MEM)
  
  # pdf(paste0(figure_folder,"/FigS_8.most.selected.surf.MEM_indiv.axes.HB.correction_marmap.pdf"),
  pdf(paste0(figure_folder,"/FigS_8.most.often.selected.surf.MEM_indiv.axes.HB.correction_marmap.pdf"),
      width=7.5*2.2/1.2*2,height=12/3*4/1.2/2*4,onefile=T)
  layout(rbind(c(1,2),c(3,4),c(5,6),c(7,8)))
  par(mar = c(3,1,4,1))
  kk=0
  # for (k in c(5,6,1,3,2,10,8,7))
  # for (k in sort.int(median.adjR2.per.MEM,decreasing=T,index.return = T)$ix)
  for (k in sort.int(nb.selecting.groups.per.MEM,decreasing=T,index.return = T)$ix)
  {
    # kk=kk+1
    # Plotting Surface:
    plot(bat, image = T, land = T, lty = 0, bty = "n", xaxt = "n", yaxt = "n", ann = F, bpal = list(c(0, max(bat), greys[1]), c(min(bat), 0, blues[3]))) #plot map without isobaths
    plot(bat, lwd = 0.2, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline
    z = as.data.frame(SUR_MEM_restricted)[rownames(SUR_MEM_restricted) %in% rownames(pie_positions.SUR),k]
    col = rgb(colorRamp(c("blue","grey","red"),space = "Lab")(rescale(z,to=c(0,1),from=c(-max(range(abs(z))),max(range(abs(z)))))),maxColorValue = 255)
    pies.coord = pie_positions.SUR[rownames(pie_positions.SUR) %in% rownames(SUR_MEM_restricted),]
    stations.coord = coord[rownames(coord) %in% rownames(pie_positions.SUR) & rownames(coord) %in% rownames(SUR_MEM_restricted),]
    segments(x0=stations.coord[,2],y0=stations.coord[,1],x1=pies.coord[,1],y1=pies.coord[,2])
    points(pies.coord,
           pch = 21, cex=7, bg = col)
    title(paste("MEM",k,"- Charac. scale =",
                # signif(charac_tmin.SUR_scale[k], digits = 3),"y"), cex.main = 4)
                signif(charac_dist.SUR_scale[k], digits = 3),"km"), cex.main = 4)
    # text(x=0,y=0,labels=LETTERS[kk])
  }
  dev.off()
  
  ###
  
  pdf(paste0(figure_folder,"/FigS_14first.dist.surf.MEM_marmap.pdf"),
      width=7.5*2.2/1.2*2,height=12/3*4/1.2/2*4,onefile=T)
  layout(rbind(c(1,2),c(3,4),c(5,6),c(7,8)))
  par(mar = c(3,1,4,1))
  kk=0
  # for (k in c(5,6,1,3,2,10,8,7))
  for (k in 1:14)
  {
    # kk=kk+1
    # Plotting Surface:
    plot(bat, image = T, land = T, lty = 0, bty = "n", xaxt = "n", yaxt = "n", ann = F, bpal = list(c(0, max(bat), greys[1]), c(min(bat), 0, blues[3]))) #plot map without isobaths
    plot(bat, lwd = 0.2, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline
    z = as.data.frame(SUR_dis_MEM_restricted)[rownames(SUR_dis_MEM_restricted) %in% rownames(pie_positions.SUR),k]
    col = rgb(colorRamp(c("blue","grey","red"),space = "Lab")(rescale(z,to=c(0,1),from=c(-max(range(abs(z))),max(range(abs(z)))))),maxColorValue = 255)
    pies.coord = pie_positions.SUR[rownames(pie_positions.SUR) %in% rownames(SUR_dis_MEM_restricted),]
    stations.coord = coord[rownames(coord) %in% rownames(pie_positions.SUR) & rownames(coord) %in% rownames(SUR_dis_MEM_restricted),]
    segments(x0=stations.coord[,2],y0=stations.coord[,1],x1=pies.coord[,1],y1=pies.coord[,2])
    points(pies.coord,
           pch = 21, cex=7, bg = col)
    title(paste("MEM",k), cex.main = 4)
    # signif(charac_dist.SUR_scale[k], digits = 3),"km"), cex.main = 4)
    # text(x=0,y=0,labels=LETTERS[kk])
  }
  dev.off()
  
}

# Table S. - groups carach.
{
  table = cbind(taxo_groups_unmodified,diversity,format(size_relativeAbund,digits=0),format(mean_sim,digits=2),"NA")[selected_groups & diversity>div_threshold,]
  table[!taxo_groups[selected_groups & diversity>div_threshold] %in% c("Dinophyceae","Collodaria"),5] = as.vector(dominant_function1)
  table = cbind(table,round(NormalizedVI_pcoa[[3]][,1:2],digits=2))
  write.table(table, file = paste0(figure_folder,"/TableS1.csv"), quote = F, sep = ";", col.names = F, row.names = F)
  
  # Phototrophs only:
  table1 = cbind(round(NormalizedVI_pcoa[[3]][dominant_function0 %in% c("Other phototrophs","Bacillariophyta"),1:2],digits=3),
                cbind(diversity,round(size_relativeAbund,digits=1))[selected_groups & diversity>100,][dominant_function0 %in% c("Other phototrophs","Bacillariophyta"),])
  colnames(table1) = c("Axis.1","Axis.2","diversity","body size")
  rownames(table1) = c(taxo_groups_unmodified[selected_groups & diversity>100][dominant_function0 %in% c("Other phototrophs","Bacillariophyta")])
  write.table(table1, file = paste0(figure_folder,"/Phototrophs.csv"), quote = F, sep = ";", col.names = T, row.names = T)
  
  table2 = cbind(round(NormalizedVI_pcoa[[3]][dominant_function0 %in% c("Gel. carn. filterers","Other metazoa","Pteropoda","Copepoda") | taxo_groups[selected_groups & diversity>100] == "Collodaria",1:2],digits=3),
    cbind(diversity,round(size_relativeAbund,digits=1))[selected_groups & diversity>100,][dominant_function0 %in% c("Gel. carn. filterers","Other metazoa","Pteropoda","Copepoda") | taxo_groups[selected_groups & diversity>100] == "Collodaria",])
  colnames(table2) = c("Axis.1","Axis.2","diversity","body size")
  rownames(table2) = c(taxo_groups_unmodified[selected_groups & diversity>100][dominant_function0 %in% c("Gel. carn. filterers","Other metazoa","Pteropoda","Copepoda") | taxo_groups[selected_groups & diversity>100] == "Collodaria"])
  write.table(table2, file = paste0(figure_folder,"/Metazoans_Collodaria.csv"), quote = F, sep = ";", col.names = T, row.names = T)
  
}

# Fig. S. adj. R2 hist per MEM
{
  adjr2_individualAxes_MEM = readRDS(paste0(results_folder,"/",
                                            "RDA_Gibbs.prevalence.min.crossValid10sampleFolds_separate.SUR.DCM_individualAxes",
                                            "_abioticPCA_bioticPCA_eigenvalueThres0.8.rds"))[[5]]
  envSelected = readRDS(paste0(results_folder,"/",
                               "RDA_Gibbs.prevalence.min.crossValid10sampleFolds_separate.SUR.DCM_",
                               # "noSelectedAxes-envSelected-significantGlobalModel-nbSites_",
                               "noSelectedAxes-envSelected-nbSites_",
                               # "both.directions.independent.selection",
                               "BH.correction",
                               "_abioticPCA_bioticPCA_eigenvalueThres0.8.rds"))[[2]]
  
  charac_dist.SUR_scale = readRDS(paste0(results_folder,"/MEM.charac.dist.scale_SUR.rds"))
  charac_tmin.SUR_scale = readRDS(paste0(results_folder,"/MEM.charac.tmin.scale_SUR.rds"))
  
  # max.selected.axes = vector(length = length(taxo_groups), mode = "numeric")
  # envSelected.SUR.MEM = list()
  # eigenvalues = list()
  # ii = 0
  # for (i in 1:length(taxo_groups))
  # {
  #   if (i %in% which(selected_groups))
  #   {
  #     ii = ii+1
  #     max.selected.axes[i] = max(which(envSelected$SUR[[i]][[3]]))
  #     envSelected.SUR.MEM[[ii]] = which(envSelected$SUR[[i]][[3]])
  #     eigenvalues[[ii]] = truncated_Tmin_SUR_pcoa$values[envSelected$SUR[[i]][[3]]]
  #   }
  # }
  
  # hist(unlist(envSelected.SUR.MEM),breaks=50)
  # hist(unlist(eigenvalues))
  
  # hist(charac_dist.SUR_scale[1:14])
  # sort(table(as.factor(unlist(envSelected.SUR.MEM))),decreasing=T)
  
  median.adjR2.per.MEM = nb.selecting.groups.per.MEM = vector(length = 14, mode = "numeric")
  adjR2.per.MEM = list()
  for (k in 1:14)
  {
    envSelected_k = rep(F,length(taxo_groups))
    for (i in which(selected_groups))
    {
      envSelected_k[i] = envSelected$SUR[[i]][[3]][k]
    }    
    median.adjR2.per.MEM[k] = median(adjr2_individualAxes_MEM$SUR[envSelected_k,k])
    adjR2.per.MEM[[k]] = adjr2_individualAxes_MEM$SUR[envSelected_k,k]
    nb.selecting.groups.per.MEM[k] = length(adjR2.per.MEM[[k]])
  }
  
  cor.test.MEM.adjR2 = matrix(nrow = 14, ncol = 4, dimnames = list(1:14,c("rho.body.size","p.val.body.size",
                                                                          "p.val.F-test","p.val.t-test.MetaVsPhoto")), data = 0)
  for (k in 1:14)
  {
    cor.test.res = cor.test(adjR2.per.MEM[[k]],log10(size_relativeAbund)[taxo_groups %in% names(adjR2.per.MEM[[k]])])
    cor.test.MEM.adjR2[k,1:2] = c(cor.test.res$estimate,cor.test.res$p.value)
    # aov.res = aov(adjR2.per.MEM[[k]][names(adjR2.per.MEM[[k]]) %in% names(dominant_function1)] ~
    #                 factor(dominant_function1[names(dominant_function1) %in% names(adjR2.per.MEM[[k]])],
    #                        c("Phagotrophs","Phototrophs","Metazoans","Parasites")))
    # cor.test.MEM.adjR2[k,3] = summary(aov.res)[[1]]$`Pr(>F)`[1]
    y = adjR2.per.MEM[[k]][names(adjR2.per.MEM[[k]]) %in% names(dominant_function1)]
    g = factor(dominant_function1[names(dominant_function1) %in% names(adjR2.per.MEM[[k]])],
               c("Phototrophs","Metazoans","Phagotrophs","Parasites"))
    lm.res = lm(y ~ g)
    cor.test.MEM.adjR2[k,3] = anova(lm.res)$`Pr(>F)`[1]
    cor.test.MEM.adjR2[k,4] = coefficients(summary(lm.res))[2,4]
  }
   
  function_colors = c("darkviolet","darkblue","firebrick","firebrick2",
    colorRampPalette(c("darkgreen","chartreuse2"),space = "Lab")(4)[c(1,3,4)],
    "darkgoldenrod1","grey")[c(1,3,6,8,9)]# [c(6,9,8,3,1)]
  
  adjR2.vect.per.MEM = list()
  adj.r2.hist.per.MEM = list()
  kk = 0
  # for (k in sort.int(median.adjR2.per.MEM,decreasing=T,index.return = T)$ix)
  for (k in sort.int(nb.selecting.groups.per.MEM,decreasing=T,index.return = T)$ix)
  {
    kk = kk+1
    
    envSelected_k = rep(F,length(taxo_groups))
    for (i in which(selected_groups))
    {
      envSelected_k[i] = envSelected$SUR[[i]][[3]][k]
    }    
    adjR2.vect.per.MEM[[k]] =  adjr2_individualAxes_MEM$SUR[envSelected_k,k]
    
    function.grouping = vector(length = length(taxo_groups), mode = "character")
    names(function.grouping) = taxo_groups
    function.grouping[names(dominant_function1)] = as.vector(dominant_function1)
    function.grouping[c("Dinophyceae","Collodaria")] = rep("Others",2)
    
    adj.r2.hist.per.MEM[[kk]] = ggplot(data = data.frame(x = adjR2.vect.per.MEM[[k]], 
                                                         group = factor(function.grouping[envSelected_k],
                                                                        levels = c("Phototrophs","Phagotrophs",
                                                                                   "Metazoans","Parasites","Others")))) +
      geom_histogram(aes(x, fill = group),
                     position = "stack",
                     # position = "identity", alpha = 0.5,
                     binwidth = 0.005, center = 0.0025) +
      scale_x_continuous(limits = c(0, 0.16)) +
      # scale_y_discrete(limits = c(0, 18)) +
      scale_fill_manual(values = function_colors) +
      ylim(c(0,16)) +
      # ylim(c(0,8)) +
      labs(y="Number of plankton groups", x="Explained variance") +
      theme_bw() +
      # ggtitle(paste("MEM",k)) +
      ggtitle(paste("MEM",k,"- Char.sc. =",
      signif(charac_dist.SUR_scale[k], digits = 3),"km")) +
                    # signif(charac_tmin.SUR_scale[k], digits = 3),"y")) +
      theme(axis.title=element_text(size=15),
            axis.text = element_text(size=15),
            plot.title=element_text(hjust=0, size=17),
            legend.position = "none",
            plot.margin=unit(c(1,1,1,0.5),"mm")) +
      # Add median adj. R2 annotation:
      annotate(geom="text",
               x=median.adjR2.per.MEM[k] + 0.16*0.05,
               y=0.9*16,
               # label=bquote(atop("Ch.s.s."==.(format(charac_dist.SUR_scale[k],digits=0,nsmall=3))*.km,
               #                   "Median"==.(format(median.adjR2.per.MEM[k],digits=2)))),
               label = bquote("Median"==.(format(median.adjR2.per.MEM[k],digits=2))),
               hjust = 0,
               size=6) +
      # Adding letter labels:
      annotate(geom="text",
               x = 0.16*0.05,
               y=0.9*16,
               label=LETTERS[kk],
               # label = bquote("Median"==.(format(median.adjR2.per.MEM[k],digits=2))),
               hjust = 1,
               size=9) +
      # # Adding charac. scale annotations:
      # annotate(geom="text",
      #          x=0.16*0.5,
      #          y=0.6*16,
      #          label=paste("Charac.scales:\n",
      #                      "+",signif(charac_tmin.SUR_scale[k], digits = 3),"y\n",
      #                      "+",signif(charac_dist.SUR_scale[k], digits = 3),"km"),
      #          # label = bquote("Median"==.(format(median.adjR2.per.MEM[k],digits=2))),
      #          hjust = 0,
      #          size=6) +
      geom_vline(xintercept = median.adjR2.per.MEM[k], linetype = "dashed", size = 0.8)
  }
  
  spl = split(adj.r2.hist.per.MEM[!is.na(adj.r2.hist.per.MEM)], 
              (seq_along(adj.r2.hist.per.MEM[!is.na(adj.r2.hist.per.MEM)])-1) %/% 8)
  ppl = lapply(spl, function(g) marrangeGrob(grobs = g, nrow = 3, ncol = 3, layout_matrix = matrix(data=1:9,nrow=3,byrow=T)))
  pdf(paste0(figure_folder,"/Adj.r2.hist.14.first.MEM_ind.axes.BH.correction_most.often.selected.pdf"),height = 1.5*10*3/4, width = 1.5*10*3/4)
  print(ppl)
  dev.off()
  
  # Extracting the legend manually:
  pdf(paste0(figure_folder,"/Adj.r2.hist.14.first.MEM_get-legend.pdf"),height = 1.5*10/4, width = 1.5*10/4)
  print(adj.r2.hist.per.MEM[[1]] + theme(legend.position = "right",
                                         legend.text = element_text(size=14),
                                         legend.title = element_text(size=16, hjust = 10)) +
          guides(fill = guide_legend(title = "Ecological categories")))
  dev.off()
  
  # Weighted mean adj.R2
  # weighted_mean_MEM_dist_scale = matrix(rep(NA,2*length(taxo_groups)),ncol=2,dimnames=list(taxo_groups,c("SUR","DCM")))
  weighted_mean_MEM_dist_scale = rep(NA,length(taxo_groups))
  adj.r2.vs.charac.dist.plot = list()
  for (i in 1:length(taxo_groups))
  {
    if (i %in% which(selected_groups))
    {
      envSelected.SUR_pos.autocorr = c(envSelected$SUR[[i]][[3]][1:14],rep(F,length(envSelected$SUR[[i]][[3]])-14))
      weighted_mean_MEM_dist_scale[i] = weighted.mean(x = charac_dist.SUR_scale[envSelected.SUR_pos.autocorr],
                                                       w = adjr2_individualAxes_MEM$SUR[i,envSelected.SUR_pos.autocorr])
      # weighted_mean_MEM_dist_scale[i,2] = weighted.mean(x = charac_dist.SUR_scale[envSelected.SUR_pos.autocorr],
      #                                                    w = adjr2_individualAxes_MEM$DCM[i,envSelected.SUR_pos.autocorr])
      # if (i == 1)
      # {
      #   plot(charac_dist.SUR_scale[envSelected.SUR_pos.autocorr],
      #        adjr2_individualAxes_MEM$SUR[i,envSelected.SUR_pos.autocorr])
      # } else
      #   lines(charac_dist.SUR_scale[envSelected.SUR_pos.autocorr],
      #         adjr2_individualAxes_MEM$SUR[i,envSelected.SUR_pos.autocorr],type="p")
      
      adj.r2.vs.charac.dist.plot[[i]] = ggplot(data = data.frame(x = charac_dist.SUR_scale[envSelected.SUR_pos.autocorr],
                                                       y = adjr2_individualAxes_MEM$SUR[i,envSelected.SUR_pos.autocorr])) +
        labs(x="Charac. spatial scale (km)", y="Surf. adj. R2 with MEM") +
        geom_point(aes(x,y), size = 0.6) +
        theme_bw() +
        ggtitle(taxo_groups[i]) +
        theme(axis.title=element_text(size=9),
              plot.title=element_text(hjust=0, size=15),
              plot.margin=unit(c(1,1,1,0.5),"mm"))
    } else
      adj.r2.vs.charac.dist.plot[[i]] = NA
  }
  
  spl = split(adj.r2.vs.charac.dist.plot[!is.na(adj.r2.vs.charac.dist.plot)], 
              (seq_along(adj.r2.vs.charac.dist.plot[!is.na(adj.r2.vs.charac.dist.plot)])-1) %/% 20)
  ppl = lapply(spl, function(g) marrangeGrob(grobs = g, nrow = 4, ncol = 5, layout_matrix = matrix(data=1:20,nrow=4,byrow=T)))
  pdf(paste0(figure_folder,"/Adj.r2.MEM.vs.charac.dist.SUR.pdf"),height = 1.5*10, width = 1.5*10)
  print(ppl)
  dev.off()

}

# Fig. S. adj. R2 hist per abiotic axis:
{
  adjr2_individualAxes_abiotic = readRDS(paste0(results_folder,"/",
                                            "RDA_Gibbs.prevalence.min.crossValid10sampleFolds_separate.SUR.DCM_individualAxes_",
                                            "independent.selection_abioticPCA_bioticPCA_eigenvalueThres0.8_noLagoon.rds"))[[1]]
  envSelected = readRDS(paste0(results_folder,"/",
                               "RDA_Gibbs.prevalence.min.crossValid10sampleFolds_separate.SUR.DCM_",
                               "noSelectedAxes-envSelected-significantGlobalModel-nbSites_",
                               "both.directions.independent.selection_abioticPCA_bioticPCA_eigenvalueThres0.8_noLagoon.rds"))[[2]]
  
  median.adjR2.per.axis = vector(length = 3, mode = "numeric")
  adjR2.per.axis = list()
  for (k in 1:3)
  {
    envSelected_k = rep(F,length(taxo_groups))
    for (i in which(selected_groups))
    {
      envSelected_k[i] = envSelected$SUR[[i]][[1]][k]
    }    
    median.adjR2.per.axis[k] = median(adjr2_individualAxes_abiotic$SUR[envSelected_k,k])
    adjR2.per.axis[[k]] = adjr2_individualAxes_abiotic$SUR[envSelected_k,k]
  }
  
  ###########################################
  # Meta.vs.Photo = c(0,-1,1,0)
  # Phago.vs.Para = c(1,0,0,-1)
  # Phago.Photo.vs.Meta.Para = c(1,1,-1,-1)
  # Phago.vs.Meta = c(1,0,-1,0)
  # Phago.vs.Photo = c(1,-1,0,0)
  # Para.vs.Photo = c(0,-1,0,1)
  cor.test.abiotic.adjR2 = matrix(nrow = 3, ncol = 4, dimnames = list(1:3,c("rho.body.size","p.val.body.size","p.val.F-test","p.val.t-test.MetaVsPhoto")), data = 0)
  for (k in 1:3)
  {
    cor.test.res = cor.test(adjR2.per.axis[[k]],log10(size_relativeAbund)[taxo_groups %in% names(adjR2.per.axis[[k]])])
    cor.test.abiotic.adjR2[k,1:2] = c(cor.test.res$estimate,cor.test.res$p.value)
    y = adjR2.per.axis[[k]][names(adjR2.per.axis[[k]]) %in% names(dominant_function1)]
    
    # g = dominant_function1[names(dominant_function1) %in% names(adjR2.per.axis[[k]])]
    # g1 = g2 = g3 = g4 = vector(length = length(g), mode = "numeric")
    # g1[g == "Phototrophs"] = 1
    # g2[g == "Phagotrophs"] = 1
    # g3[g == "Metazoans"] = 1
    # g4[g == "Parasites"] = 1
    # lm.res = lm(y ~ 1 + g2 + g3 + g4)
    # --> Equivalent to lm(y ~ g) with levels in g ordered as c("Phototrophs","Phagotrophs","Metazoans","Parasites")
    
    g = factor(dominant_function1[names(dominant_function1) %in% names(adjR2.per.axis[[k]])],
               c("Phototrophs","Metazoans","Phagotrophs","Parasites"))
    # contrasts(g) = cbind(Meta.vs.Photo,Phago.vs.Photo,Phago.vs.Meta)
    # contrasts(g) = cbind(Phago.Photo.vs.Meta.Para,Meta.vs.Photo,Phago.vs.Para)
    # contrasts(g) = cbind(Meta.vs.Photo,Phago.vs.Photo,Para.vs.Photo)
    # aov.res = car::Anova(aov(y ~ g), type = "II")
    # cor.test.abiotic.adjR2[k,3] = aov.res$`Pr(>F)`[1]
    lm.res = lm(y ~ g)
    cor.test.abiotic.adjR2[k,3] = anova(lm.res)$`Pr(>F)`[1]
    cor.test.abiotic.adjR2[k,4] = coefficients(summary(lm.res))[2,4]
  }
  #######################################
  
  # y = adjR2.per.axis[[1]][names(adjR2.per.axis[[1]]) %in% names(dominant_function1)]
  # g = factor(dominant_function1[names(dominant_function1) %in% names(adjR2.per.axis[[1]])],
  #            c("Phagotrophs","Phototrophs","Metazoans","Parasites"))
  # Meta.vs.Photo = c(0,1,-1,0)
  # Phago.vs.Meta = c(1,0,-1,0)
  # Phago.vs.Photo = c(1,-1,0,0)
  # contrasts(g) = cbind(Meta.vs.Photo,Phago.vs.Photo,Phago.vs.Meta)
  # res = aov(y ~ g)
  # CList = list("Meta.vs.Photo" = 1,
  #              "Phago.vs.Photo" = 2,
  #              "Phago.vs.Meta" = 3)
  # summary(res,
  #         split=list(g=CList))
  
  #####################################
  function_colors = c("darkviolet","darkblue","firebrick","firebrick2",
                      colorRampPalette(c("darkgreen","chartreuse2"),space = "Lab")(4)[c(1,3,4)],
                      "darkgoldenrod1","grey")[c(1,3,6,8,9)]# [c(6,9,8,3,1)]
  
  adjR2.vect.per.axis = list()
  adj.r2.hist.per.axis = list()
  plot.adj.r2.per.axis.vs.body.size = list()
  boxplot.adj.r2.per.axis = list()
  kk = 0
  for (k in sort.int(median.adjR2.per.axis,decreasing=T,index.return = T)$ix)
  {
    kk = kk+1
    
    envSelected_k = rep(F,length(taxo_groups))
    for (i in which(selected_groups))
    {
      envSelected_k[i] = envSelected$SUR[[i]][[1]][k]
    }    
    adjR2.vect.per.axis[[k]] =  adjr2_individualAxes_abiotic$SUR[envSelected_k,k]
    
    function.grouping = vector(length = length(taxo_groups), mode = "character")
    names(function.grouping) = taxo_groups
    function.grouping[names(dominant_function1)] = as.vector(dominant_function1)
    function.grouping[c("Dinophyceae","Collodaria")] = rep("Others",2)
    
    ###########################
    adj.r2.hist.per.axis[[kk]] = ggplot(data = data.frame(x = adjR2.vect.per.axis[[k]], 
                                                         group = factor(function.grouping[envSelected_k],
                                                                        levels = c("Phototrophs","Phagotrophs",
                                                                                   "Metazoans","Parasites","Others")))) +
      geom_histogram(aes(x, fill = group),
                     position = "stack",
                     # position = "identity", alpha = 0.5,
                     binwidth = 0.005, center = 0.0025) +
      scale_x_continuous(limits = c(0, 0.16)) +
      # scale_y_discrete(limits = c(0, 18)) +
      scale_fill_manual(values = function_colors) +
      ylim(c(0,16)) +
      # ylim(c(0,8)) +
      labs(y="Number of plankton groups", x="Explained variance") +
      theme_bw() +
      ggtitle(paste("Abiotic axis",k)) +
      # ggtitle(paste("MEM",k,"- Char.sc. =",
      #               signif(charac_dist.SUR_scale[k], digits = 3),"km")) +
      # signif(charac_tmin.SUR_scale[k], digits = 3),"y")) +
      theme(axis.title=element_text(size=15),
            axis.text = element_text(size=15),
            plot.title=element_text(hjust=0, size=17),
            legend.position = "none",
            plot.margin=unit(c(1,1,1,0.5),"mm")) +
      # Add median adj. R2 annotation:
      annotate(geom="text",
               x=median.adjR2.per.axis[k] + 0.16*0.05,
               y=0.9*16,
               # label=bquote(atop("Ch.s.s."==.(format(charac_dist.SUR_scale[k],digits=0,nsmall=3))*.km,
               #                   "Median"==.(format(median.adjR2.per.MEM[k],digits=2)))),
               label = bquote("Median"==.(format(median.adjR2.per.axis[k],digits=2))),
               hjust = 0,
               size=6) +
      # Adding letter labels:
      annotate(geom="text",
               x = 0.16*0.05,
               y=0.9*16,
               label=LETTERS[kk],
               # label = bquote("Median"==.(format(median.adjR2.per.MEM[k],digits=2))),
               hjust = 1,
               size=9) +
      # # Adding charac. scale annotations:
      # annotate(geom="text",
      #          x=0.16*0.5,
      #          y=0.6*16,
      #          label=paste("Charac.scales:\n",
      #                      "+",signif(charac_tmin.SUR_scale[k], digits = 3),"y\n",
      #                      "+",signif(charac_dist.SUR_scale[k], digits = 3),"km"),
      #          # label = bquote("Median"==.(format(median.adjR2.per.MEM[k],digits=2))),
      #          hjust = 0,
      #          size=6) +
      geom_vline(xintercept = median.adjR2.per.axis[k], linetype = "dashed", size = 0.8)
    
    ##########################
    plot.adj.r2.per.axis.vs.body.size[[kk]] = cor.plot(x = size_relativeAbund[envSelected_k],
                                                       y = adjR2.vect.per.axis[[k]],
                                                       x.lab = expression("Body size ("*mu*"m)"),
                                                       y.lab = "Explained variance",
                                                       y.lab.hjust = 0.5,
                                                       x.log = T,
                                                       y.log = F,
                                                       fit = T,
                                                       x.cor.pos = 0.8,
                                                       y.cor.pos = 0.8,
                                                       mar.vect = c(5,5,1,5)) + 
      ggtitle("")
    
    boxplot.adj.r2.per.axis[[kk]] = functional.boxplot(classes = factor(dominant_function1[names(dominant_function1) %in% names(adjR2.vect.per.axis[[k]])], 
                                                                        c("Phototrophs","Phagotrophs","Metazoans","Parasites")),
                                                       values = adjR2.vect.per.axis[[k]][!names(adjR2.vect.per.axis[[k]]) %in% c("Dinophyceae","Collodaria")],
                                                       y.lab = "Explained variance",
                                                       y.lab.hjust = 0.39,
                                                       mar.vect = c(5,7,-7,5)) +
      ggtitle(paste("Abiotic axis",k))
    
  }
  
  #####################################
  pdf(paste0(figure_folder,"/Adj.r2.hist.abiotic.pdf"),height = 1.5*10*2/4, width = 1.5*10*2/4)
  grid.arrange(grobs = adj.r2.hist.per.axis, 
               layout_matrix = matrix(c(1:3,NA), nrow = 2, byrow = T))
  dev.off()
  
  #####################################
  pdf(paste0(figure_folder,"/Adj.r2.abiotic.vs.body.size.function.pdf"),
      width=7.5*2.2/1.2/2*2,height=12/3*4/1.2/2*3,onefile=F)
  grid.arrange(grobs = c(boxplot.adj.r2.per.axis,
                         plot.adj.r2.per.axis.vs.body.size), 
               layout_matrix = matrix(1:6, nrow = 3, ncol = 2, byrow = F))
  dev.off()

}

####################################### Marker comparison:

# Fig. S. biogeo. metrics V9 vs psbO
{
  i_case = 1
  
  y = I_square.observed_w.mean_psbO[selected_groups_psbO.4nd,i_case]
  x = I_square.observed_w.mean_V9.psbO[taxo_groups %in% taxo_groups_psbO[selected_groups_psbO.4nd],i_case][order(factor(taxo_groups[taxo_groups %in% taxo_groups_psbO[selected_groups_psbO.4nd]],
                                                                                                     levels=taxo_groups_psbO[selected_groups_psbO.4nd]))]
  plot.autocorr.psbO.vs.V9 = cor.plot(x = x,
                                      y = y,
                                      # y.lab = paste0("Amount of biogeographic structure",ifelse(i_case == 1,""," (DCM)")),
                                      y.lab = "Short-distance spatial autocorr. psbO",
                                      x.lab = "Short-distance spatial autocorr. 18S-V9",
                                      x.log = F,
                                      y.log = F,
                                      # fit = T,
                                      fit = F,
                                      x.cor.pos = 0.8,
                                      y.cor.pos = if (i_case == 1) -0.2 else 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype="dashed")
  
  y = charac_scale_psbO[selected_groups_psbO.4nd,i_case]
  x = charac_scale_V9.psbO[taxo_groups %in% taxo_groups_psbO[selected_groups_psbO.4nd],i_case][order(factor(taxo_groups[taxo_groups %in% taxo_groups_psbO[selected_groups_psbO.4nd]],
                                                                                                levels=taxo_groups_psbO[selected_groups_psbO.4nd]))]
  plot.autocorr.scale.psbO.vs.V9 = cor.plot(x = x,
                                            y = y,
                                            # y.lab = paste(ifelse(i_case == 1,"Scale of biogeographic organization (km)",
                                            #                      "Scale of biogeo. organization(km; DCM)")),
                                            y.lab = "Scale of biogeo. organization psbO (km)",
                                            x.lab = "Scale of biogeo. organization 18S-V9 (km)",
                                            x.log = F,
                                            y.log = F,
                                            # fit = T,
                                            fit = F,
                                            x.cor.pos = 0.8,
                                            y.cor.pos = if (i_case == 1) -0.3 else 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype="dashed")
  
  y = basin_I_psbO[selected_groups_psbO.4nd,i_case]/I_square.observed_w.mean_psbO[selected_groups_psbO.4nd,i_case]
  # x = basin_I[selected_groups & diversity>div_threshold,1]
  x = (basin_I_V9.psbO[taxo_groups %in% taxo_groups_psbO[selected_groups_psbO.4nd],i_case]/
         I_square.observed_w.mean_V9.psbO[taxo_groups %in% taxo_groups_psbO[selected_groups_psbO.4nd],i_case])[order(factor(taxo_groups[taxo_groups %in% taxo_groups_psbO[selected_groups_psbO.4nd]],
                                                                                                                 levels=taxo_groups_psbO[selected_groups_psbO.4nd]))]
  plot.basin.I.within.psbO.vs.V9 = cor.plot(x = x,
                                            y = y,
                                            # y.lab = paste0("Within-basin homogeneity",ifelse(i_case == 1,""," (DCM)")),
                                            y.lab = "Within-basin homogeneity psbO",
                                            x.lab = "Within-basin homogeneity 18S-V9",
                                            x.log = F,
                                            y.log = F,
                                            # fit = T,
                                            fit = F,
                                            x.cor.pos = 0.8,
                                            y.cor.pos = -0.2) +
    geom_abline(slope = 1, intercept = 0, linetype="dashed")
  
  y = lat_I_psbO[selected_groups_psbO.4nd,i_case]/I_square.observed_w.mean_psbO[selected_groups_psbO.4nd,i_case]
  # x = lat_I_sym[selected_groups & diversity>div_threshold,1]
  x = (lat_I_V9.psbO[taxo_groups %in% taxo_groups_psbO[selected_groups_psbO.4nd],i_case]/
         I_square.observed_w.mean_V9.psbO[taxo_groups %in% taxo_groups_psbO[selected_groups_psbO.4nd],i_case])[order(factor(taxo_groups[taxo_groups %in% taxo_groups_psbO[selected_groups_psbO.4nd]],
                                                                                                                levels=taxo_groups_psbO[selected_groups_psbO.4nd]))]
  plot.lat.I.sym.vs.psbO.vs.V9 = cor.plot(x = x,
                                          y = y,
                                          # y.lab = paste0("Latitudinal symmetry",ifelse(i_case == 1,""," (DCM)")),
                                          y.lab = "Latitudinal symmetry psbO",
                                          x.lab = "Latitudinal symmetry 18S-V9",
                                          x.log = F,
                                          y.log = F,
                                          # fit = T,
                                          fit = F,
                                          x.cor.pos = 0.8,
                                          y.cor.pos = -0.1) +
    geom_abline(slope = 1, intercept = 0, linetype="dashed")
  
  pdf(paste0(figure_folder,"/Fig.SI_psbO.4noDino.vs.V9.interpretation_psbO.stations_no.fit_",ifelse(i_case==1,"","DCM_"),div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,".pdf"),
      width=7.5*2.2/1.2*2/2,height=12/3*4/1.2,onefile=F)
  g1 = ggplotGrob(plot.autocorr.psbO.vs.V9)
  g2 = ggplotGrob(plot.autocorr.scale.psbO.vs.V9)
  g3 = ggplotGrob(plot.basin.I.within.psbO.vs.V9)
  g4 = ggplotGrob(plot.lat.I.sym.vs.psbO.vs.V9)
  grid.arrange(grobs = list(g1,g2,g3,g4),
               # widths = rep(1,7),
               # heights = c(0.8,0.8,1),
               layout_matrix = rbind(c(1, 2),c(3, 4)),
               respect=T)
  # grid.draw(arrangeGrob(grobs = list(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11),
  #                       layout_matrix = rbind(c(1, 1, 2),c(3, 3, 4),c(5, 6, 7),c(8, 8, 9),c(10, 10, 11)),
  #                       respect=T))
  dev.off()
}

# Fig. S. biogeo. metrics V9 vs V4
{
  i_case = 1
  
  y = I_square.observed_w.mean_V4[selected_groups_V4,i_case]
  x = I_square.observed_w.mean_V9.V4[taxo_groups %in% taxo_groups_V4[selected_groups_V4],i_case][order(factor(taxo_groups[taxo_groups %in% taxo_groups_V4[selected_groups_V4]],
                                                                                                            levels=taxo_groups_V4[selected_groups_V4]))]
  plot.autocorr.V4.vs.V9 = cor.plot(x = x,
                                      y = y,
                                      # y.lab = paste0("Amount of biogeographic structure",ifelse(i_case == 1,""," (DCM)")),
                                      y.lab = "Short-distance spatial autocorr. 18S-V4",
                                      x.lab = "Short-distance spatial autocorr. 18S-V9",
                                      x.log = F,
                                      y.log = F,
                                      # fit = T,
                                      fit = F,
                                      x.cor.pos = if (i_case == 1) 0.8 else 0.15,
                                      y.cor.pos = if (i_case == 1) 0.1 else 0.9) +
    geom_abline(slope = 1, intercept = 0, linetype="dashed")
  
  y = charac_scale_V4[selected_groups_V4,i_case]
  x = charac_scale_V9.V4[taxo_groups %in% taxo_groups_V4[selected_groups_V4],i_case][order(factor(taxo_groups[taxo_groups %in% taxo_groups_V4[selected_groups_V4]],
                                                                                                levels=taxo_groups_V4[selected_groups_V4]))]
  plot.autocorr.scale.V4.vs.V9 = cor.plot(x = x,
                                            y = y,
                                            # y.lab = paste(ifelse(i_case == 1,"Scale of biogeographic organization (km)",
                                            #                      "Scale of biogeo. organization(km; DCM)")),
                                            y.lab = "Scale of biogeo. organization 18S-V4 (km)",
                                            x.lab = "Scale of biogeo. organization 18S-V9 (km)",
                                            x.log = F,
                                            y.log = F,
                                            # fit = T,
                                            fit = F,
                                            x.cor.pos = if (i_case == 1) 0.2 else 0.85,
                                            y.cor.pos = if (i_case == 1) 0.8 else 0.1) +
    geom_abline(slope = 1, intercept = 0, linetype="dashed")
  
  y = basin_I_V4[selected_groups_V4,i_case]/I_square.observed_w.mean_V4[selected_groups_V4,i_case]
  # x = basin_I[selected_groups & diversity>div_threshold,1]
  x = (basin_I_V9.V4[taxo_groups %in% taxo_groups_V4[selected_groups_V4],i_case]/
         I_square.observed_w.mean_V9.V4[taxo_groups %in% taxo_groups_V4[selected_groups_V4],i_case])[order(factor(taxo_groups[taxo_groups %in% taxo_groups_V4[selected_groups_V4]],
                                                                                                                levels=taxo_groups_V4[selected_groups_V4]))]
  plot.basin.I.within.V4.vs.V9 = cor.plot(x = x,
                                            y = y,
                                            # y.lab = paste0("Within-basin homogeneity",ifelse(i_case == 1,""," (DCM)")),
                                            y.lab = "Within-basin homogeneity 18S-V4",
                                            x.lab = "Within-basin homogeneity 18S-V9",
                                            x.log = F,
                                            y.log = F,
                                            # fit = T,
                                            fit = F,
                                            x.cor.pos = if (i_case == 1) 0.8 else 0.15,
                                            y.cor.pos = if (i_case == 1) 0.8 else 0.9) +
    geom_abline(slope = 1, intercept = 0, linetype="dashed")
  
  y = lat_I_V4[selected_groups_V4,i_case]/I_square.observed_w.mean_V4[selected_groups_V4,i_case]
  # x = lat_I_sym[selected_groups & diversity>div_threshold,1]
  x = (lat_I_V9.V4[taxo_groups %in% taxo_groups_V4[selected_groups_V4],i_case]/
         I_square.observed_w.mean_V9.V4[taxo_groups %in% taxo_groups_V4[selected_groups_V4],i_case])[order(factor(taxo_groups[taxo_groups %in% taxo_groups_V4[selected_groups_V4]],
                                                                                                                levels=taxo_groups_V4[selected_groups_V4]))]
  plot.lat.I.sym.vs.V4.vs.V9 = cor.plot(x = x,
                                        y = y,
                                        # y.lab = paste0("Latitudinal symmetry",ifelse(i_case == 1,""," (DCM)")),
                                        y.lab = "Latitudinal symmetry 18S-V4",
                                        x.lab = "Latitudinal symmetry 18S-V9",
                                        x.log = F,
                                        y.log = F,
                                        # fit = T,
                                        fit = F,
                                        x.cor.pos = if (i_case == 1) 0.9 else 0.15,
                                        y.cor.pos = if (i_case == 1) 0.1 else 0.1) +
    geom_abline(slope = 1, intercept = 0, linetype="dashed")
  
  pdf(paste0(figure_folder,"/Fig.SI_V4.vs.V9.interpretation_V4.stations_no.fit_",ifelse(i_case==1,"","DCM_"),div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,".pdf"),
      width=7.5*2.2/1.2*2/2,height=12/3*4/1.2,onefile=F)
  g1 = ggplotGrob(plot.autocorr.V4.vs.V9)
  g2 = ggplotGrob(plot.autocorr.scale.V4.vs.V9)
  g3 = ggplotGrob(plot.basin.I.within.V4.vs.V9)
  g4 = ggplotGrob(plot.lat.I.sym.vs.V4.vs.V9)
  grid.arrange(grobs = list(g1,g2,g3,g4),
               # widths = rep(1,7),
               # heights = c(0.8,0.8,1),
               layout_matrix = rbind(c(1, 2),c(3, 4)),
               respect=T)
  # grid.draw(arrangeGrob(grobs = list(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11),
  #                       layout_matrix = rbind(c(1, 1, 2),c(3, 3, 4),c(5, 6, 7),c(8, 8, 9),c(10, 10, 11)),
  #                       respect=T))
  dev.off()
}

# Fig. S. biogeo. metrics V9 vs V4 and psbO
{
  i_case = 1
  
  y = I_square.observed_w.mean_V4[selected_groups_V4,i_case]
  x = I_square.observed_w.mean_V9.V4[taxo_groups %in% taxo_groups_V4[selected_groups_V4],i_case][order(factor(taxo_groups[taxo_groups %in% taxo_groups_V4[selected_groups_V4]],
                                                                                                              levels=taxo_groups_V4[selected_groups_V4]))]
  x.range = range(x)
  y.range = range(y)
  plot.autocorr.V4.vs.V9 = cor.plot(x = x,
                                    y = y,
                                    # y.lab = paste0("Amount of biogeographic structure",ifelse(i_case == 1,""," (DCM)")),
                                    y.lab = "Short-distance spatial autocorr. - V4",
                                    x.lab = "Short-distance spatial autocorr. - V9",
                                    x.log = F,
                                    y.log = F,
                                    # fit = T,
                                    fit = F,
                                    x.cor.pos = if (i_case == 1) 0.8 else 0.15,
                                    y.cor.pos = if (i_case == 1) 0.1 else 0.9,
                                    mar.vect=c(5,15,5,5)) +
                                    # mar.vect=c(15,10,5,5)) +
    geom_abline(slope = 1, intercept = 0, linetype="dashed")
  
  y = I_square.observed_w.mean_psbO[selected_groups_psbO.4nd,i_case]
  x = I_square.observed_w.mean_V9.psbO[taxo_groups %in% taxo_groups_psbO[selected_groups_psbO.4nd],i_case][order(factor(taxo_groups[taxo_groups %in% taxo_groups_psbO[selected_groups_psbO.4nd]],
                                                                                                                        levels=taxo_groups_psbO[selected_groups_psbO.4nd]))]
  plot.autocorr.psbO.vs.V9 = cor.plot(x = x,
                                      y = y,
                                      # y.lab = paste0("Amount of biogeographic structure",ifelse(i_case == 1,""," (DCM)")),
                                      y.lab = "Short-distance spatial autocorr. - psbO",
                                      x.lab = "Short-distance spatial autocorr. - V9",
                                      x.log = F,
                                      y.log = F,
                                      # fit = T,
                                      fit = F,
                                      x.cor.pos = 0.8,
                                      y.cor.pos = if (i_case == 1) -0.2 else 1.2,
                                      mar.vect=c(5,15,5,5)) +
                                      # mar.vect=c(15,10,5,5)) +
    xlim(x.range[1],max(x.range[2],max(x))) +
    ylim(y.range[1:2]) +
    geom_abline(slope = 1, intercept = 0, linetype="dashed")
  
  y = charac_scale_V4[selected_groups_V4,i_case]
  x = charac_scale_V9.V4[taxo_groups %in% taxo_groups_V4[selected_groups_V4],i_case][order(factor(taxo_groups[taxo_groups %in% taxo_groups_V4[selected_groups_V4]],
                                                                                                  levels=taxo_groups_V4[selected_groups_V4]))]
  x.range = range(x)
  y.range = range(y)
  plot.autocorr.scale.V4.vs.V9 = cor.plot(x = x,
                                          y = y,
                                          # y.lab = paste(ifelse(i_case == 1,"Scale of biogeographic organization (km)",
                                          #                      "Scale of biogeo. organization(km; DCM)")),
                                          y.lab = "Scale of biogeo. organization (km) - V4 ",
                                          x.lab = "Scale of biogeo. organization (km) - V9 ",
                                          x.log = F,
                                          y.log = F,
                                          # fit = T,
                                          fit = F,
                                          x.cor.pos = if (i_case == 1) 0.2 else 0.85,
                                          y.cor.pos = if (i_case == 1) 0.8 else 0.1,
                                          mar.vect=c(5,15,5,5)) +
                                          # mar.vect=c(15,10,5,5)) +
    geom_abline(slope = 1, intercept = 0, linetype="dashed")
  
  y = charac_scale_psbO[selected_groups_psbO.4nd,i_case]
  x = charac_scale_V9.psbO[taxo_groups %in% taxo_groups_psbO[selected_groups_psbO.4nd],i_case][order(factor(taxo_groups[taxo_groups %in% taxo_groups_psbO[selected_groups_psbO.4nd]],
                                                                                                            levels=taxo_groups_psbO[selected_groups_psbO.4nd]))]
  plot.autocorr.scale.psbO.vs.V9 = cor.plot(x = x,
                                            y = y,
                                            # y.lab = paste(ifelse(i_case == 1,"Scale of biogeographic organization (km)",
                                            #                      "Scale of biogeo. organization(km; DCM)")),
                                            y.lab = "Scale of biogeo. organization (km) - psbO",
                                            x.lab = "Scale of biogeo. organization (km) - V9",
                                            x.log = F,
                                            y.log = F,
                                            # fit = T,
                                            fit = F,
                                            x.cor.pos = 0.8,
                                            y.cor.pos = if (i_case == 1) -0.3 else 1.2,
                                            mar.vect=c(5,15,5,5)) +
    # mar.vect=c(15,10,5,5)) +
    xlim(x.range[1:2]) +
    ylim(y.range[1:2]) +
    geom_abline(slope = 1, intercept = 0, linetype="dashed")
  
  y = basin_I_V4[selected_groups_V4,i_case]/I_square.observed_w.mean_V4[selected_groups_V4,i_case]
  # x = basin_I[selected_groups & diversity>div_threshold,1]
  x = (basin_I_V9.V4[taxo_groups %in% taxo_groups_V4[selected_groups_V4],i_case]/
         I_square.observed_w.mean_V9.V4[taxo_groups %in% taxo_groups_V4[selected_groups_V4],i_case])[order(factor(taxo_groups[taxo_groups %in% taxo_groups_V4[selected_groups_V4]],
                                                                                                                  levels=taxo_groups_V4[selected_groups_V4]))]
  x.range = range(x)
  y.range = range(y)
  plot.basin.I.within.V4.vs.V9 = cor.plot(x = x,
                                          y = y,
                                          # y.lab = paste0("Within-basin homogeneity",ifelse(i_case == 1,""," (DCM)")),
                                          y.lab = "Within-basin homogeneity - V4",
                                          x.lab = "Within-basin homogeneity - V9",
                                          x.log = F,
                                          y.log = F,
                                          # fit = T,
                                          fit = F,
                                          x.cor.pos = if (i_case == 1) 0.8 else 0.15,
                                          y.cor.pos = if (i_case == 1) 0.8 else 0.9,
                                          mar.vect=c(5,15,5,5)) +
                                          # mar.vect=c(15,10,5,5)) +
    geom_abline(slope = 1, intercept = 0, linetype="dashed")
  
  y = basin_I_psbO[selected_groups_psbO.4nd,i_case]/I_square.observed_w.mean_psbO[selected_groups_psbO.4nd,i_case]
  # x = basin_I[selected_groups & diversity>div_threshold,1]
  x = (basin_I_V9.psbO[taxo_groups %in% taxo_groups_psbO[selected_groups_psbO.4nd],i_case]/
         I_square.observed_w.mean_V9.psbO[taxo_groups %in% taxo_groups_psbO[selected_groups_psbO.4nd],i_case])[order(factor(taxo_groups[taxo_groups %in% taxo_groups_psbO[selected_groups_psbO.4nd]],
                                                                                                                            levels=taxo_groups_psbO[selected_groups_psbO.4nd]))]
  plot.basin.I.within.psbO.vs.V9 = cor.plot(x = x,
                                            y = y,
                                            # y.lab = paste0("Within-basin homogeneity",ifelse(i_case == 1,""," (DCM)")),
                                            y.lab = "Within-basin homogeneity - psbO",
                                            x.lab = "Within-basin homogeneity - V9",
                                            x.log = F,
                                            y.log = F,
                                            # fit = T,
                                            fit = F,
                                            x.cor.pos = 0.8,
                                            y.cor.pos = -0.2,
                                            mar.vect=c(5,15,5,5)) +
                                            # mar.vect=c(15,10,5,5)) +
    xlim(x.range[1:2]) +
    ylim(y.range[1:2]) +
    geom_abline(slope = 1, intercept = 0, linetype="dashed")
  
  y = lat_I_V4[selected_groups_V4,i_case]/I_square.observed_w.mean_V4[selected_groups_V4,i_case]
  # x = lat_I_sym[selected_groups & diversity>div_threshold,1]
  x = (lat_I_V9.V4[taxo_groups %in% taxo_groups_V4[selected_groups_V4],i_case]/
         I_square.observed_w.mean_V9.V4[taxo_groups %in% taxo_groups_V4[selected_groups_V4],i_case])[order(factor(taxo_groups[taxo_groups %in% taxo_groups_V4[selected_groups_V4]],
                                                                                                                  levels=taxo_groups_V4[selected_groups_V4]))]
  x.range = range(x)
  y.range = range(y)
  plot.lat.I.sym.vs.V4.vs.V9 = cor.plot(x = x,
                                        y = y,
                                        # y.lab = paste0("Latitudinal symmetry",ifelse(i_case == 1,""," (DCM)")),
                                        y.lab = "Latitudinal symmetry - V4",
                                        x.lab = "Latitudinal symmetry - V9",
                                        x.log = F,
                                        y.log = F,
                                        # fit = T,
                                        fit = F,
                                        x.cor.pos = if (i_case == 1) 0.9 else 0.15,
                                        y.cor.pos = if (i_case == 1) 0.1 else 0.1,
                                        mar.vect=c(5,15,5,5)) +
    geom_abline(slope = 1, intercept = 0, linetype="dashed")
  
  y = lat_I_psbO[selected_groups_psbO.4nd,i_case]/I_square.observed_w.mean_psbO[selected_groups_psbO.4nd,i_case]
  # x = lat_I_sym[selected_groups & diversity>div_threshold,1]
  x = (lat_I_V9.psbO[taxo_groups %in% taxo_groups_psbO[selected_groups_psbO.4nd],i_case]/
         I_square.observed_w.mean_V9.psbO[taxo_groups %in% taxo_groups_psbO[selected_groups_psbO.4nd],i_case])[order(factor(taxo_groups[taxo_groups %in% taxo_groups_psbO[selected_groups_psbO.4nd]],
                                                                                                                            levels=taxo_groups_psbO[selected_groups_psbO.4nd]))]
  plot.lat.I.sym.vs.psbO.vs.V9 = cor.plot(x = x,
                                          y = y,
                                          # y.lab = paste0("Latitudinal symmetry",ifelse(i_case == 1,""," (DCM)")),
                                          y.lab = "Latitudinal symmetry - psbO",
                                          x.lab = "Latitudinal symmetry - V9",
                                          x.log = F,
                                          y.log = F,
                                          # fit = T,
                                          fit = F,
                                          x.cor.pos = 0.8,
                                          y.cor.pos = -0.1,
                                          mar.vect=c(5,15,5,5)) +
    xlim(x.range[1:2]) +
    ylim(y.range[1:2]) +
    geom_abline(slope = 1, intercept = 0, linetype="dashed")
  
  pdf(paste0(figure_folder,"/Fig.SI_V9.psbO.4noDino_V9.V4_interpretation_no.fit_",ifelse(i_case==1,"","DCM_"),div_threshold,"plusOTUs",noArcticNoBiomark_insert,noLagoon_insert,".pdf"),
      width=7.5*2.2/1.3/2*2,height=12/3*2/1.2*4,onefile=F)
      # width=7.5*2.2/1.2/2*2,height=12/3*2/1.2*4,onefile=F)
  g1 = ggplotGrob(plot.autocorr.V4.vs.V9)
  g2 = ggplotGrob(plot.autocorr.psbO.vs.V9)
  g3 = ggplotGrob(plot.autocorr.scale.V4.vs.V9)
  g4 = ggplotGrob(plot.autocorr.scale.psbO.vs.V9)
  g5 = ggplotGrob(plot.basin.I.within.V4.vs.V9)
  g6 = ggplotGrob(plot.basin.I.within.psbO.vs.V9)
  g7 = ggplotGrob(plot.lat.I.sym.vs.V4.vs.V9)
  g8 = ggplotGrob(plot.lat.I.sym.vs.psbO.vs.V9)
  grid.arrange(grobs = list(g1,g2,g3,g4,g5,g6,g7,g8),
               # widths = rep(1,7),
               # heights = c(0.8,0.8,1),
               layout_matrix = rbind(c(1, 2),c(3, 4),c(5, 6),c(7, 8)),
               respect=F)
  # grid.draw(arrangeGrob(grobs = list(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11),
  #                       layout_matrix = rbind(c(1, 1, 2),c(3, 3, 4),c(5, 6, 7),c(8, 8, 9),c(10, 10, 11)),
  #                       respect=T))
  dev.off()
}

# Marker comparison
{
  # x = I_square.observed_w.mean_V4[selected_groups_V4,1]
  x = (I_square.observed_w.mean_V4[selected_groups_V4,1] + I_square.observed_w.mean[taxo_groups_V4[selected_groups_V4],1])/2
  V4.V9_dis = Normalized_VI_V9.V4[selected_groups_V4,selected_groups_V4]
  # y = diag(V4.V9_dis)
  # null.dis = mean(mean(V4.V9_dis[lower.tri(V4.V9_dis)],na.rm=T),mean(V4.V9_dis[upper.tri(V4.V9_dis)],na.rm=T))
  V4.V9_dis_no.diag = V4.V9_dis
  diag(V4.V9_dis_no.diag) = NA
  null.dis = colMeans(V4.V9_dis_no.diag,na.rm=T)
  y = (null.dis - diag(V4.V9_dis))/(1-null.dis)
  plot.VI.V9.V4.vs.Moran.I.V9.V4 = cor.plot(x = x,
                                            y = y,
                                            x.lab = "Mean spatial autocorrelation of V4 and V9",
                                            # y.lab = "V9-V4 biogeographic dissimilarity",
                                            y.lab = "V9-V4 similarity relative to null expectation",
                                            y.lab.hjust = 0.5,
                                            x.log = F,
                                            y.log = F,
                                            fit = F,
                                            # fit.display = "pearson.spearman",
                                            mar.vect = c(5,10,1,5)) +
    # geom_segment(aes(x = min(x), xend = max(x), y = null.dis, yend = null.dis), linetype="dashed")
    geom_segment(aes(x = min(x), xend = max(x), y = 0, yend = 0), linetype="dashed")
  # pdf(paste0(figure_folder,"/FigS_V9-V4.VI.dissimilarity_vs_V9-V4.average.Moran.I.pdf"),
  pdf(paste0(figure_folder,"/FigS_V9-V4.VI.similarity.relative.to.null.expectation_vs_V9-V4.average.Moran.I.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.VI.V9.V4.vs.Moran.I.V9.V4)
  dev.off()
  
  x = (I_square.observed_w.mean_V4[selected_groups_V4,1] + I_square.observed_w.mean[taxo_groups_V4[selected_groups_V4],1])/2
  V4.V9_dis = Normalized_VI_V9.V4[selected_groups_V4,selected_groups_V4]
  y = 1 - diag(V4.V9_dis)
  y.range_V9.V4 = range(y)
  x.range_V9.V4 = range(x)
  col.factor = rep("black",length(x))
  col.factor[taxo_groups_V4[selected_groups_V4] %in% c("Collodaria","Bacillariophyta","MAST-3,_12","Dinophyceae","Haptophyta")] = "blue"
  plot.VI.V9.V4.vs.Moran.I.V9.V4 = cor.plot(x = x,
                                            y = y,
                                            col.factor = col.factor,
                                            col.vect = c("black","blue"),
                                            leg.name = NULL,
                                            x.lab = "Mean spatial autocorrelation of V4 and V9",
                                            # y.lab = "V9-V4 biogeographic dissimilarity",
                                            y.lab = "V9-V4 similarity",
                                            y.lab.hjust = 0.5,
                                            x.log = F,
                                            y.log = F,
                                            fit = F,
                                            # fit.display = "pearson.spearman",
                                            mar.vect = c(5,10,1,5))
    # geom_segment(aes(x = min(x), xend = max(x), y = null.dis, yend = null.dis), linetype="dashed")
    # geom_segment(aes(x = min(x), xend = max(x), y = 0, yend = 0), linetype="dashed")
  # pdf(paste0(figure_folder,"/FigS_V9-V4.VI.dissimilarity_vs_V9-V4.average.Moran.I.pdf"),
  pdf(paste0(figure_folder,"/FigS_V9-V4.VI.similarity_vs_V9-V4.average.Moran.I.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.VI.V9.V4.vs.Moran.I.V9.V4)
  dev.off()
  
  # x = I_square.observed_w.mean_psbO[selected_groups_psbO,1]
  x = (I_square.observed_w.mean_psbO[selected_groups_psbO.4nd,1] + I_square.observed_w.mean[taxo_groups_psbO[selected_groups_psbO.4nd],1])/2
  psbO.V9_dis = Normalized_VI_V9.psbO[selected_groups_psbO.4nd,selected_groups_psbO.4nd]
  # y = diag(psbO.V9_dis)
  # null.dis = mean(mean(psbO.V9_dis[lower.tri(psbO.V9_dis)],na.rm=T),mean(psbO.V9_dis[upper.tri(psbO.V9_dis)],na.rm=T))
  psbO.V9_dis_no.diag = psbO.V9_dis
  diag(psbO.V9_dis_no.diag) = NA
  null.dis = colMeans(psbO.V9_dis_no.diag,na.rm=T)
  y = (null.dis - diag(psbO.V9_dis))/(1-null.dis)
  plot.VI.V9.psbO.vs.Moran.I.V9.psbO = cor.plot(x = x,
                                                y = y,
                                                col.vect = col.vect,
                                                x.lab = "Mean spatial autocorrelation of psbO and V9",
                                                # y.lab = "V9-psbO biogeographic dissimilarity",
                                                y.lab = "V9-psbO similarity relative to null expectation",
                                                y.lab.hjust = 0.5,
                                                x.log = F,
                                                y.log = F,
                                                fit = F,
                                                # fit.display = "pearson.spearman",
                                                mar.vect = c(5,15,1,5)) +
    geom_segment(aes(x = min(x), xend = max(x), y = 0, yend = 0), linetype="dashed")
  # geom_segment(aes(x = min(x), xend = max(x), y = null.dis, yend = null.dis), linetype="dashed")
  # pdf(paste0(figure_folder,"/FigS_v9-psbO.VI.dissimilarity_vs_V9-psbO.average.Moran.I.pdf"),
  pdf(paste0(figure_folder,"/FigS_v9-psbO.4noDino.VI.similarity.relative.to.null.expectation_vs_V9-psbO.4noDino.average.Moran.I.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.VI.V9.psbO.vs.Moran.I.V9.psbO)
  dev.off()
  
  x = (I_square.observed_w.mean_psbO[selected_groups_psbO.4nd,1] + I_square.observed_w.mean[taxo_groups_psbO[selected_groups_psbO.4nd],1])/2
  psbO.V9_dis = Normalized_VI_V9.psbO[selected_groups_psbO.4nd,selected_groups_psbO.4nd]
  y = 1 - diag(psbO.V9_dis)
  col.factor = rep("black",length(x))
  col.factor[taxo_groups_psbO[selected_groups_psbO.4nd] %in% c("Bacillariophyta","Haptophyta")] = "blue"
  plot.VI.V9.psbO.vs.Moran.I.V9.psbO = cor.plot(x = x,
                                                y = y,
                                                col.factor = col.factor,
                                                col.vect = c("black","blue"),
                                                leg.name = NULL,
                                                x.lab = "Mean spatial autocorrelation of psbO and V9",
                                                # y.lab = "V9-psbO biogeographic dissimilarity",
                                                y.lab = "V9-psbO similarity",
                                                y.lab.hjust = 0.5,
                                                x.log = F,
                                                y.log = F,
                                                fit = F,
                                                # fit.display = "pearson.spearman",
                                                mar.vect = c(5,15,1,5)) +
    ylim(y.range_V9.V4[1:2]) +
    xlim(x.range_V9.V4[1:2])
  # geom_segment(aes(x = min(x), xend = max(x), y = null.dis, yend = null.dis), linetype="dashed")
  # pdf(paste0(figure_folder,"/FigS_v9-psbO.VI.dissimilarity_vs_V9-psbO.average.Moran.I.pdf"),
  pdf(paste0(figure_folder,"/FigS_v9-psbO.4noDino.VI.similarity_vs_V9-psbO.4noDino.average.Moran.I.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.VI.V9.psbO.vs.Moran.I.V9.psbO)
  dev.off()
}

# Group dissimilarity within marker
{
  x = Normalized_VI[[1]]
  x[upper.tri(x)] = t(x)[upper.tri(t(x))]
  x = x[taxo_groups_V4[selected_groups_V4],taxo_groups_V4[selected_groups_V4]]
  x = 1 - x[lower.tri(x)]
  x.range_V9.V4 = range(x)
  y = Normalized_VI_V4[selected_groups_V4,selected_groups_V4]
  y = 1 - y[lower.tri(y)]
  y.range_V9.V4 = range(y)
  plot.VI.V4.vs.VI.V9 = cor.plot(x = x,
                                 y = y,
                                 col.vect = col.vect,
                                 x.lab = "Between-group similarity - V9",
                                 # y.lab = "V9-V4 biogeographic dissimilarity",
                                 y.lab = "Between-group similarity - V4",
                                 y.lab.hjust = 0.5,
                                 x.log = F,
                                 y.log = F,
                                 fit = F,
                                 # fit.display = "pearson.spearman",
                                 mar.vect = c(5,10,1,5)) +
    geom_abline(slope = 1, intercept = 0, linetype="dashed")
    
  pdf(paste0(figure_folder,"/FigS_group.SUR.similarity_V4_vs_V9.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.VI.V4.vs.VI.V9)
  dev.off()
  
  x = Normalized_VI[[1]]
  x[upper.tri(x)] = t(x)[upper.tri(t(x))]
  x = x[taxo_groups_psbO[selected_groups_psbO.4nd],taxo_groups_psbO[selected_groups_psbO.4nd]]
  x = 1 - x[lower.tri(x)]
  y = Normalized_VI_psbO[selected_groups_psbO.4nd,selected_groups_psbO.4nd]
  y = 1 - y[lower.tri(y)]
  plot.VI.psbO.vs.VI.V9 = cor.plot(x = x,
                                   y = y,
                                   col.vect = col.vect,
                                   x.lab = "Between-group similarity - V9",
                                   # y.lab = "V9-V4 biogeographic dissimilarity",
                                   y.lab = "Between-group similarity - psbO",
                                   y.lab.hjust = 0.5,
                                   x.log = F,
                                   y.log = F,
                                   fit = F,
                                   # fit.display = "pearson.spearman",
                                   mar.vect = c(5,10,1,5)) +
    xlim(x.range_V9.V4[1:2]) +
    ylim(y.range_V9.V4[1:2]) +
    geom_abline(slope = 1, intercept = 0, linetype="dashed") 
  
  pdf(paste0(figure_folder,"/FigS_group.SUR.similarity_psbO.4.noDino_vs_V9.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.VI.psbO.vs.VI.V9)
  dev.off()
}

######################################### PCoA axes justification:

# Computing distances back from PCoA axes
{
  original.distances = Normalized_VI[[3]][selected_groups,selected_groups]
  cor_pearson_distances = cor_spearman_distances = matrix(nrow = ncol(NormalizedVI_pcoa[[3]]), ncol = 2, data = 0, 
                                                          dimnames = list(colnames(NormalizedVI_pcoa[[3]]),c("rho","p.val")))
  for (i in 1:ncol(NormalizedVI_pcoa[[3]]))
  {
    distances.from.axes_i = as.matrix(dist(NormalizedVI_pcoa[[3]][,1:i]))
    test.p = cor.test(distances.from.axes_i[lower.tri(distances.from.axes_i)], original.distances[lower.tri(original.distances)])
    cor_pearson_distances[i,] = c(test.p$estimate, test.p$p.value)
    test.s = cor.test(distances.from.axes_i[lower.tri(distances.from.axes_i)], original.distances[lower.tri(original.distances)], method = "spearman")
    cor_spearman_distances[i,] = c(test.s$estimate, test.s$p.value)
  }
  
  x = 1:ncol(NormalizedVI_pcoa[[3]])
  y = cor_pearson_distances[,1]
  col.vect = rep("black",length(x))
  col.vect[cor_pearson_distances[,2]>0.05] = "red"
  plot.cor.pearson.distances.vs.nb.axes = cor.plot(x = x,
                                                   y = y,
                                                   col.vect = col.vect,
                                                   x.lab = "Number of PCoA axes retained",
                                                   y.lab = "Pearson correlation\n between axes-based and original distances",
                                                   y.lab.hjust = 0.5,
                                                   x.log = F,
                                                   y.log = F,
                                                   fit = F,
                                                   # fit.display = "pearson.spearman",
                                                   mar.vect = c(5,5,1,5))
  
  x = VI.pcoa$values$Cumul_eig
  y = cor_pearson_distances[,1]
  # col.vect = rep("black",length(x))
  # col.vect[cor_pearson_distances[,2]>0.05] = "red"
  col.factor = c(rep("blue",2),rep("black",length(x)-2))
  plot.cor.pearson.distances.vs.explained.var = cor.plot(x = x,
                                                         y = y,
                                                         col.factor = col.factor,
                                                         col.vect = c("black","blue"),
                                                         leg.name=NULL,
                                                         x.lab = "PCoA axes cumulative variance fraction",
                                                         y.lab = "Pearson correlation\n between axes-based and original distances",
                                                         y.lab.hjust = 0.5,
                                                         x.log = F,
                                                         y.log = F,
                                                         fit = F,
                                                         # fit.display = "pearson.spearman",
                                                         mar.vect = c(5,8,1,5)) + ylim(0,1)
  
  # x = 1:ncol(NormalizedVI_pcoa[[3]])
  # y = cor_spearman_distances[,1]
  # col.vect = rep("black",length(x))
  # col.vect[cor_spearman_distances[,2]>0.05] = "red"
  # plot.cor.spearman.distances.vs.nb.axes = cor.plot(x = x,
  #                                          y = y,
  #                                          col.vect = col.vect,
  #                                          x.lab = "Number of PCoA axes retained",
  #                                          y.lab = "Spearman correlation between distances",
  #                                          y.lab.hjust = 0.5,
  #                                          x.log = F,
  #                                          y.log = F,
  #                                          fit = F,
  #                                          # fit.display = "pearson.spearman",
  #                                          mar.vect = c(5,5,1,5))
  
  x = original.distances#[lower.tri(original.distances)]
  distances.from.axes = as.matrix(dist(NormalizedVI_pcoa[[3]][,1:2]))
  y =  distances.from.axes#[lower.tri(distances.from.axes)]
  plot.distances.first.2.axes.vs.original.distances = cor.plot(x = x,
                                                               y = y,
                                                               x.lab = "Original distances between groups",
                                                               y.lab = "Distances based on first 2 PCoA axes",
                                                               y.lab.hjust = 0.5,
                                                               x.log = F,
                                                               y.log = F,
                                                               fit = F,
                                                               fit.display = "pearson.spearman",
                                                               x.cor.pos=0.2,y.cor.pos=0.8,
                                                               dis = T,
                                                               mar.vect = c(5,8,1,5))
  
  x = as.matrix(dist(NormalizedVI_pcoa[[3]]))
  y = as.matrix(dist(NormalizedVI_pcoa[[3]][,1:2]))
  plot.distances.first.2.axes.vs.all.axes = cor.plot(x = x,
                                                     y = y,
                                                     x.lab = "Distances between groups\n based on all axes",
                                                     y.lab = "Distances between groups\n based on the first 2 PCoA axes",
                                                     y.lab.hjust = 0.5,
                                                     x.log = F,
                                                     y.log = F,
                                                     fit = F,
                                                     # fit.display = "pearson.spearman",
                                                     x.cor.pos=0.2,y.cor.pos=0.8,
                                                     dis = T,
                                                     mar.vect = c(5,8,1,5))
  
  original.distances.200 = original.distances[diversity[selected_groups] > 200,diversity[selected_groups] > 200]
  x = original.distances.200[lower.tri(original.distances.200)]
  distances.from.axes.200 = as.matrix(dist(NormalizedVI_pcoa[[3]][diversity[selected_groups] > 200,1:2]))
  y =  distances.from.axes.200[lower.tri(distances.from.axes.200)]
  plot.distances.200.first.2.axes.vs.original.distances = cor.plot(x = x,
                                                                   y = y,
                                                                   x.lab = "Original VI distances between >200-OTU groups",
                                                                   y.lab = "Distances based on first 2 PCoA axes",
                                                                   y.lab.hjust = 0.5,
                                                                   x.log = F,
                                                                   y.log = F,
                                                                   fit = T,
                                                                   fit.display = "pearson.spearman",
                                                                   x.cor.pos=0.2,y.cor.pos=0.8,
                                                                   mar.vect = c(5,5,1,5))
  
  # x = original.distances[lower.tri(original.distances)]
  # distances.from.axes = as.matrix(dist(NormalizedVI_pcoa[[3]][,1:9]))
  # y =  distances.from.axes[lower.tri(distances.from.axes)]
  # plot.distances.first.9.axes.vs.original.distances = cor.plot(x = x,
  #                                                              y = y,
  #                                                              x.lab = "Original VI distances between groups",
  #                                                              y.lab = "Distances based on first 9 PCoA axes",
  #                                                              y.lab.hjust = 0.5,
  #                                                              x.log = F,
  #                                                              y.log = F,
  #                                                              fit = T,
  #                                                              fit.display = "pearson.spearman",
  #                                                              x.cor.pos=0.2,y.cor.pos=0.8,
  #                                                              mar.vect = c(5,5,1,5))
  
  original.distances.500 = original.distances[diversity[selected_groups] > 500,diversity[selected_groups] > 500]
  x = original.distances.500[lower.tri(original.distances.500)]
  distances.from.axes.500 = as.matrix(dist(NormalizedVI_pcoa[[3]][diversity[selected_groups] > 500,1:2]))
  y =  distances.from.axes.500[lower.tri(distances.from.axes.500)]
  plot.distances.500.first.2.axes.vs.original.distances = cor.plot(x = x,
                                                                   y = y,
                                                                   x.lab = "Original distances between >500-OTU groups",
                                                                   y.lab = "Distances based on first 2 PCoA axes",
                                                                   y.lab.hjust = 0.5,
                                                                   x.log = F,
                                                                   y.log = F,
                                                                   fit = T,
                                                                   fit.display = "pearson.spearman",
                                                                   x.cor.pos=0.2,y.cor.pos=0.8,
                                                                   mar.vect = c(5,5,1,5))
  
  original.distances.1000 = original.distances[diversity[selected_groups] > 1000,diversity[selected_groups] > 1000]
  x = original.distances.1000#[lower.tri(original.distances.1000)]
  distances.from.axes.1000 = as.matrix(dist(NormalizedVI_pcoa[[3]][diversity[selected_groups] > 1000,1:2]))
  y =  distances.from.axes.1000#[lower.tri(distances.from.axes.1000)]
  plot.distances.1000.first.2.axes.vs.original.distances = cor.plot(x = x,
                                                                    y = y,
                                                                    x.lab = "Original distances between >1000-OTU groups     ",
                                                                    y.lab = "Distances based on first 2 PCoA axes",
                                                                    x.lab.hjust = 0.5,
                                                                    y.lab.hjust = 0.5,
                                                                    x.log = F,
                                                                    y.log = F,
                                                                    fit = T,
                                                                    fit.display = "pearson.p",
                                                                    x.cor.pos=0.2,y.cor.pos=1,
                                                                    dis=T,
                                                                    mar.vect = c(5,8,1,5)) + 
    ylim(range(as.matrix(dist(NormalizedVI_pcoa[[3]][,1:2])))) +
    xlim(range(original.distances,na.rm=T))
  
  x = as.matrix(dist(NormalizedVI_pcoa[[3]][diversity[selected_groups] > 1000,]))
  y = as.matrix(dist(NormalizedVI_pcoa[[3]][diversity[selected_groups] > 1000,1:2]))
  plot.distances.1000.first.2.axes.vs.all.axes = cor.plot(x = x,
                                                          y = y,
                                                          x.lab = "Distances between >1000-OTU groups\n based on all axes",
                                                          y.lab = "Distances between >1000-OTU groups\n based on the first 2 PCoA axes",
                                                          x.lab.hjust = 0.5,
                                                          y.lab.hjust = 0.5,
                                                          x.log = F,
                                                          y.log = F,
                                                          fit = F,
                                                          fit.display = "pearson.spearman",
                                                          x.cor.pos=0.2,y.cor.pos=0.8,
                                                          dis=T,
                                                          mar.vect = c(5,5,1,5))
  
  # original.distances.1000 = original.distances[diversity[selected_groups] > 1000,diversity[selected_groups] > 1000]
  # x = original.distances.1000[lower.tri(original.distances.1000)]
  # distances.from.axes.1000 = as.matrix(dist(NormalizedVI_pcoa[[3]][diversity[selected_groups] > 1000,1:9]))
  # y =  distances.from.axes.1000[lower.tri(distances.from.axes.1000)]
  # plot.distances.1000.first.9.axes.vs.original.distances = cor.plot(x = x,
  #                                                              y = y,
  #                                                              x.lab = "Original VI distances between >1000-OTU groups",
  #                                                              y.lab = "Distances based on first 9 PCoA axes",
  #                                                              y.lab.hjust = 0.5,
  #                                                              x.log = F,
  #                                                              y.log = F,
  #                                                              fit = T,
  #                                                              fit.display = "pearson.spearman",
  #                                                              x.cor.pos=0.2,y.cor.pos=0.8,
  #                                                              mar.vect = c(5,5,1,5))
  
  pdf(paste0(figure_folder,"/FigS_VI.distances.computed.from.PCoA.axes_vs_original.VI.distances.pdf"),
      width=7.5*2.2/1.2/2*2,height=12/3*4/1.2/2*3,onefile=F)
  g1 = ggplotGrob(plot.distances.first.2.axes.vs.original.distances)
  g2 = ggplotGrob(plot.distances.200.first.2.axes.vs.original.distances)
  g3 = ggplotGrob(plot.distances.500.first.2.axes.vs.original.distances)
  g4 = ggplotGrob(plot.distances.1000.first.2.axes.vs.original.distances)
  g5 = ggplotGrob(plot.cor.pearson.distances.vs.nb.axes)
  g6 = ggplotGrob(plot.cor.pearson.distances.vs.explained.var)
  grid.arrange(grobs = list(g1, g2, g3, g4, g5, g6), nrow = 3, ncol = 2)
  dev.off()
  
  ############ Axes correlation to explanatory variables
  
  descriptors = list()
  for (i_case in 1:2)
  {
    descriptors[[i_case]] = cbind(`Spatial.autocorrelation` = I_square.observed_w.mean[,i_case],
                                  `Charac.spatial.scale` = charac_scale[,i_case],
                                  `Latitudinal.symmetry.I` = lat_I[,i_case]/I_square.observed_w.mean[,i_case],
                                  `Basin.scale.similarity.I` = basin_I[,i_case]/I_square.observed_w.mean[,i_case],
                                  # `OTU.richness.log10` = log10(as.vector(diversity)),
                                  # `Body.size.log10` = log10(size_relativeAbund))
                                  `mean.Shannon.entropy` = if (i_case == 1) shannon_SUR else shannon_DCM,
                                  `total.Shannon` = if (i_case == 1) shannon_total_SUR else shannon_total_DCM)
  }
  
  all.PCoA.axes.pval = list()
  # All axes vs. each variable:
  for (i_case in 1:2)
  {
    all.PCoA.axes.pval[[i_case]] = matrix(nrow = ncol(descriptors[[i_case]]), ncol = 2, data = 0,
                                          dimnames = list(colnames(descriptors[[i_case]]),c("Adj.R2","p.val")))
    for (i in 1:ncol(descriptors[[i_case]]))
    {
      RDA = rda(NormalizedVI_pcoa[[3]] ~ descriptors[[i_case]][selected_groups,i], na.action = na.exclude)
      all.PCoA.axes.pval[[i_case]][i,1] = as.numeric(RsquareAdj(RDA)[2])
      all.PCoA.axes.pval[[i_case]][i,2] = as.numeric(anova(RDA)$'Pr(>F)'[1])
    }
  }
  
  # Each axis vs all variables (separately for Sur and DCM):
  PCoA.axis.pval = list()
  for (i_case in 1:2)
  {
    PCoA.axis.pval[[i_case]] = matrix(nrow = ncol(NormalizedVI_pcoa[[3]]), ncol = 3, data = 0,
                                      dimnames = list(colnames(NormalizedVI_pcoa[[3]]),c("R2","Adj.R2","p.val")))
    for (i in 1:ncol(NormalizedVI_pcoa[[3]]))
    {
      LM = lm(as.vector(NormalizedVI_pcoa[[3]][,i]) ~ descriptors[[i_case]][selected_groups,], na.action = na.exclude)
      PCoA.axis.pval[[i_case]][i,1] = summary(LM)$r.squared
      PCoA.axis.pval[[i_case]][i,2] = summary(LM)$adj.r.squared #as.numeric(RsquareAdj(RDA)[2])
      PCoA.axis.pval[[i_case]][i,3] = as.numeric(anova(LM)$'Pr(>F)'[1]) #as.numeric(anova(RDA)$'Pr(>F)'[1])
    }
  }
  
  # Each axis vs all variables:
  PCoA.axis.pval.SUR.DCM = matrix(nrow = ncol(NormalizedVI_pcoa[[3]]), ncol = 3, data = 0,
                                  dimnames = list(colnames(NormalizedVI_pcoa[[3]]),c("R2","Adj.R2","p.val")))
  for (i in 1:ncol(NormalizedVI_pcoa[[3]]))
  {
    LM = lm(as.vector(NormalizedVI_pcoa[[3]][,i]) ~ cbind(descriptors[[1]][selected_groups,],descriptors[[2]][selected_groups,]), na.action = na.exclude)
    PCoA.axis.pval.SUR.DCM[i,1] = summary(LM)$r.squared
    PCoA.axis.pval.SUR.DCM[i,2] = summary(LM)$adj.r.squared #as.numeric(RsquareAdj(RDA)[2])
    PCoA.axis.pval.SUR.DCM[i,3] = as.numeric(anova(LM)$'Pr(>F)'[1]) #as.numeric(anova(RDA)$'Pr(>F)'[1])
  }
  
  x = 1:ncol(NormalizedVI_pcoa[[3]])
  # y = PCoA.axis.pval[[1]][,2]
  y = PCoA.axis.pval.SUR.DCM[,2]
  # col.factor = rep("signif",length(x))
  # col.factor[PCoA.axis.pval[[1]][,3]>0.05] = "non-signif"
  col.factor = c(rep("blue",2),rep("black",length(x)-2))
  plot.var.explained.by.descriptors.vs.explained.var = cor.plot(x = x,
                                                                y = y,
                                                                col.factor = col.factor,
                                                                # col.vect = c("black","red"),
                                                                col.vect = c("black","blue"),
                                                                leg.name = NULL,
                                                                x.lab = "Rank of PCoA axis",
                                                                # y.lab = "Variance along PCoA axis linearly explained\n by biogeographic and ecological descriptors",
                                                                y.lab = "Variance along PCoA axis linearly explained\n by biogeographic descriptors",
                                                                y.lab.hjust = 0.5,
                                                                x.log = F,
                                                                y.log = F,
                                                                fit = F,
                                                                # fit.display = "pearson.spearman",
                                                                mar.vect = c(5,8,1,5)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylim(min(PCoA.axis.pval.SUR.DCM[,2]),1)
  
  pdf(paste0(figure_folder,"/FigS_PCoA.axes_explained.variance_normalized.descriptors.SUR.DCM.pdf"),
  # pdf(paste0(figure_folder,"/FigS_PCoA.axes_explained.variance_distances.all.axes.pdf"),
      width=7.5*2.2/1.2/2*2,height=12/3*4/1.2/2*2,onefile=F)
  g1 = ggplotGrob(plot.var.explained.by.descriptors.vs.explained.var)
  g2 = ggplotGrob(plot.cor.pearson.distances.vs.explained.var)
  g3 = ggplotGrob(plot.distances.first.2.axes.vs.original.distances)
  # g3 = ggplotGrob(plot.distances.first.2.axes.vs.all.axes)
  g4 = ggplotGrob(plot.distances.1000.first.2.axes.vs.original.distances)
  # g4 = ggplotGrob(plot.distances.1000.first.2.axes.vs.all.axes)
  grid.arrange(grobs = list(g1, g2, g3, g4), nrow = 2, ncol = 2)
  dev.off()
  
  ############ First PCoA axis only:
  
  x = original.distances[lower.tri(original.distances)]
  distances.from.axes = as.matrix(dist(NormalizedVI_pcoa[[3]][,1]))
  y =  distances.from.axes[lower.tri(distances.from.axes)]
  plot.distances.first.axis.vs.original.distances = cor.plot(x = x,
                                                             y = y,
                                                             x.lab = "Original VI distances between groups",
                                                             y.lab = "Distances based on the first PCoA axis",
                                                             y.lab.hjust = 0.5,
                                                             x.log = F,
                                                             y.log = F,
                                                             fit = T,
                                                             fit.display = "pearson.spearman",
                                                             x.cor.pos=0.2,y.cor.pos=0.8,
                                                             mar.vect = c(5,5,1,5))
  
  original.distances.200 = original.distances[diversity[selected_groups] > 200,diversity[selected_groups] > 200]
  x = original.distances.200[lower.tri(original.distances.200)]
  distances.from.axes.200 = as.matrix(dist(NormalizedVI_pcoa[[3]][diversity[selected_groups] > 200,1]))
  y =  distances.from.axes.200[lower.tri(distances.from.axes.200)]
  plot.distances.200.first.axis.vs.original.distances = cor.plot(x = x,
                                                                 y = y,
                                                                 x.lab = "Original VI distances between >200-OTU groups",
                                                                 y.lab = "Distances based on the first PCoA axis",
                                                                 y.lab.hjust = 0.5,
                                                                 x.log = F,
                                                                 y.log = F,
                                                                 fit = T,
                                                                 fit.display = "pearson.spearman",
                                                                 x.cor.pos=0.2,y.cor.pos=0.8,
                                                                 mar.vect = c(5,5,1,5))
  
  original.distances.500 = original.distances[diversity[selected_groups] > 500,diversity[selected_groups] > 500]
  x = original.distances.500[lower.tri(original.distances.500)]
  distances.from.axes.500 = as.matrix(dist(NormalizedVI_pcoa[[3]][diversity[selected_groups] > 500,1]))
  y =  distances.from.axes.500[lower.tri(distances.from.axes.500)]
  plot.distances.500.first.axis.vs.original.distances = cor.plot(x = x,
                                                                 y = y,
                                                                 x.lab = "Original VI distances between >500-OTU groups",
                                                                 y.lab = "Distances based on the first PCoA axis",
                                                                 y.lab.hjust = 0.5,
                                                                 x.log = F,
                                                                 y.log = F,
                                                                 fit = T,
                                                                 fit.display = "pearson.spearman",
                                                                 x.cor.pos=0.2,y.cor.pos=0.8,
                                                                 mar.vect = c(5,5,1,5))
  
  original.distances.1000 = original.distances[diversity[selected_groups] > 1000,diversity[selected_groups] > 1000]
  x = original.distances.1000[lower.tri(original.distances.1000)]
  distances.from.axes.1000 = as.matrix(dist(NormalizedVI_pcoa[[3]][diversity[selected_groups] > 1000,1]))
  y =  distances.from.axes.1000[lower.tri(distances.from.axes.1000)]
  plot.distances.1000.first.axis.vs.original.distances = cor.plot(x = x,
                                                                  y = y,
                                                                  x.lab = "Original VI distances between >1000-OTU groups",
                                                                  y.lab = "Distances based on the first PCoA axis",
                                                                  y.lab.hjust = 0.5,
                                                                  x.log = F,
                                                                  y.log = F,
                                                                  fit = T,
                                                                  fit.display = "pearson.spearman",
                                                                  x.cor.pos=0.2,y.cor.pos=0.8,
                                                                  mar.vect = c(5,5,1,5))
  
  pdf(paste0(figure_folder,"/FigS_VI.distances.computed.from.first.PCoA.axis_vs_original.VI.distances.pdf"),
      width=7.5*2.2/1.2/2*2,height=12/3*4/1.2/2*2,onefile=F)
  g1 = ggplotGrob(plot.distances.first.axis.vs.original.distances)
  g2 = ggplotGrob(plot.distances.200.first.axis.vs.original.distances)
  g3 = ggplotGrob(plot.distances.500.first.axis.vs.original.distances)
  g4 = ggplotGrob(plot.distances.1000.first.axis.vs.original.distances)
  grid.arrange(grobs = list(g1, g2, g3, g4), nrow = 2, ncol = 2)
  dev.off()
  
  ############ Second PCoA axis only:
  
  x = original.distances[lower.tri(original.distances)]
  distances.from.axes = as.matrix(dist(NormalizedVI_pcoa[[3]][,2]))
  y =  distances.from.axes[lower.tri(distances.from.axes)]
  plot.distances.second.axis.vs.original.distances = cor.plot(x = x,
                                                              y = y,
                                                              x.lab = "Original VI distances between groups",
                                                              y.lab = "Distances based on the second PCoA axis",
                                                              y.lab.hjust = 0.5,
                                                              x.log = F,
                                                              y.log = F,
                                                              fit = T,
                                                              fit.display = "pearson.spearman",
                                                              x.cor.pos=0.2,y.cor.pos=0.8,
                                                              mar.vect = c(5,5,1,5))
  
  original.distances.200 = original.distances[diversity[selected_groups] > 200,diversity[selected_groups] > 200]
  x = original.distances.200[lower.tri(original.distances.200)]
  distances.from.axes.200 = as.matrix(dist(NormalizedVI_pcoa[[3]][diversity[selected_groups] > 200,2]))
  y =  distances.from.axes.200[lower.tri(distances.from.axes.200)]
  plot.distances.200.second.axis.vs.original.distances = cor.plot(x = x,
                                                                  y = y,
                                                                  x.lab = "Original VI distances between >200-OTU groups",
                                                                  y.lab = "Distances based on the second PCoA axis",
                                                                  y.lab.hjust = 0.5,
                                                                  x.log = F,
                                                                  y.log = F,
                                                                  fit = T,
                                                                  fit.display = "pearson.spearman",
                                                                  x.cor.pos=0.2,y.cor.pos=0.8,
                                                                  mar.vect = c(5,5,1,5))
  
  original.distances.500 = original.distances[diversity[selected_groups] > 500,diversity[selected_groups] > 500]
  x = original.distances.500[lower.tri(original.distances.500)]
  distances.from.axes.500 = as.matrix(dist(NormalizedVI_pcoa[[3]][diversity[selected_groups] > 500,2]))
  y =  distances.from.axes.500[lower.tri(distances.from.axes.500)]
  plot.distances.500.second.axis.vs.original.distances = cor.plot(x = x,
                                                                  y = y,
                                                                  x.lab = "Original VI distances between >500-OTU groups",
                                                                  y.lab = "Distances based on the second PCoA axis",
                                                                  y.lab.hjust = 0.5,
                                                                  x.log = F,
                                                                  y.log = F,
                                                                  fit = T,
                                                                  fit.display = "pearson.spearman",
                                                                  x.cor.pos=0.2,y.cor.pos=0.8,
                                                                  mar.vect = c(5,5,1,5))
  
  original.distances.1000 = original.distances[diversity[selected_groups] > 1000,diversity[selected_groups] > 1000]
  x = original.distances.1000[lower.tri(original.distances.1000)]
  distances.from.axes.1000 = as.matrix(dist(NormalizedVI_pcoa[[3]][diversity[selected_groups] > 1000,2]))
  y =  distances.from.axes.1000[lower.tri(distances.from.axes.1000)]
  plot.distances.1000.second.axis.vs.original.distances = cor.plot(x = x,
                                                                   y = y,
                                                                   x.lab = "Original VI distances between >1000-OTU groups",
                                                                   y.lab = "Distances based on the second PCoA axis",
                                                                   y.lab.hjust = 0.5,
                                                                   x.log = F,
                                                                   y.log = F,
                                                                   fit = T,
                                                                   fit.display = "pearson.spearman",
                                                                   x.cor.pos=0.2,y.cor.pos=0.8,
                                                                   mar.vect = c(5,5,1,5))
  
  pdf(paste0(figure_folder,"/FigS_VI.distances.computed.from.second.PCoA.axis_vs_original.VI.distances.pdf"),
      width=7.5*2.2/1.2/2*2,height=12/3*4/1.2/2*2,onefile=F)
  g1 = ggplotGrob(plot.distances.second.axis.vs.original.distances)
  g2 = ggplotGrob(plot.distances.200.second.axis.vs.original.distances)
  g3 = ggplotGrob(plot.distances.500.second.axis.vs.original.distances)
  g4 = ggplotGrob(plot.distances.1000.second.axis.vs.original.distances)
  grid.arrange(grobs = list(g1, g2, g3, g4), nrow = 2, ncol = 2)
  dev.off()
  
  ############ First 9 PCoA axes:
  
  x = original.distances[lower.tri(original.distances)]
  distances.from.axes = as.matrix(dist(NormalizedVI_pcoa[[3]][,1:9]))
  y =  distances.from.axes[lower.tri(distances.from.axes)]
  plot.distances.first.9.axes.vs.original.distances = cor.plot(x = x,
                                                               y = y,
                                                               x.lab = "Original VI distances between groups",
                                                               y.lab = "Distances based on the first 9 PCoA axes",
                                                               y.lab.hjust = 0.5,
                                                               x.log = F,
                                                               y.log = F,
                                                               fit = T,
                                                               fit.display = "pearson.spearman",
                                                               x.cor.pos=0.2,y.cor.pos=0.8,
                                                               mar.vect = c(5,5,1,5))
  
  original.distances.200 = original.distances[diversity[selected_groups] > 200,diversity[selected_groups] > 200]
  x = original.distances.200[lower.tri(original.distances.200)]
  distances.from.axes.200 = as.matrix(dist(NormalizedVI_pcoa[[3]][diversity[selected_groups] > 200,1:9]))
  y =  distances.from.axes.200[lower.tri(distances.from.axes.200)]
  plot.distances.200.first.9.axes.vs.original.distances = cor.plot(x = x,
                                                                   y = y,
                                                                   x.lab = "Original VI distances between >200-OTU groups",
                                                                   y.lab = "Distances based on the first 9 PCoA axes",
                                                                   y.lab.hjust = 0.5,
                                                                   x.log = F,
                                                                   y.log = F,
                                                                   fit = T,
                                                                   fit.display = "pearson.spearman",
                                                                   x.cor.pos=0.2,y.cor.pos=0.8,
                                                                   mar.vect = c(5,5,1,5))
  
  original.distances.500 = original.distances[diversity[selected_groups] > 500,diversity[selected_groups] > 500]
  x = original.distances.500[lower.tri(original.distances.500)]
  distances.from.axes.500 = as.matrix(dist(NormalizedVI_pcoa[[3]][diversity[selected_groups] > 500,1:9]))
  y =  distances.from.axes.500[lower.tri(distances.from.axes.500)]
  plot.distances.500.first.9.axes.vs.original.distances = cor.plot(x = x,
                                                                   y = y,
                                                                   x.lab = "Original VI distances between >500-OTU groups",
                                                                   y.lab = "Distances based on the first 9 PCoA axes",
                                                                   y.lab.hjust = 0.5,
                                                                   x.log = F,
                                                                   y.log = F,
                                                                   fit = T,
                                                                   fit.display = "pearson.spearman",
                                                                   x.cor.pos=0.2,y.cor.pos=0.8,
                                                                   mar.vect = c(5,5,1,5))
  
  original.distances.1000 = original.distances[diversity[selected_groups] > 1000,diversity[selected_groups] > 1000]
  x = original.distances.1000[lower.tri(original.distances.1000)]
  distances.from.axes.1000 = as.matrix(dist(NormalizedVI_pcoa[[3]][diversity[selected_groups] > 1000,1:9]))
  y =  distances.from.axes.1000[lower.tri(distances.from.axes.1000)]
  plot.distances.1000.first.9.axes.vs.original.distances = cor.plot(x = x,
                                                                    y = y,
                                                                    x.lab = "Original VI distances between >1000-OTU groups",
                                                                    y.lab = "Distances based on the first 9 PCoA axes",
                                                                    y.lab.hjust = 0.5,
                                                                    x.log = F,
                                                                    y.log = F,
                                                                    fit = T,
                                                                    fit.display = "pearson.spearman",
                                                                    x.cor.pos=0.2,y.cor.pos=0.8,
                                                                    mar.vect = c(5,5,1,5))
  
  pdf(paste0(figure_folder,"/FigS_VI.distances.computed.from.first.9.PCoA.axes_vs_original.VI.distances.pdf"),
      width=7.5*2.2/1.2/2*2,height=12/3*4/1.2/2*2,onefile=F)
  g1 = ggplotGrob(plot.distances.first.9.axes.vs.original.distances)
  g2 = ggplotGrob(plot.distances.200.first.9.axes.vs.original.distances)
  g3 = ggplotGrob(plot.distances.500.first.9.axes.vs.original.distances)
  g4 = ggplotGrob(plot.distances.1000.first.9.axes.vs.original.distances)
  grid.arrange(grobs = list(g1, g2, g3, g4), nrow = 2, ncol = 2)
  dev.off()
}

# Original inter-group distances vs. log body size and log diversity distances
{
  original.distances = Normalized_VI[[3]][selected_groups,selected_groups]
  log.body.size.distances = as.matrix(dist(log(size_relativeAbund[selected_groups])))
  log.diversity.distances = as.matrix(dist(log(as.vector(diversity)[selected_groups])))
  body.size.distances = as.matrix(dist(size_relativeAbund[selected_groups]))
  diversity.distances = as.matrix(dist(as.vector(diversity)[selected_groups]))
  plot(original.distances[diversity[selected_groups]>1000,diversity[selected_groups]>1000],
       log.diversity.distances[diversity[selected_groups]>1000,diversity[selected_groups]>1000])
  plot(original.distances[diversity[selected_groups]>1000,diversity[selected_groups]>1000],
       log.body.size.distances[diversity[selected_groups]>1000,diversity[selected_groups]>1000])
  plot(original.distances[diversity[selected_groups]>1000,diversity[selected_groups]>1000],
       diversity.distances[diversity[selected_groups]>1000,diversity[selected_groups]>1000])
  plot(original.distances[diversity[selected_groups]>1000,diversity[selected_groups]>1000],
       body.size.distances[diversity[selected_groups]>1000,diversity[selected_groups]>1000])
}

################################################ PCoA1 biogeographic descriptors:

# Fig. S Simpson-Shannon vs. Moran I across groups
{
  div_threshold = 100
  
  plot.shannon.vs.Moran.I = cor.plot(x = I_square.observed_w.mean[selected_groups,1],
                                     y = shannon_SUR[selected_groups],
                                     excluded.points = which(taxo_groups[selected_groups] == "Porifera"),
                                     x.lab = "Short-distance spatial autocorrelation",
                                     y.lab = "Mean Shannon diversity across stations",
                                     y.lab.hjust = 0.5,
                                     x.log = F,
                                     y.log = F,
                                     fit = T,
                                     x.cor.pos = 0.1,
                                     y.cor.pos = 0.1,
                                     mar.vect = c(5,5,1,5))
  
  plot.simpson.vs.Moran.I = cor.plot(x = I_square.observed_w.mean[selected_groups,1],
                                     y = simpson_SUR[selected_groups],
                                     excluded.points = which(taxo_groups[selected_groups] == "Porifera"),
                                     x.lab = "Short-distance spatial autocorrelation",
                                     y.lab = "Mean Simpson diversity across stations",
                                     y.lab.hjust = 0.5,
                                     x.log = F,
                                     y.log = F,
                                     fit = T,
                                     x.cor.pos = 0.1,
                                     y.cor.pos = 0.1,
                                     mar.vect = c(5,5,1,5))
  
  plot.invsimpson.vs.Moran.I = cor.plot(x = I_square.observed_w.mean[selected_groups,1],
                                     y = invsimpson_SUR[selected_groups],
                                     excluded.points = which(taxo_groups[selected_groups] == "Porifera"),
                                     x.lab = "Short-distance spatial autocorrelation",
                                     y.lab = "Mean inverse Simpson diversity across stations",
                                     y.lab.hjust = 0.5,
                                     x.log = F,
                                     y.log = F,
                                     fit = T,
                                     x.cor.pos = 0.1,
                                     y.cor.pos = 0.1,
                                     mar.vect = c(5,5,1,5))
  
  plot.shannon.total.vs.Moran.I = cor.plot(x = I_square.observed_w.mean[selected_groups,1],
                                     y = shannon_total_SUR[selected_groups],
                                     x.lab = "Short-distance spatial autocorrelation",
                                     y.lab = "Shannon diversity of pooled stations",
                                     y.lab.hjust = 0.5,
                                     x.log = F,
                                     y.log = F,
                                     fit = T,
                                     x.cor.pos = 0.1,
                                     y.cor.pos = 0.1,
                                     mar.vect = c(5,5,1,5))
  
  plot.simpson.total.vs.Moran.I = cor.plot(x = I_square.observed_w.mean[selected_groups,1],
                                     y = simpson_total_SUR[selected_groups],
                                     x.lab = "Short-distance spatial autocorrelation",
                                     y.lab = "Simpson diversity of pooled stations",
                                     y.lab.hjust = 0.5,
                                     x.log = F,
                                     y.log = F,
                                     fit = T,
                                     x.cor.pos = 0.1,
                                     y.cor.pos = 0.1,
                                     mar.vect = c(5,5,1,5))
  
  plot.invsimpson.total.vs.Moran.I = cor.plot(x = I_square.observed_w.mean[selected_groups,1],
                                           y = invsimpson_total_SUR[selected_groups],
                                           x.lab = "Short-distance spatial autocorrelation",
                                           y.lab = "Inverse Simpson diversity of pooled stations",
                                           y.lab.hjust = 0.5,
                                           x.log = F,
                                           y.log = F,
                                           fit = T,
                                           x.cor.pos = 0.1,
                                           y.cor.pos = 0.1,
                                           mar.vect = c(5,5,1,5))
  
  pdf(paste0(figure_folder,"/FigS_Shannon-Simpson.vs.Moran.I",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1_(1).pdf"),
      width=7.5*2.2/1.2/2*3,height=12/3*4/1.2,onefile=F)
  g1 = ggplotGrob(plot.shannon.vs.Moran.I)
  g2 = ggplotGrob(plot.simpson.vs.Moran.I)
  g3 = ggplotGrob(plot.invsimpson.vs.Moran.I)
  g4 = ggplotGrob(plot.shannon.total.vs.Moran.I)
  g5 = ggplotGrob(plot.simpson.total.vs.Moran.I)
  g6 = ggplotGrob(plot.invsimpson.total.vs.Moran.I)
  grid.arrange(grobs = list(g1, g2, g3, g4, g5, g6), nrow = 2, ncol = 3)
  dev.off()
}

# Fig. S Simpson-Shannon vs. diversity across groups
{
  div_threshold = 100
  
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera"]
  y = shannon_SUR[selected_groups & taxo_groups != "Porifera"]
  x.2000 = x[x<2000]
  y.2000 = y[x<2000]
  x.cor.pos = 0.8
  y.cor.pos = 0.8
  cor.test = cor.test(log10(x.2000),y.2000,na.rm=T)
  plot.shannon.vs.diversity = cor.plot(x = x,
                                     y = y,
                                     x.lab = "OTU richness",
                                     y.lab = "Mean Shannon diversity across stations",
                                     y.lab.hjust = 0.5,
                                     x.log = T,
                                     y.log = F,
                                     fit = F,
                                     mar.vect = c(5,5,1,5)) +
    geom_smooth(data = data.frame(x = x.2000, y = y.2000),aes(x,y),method='lm',col="black",se=T) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    annotate(geom="text", 
             x=(max(x)/min(x))^x.cor.pos*min(x),
             y=y.cor.pos*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=8)
    
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera"]
  y = simpson_SUR[selected_groups & taxo_groups != "Porifera"]
  x.2000 = x[x<2000]
  y.2000 = y[x<2000]
  x.cor.pos = 0.8
  y.cor.pos = 0.8
  cor.test = cor.test(log10(x.2000),y.2000,na.rm=T)
  plot.simpson.vs.diversity = cor.plot(x = x,
                                       y = y,
                                       x.lab = "OTU richness",
                                       y.lab = "Mean Simpson diversity across stations",
                                       y.lab.hjust = 0.5,
                                       x.log = T,
                                       y.log = F,
                                       fit = F,
                                       mar.vect = c(5,5,1,5)) +
    geom_smooth(data = data.frame(x = x.2000, y = y.2000),aes(x,y),method='lm',col="black",se=T) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    annotate(geom="text", 
             x=(max(x)/min(x))^x.cor.pos*min(x),
             y=y.cor.pos*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=8)
  
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera"]
  y = invsimpson_SUR[selected_groups & taxo_groups != "Porifera"]
  x.2000 = x[x<2000]
  y.2000 = y[x<2000]
  x.cor.pos = 0.8
  y.cor.pos = 0.8
  cor.test = cor.test(log10(x.2000),y.2000,na.rm=T)
  plot.invsimpson.vs.diversity = cor.plot(x = x,
                                       y = y,
                                       x.lab = "OTU richness",
                                       y.lab = "Mean inverse Simpson diversity across stations",
                                       y.lab.hjust = 0.5,
                                       x.log = T,
                                       y.log = F,
                                       fit = F,
                                       mar.vect = c(5,5,1,5)) +
    geom_smooth(data = data.frame(x = x.2000, y = y.2000),aes(x,y),method='lm',col="black",se=T) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    annotate(geom="text", 
             x=(max(x)/min(x))^x.cor.pos*min(x),
             y=y.cor.pos*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=8)
  
  x = as.vector(diversity)[selected_groups]
  y = shannon_total_SUR[selected_groups]
  x.2000 = x[x<2000]
  y.2000 = y[x<2000]
  x.cor.pos = 0.8
  y.cor.pos = 0.8
  cor.test = cor.test(log10(x.2000),y.2000,na.rm=T)
  plot.shannon.total.vs.diversity = cor.plot(x = x,
                                          y = y,
                                          x.lab = "OTU richness",
                                          y.lab = "Shannon diversity of pooled stations",
                                          y.lab.hjust = 0.5,
                                          x.log = T,
                                          y.log = F,
                                          fit = F,
                                          mar.vect = c(5,5,1,5)) +
    geom_smooth(data = data.frame(x = x.2000, y = y.2000),aes(x,y),method='lm',col="black",se=T) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    annotate(geom="text", 
             x=(max(x)/min(x))^x.cor.pos*min(x),
             y=y.cor.pos*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=8)
  
  x = as.vector(diversity)[selected_groups]
  y = simpson_total_SUR[selected_groups]
  x.2000 = x[x<2000]
  y.2000 = y[x<2000]
  x.cor.pos = 0.8
  y.cor.pos = 0.8
  cor.test = cor.test(log10(x.2000),y.2000,na.rm=T)
  plot.simpson.total.vs.diversity = cor.plot(x = x,
                                             y = y,
                                             x.lab = "OTU richness",
                                             y.lab = "Simpson diversity of pooled stations",
                                             y.lab.hjust = 0.5,
                                             x.log = T,
                                             y.log = F,
                                             fit = F,
                                             mar.vect = c(5,5,1,5)) +
    geom_smooth(data = data.frame(x = x.2000, y = y.2000),aes(x,y),method='lm',col="black",se=T) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    annotate(geom="text", 
             x=(max(x)/min(x))^x.cor.pos*min(x),
             y=y.cor.pos*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=8)
  
  x = as.vector(diversity)[selected_groups]
  y = invsimpson_total_SUR[selected_groups]
  x.2000 = x[x<2000]
  y.2000 = y[x<2000]
  x.cor.pos = 0.8
  y.cor.pos = 0.8
  cor.test = cor.test(log10(x.2000),y.2000,na.rm=T)
  plot.invsimpson.total.vs.diversity = cor.plot(x = x,
                                             y = y,
                                             x.lab = "OTU richness",
                                             y.lab = "Inverse Simpson diversity of pooled stations",
                                             y.lab.hjust = 0.5,
                                             x.log = T,
                                             y.log = F,
                                             fit = F,
                                             mar.vect = c(5,5,1,5)) +
    geom_smooth(data = data.frame(x = x.2000, y = y.2000),aes(x,y),method='lm',col="black",se=T) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    annotate(geom="text", 
             x=(max(x)/min(x))^x.cor.pos*min(x),
             y=y.cor.pos*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=8)
  
  pdf(paste0(figure_folder,"/FigS_Shannon-Simpson.vs.diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1_(1)_noPorifera.pdf"),
      width=7.5*2.2/1.2/2*3,height=12/3*4/1.2,onefile=F)
  g1 = ggplotGrob(plot.shannon.vs.diversity)
  g2 = ggplotGrob(plot.simpson.vs.diversity)
  g3 = ggplotGrob(plot.invsimpson.vs.diversity)
  g4 = ggplotGrob(plot.shannon.total.vs.diversity)
  g5 = ggplotGrob(plot.simpson.total.vs.diversity)
  g6 = ggplotGrob(plot.invsimpson.total.vs.diversity)
  grid.arrange(grobs = list(g1, g2, g3, g4, g5, g6), nrow = 2, ncol = 3)
  dev.off()
}

# Fig. S total explained variance vs. Simpson-Shannon across groups
{
  div_threshold = 100
  
  plot.tot.var.vs.shannon = cor.plot(x = shannon_SUR[selected_groups & !is.na(colSums(varpart.env.spatial[[1]]))],
                                     y = colSums(varpart.env.spatial[[1]])[selected_groups & !is.na(colSums(varpart.env.spatial[[1]]))],
                                     x.lab = "Mean Shannon diversity across stations",
                                     y.lab = "Total variance explained\n by connectivity and local environment",
                                     y.lab.hjust = 0.5,
                                     x.log = F,
                                     y.log = F,
                                     fit = T,
                                     x.cor.pos = 0.8,
                                     y.cor.pos = 0.8,
                                     mar.vect = c(5,5,1,5))
  
  plot.tot.var.vs.simpson = cor.plot(x = simpson_SUR[selected_groups & !is.na(colSums(varpart.env.spatial[[1]]))],
                                     y = colSums(varpart.env.spatial[[1]])[selected_groups & !is.na(colSums(varpart.env.spatial[[1]]))],
                                     x.lab = "Mean Simpson diversity across stations",
                                     y.lab = "Total variance explained\n by connectivity and local environment",
                                     y.lab.hjust = 0.5,
                                     x.log = F,
                                     y.log = F,
                                     fit = T,
                                     x.cor.pos = 0.8,
                                     y.cor.pos = 0.8,
                                     mar.vect = c(5,5,1,5))
  
  plot.tot.var.vs.invsimpson = cor.plot(x = invsimpson_SUR[selected_groups & !is.na(colSums(varpart.env.spatial[[1]]))],
                                        y = colSums(varpart.env.spatial[[1]])[selected_groups & !is.na(colSums(varpart.env.spatial[[1]]))],
                                        x.lab = "Mean inverse Simpson diversity across stations",
                                        y.lab = "Total variance explained\n by connectivity and local environment",
                                        y.lab.hjust = 0.5,
                                        x.log = F,
                                        y.log = F,
                                        fit = T,
                                        x.cor.pos = 0.8,
                                        y.cor.pos = 0.8,
                                        mar.vect = c(5,5,1,5))
  
  plot.tot.var.vs.shannon.total = cor.plot(x = shannon_total_SUR[selected_groups & !is.na(colSums(varpart.env.spatial[[1]]))],
                                           y = colSums(varpart.env.spatial[[1]])[selected_groups & !is.na(colSums(varpart.env.spatial[[1]]))],
                                           x.lab = "Shannon diversity of pooled stations",
                                           y.lab = "Total variance explained\n by connectivity and local environment",
                                           y.lab.hjust = 0.5,
                                           x.log = F,
                                           y.log = F,
                                           fit = T,
                                           x.cor.pos = 0.8,
                                           y.cor.pos = 0.8,
                                           mar.vect = c(5,5,1,5))
  
  plot.tot.var.vs.simpson.total = cor.plot(x = simpson_total_SUR[selected_groups & !is.na(colSums(varpart.env.spatial[[1]]))],
                                           y = colSums(varpart.env.spatial[[1]])[selected_groups & !is.na(colSums(varpart.env.spatial[[1]]))],
                                           x.lab = "Simpson diversity of pooled stations",
                                           y.lab = "Total variance explained\n by connectivity and local environment",
                                           y.lab.hjust = 0.5,
                                           x.log = F,
                                           y.log = F,
                                           fit = T,
                                           x.cor.pos = 0.8,
                                           y.cor.pos = 0.8,
                                           mar.vect = c(5,5,1,5))
  
  plot.tot.var.vs.invsimpson.total = cor.plot(x = invsimpson_total_SUR[selected_groups & !is.na(colSums(varpart.env.spatial[[1]]))],
                                              y = colSums(varpart.env.spatial[[1]])[selected_groups & !is.na(colSums(varpart.env.spatial[[1]]))],
                                              x.lab = "Inverse Simpson diversity of pooled stations",
                                              y.lab = "Total variance explained\n by connectivity and local environment",
                                              y.lab.hjust = 0.5,
                                              x.log = F,
                                              y.log = F,
                                              fit = T,
                                              x.cor.pos = 0.8,
                                              y.cor.pos = 0.8,
                                              mar.vect = c(5,5,1,5))
  
  pdf(paste0(figure_folder,"/FigS_tot.explained.var.vs.Shannon-Simpson",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2*3,height=12/3*4/1.2,onefile=F)
  g1 = ggplotGrob(plot.tot.var.vs.shannon)
  g2 = ggplotGrob(plot.tot.var.vs.simpson)
  g3 = ggplotGrob(plot.tot.var.vs.invsimpson)
  g4 = ggplotGrob(plot.tot.var.vs.shannon.total)
  g5 = ggplotGrob(plot.tot.var.vs.simpson.total)
  g6 = ggplotGrob(plot.tot.var.vs.invsimpson.total)
  grid.arrange(grobs = list(g1, g2, g3, g4, g5, g6), nrow = 2, ncol = 3)
  dev.off()
}

# Fig. S Simpson-Shannon vs. PCoA1 across groups
{
  div_threshold = 100
  
  plot.pcoa1.vs.shannon = cor.plot(x = NormalizedVI_pcoa[[3]][,1],
                                   y = shannon_SUR[selected_groups],
                                   excluded.points = which(taxo_groups[selected_groups] == "Porifera"),
                                   x.lab = "PCoA axis 1",
                                   y.lab = "Mean Shannon diversity across stations",
                                   y.lab.hjust = 0.5,
                                   x.log = F,
                                   y.log = F,
                                   fit = T,
                                   x.cor.pos = 0.8,
                                   y.cor.pos = 0.8,
                                   mar.vect = c(5,5,1,5))
  
  plot.pcoa1.vs.simpson = cor.plot(x = NormalizedVI_pcoa[[3]][,1],
                                   y = simpson_SUR[selected_groups],
                                   excluded.points = which(taxo_groups[selected_groups] == "Porifera"),
                                   x.lab = "PCoA axis 1",
                                   y.lab = "Mean Simpson diversity across stations",
                                   y.lab.hjust = 0.5,
                                   x.log = F,
                                   y.log = F,
                                   fit = T,
                                   x.cor.pos = 0.8,
                                   y.cor.pos = 0.8,
                                   mar.vect = c(5,5,1,5))
  
  plot.pcoa1.vs.invsimpson = cor.plot(x = NormalizedVI_pcoa[[3]][,1],
                                      y = invsimpson_SUR[selected_groups],
                                      excluded.points = which(taxo_groups[selected_groups] == "Porifera"),
                                      x.lab = "PCoA axis 1",
                                      y.lab = "Mean inverse Simpson diversity across stations",
                                      y.lab.hjust = 0.5,
                                      x.log = F,
                                      y.log = F,
                                      fit = T,
                                      x.cor.pos = 0.8,
                                      y.cor.pos = 0.8,
                                      mar.vect = c(5,5,1,5))
  
  plot.pcoa1.vs.shannon.total = cor.plot(x = NormalizedVI_pcoa[[3]][,1],
                                         y = shannon_total_SUR[selected_groups],
                                         x.lab = "PCoA axis 1",
                                         y.lab = "Shannon diversity of pooled stations",
                                         y.lab.hjust = 0.5,
                                         x.log = F,
                                         y.log = F,
                                         fit = T,
                                         x.cor.pos = 0.8,
                                         y.cor.pos = 0.8,
                                         mar.vect = c(5,5,1,5))
  
  plot.pcoa1.vs.simpson.total = cor.plot(x = NormalizedVI_pcoa[[3]][,1],
                                         y = simpson_total_SUR[selected_groups],
                                         x.lab = "PCoA axis 1",
                                         y.lab = "Simpson diversity of pooled stations",
                                         y.lab.hjust = 0.5,
                                         x.log = F,
                                         y.log = F,
                                         fit = T,
                                         x.cor.pos = 0.8,
                                         y.cor.pos = 0.8,
                                         mar.vect = c(5,5,1,5))
  
  plot.pcoa1.vs.invsimpson.total = cor.plot(x = NormalizedVI_pcoa[[3]][,1],
                                            y = invsimpson_total_SUR[selected_groups],
                                            x.lab = "PCoA axis 1",
                                            y.lab = "Inverse Simpson diversity of pooled stations",
                                            y.lab.hjust = 0.5,
                                            x.log = F,
                                            y.log = F,
                                            fit = T,
                                            x.cor.pos = 0.8,
                                            y.cor.pos = 0.8,
                                            mar.vect = c(5,5,1,5))
  
  pdf(paste0(figure_folder,"/FigS_PCoA1.vs.Shannon-Simpson",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2*3,height=12/3*4/1.2,onefile=F)
  g1 = ggplotGrob(plot.pcoa1.vs.shannon)
  g2 = ggplotGrob(plot.pcoa1.vs.simpson)
  g3 = ggplotGrob(plot.pcoa1.vs.invsimpson)
  g4 = ggplotGrob(plot.pcoa1.vs.shannon.total)
  g5 = ggplotGrob(plot.pcoa1.vs.simpson.total)
  g6 = ggplotGrob(plot.pcoa1.vs.invsimpson.total)
  grid.arrange(grobs = list(g1, g2, g3, g4, g5, g6), nrow = 2, ncol = 3)
  dev.off()
}

# Fig. S comparison mean Shannon across stations, Moran I and PCoA1 vs. tot. explained variance and diversity:
{
  div_threshold = 100
  
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera"]
  y = shannon_SUR[selected_groups & taxo_groups != "Porifera"]
  x.2000 = x[x<2000]
  y.2000 = y[x<2000]
  x.cor.pos = 0.8
  y.cor.pos = 0.8
  cor.test = cor.test(log10(x.2000),y.2000,na.rm=T)
  plot.shannon.vs.diversity = cor.plot(x = x,
                                       y = y,
                                       x.lab = "OTU richness",
                                       y.lab = "Mean Shannon entropy across stations",
                                       y.lab.hjust = 0.5,
                                       x.log = T,
                                       y.log = F,
                                       fit = F,
                                       mar.vect = c(5,5,1,5)) +
    geom_smooth(data = data.frame(x = x.2000, y = y.2000),aes(x,y),method='lm',col="black",se=T) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    annotate(geom="text", 
             x=(max(x)/min(x))^x.cor.pos*min(x),
             y=y.cor.pos*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=8)
  
  plot.tot.var.vs.shannon = cor.plot(x = shannon_SUR[selected_groups & !is.na(colSums(varpart.env.spatial[[1]]))],
                                     y = colSums(varpart.env.spatial[[1]])[selected_groups & !is.na(colSums(varpart.env.spatial[[1]]))],
                                     x.lab = "Mean Shannon entropy across stations",
                                     y.lab = "Total variance explained\n by connectivity and local environment",
                                     y.lab.hjust = 0.5,
                                     x.log = F,
                                     y.log = F,
                                     fit = T,
                                     x.cor.pos = 0.8,
                                     y.cor.pos = 0.8,
                                     mar.vect = c(5,5,1,5))
  
  x = as.vector(diversity)[selected_groups]
  y = I_square.observed_w.mean[selected_groups,1]/optimalK_prevalence.min.crossValid[selected_groups]
  x.2000 = x[x<2000]
  y.2000 = y[x<2000]
  x.cor.pos = 0.8
  y.cor.pos = 0.1
  cor.test = cor.test(log10(x.2000),y.2000,na.rm=T)
  plot.Moran.I.vs.diversity = cor.plot(x = x,
                                       y = y,
                                       x.lab = "OTU richness",
                                       y.lab = "Short-distance spatial autocorrelation",
                                       y.lab.hjust = 0.5,
                                       x.log = T,
                                       y.log = F,
                                       fit = F,
                                       mar.vect = c(5,5,1,5)) +
    geom_smooth(data = data.frame(x = x.2000, y = y.2000),aes(x,y),method='lm',col="black",se=T) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    annotate(geom="text", 
             x=(max(x)/min(x))^x.cor.pos*min(x),
             y=y.cor.pos*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=8)
  
  plot.tot.var.vs.Moran.I = cor.plot(x = I_square.observed_w.mean[selected_groups,1]/optimalK_prevalence.min.crossValid[selected_groups],
                                     y = colSums(varpart.env.spatial[[1]])[selected_groups],
                                     x.lab = "Short-distance spatial autocorrelation",
                                     y.lab = "Total variance explained\n by connectivity and local environment",
                                     y.lab.hjust = 0.5,
                                     x.log = F,
                                     y.log = F,
                                     fit = T,
                                     x.cor.pos = 0.8,
                                     y.cor.pos = 0.1,
                                     mar.vect = c(5,5,1,5))
  
  x = as.vector(diversity)[selected_groups]
  y = NormalizedVI_pcoa[[3]][,1]
  x.2000 = x[x<2000]
  y.2000 = y[x<2000]
  x.cor.pos = 0.8
  y.cor.pos = 0.1
  cor.test = cor.test(log10(x.2000),y.2000,na.rm=T)
  plot.PCoA1.vs.diversity = cor.plot(x = x,
                                     y = y,
                                     x.lab = "OTU richness",
                                     y.lab = "PCoA axis 1",
                                     y.lab.hjust = 0.5,
                                     x.log = T,
                                     y.log = F,
                                     fit = F,
                                     mar.vect = c(5,5,1,5)) +
    geom_smooth(data = data.frame(x = x.2000, y = y.2000),aes(x,y),method='lm',col="black",se=T) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    annotate(geom="text", 
             x=(max(x)/min(x))^x.cor.pos*min(x),
             y=y.cor.pos*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=8)
  
  plot.tot.var.vs.PCoA1 = cor.plot(x = NormalizedVI_pcoa[[3]][,1],
                                   y = colSums(varpart.env.spatial[[1]])[selected_groups],
                                   x.lab = "PCoA axis 1",
                                   y.lab = "Total variance explained\n by connectivity and local environment",
                                   y.lab.hjust = 0.5,
                                   x.log = F,
                                   y.log = F,
                                   fit = T,
                                   x.cor.pos = 0.8,
                                   y.cor.pos = 0.1,
                                   mar.vect = c(5,5,1,5))
  
  pdf(paste0(figure_folder,"/FigS_mean.Shannon.normalized.Moran.I.PCoA1_vs_diversity.tot.explained.var",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2*2,height=12/3*4/1.2/2*3,onefile=F)
  g1 = ggplotGrob(plot.shannon.vs.diversity)
  g2 = ggplotGrob(plot.tot.var.vs.shannon)
  g3 = ggplotGrob(plot.Moran.I.vs.diversity)
  g4 = ggplotGrob(plot.tot.var.vs.Moran.I)
  g5 = ggplotGrob(plot.PCoA1.vs.diversity)
  g6 = ggplotGrob(plot.tot.var.vs.PCoA1)
  grid.arrange(grobs = list(g1, g2, g3, g4, g5, g6), nrow = 3, ncol = 2)
  dev.off()
}

# Fig. S Shannon - Moran I - PCoA1 correlations:
{
  plot.shannon.vs.Moran.I = cor.plot(x = I_square.observed_w.mean[selected_groups,1]/optimalK_prevalence.min.crossValid[selected_groups],
                                     y = shannon_SUR[selected_groups],
                                     excluded.points = which(taxo_groups[selected_groups] == "Porifera"),
                                     x.lab = "Short-distance spatial autocorrelation",
                                     y.lab = "Mean Shannon diversity across stations",
                                     y.lab.hjust = 0.5,
                                     x.log = F,
                                     y.log = F,
                                     fit = T,
                                     x.cor.pos = 0.1,
                                     y.cor.pos = 0.1,
                                     mar.vect = c(5,5,1,5))
  
  plot.pcoa1.vs.shannon = cor.plot(x = shannon_SUR[selected_groups],
                                   y = NormalizedVI_pcoa[[3]][,1],
                                   excluded.points = which(taxo_groups[selected_groups] == "Porifera"),
                                   x.lab = "Mean Shannon diversity across stations",
                                   y.lab = "PCoA axis 1",
                                   y.lab.hjust = 0.5,
                                   x.log = F,
                                   y.log = F,
                                   fit = T,
                                   x.cor.pos = 0.8,
                                   y.cor.pos = 0.8,
                                   mar.vect = c(5,5,1,5))
  
  plot.pcoa1.vs.Moran.I = cor.plot(x = I_square.observed_w.mean[selected_groups,1]/optimalK_prevalence.min.crossValid[selected_groups],
                                   y = NormalizedVI_pcoa[[3]][,1],
                                   x.lab = "Short-distance spatial autocorrelation",
                                   y.lab = "PCoA axis 1",
                                   y.lab.hjust = 0.5,
                                   x.log = F,
                                   y.log = F,
                                   fit = T,
                                   x.cor.pos = 0.8,
                                   y.cor.pos = 0.1,
                                   mar.vect = c(5,5,1,5))
  
  pdf(paste0(figure_folder,"/FigS_Shannon.normalized.Moran.I.PCoA1.correlations",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2*3,height=12/3*4/1.2/2,onefile=F)
  g1 = ggplotGrob(plot.shannon.vs.Moran.I)
  g2 = ggplotGrob(plot.pcoa1.vs.shannon)
  g3 = ggplotGrob(plot.pcoa1.vs.Moran.I)
  grid.arrange(grobs = list(g1, g2, g3), nrow = 1, ncol = 3)
  dev.off()
  
  #######
  y1 = as.vector(NormalizedVI_pcoa[[3]][taxo_groups[selected_groups] != "Porifera",1])
  x1.1 = shannon_SUR[selected_groups & taxo_groups != "Porifera"]
  x1.2 = shannon_total_SUR[selected_groups & taxo_groups != "Porifera"]
  x2 = as.vector(I_square.observed_w.mean[selected_groups & taxo_groups != "Porifera",1])
  x3 = nb_absolute_dominants[selected_groups & taxo_groups != "Porifera",1]/optimalK_prevalence.min.crossValid[selected_groups & taxo_groups != "Porifera"]
  ###
  model1 = lm(y1 ~ x1.1 + x2)
  model2 = lm(y1 ~ x1.1)
  model3 = lm(y1 ~ x2)
  varpart1 = varpart(y1,x1.1,x2)
  #varpart1$part$indfract
  anova(rda(y1,x1.1,x2))
  anova(rda(y1,x2,x1.1))
  ###
  # model4 = lm(y1 ~ x3 + x2)
  # model5 = lm(y1 ~ x3 + x2 + x1.1)
  model5 = lm(y1 ~ x1.2 + x2 + x1.1)
  model6 = lm(y1 ~ x3 + x2 + x1.1 + x1.2)
  ###
  model1 = lm(y1 ~ x1.1 + x1.2)
  model2 = lm(y1 ~ x1.1)
  anova(model1,model2)
  model3 = lm(y1 ~ x1.2)
  #
  model1 = lm(y1 ~ x1.1 + x1.2 + x2)
  model2 = lm(y1 ~ x1.1 + x1.2)
  anova(model1,model2)
  model3 = lm(y1 ~ x1.1 + x2)
  anova(model1,model3)
  varpart1 = varpart(y1,x1.1,x1.2,x2)
  #varpart1$part$indfract
  varpart2 = varpart(y1,x1.1,x1.2)
  #varpart2$part$indfract
  ###
  model1 = lm(y1 ~ x3)
  model2 = lm(y1 ~ x1.1 + x1.2)
  varpart1 = varpart(y1,cbind(x1.1,x1.2),x2,x3)
}

# Fig. S autocorrelation slope - intercept - Moran I - charac scale correlations:
{
  plot.intercept.vs.Moran.I = cor.plot(x = I_square.observed_w.mean[selected_groups,1],
                                       y = intercept[selected_groups,1],
                                       # excluded.points = which(taxo_groups[selected_groups] == "Porifera"),
                                       x.lab = "Short-distance spatial autocorrelation",
                                       y.lab = "Spatial autocorrelation intercept",
                                       y.lab.hjust = 0.5,
                                       x.log = F,
                                       y.log = F,
                                       fit = T,
                                       x.cor.pos = 0.2,
                                       y.cor.pos = 0.8,
                                       mar.vect = c(5,5,1,5))
  
  plot.slope.vs.charac.scale = cor.plot(x = charac_scale[selected_groups,1],
                                        y = slope[selected_groups,1],
                                        # excluded.points = which(taxo_groups[selected_groups] == "Porifera"),
                                        x.lab = "Characteristic scale of spatial autocorrelation (km)",
                                        y.lab = "Spatial autocorrelation slope",
                                        y.lab.hjust = 0.5,
                                        x.log = F,
                                        y.log = F,
                                        fit = T,
                                        x.cor.pos = 0.8,
                                        y.cor.pos = 0.8,
                                        mar.vect = c(5,5,1,5))
  
  plot.slope.vs.Moran.I = cor.plot(x = I_square.observed_w.mean[selected_groups,1],
                                   y = slope[selected_groups,1],
                                   # excluded.points = which(taxo_groups[selected_groups] == "Porifera"),
                                   x.lab = "Short-distance spatial autocorrelation",
                                   y.lab = "Spatial autocorrelation slope",
                                   y.lab.hjust = 0.5,
                                   x.log = F,
                                   y.log = F,
                                   fit = T,
                                   x.cor.pos = 0.8,
                                   y.cor.pos = 0.8,
                                   mar.vect = c(5,5,1,5))
  
  plot.intercept.vs.charac.scale = cor.plot(x = charac_scale[selected_groups,1],
                                            y = intercept[selected_groups,1],
                                            # excluded.points = which(taxo_groups[selected_groups] == "Porifera"),
                                            x.lab = "Characteristic scale of spatial autocorrelation (km)",
                                            y.lab = "Spatial autocorrelation intercept",
                                            y.lab.hjust = 0.5,
                                            x.log = F,
                                            y.log = F,
                                            fit = T,
                                            x.cor.pos = 0.2,
                                            y.cor.pos = 0.8,
                                            mar.vect = c(5,5,1,5))
  
  plot.slope.vs.intercept = cor.plot(x = intercept[selected_groups,1],
                                     y = slope[selected_groups,1],
                                     # excluded.points = which(taxo_groups[selected_groups] == "Porifera"),
                                     x.lab = "Spatial autocorrelation intercept",
                                     y.lab = "Spatial autocorrelation slope",
                                     y.lab.hjust = 0.5,
                                     x.log = F,
                                     y.log = F,
                                     fit = T,
                                     x.cor.pos = 0.8,
                                     y.cor.pos = 0.8,
                                     mar.vect = c(5,5,1,5))
  
  plot.Moran.I.vs.charac.scale = cor.plot(x = charac_scale[selected_groups,1],
                                          y = I_square.observed_w.mean[selected_groups,1],
                                          # excluded.points = which(taxo_groups[selected_groups] == "Porifera"),
                                          x.lab = "Characteristic scale of spatial autocorrelation (km)",
                                          y.lab = "Short-distance spatial autocorrelation",
                                          y.lab.hjust = 0.5,
                                          x.log = F,
                                          y.log = F,
                                          fit = T,
                                          x.cor.pos = 0.2,
                                          y.cor.pos = 0.8,
                                          mar.vect = c(5,5,1,5))
  
  pdf(paste0(figure_folder,"/FigS_intercept.slope.Moran.I.charac.scale.correlations",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2*2,height=12/3*4/1.2/2*3,onefile=F)
  g1 = ggplotGrob(plot.intercept.vs.Moran.I)
  g2 = ggplotGrob(plot.slope.vs.charac.scale)
  g3 = ggplotGrob(plot.slope.vs.Moran.I)
  g4 = ggplotGrob(plot.intercept.vs.charac.scale)
  g5 = ggplotGrob(plot.slope.vs.intercept)
  g6 = ggplotGrob(plot.Moran.I.vs.charac.scale)
  grid.arrange(grobs = list(g1, g2, g3, g4, g5, g6), nrow = 3, ncol = 2)
  dev.off()
}

# Fig. S autocorrelation slope vs. PCoA2 - body size:
{
  plot.slope.vs.PCoA2 = cor.plot(x = NormalizedVI_pcoa[[3]][,2],
                                 y = slope[selected_groups,1],
                                 # excluded.points = which(taxo_groups[selected_groups] == "Porifera"),
                                 x.lab = "PCoA axis 2",
                                 y.lab = "Spatial autocorrelation slope",
                                 y.lab.hjust = 0.5,
                                 x.log = F,
                                 y.log = F,
                                 fit = T,
                                 x.cor.pos = 0.8,
                                 y.cor.pos = 0.8,
                                 mar.vect = c(5,5,1,5))
  
  plot.slope.vs.body.size = cor.plot(x = size_relativeAbund[selected_groups],
                                     y = slope[selected_groups,1],
                                     # excluded.points = which(taxo_groups[selected_groups] == "Porifera"),
                                     x.lab = expression("Mean body size ("*mu*"m)"),
                                     y.lab = "Spatial autocorrelation slope",
                                     y.lab.hjust = 0.5,
                                     x.log = T,
                                     y.log = F,
                                     fit = T,
                                     x.cor.pos = 0.8,
                                     y.cor.pos = 0.8,
                                     mar.vect = c(5,5,1,5))
  
  pdf(paste0(figure_folder,"/FigS_slope_vs_PCoA2.body.size",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2*2,height=12/3*4/1.2/2,onefile=F)
  g1 = ggplotGrob(plot.slope.vs.PCoA2)
  g2 = ggplotGrob(plot.slope.vs.body.size)
  grid.arrange(grobs = list(g1, g2), nrow = 1, ncol = 2)
  dev.off()
}

# Fig. S Shannon - Moran I - PCoA1 vs. K:
{
  plot.shannon.vs.K = cor.plot(x = optimalK_prevalence.min.crossValid[selected_groups],
                               y = shannon_SUR[selected_groups],
                               excluded.points = which(taxo_groups[selected_groups] == "Porifera"),
                               x.lab = "Number of assemblages",
                               y.lab = "Mean Shannon diversity across stations",
                               y.lab.hjust = 0.5,
                               x.log = F,
                               y.log = F,
                               fit = T,
                               x.cor.pos = 0.1,
                               y.cor.pos = 0.1,
                               mar.vect = c(5,5,1,5))
  
  plot.pcoa1.vs.K = cor.plot(x = optimalK_prevalence.min.crossValid[selected_groups],
                             y = NormalizedVI_pcoa[[3]][,1],
                             x.lab = "Number of assemblages",
                             y.lab = "PCoA axis 1",
                             y.lab.hjust = 0.5,
                             x.log = F,
                             y.log = F,
                             fit = T,
                             x.cor.pos = 0.8,
                             y.cor.pos = 0.8,
                             mar.vect = c(5,5,1,5))
  
  plot.Moran.I.vs.K = cor.plot(x = optimalK_prevalence.min.crossValid[selected_groups],
                               y = I_square.observed_w.mean[selected_groups,1],
                               x.lab = "Number of assemblages",
                               y.lab = "Short-distance spatial autocorrelation",
                               y.lab.hjust = 0.5,
                               x.log = F,
                               y.log = F,
                               fit = T,
                               x.cor.pos = 0.8,
                               y.cor.pos = 0.1,
                               mar.vect = c(5,5,1,5))
  
  pdf(paste0(figure_folder,"/FigS_Shannon.Moran.I.PCoA1_vs_K",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2*3,height=12/3*4/1.2/2,onefile=F)
  g1 = ggplotGrob(plot.shannon.vs.K)
  g2 = ggplotGrob(plot.pcoa1.vs.K)
  g3 = ggplotGrob(plot.Moran.I.vs.K)
  grid.arrange(grobs = list(g1, g2, g3), nrow = 1, ncol = 3)
  dev.off()
}

# Fig. S Shannon - Moran I - PCoA1 vs. nb. dominant assemblages:
{
  plot.shannon.vs.nb.dominants = cor.plot(x = nb_dominants[selected_groups,1],
                                          y = shannon_SUR[selected_groups],
                                          excluded.points = which(taxo_groups[selected_groups] == "Porifera"),
                                          x.lab = "Number of dominant assemblages",
                                          y.lab = "Mean Shannon diversity across stations",
                                          y.lab.hjust = 0.5,
                                          x.log = F,
                                          y.log = F,
                                          fit = T,
                                          x.cor.pos = 0.1,
                                          y.cor.pos = 0.8,
                                          mar.vect = c(5,5,1,5))
  
  plot.pcoa1.vs.nb.dominants = cor.plot(x = nb_dominants[selected_groups,1],
                                        y = NormalizedVI_pcoa[[3]][,1],
                                        x.lab = "Number of dominant assemblages",
                                        y.lab = "PCoA axis 1",
                                        y.lab.hjust = 0.5,
                                        x.log = F,
                                        y.log = F,
                                        fit = T,
                                        x.cor.pos = 0.8,
                                        y.cor.pos = 0.8,
                                        mar.vect = c(5,5,1,5))
  
  plot.Moran.I.vs.nb.dominants = cor.plot(x = nb_dominants[selected_groups,1],
                                          y = I_square.observed_w.mean[selected_groups,1],
                                          x.lab = "Number of dominant assemblages",
                                          y.lab = "Short-distance spatial autocorrelation",
                                          y.lab.hjust = 0.5,
                                          x.log = F,
                                          y.log = F,
                                          fit = T,
                                          x.cor.pos = 0.8,
                                          y.cor.pos = 0.8,
                                          mar.vect = c(5,5,1,5))
  
  pdf(paste0(figure_folder,"/FigS_Shannon.Moran.I.PCoA1_vs_nb.dominant",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2*3,height=12/3*4/1.2/2,onefile=F)
  g1 = ggplotGrob(plot.shannon.vs.nb.dominants)
  g2 = ggplotGrob(plot.pcoa1.vs.nb.dominants)
  g3 = ggplotGrob(plot.Moran.I.vs.nb.dominants)
  grid.arrange(grobs = list(g1, g2, g3), nrow = 1, ncol = 3)
  dev.off()
}

# Fig. S Shannon - Moran I - PCoA1 vs. nb. absolute dominant assemblages:
{
  plot.shannon.vs.nb.absolute.dominants = cor.plot(x = nb_absolute_dominants[selected_groups,1],
                                                   y = shannon_SUR[selected_groups],
                                                   excluded.points = which(taxo_groups[selected_groups] == "Porifera"),
                                                   x.lab = "Number of dominant (>50%) assemblages",
                                                   y.lab = "Mean Shannon diversity across stations",
                                                   y.lab.hjust = 0.5,
                                                   x.log = F,
                                                   y.log = F,
                                                   fit = T,
                                                   x.cor.pos = 0.1,
                                                   y.cor.pos = 0.8,
                                                   mar.vect = c(5,5,1,5))
  
  plot.pcoa1.vs.nb.absolute.dominants = cor.plot(x = nb_absolute_dominants[selected_groups,1],
                                                 y = NormalizedVI_pcoa[[3]][,1],
                                                 x.lab = "Number of dominant (>50%) assemblages",
                                                 y.lab = "PCoA axis 1",
                                                 y.lab.hjust = 0.5,
                                                 x.log = F,
                                                 y.log = F,
                                                 fit = T,
                                                 x.cor.pos = 0.8,
                                                 y.cor.pos = 0.8,
                                                 mar.vect = c(5,5,1,5))
  
  plot.Moran.I.vs.nb.absolute.dominants = cor.plot(x = nb_absolute_dominants[selected_groups,1],
                                                   y = I_square.observed_w.mean[selected_groups,1],
                                                   x.lab = "Number of dominant (>50%) assemblages",
                                                   y.lab = "Short-distance spatial autocorrelation",
                                                   y.lab.hjust = 0.5,
                                                   x.log = F,
                                                   y.log = F,
                                                   fit = T,
                                                   x.cor.pos = 0.8,
                                                   y.cor.pos = 0.8,
                                                   mar.vect = c(5,5,1,5))
  
  pdf(paste0(figure_folder,"/FigS_Shannon.Moran.I.PCoA1_vs_nb.absolute.dominant",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2*3,height=12/3*4/1.2/2,onefile=F)
  g1 = ggplotGrob(plot.shannon.vs.nb.absolute.dominants)
  g2 = ggplotGrob(plot.pcoa1.vs.nb.absolute.dominants)
  g3 = ggplotGrob(plot.Moran.I.vs.nb.absolute.dominants)
  grid.arrange(grobs = list(g1, g2, g3), nrow = 1, ncol = 3)
  dev.off()
}

# Fig. S OTU richness vs abundance
{
  plot.diversity.vs.abundance = cor.plot(x = colMeans(relativeAbund[,selected_groups],na.rm = T),
                                         y = as.vector(diversity)[selected_groups],
                                         x.lab = "Mean relative abundance across stations",
                                         y.lab = "Diversity (#OTU)",
                                         y.lab.hjust = 0.5,
                                         x.log = T,
                                         y.log = T,
                                         fit = T,
                                         x.cor.pos = 0.2,
                                         y.cor.pos = 0.8,
                                         mar.vect = c(5,5,1,5))
  
  pdf(paste0(figure_folder,"/FigS_diversity_vs_mean.rel.abundance",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.diversity.vs.abundance)
  dev.off()
  
  plot.diversity.vs.tot_reads = cor.plot(x = tot_reads[selected_groups],
                                         y = as.vector(diversity)[selected_groups],
                                         x.lab = "Total number of reads across stations and fractions",
                                         y.lab = "Diversity (#OTU)",
                                         y.lab.hjust = 0.5,
                                         x.log = T,
                                         y.log = T,
                                         fit = T,
                                         x.cor.pos = 0.2,
                                         y.cor.pos = 0.8,
                                         mar.vect = c(5,5,1,5))
  
  pdf(paste0(figure_folder,"/FigS_diversity_vs_tot.reads",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.diversity.vs.tot_reads)
  dev.off()
}

# Fig. S OTU richness vs diversity divided by relative abundance
{
  plot.norm.diversity.vs.diversity = cor.plot(x = as.vector(diversity)[selected_groups],
                                              y = as.vector(diversity)[selected_groups]/colMeans(relativeAbund[,selected_groups],na.rm = T),
                                              x.lab = "Diversity (#OTU)",
                                              y.lab = "Diversity normalized by relative abundance",
                                              y.lab.hjust = 0.5,
                                              x.log = T,
                                              y.log = T,
                                              fit = T,
                                              x.cor.pos = 0.8,
                                              y.cor.pos = 0.8,
                                              mar.vect = c(5,5,1,5))
  
  pdf(paste0(figure_folder,"/FigS_diversity_vs_diversity.normalized.by.rel.abundance",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.norm.diversity.vs.diversity)
  dev.off()
}

# Fig. S comparison mean Shannon across stations, Moran I and PCoA1 vs. diversity divided by relative abundance:
{
  div_threshold = 100
  
  x = (as.vector(diversity)/colMeans(relativeAbund,na.rm = T))[selected_groups & taxo_groups != "Porifera"]
  y = shannon_SUR[selected_groups & taxo_groups != "Porifera"]
  x.2000 = x[as.vector(diversity)[selected_groups & taxo_groups != "Porifera"] < 2000]
  y.2000 = y[as.vector(diversity)[selected_groups & taxo_groups != "Porifera"] < 2000]
  x.cor.pos = 0.1
  y.cor.pos = 0.8
  cor.test = cor.test(log10(x.2000),y.2000,na.rm=T)
  plot.shannon.vs.norm.diversity = cor.plot(x = x,
                                       y = y,
                                       x.lab = "Diversity normalized by relative abundance",
                                       y.lab = "Mean Shannon entropy across stations",
                                       y.lab.hjust = 0.5,
                                       x.log = T,
                                       y.log = F,
                                       fit = F,
                                       mar.vect = c(5,5,1,5)) +
    geom_smooth(data = data.frame(x = x.2000, y = y.2000),aes(x,y),method='lm',col="black",se=T) +
    # geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    annotate(geom="text", 
             x=(max(x)/min(x))^x.cor.pos*min(x),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=8)
  
  x = (as.vector(diversity)/colMeans(relativeAbund,na.rm = T))[selected_groups]
  y = I_square.observed_w.mean[selected_groups,1]/optimalK_prevalence.min.crossValid[selected_groups]
  x.2000 = x[as.vector(diversity)[selected_groups] < 2000]
  y.2000 = y[as.vector(diversity)[selected_groups] < 2000]
  x.cor.pos = 0.8
  y.cor.pos = 0.8
  cor.test = cor.test(log10(x.2000),y.2000,na.rm=T)
  plot.Moran.I.vs.norm.diversity = cor.plot(x = x,
                                       y = y,
                                       x.lab = "Diversity normalized by relative abundance",
                                       y.lab = "Short-distance spatial autocorrelation",
                                       y.lab.hjust = 0.5,
                                       x.log = T,
                                       y.log = F,
                                       fit = F,
                                       mar.vect = c(5,5,1,5)) +
    geom_smooth(data = data.frame(x = x.2000, y = y.2000),aes(x,y),method='lm',col="black",se=T) +
    # geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    annotate(geom="text", 
             x=(max(x)/min(x))^x.cor.pos*min(x),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=8)
  
  x = (as.vector(diversity)/colMeans(relativeAbund,na.rm = T))[selected_groups]
  y = NormalizedVI_pcoa[[3]][,1]
  x.2000 = x[as.vector(diversity)[selected_groups] < 2000]
  y.2000 = y[as.vector(diversity)[selected_groups] < 2000]
  x.cor.pos = 0.8
  y.cor.pos = 0.8
  cor.test = cor.test(log10(x.2000),y.2000,na.rm=T)
  plot.PCoA1.vs.norm.diversity = cor.plot(x = x,
                                     y = y,
                                     x.lab = "Diversity normalized by relative abundance",
                                     y.lab = "PCoA axis 1",
                                     y.lab.hjust = 0.5,
                                     x.log = T,
                                     y.log = F,
                                     fit = F,
                                     mar.vect = c(5,5,1,5)) +
    geom_smooth(data = data.frame(x = x.2000, y = y.2000),aes(x,y),method='lm',col="black",se=T) +
    # geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    annotate(geom="text", 
             x=(max(x)/min(x))^x.cor.pos*min(x),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos*(max(y)-min(y))+min(y),
             label=bquote(atop(rho==.(format(cor.test$estimate,digits=2,nsmall=2)),p==.(format(cor.test$p.value,digits=1)))),
             size=8)
  
  pdf(paste0(figure_folder,"/FigS_mean.Shannon.normalized.Moran.I.PCoA1_vs_diversity.normalized.by.rel.abundance",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2*3,height=12/3*4/1.2/2,onefile=F)
  g1 = ggplotGrob(plot.shannon.vs.norm.diversity)
  g2 = ggplotGrob(plot.Moran.I.vs.norm.diversity)
  g3 = ggplotGrob(plot.PCoA1.vs.norm.diversity)
  grid.arrange(grobs = list(g1, g2, g3), nrow = 1, ncol = 3)
  dev.off()
}

#############################################

# Fig. S Moran I random groups vs. diversity
{
  # x.cor.pos = 0.15
  # y.cor.pos = 0.95
  plot.Moran.I.vs.log.diversity = list()
  plot.random.groups.Moran.I.vs.log.diversity = list()
  plot.Moran.I.vs.random.groups.Moran.I = list()
  for (i_case in 1:2)
  {
    y = I_square.observed_w.mean[selected_groups,i_case][!taxo_groups[selected_groups] %in% c("Dinophyceae","Diplonemida")]
    x = as.vector(diversity)[selected_groups][!taxo_groups[selected_groups] %in% c("Dinophyceae","Diplonemida")]
    cor.test = cor.test(log10(x[x<2000]),y[x<2000],na.rm=T,method="spearman")
    x.cor.pos = 0.8 
    y.cor.pos = if (i_case == 1) 0.4 else 0.5
    plot.Moran.I.vs.log.diversity[[i_case]] = cor.plot(y = y,
                                                       x = x,
                                                       # excluded.points = which(taxo_groups[selected_groups] %in% "RAD-A"),
                                                       y.lab = "Short-distance spatial autocorrelation\n in major plankton groups",
                                                       x.lab = "Diversity (#OTUs), log scale",
                                                       x.log = T,
                                                       y.log = F,
                                                       fit = F,
                                                       mar.vect=c(5,5,1,5)) +
      geom_smooth(data = data.frame(x = x[x<2000], y = y[x<2000]), aes(x,y),method='lm',col="black",se=T) +
      # geom_smooth(aes(x[x<2000], y[x<2000]),method='lm',col="black",se=T, inherit.aes = T) +
      geom_vline(xintercept = 2000, linetype="dashed") +
      annotate(geom="text", 
               x=(max(x)/min(x))^x.cor.pos*min(x),
               y=y.cor.pos*(max(y)-min(y))+min(y),
               label=bquote(atop(rho[S]==.(format(cor.test$estimate,digits=2,nsmall=2)),P==.(format(cor.test$p.value,digits=1)))),
               size=7.5) +
      ylim(0,0.8)
    
    y = I_square.observed_w.mean.random.groups[selected_groups,i_case][!taxo_groups[selected_groups] %in% c("Dinophyceae","Diplonemida")]
    x = as.vector(diversity)[selected_groups][!taxo_groups[selected_groups] %in% c("Dinophyceae","Diplonemida")]
    cor.test = cor.test(log10(x[x<2000]),y[x<2000],na.rm=T,method="spearman")
    x.cor.pos = 0.8
    y.cor.pos = if (i_case == 1) 0.35 else 0.6
    plot.random.groups.Moran.I.vs.log.diversity[[i_case]] = cor.plot(x = x,
                                                                     y = y,
                                                                     x.lab = "Diversity (#OTUs), log scale",
                                                                     y.lab = "Short-distance spatial autocorrelation\n in random mock groups",
                                                                     y.lab.hjust = 0.5,
                                                                     x.log = T,
                                                                     y.log = F,
                                                                     fit = F,
                                                                     fit.display="pearson.spearman",
                                                                     mar.vect = c(5,5,1,5)) +
      geom_smooth(data = data.frame(x = x[x<2000], y = y[x<2000]), aes(x,y),method='lm',col="black",se=T) +
      # geom_smooth(aes(x[x<2000], y[x<2000]),method='lm',col="black",se=T, inherit.aes = T) +
      geom_vline(xintercept = 2000, linetype="dashed") +
      annotate(geom="text", 
               x=(max(x)/min(x))^x.cor.pos*min(x),
               y=y.cor.pos*(max(y)-min(y))+min(y),
               label=bquote(atop(rho[S]==.(format(cor.test$estimate,digits=2,nsmall=2)),P==.(format(cor.test$p.value,digits=1)))),
               size=7.5) +
      ylim(0,0.8)
    
    x = I_square.observed_w.mean.random.groups[selected_groups,i_case][!taxo_groups[selected_groups] %in% c("Dinophyceae","Diplonemida")]
    y = I_square.observed_w.mean[selected_groups,i_case][!taxo_groups[selected_groups] %in% c("Dinophyceae","Diplonemida")]
    plot.Moran.I.vs.random.groups.Moran.I[[i_case]] = cor.plot(x = x,
                                                                     y = y,
                                                                     x.lab = "Short-distance spatial autocorrelation\n in random mock groups",
                                                                     y.lab = "Short-distance spatial autocorrelation\n in major plankton groups",
                                                                     y.lab.hjust = 0.5,
                                                                     x.log = F,
                                                                     y.log = F,
                                                                     fit = T,
                                                                     fit.display="pearson.p",
                                                                     x.cor.pos=0.35,
                                                                     y.cor.pos=0.85,
                                                                     mar.vect = c(5,5,1,5))
  }
  
  pdf(paste0(figure_folder,"/FigS_random.groups.Moran.I_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2*2,height=12/3*4/1.2/2*2,onefile=F)
  g1 = ggplotGrob(plot.Moran.I.vs.log.diversity[[1]])
  g2 = ggplotGrob(plot.random.groups.Moran.I.vs.log.diversity[[1]])
  g3 = ggplotGrob(plot.Moran.I.vs.log.diversity[[2]])
  g4 = ggplotGrob(plot.random.groups.Moran.I.vs.log.diversity[[2]])
  # g5 = ggplotGrob(plot.Moran.I.vs.random.groups.Moran.I[[1]])
  # g6 = ggplotGrob(plot.Moran.I.vs.random.groups.Moran.I[[2]])
  grid.arrange(grobs = list(g1,g2,g3,g4),
               layout_matrix = rbind(c(1, 2),c(3, 4)))
  # print(plot.random.groups.Moran.I.vs.diversity)
  dev.off()
}

# Fig. S Moran I rarefied (equal read number: 10^4) and non-rarefied vs. diversity
{
  x = as.vector(diversity)[selected_groups & tot_reads > 10^4 & taxo_groups != "Porifera"]
  y = I_square.observed_w.mean[selected_groups & tot_reads > 10^4 & taxo_groups != "Porifera",1]
  y.104 = I_square.observed_w.mean.104[selected_groups & tot_reads > 10^4 & taxo_groups != "Porifera",1]
  x.2000 = x[as.vector(diversity)[selected_groups & tot_reads > 10^4 & taxo_groups != "Porifera"] < 2000]
  y.2000 = y[as.vector(diversity)[selected_groups & tot_reads > 10^4 & taxo_groups != "Porifera"] < 2000]
  y.104.2000 = y.104[as.vector(diversity)[selected_groups & tot_reads > 10^4 & taxo_groups != "Porifera"] < 2000]
  x.cor.pos = 0.1
  y.cor.pos = 0.9
  y.cor.pos.104 = 0.8
  cor.test = cor.test(log10(x.2000),y.2000,na.rm=T)
  if (cor.test$p.value > 0.05)
  {
    stars = ""
  } else if (cor.test$p.value < 0.05 & cor.test$p.value > 0.01)
  {
    stars = "*"
  } else if (cor.test$p.value < 0.01 & cor.test$p.value > 0.001)
  {
    stars = "**"
  } else if (cor.test$p.value < 0.001)
    stars = "***"
  cor.test.104 = cor.test(log10(x.2000),y.104.2000,na.rm=T)
  if (cor.test.104$p.value > 0.05)
  {
    stars.104 = ""
  } else if (cor.test.104$p.value < 0.05 & cor.test.104$p.value > 0.01)
  {
    stars.104 = "*"
  } else if (cor.test.104$p.value < 0.01 & cor.test.104$p.value > 0.001)
  {
    stars.104 = "**"
  } else if (cor.test.104$p.value < 0.001)
    stars.104 = "***"
  plot.raref.Moran.I.vs.diversity = cor.plot(x = x,
                                            y = y,
                                            x.lab = "Diversity (#OTUs)",
                                            y.lab = "Short-distance spatial autocorrelation",
                                            y.lab.hjust = 0.5,
                                            x.log = T,
                                            y.log = F,
                                            fit = F,
                                            mar.vect = c(5,5,1,5)) +
    geom_smooth(data = data.frame(x = x.2000, y = y.2000),aes(x,y),method='lm',col="black",se=T) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    annotate(geom="text", 
             x=(max(x)/min(x))^x.cor.pos*min(x),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos*(max(y)-min(y))+min(y),
             label=bquote(rho==.(format(cor.test$estimate,digits=2,nsmall=2))^.(stars)),
             size=8) +
    geom_point(data = data.frame(x = x, y = y.104), aes(x,y.104), colour = "blue", size = 1, alpha=0.5) +
    geom_smooth(data = data.frame(x = x.2000, y = y.104.2000),aes(x,y),method='lm',col="blue",se=T, alpha=0.5) +
    annotate(geom="text", 
             x=(max(x)/min(x))^x.cor.pos*min(x),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos.104*(max(y)-min(y))+min(y),
             label=bquote(rho==.(format(cor.test.104$estimate,digits=2,nsmall=2))^.(stars.104)),
             colour="blue",
             size=8)
  pdf(paste0(figure_folder,"/FigS_rarefied.Moran.I_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.Moran.I.vs.diversity)
  dev.off()
  
  #####################
  x.cor.pos = 0.15
  y.cor.pos.104 = 0.95
  y.104 = I_square.observed_w.mean.104[selected_groups & tot_reads > 10^4,1]
  x.104 = diversity[selected_groups & tot_reads > 10^4]
  y.104.2000 = y.104[x.104 < 2000 & 
                       !names(x.104) %in% c("Porifera")]
  x.104.2000 = x.104[x.104 < 2000 & 
                       !names(x.104) %in% c("Porifera")]
  cor.test.104 = cor.test(log10(x.104.2000),y.104.2000,na.rm=T)
  if (cor.test.104$p.value > 0.05)
  {
    stars.104 = ""
  } else if (cor.test.104$p.value < 0.05 & cor.test.104$p.value > 0.01)
  {
    stars.104 = "*"
  } else if (cor.test.104$p.value < 0.01 & cor.test.104$p.value > 0.001)
  {
    stars.104 = "**"
  } else if (cor.test.104$p.value < 0.001)
    stars.104 = "***"
  excluded_groups = "Porifera"
  plot.raref.104.Moran.I.only.vs.diversity = cor.plot(x = x.104,
                                                      y = y.104,
                                                      excluded.points = which(taxo_groups[selected_groups & tot_reads > 10^4] %in% excluded_groups),
                                                      x.lab = "Diversity (#OTUs)",
                                                      y.lab = "Short-distance spatial autocorrelation",
                                                      y.lab.hjust = 0.5,
                                                      x.log = T,
                                                      y.log = F,
                                                      fit = F,
                                                      mar.vect = c(5,5,1,5)) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    geom_smooth(data = data.frame(x = as.vector(x.104.2000), y = y.104.2000),aes(x,y),method='lm',se=T,colour="black") +
    annotate(geom="text", 
             x=(max(x.104)/min(x.104))^x.cor.pos*min(x.104),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos.104*(max(y.104)-min(y.104))+min(y.104),
             label=bquote(rho==.(format(cor.test.104$estimate,digits=2,nsmall=2))^.(stars.104)),
             size=8)
  pdf(paste0(figure_folder,"/FigS_rarefied.104.Moran.I.only_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.104.Moran.I.only.vs.diversity)
  dev.off()
  
  ####### Diversity after and before rarefaction
  x = as.vector(diversity)[selected_groups & tot_reads>10^4 & taxo_groups != "Porifera"]
  x.104 = diversity.104[selected_groups & tot_reads>10^4 & taxo_groups != "Porifera"]
  plot.raref.diversity.vs.diversity = cor.plot(x = x,
                                               y = x.104,
                                               x.lab = "Diversity (#OTUs)",
                                               y.lab = "Diversity after 10^4-read rarefaction (#OTUs)",
                                               y.lab.hjust = 0.5,
                                               x.log = T,
                                               y.log = T,
                                               fit = T,
                                               fit.display = "pearson.spearman",
                                               mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.104.diversity_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.diversity.vs.diversity)
  dev.off()
  
  ####### Moran I vs rarefied diversity:
  x.104 = diversity.104[selected_groups & tot_reads>10^4 & taxo_groups != "Porifera"]
  y.104 = I_square.observed_w.mean.104[selected_groups & tot_reads>10^4 & taxo_groups != "Porifera",1]
  plot.raref.Moran.I.vs.raref.diversity = cor.plot(x = x.104,
                                                   y = y.104,
                                                   x.lab = "Diversity after 10^4-read rarefaction (#OTUs)",
                                                   y.lab = "Short-distance spatial autocorrelation",
                                                   y.lab.hjust = 0.5,
                                                   x.log = T,
                                                   y.log = F,
                                                   fit = T,
                                                   fit.display = "pearson.spearman",
                                                   mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.104.Moran.I_vs_raref.diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.Moran.I.vs.raref.diversity)
  dev.off()
}

# Fig. S Moran I rarefied (equal read number: 10^5) vs. diversity
{
  x.cor.pos = 0.15
  y.cor.pos.105 = 0.95
  y.105 = I_square.observed_w.mean.105[selected_groups & tot_reads > 10^5,1]
  x.105 = diversity[selected_groups & tot_reads > 10^5]
  y.105.2000 = y.105[x.105 < 2000 & 
                       !names(x.105) %in% c("Porifera")]
  x.105.2000 = x.105[x.105 < 2000 & 
                       !names(x.105) %in% c("Porifera")]
  cor.test.105 = cor.test(log10(x.105.2000),y.105.2000,na.rm=T)
  if (cor.test.105$p.value > 0.05)
  {
    stars.105 = ""
  } else if (cor.test.105$p.value < 0.05 & cor.test.105$p.value > 0.01)
  {
    stars.105 = "*"
  } else if (cor.test.105$p.value < 0.01 & cor.test.105$p.value > 0.001)
  {
    stars.105 = "**"
  } else if (cor.test.105$p.value < 0.001)
    stars.105 = "***"
  excluded_groups = "Porifera"
  plot.raref.105.Moran.I.only.vs.diversity = cor.plot(x = x.105,
                                                      y = y.105,
                                                      excluded.points = which(taxo_groups[selected_groups & tot_reads > 10^5] %in% excluded_groups),
                                                      x.lab = "Diversity (#OTUs)",
                                                      y.lab = "Short-distance spatial autocorrelation",
                                                      y.lab.hjust = 0.5,
                                                      x.log = T,
                                                      y.log = F,
                                                      fit = F,
                                                      mar.vect = c(5,5,1,5)) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    geom_smooth(data = data.frame(x = as.vector(x.105.2000), y = y.105.2000), aes(x,y), method='lm',se=T,colour="black") +
    annotate(geom="text", 
             x=(max(x.105)/min(x.105))^x.cor.pos*min(x.105),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos.105*(max(y.105)-min(y.105))+min(y.105),
             label=bquote(rho==.(format(cor.test.105$estimate,digits=2,nsmall=2))^.(stars.105)),
             size=8)
  pdf(paste0(figure_folder,"/FigS_rarefied.105.Moran.I.only_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.105.Moran.I.only.vs.diversity)
  dev.off()
  
  ####### Moran I vs rarefied 10^5 diversity
  x = diversity[selected_groups & tot_reads>10^5 & taxo_groups != "Porifera"]
  y.105 = diversity.105[selected_groups & tot_reads>10^5 & taxo_groups != "Porifera"]
  plot.raref.105.diversity.vs.diversity = cor.plot(x = x,
                                                   y = y.105,
                                                   x.lab = "Diversity",
                                                   y.lab = "Diversity after 10^5-read rarefaction (#OTUs)",
                                                   y.lab.hjust = 0.5,
                                                   x.log = T,
                                                   y.log = T,
                                                   fit = T,
                                                   fit.display = "pearson.spearman",
                                                   mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.105.diversity_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.105.diversity.vs.diversity)
  dev.off()
  
  #######  ####### Moran I vs rarefied 105 nb of reads
  x.105 = diversity.105[selected_groups & tot_reads>10^5 & taxo_groups != "Porifera"]
  y.105 = I_square.observed_w.mean.105[selected_groups & tot_reads>10^5 & taxo_groups != "Porifera",1]
  plot.raref.Moran.I.105.vs.raref.diversity = cor.plot(x = x.105,
                                                       y = y.105,
                                                       x.lab = "Diversity after 10^5-read rarefaction (#OTUs)",
                                                       y.lab = "Short-distance spatial autocorrelation",
                                                       y.lab.hjust = 0.5,
                                                       x.log = T,
                                                       y.log = F,
                                                       fit = T,
                                                       fit.display = "pearson.spearman",
                                                       mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.105.Moran.I_vs_raref.105.diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.Moran.I.105.vs.raref.diversity)
  dev.off()
}

# Fig. S Moran I rarefied (equal diversity: 200) and non-rarefied vs. diversity
{
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera"]
  y = I_square.observed_w.mean[selected_groups & taxo_groups != "Porifera",1]
  y.200 = I_square.observed_w.mean.200[selected_groups & diversity > 200 & taxo_groups != "Porifera",1]
  x.200 = diversity[selected_groups & diversity > 200 & taxo_groups != "Porifera"]
  x.2000 = x[x < 2000]
  y.2000 = y[x < 2000]
  y.200.2000 = y.200[x.200 < 2000 & 
                    !names(x.200) %in% c("MAST-4,_6,_7,_8,_9,_10,_11","Phaeodaria","Nassellaria_and_Eucyrtidium_*")]
  x.200.2000 = x.200[x.200 < 2000 & 
                    !names(x.200) %in% c("MAST-4,_6,_7,_8,_9,_10,_11","Phaeodaria","Nassellaria_and_Eucyrtidium_*")]
  x.cor.pos = 0.1
  y.cor.pos = 0.9
  y.cor.pos.200 = 0.8
  cor.test = cor.test(log10(x.2000),y.2000,na.rm=T)
  if (cor.test$p.value > 0.05)
  {
    stars = ""
  } else if (cor.test$p.value < 0.05 & cor.test$p.value > 0.01)
  {
    stars = "*"
  } else if (cor.test$p.value < 0.01 & cor.test$p.value > 0.001)
  {
    stars = "**"
  } else if (cor.test$p.value < 0.001)
    stars = "***"
  cor.test.200 = cor.test(log10(x.200.2000),y.200.2000,na.rm=T)
  if (cor.test.200$p.value > 0.05)
  {
    stars.200 = ""
  } else if (cor.test.200$p.value < 0.05 & cor.test.200$p.value > 0.01)
  {
    stars.200 = "*"
  } else if (cor.test.200$p.value < 0.01 & cor.test.200$p.value > 0.001)
  {
    stars.200 = "**"
  } else if (cor.test.200$p.value < 0.001)
    stars.200 = "***"
  plot.raref.200.Moran.I.vs.diversity = cor.plot(x = x,
                                             y = y,
                                             x.lab = "Diversity (#OTUs)",
                                             y.lab = "Short-distance spatial autocorrelation",
                                             y.lab.hjust = 0.5,
                                             x.log = T,
                                             y.log = F,
                                             fit = F,
                                             mar.vect = c(5,5,1,5)) +
    geom_smooth(data = data.frame(x = x.2000, y = y.2000),aes(x,y),method='lm',col="black",se=T) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    annotate(geom="text", 
             x=(max(x)/min(x))^x.cor.pos*min(x),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos*(max(y)-min(y))+min(y),
             label=bquote(rho==.(format(cor.test$estimate,digits=2,nsmall=2))^.(stars)),
             size=8) +
    geom_point(data = data.frame(x = as.vector(x.200), y = y.200), aes(x,y), colour = "blue", size = 1, alpha=0.5) +
    geom_smooth(data = data.frame(x = as.vector(x.200.2000), y = y.200.2000),aes(x,y),method='lm',col="blue",se=T, alpha=0.5) +
    annotate(geom="text", 
             x=(max(x)/min(x))^x.cor.pos*min(x),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos.200*(max(y)-min(y))+min(y),
             label=bquote(rho==.(format(cor.test.200$estimate,digits=2,nsmall=2))^.(stars.200)),
             colour="blue",
             size=8)
  pdf(paste0(figure_folder,"/FigS_rarefied.200.Moran.I_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.200.Moran.I.vs.diversity)
  dev.off()
  
  ##########
  x.cor.pos = 0.15
  y.cor.pos.200 = 0.95
  y.200 = I_square.observed_w.mean.200[selected_groups & diversity > 200,1]
  x.200 = diversity[selected_groups & diversity > 200]
  y.200.2000 = y.200[x.200 < 2000 & 
                       !names(x.200) %in% c("Porifera","MAST-4,_6,_7,_8,_9,_10,_11","Phaeodaria","Nassellaria_and_Eucyrtidium_*")]
  x.200.2000 = x.200[x.200 < 2000 & 
                       !names(x.200) %in% c("Porifera","MAST-4,_6,_7,_8,_9,_10,_11","Phaeodaria","Nassellaria_and_Eucyrtidium_*")]
  excluded_groups = c("Porifera","MAST-4,_6,_7,_8,_9,_10,_11","Phaeodaria","Nassellaria_and_Eucyrtidium_*")
  plot.raref.200.Moran.I.only.vs.diversity = cor.plot(x = x.200,
                                                 y = y.200,
                                                 excluded.points = which(taxo_groups[selected_groups & diversity > 200] %in% excluded_groups),
                                                 x.lab = "Diversity (#OTUs)",
                                                 y.lab = "Short-distance spatial autocorrelation",
                                                 y.lab.hjust = 0.5,
                                                 x.log = T,
                                                 y.log = F,
                                                 fit = F,
                                                 mar.vect = c(5,5,1,5)) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    geom_smooth(data = data.frame(x = as.vector(x.200.2000), y = y.200.2000),aes(x,y),method='lm',se=T,colour="black") +
    annotate(geom="text", 
             x=(max(x.200)/min(x.200))^x.cor.pos*min(x.200),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos.200*(max(y.200)-min(y.200))+min(y.200),
             # label=bquote(rho==.(format(cor.test.200$estimate,digits=2,nsmall=2))^.(stars.200)),
             label=bquote(atop(rho==.(format(cor.test.200$estimate,digits=2,nsmall=2)),p==.(format(cor.test.200$p.value,digits=1)))),
             size=8)
  pdf(paste0(figure_folder,"/FigS_rarefied.200.Moran.I.only_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.200.Moran.I.only.vs.diversity)
  dev.off()
  
  ####### Moran I vs rarefied 200 diversity (residual variation around 200)
  x.200 = diversity.200[selected_groups & taxo_groups != "Porifera"]
  y.200 = I_square.observed_w.mean.200[selected_groups & taxo_groups != "Porifera",1]
  plot.raref.200.Moran.I.vs.raref.200.diversity = cor.plot(x = x.200,
                                                   y = y.200,
                                                   x.lab = "Diversity after 200-OTU rarefaction (#OTUs)",
                                                   y.lab = "Short-distance spatial autocorrelation",
                                                   y.lab.hjust = 0.5,
                                                   x.log = T,
                                                   y.log = F,
                                                   fit = T,
                                                   # fit.display = "pearson.spearman",
                                                   mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.200.Moran.I_vs_raref.200.diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.200.Moran.I.vs.raref.200.diversity)
  dev.off()
  
  #######  ####### Moran I vs rarefied 200 nb of reads
  x.200 = tot_reads.200[selected_groups & diversity>200 & taxo_groups != "Porifera"]
  y.200 = I_square.observed_w.mean.200[selected_groups  & diversity>200 & taxo_groups != "Porifera",1]
  plot.raref.Moran.I.200.vs.raref.nb.reads = cor.plot(x = x.200,
                                                      y = y.200,
                                                      x.lab = "Number of reads after 200-OTU rarefaction",
                                                      y.lab = "Short-distance spatial autocorrelation",
                                                      y.lab.hjust = 0.5,
                                                      x.log = T,
                                                      y.log = F,
                                                      fit = T,
                                                      # fit.display = "pearson.spearman",
                                                      mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.200.Moran.I_vs_raref.200.nb.reads",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.Moran.I.200.vs.raref.nb.reads)
  dev.off()
}

# Fig. S Moran I randomly rarefied (equal diversity: 200) vs. diversity
{
  data.lme = data.frame(Groups = taxo_groups[selected_groups & diversity > 200 & diversity < 2000 & taxo_groups != "Porifera"],
                        Diversity = as.vector(diversity)[selected_groups & diversity > 200 & diversity < 2000 & taxo_groups != "Porifera"],
                        I_square.observed_w.mean.200.random.5reals[[1]][selected_groups & diversity > 200 & diversity < 2000 & taxo_groups != "Porifera",])
  melted.data.lme = reshape2::melt(data.lme,id.vars=1:2,
                                   measure.vars=3:(5+2),
                                   variable.name="Real",
                                   value.name="Moran.I")
  melted.data.lme[,"Diversity"] = log10(melted.data.lme[,"Diversity"])
  
  test = lme(Diversity ~ Moran.I, data = melted.data.lme, random = ~1|Groups)
  
  y.200 = I_square.observed_w.mean.200.random.5reals[[1]][selected_groups & diversity > 200 & taxo_groups != "Porifera",]
  melted.y.200 = reshape2::melt(y.200)$value
  x.200 = as.vector(diversity)[selected_groups & diversity > 200 & taxo_groups != "Porifera"]
  melted.x.200 = rep(x.200,5)
  y.200.2000 = y.200[x.200 < 2000,]
  melted.y.200.2000 = reshape2::melt(y.200.2000)$value
  x.200.2000 = x.200[x.200 < 2000]
  melted.x.200.2000 = rep(x.200.2000,5)
  #
  x.cor.pos = 0.1
  y.cor.pos = 0.95
  cor.test.200 = cor.test(log10(melted.x.200.2000),melted.y.200.2000,na.rm=T)
  if (cor.test.200$p.value > 0.05)
  {
    stars.200 = ""
  } else if (cor.test.200$p.value < 0.05 & cor.test.200$p.value > 0.01)
  {
    stars.200 = "*"
  } else if (cor.test.200$p.value < 0.01 & cor.test.200$p.value > 0.001)
  {
    stars.200 = "**"
  } else if (cor.test.200$p.value < 0.001)
    stars.200 = "***"
  #
  plot.raref.200.random.5real.Moran.I.vs.diversity = cor.plot(x = melted.x.200,
                                                              y = melted.y.200,
                                                              x.lab = "Initial diversity (#OTUs)",
                                                              y.lab = "Short-distance spatial autocorrelation",
                                                              y.lab.hjust = 0.5,
                                                              x.log = T,
                                                              y.log = F,
                                                              fit = F,
                                                              mar.vect = c(5,5,1,5)) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(melted.y.200), yend = max(melted.y.200)), linetype="dashed", inherit.aes = T) +
    geom_smooth(data = data.frame(x = melted.x.200.2000, y = melted.y.200.2000),aes(x,y),method='lm',se=T) +
    annotate(geom="text",
             x=(max(melted.x.200)/min(melted.x.200))^x.cor.pos*min(melted.x.200),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos*(max(melted.y.200)-min(melted.y.200))+min(melted.y.200),
             label=bquote(rho==.(format(cor.test.200$estimate,digits=2,nsmall=2))^.(stars.200)),
             size=8)
  pdf(paste0(figure_folder,"/FigS_rarefied.200.random.reals.Moran.I_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.200.random.5real.Moran.I.vs.diversity)
  dev.off()
}

# Fig. S Moran I rarefied (equal diversity: 250) vs. diversity
{
  ##########
  x.cor.pos = 0.15
  y.cor.pos.250 = 0.95
  y.250 = I_square.observed_w.mean.250[selected_groups & diversity > 250,1]
  x.250 = diversity[selected_groups & diversity > 250]
  y.250.2000 = y.250[x.250 < 2000 & 
                       !names(x.250) %in% c("Porifera")]#,"Phaeodaria","Basidiomycota")]
  x.250.2000 = x.250[x.250 < 2000 & 
                       !names(x.250) %in% c("Porifera")]#,"Phaeodaria","Basidiomycota")]
  cor.test.250 = cor.test(log10(x.250.2000),y.250.2000,na.rm=T)
  if (cor.test.250$p.value > 0.05)
  {
    stars.250 = ""
  } else if (cor.test.250$p.value < 0.05 & cor.test.250$p.value > 0.01)
  {
    stars.250 = "*"
  } else if (cor.test.250$p.value < 0.01 & cor.test.250$p.value > 0.001)
  {
    stars.250 = "**"
  } else if (cor.test.250$p.value < 0.001)
    stars.250 = "***"
  excluded_groups = c("Porifera")#,"Phaeodaria","Basidiomycota")
  plot.raref.250.Moran.I.only.vs.diversity = cor.plot(x = x.250,
                                                      y = y.250,
                                                      excluded.points = which(taxo_groups[selected_groups & diversity > 250] %in% excluded_groups),
                                                      x.lab = "Diversity (#OTUs)",
                                                      y.lab = "Short-distance spatial autocorrelation",
                                                      y.lab.hjust = 0.5,
                                                      x.log = T,
                                                      y.log = F,
                                                      fit = F,
                                                      mar.vect = c(5,5,1,5)) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    geom_smooth(data = data.frame(x = as.vector(x.250.2000), y = y.250.2000),aes(x,y),method='lm',se=T,colour="black") +
    annotate(geom="text", 
             x=(max(x.250)/min(x.250))^x.cor.pos*min(x.250),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos.250*(max(y.250)-min(y.250))+min(y.250),
             # label=bquote(rho==.(format(cor.test.250$estimate,digits=2,nsmall=2))^.(stars.250)),
             label=bquote(atop(rho==.(format(cor.test.250$estimate,digits=2,nsmall=2)),p==.(format(cor.test.250$p.value,digits=1)))),
             size=8)
  pdf(paste0(figure_folder,"/FigS_rarefied.250.Moran.I.only_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.250.Moran.I.only.vs.diversity)
  dev.off()
  
  ####### Moran I vs rarefied 200 diversity (residual variation around 200)
  x.250 = diversity.250[selected_groups & taxo_groups != "Porifera"]
  y.250 = I_square.observed_w.mean.250[selected_groups & taxo_groups != "Porifera",1]
  plot.raref.250.Moran.I.vs.raref.250.diversity = cor.plot(x = x.250,
                                                           y = y.250,
                                                           x.lab = "Diversity after 250-OTU rarefaction (#OTUs)",
                                                           y.lab = "Short-distance spatial autocorrelation",
                                                           y.lab.hjust = 0.5,
                                                           x.log = T,
                                                           y.log = F,
                                                           fit = T,
                                                           # fit.display = "pearson.spearman",
                                                           mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.250.Moran.I_vs_raref.250.diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.250.Moran.I.vs.raref.250.diversity)
  dev.off()
  
  #######  ####### Moran I vs rarefied 250 nb of reads
  x.250 = tot_reads.250[selected_groups & diversity>250 & taxo_groups != "Porifera"]
  y.250 = I_square.observed_w.mean.250[selected_groups  & diversity>250 & taxo_groups != "Porifera",1]
  plot.raref.Moran.I.250.vs.raref.nb.reads = cor.plot(x = x.250,
                                                      y = y.250,
                                                      x.lab = "Number of reads after 250-OTU rarefaction",
                                                      y.lab = "Short-distance spatial autocorrelation",
                                                      y.lab.hjust = 0.5,
                                                      x.log = T,
                                                      y.log = F,
                                                      fit = T,
                                                      # fit.display = "pearson.spearman",
                                                      mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.250.Moran.I_vs_raref.250.nb.reads",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.Moran.I.250.vs.raref.nb.reads)
  dev.off()
}

# Fig. S Moran I rarefied (equal diversity: 300) vs. diversity
{
  ##########
  x.cor.pos = 0.15
  y.cor.pos.300 = 0.95
  y.300 = I_square.observed_w.mean.300[selected_groups & diversity > 300,1]
  x.300 = diversity[selected_groups & diversity > 300]
  y.300.2000 = y.300[x.300 < 2000 & 
                       !names(x.300) %in% c("Porifera")]#,"Phaeodaria","Basidiomycota")]
  x.300.2000 = x.300[x.300 < 2000 & 
                       !names(x.300) %in% c("Porifera")]#,"Phaeodaria","Basidiomycota")]
  cor.test.300 = cor.test(log10(x.300.2000),y.300.2000,na.rm=T)
  if (cor.test.300$p.value > 0.05)
  {
    stars.300 = ""
  } else if (cor.test.300$p.value < 0.05 & cor.test.300$p.value > 0.01)
  {
    stars.300 = "*"
  } else if (cor.test.300$p.value < 0.01 & cor.test.300$p.value > 0.001)
  {
    stars.300 = "**"
  } else if (cor.test.300$p.value < 0.001)
    stars.300 = "***"
  excluded_groups = c("Porifera")#,"Phaeodaria","Basidiomycota")
  plot.raref.300.Moran.I.only.vs.diversity = cor.plot(x = x.300,
                                                      y = y.300,
                                                      excluded.points = which(taxo_groups[selected_groups & diversity > 300] %in% excluded_groups),
                                                      x.lab = "Diversity (#OTUs)",
                                                      y.lab = "Short-distance spatial autocorrelation",
                                                      y.lab.hjust = 0.5,
                                                      x.log = T,
                                                      y.log = F,
                                                      fit = F,
                                                      mar.vect = c(5,5,1,5)) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    geom_smooth(data = data.frame(x = as.vector(x.300.2000), y = y.300.2000),aes(x,y),method='lm',se=T,colour="black") +
    annotate(geom="text", 
             x=(max(x.300)/min(x.300))^x.cor.pos*min(x.300),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos.300*(max(y.300)-min(y.300))+min(y.300),
             # label=bquote(rho==.(format(cor.test.300$estimate,digits=2,nsmall=2))^.(stars.300)),
             label=bquote(atop(rho==.(format(cor.test.300$estimate,digits=2,nsmall=2)),p==.(format(cor.test.300$p.value,digits=1)))),
             size=8)
  pdf(paste0(figure_folder,"/FigS_rarefied.300.Moran.I.only_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.300.Moran.I.only.vs.diversity)
  dev.off()
  
  ####### Moran I vs rarefied 200 diversity (residual variation around 200)
  x.300 = diversity.300[selected_groups & taxo_groups != "Porifera"]
  y.300 = I_square.observed_w.mean.300[selected_groups & taxo_groups != "Porifera",1]
  plot.raref.300.Moran.I.vs.raref.300.diversity = cor.plot(x = x.300,
                                                           y = y.300,
                                                           x.lab = "Diversity after 300-OTU rarefaction (#OTUs)",
                                                           y.lab = "Short-distance spatial autocorrelation",
                                                           y.lab.hjust = 0.5,
                                                           x.log = T,
                                                           y.log = F,
                                                           fit = T,
                                                           # fit.display = "pearson.spearman",
                                                           mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.300.Moran.I_vs_raref.300.diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.300.Moran.I.vs.raref.300.diversity)
  dev.off()
  
  #######  ####### Moran I vs rarefied 300 nb of reads
  x.300 = tot_reads.300[selected_groups & diversity>300 & taxo_groups != "Porifera"]
  y.300 = I_square.observed_w.mean.300[selected_groups  & diversity>300 & taxo_groups != "Porifera",1]
  plot.raref.Moran.I.300.vs.raref.nb.reads = cor.plot(x = x.300,
                                                      y = y.300,
                                                      x.lab = "Number of reads after 300-OTU rarefaction",
                                                      y.lab = "Short-distance spatial autocorrelation",
                                                      y.lab.hjust = 0.5,
                                                      x.log = T,
                                                      y.log = F,
                                                      fit = T,
                                                      # fit.display = "pearson.spearman",
                                                      mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.300.Moran.I_vs_raref.300.nb.reads",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.Moran.I.300.vs.raref.nb.reads)
  dev.off()
}

# Fig. S Moran I rarefied (equal diversity: 350) vs. diversity
{
  ##########
  x.cor.pos = 0.15
  y.cor.pos.350 = 0.95
  y.350 = I_square.observed_w.mean.350[selected_groups & diversity > 350,1]
  x.350 = diversity[selected_groups & diversity > 350]
  y.350.2000 = y.350[x.350 < 2000 & 
                       !names(x.350) %in% c("Porifera")]#,"MAST-4,_6,_7,_8,_9,_10,_11","Phaeodaria","Nassellaria_and_Eucyrtidium_*")]
  x.350.2000 = x.350[x.350 < 2000 & 
                       !names(x.350) %in% c("Porifera")]#,"MAST-4,_6,_7,_8,_9,_10,_11","Phaeodaria","Nassellaria_and_Eucyrtidium_*")]
  cor.test.350 = cor.test(log10(x.350.2000),y.350.2000,na.rm=T)
  if (cor.test.350$p.value > 0.05)
  {
    stars.350 = ""
  } else if (cor.test.350$p.value < 0.05 & cor.test.350$p.value > 0.01)
  {
    stars.350 = "*"
  } else if (cor.test.350$p.value < 0.01 & cor.test.350$p.value > 0.001)
  {
    stars.350 = "**"
  } else if (cor.test.350$p.value < 0.001)
    stars.350 = "***"
  excluded_groups = c("Porifera")#,"MAST-4,_6,_7,_8,_9,_10,_11","Phaeodaria","Nassellaria_and_Eucyrtidium_*")
  plot.raref.350.Moran.I.only.vs.diversity = cor.plot(x = x.350,
                                                      y = y.350,
                                                      excluded.points = which(taxo_groups[selected_groups & diversity > 350] %in% excluded_groups),
                                                      x.lab = "Diversity (#OTUs)",
                                                      y.lab = "Short-distance spatial autocorrelation",
                                                      y.lab.hjust = 0.5,
                                                      x.log = T,
                                                      y.log = F,
                                                      fit = F,
                                                      mar.vect = c(5,5,1,5)) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    geom_smooth(data = data.frame(x = as.vector(x.350.2000), y = y.350.2000),aes(x,y),method='lm',se=T,colour="black") +
    annotate(geom="text", 
             x=(max(x.350)/min(x.350))^x.cor.pos*min(x.350),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos.350*(max(y.350)-min(y.350))+min(y.350),
             # label=bquote(rho==.(format(cor.test.350$estimate,digits=2,nsmall=2))^.(stars.350)),
             label=bquote(atop(rho==.(format(cor.test.350$estimate,digits=2,nsmall=2)),p==.(format(cor.test.350$p.value,digits=1)))),
             size=8)
  pdf(paste0(figure_folder,"/FigS_rarefied.350.Moran.I.only_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.350.Moran.I.only.vs.diversity)
  dev.off()
  
  ####### Moran I vs rarefied 200 diversity (residual variation around 200)
  x.350 = diversity.350[selected_groups & taxo_groups != "Porifera"]
  y.350 = I_square.observed_w.mean.350[selected_groups & taxo_groups != "Porifera",1]
  plot.raref.350.Moran.I.vs.raref.350.diversity = cor.plot(x = x.350,
                                                           y = y.350,
                                                           x.lab = "Diversity after 350-OTU rarefaction (#OTUs)",
                                                           y.lab = "Short-distance spatial autocorrelation",
                                                           y.lab.hjust = 0.5,
                                                           x.log = T,
                                                           y.log = F,
                                                           fit = T,
                                                           # fit.display = "pearson.spearman",
                                                           mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.350.Moran.I_vs_raref.350.diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.350.Moran.I.vs.raref.350.diversity)
  dev.off()
  
  #######  ####### Moran I vs rarefied 350 nb of reads
  x.350 = tot_reads.350[selected_groups & diversity>350 & taxo_groups != "Porifera"]
  y.350 = I_square.observed_w.mean.350[selected_groups  & diversity>350 & taxo_groups != "Porifera",1]
  plot.raref.Moran.I.350.vs.raref.nb.reads = cor.plot(x = x.350,
                                                      y = y.350,
                                                      x.lab = "Number of reads after 350-OTU rarefaction",
                                                      y.lab = "Short-distance spatial autocorrelation",
                                                      y.lab.hjust = 0.5,
                                                      x.log = T,
                                                      y.log = F,
                                                      fit = T,
                                                      # fit.display = "pearson.spearman",
                                                      mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.350.Moran.I_vs_raref.350.nb.reads",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.Moran.I.350.vs.raref.nb.reads)
  dev.off()
}

# Fig. S Moran I randomly rarefied (equal diversity: 400) vs. diversity
{
  data.lme = data.frame(Groups = taxo_groups[selected_groups & diversity > 400 & diversity < 2000 & taxo_groups != "Porifera"],
                        Diversity = as.vector(diversity)[selected_groups & diversity > 400 & diversity < 2000 & taxo_groups != "Porifera"],
                        I_square.observed_w.mean.400.random.5reals[[1]][selected_groups & diversity > 400 & diversity < 2000 & taxo_groups != "Porifera",])
  melted.data.lme = reshape2::melt(data.lme,id.vars=1:2,
                                   measure.vars=3:(5+2),
                                   variable.name="Real",
                                   value.name="Moran.I")
  melted.data.lme[,"Diversity"] = log10(melted.data.lme[,"Diversity"])
  
  test = lme(Diversity ~ Moran.I, data = melted.data.lme, random = ~1|Groups)
  ######
  
  diversity_reals = as.vector(diversity)[selected_groups & diversity > 400 & diversity < 2000 & taxo_groups != "Porifera"]
  cor_reals.400 = vector(length = 5, mode = "numeric")
  for (i_real in 1:5)
  {
    Moran.I_real = I_square.observed_w.mean.400.random.5reals[[1]][selected_groups & diversity > 400 & diversity < 2000 & taxo_groups != "Porifera",i_real]
    cor_reals.400[i_real] = cor(Moran.I_real,log10(diversity_reals))
  }
    
  ######
  y.400 = I_square.observed_w.mean.400.random.5reals[[1]][selected_groups & diversity > 400 & taxo_groups != "Porifera",]
  melted.y.400 = reshape2::melt(y.400)$value
  x.400 = as.vector(diversity)[selected_groups & diversity > 400 & taxo_groups != "Porifera"]
  melted.x.400 = rep(x.400,5)
  y.400.2000 = y.400[x.400 < 2000,]
  melted.y.400.2000 = reshape2::melt(y.400.2000)$value
  x.400.2000 = x.400[x.400 < 2000]
  melted.x.400.2000 = rep(x.400.2000,5)
  #
  x.cor.pos = 0.1
  y.cor.pos = 0.95
  cor.test.400 = cor.test(log10(melted.x.400.2000),melted.y.400.2000,na.rm=T)
  if (cor.test.400$p.value > 0.05)
  {
    stars.400 = ""
  } else if (cor.test.400$p.value < 0.05 & cor.test.400$p.value > 0.01)
  {
    stars.400 = "*"
  } else if (cor.test.400$p.value < 0.01 & cor.test.400$p.value > 0.001)
  {
    stars.400 = "**"
  } else if (cor.test.400$p.value < 0.001)
    stars.400 = "***"
  #
  plot.raref.400.random.5real.Moran.I.vs.diversity = cor.plot(x = melted.x.400,
                                                              y = melted.y.400,
                                                              x.lab = "Initial diversity (#OTUs)",
                                                              y.lab = "Short-distance spatial autocorrelation",
                                                              y.lab.hjust = 0.5,
                                                              x.log = T,
                                                              y.log = F,
                                                              fit = F,
                                                              mar.vect = c(5,5,1,5)) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(melted.y.400), yend = max(melted.y.400)), linetype="dashed", inherit.aes = T) +
    geom_smooth(data = data.frame(x = melted.x.400.2000, y = melted.y.400.2000),aes(x,y),method='lm',se=T) +
    annotate(geom="text",
             x=(max(melted.x.400)/min(melted.x.400))^x.cor.pos*min(melted.x.400),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos*(max(melted.y.400)-min(melted.y.400))+min(melted.y.400),
             label=bquote(rho==.(format(cor.test.400$estimate,digits=2,nsmall=2))^.(stars.400)),
             size=8)
  pdf(paste0(figure_folder,"/FigS_rarefied.400.random.reals.Moran.I_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.400.random.5real.Moran.I.vs.diversity)
  dev.off()
}

# Fig. S Moran I rarefied (equal diversity: 500) vs. diversity
{
  ##########
  x.cor.pos = 0.15
  y.cor.pos.500 = 0.95
  y.500 = I_square.observed_w.mean.500[selected_groups & diversity > 500,1]
  x.500 = diversity[selected_groups & diversity > 500]
  y.500.2000 = y.500[x.500 < 2000 & 
                       !names(x.500) %in% c("Porifera")]#,"MAST-4,_6,_7,_8,_9,_10,_11","Phaeodaria","Nassellaria_and_Eucyrtidium_*")]
  x.500.2000 = x.500[x.500 < 2000 & 
                       !names(x.500) %in% c("Porifera")]#,"MAST-4,_6,_7,_8,_9,_10,_11","Phaeodaria","Nassellaria_and_Eucyrtidium_*")]
  cor.test.500 = cor.test(log10(x.500.2000),y.500.2000,na.rm=T)
  if (cor.test.500$p.value > 0.05)
  {
    stars.500 = ""
  } else if (cor.test.500$p.value < 0.05 & cor.test.500$p.value > 0.01)
  {
    stars.500 = "*"
  } else if (cor.test.500$p.value < 0.01 & cor.test.500$p.value > 0.001)
  {
    stars.500 = "**"
  } else if (cor.test.500$p.value < 0.001)
    stars.500 = "***"
  excluded_groups = c("Porifera")#,"MAST-4,_6,_7,_8,_9,_10,_11","Phaeodaria","Nassellaria_and_Eucyrtidium_*")
  plot.raref.500.Moran.I.only.vs.diversity = cor.plot(x = x.500,
                                                      y = y.500,
                                                      excluded.points = which(taxo_groups[selected_groups & diversity > 500] %in% excluded_groups),
                                                      x.lab = "Diversity (#OTUs)",
                                                      y.lab = "Short-distance spatial autocorrelation",
                                                      y.lab.hjust = 0.5,
                                                      x.log = T,
                                                      y.log = F,
                                                      fit = F,
                                                      mar.vect = c(5,5,1,5)) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    geom_smooth(data = data.frame(x = as.vector(x.500.2000), y = y.500.2000),aes(x,y),method='lm',se=T,colour="black") +
    annotate(geom="text", 
             x=(max(x.500)/min(x.500))^x.cor.pos*min(x.500),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos.500*(max(y.500)-min(y.500))+min(y.500),
             # label=bquote(rho==.(format(cor.test.500$estimate,digits=2,nsmall=2))^.(stars.500)),
             label=bquote(atop(rho==.(format(cor.test.500$estimate,digits=2,nsmall=2)),p==.(format(cor.test.500$p.value,digits=1)))),
             size=8)
  pdf(paste0(figure_folder,"/FigS_rarefied.500.Moran.I.only_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.500.Moran.I.only.vs.diversity)
  dev.off()
  
  ####### Moran I vs rarefied 500 diversity (residual variation around 500)
  x.500 = diversity.500[selected_groups & taxo_groups != "Porifera"]
  y.500 = I_square.observed_w.mean.500[selected_groups & taxo_groups != "Porifera",1]
  plot.raref.500.Moran.I.vs.raref.500.diversity = cor.plot(x = x.500,
                                                           y = y.500,
                                                           x.lab = "Diversity after 500-OTU rarefaction (#OTUs)",
                                                           y.lab = "Short-distance spatial autocorrelation",
                                                           y.lab.hjust = 0.5,
                                                           x.log = T,
                                                           y.log = F,
                                                           fit = T,
                                                           # fit.display = "pearson.spearman",
                                                           mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.500.Moran.I_vs_raref.500.diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.500.Moran.I.vs.raref.500.diversity)
  dev.off()
  
  #######  ####### Moran I vs rarefied 500 nb of reads
  x.500 = tot_reads.500[selected_groups & diversity>500 & taxo_groups != "Porifera"]
  y.500 = I_square.observed_w.mean.500[selected_groups  & diversity>500 & taxo_groups != "Porifera",1]
  plot.raref.Moran.I.500.vs.raref.nb.reads = cor.plot(x = x.500,
                                                      y = y.500,
                                                      x.lab = "Number of reads after 500-OTU rarefaction",
                                                      y.lab = "Short-distance spatial autocorrelation",
                                                      y.lab.hjust = 0.5,
                                                      x.log = T,
                                                      y.log = F,
                                                      fit = T,
                                                      # fit.display = "pearson.spearman",
                                                      mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.500.Moran.I_vs_raref.500.nb.reads",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.Moran.I.500.vs.raref.nb.reads)
  dev.off()
}

# Fig. S Moran I randomly rarefied (equal diversity: 250) vs. diversity
{
  ##########
  x.cor.pos = 0.15
  y.cor.pos.250.random = 0.95
  y.250.random = I_square.observed_w.mean.250.random[selected_groups & diversity > 250,1]
  x.250.random = diversity[selected_groups & diversity > 250]
  y.250.random.2000 = y.250.random[x.250.random < 2000 & 
                       !names(x.250.random) %in% c("Porifera")]#,"MAST-4,_6,_7,_8,_9,_10,_11","Phaeodaria","Nassellaria_and_Eucyrtidium_*")]
  x.250.random.2000 = x.250.random[x.250.random < 2000 & 
                       !names(x.250.random) %in% c("Porifera")]#,"MAST-4,_6,_7,_8,_9,_10,_11","Phaeodaria","Nassellaria_and_Eucyrtidium_*")]
  cor.test.250.random = cor.test(log10(x.250.random.2000),y.250.random.2000,na.rm=T)
  if (cor.test.250.random$p.value > 0.05)
  {
    stars.250.random = ""
  } else if (cor.test.250.random$p.value < 0.05 & cor.test.250.random$p.value > 0.01)
  {
    stars.250.random = "*"
  } else if (cor.test.250.random$p.value < 0.01 & cor.test.250.random$p.value > 0.001)
  {
    stars.250.random = "**"
  } else if (cor.test.250.random$p.value < 0.001)
    stars.250.random = "***"
  excluded_groups = c("Porifera")#,"MAST-4,_6,_7,_8,_9,_10,_11","Phaeodaria","Nassellaria_and_Eucyrtidium_*")
  plot.raref.250.random.Moran.I.only.vs.diversity = cor.plot(x = x.250.random,
                                                      y = y.250.random,
                                                      excluded.points = which(taxo_groups[selected_groups & diversity > 250] %in% excluded_groups),
                                                      x.lab = "Diversity (#OTUs)",
                                                      y.lab = "Short-distance spatial autocorrelation",
                                                      y.lab.hjust = 0.5,
                                                      x.log = T,
                                                      y.log = F,
                                                      fit = F,
                                                      mar.vect = c(5,5,1,5)) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    geom_smooth(data = data.frame(x = as.vector(x.250.random.2000), y = y.250.random.2000),aes(x,y),method='lm',se=T,colour="black") +
    annotate(geom="text", 
             x=(max(x.250.random)/min(x.250.random))^x.cor.pos*min(x.250.random),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos.250.random*(max(y.250.random)-min(y.250.random))+min(y.250.random),
             # label=bquote(rho==.(format(cor.test.250.random$estimate,digits=2,nsmall=2))^.(stars.250.random)),
             label=bquote(atop(rho==.(format(cor.test.250.random$estimate,digits=2,nsmall=2)),p==.(format(cor.test.250.random$p.value,digits=1)))),
             size=8)
  pdf(paste0(figure_folder,"/FigS_rarefied.250.random.Moran.I.only_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.250.random.Moran.I.only.vs.diversity)
  dev.off()
}

# Fig. S Moran I randomly rarefied (equal diversity: 350) vs. diversity
{
  ##########
  x.cor.pos = 0.15
  y.cor.pos.350.random = 0.95
  y.350.random = I_square.observed_w.mean.350.random[selected_groups & diversity > 350,1]
  x.350.random = diversity[selected_groups & diversity > 350]
  y.350.random.2000 = y.350.random[x.350.random < 2000 & 
                                     !names(x.350.random) %in% c("Porifera")]#,"MAST-4,_6,_7,_8,_9,_10,_11","Phaeodaria","Nassellaria_and_Eucyrtidium_*")]
  x.350.random.2000 = x.350.random[x.350.random < 2000 & 
                                     !names(x.350.random) %in% c("Porifera")]#,"MAST-4,_6,_7,_8,_9,_10,_11","Phaeodaria","Nassellaria_and_Eucyrtidium_*")]
  cor.test.350.random = cor.test(log10(x.350.random.2000),y.350.random.2000,na.rm=T)
  if (cor.test.350.random$p.value > 0.05)
  {
    stars.350.random = ""
  } else if (cor.test.350.random$p.value < 0.05 & cor.test.350.random$p.value > 0.01)
  {
    stars.350.random = "*"
  } else if (cor.test.350.random$p.value < 0.01 & cor.test.350.random$p.value > 0.001)
  {
    stars.350.random = "**"
  } else if (cor.test.350.random$p.value < 0.001)
    stars.350.random = "***"
  excluded_groups = c("Porifera")#,"MAST-4,_6,_7,_8,_9,_10,_11","Phaeodaria","Nassellaria_and_Eucyrtidium_*")
  plot.raref.350.random.Moran.I.only.vs.diversity = cor.plot(x = x.350.random,
                                                             y = y.350.random,
                                                             excluded.points = which(taxo_groups[selected_groups & diversity > 350] %in% excluded_groups),
                                                             x.lab = "Diversity (#OTUs)",
                                                             y.lab = "Short-distance spatial autocorrelation",
                                                             y.lab.hjust = 0.5,
                                                             x.log = T,
                                                             y.log = F,
                                                             fit = F,
                                                             mar.vect = c(5,5,1,5)) +
    geom_segment(aes(x = 2000, xend = 2000, y = min(y), yend = max(y)), linetype="dashed", inherit.aes = T) +
    geom_smooth(data = data.frame(x = as.vector(x.350.random.2000), y = y.350.random.2000),aes(x,y),method='lm',se=T,colour="black") +
    annotate(geom="text", 
             x=(max(x.350.random)/min(x.350.random))^x.cor.pos*min(x.350.random),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos.350.random*(max(y.350.random)-min(y.350.random))+min(y.350.random),
             # label=bquote(rho==.(format(cor.test.350.random$estimate,digits=2,nsmall=2))^.(stars.350.random)),
             label=bquote(atop(rho==.(format(cor.test.350.random$estimate,digits=2,nsmall=2)),p==.(format(cor.test.350.random$p.value,digits=1)))),
             size=8)
  pdf(paste0(figure_folder,"/FigS_rarefied.350.random.Moran.I.only_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.350.random.Moran.I.only.vs.diversity)
  dev.off()
}

# Fig. S mean OTU prevalence in rarefied data (equal diversity: 200-250-500-1000) vs. initial diversity
{
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera" & diversity > 200]
  y = mean_prevalence.200[selected_groups & taxo_groups != "Porifera" & diversity > 200]
  plot.raref.200.mean.prevalence.vs.diversity = cor.plot(x = x,
                                                         y = y,
                                                         x.lab = "Diversity (#OTUs)",
                                                         y.lab = "Mean OTU prevalence",
                                                         y.lab.hjust = 0.5,
                                                         x.log = T,
                                                         y.log = F,
                                                         fit = T,
                                                         fit.display = "pearson.spearman",
                                                         mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.200.mean.prevalence_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.200.mean.prevalence.vs.diversity)
  dev.off()
  
  ###
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera" & diversity > 250]
  y = mean_prevalence.250[selected_groups & taxo_groups != "Porifera" & diversity > 250]
  plot.raref.250.mean.prevalence.vs.diversity = cor.plot(x = x,
                                                         y = y,
                                                         x.lab = "Diversity (#OTUs)",
                                                         y.lab = "Mean OTU prevalence",
                                                         y.lab.hjust = 0.5,
                                                         x.log = T,
                                                         y.log = F,
                                                         fit = T,
                                                         fit.display = "pearson.spearman",
                                                         mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.250.mean.prevalence_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.250.mean.prevalence.vs.diversity)
  dev.off()
  
  ###
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera" & diversity > 250]
  y = mean_prevalence.250.random[selected_groups & taxo_groups != "Porifera" & diversity > 250]
  plot.raref.250.random.mean.prevalence.vs.diversity = cor.plot(x = x,
                                                         y = y,
                                                         x.lab = "Diversity (#OTUs)",
                                                         y.lab = "Mean OTU prevalence",
                                                         y.lab.hjust = 0.5,
                                                         x.log = T,
                                                         y.log = F,
                                                         fit = T,
                                                         fit.display = "pearson.spearman",
                                                         mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.250.random.mean.prevalence_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.250.random.mean.prevalence.vs.diversity)
  dev.off()
  
  ###
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera" & diversity > 350]
  y = mean_prevalence.350[selected_groups & taxo_groups != "Porifera" & diversity > 350]
  plot.raref.350.mean.prevalence.vs.diversity = cor.plot(x = x,
                                                         y = y,
                                                         x.lab = "Diversity (#OTUs)",
                                                         y.lab = "Mean OTU prevalence",
                                                         y.lab.hjust = 0.5,
                                                         x.log = T,
                                                         y.log = F,
                                                         fit = T,
                                                         fit.display = "pearson.spearman",
                                                         mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.350.mean.prevalence_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.350.mean.prevalence.vs.diversity)
  dev.off()
  
  ###
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera" & diversity > 350]
  y = mean_prevalence.350.random[selected_groups & taxo_groups != "Porifera" & diversity > 350]
  plot.raref.350.random.mean.prevalence.vs.diversity = cor.plot(x = x,
                                                                y = y,
                                                                x.lab = "Diversity (#OTUs)",
                                                                y.lab = "Mean OTU prevalence",
                                                                y.lab.hjust = 0.5,
                                                                x.log = T,
                                                                y.log = F,
                                                                fit = T,
                                                                fit.display = "pearson.spearman",
                                                                mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.350.random.mean.prevalence_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.350.random.mean.prevalence.vs.diversity)
  dev.off()
  
  ###
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera" & diversity > 500]
  y = mean_prevalence.500[selected_groups & taxo_groups != "Porifera" & diversity > 500]
  plot.raref.500.mean.prevalence.vs.diversity = cor.plot(x = x,
                                                         y = y,
                                                         x.lab = "Diversity (#OTUs)",
                                                         y.lab = "Mean OTU prevalence",
                                                         y.lab.hjust = 0.5,
                                                         x.log = T,
                                                         y.log = F,
                                                         fit = T,
                                                         fit.display = "pearson.spearman",
                                                         mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.500.mean.prevalence_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.500.mean.prevalence.vs.diversity)
  dev.off()
  
  ###
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera" & diversity > 1000]
  y = mean_prevalence.1000[selected_groups & taxo_groups != "Porifera" & diversity > 1000]
  plot.raref.1000.mean.prevalence.vs.diversity = cor.plot(x = x,
                                                         y = y,
                                                         x.lab = "Diversity (#OTUs)",
                                                         y.lab = "Mean OTU prevalence",
                                                         y.lab.hjust = 0.5,
                                                         x.log = T,
                                                         y.log = F,
                                                         fit = T,
                                                         fit.display = "pearson.spearman",
                                                         mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.1000.mean.prevalence_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.1000.mean.prevalence.vs.diversity)
  dev.off()
  
  ###
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera"]
  y = mean_prevalence[selected_groups & taxo_groups != "Porifera"]
  plot.mean.prevalence.vs.diversity = cor.plot(x = x,
                                               y = y,
                                               x.lab = "Diversity (#OTUs)",
                                               y.lab = "Mean OTU prevalence",
                                               y.lab.hjust = 0.5,
                                               x.log = T,
                                               y.log = F,
                                               fit = T,
                                               fit.display = "pearson.spearman",
                                               mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_mean.prevalence_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.mean.prevalence.vs.diversity)
  dev.off()
}

# Fig. S median OTU prevalence in rarefied data (equal diversity: 200-250-500-1000) vs. initial diversity
{
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera" & diversity > 200]
  y = median_prevalence.200[selected_groups & taxo_groups != "Porifera" & diversity > 200]
  plot.raref.200.median.prevalence.vs.diversity = cor.plot(x = x,
                                                         y = y,
                                                         x.lab = "Diversity (#OTUs)",
                                                         y.lab = "Median OTU prevalence",
                                                         y.lab.hjust = 0.5,
                                                         x.log = T,
                                                         y.log = F,
                                                         fit = T,
                                                         fit.display = "pearson.spearman",
                                                         mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.200.median.prevalence_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.200.median.prevalence.vs.diversity)
  dev.off()
  
  ###
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera" & diversity > 250]
  y = median_prevalence.250[selected_groups & taxo_groups != "Porifera" & diversity > 250]
  plot.raref.250.median.prevalence.vs.diversity = cor.plot(x = x,
                                                           y = y,
                                                           x.lab = "Diversity (#OTUs)",
                                                           y.lab = "Median OTU prevalence",
                                                           y.lab.hjust = 0.5,
                                                           x.log = T,
                                                           y.log = F,
                                                           fit = T,
                                                           fit.display = "pearson.spearman",
                                                           mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.250.median.prevalence_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.250.median.prevalence.vs.diversity)
  dev.off()
  
  ###
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera" & diversity > 250]
  y = median_prevalence.250.random[selected_groups & taxo_groups != "Porifera" & diversity > 250]
  plot.raref.250.random.median.prevalence.vs.diversity = cor.plot(x = x,
                                                                  y = y,
                                                                  x.lab = "Diversity (#OTUs)",
                                                                  y.lab = "Median OTU prevalence",
                                                                  y.lab.hjust = 0.5,
                                                                  x.log = T,
                                                                  y.log = F,
                                                                  fit = T,
                                                                  fit.display = "pearson.spearman",
                                                                  mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.250.random.median.prevalence_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.250.random.median.prevalence.vs.diversity)
  dev.off()
  
  ###
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera" & diversity > 350]
  y = median_prevalence.350[selected_groups & taxo_groups != "Porifera" & diversity > 350]
  plot.raref.350.median.prevalence.vs.diversity = cor.plot(x = x,
                                                           y = y,
                                                           x.lab = "Diversity (#OTUs)",
                                                           y.lab = "Median OTU prevalence",
                                                           y.lab.hjust = 0.5,
                                                           x.log = T,
                                                           y.log = F,
                                                           fit = T,
                                                           fit.display = "pearson.spearman",
                                                           mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.350.median.prevalence_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.350.median.prevalence.vs.diversity)
  dev.off()
  
  ###
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera" & diversity > 350]
  y = median_prevalence.350.random[selected_groups & taxo_groups != "Porifera" & diversity > 350]
  plot.raref.350.random.median.prevalence.vs.diversity = cor.plot(x = x,
                                                                y = y,
                                                                x.lab = "Diversity (#OTUs)",
                                                                y.lab = "Median OTU prevalence",
                                                                y.lab.hjust = 0.5,
                                                                x.log = T,
                                                                y.log = F,
                                                                fit = T,
                                                                fit.display = "pearson.spearman",
                                                                mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.350.random.median.prevalence_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.350.random.median.prevalence.vs.diversity)
  dev.off()
  
  ###
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera" & diversity > 500]
  y = median_prevalence.500[selected_groups & taxo_groups != "Porifera" & diversity > 500]
  plot.raref.500.median.prevalence.vs.diversity = cor.plot(x = x,
                                                         y = y,
                                                         x.lab = "Diversity (#OTUs)",
                                                         y.lab = "Median OTU prevalence",
                                                         y.lab.hjust = 0.5,
                                                         x.log = T,
                                                         y.log = F,
                                                         fit = T,
                                                         fit.display = "pearson.spearman",
                                                         mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.500.mean.prevalence_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.500.median.prevalence.vs.diversity)
  dev.off()
  
  ###
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera" & diversity > 1000]
  y = mean_prevalence.1000[selected_groups & taxo_groups != "Porifera" & diversity > 1000]
  plot.raref.1000.median.prevalence.vs.diversity = cor.plot(x = x,
                                                          y = y,
                                                          x.lab = "Diversity (#OTUs)",
                                                          y.lab = "Median OTU prevalence",
                                                          y.lab.hjust = 0.5,
                                                          x.log = T,
                                                          y.log = F,
                                                          fit = T,
                                                          fit.display = "pearson.spearman",
                                                          mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_rarefied.1000.median.prevalence_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.raref.1000.median.prevalence.vs.diversity)
  dev.off()
  
  ###
  x = as.vector(diversity)[selected_groups & taxo_groups != "Porifera"]
  y = median_prevalence[selected_groups & taxo_groups != "Porifera"]
  plot.median.prevalence.vs.diversity = cor.plot(x = x,
                                               y = y,
                                               x.lab = "Diversity (#OTUs)",
                                               y.lab = "Median OTU prevalence",
                                               y.lab.hjust = 0.5,
                                               x.log = T,
                                               y.log = F,
                                               fit = T,
                                               fit.display = "pearson.spearman",
                                               mar.vect = c(5,5,1,5))
  pdf(paste0(figure_folder,"/FigS_median.prevalence_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.median.prevalence.vs.diversity)
  dev.off()
}

# Fig. S random raref. summary
{
  data.lme.200 = data.frame(Groups = taxo_groups[selected_groups & diversity > 200 & diversity < 2000],
                            Diversity = as.vector(diversity)[selected_groups & diversity > 200 & diversity < 2000],
                            I_square.observed_w.mean.200.random.10reals[[1]][selected_groups & diversity > 200 & diversity < 2000,])
  melted.data.lme.200 = reshape2::melt(data.lme.200,id.vars=1:2,
                                       measure.vars=3:(10+2),
                                       variable.name="Real",
                                       value.name="Moran.I")
  melted.data.lme.200[,"Diversity"] = log10(melted.data.lme.200[,"Diversity"])
  test.200 = lme(Diversity ~ Moran.I, data = melted.data.lme.200, random = ~1|Groups)
  ######
  data.lme.400 = data.frame(Groups = taxo_groups[selected_groups & diversity > 400 & diversity < 2000],
                            Diversity = as.vector(diversity)[selected_groups & diversity > 400 & diversity < 2000],
                            I_square.observed_w.mean.400.random.10reals[[1]][selected_groups & diversity > 400 & diversity < 2000,])
  melted.data.lme.400 = reshape2::melt(data.lme.400,id.vars=1:2,
                                       measure.vars=3:(10+2),
                                       variable.name="Real",
                                       value.name="Moran.I")
  melted.data.lme.400[,"Diversity"] = log10(melted.data.lme.400[,"Diversity"])
  test.400 = lme(Diversity ~ Moran.I, data = melted.data.lme.400, random = ~1|Groups)
  ######
  data.lme.800 = data.frame(Groups = taxo_groups[selected_groups & diversity > 800 & diversity < 2000],
                            Diversity = as.vector(diversity)[selected_groups & diversity > 800 & diversity < 2000],
                            I_square.observed_w.mean.800.random.10reals[[1]][selected_groups & diversity > 800 & diversity < 2000,])
  melted.data.lme.800 = reshape2::melt(data.lme.800,id.vars=1:2,
                                       measure.vars=3:(10+2),
                                       variable.name="Real",
                                       value.name="Moran.I")
  melted.data.lme.800[,"Diversity"] = log10(melted.data.lme.800[,"Diversity"])
  test.800 = lme(Diversity ~ Moran.I, data = melted.data.lme.800, random = ~1|Groups)
  
  ##############
  y.200 = I_square.observed_w.mean.200.random.10reals[[1]][selected_groups & diversity > 200] #& taxo_groups != "Porifera",]
  melted.y.200 = reshape2::melt(y.200)$value
  x.200 = as.vector(diversity)[selected_groups & diversity > 200] #& taxo_groups != "Porifera"]
  melted.x.200 = rep(x.200,10)
  y.200.2000 = y.200[x.200 < 2000,]
  melted.y.200.2000 = reshape2::melt(y.200.2000)$value
  x.200.2000 = x.200[x.200 < 2000]
  melted.x.200.2000 = rep(x.200.2000,10)
  #
  x.cor.pos = 0.55
  y.cor.pos = 0.1
  cor.test.200 = cor.test(log10(melted.x.200.2000),melted.y.200.2000,na.rm=T,method="spearman")
  plot.raref.200.random.10real.Moran.I.vs.diversity = cor.plot(x = melted.x.200,
                                                               y = melted.y.200,
                                                               x.lab = "Original diversity (#OTUs)",
                                                               y.lab = "Surface short-distance spatial autocorr.",
                                                               y.lab.hjust = 0.5,
                                                               x.log = T,
                                                               y.log = F,
                                                               fit = F,
                                                               mar.vect = c(5,5,1,5)) +
    geom_vline(xintercept = 2000, linetype="dashed") +
    geom_smooth(data = data.frame(x = melted.x.200.2000, y = melted.y.200.2000),aes(x,y),method='lm',color="black",se=T) +
    annotate(geom="text",
             x=(max(melted.x.200)/min(melted.x.200))^x.cor.pos*min(melted.x.200),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos*(max(melted.y.200)-min(melted.y.200))+min(melted.y.200),
             label=bquote(rho[S]==.(format(cor.test.200$estimate,digits=2,nsmall=2))),#,P==.(format(cor.test.200$p.value,digits=1)))),
             size=7.5)
  # pdf(paste0(figure_folder,"/FigS_rarefied.200.random.reals.Moran.I_vs_diversity",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
  #     width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  # print(plot.raref.200.random.10real.Moran.I.vs.diversity)
  # dev.off()
  
  ##############
  y.400 = I_square.observed_w.mean.400.random.10reals[[1]][selected_groups & diversity > 400] #& taxo_groups != "Porifera",]
  melted.y.400 = reshape2::melt(y.400)$value
  x.400 = as.vector(diversity)[selected_groups & diversity > 400] #& taxo_groups != "Porifera"]
  melted.x.400 = rep(x.400,10)
  y.400.2000 = y.400[x.400 < 2000,]
  melted.y.400.2000 = reshape2::melt(y.400.2000)$value
  x.400.2000 = x.400[x.400 < 2000]
  melted.x.400.2000 = rep(x.400.2000,10)
  #
  x.cor.pos = 0.5
  y.cor.pos = 0.1
  cor.test.400 = cor.test(log10(melted.x.400.2000),melted.y.400.2000,na.rm=T,method="spearman")
  plot.raref.400.random.10real.Moran.I.vs.diversity = cor.plot(x = melted.x.400,
                                                               y = melted.y.400,
                                                               x.lab = "Original diversity (#OTUs)",
                                                               y.lab = "Surface short-distance spatial autocorr.",
                                                               y.lab.hjust = 0.5,
                                                               x.log = T,
                                                               y.log = F,
                                                               fit = F,
                                                               mar.vect = c(5,5,1,5)) +
    geom_vline(xintercept = 2000, linetype="dashed") +
    geom_smooth(data = data.frame(x = melted.x.400.2000, y = melted.y.400.2000),aes(x,y),method='lm',color="black",se=T) +
    annotate(geom="text",
             x=(max(melted.x.400)/min(melted.x.400))^x.cor.pos*min(melted.x.400),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos*(max(melted.y.400)-min(melted.y.400))+min(melted.y.400),
             label=bquote(rho[S]==.(format(cor.test.400$estimate,digits=2,nsmall=2))),#,P==.(format(cor.test.400$p.value,digits=1)))),
             size=7.5)
  
  ##############
  y.800 = I_square.observed_w.mean.800.random.10reals[[1]][selected_groups & diversity > 800] #& taxo_groups != "Porifera",]
  melted.y.800 = reshape2::melt(y.800)$value
  x.800 = as.vector(diversity)[selected_groups & diversity > 800] #& taxo_groups != "Porifera"]
  melted.x.800 = rep(x.800,10)
  y.800.2000 = y.800[x.800 < 2000,]
  melted.y.800.2000 = reshape2::melt(y.800.2000)$value
  x.800.2000 = x.800[x.800 < 2000]
  melted.x.800.2000 = rep(x.800.2000,10)
  #
  x.cor.pos = 0.4
  y.cor.pos = 0.1
  cor.test.800 = cor.test(log10(melted.x.800.2000),melted.y.800.2000,na.rm=T,method="spearman")
  plot.raref.800.random.10real.Moran.I.vs.diversity = cor.plot(x = melted.x.800,
                                                               y = melted.y.800,
                                                               x.lab = "Original diversity (#OTUs)",
                                                               y.lab = "Surface short-distance spatial autocorr.",
                                                               y.lab.hjust = 0.5,
                                                               x.log = T,
                                                               y.log = F,
                                                               fit = F,
                                                               mar.vect = c(5,5,1,5)) +
    geom_vline(xintercept = 2000, linetype="dashed") +
    geom_smooth(data = data.frame(x = melted.x.800.2000, y = melted.y.800.2000),aes(x,y),method='lm',color="black",se=T) +
    annotate(geom="text",
             x=(max(melted.x.800)/min(melted.x.800))^x.cor.pos*min(melted.x.800),
             # x=x.cor.pos*(max(x)-min(x))+min(x),
             y=y.cor.pos*(max(melted.y.800)-min(melted.y.800))+min(melted.y.800),
             label=bquote(rho[S]==.(format(cor.test.800$estimate,digits=2,nsmall=2))),#,P==.(format(cor.test.800$p.value,digits=1)))),
             size=7.5)
  
  ########
  # varpart.vs.bodysize.reals.plot = ggplot(data = data.frame(x = reshape2::melt(body.size.reals[selected_groups,])$value,
  #                                                           y = rep((colSums(varpart.env.spatial[[1]][2:3,])/(colSums(varpart.env.spatial[[1]][1:2,])))[selected_groups],n_real))) + 
  #   # geom_bin2d(aes(x,y)) +
  #   stat_binhex(aes(x,y), show.legend=F) +
  #   scale_fill_viridis() +
  #   # scale_fill_gradient(name = "count", trans = "log") +
  #   labs(y = "Surface connectivity-environment\n explained variance ratio", 
  #        x = expression("Sampled body size ("*mu*"m), log scale")) +
  #   theme_bw() +
  #   theme(axis.title=element_text(size=22),
  #         axis.text = element_text(size=22),
  #         plot.title=element_text(hjust=0, size=24),
  #         plot.margin=unit(c(5,5,5,5),"mm")) +
  #   scale_x_log10() +
  #   geom_smooth(aes(x,y), method='lm', col = "black", se=F)
  
  ##########
  diversity_reals.200 = as.vector(diversity)[selected_groups & diversity > 200 & diversity < 2000]# & taxo_groups != "Porifera"]
  diversity_reals.400 = as.vector(diversity)[selected_groups & diversity > 400 & diversity < 2000]# & taxo_groups != "Porifera"]
  diversity_reals.800 = as.vector(diversity)[selected_groups & diversity > 800 & diversity < 2000]# & taxo_groups != "Porifera"]
  cor_reals.200 = cor_reals.400 = cor_reals.800 = vector(length = 10, mode = "numeric")
  for (i_real in 1:10)
  {
    Moran.I_real = I_square.observed_w.mean.200.random.10reals[[1]][selected_groups & diversity > 200 & diversity < 2000,i_real] #& taxo_groups != "Porifera",i_real]
    cor_reals.200[i_real] = cor(Moran.I_real,log10(diversity_reals.200),method = "spearman")
    Moran.I_real = I_square.observed_w.mean.400.random.10reals[[1]][selected_groups & diversity > 400 & diversity < 2000,i_real] #& taxo_groups != "Porifera",i_real]
    cor_reals.400[i_real] = cor(Moran.I_real,log10(diversity_reals.400),method = "spearman")
    Moran.I_real = I_square.observed_w.mean.800.random.10reals[[1]][selected_groups & diversity > 800 & diversity < 2000,i_real] #& taxo_groups != "Porifera",i_real]
    cor_reals.800[i_real] = cor(Moran.I_real,log10(diversity_reals.800),method = "spearman")
  }
  
  cor.spear.raref.div.hist.plot = ggplot(data = data.frame(x = c(cor_reals.200, cor_reals.400, cor_reals.800),
                                                           raref = rep(c("200","400","800"),each=10))) +
    geom_histogram(aes(x, fill = factor(raref)),
                   position = "stack",
                   binwidth = 0.05, center = 0.025) +
    scale_fill_manual(values = muted(c("red","blue","green"),l=50,c=80),
                      name = "Rarefied to:",
                      labels = c("200 OTUs","400 OTUs","800 OTUs")) +
    labs(y="Number of independent rarefactions", x="Spearman corr. btw. surf. spatial autocorr.\nand original log diversity") +
    theme_bw() +
    theme(axis.title=element_text(size=22,hjust=0.5),
          axis.text = element_text(size=22),
          legend.text = element_text(size=20),
          legend.title = element_text(size=20),
          legend.position = c(0.70,0.80),
          legend.box.background = element_rect(color="black", size=1.5),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(5,5,5,5),"mm")) +
    geom_vline(xintercept = 0, linetype = "dashed")
  
  # pdf(paste0(figure_folder,"/FigS_rarefs.histogram.pdf"),
  #     width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  # print(cor.spear.raref.div.hist.plot)
  # dev.off()
  
  ######################
  plot.Moran.I.raref.vs.non.raref = list()
  for (i_case in 1:2)
  {
    y.200 = I_square.observed_w.mean.200.random.10reals[[i_case]][selected_groups & diversity > 200,]
    melted.y.200 = reshape2::melt(y.200)$value
    x.200 = I_square.observed_w.mean[selected_groups & diversity > 200,i_case]
    melted.x.200 = rep(x.200,10)
    y.400 = I_square.observed_w.mean.400.random.10reals[[i_case]][selected_groups & diversity > 400,]
    melted.y.400 = reshape2::melt(y.400)$value
    x.400 = I_square.observed_w.mean[selected_groups & diversity > 400,i_case]
    melted.x.400 = rep(x.400,10)
    y.800 = I_square.observed_w.mean.800.random.10reals[[i_case]][selected_groups & diversity > 800,]
    melted.y.800 = reshape2::melt(y.800)$value
    x.800 = I_square.observed_w.mean[selected_groups & diversity > 800,i_case]
    melted.x.800 = rep(x.800,10)
    y = c(melted.y.200,melted.y.400,melted.y.800)
    x = c(melted.x.200,melted.x.400+0.006,melted.x.800+0.012)
    raref = c(rep("200 OTUs",length(melted.x.200)),rep("400 OTUs",length(melted.x.400)),rep("800 OTUs",length(melted.x.800)))
    col.vect = muted(c("red","blue","green"),l=50,c=80)
    
    cor.test.200 = cor.test(melted.x.200,melted.y.200)
    cor.test.400 = cor.test(melted.x.400,melted.y.400)
    cor.test.800 = cor.test(melted.x.800,melted.y.800)
    plot.Moran.I.raref.vs.non.raref[[i_case]] = cor.plot(x = x,
                                                         y = y,
                                                         col.factor = raref,
                                                         col.vect = col.vect,
                                                         leg.name = NULL,
                                                         x.lab = "Short-distance spatial autocorrelation\n before rarefaction",
                                                         y.lab = "Short-distance spatial autocorrelation\n after rarefaction",
                                                         y.lab.hjust = 0.5,
                                                         x.log = F,
                                                         y.log = F,
                                                         fit = F,
                                                         fit.display="pearson.p",
                                                         x.cor.pos=0.35,
                                                         y.cor.pos=0.85,
                                                         mar.vect = c(5,5,1,5)) +
      # geom_smooth(aes(x,y),method='lm',color="black",se=F,size=0.7) +
      geom_smooth(data=data.frame(melted.x.200,melted.y.200), aes(melted.x.200,melted.y.200),method='lm',color=col.vect[1],se=T,size=0.5) +
      geom_smooth(data=data.frame(melted.x.400,melted.y.400), aes(melted.x.400,melted.y.400),method='lm',color=col.vect[2],se=T,size=0.5) +
      geom_smooth(data=data.frame(melted.x.800,melted.y.800), aes(melted.x.800,melted.y.800),method='lm',color=col.vect[3],se=T,size=0.5) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      annotate(geom="text",
               x=0.4*(max(x)-min(x))+min(x),
               y=0.87*(max(y)-min(y))+min(y),
               label=bquote(rho[P]==.(format(cor.test.200$estimate,digits=2,nsmall=2))),#,P==.(format(cor.test.800$p.value,digits=1)))),
               color = col.vect[1],
               size=7.5) +
      annotate(geom="text",
               x=0.4*(max(x)-min(x))+min(x),
               y=0.93*(max(y)-min(y))+min(y),
               label=bquote(rho[P]==.(format(cor.test.400$estimate,digits=2,nsmall=2))),#,P==.(format(cor.test.800$p.value,digits=1)))),
               color = col.vect[2],
               size=7.5) +
    annotate(geom="text",
             x=0.4*(max(x)-min(x))+min(x),
             y=0.99*(max(y)-min(y))+min(y),
             label=bquote(rho[P]==.(format(cor.test.800$estimate,digits=2,nsmall=2))),#,P==.(format(cor.test.800$p.value,digits=1)))),
             color = col.vect[3],
             size=7.5)
  }
  # pdf(paste0(figure_folder,"/FigS_Moran.I.raref.vs.nonraref.pdf"),
  #     width=7.5*2.2/1.2/2*2,height=12/3*4/1.2/2,onefile=F)
  # g1 = ggplotGrob(plot.Moran.I.raref.vs.non.raref[[1]])
  # g2 = ggplotGrob(plot.Moran.I.raref.vs.non.raref[[2]])
  # grid.arrange(grobs = list(g1,g2),
  #              nrow = 1)
  # dev.off()
  ######################
  
  pdf(paste0(figure_folder,"/FigS_Moran.diversity.rarefs.sensitivity.pdf"),
      width=7.5*2.2/1.2/2*2,height=12/3*2/1.2*3,onefile=F)
  # g1 = ggplotGrob(plot.random.groups.Moran.I.vs.diversity)
  g1 = ggplotGrob(plot.raref.200.random.10real.Moran.I.vs.diversity)
  g2 = ggplotGrob(plot.raref.400.random.10real.Moran.I.vs.diversity + theme(axis.title.y=element_blank()))
  g3 = ggplotGrob(plot.raref.800.random.10real.Moran.I.vs.diversity)
  g4 = ggplotGrob(cor.spear.raref.div.hist.plot)
  g5 = ggplotGrob(plot.Moran.I.raref.vs.non.raref[[1]])
  g6 = ggplotGrob(plot.Moran.I.raref.vs.non.raref[[2]] + theme(axis.title.y=element_blank()))
  grid.arrange(grobs = list(g1,g2,g3,g4,g5,g6),
               layout_matrix = rbind(c(1, 2),c(3, 4),c(5, 6)))
  dev.off()
}

######################################## Sensitivity to body size and inference settings:

# Body size sensitivity
{
  mean_sizes = c(mean(c(mean(c(0.8,3)),mean(c(0.8,5)))),mean(c(mean(c(3,20)),mean(c(5,20)))),mean(c(20,180)),mean(c(180,2000)))
  fraction.sizes = matrix(nrow = 4, ncol = 2, data = c(0.8,5,5,20,20,180,180,2000), byrow = T)
  n_real = 1000
  
  # Log-normal:
  # body.size.reals = matrix(nrow=length(taxo_groups),ncol=n_real,dimnames=list(taxo_groups,1:n_real),data=0)
  # for (taxon in taxo_groups[selected_groups])
  # {
  #   i_taxon = which(taxon == taxo_groups)
  #   body.size.reals[i_taxon,] = rlnorm(n_real,
  #                                      meanlog = 2*log(size_relativeAbund[i_taxon]) - 1/2*log(sd_size_relativeAbund[i_taxon]^2+size_relativeAbund[i_taxon]^2),
  #                                      # meanlog = log(size_relativeAbund[i_taxon]) - 1/2*log((sd_size_relativeAbund[i_taxon]/size_relativeAbund[i_taxon])^2+1),
  #                                      sdlog = sqrt(log((sd_size_relativeAbund[i_taxon]/size_relativeAbund[i_taxon])^2+1)))
  # }
  # hist(as.numeric(body.size.reals[which(taxo_groups == "Bacillariophyta"),]))

  # Gamma:
  # body.size.reals = matrix(nrow=length(taxo_groups),ncol=n_real,dimnames=list(taxo_groups,1:n_real),data=0)
  # for (taxon in taxo_groups[selected_groups])
  # {
  #   i_taxon = which(taxon == taxo_groups)
  #   body.size.reals[i_taxon,] = rgamma(n_real,
  #                                      shape = size_relativeAbund[i_taxon]^2/sd_size_relativeAbund[i_taxon]^2,
  #                                      rate = size_relativeAbund[i_taxon]/sd_size_relativeAbund[i_taxon]^2)
  # }
  
  # Four categories:
  # body.size.reals = matrix(nrow=length(taxo_groups),ncol=n_real,dimnames=list(taxo_groups,1:n_real),data=0)
  # for (taxon in taxo_groups[selected_groups])
  # {
  #   i_taxon = which(taxon == taxo_groups)
  #   body.size.reals[i_taxon,] = sample(x=mean_sizes,size=n_real,prob=p.four.fractions[i_taxon,],replace=T)
  # }
  
  # Uniformly distributed continuous body size values:
  body.size.reals = matrix(nrow=length(taxo_groups),ncol=n_real,dimnames=list(taxo_groups,1:n_real),data=0)
  for (taxon in taxo_groups[selected_groups])
  {
    i_taxon = which(taxon == taxo_groups)
    fraction.counts = sample.int(4,size=n_real,prob=p.four.fractions[i_taxon,],replace=T)
    high.bound = 0
    for (i in 1:4)
    {
      low.bound = high.bound + 1
      high.bound = as.numeric(table(fraction.counts)[i]) + high.bound
      body.size.reals[i_taxon,low.bound:high.bound] = runif(n = as.numeric(table(fraction.counts)[i]), 
                                                            min = fraction.sizes[i,1],
                                                            max = fraction.sizes[i,2])
    }
  }
  
  #######
  # Mixed model across reals:
  data.lme = data.frame(Groups = taxo_groups[selected_groups],PCoA2 = NormalizedVI_pcoa[[3]][,2],body.size.reals[selected_groups,])
  melted.data.lme = reshape2::melt(data.lme,id.vars=1:2,
                                   measure.vars=3:(n_real+2),
                                   variable.name="Real",
                                   value.name="Body.size")
  melted.data.lme[,"Body.size"] = log10(melted.data.lme[,"Body.size"])
  test = lme(PCoA2 ~ Body.size, data = melted.data.lme, random = ~1|Groups)
  
  # Simple t-test across reals:
  melted_PCoA2 = rep(NormalizedVI_pcoa[[3]][,2],n_real)
  melted_body.size.reals = reshape2::melt(body.size.reals[selected_groups,])$value
  cor.test(as.vector(melted_PCoA2),log10(melted_body.size.reals),method = "spearman")
  
  ########
  # Correlation for each real:
  
  cor.spear.pcoa2.bodysize = cor.pears.pcoa2.bodysize = vector(length = n_real, mode = "numeric")
  p.spear.pcoa2.bodysize = p.pears.pcoa2.bodysize = vector(length = n_real, mode = "numeric")
  for (i_real in 1:n_real)
  {
    test = cor.test(NormalizedVI_pcoa[[3]][,2], log10(body.size.reals[selected_groups,i_real]),method = "spearman")
    cor.spear.pcoa2.bodysize[i_real] = test$estimate
    p.spear.pcoa2.bodysize[i_real] = test$p.value
    test = cor.test(NormalizedVI_pcoa[[3]][,2], log10(body.size.reals[selected_groups,i_real]),method = "pearson")
    cor.pears.pcoa2.bodysize[i_real] = test$estimate
    p.pears.pcoa2.bodysize[i_real] = test$p.value
  }
  # hist(cor.spear.pcoa2.bodysize)
  # hist(cor.pears.pcoa2.bodysize)
  # hist(p.spear.pcoa2.bodysize,breaks=100)
  # hist(p.pears.pcoa2.bodysize,breaks=100)
  
  cor.spear.pcoa2.varpart = cor.pears.pcoa2.varpart = vector(length = n_real, mode = "numeric")
  p.spear.pcoa2.varpart = p.pears.pcoa2.varpart = vector(length = n_real, mode = "numeric")
  for (i_real in 1:n_real)
  {
    test = cor.test((colSums(varpart.env.spatial[[1]][2:3,])/(colSums(varpart.env.spatial[[1]][1:2,])))[selected_groups],
                    log10(body.size.reals[selected_groups,i_real]), method = "spearman")
    cor.spear.pcoa2.varpart[i_real] = test$estimate
    p.spear.pcoa2.varpart[i_real] = test$p.value
  }
  # hist(cor.spear.pcoa2.varpart)
  # hist(p.spear.pcoa2.varpart,breaks=100)
  
  ############ Figure:
  library(hexbin)
  library(viridis)
  
  ########################################################
  # pcoa2.vs.bodysize.reals.plot = cor.plot(x = reshape2::melt(body.size.reals[selected_groups,])$value,
  #                                         y = rep(NormalizedVI_pcoa[[3]][,2],n_real),
  #                                         x.lab = expression("Sampled body size ("*mu*"m), log scale"),
  #                                         y.lab = "PCoA axis 2",
  #                                         x.log = T,
  #                                         y.log = F,
  #                                         fit = T)
  pcoa2.vs.bodysize.reals.plot = ggplot(data = data.frame(x = reshape2::melt(body.size.reals[selected_groups,])$value,
                                                          y = rep(NormalizedVI_pcoa[[3]][,2],n_real))) + 
    # geom_bin2d(aes(x,y)) +
    stat_binhex(aes(x,y), show.legend=T) +
    scale_fill_viridis(name = "Count") +
    # scale_fill_gradient(name = "count", trans = "log") +
                        # breaks = my_breaks, labels = my_breaks) +
    labs(y = "Biogeographic axis 2", 
         x = expression("Sampled body size ("*mu*"m), log scale")) +
    theme_bw() +
    theme(axis.title=element_text(size=22),
          axis.text = element_text(size=22),
          plot.title=element_text(hjust=0, size=24),
          legend.title = element_text(size=22),
          legend.text = element_text(size=20),
          legend.margin = unit(c(1,1,1,1),"mm"),
          plot.margin=unit(c(5,1,1,5),"mm")) +
    scale_x_log10() +
    geom_smooth(aes(x,y), method='lm', col = "black", se=F)
  
  ########################################################
  # varpart.vs.bodysize.reals.plot = cor.plot(x = reshape2::melt(body.size.reals[selected_groups,])$value,
  #                                           y = rep((colSums(varpart.env.spatial[[1]][2:3,])/(colSums(varpart.env.spatial[[1]][1:2,])))[selected_groups],n_real),
  #                                           x.lab = expression("Sampled body size ("*mu*"m), log scale"),
  #                                           y.lab = "Surface connectivity-environment\n explained variance ratio",
  #                                           x.log = T,
  #                                           y.log = F,
  #                                           fit = T)
  varpart.vs.bodysize.reals.plot = ggplot(data = data.frame(x = reshape2::melt(body.size.reals[selected_groups,])$value,
                                                            y = rep((colSums(varpart.env.spatial[[1]][2:3,])/(colSums(varpart.env.spatial[[1]][1:2,])))[selected_groups],n_real))) + 
    # geom_bin2d(aes(x,y)) +
    stat_binhex(aes(x,y), show.legend=T) +
    scale_fill_viridis(name = "Count") +
    # scale_fill_gradient(name = "count", trans = "log") +
    labs(y = "Surface Moran maps-environment\n explained variance ratio", 
         x = expression("Sampled body size ("*mu*"m), log scale")) +
    theme_bw() +
    theme(axis.title=element_text(size=22),
          axis.text = element_text(size=22),
          plot.title=element_text(hjust=0, size=24),
          legend.title = element_text(size=22),
          legend.text = element_text(size=20),
          legend.margin = unit(c(1,1,1,1),"mm"),
          plot.margin=unit(c(5,1,1,3),"mm")) +
    scale_x_log10() +
    geom_smooth(aes(x,y), method='lm', col = "black", se=F)
  
  ########################################################
  # pcoa2.bodysize.hist.plot = ggplot(data = data.frame(x = c(cor.spear.distrib.pcoa2,cor.pears.distrib.pcoa2),
  #                                                     corr.coef = rep(c("spear","pears"),each = n_real))) +
  cor.spear.pcoa2.bodysize.hist.plot = ggplot(data = data.frame(x = cor.spear.pcoa2.bodysize)) +
    # geom_density(aes(x,colour=corr.coef)) +
    # geom_density(aes(x,fill=corr.coef), alpha = 0.2) +
    geom_histogram(aes(x),
                   # position = "stack",
                   # position = "identity", alpha = 0.5,
                   binwidth = 0.01, center = 0.005) +
    # scale_x_continuous(limits = c(0, 0.16)) +
    # scale_y_discrete(limits = c(0, 18)) +
    # scale_fill_manual(values = function_colors) +
    # ylim(c(0,16)) +
    # ylim(c(0,8)) +
    labs(y="Randomly sampled sets of group body sizes", x="Spearman correlation between biogeo. axis 2\n and log body size") +
    theme_bw() +
    # ggtitle(paste("MEM",k)) +
    # ggtitle(paste("MEM",k,"- Char.sc. =",
    #               signif(charac_dist.SUR_scale[k], digits = 3),"km")) +
    # signif(charac_tmin.SUR_scale[k], digits = 3),"y")) +
    theme(axis.title=element_text(size=22),
          axis.text = element_text(size=22),
          plot.title=element_text(hjust=0, size=24),
          # legend.position = "none",
          plot.margin=unit(c(5,5,5,5),"mm")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    xlim(min(cor.spear.pcoa2.bodysize,cor.spear.pcoa2.varpart),
         max(cor.spear.pcoa2.bodysize,cor.spear.pcoa2.varpart))
  
  ########################################################
  cor.spear.pcoa2.varpart.hist.plot = ggplot(data = data.frame(x = cor.spear.pcoa2.varpart)) +
    geom_histogram(aes(x),
                   binwidth = 0.01, center = 0.005) +
    labs(y="Randomly sampled sets of group body sizes", x="Spearman correlation btw. biogeo. axis 2 and\n surf. Moran maps-environment expl. var. ratio") +
    theme_bw() +
    theme(axis.title=element_text(size=22),
          axis.text = element_text(size=22),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(5,5,5,5),"mm")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    xlim(min(cor.spear.pcoa2.bodysize,cor.spear.pcoa2.varpart),
         max(cor.spear.pcoa2.bodysize,cor.spear.pcoa2.varpart))

  ########################################################
  p.spear.pcoa2.bodysize.hist.plot = ggplot(data = data.frame(x = p.spear.pcoa2.bodysize)) +
    geom_histogram(aes(x),
                   bins = 100) +
                   # binwidth = 0.01, center = 0.005) +
    labs(y="Randomly sampled sets of group body sizes", x="Spearman p-val. (log scale) btw. PCoA axis 2  \nand log body size") +
    theme_bw() +
    theme(axis.title=element_text(size=22),
          axis.text = element_text(size=22),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(5,5,5,5),"mm")) +
    scale_x_log10() +
    geom_vline(xintercept = 0.05, linetype = "dashed")
  
  ########################################################
  p.spear.pcoa2.varpart.hist.plot = ggplot(data = data.frame(x = p.spear.pcoa2.varpart)) +
    geom_histogram(aes(x),
                   bins = 100) +
                   # binwidth = 0.01, center = 0.005) +
    labs(y="Randomly sampled sets of group body sizes", x="Spearman p-val. (log scale) btw. PCoA axis 2\nand surf. Moran maps-environment expl. var. ratio   ") +
    theme_bw() +
    theme(axis.title=element_text(size=22),
          axis.text = element_text(size=22),
          plot.title=element_text(hjust=0, size=24),
          plot.margin=unit(c(5,5,5,5),"mm")) +
    scale_x_log10(limits = c(4*10^(-5), max(p.spear.pcoa2.varpart))) +
    geom_vline(xintercept = 0.05, linetype = "dashed")
  
  pdf(paste0(figure_folder,"/FigS_body.size.sensitivity.pdf"),
      width=7.5*2.2/1.2/2*2,height=12/3*2/1.2*2,onefile=F)
  g1 = ggplotGrob(pcoa2.vs.bodysize.reals.plot)
  g2 = ggplotGrob(varpart.vs.bodysize.reals.plot)
  g3 = ggplotGrob(cor.spear.pcoa2.bodysize.hist.plot)
  g4 = ggplotGrob(cor.spear.pcoa2.varpart.hist.plot + theme(axis.title.y=element_blank()))
  # g5 = ggplotGrob(p.spear.pcoa2.bodysize.hist.plot)
  # g6 = ggplotGrob(p.spear.pcoa2.varpart.hist.plot + theme(axis.title.y=element_blank()))
  grid.arrange(grobs = list(g1,g2,g3,g4),
               layout_matrix = rbind(c(1, 2),c(3, 4)))
  dev.off()
  ########
}

# Comparison between groups across reals
{
  ####
  plot.mean.norm.VI.within.group.vs.Moran.I = cor.plot(x = c(I_square.observed_w.mean[selected_groups,3],I_square.observed_w.mean_allTaxa[3]),
                                                       y = c(1 - mean.normalized.VI.within.group[selected_groups], 1 - mean.normalized.VI.allTaxa),
                                                       col.factor = factor(c(rep("group",length(taxo_groups[selected_groups])),"AllTaxa"),levels=c("group","AllTaxa")),
                                                       col.vect = c("black","red"),
                                                       leg.name = NULL,
                                                       x.lab = "Surface short-distance spatial autocorr.",
                                                       y.lab = "Mean within-group posterior sample similarity",
                                                       y.lab.hjust = 0.5,
                                                       x.log = F,
                                                       y.log = F,
                                                       fit = T,
                                                       fit.display="pearson.p",
                                                       mar.vect = c(10,5,5,5))
  # pdf(paste0(figure_folder,"/FigS_mean.within.group.Norm.VI.vs.Moran.I",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
  #     width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  # print(plot.mean.norm.VI.within.group.vs.Moran.I)
  # dev.off()
  
  ######
  plot.mean.sim.vs.Moran.I = cor.plot(x = c(I_square.observed_w.mean[selected_groups,3],I_square.observed_w.mean_allTaxa[3]),
                                      y = c(mean_sim[selected_groups],mean_sim.allTaxa),
                                      col.factor = factor(c(rep("group",length(taxo_groups[selected_groups])),"AllTaxa"),levels=c("group","AllTaxa")),
                                      col.vect = c("black","red"),
                                      leg.name = NULL,
                                      x.lab = "Surface short-distance spatial autocorr.",
                                      y.lab = "Mean within-group posterior sample similarity\n relative to null expectation",
                                      y.lab.hjust = 0.5,
                                      x.log = F,
                                      y.log = F,
                                      fit = T,
                                      fit.display="pearson.p",
                                      mar.vect = c(10,5,5,5)) +
    ylim(0,1)
  # pdf(paste0(figure_folder,"/FigS_mean.sim.vs.Moran.I",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
  #     width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  # print(plot.mean.sim.vs.Moran.I)
  # dev.off()
  
  ####
  plot.R.hat.1000.vs.mean.sim = cor.plot(x = c(mean_sim[selected_groups],mean_sim.allTaxa),
                                         y = c(R.hat.1000[selected_groups],R.hat.1000.allTaxa),
                                         col.factor = factor(c(rep("group",length(taxo_groups[selected_groups])),"AllTaxa"),levels=c("group","AllTaxa")),
                                         col.vect = c("black","red"),
                                         leg.name = NULL,
                                         x.lab = "Mean within-group posterior sample similarity      \n relative to null expectation",
                                         y.lab = bquote(hat(R)),
                                         y.lab.hjust = 0.5,
                                         x.log = F,
                                         y.log = T,
                                         fit = F,
                                         fit.display="pearson.spearman",
                                         mar.vect = c(10,5,5,5)) +
    geom_hline(yintercept = 1, linetype="dashed")
  # pdf(paste0(figure_folder,"/FigS_R.hat.1000.vs.mean.sim",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
  #     width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  # print(plot.R.hat.1000.vs.mean.sim)
  # dev.off()
  
  ####
  x = normalized.VI.best.real[selected_groups,selected_groups]
  x = 1 - x[lower.tri(x)]
  x.melted = rep(x,10)
  y.mat = matrix(nrow = length(x), ncol = 10, data = 0)
  for (i_real in 1:10)
  {
    y = normalized.VI.10.reals[[i_real]][selected_groups,selected_groups]
    y.mat[,i_real] = 1 - y[lower.tri(y)]
  }
  y.melted = reshape2::melt(y.mat)$value
  plot.normalized.VI.10real.vs.bestReal = cor.plot(x = x.melted,
                                                   y = y.melted,
                                                   x.lab = "Similarity between groups",
                                                   y.lab = "Similarity between groups\n for 10 posterior samples",
                                                   y.lab.hjust = 0.5,
                                                   x.log = F,
                                                   y.log = F,
                                                   fit = F,
                                                   fit.display="pearson.spearman",
                                                   mar.vect = c(5,5,5,5)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")
  # pdf(paste0(figure_folder,"/FigS_Normalized.VI.10reals.vs.best.real",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
  #     width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  # print(plot.normalized.VI.10real.vs.bestReal)
  # dev.off()
  
  pdf(paste0(figure_folder,"/FigS_between.groups.MCMC.convergence",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2*2,height=12/3*4/1.2/2*2,onefile=F)
  g1 = ggplotGrob(plot.mean.norm.VI.within.group.vs.Moran.I)
  g2 = ggplotGrob(plot.mean.sim.vs.Moran.I)
  g3 = ggplotGrob(plot.R.hat.1000.vs.mean.sim)
  g4 = ggplotGrob(plot.normalized.VI.10real.vs.bestReal)
  grid.arrange(grobs = list(g1, g2, g3, g4), nrow = 2, ncol = 2)
  dev.off()
  
  # R.hat histogram:
  ####
  plot.R.hat.1000 = ggplot(data = data.frame(x = R.hat.1000[selected_groups])) +
    geom_histogram(aes(x), size = 1) +
                   # position = "stack",
                   # position = "identity", alpha = 0.5,
                   # binwidth = 0.005, center = 0.0025) +
    # scale_x_continuous(limits = c(0, 0.16)) +
    # ylim(c(0,16)) +
    # ylim(c(0,8)) +
    labs(y="Number of plankton groups", x=bquote(hat(R))) +
    theme_bw() +
    # ggtitle(paste("MEM",k,"- Char.sc. =",
    #               signif(charac_dist.SUR_scale[k], digits = 3),"km")) +
    theme(axis.title=element_text(size=15),
          axis.text = element_text(size=15),
          plot.title=element_text(hjust=0, size=17),
          legend.position = "none",
          plot.margin=unit(c(1,1,1,0.5),"mm")) +
    # Add median adj. R2 annotation:
    # annotate(geom="text",
    #          x=median.adjR2.per.MEM[k] + 0.16*0.05,
    #          y=0.9*16,
    #          # label=bquote(atop("Ch.s.s."==.(format(charac_dist.SUR_scale[k],digits=0,nsmall=3))*.km,
    #          #                   "Median"==.(format(median.adjR2.per.MEM[k],digits=2)))),
    #          label = bquote("Median"==.(format(median.adjR2.per.MEM[k],digits=2))),
    #          hjust = 0,
    #          size=6) +
    # Adding letter labels:
    # annotate(geom="text",
    #          x = 0.16*0.05,
    #          y=0.9*16,
    #          label=LETTERS[kk],
    #          # label = bquote("Median"==.(format(median.adjR2.per.MEM[k],digits=2))),
    #          hjust = 1,
    #          size=9) +
    # # Adding charac. scale annotations:
    # annotate(geom="text",
    #          x=0.16*0.5,
    #          y=0.6*16,
    #          label=paste("Charac.scales:\n",
    #                      "+",signif(charac_tmin.SUR_scale[k], digits = 3),"y\n",
    #                      "+",signif(charac_dist.SUR_scale[k], digits = 3),"km"),
    #          # label = bquote("Median"==.(format(median.adjR2.per.MEM[k],digits=2))),
    #          hjust = 0,
    #          size=6) +
    geom_vline(xintercept = R.hat.1000.allTaxa, linetype = "dashed", size = 0.8)
  
  pdf(paste0(figure_folder,"/FigS_R.hat.1000.hist",noArcticNoBiomark_insert,noLagoon_insert,"_selected100+1.pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(plot.R.hat.1000)
  dev.off()
}

# Comparison across reals AllTaxa
{
  sampled_logLiks = readRDS(paste0(data_folder,"/18S_V9_TARA_CompleteSizeRange_byStationByDepth_AllTaxa_noLagoon/Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics16_nb_iter1000_nb_real100_meanPosteriorDistributedLlh_thin25_burnin2000_occurrence/posterior_sampled_logLiks.rds"))
  chains.allTaxa = readRDS(paste0(results_folder,"/Chains.allTaxa_nb_iter3000_nb_real100_occurrence.rds"))
 
  x.samples = y.samples = vector(length = 100, mode = "numeric")
  for (i_chain in 1:100)
  {
    selected.sample.indices = sort.int(abs(sampled_logLiks[i_chain,] - mean(sampled_logLiks[i_chain,])),index.return = T)$ix[1]
    x.samples[i_chain] = (selected.sample.indices-1)*25+2000
    y.samples[i_chain] = chains.allTaxa[x.samples[i_chain],i_chain]
  }
  
  data.folder_name = paste0(data_folder,"/18S_V9_TARA_CompleteSizeRange_byStationByDepth_AllTaxa_noLagoon")
  param.folder_name = paste0("Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics16_nb_iter1000_nb_real100_meanPosteriorDistributedLlh_thin25_burnin2000_occurrence")
  Ordered_realizations = readRDS(paste0(data.folder_name,"/",param.folder_name,"/Ordered_realizations.rds"))
  
  # random.reals = sample(2:100,2)
  random.reals = c(37,64)
  
  col.vect = muted(c("red","blue","green"),l=60,c=150)
  # show_col(col.vect)
  start.chain = 7
  col.factor = rep(col.vect,each=length(seq(from=start.chain,to=3000,by=5)))
  chain.curv.plot = curv.plot(x=seq(from=start.chain,to=3000,by=5),
                              y=chains.allTaxa[seq(from=start.chain,to=3000,by=5),-Ordered_realizations$ix[c(1,random.reals)]],
                              # col.vect = col.vect,
                              # leg.name = NULL,
                              x.lab="Iterations",
                              y.lab="Log-likelihood",
                              size.factor = 0.7,
                              line.size = 0.1) +
    geom_line(data = data.frame(x1 = rep(seq(from=start.chain,to=3000,by=5),3),
                                y1 = reshape2::melt(chains.allTaxa[seq(from=start.chain,to=3000,by=5),Ordered_realizations$ix[c(1,random.reals)]])$value,
                                group=col.factor, colour=col.factor),
              aes(x1,y1,group=factor(group,levels=col.vect),colour=factor(group,levels=col.vect)), size = 0.01,
              show.legend = F) +
    scale_colour_manual(values = col.vect) +
    geom_point(data = data.frame(x2 = x.samples[Ordered_realizations$ix[c(1,random.reals)]], 
                                 y2 = y.samples[Ordered_realizations$ix[c(1,random.reals)]], 
                                 colour = col.vect),
               aes(x2,y2,colour=factor(colour,levels=col.vect)),
               size = 3,
               show.legend = F) +
    geom_vline(xintercept = 2000, linetype = "dashed") +
    ylim(-18200000,
         max(reshape2::melt(chains.allTaxa[seq(from=start.chain,to=3000,by=5),])$value))
    # scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits=4))
  pdf(paste0(figure_folder,"/FigS_100.chains.AllTaxa",noArcticNoBiomark_insert,noLagoon_insert,".pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2,onefile=F)
  print(chain.curv.plot)
  dev.off()
}

# Comparing alpha-delta AllTaxa
{
  data.folder_name = paste0(data_folder_workspace3,"/18S_V9_TARA_CompleteSizeRange_byStationByDepth_AllTaxa_noLagoon")
  
  nb_topics_range100r = 3:20
  K.folder_name_cases = c(paste0(data.folder_name,"/Rtopicmodels_LDA_VEM_nb_topics",nb_topics_range100r,"_nb_real100_em_tol1e-06_var_tol1e-08_best_keep_occurrence/"),
                          paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha50:K_delta0.1_nb_topics",nb_topics_range100r,"_nb_iter2000_nb_real100_occurrence/"),
                          paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics",nb_topics_range100r,"_nb_iter2000_nb_real100_occurrence/"),
                          # paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics",nb_topics_range100r,"_nb_iter1000_nb_real100_meanPosteriorDistributedLlh_thin25_burnin2000_occurrence/"),
                          paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha0.05:K_delta0.1_nb_topics",nb_topics_range100r,"_nb_iter2000_nb_real100_occurrence/"),
                          paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha1_delta1_nb_topics",nb_topics_range100r,"_nb_iter2000_nb_real100_occurrence/"),
                          paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha0.1_delta1_nb_topics",nb_topics_range100r,"_nb_iter2000_nb_real100_occurrence/"),
                          paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha0.05:K_delta",9.4-nb_topics_range100r/5,"_nb_topics",nb_topics_range100r,"_nb_iter2000_nb_real100_occurrence/"))
                          # paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha0.05:K_delta0.5_nb_topics",nb_topics_range100r,"_nb_iter2000_nb_real100_occurrence/"))
  # paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha50:K_delta0.1_nb_topics",nb_topics_range100r,"_nb_iter250_nb_real100_best_thin2_burnin0_occurrence/"),
  # paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics",nb_topics_range100r,"_nb_iter250_nb_real100_best_thin2_burnin0_occurrence/"),
  # paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha0.05:K_delta0.1_nb_topics",nb_topics_range100r,"_nb_iter250_nb_real100_best_thin2_burnin0_occurrence/"))
  
  K.folder_name_cases.perplexity = c(paste0(data.folder_name,"/Rtopicmodels_LDA_VEM_alpha0.1_nb_topics2-18_post-predictive-cross-valid_fold_size10_em_tol1e-06_var_tol1e-08_occurrence/"),
                                     paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha50:K_delta0.1_nb_topics2-20_post-predictive-cross-valid_nb_iter1000_fold_size10_occurrence/"),
                                     paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics2-20_post-predictive-cross-valid_nb_iter2000_fold_size10_occurrence/"),
                                     paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha0.05:K_delta0.1_nb_topics2-20_post-predictive-cross-valid_nb_iter1000_fold_size10_occurrence/"),
                                     paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha1_delta1_nb_topics2-20_post-predictive-cross-valid_nb_iter1000_fold_size10_occurrence/"),
                                     paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha0.1_delta1_nb_topics2-20_post-predictive-cross-valid_nb_iter1000_fold_size10_occurrence/"),
                                     paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha0.05:K_delta9.4-K:5_nb_topics2-20_nb_iter2000_nb_real100_occurrence/"))
                                     # paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha0.05:K_delta0.5_nb_topics2-20_nb_iter2000_nb_real100_occurrence/"))
  # paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha50:K_delta0.1_nb_topics",nb_topics_range100r,"_nb_iter250_nb_real100_best_thin2_burnin0_occurrence/"),
  # paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics",nb_topics_range100r,"_nb_iter250_nb_real100_best_thin2_burnin0_occurrence/"),
  # paste0(data.folder_name,"/Rtopicmodels_LDA_Gibbs_alpha0.05:K_delta0.1_nb_topics",nb_topics_range100r,"_nb_iter250_nb_real100_best_thin2_burnin0_occurrence/"))
  
  inference_method = c("VEM",
                       "Gibbs\n chain end sample\ndelta = 0.1, alpha = 50/K",
                       "Gibbs\n chain end sample\ndelta = 0.1, alpha = 0.1",
                       "Gibbs\n chain end sample\ndelta = 0.1, alpha = 0.05/K",
                       "Gibbs\n chain end sample\ndelta = 1, alpha = 1",
                       "Gibbs\n chain end sample\ndelta = 1, alpha = 0.1",
                       "Gibbs\n chain end sample\ndelta = 9.4-K/5, alpha = 0.05/K")
  cases = 1:length(inference_method)
  
  # "Gibbs\n best sample\n alpha = 50/K","Gibbs\n best sample\n alpha = 0.1","Gibbs\n best sample\n alpha = 0.05/K")
  mean_sim_100r = matrix(nrow = length(nb_topics_range100r), ncol = length(inference_method), dimnames = list(nb_topics_range100r,inference_method), data=0)
  sim_intercept_100r = matrix(nrow = length(nb_topics_range100r), ncol = length(inference_method), dimnames = list(nb_topics_range100r,inference_method), data=0)
  max_llh100r = matrix(nrow = length(nb_topics_range100r), ncol = length(inference_method), dimnames = list(nb_topics_range100r,inference_method), data=0)
  llh100r = matrix(nrow = length(nb_topics_range100r), ncol = length(inference_method), dimnames = list(nb_topics_range100r,inference_method), data=0)
  AIC100r = matrix(nrow = length(nb_topics_range100r), ncol = length(inference_method), dimnames = list(nb_topics_range100r,inference_method), data=0)
  sd_llh100r = matrix(nrow = length(nb_topics_range100r), ncol = length(inference_method), dimnames = list(nb_topics_range100r,inference_method), data=0)
  est_delta = list()
  # load(paste0(data.folder_name,"/data2m.Rdata"))
  # data2m[data2m > 0] = 1
  # sum_data2m = sum(data2m)
  # nrow_data2m = nrow(data2m)
  # remove(data2m)
  llh_allpoints100r = list()
  AIC_allpoints100r = list()
  for (i_case in cases)
  {
    cat(i_case,"/",length(cases),"\n")
    llh_allpoints100r[[i_case]] = list()
    AIC_allpoints100r[[i_case]] = list()
    for (K in nb_topics_range100r)
    {
      K.folder_name = K.folder_name_cases[which(nb_topics_range100r == K) + (i_case-1)*length(nb_topics_range100r)]
      if (file.exists(K.folder_name))
      {
        if (i_case == 1)
        {
          if (file.exists(paste0(K.folder_name,"1st_best_realization/assemblage_composition.rds")))
          {
            assemblage_composition = readRDS(paste0(K.folder_name,"1st_best_realization/assemblage_composition.rds"))
            # assemblage_composition[,1:K][assemblage_composition[,1:K] < 1/sum_data2m] = 0
            # assemblage_composition = assemblage_composition[1:max(apply(assemblage_composition,2,function(g) length(which(g > 0))))]
            # est_delta[K-2,1:K] = Compositional::diri.est(as.matrix(assemblage_composition[,1:K]), type = "mle")
            est_delta[[K-2]] = sirt::dirichlet.mle(assemblage_composition[,1:K])$alpha
          } else 
            est_delta[[K-2]] = NA
        }
        
        stability_file = paste0(K.folder_name,"Stability_assessment_samplewise_maxmatching/Stability_samplewise.rds")
        if (file.exists(stability_file))
        {
          stability_data.frame = readRDS(stability_file)
          mean_sim_100r[K-2,i_case] = stability_data.frame$`Normalized ES`[1]
          sim_intercept_100r[K-2,i_case] = stability_data.frame$`Normalized ES=f(llh) intercept`[1]
        } else 
        {
          mean_sim_100r[K-2,i_case] = NA
          sim_intercept_100r[K-2,i_case] = NA
        }
        
        ordered_real_file = paste0(K.folder_name,"Ordered_realizations.rds")
        if (file.exists(ordered_real_file))
        {
          ordered_real = readRDS(ordered_real_file)
          max_llh100r[K-2,i_case] = ordered_real$x[1]
          llh100r[K-2,i_case] = mean(ordered_real$x)
          sd_llh100r[K-2,i_case] = sd(ordered_real$x)
          llh_allpoints100r[[i_case]][[K-2]] = ordered_real$x
          if (i_case == 1)
          {
            AIC_allpoints100r[[i_case]][[K-2]] = 2*(K*(nrow_data2m-1) + 1 - ordered_real$x)
          } else
            AIC_allpoints100r[[i_case]][[K-2]] = -2*ordered_real$x
          AIC100r[K-2,i_case] = mean(AIC_allpoints100r[[i_case]][[K-2]])
        } else
        {
          max_llh100r[K-2,i_case] = NA
          llh100r[K-2,i_case] = NA
          sd_llh100r[K-2,i_case] = NA
          llh_allpoints100r[[i_case]][[K-2]] = rep(NA,100)
        }
      } else
      {
        mean_sim_100r[K-2,i_case] = NA
        sim_intercept_100r[K-2,i_case] = NA
        max_llh100r[K-2,i_case] = NA
        llh100r[K-2,i_case] = NA
        sd_llh100r[K-2,i_case] = NA
        AIC100r[K-2,i_case] = NA
        llh_allpoints100r[[i_case]][[K-2]] = rep(NA,100)
        AIC_allpoints100r[[i_case]][[K-2]] = rep(NA,100)
      }
    }
  }
  
  perplexity = matrix(nrow = 20, ncol = length(cases), data=NA)
  perplexity_allpoints = list()
  for (i_case in cases)
  {
    perplexity_allpoints[[i_case]] = matrix(nrow = 17, ncol = 20, data=NA)
    K.folder_name = K.folder_name_cases.perplexity[i_case]
    if (file.exists(K.folder_name) && i_case != 1)
    {
      perplexity_file = paste0(K.folder_name,"perplexity.Rdata")
      if (file.exists(perplexity_file))
      {
        load(perplexity_file)
        i_K = 0
        for (K in nb_topics_range)
        {
          i_K = i_K+1
          perplexity_allpoints[[i_case]][,K] = perplexity_mpar[,i_K]
          perplexity[K,i_case] = median(perplexity_mpar[,i_K])
        }
      }
    }
  }
  
  #########
  plot.mean.sim.vs.K = ggplot(data = data.frame(K = rep(nb_topics_range100r,length(cases)), 
                                                y = reshape2::melt(mean_sim_100r)$value, 
                                                colour = rep(inference_method,each=length(nb_topics_range100r)))) +
    geom_point(aes(K,y,
                   colour=factor(colour,levels=inference_method[c(3,2,4,6,5,7,1)]))) + 
    geom_line(aes(K,y,
                  group=factor(colour,levels=inference_method[c(3,2,4,6,5,7,1)]),
                  colour=factor(colour,levels=inference_method[c(3,2,4,6,5,7,1)]))) +
    theme_bw() +
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          axis.title.y = element_text(hjust = 0.5),
          plot.title=element_text(hjust=0, size=24),
          legend.text = element_text(size=15),
          legend.title = element_text(size=20),
          # legend.position = c(0.7,0.35),
          legend.position = c(0.65,0.48),
          legend.box.background = element_rect(color="black", size=1.5),
          plot.margin=unit(c(1,1,1,0.5),"mm")) +
    labs(x="Number K of assemblages", y="Mean similarity between 100 realizations\n relative to null expectation") +
    scale_colour_discrete(name = "Inference",
                          labels = c("VEM",
                                     bquote("MCMC"~delta==0.1~","~alpha==50/K),
                                     bquote("MCMC"~delta==0.1~","~alpha==0.1),
                                     bquote("MCMC"~delta==0.1~","~alpha==0.05/K),
                                     bquote("MCMC"~delta==1~","~alpha==1),
                                     bquote("MCMC"~delta==1~","~alpha==0.1),#)[c(3,4,2,6,5,1)])
                                     bquote("MCMC"~delta==9.4-K/5~","~alpha==0.05/K))[c(3,2,4,6,5,7,1)])
    # guides(colour = guide_legend(override.aes = list(size=3)))
  
  ###########
  # plot.perplexity = ggplot(data = data.frame(K = rep(rep(1:20,each=17),length(cases)), 
  #                                            llh = unlist(perplexity_allpoints[cases]),
  #                                            colour = rep(inference_method,each=17*20))) +
    # geom_boxplot(aes(x=factor(K), 
    #                  y=llh, 
    #                  colour=factor(colour))) + 
  plot.perplexity = ggplot(data = data.frame(K = rep(1:20,length(cases)), 
                                             llh = reshape2::melt(perplexity[,cases])$value,
                                             colour = rep(inference_method,each=20))) +
      geom_line(aes(x = K, 
                    y = llh, 
                    colour=factor(colour,levels=inference_method[c(3,4,2,6,5,7,1)]), 
                    group = factor(colour,levels=inference_method[c(3,4,2,6,5,7,1)]))) +
    geom_point(aes(x = K, 
                  y = llh, 
                  colour=factor(colour,levels=inference_method[c(3,4,2,6,5,7,1)]))) +
    theme_bw() + 
    theme(axis.text = element_text(size=22),
          axis.title=element_text(size=22),
          axis.title.y = element_text(hjust = 0.5),
          plot.title=element_text(hjust=0, size=24),
          legend.position = "none",
          plot.margin=unit(c(1,1,1,0.5),"mm")) +
    labs(x="Number K of assemblages", y="Mean model perplexity on 10-sample folds")
  
  pdf(paste0(figure_folder,"/FigS_AllTaxa_stability-perplexity_vs_K_alpha.delta.comparison",noArcticNoBiomark_insert,noLagoon_insert,".pdf"),
      width=7.5*2.2/1.2/2,height=12/3*4/1.2/2*2,onefile=F)
  g1 = ggplotGrob(plot.mean.sim.vs.K)
  g2 = ggplotGrob(plot.perplexity)
  grid.arrange(grobs = list(g1, g2), nrow = 2, ncol = 1)
  dev.off()
}

