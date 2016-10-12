#!/usr/bin/env Rscript

# Exécution lda :
#./lda est 0.1 20 settings.txt /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/H20_GH_sequences_lda_norep.txt random /Users/guilhemsommeria-klein/Desktop/These/Données_H20/algo_result1
#Exécution hdp-faster :
#./hdp --train_data /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/H20_GH_sequences_lda_norep.txt --directory /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/hdp-faster_alpha1_gamma1_eta0.01_maxiter100_samplehyper_norep_1/ --max_time -1 --max_iter 100 --save_lag 10 --verbose --sample_hyper yes/no
# Exécution hdp :
#./hdp --algorithm train --data /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/H20_GH_sequences_lda_norep.txt --directory /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/hdp_alpha1_gamma1_eta0.01_maxiter100_inittopic100_norep_1/ --max_iter 100 --save_lag 10 --eta 0.01 --split_merge yes/no --sample_hyper yes/no --init_topics 100

# to know where packages are stored : '.libPaths' in the R console, then access variable .Library

# Connexion ssh EDB :
# ssh gsommeriaklein@10.2.3.11
# sshfs gsommeriaklein@10.2.3.11:/home/gsommeriaklein/ /Users/guilhemsommeria-klein/Desktop/Serveur_EDB/ 
# Connexion ssh Genotoul :
# ssh gsommeria@genotoul.toulouse.inra.fr
# sshfs gsommeria@genotoul.toulouse.inra.fr:/home/gsommeria/ /Users/guilhemsommeria-klein/Desktop/Serveur_Genotoul.home/ 
# sshfs gsommeria@genotoul.toulouse.inra.fr:/home/gsommeria/work/ /Users/guilhemsommeria-klein/Desktop/Serveur_Genotoul.work/
# Connexion bioclust03 ENS:
# sshfs sommeria@bioclusts03.bioclust.biologie.ens.fr:/users/biodiv/sommeria /Users/guilhemsommeria-klein/Desktop/Serveur_Bioclust

code_insert = "/Users/guilhemsommeria-klein/Desktop/These/Code/R/LDA_project/"

#######################
# Loading parameters: #
#######################

source(paste0(code_insert,"LDA_parameters_file.R"))  

if (locally_based)
{
  lib_path = "/Library/Frameworks/R.framework/Resources/library"
} else if (genotoul_cluster_based)
{
  lib_path = "/home/gsommeria/packages_R/"
} else if (EDB_cluster_based)
{
  lib_path = "/home/gsommeriaklein/packages_R/"
  #lib_path = "/usr/lib64/R/library/"    
}

# To use command-line arguments:
#library(optparse)
# To use color2D.matplot function :
library(plotrix)
# To use errbar function :
library(Hmisc)
# To use dudi.pca function :
library(ade4)
# To use rainbow_hcl function :
library(colorspace)
#library(lattice)
# To compute KL divergence using KL.plugin:
library(entropy,lib.loc=lib_path)
# To perform kriging:
library(kriging,lib.loc=lib_path)
# To plot the kriged maps:
library(ggplot2)
library(gridExtra)
library(scales)
if (genotoul_cluster_based)
{
  library(gstat,lib.loc=lib_path)
} else
{
  library(gstat)
}
# To call grid.newpage()
library(grid)
# To plot shapes using readOGR (not available on cluster):
if (locally_based)
{
  library(rgdal)
  library(sp)
  library(rgeos)
  library(maptools)
}
# To use raster manipulating functions:
if (genotoul_cluster_based)
{
  library(raster,lib.loc=lib_path)
} else
{
  library(raster)
}
# To use coltorgb:
library(grDevices)
# To load results:
if (genotoul_cluster_based)
{
  library(topicmodels,lib.loc=lib_path)
} else
{
  library(topicmodels)
}
# To use the agnes function (hierarchical clustering):
library(cluster)
# To use the cophenetic function (cophenetic distance on a clustering tree):
library(stats)

############################
# Parsing command-line options
# option_list = list(
#                   make_option(c("-f", "--file"), type="character", default=NULL, 
#                               help="dataset file name", metavar="character"),
#                   make_option(c("-o", "--out"), type="character", default="out.txt", 
#                               help="output file name [default= %default]", metavar="character")
#                   ) 
# 
# opt_parser = OptionParser(option_list=option_list)
# opt = parse_args(opt_parser)
# 
# if (is.null(opt$file))
# {
#   print_help(opt_parser)
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# }
############################

if (length(Selected_real)==1 && realization_comparison)
  stop("Needs more than one realization to compare realizations")

if (length(which(c(barcode_gh,barcode_ghassigned,barcode_itsfungi,barcode_itsasco,barcode_itsbasidio,barcode_itsglomero,bacteriafungiarchaea,barcode_16sbact,
                   barcode_18s,barcode_18sfungi,barcode_18smetazoa,barcode_18sannelids,barcode_18sarthropods,barcode_18snematodes,barcode_18sprotists,
                   barcode_18splatyhelminthes,barcode_18splants,barcode_16sarch,barcode_itsplant,barcode_16sins,plantfungi,frogs_iucn_all,frogs_iucn_forest,
                   frogs_iucn_open,frogs_poly1_all,frogs_poly1_forest,frogs_poly1_open,frogs_poly2_all,frogs_poly2_forest,frogs_poly2_open)==1))>1)
  stop("More than one barcode provided")

##################################
# DATA LOADING AND PREPROCESSING #
##################################

if (realization_comparison)
  source(paste0(code_insert,"LDA_realization_comparison_fun.R"))

# color2D.matplot colors
col_values = as.list(as.data.frame(t(col2rgb(terrain.colors(255)))))

if (data_h20)
{
  data_insert = "Donnees_H20"
  short_data_insert = "H20"
} else if (data_pp)
{
  data_insert = "Donnees_PetitPlateau"
  short_data_insert = "PP"
} else if (data_betadiv || data_betadiv_pooled)
{
  data_insert = "Données_betadiv"
  short_data_insert = "betadiv"
} else if (data_gs)
{
  data_insert = "Donnees_PlateauDesGuyanes"
  short_data_insert = "GS"
}

if (barcode_gh) {
  barcode_insert = "Plantes_GH"
  short_barcode_insert = "GH"
} else if (barcode_ghassigned)
{
  barcode_insert = "Plantes_GH_assigned"
} else if (barcode_itsfungi)
{
  barcode_insert = "Champignons_ITS"
  short_barcode_insert = "ITSfungi"
} else if (barcode_itsasco)
{
  barcode_insert = "Ascomycetes_ITS"
} else if (barcode_itsbasidio)
{
  barcode_insert = "Basidiomycetes_ITS"
} else if (barcode_itsglomero)
{
  barcode_insert = "Glomeromycetes_ITS"
} else if (barcode_16sbact)
{
  barcode_insert = "Bacteries_16S"
  short_barcode_insert = "16Sbact"
} else if (barcode_18s) 
{
  barcode_insert = "Eukaryotes_18S"
  short_barcode_insert = "18Seuka"
} else if (barcode_18sfungi) 
{
  barcode_insert = "Champignons_18S"
} else if (barcode_18smetazoa) 
{
  barcode_insert = "Metazoaires_18S"
} else if (barcode_18sannelids)
{
  barcode_insert = "Annelides_18S"
} else if (barcode_18sarthropods)
{
  barcode_insert = "Arthropodes_18S"
} else if (barcode_18snematodes)
{
  barcode_insert = "Nematodes_18S"
} else if (barcode_18sprotists)
{
  barcode_insert = "Protistes_18S"
} else if (barcode_18splatyhelminthes)
{
  barcode_insert = "Platyhelminthes_18S"
} else if (barcode_18splants) 
{
  barcode_insert = "Plantes_18S"
} else if (barcode_16sarch)
{
  barcode_insert = "Archees_16S"
  short_barcode_insert = "16Sarch"
} else if (barcode_itsplant)
{
  barcode_insert = "Plantes_ITS"
  short_barcode_insert = "ITSplant"
} else if (barcode_16sins)
{
  barcode_insert = "Insectes_16S"
  short_barcode_insert = "16Sins"
} else if (plantfungi)
{
  barcode_insert = "Plantes_GH-Champignons_ITS"
} else if (bacteriafungiarchaea)
{
  barcode_insert = "Bactéries-Champignons-Archées"
} else if (testdata)
{
  barcode_insert = paste("Test_data/",testdata_dir,sep="")
} else if (frogs_poly2_all)
{
  barcode_insert = "Anoures_poly2_all"
  short_barcode_insert = "poly2_all"
} else if (frogs_poly2_forest)
{
  barcode_insert = "Anoures_poly2_forest"
  short_barcode_insert = "poly2_forest"
} else if (frogs_poly2_open)
{
  barcode_insert = "Anoures_poly2_open"
  short_barcode_insert = "poly2_open"
} else if (frogs_poly1_all)
{
  barcode_insert = "Anoures_poly1_all"
  short_barcode_insert = "poly1_all"
} else if (frogs_poly1_forest)
{
  barcode_insert = "Anoures_poly1_forest"
  short_barcode_insert = "poly1_forest"
} else if (frogs_poly1_open)
{
  barcode_insert = "Anoures_poly1_open"
  short_barcode_insert = "poly1_open"
} else if (frogs_iucn_all)
{
  barcode_insert = "Anoures_IUCN_all"
  short_barcode_insert = "IUCN_all"
} else if (frogs_iucn_open)
{
  barcode_insert = "Anoures_IUCN_open"
  short_barcode_insert = "IUCN_open"
} else if (frogs_iucn_forest)
{
  barcode_insert = "Anoures_IUCN_forest"
  short_barcode_insert = "IUCN_forest"
}

if (blei) {
  filename_insert = "Blei_LDA"
} else if (Rtopicmodels_Gibbs) {
  filename_insert = "Rtopicmodels_LDA_Gibbs"
} else if (Rtopicmodels_VEM) {
  filename_insert = "Rtopicmodels_LDA_VEM"
} else if (hdp_faster)
  filename_insert = "HDP-faster"

#cluster_prefix = "/Users/guilhemsommeria-klein/Desktop/Serveur_EDB/"
cluster_prefix = "/Users/guilhemsommeria-klein/Desktop/Serveur_Genotoul.work/"
if (genotoul_cluster_based)
{
  local_prefix = "/home/gsommeria/work/"
} else if (EDB_cluster_based)
{
  local_prefix = "/home/gsommeriaklein/"
} else if (locally_based)
{
  local_prefix = "/Users/guilhemsommeria-klein/Desktop/These/" 
  cluster_prefix = "/Users/guilhemsommeria-klein/Desktop/Serveur_Genotoul.work/"
}

if (occurrence)
{
  occurrence_insert = "_occurrence"
} else 
{
  occurrence_insert = ""
}

if (nb_occupied_sites_threshold > 1)
{
  remove_single_sites_insert = paste0("_woOTUsWithLessThan",nb_occupied_sites_threshold,"Sites")
} else 
{
  remove_single_sites_insert = ""
}

if (no_rare)
{
  no_rare_insert = "_noRareOTU"
} else 
{
  no_rare_insert = ""
}

if (samplewise) {
  MOTU_sample_insert = "samplewise"
} else if (MOTUwise)
  MOTU_sample_insert = "MOTUwise"

if (bij) {
  bij_insert = "bijective_correspondence"
} else if (!bij)
  bij_insert = "non-bijective_correspondence"

# H20_GH.R (retrieve general information on the original data) :
###############
setwd(paste(local_prefix,data_insert,"/",barcode_insert,"/",sep=""))
# taxo_ref <- read.table(paste(short_data_insert,"_",short_barcode_insert,"_taxo_ref_table.txt",sep=""), sep=" ")
if (plantfungi)
{
  load("taxo_ref_plantfungi.Rdata")
} else if (bacteriafungiarchaea)
{
  load("taxo_ref_bacteriafungiarchaea.Rdata")
} else if (barcode_ghassigned)
{
  load("taxo_ref_assigned.Rdata")
  taxo_ref = taxo_ref_assigned
} else if (barcode_itsasco)
{
  load("taxo_ref_asco.Rdata")
  taxo_ref = taxo_ref_asco
} else if (barcode_itsbasidio)
{
  load("taxo_ref_basidio.Rdata")
  taxo_ref = taxo_ref_basidio
} else if (barcode_itsglomero)
{
  load("taxo_ref_glomero.Rdata")
  taxo_ref = taxo_ref_glomero
} else if (barcode_18sfungi)
{
  load("taxo_ref_euka_fungi.Rdata")
  taxo_ref = taxo_ref_euka_fungi
} else if (barcode_18splants)
{
  load("taxo_ref_euka_plants.Rdata")
  taxo_ref = taxo_ref_euka_plants
} else if (barcode_18smetazoa)
{
  load("taxo_ref_euka_metazoa.Rdata")
  taxo_ref = taxo_ref_euka_metazoa
} else if (barcode_18sannelids)
{
  load("taxo_ref_euka_annelids.Rdata")
  taxo_ref = taxo_ref_euka_annelids
} else if (barcode_18sarthropods)
{
  load("taxo_ref_euka_arthropods.Rdata")
  taxo_ref = taxo_ref_euka_arthropods
} else if (barcode_18snematodes)
{
  load("taxo_ref_euka_nematodes.Rdata")
  taxo_ref = taxo_ref_euka_nematodes
} else if (barcode_18sprotists)
{
  load("taxo_ref_euka_protists.Rdata")
  taxo_ref = taxo_ref_euka_protists
} else if (barcode_18splatyhelminthes)
{
  load("taxo_ref_euka_platyhelminthes.Rdata")
  taxo_ref = taxo_ref_euka_platyhelminthes
} else if (testdata)
{
  load("taxo_ref_testdata.Rdata")
} else if (!data_gs)
  taxo_ref = read.table(paste(short_data_insert,"_",short_barcode_insert,"_taxo_ref_table_assignScore.txt",sep=""), colClasses="vector", sep=" ")

# sample counts without replicates and controles - larger number of reads per sample than in _sequences_counts_norep.txt
# sample_counts <- as.vector(t(read.table(paste(short_data_insert,"_",short_barcode_insert,"_sample_counts.txt",sep=""))))

# Originally only computed in "H20_GH.R", now computed here too
# normal_ordered_seq <- as.vector(unlist(read.table("Seq_ordered_according_to_sitenormalized_abund.txt")))
# normal_data2m <- as.matrix(read.table("Site_composition_in_sequences.txt"))
# KL_normal_data2m1 <- as.matrix(read.table("Sequence_composition_in_sites.txt"))

# if (!existingresult || dominance_calculation || data_betadiv || data_betadiv_pooled || topic_read_proportion_comput)
# data2m should be loaded even in the case of exisitng reasults so as to be able to modify taxo_ref in case of empty rows
if (plantfungi)
{
  if (filled)
    load("data2m_filled.Rdata")
  else 
    load("data2m_plantfungi.Rdata")
} else if (bacteriafungiarchaea)
{
  if (filled)
    load("data2m_filled.Rdata")
  else 
    load("data2m_bacteriafungiarchaea.Rdata")
} else if (barcode_ghassigned)
{
  if (filled)
    load("data2m_filled.Rdata")
  else 
  {
    load("data2m_assigned.Rdata")
    data2m = data2m_assigned
  }
} else if (barcode_itsasco)
{
  if (filled)
    load("data2m_filled.Rdata")
  else
  {
    load("data2m_asco.Rdata")
    data2m = data2m_asco
  }
} else if (barcode_itsbasidio)
{
  if (filled)
    load("data2m_filled.Rdata")
  else
  {
    load("data2m_basidio.Rdata")
    data2m = data2m_basidio
  }
} else if (barcode_itsglomero)
{
  if (filled)
    load("data2m_filled.Rdata")
  else
  {
    load("data2m_glomero.Rdata")
    data2m = data2m_glomero
  }
} else if (barcode_18sfungi)
{
  if (filled)
    load("data2m_filled.Rdata")
  else
  {
    load("data2m_euka_fungi.Rdata")
    data2m = data2m_euka_fungi
  }
} else if (barcode_18splants)
{
  if (filled)
    load("data2m_filled.Rdata")
  else 
  {
    load("data2m_euka_plants.Rdata")
    data2m = data2m_euka_plants
  }
} else if (barcode_18smetazoa)
{
  if (filled)
    load("data2m_filled.Rdata")
  else 
  {
    load("data2m_euka_metazoa.Rdata")
    data2m = data2m_euka_metazoa
  }
} else if (barcode_18sannelids)
{
  if (filled)
    load("data2m_filled.Rdata")
  else 
  {
    load("data2m_euka_annelids.Rdata")
    data2m = data2m_euka_annelids
  }
} else if (barcode_18sarthropods)
{
  if (filled)
    load("data2m_filled.Rdata")
  else 
  {
    load("data2m_euka_arthropods.Rdata")
    data2m = data2m_euka_arthropods
  }
} else if (barcode_18snematodes)
{
  if (filled)
    load("data2m_filled.Rdata")
  else 
  {
    load("data2m_euka_nematodes.Rdata")
    data2m = data2m_euka_nematodes
  }
} else if (barcode_18sprotists)
{
  if (filled)
    load("data2m_filled.Rdata")
  else 
  {
    load("data2m_euka_protists.Rdata")
    data2m = data2m_euka_protists
  }
} else if (barcode_18splatyhelminthes)
{
  if (filled)
    load("data2m_filled.Rdata")
  else 
  {
    load("data2m_euka_platyhelminthes.Rdata")
    data2m = data2m_euka_platyhelminthes
  }
} else if (testdata)
{
  load("data2m_testdata.Rdata")
  load("true_documents.Rdata")
  true_documents = true_norm_documents
  # true_topic_compo to be loaded
} else 
{
  if (filled)
  {
    if (data_pp || data_h20)
      load("data2m_filled.Rdata")
    else if (data_gs)
    {
      data2 = read.table(paste(short_data_insert,"_",short_barcode_insert,".txt",sep=""), colClasses="vector", sep=" ")
      data2 = t(data2)
      data2m = as.matrix(data2[-c(1,2),-1])
      colnames(data2m) = vector(length=ncol(data2m),mode="character")
      for (k in 1:ncol(data2m))
        colnames(data2m)[k] = paste("x =",data2[1,k+1],"- y =",data2[2,k+1])
      coordGS = data.frame(x=as.numeric(data2[1,-1])*10,y=as.numeric(data2[2,-1])*10)
      taxo_ref = data.frame(V1=c("Sequence_number",seq(1,nrow(data2m),1)), V2=c("kingdom_name_ok",rep(NA,nrow(data2m))), V3=c("phylum_name_ok",rep(NA,nrow(data2m))), V4=c("class_name_ok",rep(NA,nrow(data2m))),
                            V5=c("order_name_ok",rep(NA,nrow(data2m))), V6=c("family_name_ok",rep(NA,nrow(data2m))), V7=c("genus_name_ok",rep(NA,nrow(data2m))), V8=c("species_name_ok",rep(NA,nrow(data2m))),
                            V9=c("scientific_name_ok",data2[-c(1,2),1]), V10=c("taxonomic_rank_ok",rep("species",nrow(data2m))), V11=c("best_identity_ok",rep(1,nrow(data2m))))
    }
  } else 
  {
    data2 = read.table(paste(short_data_insert,"_",short_barcode_insert,"_sequences_counts_norep.txt",sep=""), colClasses="vector", sep=" ")
    if (data_betadiv)
    {
      if (length(grep("^sample\\..*T$",data2[1,])) > 0)
        data2 = data2[,-grep("^sample\\..*T$",data2[1,])]
      # removing H20-Mobio samples
      if (length(grep("^sample\\.NH20.M$",data2[1,])) > 0)
        data2 = data2[,-grep("^sample\\.NH20.M$",data2[1,])]
    }
    data2m = as.matrix(data2[-1,])
    colnames(data2m) = data2[1,]
  }
}
data2m = apply(data2m,2,as.numeric)
rownames(data2m) = seq(1,nrow(data2m),1)

if (data_pp && !filled)
{
  #Reminding the spatial positions (column indices) where a blank space needs to be introduced when plotting the spatial distibution of samples
  # cat("Detecting missing samples ...\n")
  # sample_names = colnames(data2m)
  # NEXT_CHAR = function(previous_last_char_index) {
  #   if (previous_last_char_index == 29) {
  #     next_last_char_index = 1
  #   } else {
  #     next_last_char_index = previous_last_char_index+1}
  #   next_last_char_index
  # }
  # Sample_name_endings = c("7","8","9","A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")
  # # Initialization
  # Missing_positions_indices = vector(length=29*39,mode="numeric")
  # char_chain = as.character(sample_names[1])
  # last_char = substr(char_chain,nchar(char_chain),nchar(char_chain))
  # last_char_index = which(Sample_name_endings==last_char)
  # spatial_index = 1
  # if (last_char_index!=1)
  #   Missing_positions_indices[spatial_index] = 1 
  # previous_last_char_index = 1 
  # # loop over non-empty samples
  # for (j in 2:length(sample_names)) 
  # {
  #   spatial_index = spatial_index+1
  #   char_chain = as.character(sample_names[j])
  #   last_char = substr(char_chain,nchar(char_chain),nchar(char_chain))
  #   last_char_index = which(Sample_name_endings==last_char)
  #   if (last_char_index!=NEXT_CHAR(previous_last_char_index))
  #   {
  #     Missing_positions_indices[spatial_index] = 1
  #     spatial_index = spatial_index+1
  #     iter_last_char_index = NEXT_CHAR(previous_last_char_index)
  #     while (last_char_index != NEXT_CHAR(iter_last_char_index))
  #     {
  #       Missing_positions_indices[spatial_index] = 1
  #       iter_last_char_index = NEXT_CHAR(iter_last_char_index)
  #       spatial_index = spatial_index+1
  #     }
  #   }
  #   previous_last_char_index = last_char_index
  # }
  
  # Loading the precomputed Missing_positions_indices (one file for each barcode)
  load("Missing_positions_indices.Rdata")
} else if (data_betadiv)
  load("Missing_positions_indices_nomobio.Rdata")

OTUs_to_be_removed = vector(length=nrow(data2m),mode="logical")
Prop_OTU_removed = 0
Prop_reads_removed = 0
for (OTU in 1:nrow(data2m))
  # Removing all OTUs which are present in less than nb_occupied_sites_threshold
  OTUs_to_be_removed[OTU] = (length(which(data2m[OTU,]>0))<nb_occupied_sites_threshold)
if (no_rare)
{
  for (OTU in 1:nrow(data2m))
  {
    if (OTUs_to_be_removed[OTU] == 0)
      OTUs_to_be_removed[OTU] = !(sum(data2m[OTU,])>length(data2m[OTU,]))
  }
}
if (any(OTUs_to_be_removed))
{
  Prop_OTU_removed = length(which(OTUs_to_be_removed))/nrow(data2m)
  Prop_reads_removed = sum(data2m[which(OTUs_to_be_removed),])/sum(data2m)
  data2m = data2m[-which(OTUs_to_be_removed),]
  taxo_ref = taxo_ref[-(which(OTUs_to_be_removed)+1),]
  rownames(data2m) = seq(1,nrow(data2m),1)
  rownames(taxo_ref) = seq(0,nrow(data2m),1)
  taxo_ref[2:(nrow(data2m)+1),1] = seq(1,nrow(data2m),1)
  count2 = rowSums(data2m)
}

if (occurrence)
  data2m[data2m>0] = 1

# Dealing with the case of split datasets where some samples are empty despite having filling them at the level of the whole dataset
filled_with_gaps = 0
if (filled && (length(which(colSums(data2m)==0))!=0))
{
  filled_with_gaps = 1
  Missing_positions_indices = vector(length=29*39,mode="numeric")
  Missing_positions_indices[which(colSums(data2m)==0)] = 1
  if (data_gs)
    coordGS = coordGS[-which(colSums(data2m)==0),]
  data2m = data2m[,-which(colSums(data2m)==0)]
}

if (data_betadiv_pooled)
{
  data2m_pooled = matrix(nrow=nrow(data2m),ncol=19)
  data2m_pooled[,1] = rowSums(data2m[,grep("^sample\\.ARB1",colnames(data2m))])
  data2m_pooled[,2] = rowSums(data2m[,grep("^sample\\.ARB2",colnames(data2m))])
  data2m_pooled[,3] = rowSums(data2m[,grep("^sample\\.NBA1",colnames(data2m))])
  data2m_pooled[,4] = rowSums(data2m[,grep("^sample\\.NBA2",colnames(data2m))])
  data2m_pooled[,5] = rowSums(data2m[,grep("^sample\\.NBA3",colnames(data2m))])
  data2m_pooled[,6] = rowSums(data2m[,grep("^sample\\.NF21",colnames(data2m))])
  data2m_pooled[,7] = rowSums(data2m[,grep("^sample\\.NH20",colnames(data2m))])
  data2m_pooled[,8] = rowSums(data2m[,grep("^sample\\.NH21",colnames(data2m))])
  data2m_pooled[,9] = rowSums(data2m[,grep("^sample\\.NL11",colnames(data2m))])
  data2m_pooled[,10] = rowSums(data2m[,grep("^sample\\.NL12",colnames(data2m))])
  data2m_pooled[,11] = rowSums(data2m[,grep("^sample\\.NO13",colnames(data2m))])
  data2m_pooled[,12] = rowSums(data2m[,grep("^sample\\.NPA5",colnames(data2m))])
  data2m_pooled[,13] = rowSums(data2m[,grep("^sample\\.NPA6",colnames(data2m))])
  data2m_pooled[,14] = rowSums(data2m[,grep("^sample\\.NPAC",colnames(data2m))])
  data2m_pooled[,15] = rowSums(data2m[,grep("^sample\\.PA01",colnames(data2m))])
  data2m_pooled[,16] = rowSums(data2m[,grep("^sample\\.PA02",colnames(data2m))])
  data2m_pooled[,17] = rowSums(data2m[,grep("^sample\\.PA03",colnames(data2m))])
  data2m_pooled[,18] = rowSums(data2m[,grep("^sample\\.PA11",colnames(data2m))])
  data2m_pooled[,19] = rowSums(data2m[,grep("^sample\\.PA12",colnames(data2m))])
  data2m = data2m_pooled
  # colnames(norm_data2m_pooled) = c("ARB1","ARB2","NBA1","NBA2","NBA3","NF21","NH20","NH21","NL11","NL12","NO13","NPA5","NPA6","NPAC","PA01","PA02","PA03","PA11","PA12")
  # Giving names corresponding to real plots, not sample names
  colnames(data2m) = c("AR74","AR73","NBA1","NBA2","NINS","NF21","NH20","NH21","NL11","NL12","NO13","NPA5","NPA6","NL18","P063","P121","P064","P111","P122")
  #data2m = apply(data2m,MARGIN=2,FUN=as.character)
}

# sample_count2 contains the number of reads per sample, in order to normalize the number of reads per sequence
sample_count2 = colSums(data2m)

# # computing the variance of the number of reads per sample
# var_sample_count2 = 0
# for (i in 1:ncol(data2m))
#   var_sample_count2 = (sample_count2[i] - sum(sample_count2)/ncol(data2m))^2/ncol(data2m) + var_sample_count2
# std_dev_sample_count2 = sqrt(var_sample_count2)/(sum(sample_count2)/ncol(data2m))

# normal_count2 contains the proportion of each sequence average over the sites
normal_count2 = rowSums(data2m)/ncol(data2m)
# nrow(data2m) = nb_seq et ncol(data2m) = nb_sites
# computing the read proportion per site for each sequence
normal_data2m = sweep(data2m,MARGIN=2,sample_count2,'/') 

sorted_normal_abundances2 = sort.int(normal_count2,decreasing=T,index.return=T)
normal_ordered_seq = sorted_normal_abundances2$ix
####################
# end condition (!existingresult || dominance_calculation || data_betadiv || data_betadiv_pooled)
#}

######################################
# END DATA LOADING AND PREPROCESSING #
######################################


#####################
# LDA DECOMPOSITION #
#####################

LLH_final0 = matrix(nrow=2,ncol=mpar_range,data=0)
LLH_final1 = matrix(nrow=nb_real,ncol=mpar_range,data=0)
AIC0 = matrix(nrow=6,ncol=mpar_range,data=0)
rownames(AIC0) = c("AIC_phi_theta","stdev_AIC_phi_theta","AIC_zi","stdev_AIC_zi","AIC_phi_alpha","stdev_AIC_phi_alpha")
# AIC phi theta:
AIC1 = matrix(nrow=nb_real,ncol=mpar_range,data=0)
# AIC zi:
AIC2 = matrix(nrow=nb_real,ncol=mpar_range,data=0)
# AIC phi alpha:
AIC3 = matrix(nrow=nb_real,ncol=mpar_range,data=0)
if (Rtopicmodels_VEM) 
{
  alpha_est1 = matrix(nrow=nb_real,ncol=mpar_range,data=0)
  alpha_est0 = matrix(nrow=2,ncol=mpar_range,data=0)
}

# creating the directories, and downloading the .Rdata file from the cluster for (mpar)
#####################
if (data_h20) 
{
  local_dirname = paste(local_prefix,data_insert,"/plots/",barcode_insert,"/",filename_insert,"/",sep="")
  cluster_dirname = paste(cluster_prefix,data_insert,"/plots/",barcode_insert,"/",filename_insert,"/",sep="")
} else if (data_pp || data_gs || data_betadiv || data_betadiv_pooled)
{
  if (filled)
  {
    local_dirname = paste(local_prefix,data_insert,"/",barcode_insert,"/",filename_insert,"/",sep="")
    cluster_dirname = paste(cluster_prefix,data_insert,"/",barcode_insert,"/",filename_insert,"/",sep="")
  } else
  {
    local_dirname = paste(local_prefix,data_insert,"/",barcode_insert,"/",filename_insert,"/Unfilled/",sep="")
    cluster_dirname = paste(cluster_prefix,data_insert,"/",barcode_insert,"/",filename_insert,"/Unfilled/",sep="")
  }
}

if (!(file.exists(local_dirname)))
  dir.create(local_dirname)
setwd(local_dirname)  

if (mpar)
{
  if (cluster)
    cat("Loading data from cluster ...\n")
  
  if (mnb_topics)
  {
    if (Rtopicmodels_Gibbs)
    {
      local_subdirname = paste0(local_dirname,filename_insert,"_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_iter",nb_iter,"_nb_real",nb_real,occurrence_insert,remove_single_sites_insert,no_rare_insert,"/")
      cluster_subdirname = paste0(cluster_dirname,filename_insert,"_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_iter",nb_iter,"_nb_real",nb_real,occurrence_insert,remove_single_sites_insert,no_rare_insert,"/")
      filename = paste0(filename_insert,"_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,"_nb_iter",nb_iter,occurrence_insert,remove_single_sites_insert,no_rare_insert,".Rdata")  
    } else if (Rtopicmodels_VEM)
    {
      local_subdirname = paste0(local_dirname,filename_insert,"_",alpha_insert,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,occurrence_insert,remove_single_sites_insert,no_rare_insert,"/")
      cluster_subdirname = paste0(cluster_dirname,filename_insert,"_",alpha_insert,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,occurrence_insert,remove_single_sites_insert,no_rare_insert,"/")
      filename = paste0(filename_insert,"_",alpha_insert,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,occurrence_insert,remove_single_sites_insert,no_rare_insert,".Rdata")
    }
    
    #     Not erasing the folder if it already exists:
    #     if (!(file.exists(local_subdirname)))
    #     {
    #       dir.create(local_subdirname)
    #       local_subdirname = paste(local_subdirname,"/",sep="")
    #     }
    #     else if (file.exists(local_subdirname))
    #     {
    #       local_subdirname_modified = local_subdirname
    #       j_index = 1
    #       while (file.exists(local_subdirname_modified))
    #       {
    #         local_subdirname_modified = paste(local_subdirname,j_index,sep="")
    #         j_index = j_index+1
    #       }
    #       dir.create(local_subdirname_modified)
    #       local_subdirname = paste(local_subdirname_modified,"/",sep="")
    #     }
    #     setwd(local_subdirname)
    
    if (!(file.exists(local_subdirname)))
      dir.create(local_subdirname)
    setwd(local_subdirname)
    
    directory_file = "Directory.txt"
    write(paste(local_subdirname,filename,sep=""),directory_file)
    
    if (cluster)
    {  
      command = paste0("cp ",cluster_subdirname,filename," ",local_subdirname,filename)
      system(command, intern=TRUE)
      cat("Loading",filename,"\n")
      load(filename)
    }
    
    if (existingresult)
    {
      cat("Loading",filename,"\n")
      load(filename)
    }
  } 
}

for (par_index in 1:mpar_range)
{
  if (mpar)
    cat("Parameter value",par_index,"\n")  
  
  if (mpar) 
  {
    if (mnb_topics)
      nb_topics = nb_topics_range[par_index]  
  }
  if (alpha_insert == "alpha50_over_nb_topics")
    alpha = 50/nb_topics
  else if (alpha_insert == "alpha0.1")
    alpha = 0.1
  
  #   if (blei)
  #   {
  #     setwd(paste("/Users/guilhemsommeria-klein/Desktop/These/",data_insert,"/",barcode_insert,"/Blei_lda_alphaes",alpha,"_topics",nb_topics_range[i],"_norep_strset_3/"),sep="")
  #     
  #     #assignment <- read.table("word-assignments.dat", sep=" ")
  #     logbeta <- read.table("final.beta", sep=" ")
  #     logbeta = logbeta[,-c(1,2)]
  #     gamma <- read.table("final.gamma", sep=" ")
  #     alpha <- read.table("final.other", sep=" ")
  #     
  #     alpha_value = alpha[3,2]
  #     nb_terms = alpha[2,2]-1
  #     nb_doc = length(gamma[,1])
  #     
  #     documents = as.matrix(gamma-alpha_value)
  #     documents[documents<0] = 0
  #   }
  
  # creating the directories and downloading the .Rdata file from the cluster for (!mpar) condition 
  ######################
  
  if (!mpar)
  {
    if (cluster)
      cat("Loading data from cluster ...\n")
    
    if (Rtopicmodels_Gibbs || Rtopicmodels_VEM)
    {
      if (!best && !best_keep)
      {  
        
        if (Rtopicmodels_Gibbs) {
          filename = paste0(filename_insert,"_alpha",alpha,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,remove_single_sites_insert,".Rdata")
        } else if (Rtopicmodels_VEM)
          filename = paste0(filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,remove_single_sites_insert,".Rdata")
        
        directory_file = "Directory.txt"
        write(paste(local_dirname,filename,sep=""),directory_file)
        
        if (cluster)
        {  
          command = paste("cp ",cluster_dirname,filename," ",local_dirname,filename,sep="")
          system(command, intern=TRUE)
          cat("Loading",filename,"\n")
          load(filename)
        }
        
        if (existingresult)
        {
          cat("Loading",filename,"\n")
          load(filename)
        }
        
      } else if (best || best_keep)
      {
        if (best) 
        {
          if (Rtopicmodels_Gibbs)
          {
            local_subdirname = paste0(local_dirname,filename_insert,"_alpha",alpha,"_delta",delta,"_nb_topics",nb_topics,"_nb_iter",nb_iter,"_nb_real",nb_real,"_best",occurrence_insert,remove_single_sites_insert,no_rare_insert,"/")
            cluster_subdirname = paste0(cluster_dirname,filename_insert,"_alpha",alpha,"_delta",delta,"_nb_topics",nb_topics,"_nb_iter",nb_iter,"_nb_real",nb_real,"_best",occurrence_insert,remove_single_sites_insert,no_rare_insert,"/")
            filename = paste0(filename_insert,"_alpha",alpha,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,"_best",occurrence_insert,remove_single_sites_insert,no_rare_insert,".Rdata")  
          } else if (Rtopicmodels_VEM)
          {
            local_subdirname = paste0(local_dirname,filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best",occurrence_insert,remove_single_sites_insert,no_rare_insert,"/")
            cluster_subdirname = paste0(cluster_dirname,filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best",occurrence_insert,remove_single_sites_insert,no_rare_insert,"/")
            filename = paste0(filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best",occurrence_insert,remove_single_sites_insert,no_rare_insert,".Rdata")
          }
          if (!(file.exists(local_subdirname)))
            dir.create(local_subdirname)
          setwd(local_subdirname) 
          
          if (cluster)
          {  
            command = paste("cp ",cluster_subdirname,filename," ",local_subdirname,filename,sep="")
            system(command, intern=TRUE)
            cat("Loading",filename,"\n")
            load(filename)
          }
          
          if (existingresult)
          {
            cat("Loading",filename,"\n")
            load(filename)
          }
          
          directory_file = "Directory.txt"
          write(paste(local_subdirname,filename,sep=""),directory_file)
          
          if (local)
            cat(filename,"\n")
          
        } else if (best_keep)
        {
          if (Rtopicmodels_Gibbs)
          {
            local_subdirname = paste0(local_dirname,filename_insert,"_alpha",alpha,"_delta",delta,"_topics",nb_topics,"_nb_iter",nb_iter,"_nb_real",nb_real,"_best_keep",occurrence_insert,remove_single_sites_insert,no_rare_insert,"/")
            cluster_subdirname = paste0(cluster_dirname,filename_insert,"_alpha",alpha,"_delta",delta,"_topics",nb_topics,"_nb_iter",nb_iter,"_nb_real",nb_real,"_best_keep",occurrence_insert,remove_single_sites_insert,no_rare_insert,"/")
            filename = paste0(filename_insert,"_alpha",alpha,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,"_best_keep",occurrence_insert,remove_single_sites_insert,no_rare_insert,".Rdata")
          } else if (Rtopicmodels_VEM)
          {
            local_subdirname = paste0(local_dirname,filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best_keep",occurrence_insert,remove_single_sites_insert,no_rare_insert,"/")
            cluster_subdirname = paste0(cluster_dirname,filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best_keep",occurrence_insert,remove_single_sites_insert,no_rare_insert,"/")
            filename = paste0(filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best_keep",occurrence_insert,remove_single_sites_insert,no_rare_insert,".Rdata")
          }
          if (!(file.exists(local_subdirname)))
            dir.create(local_subdirname)
          setwd(local_subdirname) 
          
          if (cluster)
          {  
            command = paste("cp ",cluster_subdirname,filename," ",local_subdirname,filename,sep="")
            system(command, intern=TRUE)
            cat("Loading",filename,"\n")
            load(filename)
          }
          
          if (existingresult)
          {
            cat("Loading",filename,"\n")
            load(filename)
          }
          
          directory_file = "Directory.txt"
          write(paste(local_subdirname,filename,sep=""),directory_file)
          
          if (local) 
            cat(filename,"\n")
        }
      }
    } else if (hdp_faster)
    {
      local_subdirname = paste0(local_dirname,filename_insert,"hdp-faster_alpha",alpha,"_gamma",gamma,"_eta",eta,"_maxiter",maxiter,"_norep",occurrence_insert,remove_single_sites_insert,no_rare_insert)
      if (!(file.exists(local_subdirname)))
        dir.create(local_subdirname)
      setwd(local_subdirname) 
      
      read.table("")
    }
    
    ###### end of the !mpar condition
  } 
  
  #########################
  
  # beta = #topics x #terms
  # gamma = # documents x #topics
  if (Rtopicmodels_Gibbs)
  {
    if ((cluster || existingresult) && mpar)
    {
      for (j in 1:nb_real)
      {
        if (j==1)
          Result = Result_mpar[[(par_index-1)*nb_real+j]]
        else if (j>1)
          Result = c(Result,Result_mpar[[(par_index-1)*nb_real+j]])
      }
    } else if (local)
    {
      SEED = vector(length=nb_real,mode="integer")
      for (j in 1:nb_real)
        SEED[j] = as.integer(Sys.time()) - j*10^7
      control_LDA_Gibbs = list(alpha=alpha, estimate.beta=TRUE,
                               verbose = 100, prefix = tempfile(), save = 0, keep = llh_keep,
                               seed = SEED, nstart = nb_real, best = best,
                               delta = delta,
                               iter = nb_iter, burnin = 0, thin = nb_iter)
      
      Result = topicmodels::LDA(x=t(data2m),k=nb_topics,method = "Gibbs",control=control_LDA_Gibbs,model=NULL)
      
      if (mpar)
      {
        if (par_index == 1)
          Result_mpar = Result
        else if (par_index > 1)
          Result_mpar = c(Result_mpar,Result)
      }
    }
  } else if (Rtopicmodels_VEM)
  {
    if ((cluster || existingresult) && mpar)
    {
      for (j in 1:nb_real)
      {
        if (j==1)
          Result = Result_mpar[[(par_index-1)*nb_real+j]]
        else if (j>1)
          Result = c(Result,Result_mpar[[(par_index-1)*nb_real+j]])
      }
    } else if (local)
    {
      SEED = vector(length=nb_real,mode="integer")
      for (j in 1:nb_real)
        SEED[j] = as.integer(Sys.time()) - j*10^7
      control_LDA_VEM = list(estimate.alpha=TRUE, alpha=alpha, estimate.beta=TRUE,
                             verbose = 1, prefix = tempfile(), save = 0, keep = 1,
                             seed = SEED, nstart = nb_real, best = best,
                             var = list(iter.max = 500, tol = var_tol),
                             em = list(iter.max = 1000, tol = em_tol),
                             initialize = "random")
      
      Result = topicmodels::LDA(x=t(data2m),k=nb_topics,method = "VEM",control=control_LDA_VEM,model=NULL)
      
      #       ##########################
      #       library(topicmodels)
      #       nb_real = 5
      #       SEED = vector(length=nb_real,mode="integer")
      #       for (j in 1:nb_real)
      #         SEED[j] = as.integer(Sys.time()) - j*10^7
      #       alpha = 0.1
      #       nb_topics = 2
      #       em_tol = 10^-7
      #       var_tol = 10^-8
      #       #data2m = t(matrix(nrow=3,ncol=3,data=c(c(3,2,1),c(4,0,3),c(1,3,0))))
      #       data2m = t(matrix(nrow=3,ncol=3,data=c(c(4,2,0),c(1,2,3),c(0,1,3))))
      #       control_LDA_VEM = list(estimate.alpha=TRUE, alpha=alpha, estimate.beta=TRUE,
      #                              verbose = 1, prefix = tempfile(), save = 0, keep = 1,
      #                              seed = SEED, nstart = nb_real, best = 0,
      #                              var = list(iter.max = 500, tol = var_tol),
      #                              em = list(iter.max = 1000, tol = em_tol),
      #                              initialize = "random")
      #       illustrative_result = topicmodels::LDA(x=t(data2m),k=nb_topics,method = "VEM",control=control_LDA_VEM,model=NULL)
      #       LLH_final_real1 = vector(length=nb_real,mode="numeric")
      #       for (j in 1:nb_real)
      #         LLH_final_real1[j] = sum(illustrative_result[[j]]@loglikelihood)
      #       Ordered_realizations = sort.int(LLH_final_real1,decreasing=T,index.return=T)
      #       illustrative_documents = illustrative_result[[Ordered_realizations$ix[1]]]@gamma
      #       illustrative_beta = exp(illustrative_result[[Ordered_realizations$ix[1]]]@beta) 
      #       illustrative_alpha = illustrative_result[[Ordered_realizations$ix[1]]]@alpha
      #       #########################
      
      if (!best)
      {
        discarded_real = 1
        while (length(discarded_real)>0)
        {
          discarded_real = NULL
          for (j in 1:nb_real)
          { 
            if (max(Result[[j]]@logLiks)>=0)
            {
              if (length(discarded_real)>0)
                discarded_real = c(discarded_real,j)
              else
                discarded_real = j
              cat("Positive loglikelihood value: realization",j,"computed again\n")
            }
          }
          
          if (length(discarded_real)>0)
          { 
            SEED = vector(length=length(discarded_real),mode="integer")
            jj = 1
            for (j in discarded_real)
            {
              SEED[jj] = as.integer(Sys.time()) - j*10^7
              jj = jj+1
            }
            control_LDA_VEM = list(estimate.alpha=TRUE, alpha=alpha, estimate.beta=TRUE,
                                   verbose = 1, prefix = tempfile(), save = 0, keep = 1,
                                   seed = SEED, nstart = length(discarded_real), best = best,
                                   var = list(iter.max = 500, tol = var_tol),
                                   em = list(iter.max = 1000, tol = em_tol),
                                   initialize = "random")
            
            Result_add = topicmodels::LDA(x=t(data2m),k=nb_topics,method = "VEM",control=control_LDA_VEM,model=NULL)
            
            jj = 1
            for (j in discarded_real)
            {
              Result[[j]] = Result_add[[jj]]
              jj = jj+1
            }
          }
        }
      }
      
      if (mpar)
      {
        if (par_index == 1)
          Result_mpar = Result
        else if (par_index > 1)
          Result_mpar = c(Result_mpar,Result)
      }
    }
  }
  
  if (mpar)
  {
    if (Rtopicmodels_Gibbs)
    {
      nb_terms = Result[[1]]@fitted[[1]]@wordassignments$ncol
      nb_doc = Result[[1]]@fitted[[1]]@wordassignments$nrow
      nb_words = Result[[1]]@fitted[[1]]@n
      for (j in 1:nb_real)
      {
        LLH_final1[j,par_index] = Result[[j]]@fitted[[1]]@loglikelihood
        AIC1[j,par_index] = 2*(nb_topics*(nb_terms+nb_doc) - Result[[j]]@fitted[[1]]@loglikelihood)
        AIC2[j,par_index] = 2*(nb_words - Result[[j]]@fitted[[1]]@loglikelihood)
        AIC3[j,par_index] = 2*(nb_topics*nb_terms + 1 - Result[[j]]@fitted[[1]]@loglikelihood)
      } 
      LLH_final0[1,par_index] = mean(LLH_final1[,par_index])
      LLH_final0[2,par_index] = sd(LLH_final1[,par_index])      
      AIC0[1,par_index] = 2*(nb_topics*(nb_terms+nb_doc) - LLH_final0[1,par_index])
      AIC0[3,par_index] = 2*(nb_words - LLH_final0[1,par_index])
      AIC0[5,par_index] = 2*(nb_topics*nb_terms + 1 - LLH_final0[1,par_index])
      AIC0[2,par_index] = sd(AIC1[,par_index])
      AIC0[4,par_index] = sd(AIC2[,par_index])
      AIC0[6,par_index] = sd(AIC3[,par_index])
    } else if (Rtopicmodels_VEM)
    {
      nb_terms = Result[[1]]@wordassignments$ncol
      nb_doc = Result[[1]]@wordassignments$nrow
      nb_words = Result[[1]]@n
      for (j in 1:nb_real)
      {
        LLH_final1[j,par_index] = sum(Result[[j]]@loglikelihood)
        AIC1[j,par_index] = 2*(nb_topics*(nb_terms+nb_doc) - sum(Result[[j]]@loglikelihood))
        AIC2[j,par_index] = 2*(nb_words - sum(Result[[j]]@loglikelihood))
        AIC3[j,par_index] = 2*(nb_topics*nb_terms + 1 - sum(Result[[j]]@loglikelihood))
        alpha_est1[j,par_index] = Result[[j]]@alpha
      } 
      LLH_final0[1,par_index] = mean(LLH_final1[,par_index])
      LLH_final0[2,par_index] = sd(LLH_final1[,par_index])
      AIC0[1,par_index] = 2*(nb_topics*(nb_terms+nb_doc) - LLH_final0[1,par_index])
      AIC0[3,par_index] = 2*(nb_words - LLH_final0[1,par_index])
      AIC0[5,par_index] = 2*(nb_topics*nb_terms + 1 - LLH_final0[1,par_index])
      AIC0[2,par_index] = sd(AIC1[,par_index])
      AIC0[4,par_index] = sd(AIC2[,par_index])
      AIC0[6,par_index] = sd(AIC3[,par_index])
      alpha_est0[1,par_index] = mean(alpha_est1[,par_index])
      alpha_est0[2,par_index] = sd(alpha_est1[,par_index])
    }
    
    if (Rtopicmodels_VEM && (par_index == length(mpar_range)))
    {
      alpha_file = "estimated_alpha.txt"
      write(paste("par =",mpar_range[1],"- alpha =",alpha_est0[1,1]),alpha_file,ncolumns=1)
      if (length(mpar_range)>1)
      {
        for (par_index in 2:length(mpar_range))
          write(paste("par =",mpar_range[par_index],"- alpha =",alpha_est0[1,par_index]),alpha_file,append=T)
      }
    }
    # This (!mpar) section could be put after the the loop over mpar
  } else if (!mpar)
  {
    if (best_keep)
    {
      #AIC1 = vector(length=nb_real,mode="numeric")
      #AIC2 = vector(length=nb_real,mode="numeric")
      LLH_final_real1 = vector(length=nb_real,mode="numeric")
      #Akaike_weights_AIC1 = vector(length=nb_real,mode="numeric")
      #Akaike_weights_AIC2 = vector(length=nb_real,mode="numeric")
      Akaike_weights_llh = vector(length=nb_real,mode="numeric")
      for (j in 1:nb_real)
      {
        if (Rtopicmodels_VEM)
          LLH_final_real1[j] = sum(Result[[j]]@loglikelihood)
        else if (Rtopicmodels_Gibbs)
          LLH_final_real1[j] = Result[[j]]@fitted[[1]]@loglikelihood
        #AIC1[j] = 2*(nb_topics*(nb_terms+nb_doc) - LLH_final_real1[j])
        #AIC2[j] = 2*(nb_words - LLH_final_real1[j])
      }
      Ordered_realizations = sort.int(LLH_final_real1,decreasing=T,index.return=T)
      saveRDS(Ordered_realizations,file="Ordered_realizations.rds")
      best_real = which(LLH_final_real1 == max(LLH_final_real1))
      for (j in 1:nb_real) 
      {
        Akaike_weights_llh[j] = exp(LLH_final_real1[j] - max(LLH_final_real1))/sum(exp(LLH_final_real1 - max(LLH_final_real1)))
        #Akaike_weights_AIC1[j] = exp(1/2*(min(AIC1) - AIC1[j]))/sum(exp(1/2*(min(AIC1) - AIC1)))
        #Akaike_weights_AIC2[j] = exp(1/2*(min(AIC2) - AIC2[j]))/sum(exp(1/2*(min(AIC2) - AIC2)))
      }
    }
    if (!best)
    {
      if (Rtopicmodels_Gibbs)
      {
        LLH1 = matrix(nrow=nb_real,ncol=nb_iter/llh_keep,data=0)
        LLH0 = matrix(nrow=2,ncol=nb_iter/llh_keep,data=0)
        for (j in 1:nb_real)
        {
          #logbeta1 = Result[[j]]@fitted[[1]]@beta
          #documents1[j] = Result[[j]]@fitted[[1]]@gamma
          #nb_terms1 = Result[[j]]@fitted[[1]]@wordassignments$ncol
          #nb_doc = Result[[j]]@fitted[[1]]@wordassignments$nrow
          LLH1[j,] = Result[[j]]@fitted[[1]]@logLiks
          #           LLH0[1,] = Result[[j]]@fitted[[1]]@logLiks/nb_real + LLH0[1,]
          #           LLH0[2,] = Result[[j]]@fitted[[1]]@logLiks^2/nb_real + LLH0[2,]
          #llh[j] = Result[[j]]@fitted[[1]]@loglikelihood
        }
        LLH0[1,] = apply(LLH1,2,mean)
        LLH0[2,] = apply(LLH1,2,sd)
        #         LLH0[2,] = LLH0[2,] - LLH0[1,]^2
        #         LLH0[2,] = sqrt(LLH0[2,])
      } else if (Rtopicmodels_VEM)
      {
        # different number of iterations for each realization
        #       LLH1 = matrix(nrow=nb_real,ncol=nb_iter/llh_keep,data=0)
        #       LLH0 = matrix(nrow=2,ncol=nb_iter/llh_keep,data=0)
        #       for (j in 1:nb_real)
        #       { 
        #         #logbeta1 = Result[[j]]@fitted[[1]]@beta
        #         #documents1[j] = Result[[j]]@fitted[[1]]@gamma
        #         #nb_terms1 = Result[[j]]@fitted[[1]]@wordassignments$ncol
        #         #nb_doc = Result[[j]]@fitted[[1]]@wordassignments$nrow
        #         LLH1[j,] = Result[[j]]@fitted[[1]]@logLiks
        #         LLH0[1,] = Result[[j]]@fitted[[1]]@logLiks/nb_real + LLH0[1,]
        #         LLH0[2,] = Result[[j]]@fitted[[1]]@logLiks^2/nb_real + LLH0[2,]
        #         #llh[j] = Result[[j]]@fitted[[1]]@loglikelihood
        #       }
        #       LLH0[2,] = LLH0[2,] - LLH0[1,]^2
        #       LLH0[2,] = sqrt(LLH0[2,])
      }
    }
  }
  
  #### end of the loop over mpar_range
}
# if (Rtopicmodels_VEM)
#   if (mpar) {
#     alpha_insert = paste("alpha",format(min(alpha_est0[1,]),digits=3),"-",format(max(alpha_est0[1,]),digits=3),sep="")
#   } else if (!mpar) {
#     if (!best)
#       alpha_insert = paste("alpha",format(alpha_est0[1,],digits=3),sep="")
#   }

write(paste("Proportion of OTUs removed before LDA decomposition:",Prop_OTU_removed),file="Removed_OTUs_and_reads.txt")
write(paste("Proportion of reads removed before LDA decomposition:",Prop_reads_removed),file="Removed_OTUs_and_reads.txt",append=T)

############################
# END OF LDA DECOMPOSITION #
############################


################################################
# ############################################ #
# ######## POSTPROCESSING AND FIGURES ######## #
# ############################################ #
################################################

if (!mpar)
{ 
  if (select_real && best_keep)
  {
    length_selected_real = length(Selected_real)
    alpha_est_mean = 0
  } else if ((!select_real && best_keep) || best)
    length_selected_real = 1
  
  if (select_real && (testdata || data_betadiv_pooled) && realization_comparison)
    Correlation_old = vector(mode="list",length=length_selected_real-1)
  
  #########################################################
  # Defining variables outside the loop over realizations #
  #########################################################
  
  if (select_real && realization_comparison)
  {
    KL = vector(mode="list",length=length_selected_real-1)
    DKL100 = vector(mode="list",length=length_selected_real-1)
    Correlation = vector(mode="list",length=length_selected_real-1)
    
    KL_allRealPairs = matrix(nrow=length_selected_real,ncol=length_selected_real,data=0)
    KL_allRealPairs[lower.tri(KL_allRealPairs,diag=T)] = NA
    KL_allRealPairs_byTopic = list()
    for (k in 1:nb_topics)
      KL_allRealPairs_byTopic[[k]] = matrix(nrow=length_selected_real,ncol=length_selected_real,data=NA)
    
    p_value_allRealPairs = matrix(nrow=length_selected_real,ncol=length_selected_real,data=0)
    p_value_allRealPairs[lower.tri(p_value_allRealPairs,diag=T)] = NA
    p_value_allRealPairs_byTopic = list()
    for (k in 1:nb_topics)
      p_value_allRealPairs_byTopic[[k]] = matrix(nrow=length_selected_real,ncol=length_selected_real,data=NA)
    
    KL_allRealPairs_w_rndzations = matrix(nrow=length_selected_real,ncol=length_selected_real,data=0)
    diag(KL_allRealPairs_w_rndzations) = NA
    KL_allRealPairs_randomized_byTopic = list()
    for (k in 1:nb_topics)
      KL_allRealPairs_randomized_byTopic[[k]] = matrix(nrow=length_selected_real,ncol=length_selected_real,data=NA)
    
    DKL100_allRealPairs = matrix(nrow=length_selected_real,ncol=length_selected_real,data=0)
    DKL100_allRealPairs[lower.tri(DKL100_allRealPairs,diag=T)] = NA
    DKL100_allRealPairs_byTopic = list()
    for (k in 1:nb_topics)
      DKL100_allRealPairs_byTopic[[k]] = matrix(nrow=length_selected_real,ncol=length_selected_real,data=NA)
    
    Var_KL_topic_comparison_randomized_allRealPairs = matrix(nrow=length_selected_real,ncol=length_selected_real,data=0)
    Var_KL_topic_comparison_randomized_allRealPairs[lower.tri(Var_KL_topic_comparison_randomized_allRealPairs,diag=T)] = NA
    Var_KL_topic_comparison_randomized_allRealPairs_byTopic = list()
    for (k in 1:nb_topics)
      Var_KL_topic_comparison_randomized_allRealPairs_byTopic[[k]] = matrix(nrow=length_selected_real,ncol=length_selected_real,data=NA)
    
    SES_allRealPairs = matrix(nrow=length_selected_real,ncol=length_selected_real,data=0)
    SES_allRealPairs[lower.tri(SES_allRealPairs,diag=T)] = NA
    SES_allRealPairs_byTopic = list()
    for (k in 1:nb_topics)
      SES_allRealPairs_byTopic[[k]] = matrix(nrow=length_selected_real,ncol=length_selected_real,data=NA)
    #     SES_samplewise_allRealPairs = matrix(nrow=length_selected_real,ncol=length_selected_real,data=0)
    #     SES_samplewise_allRealPairs[lower.tri(SES_samplewise_allRealPairs,diag=T)] = NA
    #     SES_samplewise_allRealPairs_byTopic = list()
    #     for (k in 1:nb_topics)
    #       SES_samplewise_allRealPairs_byTopic[[k]] = matrix(nrow=length_selected_real,ncol=length_selected_real,data=NA)
    
    Correlation_allRealPairs = matrix(nrow=length_selected_real,ncol=length_selected_real,data=0)
    Correlation_allRealPairs_byTopic = list()
    for (k in 1:nb_topics)
      Correlation_allRealPairs_byTopic[[k]] = matrix(nrow=length_selected_real,ncol=length_selected_real,data=NA)
    
    Topic_correspondence_to_best_real = list()
    Topic_correspondence_to_best_real[[1]] = seq(1,nb_topics,1)
    
    llh_differences_allRealPairs = matrix(nrow=length_selected_real,ncol=length_selected_real,data=NA)
    
    #KL_sym = vector(mode="list",length=length_selected_real-1)
    documents_allreal = vector(mode="list",length=length_selected_real)
    KL_documents_allreal = vector(mode="list",length=length_selected_real)
    topic_compo_allreal = vector(mode="list",length=length_selected_real)
    KL_topic_compo_allreal = vector(mode="list",length=length_selected_real)
    sort_normal_topic_allreal = vector(mode="list",length=length_selected_real)
  }
  
  ##########################
  # Loading earth map data #
  ##########################
  
  if (kriging && data_gs)
  {
    borders = readOGR(dsn=paste0(local_prefix,data_insert), layer="World_car")
    borders@data$id = rownames(borders@data)
    bordersPoints = fortify(borders, region="id")
    bordersGgplot = merge(bordersPoints, borders@data, by="id")
    
    rivers = readOGR(dsn=paste0(local_prefix,data_insert,"/ne_50m_rivers_lake_centerlines"), layer="ne_50m_rivers_lake_centerlines")
    rivers@data$id = rownames(rivers@data)
    riversPoints = fortify(rivers, region="id")
    riversGgplot = merge(riversPoints, rivers@data, by="id")
    
    #       ocean = readOGR(dsn=paste0(local_prefix,data_insert,"/ne_50m_ocean"), layer="ne_50m_ocean")
    #       ocean@data$id = rownames(ocean@data)
    #       oceanPoints = fortify(ocean, region="id")
    #       oceanGgplot = merge(oceanPoints, ocean@data, by="id")  
    #     Americas = oceanGgplot[which(oceanGgplot$piece==2),]
    
    land = readOGR(dsn=paste0(local_prefix,data_insert,"/ne_10m_land"), layer="ne_10m_land")
    #land_shape = readShapePoly(paste0(local_prefix,data_insert,"/ne_10m_land/ne_10m_land"))
    
    #     land@data$id = rownames(land@data)
    #     landPoints = fortify(land, region="id")
    #     landGgplot = merge(landPoints, land@data, by="id") 
  }
  
  ##########################
  # LOOP OVER REALIZATIONS #
  ##########################
  for (j_select in 1:length_selected_real)
  {
    
    #########################################################
    # Defining variables within the loop over realizations #
    #########################################################
{
  
  if (best)
  {
    logbeta = Result@beta 
    documents = Result@gamma
    nb_terms = Result@wordassignments$ncol
    nb_doc = Result@wordassignments$nrow
    if (Rtopicmodels_Gibbs) {
      llh = Result@loglikelihood
      if (llh_keep)
        llh_iterations = Result@logLiks
    } else if (Rtopicmodels_VEM) 
    {
      llh = sum(Result@loglikelihood)
      llh_iterations = Result@logLiks
      alpha_est = Result@alpha
    }
  } else if (best_keep)
  { 
    if (select_real)
    {
      if (Rtopicmodels_VEM)
      {
        llh = LLH_final_real1[Ordered_realizations$ix[Selected_real[j_select]]] 
        llh_iterations = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@logLiks
        logbeta = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@beta 
        documents = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@gamma
        nb_terms = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@wordassignments$ncol
        nb_doc = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@wordassignments$nrow
        alpha_est = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@alpha
      } else if (Rtopicmodels_Gibbs)
      {
        llh = LLH_final_real1[Ordered_realizations$ix[Selected_real[j_select]]] 
        if (llh_keep)
          llh_iterations = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@fitted[[1]]@logLiks 
        logbeta = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@fitted[[1]]@beta 
        documents = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@fitted[[1]]@gamma
        nb_terms = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@fitted[[1]]@wordassignments$ncol
        nb_doc = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@fitted[[1]]@wordassignments$nrow
      }
    } else if (!select_real)
    {
      if (Rtopicmodels_VEM)
      {
        llh = LLH_final_real1[best_real] 
        llh_iterations = Result[[best_real]]@logLiks
        logbeta = Result[[best_real]]@beta 
        documents = Result[[best_real]]@gamma
        nb_terms = Result[[best_real]]@wordassignments$ncol
        nb_doc = Result[[best_real]]@wordassignments$nrow
        alpha_est = Result[[best_real]]@alpha
      } else if (Rtopicmodels_Gibbs)
      {
        llh = LLH_final_real1[best_real] 
        if (llh_keep)
          llh_iterations = Result[[best_real]]@fitted[[1]]@logLiks 
        logbeta = Result[[best_real]]@fitted[[1]]@beta 
        documents = Result[[best_real]]@fitted[[1]]@gamma
        nb_terms = Result[[best_real]]@fitted[[1]]@wordassignments$ncol
        nb_doc = Result[[best_real]]@fitted[[1]]@wordassignments$nrow
      }
    }
  }
  
  #### August 2014 code #####
  # nrow(data2m) = nb_seq et ncol(data2m) = nb_sites
  # nrow(documents) = nb_doc and ncol(documents) = nb_topics
  
  # saving the result of all realizations, without reordering the MOTUs as in topic_compo
  topic_compo = matrix(nrow=nb_terms,ncol=nb_topics,data=0)
  for (k in 1:nb_topics)
    topic_compo[,k] = exp(t(logbeta[k,]))
  KL_topic_compo = topic_compo
  KL_topic_compo[which(topic_compo < 1/sum(data2m))] = 1/sum(data2m)
  KL_topic_compo = sweep(KL_topic_compo,MARGIN=2,colSums(KL_topic_compo),`/`)   
  
  # setting the minimal proportion of a topic in a sample as 1/#readsInSample :
  #     for (i in 1:nrow(documents))
  #       documents[which(documents[i,]<1/sum(data2m[,i]))] = 1/sum(data2m[,i])
  # setting the minimal proportion of a topic in a sample as 1/total#reads (same as for the minimal proportion of an OTU in a sample)
  # in the case of occurrence data, it is the sum of OTU occurrences over the dataset (the total number of reads would be given by sum(count2) for occurrence data)
  KL_documents = documents
  KL_documents[which(documents < 1/sum(data2m))] = 1/sum(data2m)
  
  # documents : propotion of each topic in a site/document (sums to 1 over topics)
  # KL_documents : propotion of each site/document in a topic (sums to 1 over sites/documents)
  #     KL_documents = sweep(documents,MARGIN=2,colSums(documents),`/`)
  KL_documents = sweep(KL_documents,MARGIN=2,colSums(KL_documents),`/`)
  
  # smoothed_KL_documents = sweep(documents+1,MARGIN=2,colSums(documents)+nb_doc,`/`)
  # smoothed_documents = sweep(documents+1,MARGIN=1,rowSums(documents)+nb_topics,`/`)
  
  if (realization_comparison)
  {
    documents_allreal[[j_select]] = documents
    KL_documents_allreal[[j_select]] = KL_documents
    
    topic_compo_allreal[[j_select]] = topic_compo
    KL_topic_compo_allreal[[j_select]] = KL_topic_compo
  }
  
}

################################
# Topic dominance calculations #
################################

if (dominance_calculation)
{
  dominant_topic=vector(length=nb_doc,mode="integer")
  for (i in 1:nb_doc)
  {
    dominant_topic[i]=which(documents[i,]==max(documents[i,]))
  }
  spatial_dominant_topic = matrix(nrow=19,ncol=19)
  for (j in 1:19)
  {
    for (i in 1:19)
    {
      spatial_dominant_topic[i,j] = dominant_topic[(j-1)*19+i]
    }
  }
  
  # dominant topic per site with above-threshold proportion
  dominant_topic_thres=matrix(nrow=nb_doc,ncol=5)
  threshold=c(0.5,0.6,0.7,0.8,0.9)
  for (i in 1:nb_doc)
  {
    dom_topic=which(documents[i,]==max(documents[i,]))
    for (k in 1:5)
    {
      if (documents[i,dom_topic]>=threshold[k])
        dominant_topic_thres[i,k]=dom_topic
      else dominant_topic_thres[i,k]=0
    }
  }
  
  topic_site_nb=matrix(nrow=nb_topics,ncol=5,data=0)
  for (k in 1:5)
  {
    for (j in 1:nb_topics)
    {
      for (i in 1:nb_doc)
      {
        if (documents[i,j]>=threshold[k])
          topic_site_nb[j,k]=topic_site_nb[j,k]+1
      }
    }
  }
}

##################
# Topic ordering #
##################

if (topic_read_proportion_comput)
{
  #tot_reads=sum(sample_counts)
  tot_reads=sum(sample_count2)
  # Computing the proportion of total read number for each topic
  prop_topic = vector(length=nb_topics,mode="numeric")
  for (j in 1:nb_topics)
  { 
    for (i in 1:nb_doc)
    {
      prop_topic[j] =  prop_topic[j] + sample_count2[i]*documents[i,j]/tot_reads 
    }
  }
  sort_prop_topic = sort.int(prop_topic,index.return=T)
}

##########################################################################
# Computing the proportion of site-normalized read number for each topic #
##########################################################################

normal_topic = vector(length=nb_topics,mode="numeric")
for (j in 1:nb_topics)
{ 
  for (i in 1:nb_doc)
  {
    normal_topic[j] =  normal_topic[j] + documents[i,j]/nb_doc 
  }
}
sort_normal_topic = sort.int(normal_topic,index.return=T)
save(sort_normal_topic,file="sort_normal_topic.Rdata")
if (realization_comparison)
  sort_normal_topic_allreal[[j_select]] = sort_normal_topic

############################
# Saving topic composition #
############################
topic_compo_taxo = vector(length = nb_topics, mode = "list")
for (k in 1:nb_topics)
{
  # saving the result of all realizations, without reordering the MOTUs
  #topic_compo_allreal[[j_select]][,k] = exp(t(logbeta[k,]))
  
  # MOTUs are reordered by abundance within each site for the current realization
  vect=sort.int(t(logbeta[k,]),decreasing=T,index.return=T)
  proba=exp(vect$x)
  #mat=matrix(nrow=nb_terms,ncol=length(taxo_ref[1,]),dimnames=list(as.character(1:nb_terms),as.character(as.list(taxo_ref[1,]))))
  mat=matrix(nrow=nb_terms,ncol=length(taxo_ref[1,]))
  taxo_ref_names=as.character(taxo_ref[1,1])
  for (i in 2:length(taxo_ref[1,]))
    taxo_ref_names = c(taxo_ref_names,as.character(taxo_ref[1,i]))
  colnames(mat)=taxo_ref_names
  #topic_compo1=data.frame(vect,cbind(mat[1:(length(taxo_ref[1,])+1),]))
  #topic_compo1=data.frame(vect,mat)
  topic_compo_taxo[[k]]=data.frame(proba,mat)
  #     for (i in 1:nb_terms)
  #     {
  #     #topic_compo_index[i]=taxo_ref[i+1,1]
  #     topic_compo[[1,k]][i,2:(length(taxo_ref[1,])+1)]=as.character(taxo_ref[vect$ix[i]+1,])
  #     }
  for (i in 1:nb_terms)
  {
    #topic_compo_index[i]=taxo_ref[i+1,1]
    for (j in 1:length(taxo_ref[1,]))
      topic_compo_taxo[[k]][i,j+1]=as.character(taxo_ref[vect$ix[i]+1,j])
  }
}

#######################################
# Comparing realizations: computation #
#######################################

# Computing the correlation between topics (columns) over sites (lines):
# -> Redundant with Correlation_samplewise up to the ordering of the topics, but used for testdata=1 and data_betadiv_pooled=1
if ((j_select != 1) && (testdata || data_betadiv_pooled) && select_real && realization_comparison)
  Correlation_old[[j_select-1]] = cor(documents_allreal[[1]],documents)
# Correlation_samplewise[[j_select-1]] = cor(topic_compo_allreal[[1]],topic_compo_allreal[[j_select]][,Topic_correspondence_samplewise])

# Computing correlation between the topics of the highest-llh realization and the other
if ((j_select != 1) && realization_comparison)
{
  cat("Computing correlations between realizations 1 to",Selected_real[j_select]-1,"and realization",Selected_real[j_select],"...\n")  
  
  # Retrieving the correspondence between topics among realizations
  #     Topic_comparison_cor[j_select-1,] = apply(Correlation[[j_select-1]],1,max)
  #     Topic_correspondence_cor = vector(length=nb_topics,mode="numeric")
  #     for (k in 1:nb_topics)
  #       Topic_correspondence_cor[k] = which(Topic_comparison_cor[j_select-1,k] == Correlation[[j_select-1]][k,])
  
  # Using Hellinger dissimilarity between topics on documents-wise normalized topic decomposition (KL_norm_document sums to 1 over all documents for each topic (i.e. column))
  #     Similarity = 1 - 1/sqrt(2)*dist(t(sqrt(cbind(KL_documents_allreal[[1]][,rev(sort_normal_topic_allreal[[1]]$ix)],KL_documents[,rev(sort_normal_topic$ix)]))),
  #                                              method = "euclidean", diag = FALSE, upper = FALSE)
  #     Hellinger[[j_select-1]] = as.matrix(Similarity)[(nb_topics+1):(2*nb_topics),1:nb_topics]
  
  if (samplewise)
  {
    KL_documents_jselect_randomized = vector(length=nb_rndzations,mode="list")
  } else if (MOTUwise)
    KL_topic_compo_jselect_randomized = vector(length=nb_rndzations,mode="list")
  
  # Computing the mean KL distance between topics for all pairs of realizations
  for (j_real in 1:(j_select-1))
  {
    
    # Call to function LDA_realization_comparison_fun
    # KL_documents_jselect_randomized (if samplewise) or KL_topic_compo_jselect_randomized (if MOTUwise) are computed if j_real == 1 or just made use of is j_real>1
    LDA_realization_comparison_result = 
      LDA_realization_comparison_fun(j_select,j_real,ifelse(exists("KL_documents_jselect_randomized"),KL_documents_jselect_randomized,NA),
                                     ifelse(exists("KL_topic_compo_jselect_randomized"),KL_topic_compo_jselect_randomized,NA),
                                     nb_topics,samplewise,MOTUwise,documents_allreal,KL_documents_allreal,topic_compo_allreal,KL_topic_compo_allreal,
                                     testdata,filled,filled_with_gaps,nb_rndzations,ncol_data2m = ncol(data2m),nrow_data2m = nrow(data2m),
                                     Missing_positions_indices = ifelse(exists("Missing_positions_indices"),Missing_positions_indices,NA),bij)
    
    Topic_correspondence = LDA_realization_comparison_result$Topic_correspondence
    KL_topic_comparison = LDA_realization_comparison_result$KL_topic_comparison
    Mean_KL_topic_comparison_randomized = LDA_realization_comparison_result$Mean_KL_topic_comparison_randomized
    p_value = LDA_realization_comparison_result$p_value
    Var_KL_topic_comparison_randomized = LDA_realization_comparison_result$Var_KL_topic_comparison_randomized
    if (samplewise)
      KL_documents_jselect_randomized = LDA_realization_comparison_result$KL_documents_jselect_randomized
    else if (MOTUwise)
      KL_topic_compo_jselect_randomized = LDA_realization_comparison_result$KL_topic_compo_jselect_randomized
    
    ##############################
    # Saving the KL distance between the current and the best realizations, for all topics ordered by best KL-based matches
    if (j_real == 1)
    {
      Topic_correspondence_to_best_real[[j_select]] = Topic_correspondence
      
      if (bij)
      {
        KL[[j_select-1]] = KL_topic_comparison[,Topic_correspondence]
        
        if (samplewise)
          Correlation[[j_select-1]] = cor(documents_allreal[[1]],documents_allreal[[j_select]][,Topic_correspondence])
        else if (MOTUwise)
          Correlation[[j_select-1]] = cor(topic_compo_allreal[[1]],topic_compo_allreal[[j_select]][,Topic_correspondence])
        
        DKL100[[j_select-1]] = vector(length=nb_topics)
        for (k in 1:nb_topics)
          DKL100[[j_select-1]][k] = Mean_KL_topic_comparison_randomized[k] - KL_topic_comparison[k,Topic_correspondence[k]]
      } else if (!bij)
      {
        # ! These quantities are not really meaningful in the !bij case !
        KL[[j_select-1]] = KL_topic_comparison[Topic_correspondence[,1],Topic_correspondence[,2]]
        
        if (samplewise)
          Correlation[[j_select-1]] = cor(documents_allreal[[1]][,Topic_correspondence[,1]],documents_allreal[[j_select]][,Topic_correspondence[,2]])
        else if (MOTUwise)
          Correlation[[j_select-1]] = cor(topic_compo_allreal[[1]][,Topic_correspondence[,1]],topic_compo_allreal[[j_select]][,Topic_correspondence[,2]])
        
        DKL100[[j_select-1]] = vector(length=nb_topics)
        for (i in 1:nb_topics)
          DKL100[[j_select-1]][i] = Mean_KL_topic_comparison_randomized[i] - KL_topic_comparison[Topic_correspondence[i,1],Topic_correspondence[i,2]]
      }
    }
    
    # Saving the KL distance between the current and all previous realizations, averaged over best matching topics
    for (k in 1:nb_topics)
    {
      if (bij)
      {
        # The value of k indexes the topics of j_real in their original order
        # Topic_correspondence_samplewise_to_best_real[[j_real]][k] contains the order of the topics in j_real with respect to the topics of the best realization
        p_value_allRealPairs_byTopic[[k]][j_real,j_select] = p_value[k]
        Var_KL_topic_comparison_randomized_allRealPairs_byTopic[[k]][j_real,j_select] = Var_KL_topic_comparison_randomized
        SES_allRealPairs_byTopic[[k]][j_real,j_select] = (Mean_KL_topic_comparison_randomized[k] - KL_topic_comparison[k,Topic_correspondence[k]])/
          sqrt(Var_KL_topic_comparison_randomized[k] - Mean_KL_topic_comparison_randomized[k]^2)
        Var_KL_topic_comparison_randomized_allRealPairs[j_real,j_select] = Var_KL_topic_comparison_randomized[k]/nb_topics + 
          Var_KL_topic_comparison_randomized_allRealPairs[j_real,j_select]
        
        #KL_MOTUwise_allRealPairs[j_real,j_select] = KL_MOTUwise_topic_comparison[k,Topic_correspondence_MOTUwise[k]]/nb_topics + KL_MOTUwise_allRealPairs[j_real,j_select]
        KL_allRealPairs[j_real,j_select] = KL_topic_comparison[k,Topic_correspondence[k]]/nb_topics + KL_allRealPairs[j_real,j_select]
        KL_allRealPairs_byTopic[[k]][j_real,j_select] = KL_topic_comparison[k,Topic_correspondence[k]]
        KL_allRealPairs_randomized_byTopic[[k]][j_real,j_select] = Mean_KL_topic_comparison_randomized[k]
        
        KL_allRealPairs_w_rndzations[j_real,j_select] = KL_topic_comparison[k,Topic_correspondence[k]]/nb_topics + KL_allRealPairs_w_rndzations[j_real,j_select]
        KL_allRealPairs_w_rndzations[j_select,j_real] = Mean_KL_topic_comparison_randomized[k]/nb_topics + KL_allRealPairs_w_rndzations[j_select,j_real]
        
        DKL100_allRealPairs[j_real,j_select] = (Mean_KL_topic_comparison_randomized[k] - KL_topic_comparison[k,Topic_correspondence[k]])/nb_topics + 
          DKL100_allRealPairs[j_real,j_select]
        DKL100_allRealPairs_byTopic[[k]][j_real,j_select] = 
          Mean_KL_topic_comparison_randomized[k] - KL_topic_comparison[k,Topic_correspondence[k]]
        
        SES_allRealPairs[j_real,j_select] = SES_allRealPairs_byTopic[[k]][j_real,j_select]/nb_topics + SES_allRealPairs[j_real,j_select]
        
        #             SES_samplewise_allRealPairs_byTopic[[k]][j_real,j_select] = 
        #               DKL100_samplewise_allRealPairs_byTopic[[k]][j_real,j_select]/sqrt(Var_KL_samplewise_topic_comparison_randomized[k] - Mean_KL_samplewise_topic_comparison_randomized[k]^2)
        
        #Correlation_MOTUwise_allRealPairs[j_real,j_select] = cor(documents_allreal[[j_real]],documents[,Topic_correspondence_MOTUwise])[k,k]
        
        if (samplewise)
        {
          Correlation_allRealPairs[j_real,j_select] = cor(documents_allreal[[j_real]],documents_allreal[[j_select]][,Topic_correspondence])[k,k]/nb_topics +
            Correlation_allRealPairs[j_real,j_select]
          Correlation_allRealPairs_byTopic[[k]][j_real,j_select] = 
            cor(documents_allreal[[j_real]],documents_allreal[[j_select]][,Topic_correspondence])[k,k]
        } else if (MOTUwise)
        {
          Correlation_allRealPairs[j_real,j_select] = cor(topic_compo_allreal[[j_real]],topic_compo_allreal[[j_select]][,Topic_correspondence])[k,k]/nb_topics +
            Correlation_allRealPairs[j_real,j_select]
          Correlation_allRealPairs_byTopic[[k]][j_real,j_select] = 
            cor(topic_compo_allreal[[j_real]],topic_compo_allreal[[j_select]][,Topic_correspondence])[k,k]
        }
      } else if (!bij)
      {
        k0 = Topic_correspondence[k,1]
        k1 = Topic_correspondence[k,2]
        
        p_value_allRealPairs_byTopic[[k]][j_real,j_select] = p_value[k]
        Var_KL_topic_comparison_randomized_allRealPairs_byTopic[[k]][j_real,j_select] = Var_KL_topic_comparison_randomized
        SES_allRealPairs_byTopic[[k]][j_real,j_select] = (Mean_KL_topic_comparison_randomized[k] - KL_topic_comparison[k0,k1])/
          sqrt(Var_KL_topic_comparison_randomized - Mean_KL_topic_comparison_randomized[k]^2)
        Var_KL_topic_comparison_randomized_allRealPairs[j_real,j_select] = Var_KL_topic_comparison_randomized/nb_topics + 
          Var_KL_topic_comparison_randomized_allRealPairs[j_real,j_select]
        
        #KL_MOTUwise_allRealPairs[j_real,j_select] = KL_MOTUwise_topic_comparison[k,Topic_correspondence_MOTUwise[k]]/nb_topics + KL_MOTUwise_allRealPairs[j_real,j_select]
        KL_allRealPairs[j_real,j_select] = KL_topic_comparison[k0,k1]/nb_topics + KL_allRealPairs[j_real,j_select]
        KL_allRealPairs_byTopic[[k]][j_real,j_select] = KL_topic_comparison[k0,k1]
        KL_allRealPairs_randomized_byTopic[[k]][j_real,j_select] = Mean_KL_topic_comparison_randomized[k]
        
        KL_allRealPairs_w_rndzations[j_real,j_select] = KL_topic_comparison[k0,k1]/nb_topics + KL_allRealPairs_w_rndzations[j_real,j_select]
        KL_allRealPairs_w_rndzations[j_select,j_real] = Mean_KL_topic_comparison_randomized[k]/nb_topics + KL_allRealPairs_w_rndzations[j_select,j_real]
        
        DKL100_allRealPairs[j_real,j_select] = (Mean_KL_topic_comparison_randomized[k] - KL_topic_comparison[k0,k1])/nb_topics + 
          DKL100_allRealPairs[j_real,j_select]
        DKL100_allRealPairs_byTopic[[k]][j_real,j_select] = 
          Mean_KL_topic_comparison_randomized[k] - KL_topic_comparison[k0,k1]
        
        SES_allRealPairs[j_real,j_select] = SES_allRealPairs_byTopic[[k]][j_real,j_select]/nb_topics + SES_allRealPairs[j_real,j_select]
        
        #           SES_samplewise_allRealPairs_byTopic[[k]][j_real,j_select] = 
        #             DKL100_samplewise_allRealPairs_byTopic[[k]][j_real,j_select]/sqrt(Var_KL_samplewise_topic_comparison_randomized[k] - Mean_KL_samplewise_topic_comparison_randomized[k]^2)
        
        #Correlation_MOTUwise_allRealPairs[j_real,j_select] = cor(documents_allreal[[j_real]],documents[,Topic_correspondence_MOTUwise])[k,k]
        
        if (samplewise)
        {
          Correlation_allRealPairs[j_real,j_select] = cor(documents_allreal[[j_real]][,Topic_correspondence[,1]],documents_allreal[[j_select]][,Topic_correspondence[,2]])[k,k]/nb_topics +
            Correlation_allRealPairs[j_real,j_select]
          Correlation_allRealPairs_byTopic[[k]][j_real,j_select] = 
            cor(documents_allreal[[j_real]][,Topic_correspondence[,1]],documents_allreal[[j_select]][,Topic_correspondence[,2]])[k,k]
        } else if (MOTUwise)
        {
          Correlation_allRealPairs[j_real,j_select] = cor(topic_compo_allreal[[j_real]][,Topic_correspondence[,1]],topic_compo_allreal[[j_select]][,Topic_correspondence[,2]])[k,k]/nb_topics +
            Correlation_allRealPairs[j_real,j_select]
          Correlation_allRealPairs_byTopic[[k]][j_real,j_select] = 
            cor(topic_compo_allreal[[j_real]][,Topic_correspondence[,1]],topic_compo_allreal[[j_select]][,Topic_correspondence[,2]])[k,k]
        }
      }
    }
    
    llh_differences_allRealPairs[j_real,j_select] = Ordered_realizations$x[j_real]-Ordered_realizations$x[j_select]
    p_value_allRealPairs[j_real,j_select] = mean(p_value)
    
    #         SES_samplewise_allRealPairs[j_real,j_select] = DKL100_samplewise_allRealPairs[j_real,j_select]/
    #             sqrt(mean(Var_KL_samplewise_topic_comparison_randomized)-mean(Mean_KL_samplewise_topic_comparison_randomized)^2)
    
    # End of the j_real loop
  }   
}

##################################################
# Comparison with abiotic variables: computation #
##################################################

# Computing correlations between the best realization and abiotic variables
if ((j_select == 1) && abiotic_variables)
{
  cat("Computing correlations between the best realization and abiotic variables ...\n")
  
  source(paste0(code_insert,"LDA_abiotic_variables_fun.R"))
  
  if (data_pp && !testdata)
    case_list = c("lidar","chemi")
  else if (data_gs)
    case_list = "climate"
  
  for (case in case_list)
  {
    # Call to the LDA_abiotic_variables_fun function
    LDA_abiotic_variables_result = LDA_abiotic_variables_fun(data_pp,testdata,data_gs,case,local_prefix,data_insert,coordGS = ifelse(exists("coordGS"),coordGS,NA),
                                                             filled,filled_with_gaps,Missing_positions_indices = ifelse(exists("Missing_positions_indices"),Missing_positions_indices,NA),
                                                             ncol_data2m = ncol(data2m),sum_data2m = sum(data2m),pca_abiotic,nb_topics,nb_abiotic_rndzations,
                                                             documents,KL_documents,occurrence,case_list)
    
    #         if (case == "chemi")
    #           browser()
    
    if ((case == case_list[1]) && data_pp && (nb_topics == 3))
    {
      SKL_bacteria = LDA_abiotic_variables_result$SKL_bacteria
      p_value_SKL_bacteria = LDA_abiotic_variables_result$p_value_SKL_bacteria
    }
    
    if (case == "lidar")
    {
      Correlation_lidar = LDA_abiotic_variables_result$Correlation_abiotic
      ncol_lidar = LDA_abiotic_variables_result$ncol0
      Mean_cor_lidar_comparison_randomized = LDA_abiotic_variables_result$Mean_cor_abiotic_comparison_randomized 
      Sd_cor_lidar_comparison_randomized = LDA_abiotic_variables_result$Sd_cor_abiotic_comparison_randomized
      p_value_lidar = LDA_abiotic_variables_result$p_value_abiotic 
      colnames_lidar = LDA_abiotic_variables_result$colnames_abiotic
    } else if (case == "chemi")
    {
      Correlation_chemi = LDA_abiotic_variables_result$Correlation_abiotic
      ncol_chemi = LDA_abiotic_variables_result$ncol0
      Mean_cor_chemi_comparison_randomized = LDA_abiotic_variables_result$Mean_cor_abiotic_comparison_randomized 
      Sd_cor_chemi_comparison_randomized = LDA_abiotic_variables_result$Sd_cor_abiotic_comparison_randomized
      p_value_chemi = LDA_abiotic_variables_result$p_value_abiotic 
      colnames_chemi = LDA_abiotic_variables_result$colnames_abiotic
    } else if (case == "climate")
    {
      Correlation_climate = LDA_abiotic_variables_result$Correlation_abiotic
      ncol_climate = LDA_abiotic_variables_result$ncol0
      #           Mean_cor_climate_comparison_randomized = Mean_cor_abiotic_comparison_randomized 
      #           Var_cor_climate_comparison_randomized = Var_cor_abiotic_comparison_randomized
      #           p_value_climate = p_value_abiotic 
      colnames_climate = LDA_abiotic_variables_result$colnames_abiotic
    }
    
    # End loop over case list :
  }
  
  # Resuming to the main subdirectory
  setwd(local_subdirname)
  
}

#########################################
# Goodness_of_fit testdata: computation #
#########################################

# Computing KL distance between true documents and the LDA topics, as well as between true documents and randomizations of the LDA topics
if (best_keep && (j_select == 1) && testdata && (nb_topics == true_nb_topics))
{
  source(paste0(code_insert,"LDA_testdata_fun.R"))
  
  # Call to the LDA_testdata_fun function
  LDA_testdata_result = LDA_testdata_fun(nb_rndzations_true_documents,samplewise,MOTUwise,ncol_data2m=ncol(data2m),
                                         nrow_data2m=nrow(data2m),nb_topics,KL_documents,documents,true_documents,KL_topic_compo)
  
  Correlation_true_documents = LDA_testdata_result$Correlation_true_documents
  Mean_KL_topic_comparison_true_documents_randomized = LDA_testdata_result$Mean_KL_topic_comparison_true_documents_randomized
  p_value_true_documents = LDA_testdata_result$p_value_true_documents
  KL_topic_comparison_true_documents = LDA_testdata_result$KL_topic_comparison_true_documents
  DKL100_true_documents = LDA_testdata_result$DKL100_true_documents
  nES_true_documents = LDA_testdata_result$nES_true_documents
  
  save(Correlation_true_documents,KL_topic_comparison_true_documents,Mean_KL_topic_comparison_true_documents_randomized,
       p_value_true_documents,nES_true_documents,file="Goodness-of-fit_true_documents.Rdata")
}

#########################################################################################
# Correlation to bacterial assemblages to name the assemblages in case (nb_topics == 3) #
#########################################################################################

if (j_select == 1)
{
  # assemblage_names_vect and correlations_to_bacterial_assemblages are reused to label the assemblages
  # for realization_comparison in the (data_pp && nb_topics == 3) case
  assemblage_names_vect = vector(length=nb_topics,mode="character")
  correlations_to_bacterial_assemblages = vector(length=nb_topics,mode="numeric")
  if (data_pp && (nb_topics==3))
  {
    # Retrieving the name of the assemblage based on the correlation with the 3 bacterial assemblages:
    setwd(paste0(local_prefix,data_insert))
    documents_bacteria = readRDS("documents_bacteria_PP_3topics_best_real.rds")
    setwd(local_subdirname)
    if (!filled || filled_with_gaps)
      documents_bacteria = documents_bacteria[-which(Missing_positions_indices==1),]
    colnames(documents_bacteria) = c("Terra firme","Hydromorphic","Exposed rock")
    for (k in 1:nb_topics)
    {
      k0 = rev(sort_normal_topic$ix)[k]
      Correlation_bacteria = cor(documents[,k0],documents_bacteria)
      sorted_correlations_to_bacteria = sort.int(Correlation_bacteria,index=T,decreasing=T)
      k_bact = sorted_correlations_to_bacteria$ix[1]
      correlations_to_bacterial_assemblages[k] = Correlation_bacteria[k_bact]
      assemblage_names_vect[k] = colnames(documents_bacteria)[k_bact]
    }
  } else
    assemblage_names_vect = seq(1,nb_topics,1)
  
  diversity = vector(length = nb_topics+1, mode="numeric")
  for (k in 1:nb_topics)
  {
    k0 = rev(sort_normal_topic$ix)[k]
    diversity[k] = length(which(topic_compo[,k0] > 1/sum(data2m)))
  }
  diversity[nb_topics+1] = nrow(data2m)
  diversity = data.frame(Assemblage = c(assemblage_names_vect,"Total diversity"),Diversity = diversity)
  
}

################################################
# Topic comparison within the best realization #
################################################

if ((j_select == 1) && best_real_comparison)
{
  source(paste0(code_insert,"LDA_bestreal_fun.R"))
  
  LDA_best_real_result = LDA_best_real_fun(nb_rndzations_best_real,nrow_data2m = nrow(data2m),sort_normal_topic,KL_documents,
                                           KL_topic_compo,assemblage_names_vect,documents,topic_compo,sum_data2m = sum(data2m))
  
  Hellinger_topic_comparison_MOTUwise = LDA_best_real_result$Hellinger_topic_comparison_MOTUwise
  Jaccard_topic_comparison_MOTUwise = LDA_best_real_result$Jaccard_topic_comparison_MOTUwise
  Beta.sim_topic_comparison_MOTUwise = LDA_best_real_result$Beta.sim_topic_comparison_MOTUwise
  Corr_topic_comparison_samplewise = LDA_best_real_result$Corr_topic_comparison_samplewise
  Corr_topic_comparison_MOTUwise = LDA_best_real_result$Corr_topic_comparison_MOTUwise
  KL_topic_comparison_samplewise = LDA_best_real_result$KL_topic_comparison_samplewise
  KL_topic_comparison_MOTUwise = LDA_best_real_result$KL_topic_comparison_MOTUwise
  nES_topic_comparison_MOTUwise = LDA_best_real_result$nES_topic_comparison_MOTUwise
  
  #       cor(cophenetic(upgma_Hellinger2),Hellinger_topic_comparison_MOTUwise2)
  #       cor(cophenetic(upgma_Hellinger1),Hellinger_topic_comparison_MOTUwise1)
  #       cor(cophenetic(upgma_Hellinger),Hellinger_topic_comparison_MOTUwise)
  #       cor(cophenetic(upgma_Jaccard),Jaccard_topic_comparison_MOTUwise)
  
  #       coph_Hellinger = cophenetic(upgma)
  #       cor(spdist, coph)
  #       # 0.7432604
  # 
  #       wpgma <- agnes(spdist, diss =T, method = "weighted", keep.diss =F, keep.data =F)
  #       coph=cophenetic(wpgma)
  #       cor(spdist, coph)
  #       # 0.6194799
  # 
  #       ward <- agnes(spdist, diss =T, method = "ward", keep.diss =F, keep.data =F)
  #       coph=cophenetic(ward)
  #       cor(spdist, coph)
  #       # 0.6047672
  # 
  # sl <- agnes(spdist, diss =T, method = "single", keep.diss =F, keep.data =F)
  # coph=cophenetic(sl)
  # cor(spdist, coph)
  # # 0.1811903
  # 
  # cl<- agnes(spdist, diss =T, method = "complete", keep.diss =F, keep.data =F)
  # coph=cophenetic(cl)
  # cor(spdist, coph) 
  
  # Code Antoine :
  #       upgma_Hellinger2 <- agnes(Hellinger_topic_comparison_MOTUwise2, diss =T, method = "average", keep.diss =F, keep.data =F)
  #       upgma_Hellinger1 <- agnes(Hellinger_topic_comparison_MOTUwise1, diss =T, method = "average", keep.diss =F, keep.data =F)
  upgma_Hellinger <- agnes(Hellinger_topic_comparison_MOTUwise, diss =T, method = "average", keep.diss =F, keep.data =F)
  upgma_Jaccard <- agnes(Jaccard_topic_comparison_MOTUwise, diss =T, method = "average", keep.diss =F, keep.data =F)
  upgma_Beta.sim <- agnes(Beta.sim_topic_comparison_MOTUwise, diss =T, method = "average", keep.diss =F, keep.data =F)
  
  bestreal_dirname = paste0(local_subdirname,"Comparison_in_best_real")
  setwd(bestreal_dirname)
  
  pdf("Hierarchical_clustering_of_topics_Hellinger.pdf")
  plot(upgma_Hellinger, which.plots=2, ann=F, cex=2, cex.axis = 2, lwd=1.5)
  title("Average Hellinger")
  dev.off()
  
  pdf("Hierarchical_clustering_of_topics_Jaccard.pdf")
  plot(upgma_Jaccard, which.plots=2, ann=F, cex=2, cex.axis = 2, lwd=1.5)
  title("Average Jaccard")
  dev.off()
  
  pdf("Hierarchical_clustering_of_topics_beta.sim.pdf")
  plot(upgma_Beta.sim, which.plots=2, ann=F, cex=2, cex.axis = 2, lwd=1.5)
  title("Average beta.sim")
  dev.off()
  
  save(KL_topic_comparison_samplewise, KL_topic_comparison_MOTUwise, nES_topic_comparison_MOTUwise, Corr_topic_comparison_samplewise, Corr_topic_comparison_MOTUwise,
       Hellinger_topic_comparison_MOTUwise, Jaccard_topic_comparison_MOTUwise, Beta.sim_topic_comparison_MOTUwise, diversity, file = "Topic_comparison_inBestReal.Rdata")
  
  setwd(local_subdirname)
}

if (!best && !best_keep)
{  
  if (local)
  {
    if (Rtopicmodels_Gibbs)
    {
      save(Result,LLH0,LLH1,alpha_insert,delta,nb_topics,nb_real,nb_iter,llh_keep,Prop_OTU_removed,Prop_reads_removed,file=filename)
    } else if (Rtopicmodels_VEM)
    {
      save(Result,LLH0,LLH1,alpha_insert,nb_topics,nb_real,Prop_OTU_removed,Prop_reads_removed,file=filename)
    }
  }
  
  #########################
  # Llh convergence plots #
  #########################
  
  if (Rtopicmodels_Gibbs)
  {
    plotname = paste("llh_convergence_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,".pdf",sep="")
    pdf(plotname)
    #par(mar=c(5.1,4.1,4.1,2.1))
    par(mar=c(5.1,5.1,4.1,2.1))
    
    errbar(seq(llh_keep,nb_iter,llh_keep),LLH0[1,],LLH0[1,]+LLH0[2,],LLH0[1,]-LLH0[2,],ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
    
    legend(x="bottomright",legend=c(filename_insert,alpha_insert,paste("delta =",delta),paste("nb topics =",nb_topics),paste("nb realizations =",nb_real)),col="black",inset=0)
    title("Log-likelihood value vs number of iterations",cex.main=1.7)
    title(xlab="Number of iterations",ylab="Log-likelihood value",cex.lab=1.5)
    dev.off()
    
    plotname = paste("llh_convergence_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,"_allpoints.pdf",sep="")
    pdf(plotname)
    #par(mar=c(5.1,4.1,4.1,2.1))
    par(mar=c(5.1,5.1,4.1,2.1))
    
    plot(seq(llh_keep,nb_iter,llh_keep),LLH1[1,],ylim=range(LLH1),ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
    for (j in 2:nb_real)
      lines(seq(llh_keep,nb_iter,llh_keep),LLH1[j,],ylim=range(LLH1),lwd=2,type="p",col="black")
    
    legend(x="bottomright",legend=c(filename_insert,alpha_insert,paste("delta =",delta),paste("nb topics =",nb_topics),paste("nb realizations =",nb_real)),col="black",inset=0)
    title("Log-likelihood value vs number of iterations",cex.main=1.7)
    title(xlab="Number of iterations",ylab="Log-likelihood value",cex.lab=1.5)
    dev.off()
  }
  
} else if (best || best_keep)
{
  
  if (best && local) 
  {
    if (Rtopicmodels_Gibbs)
    {
      save(Result,alpha,delta,logbeta,documents,nb_doc,nb_terms,nb_topics,nb_real,nb_iter,llh_keep,llh,
           topic_compo_taxo,sort_normal_topic,KL_documents,documents,Prop_OTU_removed,Prop_reads_removed,file=filename)
    } else if (Rtopicmodels_VEM) {
      cat("Saving result ...\n")
      save(Result,alpha_est,logbeta,documents,nb_doc,nb_terms,nb_topics,nb_real,llh,
           topic_compo_taxo,sort_normal_topic,KL_documents,documents,Prop_OTU_removed,Prop_reads_removed,file=filename)
    }
  } else if (best_keep && (j_select == 1) && local)
  {
    if (Rtopicmodels_Gibbs)
    {
      save(Result,alpha,delta,logbeta,documents,nb_doc,nb_terms,nb_topics,nb_real,nb_iter,llh_keep,llh,AIC1,AIC2,LLH_final_real1,Akaike_weights_llh,llh_iterations,
           topic_compo_taxo,sort_normal_topic,KL_documents,documents,Prop_OTU_removed,Prop_reads_removed,file=filename)
    } else if (Rtopicmodels_VEM) {
      cat("Saving result ...\n")
      save(Result,alpha_est,logbeta,documents,nb_doc,nb_terms,nb_topics,nb_real,llh,AIC1,AIC2,LLH_final_real1,Akaike_weights_llh,llh_iterations,
           topic_compo_taxo,sort_normal_topic,KL_documents,documents,Prop_OTU_removed,Prop_reads_removed,file=filename)
    }
  }
  
  #############################################
  # Comparison with abiotic variables: output #
  #############################################
  
  if (best_keep && (j_select == 1) && ((data_pp && !testdata) || data_gs) && abiotic_variables)
  {
    abiotic_variables_dirname = paste0(local_subdirname,"Abiotic_variables")
    if (!(file.exists(abiotic_variables_dirname)))
      dir.create(abiotic_variables_dirname)
    setwd(abiotic_variables_dirname)
    
    #######################
    if (data_pp)
    {
      # saving the lidar and chemistery correlation values and associated p-values
      save(Correlation_lidar,Mean_cor_lidar_comparison_randomized,Sd_cor_lidar_comparison_randomized,p_value_lidar,file="Correlation_between_lidar_data_and_topics_for_best_real.Rdata")
      save(Correlation_chemi,Mean_cor_chemi_comparison_randomized,Sd_cor_chemi_comparison_randomized,p_value_chemi,file="Correlation_between_soil_chemistery_data_and_topics_for_best_real.Rdata")
      if (nb_topics == 3)
        save(SKL_bacteria,Mean_SKL_bacteria_comparison_randomized,Sd_SKL_bacteria_comparison_randomized,p_value_SKL_bacteria,file="SKL_samplewise_between_bacterial_topics_and_best_real.Rdata")
    } else if (data_gs)
      save(Correlation_climate,file="Correlation_between_climate_data_and_topics_for_best_real.Rdata")
    #######################
    
    ############
    # txt file #
    ############
    
    if (data_pp)
      case0_range = c("lidar","chemi")
    else if (data_gs)
      case0_range = "climate"
    
    for (case0 in case0_range)
    {
      if (case0 == "lidar")
      {
        Correlation_abiotic = Correlation_lidar
        Mean_cor_abiotic_comparison_randomized = Mean_cor_lidar_comparison_randomized
        Sd_cor_abiotic_comparison_randomized = Sd_cor_lidar_comparison_randomized
        p_value_abiotic = p_value_lidar
        ncol0 = ncol_lidar
        colnames_abiotic = colnames_lidar
      } else if (case0 == "chemi")
      {
        Correlation_abiotic = Correlation_chemi
        Mean_cor_abiotic_comparison_randomized = Mean_cor_chemi_comparison_randomized
        Sd_cor_abiotic_comparison_randomized = Sd_cor_chemi_comparison_randomized
        p_value_abiotic = p_value_chemi
        ncol0 = ncol_chemi
        colnames_abiotic = colnames_chemi
      } else if (case0 == "climate")
      {
        Correlation_climate = Correlation_climate
        #             Mean_cor_abiotic_comparison_randomized = Mean_cor_climate_comparison_randomized
        #             Var_cor_abiotic_comparison_randomized = Var_cor_climate_comparison_randomized
        #             p_value_abiotic = p_value_climate
        ncol0 = ncol_climate
        colnames_abiotic = colnames_climate
      }
      
      # Abiotic file does not report the comparison to PCA axes in the case pca_abiotic = 1
      Abiotic_file = paste0(case0,"_comparison.txt")
      if (!data_gs)
        abiotic_data.frame = as.data.frame(matrix(ncol=length(colnames_abiotic)*3,nrow=nb_topics,data=0))
      else
        abiotic_data.frame = as.data.frame(matrix(ncol=length(colnames_abiotic)*2,nrow=nb_topics,data=0))
      
      if (data_pp && (nb_topics == 3))
        SKL_bacteria_data.frame = as.data.frame(matrix(ncol=nb_topics*3,nrow=nb_topics,data=0))
      
      write(paste(barcode_insert,"-",nb_topics,"topics\n%%%%%%%%%%%%%\n"),file=Abiotic_file,append=F)
      for (k in 1:nb_topics)
      {
        k0 = rev(sort_normal_topic$ix)[k]
        if (nb_topics == 3 && data_pp)
        {
          sorted_correlations_to_bacteria = sort.int(Correlation_abiotic[[1]][k0,ncol0+(1:3)],index=T,decreasing=T)
          j0 = ncol0+sorted_correlations_to_bacteria$ix[1]
          write(paste0("\n%%%%%%%%%%%%%\nFor topic ",k," - ",colnames_abiotic[j0],
                       " assemblage (rho = ",Correlation_abiotic[[1]][k0,j0],", p = ",p_value_abiotic[k0,j0],"):\n%%%%%%%%%%%%%"),file=Abiotic_file,append=T)
          i = 0
          i_insert = ""
          while (any(rownames(abiotic_data.frame)==paste(assemblage_names_vect[k],barcode_insert,i_insert)))
          {
            i = i+1
            i_insert = i
          }
          rownames(abiotic_data.frame)[k] = paste(assemblage_names_vect[k],barcode_insert,i_insert)
        } else
        {
          write(paste0("\n%%%%%%%%%%%%%\nFor topic ",k,":\n%%%%%%%%%%%%%"),file=Abiotic_file,append=T)
          rownames(abiotic_data.frame)[k] = paste("Assemblage",k,barcode_insert)
        }
        for (j in 1:length(colnames_abiotic))
        {
          write(paste0("\n",colnames_abiotic[j],"\n%%%%%%%%%%%%%"),file=Abiotic_file,append=T)
          write(paste("Correlation =",Correlation_abiotic[[1]][k0,j]),file=Abiotic_file,append=T)
          # Comparison to randomizations is not implemented yet for GS data:
          if (!data_gs)
          {
            write(paste("Effect size with spatial randomizations =",Correlation_abiotic[[1]][k0,j] - Mean_cor_abiotic_comparison_randomized[k0,j]),file=Abiotic_file,append=T)
            write(paste("Standardized effect size with spatial randomizations =",
                        (Correlation_abiotic[[1]][k0,j] - Mean_cor_abiotic_comparison_randomized[k0,j])/Sd_cor_abiotic_comparison_randomized[k0,j]),file=Abiotic_file,append=T)
            write(paste("Normalized effect size with spatial randomizations =",
                        (Correlation_abiotic[[1]][k0,j] - Mean_cor_abiotic_comparison_randomized[k0,j])/(1-Mean_cor_abiotic_comparison_randomized[k0,j])),file=Abiotic_file,append=T)
            write(paste("p-value for spatial randomizations =",p_value_abiotic[k0,j]),file=Abiotic_file,append=T)
            
            abiotic_data.frame[k,(3*(j-1)+1):(3*(j-1)+3)] = c(Correlation_abiotic[[1]][k0,j],(Correlation_abiotic[[1]][k0,j] - Mean_cor_abiotic_comparison_randomized[k0,j])/(1-Mean_cor_abiotic_comparison_randomized[k0,j]),p_value_abiotic[k0,j])
            if (k==1)
              colnames(abiotic_data.frame)[(3*(j-1)+1):(3*(j-1)+3)] = c(paste(colnames_abiotic[j],"corr."),paste(colnames_abiotic[j],"norm. ES"),paste(colnames_abiotic[j],"p-value"))  
          } else
          {
            abiotic_data.frame[k,(2*(j-1)+1):(2*(j-1)+2)] = c(Correlation_abiotic[[1]][k0,j],p_value_abiotic[k0,j])
            if (k==1)
              colnames(abiotic_data.frame)[(2*(j-1)+1):(2*(j-1)+2)] = c(paste(colnames_abiotic[j],"corr."),paste(colnames_abiotic[j],"p-value"))
          }
          
          # write(paste("p-value for Student's t-test =",Correlation_abiotic[[2]][k0,j]),file=Abiotic_file,append=T)
        }
        
        if (data_pp && nb_topics == 3)
        {
          for (j in 1:nb_topics)
          {
            SKL_bacteria_data.frame[k,(3*(j-1)+1):(3*(j-1)+3)] = c(SKL_bacteria[[1]][k0,j],SKL_bacteria[[2]][k0,j],p_value_SKL_bacteria[k0,j])
            if (k==1)
              colnames(SKL_bacteria_data.frame)[(3*(j-1)+1):(3*(j-1)+3)] = c(paste(colnames_abiotic[ncol0+j],"SKL"),paste(colnames_abiotic[ncol0+j],"norm. ES"),paste(colnames_abiotic[ncol0+j],"p-value"))  
          }
          rownames(SKL_bacteria_data.frame) = rownames(abiotic_data.frame)
        }
      }
      
      # Comparison to PCA components is not included in the output
      saveRDS(abiotic_data.frame,file=paste0(case0,"_data.frame.rds"))
      if (data_pp && nb_topics == 3)
        saveRDS(SKL_bacteria_data.frame,file="SKL_bacteria_data.frame.rds")
      setwd(local_subdirname)
      write.csv(abiotic_data.frame,paste0(case0,"_comparison.csv"))
      if (data_pp && nb_topics == 3)
        write.csv(SKL_bacteria_data.frame,"bacteria_SKL_samplewise_comparison.csv")
      setwd(abiotic_variables_dirname)
      
      #########
      # plots #
      #########
      
      if (pca_abiotic)
        case1_range = c(1,3)
      else 
        case1_range = 1
      
      case2_range = c("_Bonferroni","")
      
      for (case1 in case1_range)
      {
        for (case2 in case2_range)
        {
          ######################
          if (case1 == 1)
            pdf(paste0("Correlation_between_",case0,"_data_and_topics_for_best_real",case2,".pdf"))
          else if (case1 == 3)
            pdf(paste0("Correlation_between_",case0,"_data_and_topics_for_best_real_pca",case2,".pdf"))
          par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
          #par(mar=c(5.1,4.1,4.1,2.1)
          #bottom left top right
          par(mar=c(15.1,10.1,4.1,4.1))
          for (k in 1:nb_topics)
          {            
            k0 = rev(sort_normal_topic$ix)[k]
            plot(Correlation_abiotic[[case1]][k0,],type="p",ann=F,yaxt="n",xaxt="n")
            for (j in 1:ncol0)
            {
              if (case2 == "")
              {
                if (Correlation_abiotic[[case1+1]][k0,j]>0.05)
                  lines(j,Correlation_abiotic[[case1]][k0,j],col="red",type="p")
              } else if (case2 == "_Bonferroni")
              {
                if (Correlation_abiotic[[case1+1]][k0,j]>0.05/ncol0)
                  lines(j,Correlation_abiotic[[case1]][k0,j],col="red",type="p")
              }
            }
            
            axis(2, ylim=range(Correlation_abiotic[[case1]][k0,]), col='black')
            if (ncol0 < 15)
            {
              axis(1, at=1:ncol0, labels = F)
              labels = vector(length=ncol0,mode="character")
              for (i in 1:ncol0)
              {
                if (case1 == 3)
                  labels[i] = paste(i,"th PCA axis",sep="")
                else if (case1 == 1)
                  labels[i] = colnames_abiotic[i]
              }
            } else 
            {
              if (case1 == 1)
              {
                axis(1, at=1:ncol0, labels = F)
                text(1:ncol0, par("usr")[3], srt = 90, pos=1, offset=1.5, labels = colnames_abiotic[1:ncol0], xpd = T, cex=0.6)
              }
              else if (case1 == 3)
              {
                axis(1, at=1:ncol0, labels = T)
                xlab("PCA axes")
              } 
            }
            
            text(1:ncol0, par("usr")[3], srt = 45, adj = c(1.05,1.75), labels = labels, xpd = TRUE, cex=1.1)
            title(ylab="Correlation value \nbetween spatial distributions")
            if (k==1)
              title(paste0("Correlation between ",case0," data and \nthe 1st most abundant assemblage"))
            else if (k==2)
              title(paste0("Correlation between ",case0," data and \nthe 2nd most abundant assemblage"))
            else if (k==3)
              title(paste0("Correlation between ",case0," data and \nthe 3rd most abundant assemblage"))
            else 
              title(paste0("Correlation between ",case0," data and the",k,"th most abundant assemblage",sep=""))
          }
          dev.off()
          ######################
        }
      }
      
    }
    setwd(local_subdirname)
  }
  
  ########################################################
  
  col_topic=c("black","blue","red","forestgreen","gold","mediumorchid2","chocolate4","cadetblue4","maroon","tan2","mediumturquoise","lightsalmon2","hotpink","greenyellow","lavender","khaki4")
  
  #########################################################
  # Creating the subdirectory specific to the realization #
  #########################################################
  
  if (select_real)
  {  
    subsubdirname = paste(local_subdirname,"ordered_realization_number_",Selected_real[j_select],"/",sep="")
    if (!(file.exists(subsubdirname)))
      dir.create(subsubdirname)
    setwd(subsubdirname)
    
    # Saving the variables specific to the realization in the subfolder, except for the best realization as it is already saved in the main folder
    if (j_select>1)
    {
      directory_file = "Directory.txt"
      write(paste(subsubdirname,filename,sep=""),directory_file)
      if (Rtopicmodels_Gibbs)
      {
        save(alpha,delta,logbeta,documents,nb_doc,nb_terms,nb_topics,nb_real,nb_iter,llh_keep,llh,AIC1,AIC2,LLH_final_real1,Akaike_weights_llh,llh_iterations,
             topic_compo_taxo,sort_normal_topic,KL_documents,documents,Ordered_realizations,file=filename)
      } else if (Rtopicmodels_VEM) 
      {
        save(alpha_est,logbeta,documents,nb_doc,nb_terms,nb_topics,nb_real,llh,AIC1,AIC2,LLH_final_real1,Akaike_weights_llh,llh_iterations,
             topic_compo_taxo,sort_normal_topic,KL_documents,documents,Ordered_realizations,file=filename)
      }
    }
    
    if (Rtopicmodels_VEM )
      alpha_est_mean = alpha_est/length_selected_real + alpha_est_mean
  }
  
  if (Rtopicmodels_VEM )
  {  
    alpha_file = "estimated_alpha.txt"
    write(alpha_est,alpha_file)
  }
  
  #############################
  # loglikelihood convergence #
  #############################
  
  if ((Rtopicmodels_VEM || (Rtopicmodels_Gibbs && llh_keep)))
  {
    plotname = paste("Loglikelihood_convergence_",alpha_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,".pdf",sep="")
    pdf(plotname)
    #par(mar=c(5.1,4.1,4.1,2.1))
    par(mar=c(5.1,5.1,4.1,2.1))
    
    plot(llh_iterations,ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
    
    if (Rtopicmodels_Gibbs)
      legend(x="bottomright",legend=c(filename_insert,alpha_insert,paste("delta =",delta),paste("nb topics =",nb_topics),paste("nb realizations =",nb_real)),col="black",inset=0)
    else if (Rtopicmodels_VEM)
      legend(x="bottomright",legend=c(filename_insert,alpha_insert,paste("nb topics =",nb_topics),paste("nb realizations =",nb_real)),col="black",inset=0)
    title("Loglikelihood convergence",cex.main=1.7)
    title(xlab="Iterations",ylab="Loglikelihood value",cex.lab=1.5)
    dev.off()
    
    ##############
    if (length(llh_iterations)>10)
    {
      plotname = paste("Loglikelihood_convergence_zoom_",alpha_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,".pdf",sep="")
      pdf(plotname)
      #par(mar=c(5.1,4.1,4.1,2.1))
      par(mar=c(5.1,5.1,4.1,2.1))
      
      plot(seq(11,length(llh_iterations),1),llh_iterations[-(1:10)],ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
      
      if (Rtopicmodels_Gibbs)
        legend(x="bottomright",legend=c(filename_insert,alpha_insert,paste("delta =",delta),paste("nb topics =",nb_topics),paste("nb realizations =",nb_real)),col="black",inset=0)
      else if (Rtopicmodels_VEM)
        legend(x="bottomright",legend=c(filename_insert,alpha_insert,paste("nb topics =",nb_topics),paste("nb realizations =",nb_real)),col="black",inset=0)
      title("Loglikelihood convergence",cex.main=1.7)
      title(xlab="Iterations",ylab="Loglikelihood value",cex.lab=1.5)
      dev.off()
    }
  }
  
  ###################################
  #                                 #
  ###############################   #
  # Topic abundance information #   #
  ###############################   #
  #                                 #
  ###################################
  
{
  
  if (select_real)
  {
    subsubsubdirname = paste(subsubdirname,"topics_abundance_info/",sep="")
    if (!(file.exists(subsubsubdirname)))
    {dir.create(subsubsubdirname)}
    setwd(subsubsubdirname) 
  } else if (!select_real)
  {
    subsubdirname = paste(local_subdirname,"topics_abundance_info/",sep="")
    if (!(file.exists(subsubdirname)))
    {dir.create(subsubdirname)}
    setwd(subsubdirname)  
  }
  
  # pdf("Global_topic_proportions.pdf")
  # plot(1:nb_topics,prop_topic,ann=FALSE)
  # title(xlab="Topic index",ylab="Global topic proportion",main="Global proportion of each topic")
  # dev.off()  
  
  if (topic_read_proportion_comput)
  {
    pdf("Sorted_global_topic_proportions.pdf")
    plot(nb_topics:1,sort_prop_topic$x,ann=FALSE,ylim=c(0,max(prop_topic)))
    lines(c(1,nb_topics),c(0,0),type="l",lty=2)
    title(ylab="Global topic proportion",main="Global proportion of each topic")
    dev.off()  
  }
  
  ########## to finish ###############
  pdf("Norm_global_topic_proportions.pdf")
  #par(mar=c(5.1,4.1,4.1,2.1))
  par(mar=c(5.1,5.1,4.1,2.1)) 
  par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
  plot(nb_topics:1,sort_normal_topic$x,ann=FALSE,ylim=c(0,max(normal_topic)))
  lines(c(1,nb_topics),c(0,0),type="l",lty=2)
  title(ylab="Relative abundance averaged over samples",main="Relative abundance of classes\n averaged over samples")
  mtext(side=1,"Classes ranked by\n relative abundances averaged over samples",line=4,cex=1.5)
  dev.off()  
  
  if (topic_read_proportion_comput)
  {
    pdf("Sorted_global_topic_log_proportions.pdf")
    plot(nb_topics:1,log(sort_prop_topic$x),ann=FALSE)
    #lines(c(1,nb_topics),c(1,1),type="l",lty=2)
    title(ylab="Global topic proportion (log)",main="Global proportion of each topic (log)")
    dev.off()  
  }
  
}

######################################################
#                                                    #
##################################################   #
# Information about topics' sequence composition #   #
##################################################   #
#                                                    #
######################################################

{
  
  if (select_real)
  {
    subsubsubdirname = paste(subsubdirname,"topics_sequence_composition_info/",sep="")
    if (!(file.exists(subsubsubdirname)))
      dir.create(subsubsubdirname)
    setwd(subsubsubdirname)
  } else if (!select_real)
  {
    subsubdirname = paste(local_subdirname,"topics_sequence_composition_info/",sep="")
    if (!(file.exists(subsubdirname)))
      dir.create(subsubdirname)
    setwd(subsubdirname)
  }
  
  if (dominance_calculation)
  {
    
    nb_bins = 10
    max_topic = vector(length=nb_topics,mode="numeric")
    for (j in 1:nb_topics)
      max_topic[j] = max(exp(logbeta)[j,]) 
    
    pdf("topic_dominance_by_sequence_hist.pdf")
    hist(max_topic, breaks = nb_bins, freq=T, xlab = "Topic's proportion taken by the most abundant sequence in topic",
         ylab = "Number of topics", main = "Number of topics\n with respect to the topic's proportion\n taken by the most abundant sequence in topic", xaxp = c(0, 1, 10))
    #ylab = "Number of sites", main = "Histogramme of the number of sites\n with respect to the proportion\n of the most abundant sequence")
    dev.off()
    
    pdf("topic_dominance_by_sequence_hist_percent.pdf")
    par(mar=c(5.1,5.1,5.1,2.1)) 
    par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
    h = hist(max_topic, breaks = nb_bins, plot=F)
    h$density = h$counts/sum(h$counts)
    plot(h, freq=F, ann=F, xaxp = c(0, 1, 10))
    title(main = "Proportion of the total number of classes\n with respect to how much dominant\n the most abundant sequence in the class is",
          ylab = "Proportion of the total number of classes")
    #ylab = "Number of sites", main = "Proportion of sites\n with respect to the proportion\n of the most abundant sequence")
    mtext("Proportion of the most abundant sequence in the class",side=1,line=4,cex=1.5) 
    dev.off()
    
    pdf("topic_dominance_by_sequence_vs_site-normalized_topic_abundance.pdf")
    par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
    par(mar=c(5.1,5.1,4.1,2.1))
    plot(max_topic[rev(sort_normal_topic$ix)],
         main="Proportion of the dominant sequence\n in each class",
         xlab="",
         ylab="Proportion of the dominant sequence in class")
    mtext("Classes ranked by\n relative abundances averaged over samples",side=1,line=4,cex=1.5)
    dev.off()
    
    nb_bins = 10
    prop_reads_seq_main_topic = vector(length=nb_terms,mode="numeric")  
    topicNormalized_prop_reads_seq_main_topic = vector(length=nb_terms,mode="numeric") 
    prop_max_topic_seq = vector(length=nb_terms,mode="numeric") 
    index_max_topic_seq = vector(length=nb_terms,mode="numeric") 
    # equality = vector(length=nb_terms,mode="numeric") 
    exp_logbeta=exp(logbeta)
    for (j in 1:nb_terms)
    { 
      count_total_seq = 0
      for (i in 1:nb_topics)
      {
        count_total_seq = exp_logbeta[i,j]*prop_topic[i]*tot_reads + count_total_seq
      }
      prop_max_topic_seq[j] = max(exp_logbeta[,j])
      # prop_max_topic_seq = sort.int(exp_logbeta[,j],decreasing=T,index.return=T)
      index_max_topic_seq = which(exp_logbeta[,j]==prop_max_topic_seq[j])
      # since 4 sequences are totally absent due to the removing of replicates and controles, they are hence present in the same default proportion in all 4 topics,
      # and index_max_topic_seq yields the indices of all 40 topics instead of one topic.
      # To reduce the dimension of index_max_topic_seq and prevent printing an error in the prompt, the first topic is arbitrarily chosen.
      if (length(index_max_topic_seq) > 1)
        index_max_topic_seq = index_max_topic_seq[1]
      # equality[j] = length(index_max_topic_seq) 
      prop_reads_seq_main_topic[j] = prop_max_topic_seq[j]*prop_topic[index_max_topic_seq]*tot_reads/count_total_seq
      topicNormalized_prop_reads_seq_main_topic[j] = prop_max_topic_seq[j]/sum(exp_logbeta[,j])
      # prop_reads_seq_main_topic[j] = sorted_seq$x[1]*prop_topic[sorted_seq$ix[1]]*tot_reads/count_total_seq
    }
    
    pdf("topic_exclusivity_hist.pdf")
    hist(prop_reads_seq_main_topic, breaks = nb_bins, freq=T, xlab = "Share of the sequence's total number of reads\n in the topic where the sequence is most dominant",
         ylab = "Number of sequences", main = "Number of sequences\n with respect to the share of their total number of reads\n in the topic where they are most dominant", xaxp = c(0, 1, 10))
    dev.off()
    
    pdf("topic_exclusivity_hist_percent.pdf")
    par(mar=c(5.1,5.1,7.1,2.1)) 
    par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
    h = hist(prop_reads_seq_main_topic, breaks = nb_bins, plot=F)
    h$density = h$counts/sum(h$counts)
    plot(h, freq=F, ann=F, xaxp = c(0, 1, 10))
    title(main = "Proportion of the number of sequences\n with respect to the share of their 
          total number of reads\n in the class where they are most dominant",
          ylab = "Proportion of the total number of sequences")
    #ylab = "Number of sites", main = "Proportion of sites\n with respect to the proportion\n of the most abundant sequence")
    mtext("Share of the sequence's total number of reads\n in the class where the sequence is most dominant",side=1,line=4,cex=1.5) 
    dev.off()
    
    pdf("topic_distinctiveness.pdf")
    par(mar=c(5.1,5.1,7.1,2.1)) 
    par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
    h = hist(topicNormalized_prop_reads_seq_main_topic, breaks = nb_bins, plot=F)
    h$density = h$counts/sum(h$counts)
    plot(h, freq=F, ann=F, xaxp = c(0, 1, 10))
    title(main = "Community distinctiveness",
          ylab = "Proportion of MOTUs")
    #ylab = "Number of sites", main = "Proportion of sites\n with respect to the proportion\n of the most abundant sequence")
    mtext("Community-normalized proportion of sequences\n in the community where the MOTU is most dominant",side=1,line=4,cex=1.5) 
    dev.off()
    
    pdf("Topic_dominance_by_sequence_wrt_site-normalized_sequence_rank.pdf")
    par(mar=c(5.1,5.1,4.1,4.1))
    plot(prop_max_topic_seq[normal_ordered_seq],
         main="Proportion of a sequence\n in the topic where it is most dominant \n- sequences ranked by site-normalized abundances",
         xlab="Sequences ranked by site-normalized abundances",
         ylab="Proportion of a sequence\n in the topic where it is most dominant")
    dev.off()
    
    pdf("Topic_dominance_by_sequence_wrt_site-normalized_sequence_rank_250firstseq.pdf")
    par(mar=c(5.1,5.1,4.1,4.1))
    par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
    plot(prop_max_topic_seq[normal_ordered_seq[1:250]],
         main="Proportion of a sequence\n in the topic where it is most dominant\n- 250 first sequences",
         xlab="Sequences ranked by relative abundances averaged over samples",
         ylab="Proportion of a sequence\n in the class where it is most dominant")
    dev.off()
    
    pdf("Total_read_share_in_main_topic_vs_site-normalized_seq_ordening.pdf")
    par(mar=c(5.1,5.1,4.1,4.1))
    par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
    plot(prop_reads_seq_main_topic[normal_ordered_seq],
         main="Share of a sequence's total number of reads\n in the class where it is most dominant",
         ylab = "Share of the sequence's total number of reads",
         xlab = "")
    mtext("Sequences ranked by\n relative abundances averaged over samples",side=1,line=4,cex=1.5)
    dev.off()
    
    pdf("Total_read_share_in_main_topic_vs_site-normalized_seq_ordening_40firstseq.pdf")
    par(mar=c(5.1,5.1,5.1,4.1))
    par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
    plot(prop_reads_seq_main_topic[normal_ordered_seq[1:40]],
         main="Share of a sequence's total number of reads\n in the class where it is most dominant\n- 40 first sequences",
         ylab = "Share of the sequence's total number of reads",
         xlab = "")
    mtext("Sequences ranked by\n relative abundances averaged over samples",side=1,line=4,cex=1.5)
    dev.off()
    
    # end dominance_calculation
  }
  
  pdf("Site-normalized_MOTU_composition_for_each_topic_20firstseq_assignPrecision.pdf")
  par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
  #par(mar=c(5.1,4.1,4.1,2.1)
  #bottom left top right
  par(mar=c(15.1,10.1,4.1,4.1))
  for (k in 1:nb_topics)
  {
    plot(topic_compo_taxo[[rev(sort_normal_topic$ix)[k]]][1:20,1],type="p",ann=F,yaxt="n",xaxt="n")
    axis(2, ylim=range(topic_compo_taxo[[rev(sort_normal_topic$ix)[k]]][1:20,1]), col='black')
    axis(1, at=1:20, labels = F)
    labels = vector(length=20,mode="character")
    for (i in 1:20)
      labels[i]=paste(topic_compo_taxo[[rev(sort_normal_topic$ix)[k]]][i,11],topic_compo_taxo[[rev(sort_normal_topic$ix)[k]]][i,10],topic_compo_taxo[[rev(sort_normal_topic$ix)[k]]][i,2],sep=" - ")
    text(1:20, par("usr")[3], srt = 45, adj = c(1.05,1.75), labels = labels, xpd = TRUE, cex=1.1)
    if (k==1)
      title(ylab="Relative abundance in the\n 1st most abundant class")
    else if (k==2)
      title(ylab="Relative abundance in the\n 2nd most abundant class")
    else 
      title(ylab=paste("Relative abundance in the\n",k,"th most abundant class",sep=""))
    title("Relative class composition\n for the 20 most abundant sequences")
  }
  dev.off()
  
  pdf("Site-normalized_MOTU_composition_for_each_topic_20firstseq_assignPrecision_noNumber.pdf")
  par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
  #par(mar=c(5.1,4.1,4.1,2.1)
  #bottom left top right
  par(mar=c(15.1,10.1,4.1,4.1))
  for (k in 1:nb_topics)
  {
    plot(topic_compo_taxo[[rev(sort_normal_topic$ix)[k]]][1:20,1],type="p",ann=F,yaxt="n",xaxt="n")
    axis(2, ylim=range(topic_compo_taxo[[rev(sort_normal_topic$ix)[k]]][1:20,1]), col='black')
    axis(1, at=1:20, labels = F)
    labels = vector(length=20,mode="character")
    for (i in 1:20)
      labels[i]=paste(topic_compo_taxo[[rev(sort_normal_topic$ix)[k]]][i,11],topic_compo_taxo[[rev(sort_normal_topic$ix)[k]]][i,10],sep=" - ")
    text(1:20, par("usr")[3], srt = 45, adj = c(1.05,1.75), labels = labels, xpd = TRUE, cex=1.1)
    if (k==1)
      title(ylab="Relative abundance in the\n 1st most abundant class")
    else if (k==2)
      title(ylab="Relative abundance in the\n 2nd most abundant class")
    else 
      title(ylab=paste("Relative abundance in the\n",k,"th most abundant class",sep=""))
    title("Relative class composition\n for the 20 most abundant sequences")
  }
  dev.off()
  
  # Writing the taxonomic composition to a .txt file
  Taxo_compo_file = "Taxonomic_composition.txt"
  Ordered_topic_compo_taxo_full = list()
  Ordered_topic_compo_taxo = list()
  write(paste(barcode_insert,"-",nb_topics,"topics\n%%%%%%%%%%%%%\n"),file=Taxo_compo_file,append=F)
  for (k in 1:nb_topics)
  {
    k0 = rev(sort_normal_topic$ix)[k]
    Ordered_topic_compo_taxo[[k]] = topic_compo_taxo[[k0]][which(topic_compo_taxo[[k0]][,1]>1/sum(data2m)),]
    Ordered_topic_compo_taxo_full[[k]] = topic_compo_taxo[[k0]]
    if (nb_topics == 3 && data_pp && abiotic_variables)
    {
      sorted_correlations_to_bacteria = sort.int(Correlation_abiotic[[1]][k0,ncol0+(1:3)],index=T,decreasing=T)
      j0 = ncol0+sorted_correlations_to_bacteria$ix[1]
      write(paste0("\n%%%%%%%%%%%%%\nFor topic ",k," - ",colnames_abiotic[j0],
                   " assemblage (rho = ",Correlation_abiotic[[1]][k0,j0],", p = ",p_value_abiotic[k0,j0],"):\n%%%%%%%%%%%%%"),file=Taxo_compo_file,append=T)
    } else
      write(paste0("\n%%%%%%%%%%%%%\nFor topic ",k,":\n%%%%%%%%%%%%%"),file=Taxo_compo_file,append=T)
    for (i in 1:min(50,length(which(topic_compo_taxo[[k0]][,1]>1/sum(data2m)))))
      write(paste(topic_compo_taxo[[k0]][i,1],"-",topic_compo_taxo[[k0]][i,11],"-",topic_compo_taxo[[k0]][i,10]),file=Taxo_compo_file,append=T)
  }
  
  save(Ordered_topic_compo_taxo,Ordered_topic_compo_taxo_full,file = "Taxonomic_composition.Rdata")
}

##################################################
#                                                #
##############################################   #
# Information about topics' site repartition #   #
##############################################   #
#                                                #
##################################################

{
  
  if (select_real)
  {
    subsubsubdirname = paste(subsubdirname,"/topics_site_repartition_info/",sep="")
    if (!(file.exists(subsubsubdirname)))
      dir.create(subsubsubdirname)
    setwd(subsubsubdirname) 
  } else if (!select_real)
  {
    subsubdirname = paste(local_subdirname,"/topics_site_repartition_info/",sep="")
    if (!(file.exists(subsubdirname)))
      dir.create(subsubdirname)
    setwd(subsubdirname) 
  }
  
  if (dominance_calculation)
  {
    
    nb_bins = 10
    topic_dominance_by_site = vector(length = nb_doc, mode = "numeric")
    for (i in 1:nb_doc)
      topic_dominance_by_site[i] = max(documents[i,])
    
    pdf("site_dominance_by_topic_hist_10bins.pdf")
    hist(topic_dominance_by_site, breaks = nb_bins, freq=T, xlab = "Site's proportion taken\n by the most abundant topic in site",
         ylab = "Number of sites", main = "Number of sites\n with respect to the site's proportion taken\n by the most abundant topic in site", xaxp = c(0, 1, 10))
    dev.off()
    
    pdf("site_dominance_by_topic_hist_percent_10bins.pdf")
    h = hist(topic_dominance_by_site, breaks = nb_bins, plot=F)
    h$density = h$counts/sum(h$counts)
    plot(h, freq=F, xlab = "Site's proportion taken\n by the most abundant topic in site",
         ylab = "Proportion of the total number of sites", main = "Proportion of the total number of sites\n with respect to the site's proportion taken\n by the most abundant topic in site", xaxp = c(0, 1, 10))
    dev.off()
    
    nb_bins = 20
    pdf("site_dominance_by_topic_hist_20bins.pdf")
    hist(topic_dominance_by_site, breaks = nb_bins, freq=T, xlab = "Site's proportion taken\n by the most abundant topic in site",
         ylab = "Number of sites", main = "Number of sites\n with respect to the site's proportion taken\n by the most abundant topic in site", xaxp = c(0, 1, 20))
    dev.off()
    
    # pdf("site_dominance_by_topic_hist_percent_20bins.pdf")
    # h = hist(topic_dominance_by_site, breaks = nb_bins, plot=F)
    # h$density = h$counts/sum(h$counts)
    # plot(h, freq=F, xlab = "Sample's proportion taken\n by the most abundant class in sample",
    #      ylab = "Proportion of the total number of samples", main = "Proportion of the total number of samples\n with respect to the sample's proportion taken\n by the most abundant class in sample", xaxp = c(0, 1, 20))
    # dev.off()
    
    pdf("site_dominance_by_topic_hist_percent_20bins.pdf")
    par(mar=c(5.1,5.1,5.1,2.1)) 
    par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
    h = hist(topic_dominance_by_site, breaks = nb_bins, plot=F)
    h$density = h$counts/sum(h$counts)
    plot(h, freq=F, ann=F, xaxp = c(0, 1, 10))
    title(main = "Proportion of the total number of samples\n with respect to the sample's proportion taken\n by the most abundant class in sample",
          ylab = "Proportion of the total number of samples")
    #ylab = "Number of sites", main = "Proportion of sites\n with respect to the proportion\n of the most abundant sequence")
    mtext("Sample's proportion taken\n by the most abundant class in sample",side=1,line=4,cex=1.5) 
    dev.off()
    
    if (topic_read_proportion_comput)
    {
      pdf("Number_of_dominant_sites_for_each_topic.pdf")
      plot(nb_topics:1,topic_site_nb[sort_prop_topic$ix,1],ann=F,type="p",ylim=range(topic_site_nb))
      for (k in 2:5)
        lines(nb_topics:1,topic_site_nb[sort_prop_topic$ix,k],type="p",col=col_topic[k])
      title(ylab="Number of dominant sites",xlab="Topics ordered according to their global read proportion",main="Number of dominant sites for each topic")
      legend(x="topleft",legend=c("Threshold 0.5","Threshold 0.6","Threshold 0.7","Threshold 0.8","Threshold 0.9"),text.col=col_topic[1:5],inset=0.1)
      dev.off()
    }
    
    # end of abundance_calculation
  }
  
  if (testdata)
  {
    pdf("Topic_ordered_by_site-normalized_abundance_sample_distribution.pdf")
    #par(mfrow=c(2,2))
    for (k in 1:nb_topics)
    {
      par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
      plot(documents[,rev(sort_normal_topic$ix)[k]],type="l",
           main=paste("Community #",rev(sort_normal_topic$ix)[k],sep=""),
           xlab="Samples",
           ylab="Proportion of the community in the sample")
    }
    dev.off()
    
    pdf("Topic_ordered_by_site-normalized_abundance_sample_distribution_oneplot.pdf")
    for (k in 1:nb_topics)
    {
      par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
      if (k==1)
      {
        plot(documents[,rev(sort_normal_topic$ix)[k]],type="l",col=rainbow_hcl(nb_topics)[k],
             main="Composition of samples",
             xlab="Samples",
             ylab="Proportion of the community in the sample")
      } else
        lines(documents[,rev(sort_normal_topic$ix)[k]],type="l",col=rainbow_hcl(nb_topics)[k])
    }
    for (k in 1:true_nb_topics)
      lines(true_documents[,k],type="l",col="black",lty=2)
    dev.off()
  }
  
  #########################################################################
  # Topic composition maps                                                #
  #########################################################################
  
  if ((data_pp && !testdata) || data_h20 || data_gs)
  {
    source(paste0(code_insert,"LDA_spatial_maps_fun.R"))
    
    LDA_spatial_maps_result = LDA_spatial_maps_fun(j_select,kriged_real,kriging,nb_topics,data_pp,data_h20,
                                                   data_gs,filled,filled_with_gaps,ifelse(exists("Missing_positions_indices"),Missing_positions_indices,NA),
                                                   documents,ifelse(exists("coordGS"),coordGS,NA),local_prefix,data_insert,sort_normal_topic,subsubsubdirname)
    
    lidar.plot = LDA_spatial_maps_result$lidar.plot
    tmp.plot = LDA_spatial_maps_result$tmp.plot
    tmp.one.plot = LDA_spatial_maps_result$tmp.one.plot
    tmp.one.plot.discrete = LDA_spatial_maps_result$tmp.one.plot.discrete
    spatial_topicmix = LDA_spatial_maps_result$spatial_topicmix
    # spatial_topicmix_kriged is NA if it had already been saved previously
    spatial_topicmix_kriged = LDA_spatial_maps_result$spatial_topicmix_kriged
    
    if (kriging && any(kriged_real == j_select))
    {
      if (!(file.exists("Spatial_topicmix_kriged.rds")))
        saveRDS(spatial_topicmix_kriged,file="Spatial_topicmix_kriged.rds")
      
      if (data_pp && (nb_topics == 3))
      {
        #ggsave(filename = "Topic_ordered_by_site-normalized_abundance_composition_maps_kriged.pdf", do.call("arrangeGrob", c(tmp.plot, ncol=nb_topics)), height = 4, width = 10)
        ggsave(filename = "Topic_ordered_by_site-normalized_abundance_composition_maps_kriged.pdf", do.call("arrangeGrob", c(tmp.plot, lidar.plot, ncol=nb_topics)), height = 10/3*4/3*2, width = 10)
        ggsave(filename = "Topic_ordered_by_site-normalized_abundance_composition_maps_kriged_oneplot.pdf",  do.call("arrangeGrob", list(tmp.one.plot, tmp.one.plot.discrete, ncol=2)), height = 10/2*4/3, width = 10)
        #ggsave(filename = "Topic_ordered_by_site-normalized_abundance_composition_maps_kriged.pdf", do.call("arrangeGrob", c(tmp.plot, lidar.plot, ncol=nb_topics)), width = 10)
        #ggsave(filename = "Test1.pdf", tmp.one.plot, width=10)
      } else
      {
        ggsave(filename = "Topic_ordered_by_site-normalized_abundance_composition_maps_kriged.pdf", do.call("arrangeGrob", c(tmp.plot, ncol=min(nb_topics,4))), height = 10/min(nb_topics,4)*((nb_topics-1)%/%4 + 1)*4/3, width = 10)
        #ggsave(filename = "Topic_ordered_by_site-normalized_abundance_composition_maps_kriged.pdf", do.call("arrangeGrob", c(tmp.plot, ncol=min(nb_topics,4))), width = 10)
        ggsave(filename = "Topic_ordered_by_site-normalized_abundance_composition_maps_kriged_oneplot.pdf",  do.call("arrangeGrob", list(tmp.one.plot, tmp.one.plot.discrete, ncol=2)), height = 10/2*4/3, width = 10)
      }
      
      pdf("Topic_ordered_by_site-normalized_abundance_composition_maps_kriged_onebyone.pdf")
      for (k in 1:nb_topics)
      {
        #grid.newpage()
        #do.call("arrangeGrob",tmp.plot[[k]])
        print(tmp.plot[[k]])
      }
      dev.off()
    }
    
    if (topic_read_proportion_comput && !data_gs)
    {
      ############################
      pdf("Topic_ordered_by_total_abundance_composition_maps.pdf")
      # 2 topics :
      #par(mfrow=c(2,1))
      # 5 topics :
      #par(mfrow=c(3,2))
      # 10 topics
      par(mfrow=c(2,2))
      # par(mar = c(bottom, left, top, right))
      par(mar = c(5, 4, 4, 5) + 0.1)
      #lattice::levelplot(abund2.pred~x+y, z2, 
      #col.regions=topo.colors(100), aspect = "iso",contour=T,main=binomial)
      # Loop over topics (one map per topic, the color stands for the proportion of the topic)
      for (k in 1:nb_topics)
      {  
        k0 = rev(sort_prop_topic$ix)[k]
        # colors Red Green Blue
        # c(226/255,2/255),c(226/255,63/255),c(226/255,165/255)
        color2D.matplot(spatial_topicmix[[k0]],c(226/255,2/255),c(226/255,63/255),c(226/255,165/255),
                        extremes=NA,cellcolors=NA,show.legend=F,nslices=10,xlab="",xrange=c(0,1),
                        ylab="",do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
                        border=NA,na.color="white",main=paste("Topic ",k," - ",k0,sep=""))
        axis(2, at=c(5,10,15,20,25), labels=c("50","100","150","200","250"))
        axis(1, at=c(5,10,15,20,25,30,35), labels=c("50","100","150","200","250","300","350"))
        color.legend(xl=39+1.392857,yb=4.53125,xr=39+2.785714,yt=13.59375,legend=c(0,1),
                     rect.col=color.scale(seq(0,1,1/255),c(226/255,2/255),c(226/255,63/255),c(226/255,165/255)),gradient="y",cex=1.5, 
                     align="rb")
      }
      dev.off()
    }
    
    if ((!kriging || !any(kriged_real == j_select)) && !data_gs)
    {
      pdf("Topic_ordered_by_site-normalized_abundance_composition_maps.pdf")
      # 2 topics :
      #par(mfrow=c(2,1))
      # 5 topics :
      #par(mfrow=c(3,2))
      # 10 topics
      par(mfrow=c(2,2))
      par(cex.lab=1.5,cex.main=1.7,cex.axis=2)
      # par(mar = c(bottom, left, top, right))
      #par(mar = c(5, 4, 4, 5) + 0.1)
      #lattice::levelplot(abund2.pred~x+y, z2, 
      #col.regions=topo.colors(100), aspect = "iso",contour=T,main=binomial)
      for (k in 1:nb_topics)
      {
        k0 = rev(sort_normal_topic$ix)[k]
        # colors Red Green Blue
        # rev(col_values$red)/255,rev(col_values$green)/255,rev(col_values$blue)/255
        # c(226/255,2/255),c(226/255,63/255),c(226/255,165/255)
        color2D.matplot(spatial_topicmix[[k0]],c(226/255,2/255),c(226/255,63/255),c(226/255,165/255),
                        extremes=NA,cellcolors=NA,show.legend=F,nslices=10,xlab="",xrange=c(0,1),
                        ylab="",do.hex=FALSE,axes=FALSE,show.values=FALSE,vcol="white",vcex=1,
                        border=NA,na.color="white",main=paste("Topic ",k0," - #",k,sep=""))
        axis(2, at=c(0,5,10,15,20,25,30), labels=c("","50","100","150","200","250","300"))
        axis(1, at=c(0,5,10,15,20,25,30,35,40), labels=c("0","","100","","200","","300","","400"))
        #       color.legend(xl=39+1.392857,yb=4.53125,xr=39+2.785714,yt=13.59375,legend=c(0,1),
        #                    rect.col=color.scale(seq(0,1,1/255),c(226/255,2/255),c(226/255,63/255),c(226/255,165/255)),gradient="y",cex=1.5, 
        #                    align="rb")
      }
      dev.off()
    }
    
    # Plots the spatial distribution of topics while removing abundance values lower than threshold_val   
    #     threshold_val = 0.1
    #     if ((!kriging || !any(kriged_real == j_select)) && !data_gs)
    #     {
    #       if (filled && !filled_with_gaps)
    #       {
    #         pdf(paste("Topic_ordered_by_site-normalized_abundance_composition_maps_onebyone_thres",threshold_val,".pdf",sep=""))
    #         # 2 topics :
    #         #par(mfrow=c(2,1))
    #         # 5 topics :
    #         #par(mfrow=c(3,2))
    #         # 10 topics
    #         #par(mfrow=c(2,2))
    #         par(cex.lab=1.5,cex.main=1.7,cex.axis=2)
    #         # par(mar = c(bottom, left, top, right))
    #         par(mar = c(5, 4, 4, 5) + 0.1)
    #         #lattice::levelplot(abund2.pred~x+y, z2, 
    #         #col.regions=topo.colors(100), aspect = "iso",contour=T,main=binomial)
    #         # Loop over topics (one map per topic, the color stands for the proportion of the topic)
    #         for (k in 1:nb_topics)
    #         {  
    #           k0 = rev(sort_normal_topic$ix)[k]
    #           # colors based on colorspace::divergence_hcl(3)
    #           # c(226/255,2/255),c(226/255,63/255),c(226/255,165/255)
    #           # rev(col_values$red)/255,rev(col_values$green)/255,rev(col_values$blue)/255
    #           # grey to green : c(242,0)/255,c(242,166)/255,c(242,0)/255
    #           # grey to blue : c(226/255,2/255),c(226/255,63/255),c(226/255,165/255)
    #           cellcol = matrix(nrow=29,ncol=39,data=0.1)
    #           # Uniform grey between 0 and threshold_val, and then a grey to blue color gradient
    #           cellcol[spatial_topicmix[[k0]]<threshold_val] = color.scale(spatial_topicmix[[k0]][spatial_topicmix[[k0]]<threshold_val],226/255,226/255,226/255)
    #           cellcol[spatial_topicmix[[k0]]>threshold_val] = color.scale(spatial_topicmix[[k0]][spatial_topicmix[[k0]]>threshold_val],c(226/255,2/255),c(226/255,63/255),c(226/255,165/255))
    #           color2D.matplot(spatial_topicmix[[k0]],
    #                           extremes=NA,cellcolors=cellcol,show.legend=F,nslices=10,xlab="",
    #                           ylab="",do.hex=FALSE,axes=F,show.values=FALSE,vcol="white",vcex=1,xrange=c(0,1),
    #                           border=NA,na.color="white",main=paste("Topic ",k0," - #",k,sep=""))
    #           axis(2, at=c(0,5,10,15,20,25,30), labels=c("","50","100","150","200","250","300"))
    #           axis(1, at=c(0,5,10,15,20,25,30,35,40), labels=c("0","","100","","200","","300","","400"))
    #           #       color.legend(xl=39+1.392857,yb=4.53125,xr=39+2.785714,yt=13.59375,legend=c(0,1),
    #           #                    rect.col=color.scale(seq(0,1,1/255),c(226/255,2/255),c(226/255,63/255),c(226/255,165/255)),gradient="y",cex=1.5, 
    #           #                    align="rb")
    #           color.legend(xl=39+1.392857,yb=4.53125,xr=39+2.785714,yt=13.59375,legend=c(0,1),
    #                        rect.col=color.scale(c(rep(threshold_val,threshold_val/(1-threshold_val)*255),
    #                                               seq(threshold_val,1,1/(255+threshold_val/(1-threshold_val)*255))),
    #                                             c(226/255,2/255),c(226/255,63/255),c(226/255,165/255)),gradient="y",cex=1.5, 
    #                        align="rb")
    #           
    #         }
    #         dev.off()
    #       }
    #     }
    
  } 
  
  else if (data_betadiv_pooled)
  {
    pdf("Samples_composition.pdf")
    par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
    #bottom left top right
    #par(mar=c(5.1,4.1,4.1,2.1)
    par(mar=c(15.1,10.1,4.1,4.1))
    # old_names = c("AR74","AR73","NBA1","NBA2","NINS","NF21","NH20","NH21","NL11","NL12","NO13","NPA5","NPA6","NL18","P063","P121","P064","P111","P122")
    new_names = c("NH20","NH21","NF21","NL11","NL12","NO13","NL18","NBA1","NBA2","NPA5","NPA6","NINS","AR74","AR73","P121","P122","P111","P063","P064")
    new_names_indices = vector(length=length(new_names),mode="numeric")
    for (ii in 1:length(new_names)) 
      new_names_indices[ii] = which(colnames(data2m) == new_names[ii])
    plot(documents[new_names_indices,rev(sort_normal_topic$ix)[1]],col=terrain.colors(nb_topics)[1],type="p",ann=F,yaxt="n",xaxt="n")
    axis(2, ylim=c(0,1), col='black')
    axis(1, at=1:19, labels = F)
    text(1:19, par("usr")[3], srt = 45, adj = c(1.05,1.75), labels = new_names, xpd = TRUE, cex=1.1)
    for (k in 2:nb_topics)
    {
      lines(documents[new_names_indices,rev(sort_normal_topic$ix)[k]],col=terrain.colors(nb_topics)[k],type="p")
    }
    title(ylab="Plots composition")
    title(paste("Distribution of the ",nb_topics,"\ncomponent communities",sep=""),cex.main=1.7)
    dev.off()
    
    pdf("Samples_composition_barplot.pdf")
    par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
    #bottom left top right
    #par(mar=c(5.1,4.1,4.1,2.1)
    par(mar=c(15.1,10.1,4.1,4.1))
    barplot(t(documents[new_names_indices,rev(sort_normal_topic$ix)]),col=terrain.colors(nb_topics),names.arg=new_names,ylim=c(0,1),ann=F,las=3)
    #axis(1, at=1:19, labels = F)
    #text(1:19, par("usr")[3], srt = 45, adj = c(1.05,1.75), labels = new_names, xpd = TRUE, cex=1.1)
    title(ylab="Plots composition",main=paste("Distribution of the ",nb_topics,"\ncomponent communities",sep=""))
    dev.off()
    
  }
  
  # end of "site repartition information" subsection
}

# end of (best || best_keep) condition
}

# end of loop over j_select
  }

##################
setwd(local_subdirname)

if (Rtopicmodels_VEM && select_real)
{  
  mean_alpha_file = "mean_estimated_alpha.txt"
  write(alpha_est_mean,mean_alpha_file)
}

##################################
# Comparing realizations: output #
##################################

if (select_real && realization_comparison)
{
  realization_comparison_dirname = paste0(local_subdirname,"Realization_comparison_",MOTU_sample_insert,"_",bij_insert)
  if (!(file.exists(realization_comparison_dirname)))
    dir.create(realization_comparison_dirname)
  setwd(realization_comparison_dirname)
  
  # Comparing realizations, averaged over topics
  #######
  # saving the correlation file
  #     save(Correlation_MOTUwise,Correlation_samplewise,KL_samplewise,KL_MOTUwise,file="Topic_KL-correlation_between_best_real_and_following_reals.Rdata")
  #     save(Correlation_samplewise_allRealPairs,Correlation_MOTUwise_allRealPairs,KL_samplewise_allRealPairs,KL_MOTUwise_allRealPairs,file="Topic_KL-correlation_matrices_between_all_real.Rdata")
  save(Correlation,KL,DKL100,file=paste0("SKL-correlation_toBestReal_",MOTU_sample_insert,"_",bij_insert,".Rdata"))
  save(Correlation_allRealPairs,KL_allRealPairs_w_rndzations,KL_allRealPairs,DKL100_allRealPairs,llh_differences_allRealPairs,p_value_allRealPairs,SES_allRealPairs,file=paste0("SKL-correlation_allRealPairs_",MOTU_sample_insert,"_",bij_insert,".Rdata"))
  save(Correlation_allRealPairs_byTopic,Correlation_allRealPairs_byTopic,KL_allRealPairs_byTopic,KL_allRealPairs_randomized_byTopic,SES_allRealPairs_byTopic,DKL100_allRealPairs_byTopic,Var_KL_topic_comparison_randomized_allRealPairs_byTopic,
       p_value_allRealPairs_byTopic,Topic_correspondence_to_best_real,file=paste0("SKL-correlation_allRealPairs_byTopic_",MOTU_sample_insert,"_",bij_insert,".Rdata"))
  
  Stability_file = paste0("Stability_",MOTU_sample_insert,"_",bij_insert,".txt")
  
  write(paste(barcode_insert,"-",nb_topics,"topics -",length_selected_real,"realizations -",MOTU_sample_insert,"-",bij_insert,"\n%%%%%%%%%%%%%\n"),file=Stability_file,append=F)
  write(paste("For all topics :\n%%%%%%%%%%%%%\n"),file=Stability_file,append=T)
  ###
  write(paste("Mean correlation =",mean(Correlation_allRealPairs[which(Correlation_allRealPairs!=0)]),"\n"),file=Stability_file,append=T)
  ###
  write(paste("Mean KL distance =",mean(KL_allRealPairs[which(!is.na(KL_allRealPairs))]),"\n"),file=Stability_file,append=T)
  ###
  # Intercept of SKL = f(llh)
  fitted_llh = Ordered_realizations$x[1]-Ordered_realizations$x[2:length_selected_real]
  fit1 = lm(KL_allRealPairs[1,-1] ~ fitted_llh)
  write(paste("Mean KL distance = f(llh difference) intercept =",fit1$coefficient[1],"\n"),file=Stability_file,append=T)
  ###
  #write(paste("Mean randomized KL distance for",nb_rndzations,"randomizations per topic =",mean(KL_allRealPairs_w_rndzations[lower.tri(KL_allRealPairs_w_rndzations,diag=F)]),"\n"),file=Stability_file,append=T)
  ###
  write(paste("Mean KL distance difference (ES) for",nb_rndzations,"randomizations per topic =",mean(DKL100_allRealPairs[which(!is.na(DKL100_allRealPairs))]),"\n"),file=Stability_file,append=T)
  ###
  nES = t(DKL100_allRealPairs)[lower.tri(DKL100_allRealPairs,diag=F)]/KL_allRealPairs_w_rndzations[lower.tri(KL_allRealPairs_w_rndzations,diag=F)]
  write(paste("Mean normalized KL distance difference (nES) for",nb_rndzations,"randomizations per topic =",mean(nES),"\n"),file=Stability_file,append=T)
  ###
  # Intercept of nES = f(llh) : 
  fit2 = lm(nES[1:length(fitted_llh)] ~ fitted_llh)
  write(paste("nES = f(llh difference) intercept =",fit2$coefficient[1],"\n"),file=Stability_file,append=T)  
  ###
  SES_denominator = sqrt(mean(Var_KL_topic_comparison_randomized_allRealPairs[which(!is.na(Var_KL_topic_comparison_randomized_allRealPairs))]) -
                           mean(KL_allRealPairs_w_rndzations[lower.tri(KL_allRealPairs_w_rndzations,diag=F)])^2)
  SES = mean(DKL100_allRealPairs[which(!is.na(DKL100_allRealPairs))])/SES_denominator
  write(paste("Mean standardized KL distance difference (SES) for",nb_rndzations,"randomizations per topic =",SES,"\n"),file=Stability_file,append=T)
  ###
  SES_topicByTopic = mean(SES_allRealPairs[which(!is.na(SES_allRealPairs))])
  write(paste("Mean standardized KL distance difference (SES) computed topic by topic for",nb_rndzations,"randomizations per topic =",SES_topicByTopic,"\n"),file=Stability_file,append=T)
  ###
  write(paste("Mean p-value for",nb_rndzations,"randomizations per topic =",mean(p_value_allRealPairs[which(!is.na(p_value_allRealPairs))]),"\n"),file=Stability_file,append=T)
  
  if (bij)
    stability_data.frame = as.data.frame(matrix(nrow=nb_topics+1,ncol=11,data=0))
  else
    stability_data.frame = as.data.frame(matrix(nrow=1,ncol=11,data=0))
  rownames(stability_data.frame)[1] = barcode_insert
  stability_data.frame[1,] = c(diversity$Diversity[nb_topics+1],
                               mean(Correlation_allRealPairs[which(Correlation_allRealPairs!=0)]),
                               mean(KL_allRealPairs[which(!is.na(KL_allRealPairs))]),
                               fit1$coefficient[1],
                               mean(DKL100_allRealPairs[which(!is.na(DKL100_allRealPairs))]),
                               mean(nES),
                               fit2$coefficient[1],
                               SES,
                               SES_topicByTopic,
                               mean(p_value_allRealPairs[which(!is.na(p_value_allRealPairs))]),
                               alpha_est_mean)
  
  if (bij)
  {   
    # Comparing realizations topic by topic
    for (k in 1:nb_topics)
    {
      k0 = rev(sort_normal_topic_allreal[[1]]$ix)[k]
      
      p_value = 0
      mean_KL = 0
      mean_DKL100 = 0
      mean_nES = 0
      mean_KL_randomized = 0
      mean_SES_denominator = 0
      mean_correlation = 0
      SES_topicByTopic = 0
      for (j_real in 1:(length_selected_real-1))
      {
        p_value = sum(p_value_allRealPairs_byTopic[[Topic_correspondence_to_best_real[[j_real]][k0]]][j_real,-(1:j_real)])/((length_selected_real-1)*length_selected_real/2) + p_value
        mean_KL = sum(KL_allRealPairs_byTopic[[Topic_correspondence_to_best_real[[j_real]][k0]]][j_real,-(1:j_real)])/((length_selected_real-1)*length_selected_real/2) + mean_KL
        mean_DKL100 = sum(DKL100_allRealPairs_byTopic[[Topic_correspondence_to_best_real[[j_real]][k0]]][j_real,-(1:j_real)])/((length_selected_real-1)*length_selected_real/2) + mean_DKL100 
        mean_nES = sum(DKL100_allRealPairs_byTopic[[Topic_correspondence_to_best_real[[j_real]][k0]]][j_real,-(1:j_real)]/KL_allRealPairs_randomized_byTopic[[Topic_correspondence_to_best_real[[j_real]][k0]]][j_real,-(1:j_real)])/((length_selected_real-1)*length_selected_real/2) + mean_nES
        mean_KL_randomized = sum(KL_allRealPairs_randomized_byTopic[[Topic_correspondence_to_best_real[[j_real]][k0]]][j_real,-(1:j_real)])/((length_selected_real-1)*length_selected_real/2) + mean_KL_randomized
        mean_SES_denominator = sum(Var_KL_topic_comparison_randomized_allRealPairs_byTopic[[Topic_correspondence_to_best_real[[j_real]][k0]]][j_real,-(1:j_real)])/((length_selected_real-1)*length_selected_real/2) + mean_SES_denominator
        SES_topicByTopic = sum(SES_allRealPairs_byTopic[[Topic_correspondence_to_best_real[[j_real]][k0]]][j_real,-(1:j_real)])/((length_selected_real-1)*length_selected_real/2) + SES_topicByTopic
        mean_correlation = sum(Correlation_allRealPairs_byTopic[[Topic_correspondence_to_best_real[[j_real]][k0]]][j_real,-(1:j_real)])/((length_selected_real-1)*length_selected_real/2) + mean_correlation
      }
      mean_SES_denominator = sqrt(mean_SES_denominator - mean_KL_randomized^2)
      SES = mean_DKL100/mean_SES_denominator
      
      # fit of SKL = f(llh) :
      fitted_llh = Ordered_realizations$x[1]-Ordered_realizations$x[2:length_selected_real]
      fit1 = lm(KL_allRealPairs_byTopic[[k0]][1,-1] ~ fitted_llh)
      # fit of nES = f(llh) :
      fit2 = lm(DKL100_allRealPairs_byTopic[[k0]][1,-1]/KL_allRealPairs_randomized_byTopic[[k0]][1,-1] ~ fitted_llh)
      
      if (data_pp && !testdata && (nb_topics == 3))
      {
        i = 0
        i_insert = ""
        while (any(rownames(stability_data.frame)==paste(assemblage_names_vect[k],barcode_insert,i_insert)))
        {
          i = i+1
          i_insert = i
        }
        rownames(stability_data.frame)[k+1] = paste(assemblage_names_vect[k],barcode_insert,i_insert) 
        write(paste0("For topic ",k," - ",assemblage_names_vect[k],
                     " assemblage (rho = ",correlations_to_bacterial_assemblages[k],"):\n%%%%%%%%%%%%%\n"),file=Stability_file,append=T)
      } else
      {
        write(paste0("For topic ",k,":\n%%%%%%%%%%%%%\n"),file=Stability_file,append=T)
        rownames(stability_data.frame)[k+1] = paste("Assemblage",k,barcode_insert)
      }
      
      write(paste("Mean correlation =",mean_correlation,"\n"),file=Stability_file,append=T)
      write(paste("Mean KL distance =",mean_KL,"\n"),file=Stability_file,append=T)
      write(paste("Mean KL distance = f(llh difference) intercept =",fit1$coefficient[1],"\n"),file=Stability_file,append=T)
      write(paste("Mean KL distance difference (ES) for",nb_rndzations,"randomizations per topic =",mean_DKL100,"\n"),file=Stability_file,append=T)
      write(paste("Mean normalized KL distance difference (nES) for",nb_rndzations,"randomizations per topic =",mean_nES,"\n"),file=Stability_file,append=T)
      write(paste("nES = f(llh difference) intercept =",fit2$coefficient[1],"\n"),file=Stability_file,append=T)
      write(paste("Mean standardized KL distance difference (SES) for",nb_rndzations,"randomizations per topic =",SES,"\n"),file=Stability_file,append=T)
      write(paste("Mean standardized KL distance difference (SES) computed topic by topic for",nb_rndzations,"randomizations per topic =",SES_topicByTopic,"\n"),file=Stability_file,append=T)
      write(paste("Mean p-value for",nb_rndzations,"randomizations per topic =",p_value,"\n"),file=Stability_file,append=T)
      
      stability_data.frame[k+1,] = c(diversity$Diversity[k],
                                     mean_correlation,
                                     mean_KL,
                                     fit1$coefficient[1],
                                     mean_DKL100,
                                     mean_nES,
                                     fit2$coefficient[1],
                                     SES,
                                     SES_topicByTopic,
                                     p_value,
                                     NA)
    }
  }
  
  colnames(stability_data.frame) = c("Diversity","Correlation","SKL","SKL=f(llh) intercept","ES","Normalized ES","Normalized ES=f(llh) intercept","Total SES","Mean SES over assemblages","p-value","alpha")
  saveRDS(stability_data.frame,file="stability_data.frame.rds")
  setwd(local_subdirname)
  write.csv(stability_data.frame,paste0("Stability_",MOTU_sample_insert,"_",bij_insert,".csv"))
  setwd(realization_comparison_dirname)
  
  #########
  # Plots #
  #########
  
  pdf(paste0("SKL_",MOTU_sample_insert,"_RealByReal_color2D.pdf"))
  # 2 topics :
  #par(mfrow=c(2,1))
  # 5 topics :
  #par(mfrow=c(3,2))
  # 10 topics
  par(mfrow=c(4,4))
  #par(mar=c(5.1,4.1,4.1,2.1)
  par(mar=c(5.1,5.1,4.1,2.1))
  for (j_select in 2:length_selected_real)
  {
    if (j_select == 2)
      tag = "nd"
    else if (j_select == 3)
      tag = "rd"
    else tag = "th"
    
    # red-grey-blue divergence_hcl(3) color scale: c(142/255,226/255,2/255),c(6/255,226/255,63/255),c(59/255,226/255,165/255)
    #       color2D.matplot(KL_samplewise[[j_select-1]],c(0,1),c(0,0),c(1,0),xrange=c(-1,1),
    #                       extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,xlab=paste("Topics sorted by decreasing abundances,\n",Selected_real[j_select],tag," best realization",sep=""),
    #                       ylab="Topics sorted by decreasing abundances,\n best realization",do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
    #                       border="black",na.color="white",main=paste(Selected_real[j_select],tag," best realization",sep=""))
    color2D.matplot(KL[[j_select-1]],c(0,1),c(0,0),c(1,0),xrange=c(-1,1),
                    extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,xlab=paste("Topics in the ",Selected_real[j_select],tag," best realization",sep=""),
                    ylab="Topics in the \n best realization",do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
                    border="black",na.color="white",main=paste(Selected_real[j_select],tag," best realization",sep=""))
    
  }
  dev.off()
  
  #   pdf(paste0("Correlation_",MOTU_sample_insert,"_RealByReal_color2D.pdf"))
  #   # 2 topics :
  #   #par(mfrow=c(2,1))
  #   # 5 topics :
  #   #par(mfrow=c(3,2))
  #   # 10 topics
  #   par(mfrow=c(4,4))
  #   #par(mar=c(5.1,4.1,4.1,2.1)
  #   par(mar=c(5.1,5.1,4.1,2.1))
  #   for (j_select in 2:length_selected_real)
  #   {
  #     if (j_select == 2)
  #       tag = "nd"
  #     else if (j_select == 3)
  #       tag = "rd"
  #     else tag = "th"
  #     
  #     # red-grey-blue divergence_hcl(3) color scale: c(142/255,226/255,2/255),c(6/255,226/255,63/255),c(59/255,226/255,165/255)
  #     color2D.matplot(Correlation[[j_select-1]],c(0,1),c(0,0),c(1,0),xrange=c(-1,1),
  #                     extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,xlab=paste("Topics in the ",Selected_real[j_select],tag," best realization",sep=""),
  #                     ylab="Topics in the best realization",do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
  #                     border="black",na.color="white",main=paste(Selected_real[j_select],tag," best realization",sep=""))
  #   }
  #   dev.off()
  
  pdf(paste0("SKL_",MOTU_sample_insert,"_allRealMatrix_color2D.pdf"))
  # 2 topics :
  #par(mfrow=c(2,1))
  # 5 topics :
  #par(mfrow=c(3,2))
  # 10 topics
  #par(mar=c(5.1,4.1,4.1,2.1)
  par(mar=c(5.1,5.1,4.1,2.1))
  color2D.matplot(KL_allRealPairs,c(0,1),c(0,0),c(1,0),
                  extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,xlab=paste("Realizations ranked by decreasing likelihood",sep=""),
                  ylab="Realizations ranked by decreasing likelihood",do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
                  border="black",na.color="white")
  dev.off()
  
  if (max(p_value_allRealPairs[!is.na(p_value_allRealPairs)])!=0)
  {
    pdf(paste0("p-value_",MOTU_sample_insert,"_allRealMatrix_color2D.pdf"))
    #par(mar=c(5.1,4.1,4.1,2.1)
    par(mar=c(5.1,5.1,4.1,2.1))
    color2D.matplot(p_value_allRealPairs,c(0,1),c(0,0),c(1,0),
                    extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,xlab=paste("Realizations ranked by decreasing likelihood",sep=""),
                    ylab="Realizations ranked by decreasing likelihood",do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
                    border="black",na.color="white")
    dev.off()
  }
  
  pdf(paste0("SKL_",MOTU_sample_insert,"_allRealMatrix_wRandomizations_color2D.pdf"))
  # 2 topics :
  #par(mfrow=c(2,1))
  # 5 topics :
  #par(mfrow=c(3,2))
  # 10 topics
  #par(mar=c(5.1,4.1,4.1,2.1)
  par(mar=c(5.1,5.1,4.1,2.1))
  color2D.matplot(KL_allRealPairs_w_rndzations,c(0,1),c(0,0),c(1,0),
                  extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,xlab=paste("Realizations ranked by decreasing likelihood",sep=""),
                  ylab="Realizations ranked by decreasing likelihood",do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
                  border="black",na.color="white")
  dev.off()
  
  #   pdf(paste0("Correlation_",MOTU_sample_insert,"_allRealMatrix_color2D.pdf"))
  #   # 2 topics :
  #   #par(mfrow=c(2,1))
  #   # 5 topics :
  #   #par(mfrow=c(3,2))
  #   # 10 topics
  #   #par(mar=c(5.1,4.1,4.1,2.1)
  #   par(mar=c(5.1,5.1,4.1,2.1))
  #   
  #   color2D.matplot(Correlation_allRealPairs[-nrow(Correlation_allRealPairs),-1],c(0,1),c(0,0),c(1,0),
  #                   extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,xlab=paste("Realizations ranked by decreasing likelihood",sep=""),
  #                   ylab="Realizations ranked by decreasing likelihood",do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
  #                   border="black",na.color="white")
  #   dev.off()
  
  pdf(paste0("SKL_",MOTU_sample_insert,"_allRealPairs_vs_llh-diff.pdf"))  
  par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
  # par(mar = c(bottom, left, top, right))
  par(mar = c(5, 5, 4, 3) + 0.1)
  #     plot(seq(2,length_selected_real,1),apply(Topic_comparison,1,mean),ylim=c(0,1),type="p", main = "Averaged Hellinger similarity between\n the component communities in different realizations",
  #          xlab = "Realizations ranked by decreasing likelihood", ylab = "Hellinger similarity")   
  plot(llh_differences_allRealPairs[!is.na(llh_differences_allRealPairs)],KL_allRealPairs[!is.na(KL_allRealPairs)],type="p", main = "Averaged KL distance between\n pairs of realizations",
       xlab = "Loglikelihood difference", ylab = "KL distance")
  dev.off()
  
  pdf(paste0("ES_",MOTU_sample_insert,"_allRealPairs_vs_llh-diff.pdf"))  
  par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
  # par(mar = c(bottom, left, top, right))
  par(mar = c(5, 5, 4, 3) + 0.1)
  #     plot(seq(2,length_selected_real,1),apply(Topic_comparison,1,mean),ylim=c(0,1),type="p", main = "Averaged Hellinger similarity between\n the component communities in different realizations",
  #          xlab = "Realizations ranked by decreasing likelihood", ylab = "Hellinger similarity")  
  if (length(which(!is.na(DKL100_allRealPairs))) != length(which(!is.na(llh_differences_allRealPairs))))
    stop("\nError in DKL100_allRealPairs values")
  plot(llh_differences_allRealPairs[!is.na(llh_differences_allRealPairs)],DKL100_allRealPairs[!is.na(DKL100_allRealPairs)],type="p", main = "Averaged KL distance difference between\n pairs of realizations",
       xlab = "Loglikelihood difference", ylab = "Difference of KL distance")
  dev.off()
  
  pdf(paste0("nES_",MOTU_sample_insert,"_allRealPairs_vs_llh-diff.pdf"))  
  par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
  # par(mar = c(bottom, left, top, right))
  par(mar = c(5, 5, 4, 3) + 0.1)
  #     plot(seq(2,length_selected_real,1),apply(Topic_comparison,1,mean),ylim=c(0,1),type="p", main = "Averaged Hellinger similarity between\n the component communities in different realizations",
  #          xlab = "Realizations ranked by decreasing likelihood", ylab = "Hellinger similarity")  
  if (length(which(!is.na(DKL100_allRealPairs))) != length(which(!is.na(llh_differences_allRealPairs))))
    stop("\nError in DKL100_allRealPairs values")
  plot(t(llh_differences_allRealPairs)[lower.tri(llh_differences_allRealPairs,diag=F)],nES,type="p", main = "Normalized KL distance difference between\n pairs of realizations",
       xlab = "Loglikelihood difference", ylab = "Normalized Effect Size")
  dev.off()
  
  if (bij)
  {
    pdf(paste0("SKL_",MOTU_sample_insert,"_toBestReal_byTopic_vs_llh-diff.pdf"))  
    par(mfrow=c(2,2))
    for (k in 1:nb_topics)
    {
      par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
      # par(mar = c(bottom, left, top, right))
      par(mar = c(5, 5, 4, 3) + 0.1)
      #     plot(seq(2,length_selected_real,1),apply(Topic_comparison,1,mean),ylim=c(0,1),type="p", main = "Averaged Hellinger similarity between\n the component communities in different realizations",
      #          xlab = "Realizations ranked by decreasing likelihood", ylab = "Hellinger similarity")   
      plot(Ordered_realizations$x[1]-Ordered_realizations$x[2:length_selected_real],KL_allRealPairs_byTopic[[rev(sort_normal_topic_allreal[[1]]$ix)[k]]][1,-1],
           type="p", main = paste("Averaged KL distance\n assemblage", k),
           xlab = "Loglikelihood difference\n with best realization", ylab = "KL distance\n to best realization")
    }
    dev.off()
    
    pdf(paste0("SKL_",MOTU_sample_insert,"_toBestReal_byTopic_vs_rank.pdf"))  
    par(mfrow=c(2,2))
    for (k in 1:nb_topics)
    {
      par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
      # par(mar = c(bottom, left, top, right))
      par(mar = c(5, 5, 4, 3) + 0.1)
      #     plot(seq(2,length_selected_real,1),apply(Topic_comparison,1,mean),ylim=c(0,1),type="p", main = "Averaged Hellinger similarity between\n the component communities in different realizations",
      #          xlab = "Realizations ranked by decreasing likelihood", ylab = "Hellinger similarity")   
      plot(seq(2,length_selected_real,1),KL_allRealPairs_byTopic[[rev(sort_normal_topic_allreal[[1]]$ix)[k]]][1,-1],
           type="p", main = paste("Averaged KL distance\n assemblage", k),
           xlab = "Realizations sorted\n by decreasing llh", ylab = "KL distance\n to best realization")
    }
    dev.off()
    
    nES_byTopic = list()
    for (k in 1:nb_topics)
      nES_byTopic[[k]] = t(DKL100_allRealPairs_byTopic[[k]])[lower.tri(DKL100_allRealPairs_byTopic[[k]],diag=F)]/t(KL_allRealPairs_randomized_byTopic[[k]])[lower.tri(KL_allRealPairs_randomized_byTopic[[k]],diag=F)]
    
    pdf(paste0("nES_",MOTU_sample_insert,"_toBestReal_byTopic_vs_llh-diff.pdf"))  
    par(mfrow=c(2,2))
    for (k in 1:nb_topics)
    {
      par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
      # par(mar = c(bottom, left, top, right))
      par(mar = c(5, 5, 4, 3) + 0.1)
      #     plot(seq(2,length_selected_real,1),apply(Topic_comparison,1,mean),ylim=c(0,1),type="p", main = "Averaged Hellinger similarity between\n the component communities in different realizations",
      #          xlab = "Realizations ranked by decreasing likelihood", ylab = "Hellinger similarity")   
      plot(Ordered_realizations$x[1]-Ordered_realizations$x[2:length_selected_real],nES_byTopic[[rev(sort_normal_topic_allreal[[1]]$ix)[k]]][1:(length_selected_real-1)],
           type="p", main = paste("Averaged KL distance\n assemblage", k),
           xlab = "Loglikelihood difference\n with best realization", ylab = "Normalized Effect Size")
    }
    dev.off()
  }
  
  pdf(paste0("SKL_",MOTU_sample_insert,"_toBestReal_vs_llh-diff.pdf")) 
  par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
  # par(mar = c(bottom, left, top, right))
  par(mar = c(5, 5, 4, 3) + 0.1)
  #     plot(seq(2,length_selected_real,1),apply(Topic_comparison,1,mean),ylim=c(0,1),type="p", main = "Averaged Hellinger similarity between\n the component communities in different realizations",
  #          xlab = "Realizations ranked by decreasing likelihood", ylab = "Hellinger similarity")
  plot(Ordered_realizations$x[1]-Ordered_realizations$x[2:length_selected_real],KL_allRealPairs[1,-1],type="p", main = "Averaged KL distance between\n the component communities in different realizations",
       xlab = "Loglikelihood difference with best realization", ylab = "KL distance to best realization")
  dev.off()
  
  pdf(paste0("SKL_",MOTU_sample_insert,"_toBestReal_vs_rank.pdf"))  
  par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
  # par(mar = c(bottom, left, top, right))
  par(mar = c(5, 5, 4, 3) + 0.1)
  #     plot(seq(2,length_selected_real,1),apply(Topic_comparison,1,mean),ylim=c(0,1),type="p", main = "Averaged Hellinger similarity between\n the component communities in different realizations",
  #          xlab = "Realizations ranked by decreasing likelihood", ylab = "Hellinger similarity")
  plot(seq(2,length_selected_real,1),KL_allRealPairs[1,-1],type="p", main = "Averaged KL distance between\n the component communities in different realizations",
       xlab = "Realizations sorted by decreasing llh", ylab = "KL distance to best realization")
  dev.off()
  
  pdf(paste0("ES_",MOTU_sample_insert,"_toBestReal_vs_llh-diff.pdf"))  
  par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
  # par(mar = c(bottom, left, top, right))
  par(mar = c(5, 5, 4, 3) + 0.1)
  plot(Ordered_realizations$x[1]-Ordered_realizations$x[2:length_selected_real],DKL100_allRealPairs[1,-1],type="p", main = "Averaged KL distance difference between\n the component communities in different realizations",
       xlab = "Loglikelihood difference with best realization", ylab = "Difference in KL distance to best realization")
  dev.off()
  
  pdf(paste0("nES_",MOTU_sample_insert,"_toBestReal_vs_llh-diff.pdf"))  
  par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
  # par(mar = c(bottom, left, top, right))
  par(mar = c(5, 5, 4, 3) + 0.1)
  plot(Ordered_realizations$x[1]-Ordered_realizations$x[2:length_selected_real],nES[1:(length_selected_real-1)],type="p", main = "Normalized veraged KL distance difference between\n the component communities in different realizations",
       xlab = "Loglikelihood difference with best realization", ylab = "Normalized Effect Size")
  dev.off()
  
  pdf(paste0("nES_",MOTU_sample_insert,"_toBestReal_vs_rank.pdf"))  
  par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
  # par(mar = c(bottom, left, top, right))
  par(mar = c(5, 5, 4, 3) + 0.1)
  #     plot(seq(2,length_selected_real,1),apply(Topic_comparison,1,mean),ylim=c(0,1),type="p", main = "Averaged Hellinger similarity between\n the component communities in different realizations",
  #          xlab = "Realizations ranked by decreasing likelihood", ylab = "Hellinger similarity")
  plot(seq(2,length_selected_real,1),nES[1:(length_selected_real-1)],type="p", main = "Normalized averaged KL distance between\n the component communities in different realizations",
       xlab = "Realizations sorted by decreasing llh", ylab = "Normalized Effect Size")
  dev.off()
  
  pdf(paste0("p-value_",MOTU_sample_insert,"_toBestReal_vs_llh-diff.pdf"))  
  par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
  # par(mar = c(bottom, left, top, right))
  par(mar = c(5, 5, 4, 3) + 0.1)
  plot(Ordered_realizations$x[1]-Ordered_realizations$x[2:length_selected_real],p_value_allRealPairs[1,-1],type="p", main = "Averaged stability p-value between the best real. and the following",
       xlab = "Llh difference with best realization", ylab = "p-value wrt best realization", xlim=range(p_value_allRealPairs[1,-1]))
  dev.off()
  
  pdf(paste0("Correlation_",MOTU_sample_insert,"_toBestReal_vs_llh-diff.pdf"))  
  par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
  plot(Ordered_realizations$x[1]-Ordered_realizations$x[2:length_selected_real],Correlation_allRealPairs[1,-1],type="p", main = "Averaged correlation between\n the component communities in different realizations",
       xlab = "Loglikelihood difference with best realization", ylab = "Correlation coefficient with best realization")
  dev.off()
  
  # Comparing llh and AIC of realizations:
  #######################################
  
  plotname = paste("Realization_comparison_llh_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,".pdf",sep="")
  pdf(plotname)
  #par(mar=c(5.1,4.1,4.1,2.1))
  par(mar=c(5.1,5.1,4.1,2.1))
  
  LLH_final_real1_sorted = sort.int(LLH_final_real1,decreasing=T,index.return=T)
  plot(seq(1,nb_real,1),LLH_final_real1_sorted$x,ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
  #lines(seq(1,nb_real,1),Akaike_weights_AIC1[Akaike_weights_llh_sorted$ix],lwd=2,type="p",col="forestgreen")
  #lines(seq(1,nb_real,1),Akaike_weights_AIC2[Akaike_weights_llh_sorted$ix],lwd=2,type="p",col="blue")
  
  legend(x="topright",legend=c(filename_insert,alpha_insert,paste("delta =",delta),paste("nb topics =",nb_topics),paste("nb realizations =",nb_real)),col="black",inset=0)
  title("Loglikelihood for different realizations",cex.main=1.7)
  title(xlab="Ranked realizations",ylab="Loglikelihood value",cex.lab=1.5)
  dev.off()
  #####################
  
  #####################
  plotname = paste("Realization_llh_histogram_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,".pdf",sep="")
  pdf(plotname)
  #par(mar=c(5.1,4.1,4.1,2.1))
  par(mar=c(5.1,5.1,4.1,2.1))
  
  hist(LLH_final_real1, breaks = 10, freq=T, ann=F, xaxp = c(0, 1, 10))      
  
  legend(x="topright",legend=c(filename_insert,alpha_insert,paste("delta =",delta),paste("nb topics =",nb_topics),paste("nb realizations =",nb_real)),col="black",inset=0)
  title("Loglikelihood histogram for different realizations",cex.main=1.7)
  title(ylab="Number of realizations",xlab="Loglikelihood value",cex.lab=1.5)
  dev.off()
  #####################
  
  # the probability weights are the same whether they are calculated with the llh or with AIC, since all realizations have the same number of parameters
  #####################
  plotname = paste("Realization_comparison_AIC_weight_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,".pdf",sep="")
  pdf(plotname)
  #par(mar=c(5.1,4.1,4.1,2.1))
  par(mar=c(5.1,5.1,4.1,2.1))
  
  Akaike_weights_llh_sorted = sort.int(Akaike_weights_llh,decreasing=T,index.return=T)
  plot(seq(1,nb_real,1),Akaike_weights_llh_sorted$x,ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
  #lines(seq(1,nb_real,1),Akaike_weights_AIC1[Akaike_weights_llh_sorted$ix],lwd=2,type="p",col="forestgreen")
  #lines(seq(1,nb_real,1),Akaike_weights_AIC2[Akaike_weights_llh_sorted$ix],lwd=2,type="p",col="blue")
  
  legend(x="bottomright",legend=c(filename_insert,alpha_insert,paste("delta =",delta),paste("nb topics =",nb_topics),paste("nb realizations =",nb_real)),col="black",inset=0)
  title("Akaike weights of different realizations",cex.main=1.7)
  title(xlab="Ranked realizations",ylab="Akaike weight value",cex.lab=1.5)
  dev.off()
  #####################
  
  #####################
  plotname = paste("Realization_comparison_log10(AIC_weight)_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,".pdf",sep="")
  pdf(plotname)
  #par(mar=c(5.1,4.1,4.1,2.1))
  par(mar=c(5.1,5.1,4.1,2.1))
  
  Akaike_weights_llh_sorted = sort.int(Akaike_weights_llh,decreasing=T,index.return=T)
  Akaike_weights_llh_log10 = vector(length=nb_real,mode="numeric")
  #Akaike_weights_AIC1_log10 = vector(length=nb_real,mode="numeric")
  #Akaike_weights_AIC2_log10 = vector(length=nb_real,mode="numeric")
  Akaike_weights_llh_log10[Akaike_weights_llh==0] = NA
  Akaike_weights_llh_log10[Akaike_weights_llh!=0] = log10(Akaike_weights_llh[Akaike_weights_llh!=0])
  #Akaike_weights_AIC1_log10[Akaike_weights_AIC1==0] = NA
  #Akaike_weights_AIC1_log10[Akaike_weights_AIC1!=0] = log10(Akaike_weights_AIC1[Akaike_weights_AIC1!=0])
  #Akaike_weights_AIC2_log10[Akaike_weights_AIC2==0] = NA
  #Akaike_weights_AIC2_log10[Akaike_weights_AIC2!=0] = log10(Akaike_weights_AIC2[Akaike_weights_AIC2!=0])
  
  plot(seq(1,nb_real,1),Akaike_weights_llh_log10[Akaike_weights_llh_sorted$ix],ann=F,cex.axis=1.5,lwd=2,type="p",col="black",
       #ylim=range(Akaike_weights_llh_log10[Akaike_weights_llh!=0],Akaike_weights_AIC1_log10[Akaike_weights_AIC1!=0],Akaike_weights_AIC2_log10[Akaike_weights_AIC2!=0]))
       ylim=range(Akaike_weights_llh_log10[Akaike_weights_llh!=0]))
  #lines(seq(1,nb_real,1),Akaike_weights_AIC1_log10[Akaike_weights_llh_sorted$ix],lwd=2,type="p",col="forestgreen")
  #lines(seq(1,nb_real,1),Akaike_weights_AIC2_log10[Akaike_weights_llh_sorted$ix],lwd=2,type="p",col="blue")
  
  legend(x="bottomright",legend=c(filename_insert,alpha_insert,paste("delta =",delta),paste("nb topics =",nb_topics),paste("nb realizations =",nb_real)),col="black",inset=0)
  title("Akaike weights of different realizations",cex.main=1.7)
  title(xlab="Ranked realizations",ylab="Akaike weight value (log10)",cex.lab=1.5)
  dev.off()
  #######################
  
  # Back to the main directory
  setwd(local_subdirname)
  
  # end condition (select_real && realization_comparison)
}

#######################################
# Comparing realizations for testdata #
#######################################

if (testdata)
{
  if (realization_comparison)
  {
    ##############
    pdf("Topic_ordered_by_site-normalized_abundance_sample_distribution_allreal.pdf")
    #par(mfrow=c(2,2))
    par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
    for (k in 1:nb_topics)
    {
      plot(documents_allreal[[1]][,k],type="l",ylim=c(0,1),
           main=paste("Community #",k,sep=""),
           xlab="Samples",
           ylab="Proportion of the community in the sample")
      #     plot(documents_allreal[[1]][,rev(sort_normal_topic_allreal[[1]]$ix)[k]],type="l",ylim=c(0,1),
      #          main=paste("Community #",rev(sort_normal_topic_allreal[[1]]$ix)[k],sep=""),
      #          xlab="Samples",
      #          ylab="Proportion of the community in the sample")
      for (j_select in 2:length_selected_real)
      {
        index0 = which(Correlation_old[[j_select-1]][k,] == max(Correlation_old[[j_select-1]][k,]))
        lines(documents_allreal[[j_select]][,index0],type="l")
        #lines(documents_allreal[[j_select]][,rev(sort_normal_topic_allreal[[j_select]]$ix)[index0]],type="l")
      }
      #index0 = which(Correlation_true_documents[k,] == max(Correlation_true_documents[k,]))
      #lines(true_documents[,index0],type="l",lty=2)
      for (k in 1:true_nb_topics)
        lines(true_documents[,k],type="l",lty=2)
    }
    dev.off()    
    
    #################
    pdf("Topic_ordered_by_site-normalized_abundance_sample_distribution_allreal_oneplot.pdf")
    #par(mfrow=c(2,2))
    par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
    for (k in 1:nb_topics)
    {
      if (k == 1)
      {
        plot(documents_allreal[[1]][,k],type="l",ylim=c(0,1),col=rainbow_hcl(nb_topics)[k],
             main="Composition of samples",
             xlab="Samples",
             ylab="Proportion of the community in the sample")
        #       plot(documents_allreal[[1]][,rev(sort_normal_topic_allreal[[1]]$ix)[k]],type="l",ylim=c(0,1),col=rainbow_hcl(nb_topics)[k],
        #            main="Composition of samples",
        #            xlab="Samples",
        #            ylab="Proportion of the community in the sample")
      } else
        lines(documents_allreal[[1]][,k],type="l",col=rainbow_hcl(nb_topics)[k])
      #       lines(documents_allreal[[1]][,rev(sort_normal_topic_allreal[[1]]$ix)[k]],type="l",col=rainbow_hcl(nb_topics)[k])
      for (j_select in 2:length_selected_real)
      {
        index0 = which(Correlation_old[[j_select-1]][k,] == max(Correlation_old[[j_select-1]][k,]))
        lines(documents_allreal[[j_select]][,index0],type="l",col=rainbow_hcl(nb_topics)[k])
        #       lines(documents_allreal[[j_select]][,rev(sort_normal_topic_allreal[[j_select]]$ix)[index0]],type="l",col=rainbow_hcl(nb_topics)[k])
      }
      #   index0 = which(Correlation_true_documents[k,] == max(Correlation_true_documents[k,]))
      #   lines(true_documents[,index0],type="l",col="black",lty=2)
    }
    for (k in 1:true_nb_topics)
      lines(true_documents[,k],type="l",lty=2)
    dev.off()   
    
    # Computing mean error
    ####################
    error_tot = 0
    error_real = vector(length=length_selected_real,mode="numeric")
    documents_mean = matrix(nrow=nb_topics,ncol=nb_doc,data=0)
    for (k in 1:nb_topics)
    {
      index0 = which(Correlation_true_documents[k,] == max(Correlation_true_documents[k,]))
      error_real[1] = mean(abs(documents_allreal[[1]][which(true_documents[,index0]>0),k] 
                               - true_documents[which(true_documents[,index0]>0),index0]))/nb_topics + error_real[1]
      #     error_real[1] = mean(abs(documents_allreal[[1]][which(true_documents[,index0]>0),rev(sort_normal_topic_allreal[[1]]$ix)[k]] 
      #                              - true_documents[which(true_documents[,index0]>0),index0]))/nb_topics + error_real[1]
      documents_mean[k,] = documents_allreal[[1]][,k]/length_selected_real 
      #     documents_mean[k,] = documents_allreal[[1]][,rev(sort_normal_topic_allreal[[1]]$ix)[k]]/length_selected_real    
      for (j_select in 2:length_selected_real)
      {
        index1 = which(Correlation_old[[j_select-1]][k,] == max(Correlation_old[[j_select-1]][k,]))
        error_real[j_select] = mean(abs(documents_allreal[[j_select]][which(true_documents[,index0]>0),index1] 
                                        - true_documents[which(true_documents[,index0]>0),index0]))/nb_topics + error_real[j_select]
        #       error_real[j_select] = mean(abs(documents_allreal[[j_select]][which(true_documents[,index0]>0),rev(sort_normal_topic_allreal[[j_select]]$ix)[index1]] 
        #                                       - true_documents[which(true_documents[,index0]>0),index0]))/nb_topics + error_real[j_select]
        documents_mean[k,] = documents_allreal[[j_select]][,index1]/length_selected_real + documents_mean[k,]
        #       documents_mean[k,] = documents_allreal[[j_select]][,rev(sort_normal_topic_allreal[[j_select]]$ix)[index1]]/length_selected_real + documents_mean[k,]
      }  
    }
    
    for (j_select in 1:length_selected_real)
    {
      error_tot = error_real[j_select]/length_selected_real + error_tot
      error_file = "Error.txt"
      if (j_select == 1)
        write(paste("Error for each run:\n",error_real[1],sep=""),file=error_file,append=F) 
      else
        write(error_real[j_select],file=error_file,append=T)
      if (j_select == length_selected_real)
        write(paste("Mean error:\n",error_tot,sep=""),file=error_file,append=T)  
    }
    
    # plotting the mean error around the mean result over the different runs 
    #################
    pdf("Topic_sample_distribution_error.pdf")
    #par(mfrow=c(2,2))
    par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
    for (k in 1:nb_topics)
    {
      if (k == 1)
      {
        plot(documents_mean[k,],type="l",ylim=c(0,1),col=rainbow_hcl(nb_topics)[k],
             main="Composition of samples",
             xlab="Samples",
             ylab="Proportion of the community in the sample")
      } else
        lines(documents_mean[k,],type="l",col=rainbow_hcl(nb_topics)[k])
      lines(documents_mean[k,]+error_tot*documents_mean[k,],type="l",col=rainbow_hcl(nb_topics)[k],lty=2)
      lines(documents_mean[k,]-error_tot*documents_mean[k,],type="l",col=rainbow_hcl(nb_topics)[k],lty=2)
      #   index0 = which(Correlation_true_documents[k,] == max(Correlation_true_documents[k,]))
      #   lines(true_documents[,index0],type="l",col="black",lty=2)
    }
    for (k in 1:true_nb_topics)
      lines(true_documents[,k],type="l",col="black",lty=2)
    dev.off()  
  }
  
  if (nb_topics == true_nb_topics)
  {
    # Saving the KL comparison between the true topics and the LDA topics:
    KL_true_topics_file = paste0("Goodness-of-fit_nb_rndzations",nb_rndzations_true_documents,".txt")
    
    write("KL distance between true topics and LDA topics:",file=KL_true_topics_file,append=F) 
    for (k in 1:nb_topics)
      write(KL_topic_comparison_true_documents[k],file=KL_true_topics_file,append=T)
    write(paste("Mean =",mean(KL_topic_comparison_true_documents)),file=KL_true_topics_file,append=T)
    
    write(paste("\nMean KL distance between true topics and",nb_rndzations_true_documents,"randomizations of the LDA topics:"),file=KL_true_topics_file,append=T)
    for (k in 1:nb_topics)
      write(Mean_KL_topic_comparison_true_documents_randomized[k],file=KL_true_topics_file,append=T)
    write(paste("Mean =",mean(Mean_KL_topic_comparison_true_documents_randomized)),file=KL_true_topics_file,append=T)
    
    write(paste("\nDifference between the KL distance between true topics and LDA topics and the mean KL distance (ES) over",nb_rndzations_true_documents,"randomizations:"),file=KL_true_topics_file,append=T)
    for (k in 1:nb_topics)
      write(DKL100_true_documents[k],file=KL_true_topics_file,append=T)
    write(paste("Mean =",mean(DKL100_true_documents)),file=KL_true_topics_file,append=T)
    
    write(paste("\nNormalized difference between the KL distance between true topics and LDA topics and the mean KL distance (nES) over",nb_rndzations_true_documents,"randomizations:"),file=KL_true_topics_file,append=T)
    for (k in 1:nb_topics)
      write(nES_true_documents[k],file=KL_true_topics_file,append=T)
    write(paste("Mean =",mean(nES_true_documents)),file=KL_true_topics_file,append=T)
    
    write(paste("\np-value over",nb_rndzations_true_documents,"randomizations of the KL distance between true topics and LDA topics:"),file=KL_true_topics_file,append=T)
    for (k in 1:nb_topics)
      write(p_value_true_documents[k],file=KL_true_topics_file,append=T)
    write(paste("Mean =",mean(p_value_true_documents)),file=KL_true_topics_file,append=T)
  }
}

######################################
# Comparing realizations for betadiv #
######################################

if (data_betadiv_pooled && realization_comparison)
{
  pdf("Samples_composition_allreal_oneplot.pdf")
  #par(mfrow=c(2,2))
  par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
  #bottom left top right
  #par(mar=c(5.1,4.1,4.1,2.1)
  par(mar=c(15.1,10.1,4.1,4.1))
  for (k in 1:nb_topics)
  {
    if (k == 1)
    {
      plot(documents_allreal[[1]][new_names_indices,rev(sort_normal_topic_allreal[[1]]$ix)[k]],type="p",ylim=c(0,1),col=terrain.colors(nb_topics)[k],
           main=paste("Distribution of the ",nb_topics,"\ncomponent communities",sep=""),yaxt="n",xaxt="n",ann=F)
      axis(2, ylim=c(0,1), col='black')
      axis(1, at=1:19, labels = F)
      text(1:19, par("usr")[3], srt = 45, adj = c(1.05,1.75), labels = new_names, xpd = TRUE, cex=1.1)
    } else
      lines(documents_allreal[[1]][new_names_indices,rev(sort_normal_topic_allreal[[1]]$ix)[k]],type="p",col=terrain.colors(nb_topics)[k])
    for (j_select in 2:length_selected_real)
    {
      index0 = which(Correlation_old[[j_select-1]][k,] == max(Correlation_old[[j_select-1]][k,]))
      lines(documents_allreal[[j_select]][new_names_indices,rev(sort_normal_topic_allreal[[j_select]]$ix)[index0]],type="p",col=terrain.colors(nb_topics)[k])
    }
  }
  title(ylab="Plots composition")
  dev.off()    
}

# Plotting stacked bar charts
# x<-data.frame(row=rep(c("A","B","C","D","E","F","G","H"),each=12), col=rep(c(1:12),8), n=100, l1=runif(96),l2=runif(96),l3=runif(96),l4=runif(96),l5=runif(96),l6=runif(96))
# library(reshape2)
# y<-melt(x,id.var=c("row","col","n"))
# library(ggplot2)
# 
# ggplot(y, aes(as.factor(col),value, fill=variable))+geom_bar(stat="identity",position="stack")+facet_grid(~row)  
# 
# Alternatively : function "barplot"












###### end of the !mpar condition
} else if (mpar)
{
  #####
  
  # subdirname = "several_topic_nb/"
  # if (!(file.exists(dirname)))
  #   {dir.create(dirname)}
  # setwd(subdirname)
  
  save(nb_topics_range,LLH_final0,AIC0,AIC1,file="AIC-llh.Rdata")
  
  if (mnb_topics)
  {
    if (local)
    {
      if (Rtopicmodels_Gibbs)
        save(Result_mpar,LLH_final0,AIC0,nb_topics_range,alpha_insert,delta,nb_real,nb_iter,Prop_OTU_removed,Prop_reads_removed,file=filename)
      else if (Rtopicmodels_VEM)
        save(Result_mpar,LLH_final0,AIC0,AIC1,AIC2,AIC3,nb_topics_range,alpha_insert,nb_real,alpha_est0,Prop_OTU_removed,Prop_reads_removed,file=filename)
    }
    ########
    plotname = paste("llh_compare_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,"_nb_iter",nb_iter,".pdf",sep="")
    pdf(plotname)
    #par(mar=c(5.1,4.1,4.1,2.1))
    par(mar=c(5.1,5.1,4.1,2.1))
    
    errbar(nb_topics_range,LLH_final0[1,],LLH_final0[1,]+LLH_final0[2,],LLH_final0[1,]-LLH_final0[2,],ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
    
    #legend(x="bottomright",legend=c(filename_insert,alpha_insert,paste("delta =",delta),paste("nb realizations =",nb_real),paste("nb iterations =",nb_iter)),col="black",inset=0)
    title("Final log-likelihood value vs number of topics",cex.main=1.7)
    title(xlab="K",ylab="log-likelihood value",cex.lab=1.5)
    dev.off()
    ########
    
    ########
    plotname = paste("llh_compare_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,"_nb_iter",nb_iter,"_allpoints.pdf",sep="")
    pdf(plotname)
    #par(mar=c(5.1,4.1,4.1,2.1))
    par(mar=c(5.1,5.1,4.1,2.1))
    
    plot(nb_topics_range,LLH_final1[1,],ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
    for (j in 1:nb_real)
      lines(nb_topics_range,LLH_final1[j,],lwd=2,type="p",col="black")
    
    #legend(x="bottomright",legend=c(filename_insert,alpha_insert,paste("delta =",delta),paste("nb realizations =",nb_real),paste("nb iterations =",nb_iter)),col="black",inset=0)
    title("Final log-likelihood value vs number of topics",cex.main=1.7)
    title(xlab="K",ylab="log-likelihood value",cex.lab=1.5)
    dev.off()
    ########
    
    if (Rtopicmodels_VEM)
    {
      #######
      plotname = paste("AIC_compare_nb_topics*nb_terms+1_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,"_nb_iter",nb_iter,".pdf",sep="")
      pdf(plotname)
      #par(mar=c(5.1,4.1,4.1,2.1))
      par(mar=c(5.1,5.1,4.1,2.1))
      
      errbar(nb_topics_range,AIC0[5,],AIC0[5,]+AIC0[6,],AIC0[5,]-AIC0[6,],ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
      
      legend(x="bottomright",legend=c(filename_insert,alpha_insert,paste("delta =",delta),paste("nb realizations =",nb_real)),col="black",inset=0)
      title("AIC value vs number of topics",cex.main=1.7)
      title(xlab="Number of topics",ylab="AIC value",cex.lab=1.5)
      dev.off()
      ########
      
      #######
      plotname = paste("AIC_compare_nb_topics*nb_terms+1_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,"_nb_iter",nb_iter,"_allpoints.pdf",sep="")
      pdf(plotname)
      #par(mar=c(5.1,4.1,4.1,2.1))
      par(mar=c(5.1,5.1,4.1,2.1))
      
      plot(nb_topics_range,AIC3[1,],ann=F,cex.axis=1.5,lwd=2,type="p",col="black",ylim=range(AIC3))
      for (j in 1:nb_real)
        lines(nb_topics_range,AIC3[j,],lwd=2,type="p",col="black")
      
      legend(x="bottomright",legend=c(filename_insert,alpha_insert,paste("delta =",delta),paste("nb realizations =",nb_real)),col="black",inset=0)
      title("AIC value vs number of topics",cex.main=1.7)
      title(xlab="Number of topics",ylab="AIC value",cex.lab=1.5)
      dev.off()
      ########
    }
    
    #######
    plotname = paste("AIC_compare_nb_topics*(nb_doc+nb_terms)_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,"_nb_iter",nb_iter,".pdf",sep="")
    pdf(plotname)
    #par(mar=c(5.1,4.1,4.1,2.1))
    par(mar=c(5.1,5.1,4.1,2.1))
    
    errbar(nb_topics_range,AIC0[1,],AIC0[1,]+AIC0[2,],AIC0[1,]-AIC0[2,],ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
    
    legend(x="bottomright",legend=c(filename_insert,alpha_insert,paste("delta =",delta),paste("nb realizations =",nb_real),paste("nb iterations =",nb_iter)),col="black",inset=0)
    title("AIC value vs number of topics",cex.main=1.7)
    title(xlab="Number of topics",ylab="AIC value",cex.lab=1.5)
    dev.off()
    ########
    
    #######
    plotname = paste("AIC_compare_nb_topics*(nb_doc+nb_terms)_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,"_nb_iter",nb_iter,"_allpoints.pdf",sep="")
    pdf(plotname)
    #par(mar=c(5.1,4.1,4.1,2.1))
    par(mar=c(5.1,5.1,4.1,2.1))
    
    plot(nb_topics_range,AIC1[1,],ann=F,cex.axis=1.5,lwd=2,type="p",col="black",ylim=range(AIC1))
    for (j in 1:nb_real)
      lines(nb_topics_range,AIC1[j,],lwd=2,type="p",col="black")
    
    legend(x="bottomright",legend=c(filename_insert,alpha_insert,paste("delta =",delta),paste("nb realizations =",nb_real),paste("nb iterations =",nb_iter)),col="black",inset=0)
    title("AIC value vs number of topics",cex.main=1.7)
    title(xlab="Number of topics",ylab="AIC value",cex.lab=1.5)
    dev.off()
    ########
    
    #######
    plotname = paste("AIC_compare_nb_words_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,"_nb_iter",nb_iter,".pdf",sep="")
    pdf(plotname)
    #par(mar=c(5.1,4.1,4.1,2.1))
    par(mar=c(5.1,5.1,4.1,2.1))
    
    errbar(nb_topics_range,AIC0[3,],AIC0[3,]+AIC0[4,],AIC0[3,]-AIC0[4,],ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
    
    legend(x="bottomright",legend=c(filename_insert,alpha_insert,paste("delta =",delta),paste("nb realizations =",nb_real),paste("nb iterations =",nb_iter)),col="black",inset=0)
    title("AIC value vs number of topics",cex.main=1.7)
    title(xlab="Number of topics",ylab="AIC value",cex.lab=1.5)
    dev.off()
    ########
    
    #######
    plotname = paste("AIC_compare_nb_words_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,"_nb_iter",nb_iter,"_allpoints.pdf",sep="")
    pdf(plotname)
    #par(mar=c(5.1,4.1,4.1,2.1))
    par(mar=c(5.1,5.1,4.1,2.1))
    
    plot(nb_topics_range,AIC2[1,],ann=F,cex.axis=1.5,lwd=2,type="p",col="black",ylim=range(AIC2))
    for (j in 1:nb_real)
      lines(nb_topics_range,AIC2[j,],lwd=2,type="p",col="black")
    
    legend(x="bottomright",legend=c(filename_insert,alpha_insert,paste("delta =",delta),paste("nb realizations =",nb_real),paste("nb iterations =",nb_iter)),col="black",inset=0)
    title("AIC value vs number of topics",cex.main=1.7)
    title(xlab="Number of topics",ylab="AIC value",cex.lab=1.5)
    dev.off()
    #######
    
    if (Rtopicmodels_VEM)
    {
      #######
      plotname = paste("estimated_alpha_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,".pdf",sep="")
      pdf(plotname)
      #par(mar=c(5.1,4.1,4.1,2.1))
      par(mar=c(5.1,5.1,4.1,2.1))
      
      errbar(nb_topics_range,alpha_est0[1,],alpha_est0[1,]+alpha_est0[2,],alpha_est0[1,]-alpha_est0[2,],ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
      
      legend(x="bottomright",legend=c(filename_insert,alpha_insert,paste("delta =",delta),paste("nb realizations =",nb_real),paste("nb iterations =",nb_iter)),col="black",inset=0)
      title("Estimated alpha value vs number of topics",cex.main=1.7)
      title(xlab="Number of topics",ylab="Estimated alpha value",cex.lab=1.5)
      dev.off()
      ########
    }
    
    # llh1=vector(length=20)
    # for (j in 1:20)
    #   llh1[j] = Result@fitted[[j]]@loglikelihood
    # filename = paste("llh_compare_alpha",alpha,"_delta",delta,"_topics",10,"_iter100-10,000.pdf")
    # pdf(filename)
    # plot(seq(1,100,1)*100,Result@fitted[[1]]@logLiks)
    # dev.off()
    
    ##### end of mnb_topics condition
  }
  
  #### end of mpar condition
}

if (local)
{
  end.time <- Sys.time()
  setwd(local_subdirname)
  time_file = "Computation_time.txt"
  write(end.time-start.time,time_file)
}


