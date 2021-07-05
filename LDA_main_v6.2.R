#!/usr/bin/env Rscript

# See LDA_parameters_file.R for options and documentation

# This software is run by calling in a terminal
# from the directory containing LDA_main and all the functions its calls:
# Rscript --vanilla LDA_main_v6.0.R /path.to.data LDA_parameters.R

# Alternatively, LDA_main_v6.0.R can be directly sourced from R or Rstudio 
# after commenting out lines 33 to 39 below and setting manually:  
# parameter_file = "LDA_parameters.R"
# local_data_path = "/path.to.data"

# Depending on the options chosen in LDA_parameters.R,
# this code calls the following functions:
# LDA_realization_comparison_fun.R
# LDA_abiotic_variables_fun.R
# LDA_best_real_fun.R
# LDA_spatial_maps_fun.R
# LDA_detecting_missing_samples_fun.R
# LDA_filling_missing_samples_fun.R
# LDA_testdata_fun.R
# or their "lowMemory" counterparts if the low_memory option is chosen.

# Code written by Guilhem Sommeria-Klein

######################
# Loading parameters #
######################

# This version asks for data_path instead of code_insert as argument and seeks parameter_file there
#######################################
args = commandArgs(trailingOnly = T)

local_data_path = args[1]
parameter_file = args[2]
if (!is.null(args[3]))
  process_index = as.numeric(args[3])+1
#########################################

code_insert = getwd()

if (strsplit(local_data_path,split="",fixed=T)[[1]][1] == ".")
  local_data_path = paste0(getwd(),paste0(strsplit(local_data_path,split="",fixed=T)[[1]][-1],collapse = ""))

source(paste0(local_data_path,"/",parameter_file))  

if (strsplit(abiotic_data_path,split="",fixed=T)[[1]][1] == ".")
  abiotic_data_path = paste0(getwd(),paste0(strsplit(abiotic_data_path,split="",fixed=T)[[1]][-1],collapse = ""))

local_computation = 1
if (cluster_existingresult || local_existingresult)
  local_computation = 0

if (mpar && !local_computation)
  stop("Loading existing results is only possible for mpar = 0.")

if (local_computation)
  start.time <- Sys.time()

mpar_range = length(nb_topics)
if (mpar)
{
  if (mpar_range == 1)
    stop("nb_topics should have length longer that one for mpar = 1.")
  nb_topics_range = nb_topics
  
  if (ValleExtendedLDA_Gibbs)
    stop("Valle's presence-absence LDA cannot be run for mpar = 1.")
    
  if (length(alpha)>1)
  {
    if (all(alpha == 50/nb_topics_range))
      alpha_insert = "_alpha50:K"
    else if (all(alpha == 0.05/nb_topics_range))
      alpha_insert = "_alpha0.05:K"
  } else
    alpha_insert = paste0("_alpha",alpha)

  random_folds_insert = ""
  partitioning_insert = ""
  if (elbow_aic)
  {
    model_selection_insert = "_elbow-AIC"
    real_fold_insert = "_nb_real"
  } else if (cross_validation)
  {
    model_selection_insert = "_post-predictive-cross-valid"
    real_fold_insert = "_fold_size"
    if (!partitioning)
      partitioning_insert = paste0("_overlapping-folds",nb_real,"real")
    if (!random_folds)
      random_folds_insert = "_blocks"
  }
} else
{
  if (alpha == 50/nb_topics)
    alpha_insert = "_alpha50:K"
  else if (alpha == 0.05/nb_topics)
    alpha_insert = "_alpha0.05:K"
  else
    alpha_insert = paste0("_alpha",alpha)
}

if (best != 2 && (Rtopicmodels_Gibbs || ValleExtendedLDA_Gibbs))
  warning("Gibbs samplings should only run with best = 2 (approximate maximum a posteriori).")

if (Rtopicmodels_Gibbs)
{
  filename_insert = "Rtopicmodels_LDA_Gibbs"
  if (!mpar)
  {
    if (best == 1)
    {
      best_insert = paste0("_thin",thin,"_burnin",burnin)
    } else if (best == 2)
    {
      best_insert = paste0("_meanPosteriorDistributedLlh_thin",thin,"_burnin",burnin)
    } else
      best_insert = ""
  } else
  {
    if (best != 0)
    {
      best_insert = paste0("_thin",thin,"_burnin",burnin)
    } else
      best_insert = ""
  }
} else if (Rtopicmodels_VEM) 
{
  filename_insert = "Rtopicmodels_LDA_VEM"
} else if (ValleExtendedLDA_Gibbs)
{
  filename_insert = "ValleExtended_LDA_Gibbs"
  best_insert = ""
} else if (hdp)
  filename_insert = "HDP-faster"

if (randomInit)
{
  init_insert = ""
} else if (seededInit)
{
  init_insert = "_seededInit"
}

if (occurrence)
{
  occurrence_insert = "_occurrence"
} else 
{
  if (ValleExtendedLDA_Gibbs)
    stop("Valle's presence-absence LDA cannot be run for abundance data.")
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
  no_rare_insert = paste0("_noRareOTU_",no_rare_factor,"factor")
} else 
{
  no_rare_insert = ""
}

if (maxmatching)
{
  matching_insert = "maxmatching"
} else if (greedymatching)
{
  matching_insert = "greedymatching"
}

if (samplewise)
{
  MOTU_sample_insert = "samplewise"
} else if (MOTUwise)
  MOTU_sample_insert = "MOTUwise"

################
# Loading data #
################

# setwd(local_data_path)

cat("Loading original data ...\n")

# The data may be provided under the form of text file or a Rdata file.
# The Rdata format may be convenient if the data are outputted from another R script,
# e.g. if they result from splitting a taxonomic group based on taxonomic assignments.
if (testdata)
{
  load(paste0(local_data_path,"/data2m_testdata.Rdata"))
  load(paste0(local_data_path,"/true_documents.Rdata"))
  true_documents = true_norm_documents
  # true_topic_compo not available
} else if (file.exists(paste0(local_data_path,"/data2m.Rdata")))
{
  load(paste0(local_data_path,"/data2m.Rdata"))
} else 
{
  data_file = paste0(local_data_path,"/",data_file)
  if (low_memory)
  {
    if (strsplit(data_file,split=".",fixed=T)[[1]][length(strsplit(data_file,split=".",fixed=T)[[1]])] == "csv")
      data2m = fread(data_file, colClasses="numeric", sep=",", header = T, row.names = 1)
    else
      data2m = fread(data_file, colClasses="numeric", sep="\t", header = T, row.names = 1)
  } else
  {
    if (strsplit(data_file,split=".",fixed=T)[[1]][length(strsplit(data_file,split=".",fixed=T)[[1]])] == "csv")
      data2m = read.table(data_file, colClasses="numeric", sep=",", header = T, row.names = 1)
    else 
      data2m = read.table(data_file, colClasses="numeric", sep="\t", header = T, row.names = 1)
  }
}

if (file.exists(paste0(local_data_path,"/taxo_ref.Rdata")))
{
  load(paste0(local_data_path,"/taxo_ref.Rdata"))
} else if (!testdata && exists("taxo_ref_file"))
{
  taxo_ref_file = paste0(local_data_path,"/",taxo_ref_file)
  if (strsplit(taxo_ref_file,split=".",fixed=T)[[1]][length(strsplit(taxo_ref_file,split=".",fixed=T)[[1]])] == "csv")
    taxo_ref = read.table(taxo_ref_file, colClasses="vector", sep=",", header = T, row.names = 1)
  else
    taxo_ref = read.table(taxo_ref_file, colClasses="vector", sep="\t", header = T, row.names = 1)
}

if (geographic)
{
  if (file.exists(paste0(local_data_path,"/coord.Rdata")))
  {
    load(paste0(local_data_path,"/coord.Rdata"))
  } else
  {
    coord_file = paste0(local_data_path,"/",coord_file)
    # coord_file contains the eastings (first column) and northings (second column) of all samples, without titles
    if (strsplit(coord_file,split=".",fixed=T)[[1]][length(strsplit(coord_file,split=".",fixed=T)[[1]])] == "csv")
      coord = read.table(coord_file, colClasses="numeric", sep=";", col.names = c("x","y"), header = T, row.names = 1)
    else
      coord = read.table(coord_file, colClasses="numeric", sep="\t", col.names = c("x","y"), header = T, row.names = 1)
  }
  # the coordinates are multiplied by the grain for spatial representation with kriging between samples
  if (grid)
    coord = coord*grain
  if (nrow(coord) != ncol(data2m))
    stop("Mismatch between spatial coordinates and biotic data.")
} else if (grid && !geographic)
{
  if (byRow)
    # Samples are listed row by row:
    coord = expand.grid(x = grain*seq(1,ncol_grid,1), y = grain*seq(1,nrow_grid,1))
  else
    # Samples are listed column by column:
    coord = expand.grid(y = grain*seq(1,nrow_grid,1), x = grain*seq(1,ncol_grid,1))[,c(2,1)]
} else 
  coord = NULL

missing = 0
if (grid && !geographic && nrow_grid*ncol_grid != ncol(data2m))
{
  missing = 1
}
missing0 = missing  

########################
# Creating directories #
########################

# Creating a subdirectory specific to the current set of paremeters in the data directory. 
if (grid && !geographic && missing && !filling)
{
  local_dirname = paste0(local_data_path,"/Unfilled/")
  if (cluster_existingresult)
    cluster_dirname = paste(cluster_data_path,"/Unfilled/")
} else 
{
  local_dirname = paste0(local_data_path,"/")
  if (cluster_existingresult)
    cluster_dirname = paste0(cluster_data_path,"/")
}

if (mpar)
{
  if (Rtopicmodels_Gibbs)
  {
    local_subdirname = paste0(local_dirname,filename_insert,alpha_insert,"_delta",format(delta,digits=3),"_nb_topics",nb_topics_range[1],
                              "-",nb_topics_range[length(nb_topics_range)],model_selection_insert,partitioning_insert,random_folds_insert,real_fold_insert,ifelse(cross_validation,fold_size,nb_real),"_nb_iter",nb_iter,best_insert,occurrence_insert,remove_single_sites_insert,no_rare_insert,init_insert,"/")
    if (cluster_existingresult)
      cluster_subdirname = paste0(cluster_dirname,filename_insert,alpha_insert,"_delta",format(delta,digits=3),"_nb_topics",nb_topics_range[1],
                                  "-",nb_topics_range[length(nb_topics_range)],model_selection_insert,partitioning_insert,random_folds_insert,real_fold_insert,ifelse(cross_validation,fold_size,nb_real),"_nb_iter",nb_iter,best_insert,occurrence_insert,remove_single_sites_insert,no_rare_insert,init_insert,"/")
    filename = paste0(filename_insert,alpha_insert,"_delta",format(delta,digits=3),"_nb_topics",nb_topics_range[1],
                      "-",nb_topics_range[length(nb_topics_range)],model_selection_insert,partitioning_insert,random_folds_insert,real_fold_insert,ifelse(cross_validation,fold_size,nb_real),"_nb_iter",nb_iter,best_insert,occurrence_insert,remove_single_sites_insert,no_rare_insert,init_insert,".Rdata")  
  } else if (Rtopicmodels_VEM)
  {
    local_subdirname = paste0(local_dirname,filename_insert,alpha_insert,"_nb_topics",nb_topics_range[1],
                              "-",nb_topics_range[length(nb_topics_range)],model_selection_insert,partitioning_insert,random_folds_insert,real_fold_insert,ifelse(cross_validation,fold_size,nb_real),"_em_tol",em_tol,"_var_tol",var_tol,occurrence_insert,remove_single_sites_insert,no_rare_insert,init_insert,"/")
    if (cluster_existingresult)
      cluster_subdirname = paste0(cluster_dirname,filename_insert,alpha_insert,"_nb_topics",nb_topics_range[1],
                                  "-",nb_topics_range[length(nb_topics_range)],model_selection_insert,partitioning_insert,random_folds_insert,real_fold_insert,ifelse(cross_validation,fold_size,nb_real),"_em_tol",em_tol,"_var_tol",var_tol,occurrence_insert,remove_single_sites_insert,no_rare_insert,init_insert,"/")
    filename = paste0(filename_insert,alpha_insert,"_nb_topics",nb_topics_range[1],
                      "-",nb_topics_range[length(nb_topics_range)],model_selection_insert,partitioning_insert,random_folds_insert,real_fold_insert,ifelse(cross_validation,fold_size,nb_real),"_em_tol",em_tol,"_var_tol",var_tol,occurrence_insert,remove_single_sites_insert,no_rare_insert,init_insert,".Rdata")
  }
} else if (!mpar)
{
  if (Rtopicmodels_Gibbs)
  {
    local_subdirname = paste0(local_dirname,filename_insert,alpha_insert,"_delta",format(delta,digits=3),"_nb_topics",nb_topics,"_nb_iter",nb_iter,"_nb_real",nb_real,best_insert,occurrence_insert,remove_single_sites_insert,no_rare_insert,init_insert,"/")
    if (cluster_existingresult)
      cluster_subdirname = paste0(cluster_dirname,filename_insert,alpha_insert,"_delta",format(delta,digits=3),"_nb_topics",nb_topics,"_nb_iter",nb_iter,"_nb_real",nb_real,best_insert,occurrence_insert,remove_single_sites_insert,no_rare_insert,init_insert,"/")
    filename = paste0(filename_insert,alpha_insert,"_delta",format(delta,digits=3),"_nb_topics",nb_topics,"_nb_iter",nb_iter,"_nb_real",nb_real,best_insert,occurrence_insert,remove_single_sites_insert,no_rare_insert,".Rdata")
  } else if (Rtopicmodels_VEM)
  {
    local_subdirname = paste0(local_dirname,filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best_keep",occurrence_insert,remove_single_sites_insert,no_rare_insert,init_insert,"/")
    if (cluster_existingresult)
      cluster_subdirname = paste0(cluster_dirname,filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best_keep",occurrence_insert,remove_single_sites_insert,no_rare_insert,init_insert,"/")
    filename = paste0(filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best_keep",occurrence_insert,remove_single_sites_insert,no_rare_insert,init_insert,".Rdata")
  } else if (ValleExtendedLDA_Gibbs)
  {
    local_subdirname = paste0(local_dirname,filename_insert,alpha_insert,"_delta",format(delta,digits=3),"_nb_topics",nb_topics,"_nb_iter",nb_iter,"_nb_real",nb_real,best_insert,remove_single_sites_insert,no_rare_insert,"/")
    if (cluster_existingresult)
      cluster_subdirname = paste0(cluster_dirname,filename_insert,alpha_insert,"_delta",format(delta,digits=3),"_nb_topics",nb_topics,"_nb_iter",nb_iter,"_nb_real",nb_real,best_insert,remove_single_sites_insert,no_rare_insert,"/")
    filename = paste0(filename_insert,alpha_insert,"_delta",format(delta,digits=3),"_nb_topics",nb_topics,"_nb_iter",nb_iter,"_nb_real",nb_real,best_insert,remove_single_sites_insert,no_rare_insert,".Rdata") 
  } else if (hdp)
  {
    local_subdirname = paste0(local_dirname,filename_insert,"hdp-faster_alpha",alpha,"_gamma",gamma,"_eta",eta,"_maxiter",maxiter,"_norep",occurrence_insert,remove_single_sites_insert,no_rare_insert,"/")
    if (!file.exists(local_subdirname))
      dir.create(local_subdirname)
  }
}

######################
# Data preprocessing #
######################

cat("Preprocessing original data ...\n")

if (grid && !geographic)
{
  if (missing)
  {
    # Detecting the positions of missing samples only if it has not been previously done
    if (!file.exists(paste0(local_data_path,"/Missing_positions_indices.Rdata")))
    {
      source(paste0(code_insert,"/LDA_detecting_missing_samples_fun.R"))
      if (byRow)
      {
        nb_repeats = nrow_grid
        if (ncol_grid != length(Sample_name_endings))
          stop("Mismatch between grid dimensions and sample_name_endings.")
      } else 
      {
        nb_repeats = ncol_grid
        if (nrow_grid != length(Sample_name_endings))
          stop("Mismatch between grid dimensions and sample_name_endings.")
      } 
      # Missing_positions_indices is a vector of the same length as the number of samples (missing or not) in the grid,
      # taking value 1 for a missing sample and 0 otherwise.
      # Missing_positions_indices describes sample positions in the same order as data2m.
      Missing_positions_indices = LDA_detecting_missing_samples_fun(Sample_name_endings,colnames(data2m),nb_repeats,byRow)
      
      save(Missing_positions_indices, file = "Missing_positions_indices.Rdata")
    } else 
      # Loading the precomputed Missing_positions_indices
      load("Missing_positions_indices.Rdata")
  } else
    Missing_positions_indices = vector(length = nrow_grid*ncol_grid, mode = "numeric")
} else 
  Missing_positions_indices = vector(length = ncol(data2m), mode = "numeric")

if (grid && !geographic && missing && filling)
{
  # Filling missing samples only if it has not been previously done
  if (!file.exists(paste0(local_data_path,"/data2m_filled.Rdata")) &&
      !(cluster_existingresult && file.exists(paste0(cluster_subdirname,"data2m_filled.Rdata"))))
  {
    if (all(Missing_positions_indices==0))
      cat("No missing sample.\n")
    else 
    {
      cat("Filling missing samples ...\n")
      
      source(paste0(code_insert,"/LDA_filling_missing_samples_fun.R"))
      
      data2m = LDA_filling_missing_samples_fun(Missing_positions_indices,data2m,nrow_grid,ncol_grid,byRow)
      # Missing_positions_indices = rep(0,length(Missing_positions_indices))
      save(data2m,file=paste0(local_data_path,"/data2m_filled.Rdata"))
    }
  } else if (!file.exists(paste0(local_data_path,"/data2m_filled.Rdata")) && cluster_existingresult && file.exists(paste0(cluster_subdirname,"data2m_filled.Rdata")))
  {
    # Because filling the missing samples involves random numbers, it should not be performed a second time independently
    # if it has already been performed on the cluster (for cluster_existingresult = 1): 
    command = paste0("cp ",cluster_subdirname,"data2m_filled.Rdata"," ",local_subdirname,"data2m_filled.Rdata")
    system(command, intern=TRUE)
    load(paste0(local_data_path,"/data2m_filled.Rdata"))
  } else if (file.exists("data2m_filled.Rdata"))
    load(paste0(local_data_path,"/data2m_filled.Rdata"))
  missing = 0
}

# Moving to the directory specfific to the current set of parameters, once all data have been loaded:
if (!file.exists(local_dirname))
  dir.create(local_dirname)
setwd(local_dirname) 

if (!file.exists(local_subdirname))
  dir.create(local_subdirname)
setwd(local_subdirname) 

write(paste0(local_subdirname,filename),file="Directory.txt")
write(paste0("0/",mpar_range),file="Progress.txt",ncolumns=1,append=F)

OTUs_to_be_removed = vector(length=nrow(data2m), mode="logical")
Prop_OTU_removed = 0
Prop_reads_removed = 0
# Removing all OTUs which are present in less than nb_occupied_sites_threshold, including OTUs with 0 total abundance counts
for (OTU in 1:nrow(data2m))
  OTUs_to_be_removed[OTU] = (length(which(data2m[OTU,]>0)) < nb_occupied_sites_threshold)
# Removing all OTUs which are present in less than no_rare_factor reads per site on average
if (no_rare)
{
  for (OTU in 1:nrow(data2m))
  {
    if (OTUs_to_be_removed[OTU] == 0)
      OTUs_to_be_removed[OTU] = !(sum(data2m[OTU,])>length(data2m[OTU,])*no_rare_factor)
  }
}
if (any(OTUs_to_be_removed))
{
  Prop_OTU_removed = length(which(OTUs_to_be_removed))/nrow(data2m)
  Prop_reads_removed = sum(data2m[which(OTUs_to_be_removed),])/sum(data2m)
  data2m = data2m[!OTUs_to_be_removed,]
  #rownames(data2m) = seq(1,nrow(data2m),1)
  if (exists("taxo_ref"))
  {
    taxo_ref = taxo_ref[!OTUs_to_be_removed,]
    # taxo_ref = taxo_ref[-(which(OTUs_to_be_removed)+1),]
    #rownames(taxo_ref) = seq(0,nrow(data2m),1)
    #taxo_ref[2:(nrow(data2m)+1),1] = seq(1,nrow(data2m),1)
  }
}

# Writing how many samples and read have been lost to preprocessing in text file:
write(paste("Proportion of OTUs removed before LDA decomposition:",Prop_OTU_removed),file="Removed_OTUs_and_reads.txt")
write(paste("Proportion of reads removed before LDA decomposition:",Prop_reads_removed),file="Removed_OTUs_and_reads.txt",append=T)

if (occurrence)
{
  if (any(data2m != 1 & data2m != 0))
    data2m[data2m > 0] = 1
  else
    occurrence_insert = ""
}

# Dealing with empty samples now that data have been preprocessed:
if (length(which(colSums(data2m)==0)) != 0)
{
  missing = 1
  if (grid && !geographic)
  {
    colSums_data2m = colSums(data2m)
    spatial_colSums_data2m = matrix(nrow=nrow_grid, ncol=ncol_grid, data=0)
    position_shift = 0
    # Representing the samples as a spatial matrix, so as to record the empty ones in Missing_positions_indices:
    if (byRow)
    {
      for (i in 1:nrow_grid)
      {
        for (j in 1:ncol_grid)
        {
          missing_index = (i-1)*ncol_grid+j
          if (Missing_positions_indices[missing_index]==0)
          {
            spatial_colSums_data2m[i,j] = colSums_data2m[missing_index-position_shift]    
          } else if (Missing_positions_indices[missing_index]==1)
          {
            spatial_colSums_data2m[i,j] = NA
            position_shift = position_shift+1
          }
        }
      }
      spatial_colSums_data2m = as.vector(t(spatial_colSums_data2m))
    } else
    {
      for (j in 1:ncol_grid)
      {
        for (i in 1:nrow_grid)
        {
          missing_index = (j-1)*nrow_grid+i
          if (Missing_positions_indices[missing_index]==0)
          {
            spatial_colSums_data2m[i,j] = colSums_data2m[missing_index-position_shift]    
          } else if (Missing_positions_indices[missing_index]==1)
          {
            spatial_colSums_data2m[i,j] = NA
            position_shift = position_shift+1
          }
        }
      }
      spatial_colSums_data2m = as.vector(spatial_colSums_data2m)
    }
    
    # Saving the new missing locations (locations that end up being empty due to data preprocessing) 
    # plus the initial ones as Missing_position_indices0 (for abiotic_comparison):
    Missing_positions_indices0 = Missing_positions_indices
    Missing_positions_indices0[which(spatial_colSums_data2m==0)] = 1
    # Saving only the new missing locations as Missing_positions_indices (for realization comparison):
    Missing_positions_indices = vector(length = length(Missing_positions_indices), mode = "numeric")
    Missing_positions_indices[which(spatial_colSums_data2m==0)] = 1
    # Removing both new and old missing locations from coord0, and only the new ones for coord (for spatial_maps):
    coord1 = coord[-which(Missing_positions_indices==1),]
    coord0 = coord[-which(Missing_positions_indices0==1),]
    coord = coord1
  } else if (geographic)
  {
    coord = coord[-which(colSums(data2m)==0),]
    if (grid)
      # Missing_positions_indices is used in data_abiotic and realization_comparison
      Missing_positions_indices[which(spatial_colSums_data2m==0)] = 1
  }
  # Removing the empty locations from data2m before applying LDA: 
  data2m = data2m[,-which(colSums(data2m)==0)]
} else 
{
  # For grid and !geographic plotting:
  Missing_positions_indices0 = Missing_positions_indices
  Missing_positions_indices = vector(length = length(Missing_positions_indices), mode = "numeric")
  coord0 = coord[-which(Missing_positions_indices0==1),]
}

# Saving useful properties of the original data and removing them from memory if not further needed
nb_terms = nrow(data2m)
nb_doc = ncol(data2m)
sum_data2m = sum(data2m)
# Saving sample names:
colnames_data2m = colnames(data2m)

if (!local_computation)
  remove(data2m)

if (mpar && cross_validation)
{
  # Generating folds of hold-out samples for cross-validation.
  # The folds are generated once for all and reused for each K value.
  if (partitioning)
  {
    # Partitioning the samples into folds of "fold_size" samples.
    # The number of folds ("nb_real") depends on the number of samples and on fold size.
    # If the number of samples cannot be exactly divided by fold size, one fold will contain a smaller or larger number of samples.
    if (fold_size < 1)
      stop("fold_size should be larger than 1 for partitioning = TRUE.")
    nb_real = nb_doc %/% fold_size
    if (nb_doc %% fold_size > fold_size/2)
    {
      split_factor = c(unlist(lapply(1:nb_real,rep,fold_size)),rep(nb_real+1,nb_doc %% fold_size))
      nb_real = nb_real+1
    } else
    {
      split_factor = c(unlist(lapply(1:nb_real,rep,fold_size)),rep(nb_real,nb_doc %% fold_size))
    }
    cat("Generating",nb_real,"folds of",fold_size,"samples out of",nb_doc,"samples by exact partitioning.\n")
    # If random_folds=F, samples are partitioned into folds following their original order in the data ("block cross-validation"),
    # otherwise they are shuffled.
    if (random_folds)
      split_factor = sample(split_factor,length(split_factor))
    # Generating the folds:
    folds = split(1:nb_doc,split_factor)
    if (nb_doc != length(split_factor) || nb_real != length(folds))
      stop("Error while generating folds.")
  } else
  {
    # Generating nb_real overlapping folds by independently partitioning the samples in two for each realization.
    if (fold_size > 0.5)
      stop("fold_size should be lower than 0.5 for partitioning = FALSE.")
    abs_fold_size = round(fold_size*nb_doc)
    cat("Generating",nb_real,"independent folds of",abs_fold_size,"samples out of",nb_doc,"samples.\n")
    folds = list()
    # If random_folds=F, samples are partitioned into folds following their original order in the data ("block cross-validation"),
    # otherwise they are shuffled.
    if (random_folds)
    {
      # Picking abs_fold_size randomly chosen samples for the fold.
      for (r in 1:nb_real)
        folds[[r]] = sample(1:nb_doc,abs_fold_size)
    } else
    {
      # Picking abs_fold_size contiguous samples for the fold, starting from a random index for each realization.
      for (r in 1:nb_real)
      {
        start_fold = sample(1:nb_doc,1)
        if (start_fold + abs_fold_size - 1 > nb_doc)
        {
          # Leaving the samples (abs_fold_size - (nb_doc - start_fold)):(start_fold-1) out of the fold:
          folds[[r]] = c(1:(abs_fold_size - (nb_doc - start_fold + 1)),
                         start_fold:nb_doc)
        } else 
        {
          folds[[r]] = start_fold:(start_fold + abs_fold_size - 1)
        }
      }
    }
  }
}

####################################
# Loading saved LDA decompositions # 
####################################

if (cluster_existingresult)
{  
  # cat("Loading data from cluster ...\n")
  # command = paste0("cp ",cluster_subdirname,filename," ",local_subdirname,filename)
  # system(command, intern=TRUE)
  cat("Loading LDA decomposition from file:\n",filename,"\n")
  load(paste0(cluster_subdirname,filename))
  if (Rtopicmodels_Gibbs && best == 2 && !mpar)
    sampled_logLiks = readRDS(paste0(cluster_subdirname,"posterior_sampled_logLiks.rds"))
} else if (local_existingresult)
{
  cat("Loading LDA decomposition from file:\n",filename,"\n")
  load(filename)
  if (Rtopicmodels_Gibbs && best == 2 && !mpar)
    sampled_logLiks = readRDS("posterior_sampled_logLiks.rds")
} else if (local_computation)
  cat("Results will be saved in file:\n",filename,"\n")

##################################
# Loop over the number of topics #
##################################

if (mpar)
{
  if (local_computation)
  {
    # Result_mpar = list()
    if (cross_validation)
      perplexity_mpar = matrix(nrow = nb_real, ncol = mpar_range, dimnames = list(paste0("real",1:nb_real),paste0(nb_topics_range,"assemblages")), data = 0)
    if (Rtopicmodels_Gibbs && best == 2 && elbow_aic)
      sampled_logLiks_mpar = list()
  }
  
  if (elbow_aic)
  {
    LLH_final0 = matrix(nrow=2,ncol=mpar_range,data=0)
    LLH_final = matrix(nrow=nb_real,ncol=mpar_range,data=0)
    # perplexity0 = matrix(nrow=2,ncol=mpar_range,data=0)
    # perplexity = matrix(nrow=nb_real,ncol=mpar_range,data=0)
    AIC0 = matrix(nrow=2,ncol=mpar_range,data=0)
    AICc0 = matrix(nrow=2,ncol=mpar_range,data=0)
    AIC = matrix(nrow=nb_real,ncol=mpar_range,data=0)
    AICc = matrix(nrow=nb_real,ncol=mpar_range,data=0)
    if (Rtopicmodels_VEM) 
    {
      alpha_est1 = matrix(nrow=nb_real,ncol=mpar_range,data=0)
      alpha_est0 = matrix(nrow=2,ncol=mpar_range,data=0)
    }
  }
}

# Start of the loop over the number of topics if mpar = 1
# (otherwise mpar_range has length 1)
for (par_index in 1:mpar_range)
{
  if (mpar)
  {
    cat("Parameter value",par_index,"\n")  
    nb_topics = nb_topics_range[par_index]  
  }
  
  ####################################################
  # Computing LDA decomposition                      #
  # using package topicmodels (Grün and Hornik 2011) #
  ####################################################
  
  if (local_computation)
  {
    if (Rtopicmodels_Gibbs)
    {
      # Topics are initalized at random, using seeds chosen so as to never be twice the same
      SEED = as.integer(Sys.time()) - (1:nb_real)*10^7
      # Inferred quantities:
      # beta : dimensions #topics x #terms
      # gamma : dimensions #documents x #topics
      # if best = 1, the chain is sampled every 'thin' iterations and the highest-llh sample is kept;
      # otherwise only the last of the 'nb_iter' iterations is kept
      if (!mpar || (mpar && elbow_aic))
      {
        if (!mpar)
          Result = list()
        if (best == 2)
          sampled_logLiks = matrix(nrow = nb_real, ncol = nb_iter/thin, data = 0)
        for (j in 1:nb_real)
        {
          control_LDA_Gibbs = list(alpha=ifelse(length(alpha)>1,alpha[par_index],alpha), estimate.beta=TRUE,
                                   verbose = 100, prefix = tempfile(), save = 0, keep = llh_keep,
                                   seed = SEED[j], nstart = 1, best = 0,
                                   delta = delta,
                                   iter = nb_iter, burnin = ifelse(best,burnin,0), thin = ifelse(best,thin,nb_iter))
          Result_j = topicmodels::LDA(x=t(data2m),k=nb_topics,method = "Gibbs",control=control_LDA_Gibbs,model=NULL)
          
          # if (mpar && par_index == 1)
          #   nb_words = Result_j@fitted[[1]]@n
          
          if (best == 2)
          {
            sampled_logLiks[j,] = unlist(lapply(Result_j@fitted, function(g) g@loglikelihood))
            if (mpar)
              llh = mean(sampled_logLiks[j,])
            else if (!mpar)
              Result[[j]] = Result_j@fitted[[sort.int(abs(sampled_logLiks[j,] - mean(sampled_logLiks[j,])),index.return = T)$ix[1]]]
          } else if (best == 1)
          {
            best_model_index = sort.int(unlist(lapply(Result_j@fitted, function(g) g@loglikelihood)),decreasing=T,index.return=T)$ix[1]
            if (mpar)
              llh = sum(Result_j@fitted[[best_model_index]]@loglikelihood)
            else if (!mpar)
              Result[[j]] = Result_j@fitted[[best_model_index]]
          } else if (!best)
          {
            if (mpar)
              llh = sum(Result_j@fitted[[1]]@loglikelihood)
            else if (!mpar)
              Result[[j]] = Result_j@fitted[[1]]
          }
          
          if (mpar)
          {
            LLH_final[j,par_index] = llh
            AIC[j,par_index] = 2*(nb_topics*(nb_terms-1) + 1 - llh)
            AICc[j,par_index] = 2*((nb_topics*(nb_terms-1) + 1)*(1+1/nb_doc) - llh)
          }
          
          if (!mpar && best)
            Result[[j]]@logLiks = unlist(lapply(Result_j@fitted, function(g) g@logLiks))
        }
        
        if (mpar && elbow_aic && best == 2)
        {
          sampled_logLiks_mpar[[par_index]] = sampled_logLiks
          # Intermediate results are saved for each parameter value:
          saveRDS(sampled_logLiks_mpar,file="posterior_sampled_logLiks.rds")
        }
      } else if (mpar && cross_validation)
      {
        for (j in 1:nb_real)
        {
          control_LDA_Gibbs = list(alpha=ifelse(length(alpha)>1,alpha[par_index],alpha), estimate.beta=TRUE,
                                   verbose = 100, prefix = tempfile(), save = 0, keep = llh_keep,
                                   seed = SEED[j], nstart = 1, best = 0,
                                   delta = delta,
                                   iter = nb_iter, burnin = ifelse(best,burnin,0), thin = ifelse(best,thin,nb_iter))
          
          Result = topicmodels::LDA(x=t(data2m[,-folds[[j]]]),k=nb_topics,method = "Gibbs",control=control_LDA_Gibbs,model=NULL)
          
          control_LDA_Gibbs = list(alpha=ifelse(length(alpha)>1,alpha[par_index],alpha), estimate.beta=F,
                                   verbose = 100, prefix = tempfile(), save = 0, keep = llh_keep,
                                   seed = SEED[j], nstart = 1, best = 1,
                                   delta = delta,
                                   iter = nb_iter, burnin = ifelse(best,burnin,0), thin = nb_iter)
          
          perplexity_mpar[j,par_index] = topicmodels::perplexity(object = Result, newdata = t(data2m[,folds[[j]]]), control = control_LDA_Gibbs)
        }
        
        # Intermediate results are saved for each parameter value:
        save(nb_topics_range,perplexity_mpar,folds,file="perplexity.Rdata")
        # Progress.text gives information about the progress of the computation:
        write(paste0(par_index,"/",mpar_range),file="Progress.txt",ncolumns=1,append=T)
      }
    } else if (Rtopicmodels_VEM)
    {
      # if (mpar && (cluster_existingresult || local_existingresult))
      # {
      #   Result = list()
      #   for (j in 1:nb_real)
      #     Result[[j]] = Result_mpar[[(par_index-1)*nb_real+j]]
      # if (local_computation)
      # {
      SEED = as.integer(Sys.time()) - (1:(2*nb_real))*10^7
      # Inferred quantities:
      # beta : dimensions #topics x #terms
      # gamma : dimensions #documents x #topics
      # alpha : scalar
      if (randomInit)
        initialize_arg = "random"
      else if (seededInit)
        initialize_arg = "seeded"
      
      if (!mpar || (mpar && elbow_aic))
      {
        control_LDA_VEM = list(estimate.alpha=TRUE, alpha=alpha, estimate.beta=TRUE,
                               verbose = 1, prefix = tempfile(), save = 0, keep = 1,
                               seed = SEED[1:nb_real], nstart = nb_real, best = 0,
                               var = list(iter.max = 500, tol = var_tol),
                               em = list(iter.max = 1000, tol = em_tol),
                               initialize = initialize_arg)
        
        Result = topicmodels::LDA(x=t(data2m),k=nb_topics,method = "VEM",control=control_LDA_VEM,model=NULL)
        
      } else if (mpar && cross_validation)
      {
        Perplexity = vector(length = nb_real, mode = "numeric")
        Result = list()
        for (j in 1:nb_real)
        {
          control_LDA_VEM = list(estimate.alpha=TRUE, alpha=alpha, estimate.beta=TRUE,
                                 verbose = 1, prefix = tempfile(), save = 0, keep = 1,
                                 seed = SEED[j], nstart = 1, best = 0,
                                 var = list(iter.max = 500, tol = var_tol),
                                 em = list(iter.max = 1000, tol = em_tol),
                                 initialize = initialize_arg)
          
          Result[[j]] = topicmodels::LDA(x=t(data2m[,-folds[[j]]]),k=nb_topics,method = "VEM",control=control_LDA_VEM,model=NULL)[[1]]
          
          control_LDA_VEM = list(estimate.alpha=TRUE, alpha=alpha, estimate.beta=F,
                                 verbose = 1, prefix = tempfile(), save = 0, keep = 1,
                                 seed = SEED[j+nb_real], nstart = 1, best = 0,
                                 var = list(iter.max = 500, tol = var_tol),
                                 em = list(iter.max = 1000, tol = em_tol),
                                 initialize = initialize_arg)
          
          Perplexity[j] = topicmodels::perplexity(object = Result[[j]], newdata = t(data2m[,folds[[j]]]), control = control_LDA_VEM)
        }
      }
      
      # Removing realizations with loglikelihood equal to 0 and replacing them with new realizations. 
      # Realizations with loglikelihood equal to 0 correspond to degenerate decompositions where all topics have equal proportions in all samples.
      discarded_real = 1
      while (length(discarded_real) > 0)
      {
        discarded_real = NULL
        for (j in 1:nb_real)
        { 
          if (max(Result[[j]]@logLiks) >= 0)
          {
            if (length(discarded_real) > 0)
              discarded_real = c(discarded_real,j)
            else
              discarded_real = j
            cat("Positive loglikelihood value: realization",j,"computed again\n")
          }
        }
        
        if (length(discarded_real) > 0)
        { 
          SEED = vector(length=2*length(discarded_real),mode="integer")
          jj = 0
          for (j in discarded_real)
          {
            jj = jj+1
            SEED[jj] = as.integer(Sys.time()) - j*10^7
            SEED[jj+length(discarded_real)] = as.integer(Sys.time()) - (j+nb_real)*10^7
          }
          
          if (!mpar || (mpar && elbow_aic))
          {
            control_LDA_VEM = list(estimate.alpha=TRUE, alpha=alpha, estimate.beta=TRUE,
                                   verbose = 1, prefix = tempfile(), save = 0, keep = 1,
                                   seed = SEED[1:length(discarded_real)], nstart = length(discarded_real), best = 0,
                                   var = list(iter.max = 500, tol = var_tol),
                                   em = list(iter.max = 1000, tol = em_tol),
                                   initialize = initialize_arg)
            
            Result_add = topicmodels::LDA(x=t(data2m),k=nb_topics,method = "VEM",control=control_LDA_VEM,model=NULL)
          } else if (mpar && cross_validation)
          {
            Result_add = list()
            Perplexity_add = vector(length = length(discarded_real), mode = "numeric")
            jj = 0
            for (j in discarded_real)
            {
              jj = jj+1
              control_LDA_VEM = list(estimate.alpha=TRUE, alpha=alpha, estimate.beta=TRUE,
                                     verbose = 1, prefix = tempfile(), save = 0, keep = 1,
                                     seed = SEED[jj], nstart = 1, best = 0,
                                     var = list(iter.max = 500, tol = var_tol),
                                     em = list(iter.max = 1000, tol = em_tol),
                                     initialize = initialize_arg)
              
              Result_add[[jj]] = topicmodels::LDA(x=t(data2m[,-folds[[j]]]),k=nb_topics,method = "VEM",control=control_LDA_VEM,model=NULL)[[1]]
              
              control_LDA_VEM = list(estimate.alpha=TRUE, alpha=alpha, estimate.beta=F,
                                     verbose = 1, prefix = tempfile(), save = 0, keep = 1,
                                     seed = SEED[jj+length(discarded_real)], nstart = 1, best = 0,
                                     var = list(iter.max = 500, tol = var_tol),
                                     em = list(iter.max = 1000, tol = em_tol),
                                     initialize = initialize_arg)
              
              Perplexity_add[jj] = topicmodels::perplexity(object = Result_add[[jj]], newdata = t(data2m[,folds[[j]]]), control = control_LDA_VEM)
            }
          }
          
          jj = 0
          for (j in discarded_real)
          {
            jj = jj+1
            Result[[j]] = Result_add[[jj]]
            if (mpar && cross_validation)
              Perplexity[j] = Perplexity_add[jj]
          }
        }
      }
      
      if (mpar)
      {
        for (j in 1:nb_real)
        {
          # Result_mpar[[((par_index-1)*nb_real+j)]] = Result[[j]]
          if (elbow_aic)
          {
            for (j in 1:nb_real)
            {
              llh = sum(Result[[j]]@loglikelihood)
              alpha_est1[j,par_index] = Result[[j]]@alpha
              LLH_final[j,par_index] = llh
              AIC[j,par_index] = 2*(nb_topics*(nb_terms-1) + 1 - llh)
              AICc[j,par_index] = 2*((nb_topics*(nb_terms-1) + 1)*(1+1/nb_doc) - llh)
            }
          } else if (cross_validation)
            perplexity_mpar[j,par_index] = Perplexity[j]
        }
        
        if (cross_validation)
        {
          save(nb_topics_range,perplexity_mpar,folds,file="perplexity.Rdata")
          write(paste0(par_index,"/",mpar_range),file="Progress.txt",ncolumns=1,append=T)
        }
      }
    } else if (ValleExtendedLDA_Gibbs)
    {
      #get functions
      source(paste0(code_insert,'/Valle_lda.presence.absence_fun.R'))
      source(paste0(code_insert,'/Valle_get.theta_fun.R'))
      sourceCpp(paste0(code_insert,'/aux1.cpp'))
      
      if (!mpar)
      {
        Result = list()
        # if (best == 2)
        # {
        #   sampled_logLiks = matrix(nrow = nb_real, ncol = (nb_iter - burnin)/thin, data = 0)
        #   evidence = vector(length = nb_real, mode = "numeric")
        # }
        for (j in 1:nb_real)
        {
          # Topics are initalized at random, using seeds chosen so as to never be twice the same
          set.seed(as.integer(Sys.time()) - j*10^7)
          #run Gibbs sampler
          Result[[j]] = lda.presence.absence(dat = t(data2m),id = 1:ncol(data2m), ncomm = nb_topics,
                                             a.phi = delta, b.phi = delta, gamma = alpha, ngibbs = nb_iter)
        }
      }
    }
  }
    
  ######################
  # Extracting outputs #
  ######################
  
  if (mpar && elbow_aic)
  {
    LLH_final0[1,par_index] = mean(LLH_final[,par_index])
    LLH_final0[2,par_index] = sd(LLH_final[,par_index])
    # perplexity0[1,par_index] = mean(perplexity[,par_index])
    # perplexity0[2,par_index] = sd(perplexity[,par_index])
    AIC0[1,par_index] = 2*(nb_topics*(nb_terms-1) + 1 - LLH_final0[1,par_index])
    AICc0[1,par_index] = 2*((nb_topics*(nb_terms-1) + 1)*(1+1/nb_doc) - LLH_final0[1,par_index])
    AIC0[2,par_index] = sd(AIC[,par_index])
    AICc0[2,par_index] = sd(AICc[,par_index])
    if (Rtopicmodels_VEM)
    {
      alpha_est0[1,par_index] = mean(alpha_est1[,par_index])
      alpha_est0[2,par_index] = sd(alpha_est1[,par_index])
      write(paste("nb topics =",nb_topics_range[par_index],"- mean alpha =",alpha_est0[1,par_index]),file="estimated_alpha.txt",ncolumns=1,append=ifelse(par_index == 1,F,T))
    }
    
    if (Rtopicmodels_VEM)
      save(nb_topics_range,LLH_final0,LLH_final,alpha_est1,alpha_est0,AIC0,AICc0,AIC,AICc,file="AIC-llh.Rdata")
    else 
      save(nb_topics_range,LLH_final0,LLH_final,AIC0,AICc0,AIC,AICc,file="AIC-llh.Rdata")
      
    # Progress.text gives information about the progress of the computation:
    write(paste0(par_index,"/",mpar_range),file="Progress.txt",ncolumns=1,append=T)
  } else if (!mpar)
  {  
    LLH_final_real1 = vector(length=nb_real,mode="numeric")
    for (j in 1:nb_real)
      LLH_final_real1[j] = sum(Result[[j]]@loglikelihood)
  
    if (ValleExtendedLDA_Gibbs)
    {
      for (j in 1:nb_real)
        LLH_final_real1[j] = Result[[j]]$llk[nb_iter]
    }
    # Realizations are ordered from the highest likelihood (best realization) to the lowest
    if (Rtopicmodels_Gibbs && best == 2)
    {
      Ordered_realizations = sort.int(abs(LLH_final_real1 - mean(LLH_final_real1)),index.return = T)
      HighestLlh_ordered_realizations = vector(length = nb_real,mode = "numeric")
      for (i in seq_along(Ordered_realizations$ix))
        HighestLlh_ordered_realizations[i] = which(sort.int(LLH_final_real1,decreasing=T,index.return=T)$ix == Ordered_realizations$ix[i])
      saveRDS(HighestLlh_ordered_realizations,file="LlhRank_of_ordered_realizations.rds")
    } else 
      Ordered_realizations = sort.int(LLH_final_real1,decreasing=T,index.return=T)
    saveRDS(Ordered_realizations,file="Ordered_realizations.rds")
    Akaike_weights_llh = exp(LLH_final_real1 - LLH_final_real1[Ordered_realizations$ix[1]])/sum(exp(LLH_final_real1 - LLH_final_real1[Ordered_realizations$ix[1]]))
  }
}
# End of the loop over the number of topics

##############################
# Saving computation results #
##############################

if (local_computation)
{
  cat("Saving result ...\n")
  if (!mpar)
  {
    save(Result,file=filename)
    if (Rtopicmodels_Gibbs && best == 2)
      saveRDS(sampled_logLiks,file="posterior_sampled_logLiks.rds")
  } 
  remove(data2m)
}

################################################
# ############################################ #
# ######## POSTPROCESSING AND FIGURES ######## #
# ############################################ #
################################################

cat("\nStarting postprocessing ...\n")

if (!mpar)
{ 
  length_selected_real = length(Selected_real)
  alpha_est_mean = 0
  
  #################################################
  # Defining variables for realization comparison #
  #################################################
  
  if (realization_comparison)
  {
    if (length_selected_real == 1)
      stop("Needs more than one realization to compare realizations")
    
    # if (testdata || !grid && !geographic)
    #   Correlation_samplewise = vector(mode="list",length=length_selected_real-1)
    
    if (low_memory)
      source(paste0(code_insert,"/LDA_realization_comparison_lowMemory_fun_v3.R"))
    else
      source(paste0(code_insert,"/LDA_realization_comparison_fun_v3.R"))
    
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
    
    Correlation_allRealPairs = matrix(nrow=length_selected_real,ncol=length_selected_real,data=0)
    Correlation_allRealPairs_byTopic = list()
    for (k in 1:nb_topics)
      Correlation_allRealPairs_byTopic[[k]] = matrix(nrow=length_selected_real,ncol=length_selected_real,data=NA)
    
    # Topic_correspondence_to_best_real = list()
    # Topic_correspondence_to_best_real[[1]] = seq(1,nb_topics,1)
    
    llh_differences_allRealPairs = matrix(nrow=length_selected_real,ncol=length_selected_real,data=NA)
    # perplexity_differences_allRealPairs = matrix(nrow=length_selected_real,ncol=length_selected_real,data=NA)
    
    documents_allreal = vector(mode="list",length=length_selected_real)
    KL_documents_allreal = vector(mode="list",length=length_selected_real)
    topic_compo_allreal = vector(mode="list",length=length_selected_real)
    KL_topic_compo_allreal = vector(mode="list",length=length_selected_real)
    #sort_normal_topic_allreal = vector(mode="list",length=length_selected_real)
  }
  
  cat("Starting the loop over realizations ...\n")
  
  ###########################################################
  # Loop over realizations ordered by decreasing likelihood #
  ###########################################################
  
  for (j_select in 1:length_selected_real)
  {
    cat("##############\nRealization",j_select,"\n")
    
    ###########################################################
    # Defining variables for the current realization j_select #
    ###########################################################
    
    # nrow(data2m) = nb_MOTUs = nb_terms ; ncol(data2m) = nb_samples = nb_doc
    # nrow(documents = gamma) = nb_samples ; ncol(documents) = nb_topics
    # nrow(beta) = nb_topics ; ncol(beta) = nb_MOTUs
    
    if (Rtopicmodels_VEM || Rtopicmodels_Gibbs)
    {
      llh = LLH_final_real1[Ordered_realizations$ix[Selected_real[j_select]]] 
      llh_iterations = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@logLiks
      logbeta = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@beta 
      # topic_compo : proportion of each MOTU in each topic (sums to 1 over MOTUs for each topic)
      topic_compo = exp(t(logbeta))
      documents = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@gamma
      # nb_terms = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@wordassignments$ncol
      # nb_doc = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@wordassignments$nrow
      if (Rtopicmodels_VEM)
        alpha_est = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@alpha
    } else if (ValleExtendedLDA_Gibbs)
    {
      llh = LLH_final_real1[Ordered_realizations$ix[Selected_real[j_select]]] 
      llh_iterations = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]$llk
      topic_compo = t(matrix(Result[[Ordered_realizations$ix[Selected_real[j_select]]]]$phi[nb_iter,],nb_topics,nb_terms))
      # topic_compo : proportion of each MOTU in each topic (sums to 1 over MOTUs for each topic)
      # Because phi_sk's obey independent beta priors in Valle's model, they do not sum to 1
      topic_compo = topic_compo/colSums(topic_compo)
      documents = matrix(Result[[Ordered_realizations$ix[Selected_real[j_select]]]]$theta[nb_iter,],nb_doc,nb_topics)
    }
    
    # KL_topic_compo is intended for computing KL divergence between topics based on their MOTU composition
    KL_topic_compo = topic_compo
    # Setting the minimal proportion of an OTU in a sample as 1/total#reads in KL_topic_compo, 
    # where total#reads is the total number of reads over the dataset
    # (in the case of occurrence data, it is 1/total#occurrences, 
    # where total#occurrences is the sum of OTU occurrences over the dataset): 
    KL_topic_compo[which(topic_compo < 1/sum_data2m)] = 1/sum_data2m
    # Renormalizing KL_topic_compo by topic:
    KL_topic_compo = sweep(KL_topic_compo,MARGIN=2,colSums(KL_topic_compo),`/`)   
    
    # documents : proportion of each topic in a site/document (sums to 1 over topics in each site)
    # KL_documents : proportion of each site/document in a topic (sums to 1 over sites/documents for each topic)
    # KL_documents is intended for computing KL divergence between topics based on their distribution over sites
    KL_documents = documents
    # Setting the minimal proportion of a topic in a sample as 1/total#reads (respectively 1/total#occurrences) in KL_documents:
    KL_documents[which(documents < 1/sum_data2m)] = 1/sum_data2m
    KL_documents = sweep(KL_documents,MARGIN=2,colSums(KL_documents),`/`)
    
    # In the best realization, topics are ordered by decreasing occupancy, while in subsequent realizations, 
    # they are ordered based on their match to topics in the best realization
    # This is the order of the topics in all outputs
    if (j_select == 1)
    {
      # Sorting topics based on their mean proportion over sites
      sort_normal_topic = sort.int(apply(documents,2,mean), index.return=T, decreasing = T)
      documents = documents[,sort_normal_topic$ix]
      KL_documents = KL_documents[,sort_normal_topic$ix]
      topic_compo = topic_compo[,sort_normal_topic$ix]
      KL_topic_compo = KL_topic_compo[,sort_normal_topic$ix]
      
      # if (realization_comparison)
      #   sort_normal_topic_allreal[[j_select]] = sort_normal_topic
      
      KL_documents_bestReal = KL_documents
      
      source(paste0(code_insert,"/LDA_topic_correspondence_fun.R"))
    } else 
    {
      KL_topic_comparison = matrix(nrow=nb_topics,ncol=nb_topics,data=0)
      for (k in 1:nb_topics)
      {
        for (k1 in 1:nb_topics)
        {
          KL_topic_comparison[k,k1] = 1/2*(KL.plugin(KL_documents_bestReal[,k],KL_documents[,k1]) +
                                             KL.plugin(KL_documents[,k1],KL_documents_bestReal[,k]))
          # else if (MOTUwise)
          #   KL_topic_comparison[k,k1] = 1/2*(KL.plugin(KL_topic_compo_allreal[[j_real]][,k],KL_topic_compo_allreal[[j_select]][,k1]) +
          #                                      KL.plugin(KL_topic_compo_allreal[[j_select]][,k1],KL_topic_compo_allreal[[j_real]][,k]))         
        }
      }
      
      Topic_correspondence = LDA_topic_correspondence_fun(KL_topic_comparison,maxmatching,greedymatching)
      
      documents = documents[,Topic_correspondence]
      KL_documents = KL_documents[,Topic_correspondence]
      topic_compo = topic_compo[,Topic_correspondence]
      KL_topic_compo = KL_topic_compo[,Topic_correspondence]
    }
    
    if (realization_comparison)
    {
      documents_allreal[[j_select]] = documents
      KL_documents_allreal[[j_select]] = KL_documents
      
      topic_compo_allreal[[j_select]] = topic_compo
      KL_topic_compo_allreal[[j_select]] = KL_topic_compo
    }
    
    #######################################
    # Comparing realizations: computation #
    #######################################
    
    # Computing correlation between the topics of the highest-llh realization and the others
    if (j_select != 1 && realization_comparison)
    {
      cat("Computing similarities between realizations 1 to",j_select-1,"and realization",j_select,"...\n")  
      
      # if (testdata || !grid && !geographic)
      #   Correlation_samplewise[[j_select-1]] = cor(documents_allreal[[1]],documents)
      
      if (!low_memory)
      {
        if (samplewise)
        {
          KL_documents_jselect_randomized = vector(length=nb_rndzations,mode="list")
          KL_topic_compo_jselect_randomized = NULL
        } else if (MOTUwise)
        {
          KL_topic_compo_jselect_randomized = vector(length=nb_rndzations,mode="list")
          KL_documents_jselect_randomized = NULL
        }
      }
      
      # Computing the mean KL distance between topics for all pairs of realizations.
      # Each realization j_select is compared to the j_select-1 previous realizations, 
      # so as to compare every pair of realizations once.
      for (j_real in 1:(j_select-1))
      {
        # If j_real = 1 and samplewise = 1, Topic_correspondence and KL_topic_comparison have already been computed between the best realization and the realization j_select
        if (j_real !=1 || MOTUwise)
        {
          KL_topic_comparison = matrix(nrow=nb_topics,ncol=nb_topics,data=0)
          for (k in 1:nb_topics)
          {
            for (k1 in 1:nb_topics)
            {
              if (samplewise)
                KL_topic_comparison[k,k1] = 1/2*(KL.plugin(KL_documents_allreal[[j_real]][,k],KL_documents_allreal[[j_select]][,k1]) +
                                                   KL.plugin(KL_documents_allreal[[j_select]][,k1],KL_documents_allreal[[j_real]][,k]))
              else if (MOTUwise)
                KL_topic_comparison[k,k1] = 1/2*(KL.plugin(KL_topic_compo_allreal[[j_real]][,k],KL_topic_compo_allreal[[j_select]][,k1]) +
                                                   KL.plugin(KL_topic_compo_allreal[[j_select]][,k1],KL_topic_compo_allreal[[j_real]][,k]))         
            }
          }
          
          # Computing a one-to-one topic correspondence between j_real and j_select: 
          # every topic in j_select is associated with a single topic in j_real,
          # the best matching among the topics that are not assigned yet. 
          Topic_correspondence = LDA_topic_correspondence_fun(KL_topic_comparison,maxmatching,greedymatching)
        } else
        {
          KL_topic_comparison = KL_topic_comparison[,Topic_correspondence]
          Topic_correspondence = 1:k
        }
        
        # Function LDA_realization_comparison_fun computes nb_rndzations randomizations of the j_select realization for j_real == 1,
        # stores them in KL_documents_jselect_randomized (if samplewise = 1) or KL_topic_compo_jselect_randomized (if MOTUwise = 1),
        # and then reuses them for j_real>1.
        # Function LDA_realization_comparison_lowMemory_fun erases each randomization immediately after it has been computed,
        # thus saving memory at the expense of computational speed, which may be useful for large datasets. 
        if (low_memory)
        {
          LDA_realization_comparison_result = 
            LDA_realization_comparison_lowMemory_fun(j_select,j_real,
                                                     nb_topics,samplewise,MOTUwise,KL_documents_allreal,KL_topic_compo_allreal,
                                                     byRow,testdata,missing,nb_rndzations,nrow_data2m=nb_terms,ncol_data2m=nb_doc,Missing_positions_indices, 
                                                     KL_topic_comparison,Topic_correspondence)
        } else if (!low_memory)
          LDA_realization_comparison_result = 
            LDA_realization_comparison_fun(j_select,j_real,KL_documents_jselect_randomized,KL_topic_compo_jselect_randomized,
                                           nb_topics,samplewise,MOTUwise,KL_documents_allreal,KL_topic_compo_allreal,nrow_grid,ncol_grid,byRow,
                                           testdata,missing,nb_rndzations,ncol_data2m=nb_doc,nrow_data2m=nb_terms,Missing_positions_indices,
                                           KL_topic_comparison,Topic_correspondence)
        
        Mean_KL_topic_comparison_randomized = LDA_realization_comparison_result$Mean_KL_topic_comparison_randomized
        p_value = LDA_realization_comparison_result$p_value
        if (!low_memory)
        {
          if (samplewise)
            KL_documents_jselect_randomized = LDA_realization_comparison_result$KL_documents_jselect_randomized
          else if (MOTUwise)
            KL_topic_compo_jselect_randomized = LDA_realization_comparison_result$KL_topic_compo_jselect_randomized
        }
        
        if (samplewise)
        {
          cor_allReal = cor(documents_allreal[[j_real]],documents_allreal[[j_select]])
        } else if (MOTUwise)
        {
          cor_allReal = cor(topic_compo_allreal[[j_real]],topic_compo_allreal[[j_select]])
        }
        # if (j_real == 1)
        #   Topic_correspondence_to_best_real[[j_select]] = Topic_correspondence
        
        # Saving the KL distance between the current and all previous realizations, averaged over best matching topics
        for (k in 1:nb_topics)
        {
          # The index k cycles through the topics of realization j_real in their original order.
          # The vector Topic_correspondence contains the correspondence between the topics 
          # of realization j_real and the topics of the best realization.
          k0 = which(Topic_correspondence == k)
          p_value_allRealPairs_byTopic[[k]][j_real,j_select] = 1/2*(p_value[k] + p_value[k0])
          
          KL_allRealPairs[j_real,j_select] = KL_topic_comparison[k,Topic_correspondence[k]]/nb_topics + KL_allRealPairs[j_real,j_select]
          KL_allRealPairs_byTopic[[k]][j_real,j_select] = 1/2*(KL_topic_comparison[k,Topic_correspondence[k]] + KL_topic_comparison[k0,k])
          KL_allRealPairs_randomized_byTopic[[k]][j_real,j_select] = 1/2*(Mean_KL_topic_comparison_randomized[k] + Mean_KL_topic_comparison_randomized[k0])
          
          KL_allRealPairs_w_rndzations[j_real,j_select] = KL_topic_comparison[k,Topic_correspondence[k]]/nb_topics + KL_allRealPairs_w_rndzations[j_real,j_select]
          KL_allRealPairs_w_rndzations[j_select,j_real] = Mean_KL_topic_comparison_randomized[k]/nb_topics + KL_allRealPairs_w_rndzations[j_select,j_real]
          
          DKL100_allRealPairs[j_real,j_select] = (Mean_KL_topic_comparison_randomized[k] - KL_topic_comparison[k,Topic_correspondence[k]])/nb_topics + DKL100_allRealPairs[j_real,j_select]
          DKL100_allRealPairs_byTopic[[k]][j_real,j_select] = 1/2*(Mean_KL_topic_comparison_randomized[k] - KL_topic_comparison[k,Topic_correspondence[k]] + Mean_KL_topic_comparison_randomized[k0] - KL_topic_comparison[k0,k])
          
          if (samplewise)
          {
            Correlation_allRealPairs[j_real,j_select] = cor_allReal[k,Topic_correspondence[k]]/nb_topics + Correlation_allRealPairs[j_real,j_select]
            Correlation_allRealPairs_byTopic[[k]][j_real,j_select] = 1/2*(cor_allReal[k,Topic_correspondence[k]] + cor_allReal[k0,k])
          } else if (MOTUwise)
          {
            Correlation_allRealPairs[j_real,j_select] = cor_allReal[k,Topic_correspondence[k]]/nb_topics + Correlation_allRealPairs[j_real,j_select]
            Correlation_allRealPairs_byTopic[[k]][j_real,j_select] = 1/2*(cor_allReal[k,Topic_correspondence[k]] + cor_allReal[k0,k])
          }
        }
        
        llh_differences_allRealPairs[j_real,j_select] = abs(Ordered_realizations$x[j_real]-Ordered_realizations$x[j_select])
        # perplexity_differences_allRealPairs[j_real,j_select] = perplexity1[Ordered_realizations$ix[j_real]]-perplexity1[Ordered_realizations$ix[j_select]]
        p_value_allRealPairs[j_real,j_select] = mean(p_value)
      }
    }
    # End of the loop over j_real
    
    ####################################
    # Assemblage diversity computation #
    ####################################
    
    if (j_select == 1)
    {
      diversity = vector(length = nb_topics+1, mode="numeric")
      for (k in 1:nb_topics)
      {
        #k0 = sort_normal_topic$ix[k]
        k0 = k
        diversity[k] = length(which(topic_compo[,k0] > 1/sum_data2m))
      }
      diversity[nb_topics+1] = nb_terms
      assemblage_names_vect = paste("Assemblage",1:nb_topics)
      diversity = data.frame(Assemblage = c(assemblage_names_vect,"Total diversity"),Diversity = diversity)
      # diversity is saved in "stability_data.frame.rds" and in "Topic_comparison_inBestReal.Rdata"
    }
    
    #############################################################
    # Comparison of the best realization with abiotic variables #
    #############################################################
    
    if ((j_select == 1) && abiotic_variables && grid && !testdata)
    {
      ###############
      # Computation #
      ###############
      
      cat("Computing correlations between the best realization and environmental variables ...\n")
      
      source(paste0(code_insert,"/LDA_environmental_variables_fun_v1.R"))
      
      setwd(abiotic_data_path)
      if (strsplit(abiotic_file,split=".",fixed=T)[[1]][length(strsplit(abiotic_file,split=".",fixed=T)[[1]])] == "csv")
        data_abiotic = read.table(abiotic_file, sep=",", colClasses="numeric", header = T)
      else
        data_abiotic = read.table(abiotic_file, sep="\t", colClasses="numeric", header = T, row.names = 1)
      colnames_abiotic = colnames(data_abiotic)
      ncol0 = ncol(data_abiotic)
      
      if (nrow(data_abiotic) != nb_doc)   
        stop("Environmental and biotic data do not match")     
      
      # data_abiotic = read.table(abiotic_file, sep=";", colClasses="character", header = T)
      # data_abiotic = data_abiotic[as.vector(t(matrix(1:(nrow_grid*ncol_grid),ncol=ncol_grid))),-1]
      # write.table(data_abiotic,file="Environmental_variables1.csv",sep=";",quote=F,row.names=F)
      
      # Call to the LDA_abiotic_variables_fun function
      LDA_abiotic_variables_result = LDA_environmental_variables_fun(grid,missing,Missing_positions_indices = Missing_positions_indices0,byRow,
                                                                     ncol_data2m=nb_doc,pca_abiotic,nb_topics,nb_abiotic_rndzations,
                                                                     documents,data_abiotic)
      
      Correlation_abiotic = LDA_abiotic_variables_result$Correlation_abiotic
      ncol_abiotic = LDA_abiotic_variables_result$ncol0
      colnames_abiotic = LDA_abiotic_variables_result$colnames_abiotic
      if (!geographic)
      {
        Mean_cor_abiotic_comparison_randomized = LDA_abiotic_variables_result$Mean_cor_abiotic_comparison_randomized 
        Sd_cor_abiotic_comparison_randomized = LDA_abiotic_variables_result$Sd_cor_abiotic_comparison_randomized
        p_value_abiotic = LDA_abiotic_variables_result$p_value_abiotic 
      }
      
      cat("Computing correlations between the best realization and environmental variables: Done.\n")
      
      ###########
      # Output #
      ##########
      
      cat("Outputting the comparison to environmental variables ...\n")
      
      abiotic_variables_dirname = paste0(local_subdirname,"Environmental_variables")
      if (!file.exists(abiotic_variables_dirname))
        dir.create(abiotic_variables_dirname)
      setwd(abiotic_variables_dirname)
      
      if (!geographic)
      {
        # Saving the abiotic correlation values and the associated p-values
        save(Correlation_abiotic,Mean_cor_abiotic_comparison_randomized,Sd_cor_abiotic_comparison_randomized,p_value_abiotic,file="Correlation_between_environmental_data_and_assemblages_for_best_real.Rdata")
      } else if (geographic)
        # Saving the abiotic correlation values
        save(Correlation_abiotic,file="Correlation_between_environmental_data_and_assemblages_for_best_real.Rdata")
      
      # Output is written to 3 different formats: .txt, .csv and .rds.
      # Output does not report the comparison to PCA axes in the case pca_abiotic = 1
      Abiotic_file = "Environmental_comparison_verbose.txt"
      if (!geographic)
        abiotic_data.frame = as.data.frame(matrix(ncol=ncol0*3,nrow=nb_topics,data=0))
      else
        abiotic_data.frame = as.data.frame(matrix(ncol=ncol0*2,nrow=nb_topics,data=0))
      rownames(abiotic_data.frame) = paste("Assemblage",1:nb_topics)
      
      if (exists("data_name"))
        write(paste(data_name,"-",nb_topics,"assemblages\n%%%%%%%%%%%%%\n"), file=Abiotic_file, append=F)
      else 
        write(paste(nb_topics,"assemblages\n%%%%%%%%%%%%%\n"), file=Abiotic_file, append=F)
      for (k in 1:nb_topics)
      {
        #k0 = sort_normal_topic$ix[k]
        k0 = k
        
        write(paste0("\n%%%%%%%%%%%%%\nFor assemblage ",k,":\n%%%%%%%%%%%%%"),file=Abiotic_file,append=T)
        for (j in 1:length(colnames_abiotic))
        {
          write(paste0("\n",colnames_abiotic[j],"\n%%%%%%%%%%%%%"),file=Abiotic_file,append=T)
          write(paste("Correlation =",Correlation_abiotic[[1]][k0,j]),file=Abiotic_file,append=T)
          # Comparison to randomizations is not implemented for geographic:
          if (!geographic)
          {
            write(paste("Similarity (with spatial randomizations) =",
                        (Correlation_abiotic[[1]][k0,j] - Mean_cor_abiotic_comparison_randomized[k0,j])/(1-Mean_cor_abiotic_comparison_randomized[k0,j])),file=Abiotic_file,append=T)
            write(paste("p-value for spatial randomizations =",p_value_abiotic[k0,j]),file=Abiotic_file,append=T)
            
            abiotic_data.frame[k,(3*(j-1)+1):(3*(j-1)+3)] = c(Correlation_abiotic[[1]][k0,j],(Correlation_abiotic[[1]][k0,j] - Mean_cor_abiotic_comparison_randomized[k0,j])/(1-Mean_cor_abiotic_comparison_randomized[k0,j]),p_value_abiotic[k0,j])
            if (k==1)
              colnames(abiotic_data.frame)[(3*(j-1)+1):(3*(j-1)+3)] = c(paste(colnames_abiotic[j],"corr."),paste(colnames_abiotic[j],"similarity"),paste(colnames_abiotic[j],"p-value"))  
          } else
          {
            abiotic_data.frame[k,(2*(j-1)+1):(2*(j-1)+2)] = c(Correlation_abiotic[[1]][k0,j],p_value_abiotic[k0,j])
            if (k==1)
              colnames(abiotic_data.frame)[(2*(j-1)+1):(2*(j-1)+2)] = c(paste(colnames_abiotic[j],"corr."),paste(colnames_abiotic[j],"p-value"))
          }
        }
      }
      saveRDS(abiotic_data.frame, file = "Environmental_comparison.rds")
      write.csv(abiotic_data.frame, file = "Environmental_comparison.csv")
      
      # Plotting the outpur for both pca_abiotic = 0 and pca_abiotic = 1
      if (pca_abiotic)
        case1_range = c(1,3)
      else 
        case1_range = 1
      
      for (case1 in case1_range)
      {
        if (case1 == 1)
          pdf("Correlation_between_environmental_data_and_assemblages_for_best_real.pdf")
        else if (case1 == 3)
          pdf("Correlation_between_environmental_data_and_assemblages_for_best_real_pca.pdf")
        par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
        #par(mar=c(5.1,4.1,4.1,2.1)
        #bottom left top right
        par(mar=c(15.1,10.1,4.1,4.1))
        for (k in 1:nb_topics)
        {            
          #k0 = sort_normal_topic$ix[k]
          k0 = k
          plot(Correlation_abiotic[[case1]][k0,],type="p",ann=F,yaxt="n",xaxt="n")
          for (j in 1:ncol0)
          {
            # Bonferroni correction is applied
            if (Correlation_abiotic[[case1+1]][k0,j]>0.05/ncol0)
              lines(j,Correlation_abiotic[[case1]][k0,j],col="red",type="p")
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
            text(1:ncol0, par("usr")[3], srt = 45, adj = c(1.05,1.75), labels = labels, xpd = TRUE, cex=1.1)
          } else 
          {
            if (case1 == 1)
            {
              axis(1, at=1:ncol0, labels = F)
              text(1:ncol0, par("usr")[3], srt = 90, pos=1, offset=1.5, labels = colnames_abiotic, xpd = T, cex=1)
            }
            else if (case1 == 3)
            {
              axis(1, at=1:ncol0, labels = T)
              xlab("PCA axes")
            } 
          }
          
          title(ylab="Correlation value \nbetween spatial distributions")
          if (k==1)
            title("Correlation between environmental data and \nthe 1st assemblage")
          else if (k==2)
            title("Correlation between environmental data and \nthe 2nd assemblage")
          else if (k==3)
            title("Correlation between environmental data and \nthe 3rd assemblage")
          else 
            title(paste0("Correlation between environmental data and the ",k,"th assemblage"))
        }
        dev.off()
      }
      setwd(local_subdirname)
    }
    
    ################################################################
    # Computing the goodness of fit to simulated data for testdata # 
    ################################################################
    
    # Computing KL distance between true documents and the LDA topics, as well as between true documents and randomizations of the LDA topics
    if ((j_select == 1) && testdata && (nb_topics == true_nb_topics))
    {
      source(paste0(code_insert,"/LDA_testdata_fun.R"))
      
      # Call to the LDA_testdata_fun function
      LDA_testdata_result = LDA_testdata_fun(nb_rndzations_true_documents,samplewise,MOTUwise,
                                             ncol_data2m=nb_doc,nrow_data2m=nb_terms,sum_data2m,nb_topics,
                                             KL_documents,documents,true_documents,KL_topic_compo)
      
      Correlation_true_documents = LDA_testdata_result$Correlation_true_documents
      p_value_true_documents = LDA_testdata_result$p_value_true_documents
      KL_topic_comparison_true_documents = LDA_testdata_result$KL_topic_comparison_true_documents
      nES_true_documents = LDA_testdata_result$nES_true_documents
      
      save(Correlation_true_documents,KL_topic_comparison_true_documents,
           p_value_true_documents,nES_true_documents,file="Goodness-of-fit_true_documents.Rdata")
      
      # Outputting the KL comparison between the simulated topics and the inferred topics:
      KL_true_topics_file = paste0("Goodness-of-fit_nb_rndzations",nb_rndzations_true_documents,".txt")
      
      write("KL distance between the simulated and LDA-retrieved assemblages:",file=KL_true_topics_file,append=F) 
      for (k in 1:nb_topics)
        write(KL_topic_comparison_true_documents[k],file=KL_true_topics_file,append=T)
      write(paste("Mean =",mean(KL_topic_comparison_true_documents)),file=KL_true_topics_file,append=T)
      
      write(paste("\nNormalized effect size (nES) based on the mean KL distance between the simulated and LDA-retrieved assemblages over",nb_rndzations_true_documents,"randomizations:"),file=KL_true_topics_file,append=T)
      for (k in 1:nb_topics)
        write(nES_true_documents[k],file=KL_true_topics_file,append=T)
      write(paste("Mean =",mean(nES_true_documents)),file=KL_true_topics_file,append=T)
      
      write(paste("\np-value over",nb_rndzations_true_documents,"randomizations of the KL distance between simulated and LDA-retrieved assemblages:"),file=KL_true_topics_file,append=T)
      for (k in 1:nb_topics)
        write(p_value_true_documents[k],file=KL_true_topics_file,append=T)
      write(paste("Mean =",mean(p_value_true_documents)),file=KL_true_topics_file,append=T)
    }
    
    ################################################
    # Assemblage comparison within the best realization #
    ################################################
    
    if ((j_select == 1) && best_real_comparison)
    {
      # sort_normal_topic
      
      if (low_memory)
      {
        source(paste0(code_insert,"/LDA_best_real_lowMemory_fun_v2.R"))
        LDA_best_real_result = LDA_best_real_lowMemory_fun(nb_rndzations_best_real,nrow_data2m=nb_terms,KL_documents,
                                                           KL_topic_compo,assemblage_names_vect,documents,topic_compo,sum_data2m)
      } else
      {
        source(paste0(code_insert,"/LDA_best_real_fun_v2.R"))
        LDA_best_real_result = LDA_best_real_fun(nb_rndzations_best_real,nrow_data2m=nb_terms,KL_documents,
                                                 KL_topic_compo,assemblage_names_vect,documents,topic_compo,sum_data2m)
      }
      
      Hellinger_topic_comparison_MOTUwise = LDA_best_real_result$Hellinger_topic_comparison_MOTUwise
      Jaccard_topic_comparison_MOTUwise = LDA_best_real_result$Jaccard_topic_comparison_MOTUwise
      Beta.sim_topic_comparison_MOTUwise = LDA_best_real_result$Beta.sim_topic_comparison_MOTUwise
      Corr_topic_comparison_samplewise = LDA_best_real_result$Corr_topic_comparison_samplewise
      Corr_topic_comparison_MOTUwise = LDA_best_real_result$Corr_topic_comparison_MOTUwise
      KL_topic_comparison_samplewise = LDA_best_real_result$KL_topic_comparison_samplewise
      KL_topic_comparison_MOTUwise = LDA_best_real_result$KL_topic_comparison_MOTUwise
      nES_topic_comparison_MOTUwise = LDA_best_real_result$nES_topic_comparison_MOTUwise
      
      cat("Assemblage comparison within the best realization: 100 %\n")
      
      upgma_Hellinger = agnes(Hellinger_topic_comparison_MOTUwise, diss =T, method = "average", keep.diss =F, keep.data =F)
      upgma_Jaccard = agnes(Jaccard_topic_comparison_MOTUwise, diss =T, method = "average", keep.diss =F, keep.data =F)
      upgma_Beta.sim = agnes(Beta.sim_topic_comparison_MOTUwise, diss =T, method = "average", keep.diss =F, keep.data =F)
      
      bestreal_dirname = paste0(local_subdirname,"Assemblage_comparison_within_best_realization")
      if (!file.exists(bestreal_dirname))
        dir.create(bestreal_dirname)
      setwd(bestreal_dirname)
      
      pdf("Hierarchical_clustering_of_assemblages_Hellinger.pdf")
      plot(upgma_Hellinger, which.plots=2, ann=F, cex=2, cex.axis = 2, lwd=1.5)
      title("Average Hellinger")
      dev.off()
      
      pdf("Hierarchical_clustering_of_assemblages_Jaccard.pdf")
      plot(upgma_Jaccard, which.plots=2, ann=F, cex=2, cex.axis = 2, lwd=1.5)
      title("Average Jaccard")
      dev.off()
      
      pdf("Hierarchical_clustering_of_assemblages_beta.sim.pdf")
      plot(upgma_Beta.sim, which.plots=2, ann=F, cex=2, cex.axis = 2, lwd=1.5)
      title("Average beta.sim")
      dev.off()
      
      save(KL_topic_comparison_samplewise, KL_topic_comparison_MOTUwise, nES_topic_comparison_MOTUwise, Corr_topic_comparison_samplewise, Corr_topic_comparison_MOTUwise,
           Hellinger_topic_comparison_MOTUwise, Jaccard_topic_comparison_MOTUwise, Beta.sim_topic_comparison_MOTUwise, diversity, file = "Assemblage_comparison_inBestReal.Rdata")
      
      setwd(local_subdirname)
    }
    
    #######################################################
    # Creating a subdirectory specific to the realization #
    #######################################################
    
    if (Rtopicmodels_Gibbs && best == 2)
      subsub_insert = "closestToMean"
    else
      subsub_insert = "best"
    if (Selected_real[j_select] == 1)
      subsubdirname = paste0(local_subdirname,Selected_real[j_select],"st_",subsub_insert,"_realization/")
    else if (Selected_real[j_select] == 2)
      subsubdirname = paste0(local_subdirname,Selected_real[j_select],"nd_",subsub_insert,"_realization_",matching_insert,"/")
    else if (Selected_real[j_select] == 3)
      subsubdirname = paste0(local_subdirname,Selected_real[j_select],"rd_",subsub_insert,"_realization_",matching_insert,"/")
    else
      subsubdirname = paste0(local_subdirname,Selected_real[j_select],"th_",subsub_insert,"_realization_",matching_insert,"/")
    
    if (!file.exists(subsubdirname))
      dir.create(subsubdirname)
    setwd(subsubdirname)
    
    if (Rtopicmodels_VEM)
    {  
      alpha_est_mean = alpha_est/length_selected_real + alpha_est_mean
      
      alpha_file = "estimated_alpha.txt"
      write(alpha_est,alpha_file)
    }
    
    ########################################
    # Outputting loglikelihood convergence #
    ########################################
    
    pdf("Loglikelihood_convergence.pdf")
    #par(mar=c(5.1,4.1,4.1,2.1))
    par(mar=c(5.1,5.1,4.1,2.1))
    plot((1:length(llh_iterations))*ifelse(Rtopicmodels_Gibbs,llh_keep,1),llh_iterations,ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
    if (Rtopicmodels_Gibbs || ValleExtendedLDA_Gibbs)
      legend(x="bottomright",legend=c(filename_insert,paste("alpha =",alpha),paste("delta =",delta),paste("nb topics =",nb_topics),paste("nb realizations =",nb_real)),col="black",inset=0)
    else if (Rtopicmodels_VEM)
      legend(x="bottomright",legend=c(filename_insert,paste("alpha =",alpha),paste("nb topics =",nb_topics),paste("nb realizations =",nb_real)),col="black",inset=0)
    title("Loglikelihood convergence",cex.main=1.7)
    title(xlab="Iterations",ylab="Loglikelihood value",cex.lab=1.5)
    dev.off()
    
    # Zooming on the end of the curve by removing the first 10 iterations
    if (length(llh_iterations) > 10)
    {
      pdf("Loglikelihood_convergence_zoom.pdf")
      #par(mar=c(5.1,4.1,4.1,2.1))
      par(mar=c(5.1,5.1,4.1,2.1))
      plot((11:length(llh_iterations))*ifelse(Rtopicmodels_Gibbs,llh_keep,1),llh_iterations[-(1:10)],ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
      if (Rtopicmodels_Gibbs || ValleExtendedLDA_Gibbs)
        legend(x="bottomright",legend=c(filename_insert,paste("alpha =",alpha),paste("delta =",delta),paste("nb topics =",nb_topics),paste("nb realizations =",nb_real)),col="black",inset=0)
      else if (Rtopicmodels_VEM)
        legend(x="bottomright",legend=c(filename_insert,paste("alpha =",alpha),paste("nb topics =",nb_topics),paste("nb realizations =",nb_real)),col="black",inset=0)
      title("Loglikelihood convergence (zoom)",cex.main=1.7)
      title(xlab="Iterations",ylab="Loglikelihood value",cex.lab=1.5)
      dev.off()
    }
    
    pdf("Assemblage_prevalence.pdf")
    plot(apply(documents,2,mean),ann=F,ylim = c(0,max(apply(documents,2,mean))),cex.axis=1.5,lwd=2,type="p",col="black")
    abline(h=1/nb_topics,lty=2)
    title(xlab="Community type",ylab="Mean relative proportion across samples",cex.lab=1.5)
    dev.off()
    
    #####################################
    # Outputting posterior distribution #
    #####################################
    
    if (Rtopicmodels_Gibbs && best == 2)
    {
      # if (local_computation)
      #   saveRDS(sampled_logLiks,file="sampled_logLiks.rds")
      # else if (local_existingresult)
      #   sampled_logLiks = readRDS("sampled_logLiks.rds")
      # else if (cluster_existingresult)
      #   sampled_logLiks = readRDS(paste0(cluster_subdirname,"sampled_logLiks.rds"))
 
      pdf("Sampled_posterior_hist.pdf")
      par(mar=c(5.1,5.1,4.1,2.1))
      hist(sampled_logLiks[Selected_real[j_select],],ann=F)
      abline(v=mean(sampled_logLiks[Selected_real[j_select],]),lty=2)
      abline(v=median(sampled_logLiks[Selected_real[j_select],]),lty=2,col="blue")
      title(ylab="Number of samples", xlab="Loglikelihood value", cex.lab=1.5)
      dev.off()
    }
    
    ##########################################
    # Outputting topic taxonomic composition #
    ##########################################
    
    if (exists("taxo_ref") && !testdata)
    {
      # Creating a subdirectory for taxonomic composition outputs 
      # subsubsubdirname = paste(subsubdirname,"assemblage_taxonomic_composition/",sep="")
      # if (!file.exists(subsubsubdirname))
      #   dir.create(subsubsubdirname)
      # setwd(subsubsubdirname)
      
      cat("Outputting assemblage taxonomic composition ...\n")
      
      # More convenient output than topic_compo_taxo (in a single data.frame, and OTUs are not reordered):
      # topic_compo is equal to exp(t(logbeta)), i.e. the proportion of each OTU in each assemblage, (1 assemblage = 1 column), 
      # where assemblages have been reordered according to their match to the assemblages in the best realization
      taxonomic_compo = data.frame(topic_compo,taxo_ref)
      colnames(taxonomic_compo)[1:nb_topics] = paste0("assemblage",1:nb_topics)
      saveRDS(taxonomic_compo,file = "assemblage_composition.rds")
    }
    
    ############################################
    # Outputting topic distribution over sites #  
    ############################################
    
    # subsubsubdirname = paste(subsubdirname,"/topics_site_repartition_info/",sep="")
    # if (!(file.exists(subsubsubdirname)))
    #   dir.create(subsubsubdirname)
    # setwd(subsubsubdirname) 
    
    if (!testdata)
    {
      if ((grid || geographic) && spatial_plot && any(plotted_real == j_select))
      {
        cat("Plotting spatial composition maps ...\n") 
        
        source(paste0(code_insert,"/LDA_spatial_maps_fun_v4.1.R"))
        
        if (geographic && j_select == plotted_real[1])
        {
          if (grid || UnlimitedPointsPerStation)
          {
            # Loading background map: 
            background = readOGR(dsn = paste0(background_map_path,"/",background_map_file), layer = background_map_file)
            background@data$id = rownames(background@data)
            backgroundPoints = fortify(background, region="id")
            backgroundGgplot = merge(backgroundPoints, background@data, by="id") 
          } else
            backgroundGgplot = NA
          
          if (grid)
          {
            if (exists("river_file") && file.exists(paste0(background_map_path,"/",river_file,"/",river_file,".shp")))
            {
              # Loading rivers: 
              rivers = readOGR(dsn = paste0(background_map_path,"/",river_file), layer = river_file)
              rivers@data$id = rownames(rivers@data)
              riversPoints = fortify(rivers, region="id")
              riversGgplot = merge(riversPoints, rivers@data, by="id")
            } else
              riversGgplot = NA
            
            if (exists("land_file") && file.exists(paste0(background_map_path,"/",land_file,"/",land_file,".shp")))
              # Loading land shape for extruding the kriged raster layer
              land = readOGR(dsn = paste0(background_map_path,"/",land_file), layer = land_file)
            else
              land = NA
          } else
          {
            riversGgplot = NA
            land = NA
          }
          
          if (!grid && SurDCMperStation)
            bat = getNOAA.bathy(-180, 180, -90, 90, res = 20, keep=F)
          else 
            bat = NA
        }
        
        LDA_spatial_maps_result = LDA_spatial_maps_fun(grid,geographic,oneplot,
                                                       documents,coord,coord0,
                                                       backgroundGgplot,riversGgplot,land,bat,
                                                       outputdirname = subsubdirname,
                                                       UnlimitedPointsPerStation,
                                                       SurDCMperStation,
                                                       legend_labels = paste("Assemblage",1:nb_topics),
                                                       col_range = "standard")
        
        spatial_topicmix_kriged = LDA_spatial_maps_result$spatial_topicmix_kriged
        
        # if (!(file.exists("Spatial_topicmix_kriged.rds")))
        saveRDS(spatial_topicmix_kriged,file="Spatial_topicmix_kriged.rds")
      
        # assemblage_proportions = do.call("rbind", lapply(spatial_topicmix_kriged, function(g) g$z.pred))
        # dimnames(assemblage_proportions) = list(paste("Assemblage",1:nb_topics),colnames_data2m)
        # saveRDS(assemblage_proportions,file = "assemblage_proportions.rds")
        
        if (grid || geographic)
        {
          tmp.plot = LDA_spatial_maps_result$tmp.plot
          
          ggsave(filename = "Spatial_distribution_maps_kriged.pdf", do.call("arrangeGrob", c(tmp.plot, ncol=min(nb_topics,4))), height = 10/min(nb_topics,4)*((nb_topics-1)%/%4 + 1)*4/3, width = 10)
          if (oneplot)
          {
            tmp.one.plot = LDA_spatial_maps_result$tmp.one.plot
            tmp.one.plot.discrete = LDA_spatial_maps_result$tmp.one.plot.discrete
            ggsave(filename = "Spatial_distribution_maps_kriged_oneplot.pdf",  do.call("arrangeGrob", list(tmp.one.plot, tmp.one.plot.discrete, ncol=2)), height = 10/2*4/3, width = 10)
          }
            
          # pdf("Spatial_distribution_maps_kriged_onebyone.pdf")
          # for (k in 1:nb_topics)
          #   print(tmp.plot[[k]])
          # dev.off()
        } 
      } else if (grid || geographic)
      {
        spatial_topicmix_kriged = vector(length=nb_topics,mode="list")
        for (k in 1:nb_topics)
        {
          spatial_topicmix_kriged[[k]] = data.frame(x=coord$x,y=coord$y,z.pred=documents[,k])
          rownames(spatial_topicmix_kriged[[k]]) = rownames(coord)
        }
        # if (!(file.exists("Spatial_topicmix_kriged.rds")))
        saveRDS(spatial_topicmix_kriged,file="Spatial_topicmix_kriged.rds")
      }
      
        # pdf("Samples_composition.pdf")
        # par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
        # #bottom left top right
        # #par(mar=c(5.1,4.1,4.1,2.1)
        # par(mar=c(15.1,10.1,4.1,4.1))
        # plot(documents[,1],col=terrain.colors(nb_topics)[1],type="p",ann=F,yaxt="n",xaxt="n")
        # axis(2, ylim=c(0,1), col='black')
        # axis(1, at=1:length(colnames_data2m), labels = F)
        # text(1:length(colnames_data2m), par("usr")[3], srt = 45, adj = c(1.05,1.75), labels = colnames_data2m, xpd = TRUE, cex=1.1)
        # for (k in 2:nb_topics)
        #   lines(documents[,k],col=terrain.colors(nb_topics)[k],type="p")
        # title(ylab="Samples composition")
        # title(paste("Distribution of the ",nb_topics,"\nassemblages",sep=""),cex.main=1.7)
        # dev.off()
        
        if (nb_topics < 6)
        {
          col = terrain.colors(nb_topics)
        } else
          col = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),space = "Lab")(nb_topics)
        
        pdf("Samples_composition_barplot.pdf")
        par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=1)
        #bottom left top right
        #par(mar=c(5.1,4.1,4.1,2.1)
        par(mar=c(15.1,10.1,4.1,4.1))
        barplot(t(documents),col=col[sample(1:nb_topics,nb_topics)],names.arg=colnames_data2m,ylim=c(0,1),ann=F,las=3,cex.names=0.5)
        title(ylab="Samples composition",main=paste("Distribution of the ",nb_topics,"\nassemblages",sep=""))
        dev.off()
  
    } else if (testdata)
    {
      pdf("Assemblage_distribution_in_samples.pdf")
      #par(mfrow=c(2,2))
      for (k in 1:nb_topics)
      {
        par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
        plot(documents[,k],type="l",
             main=paste("Assemblage #",k,sep=""),
             xlab="Samples",
             ylab="Assemblage proportion in sample")
      }
      dev.off()
      
      pdf("Assemblage_distribution_in_samples_oneplot.pdf")
      for (k in 1:nb_topics)
      {
        par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
        if (k==1)
        {
          plot(documents[,k],type="l",col=rainbow_hcl(nb_topics)[k],
               main="Sample composition",
               xlab="Samples",
               ylab="Assemblage proportion in sample")
        } else
          lines(documents[,k],type="l",col=rainbow_hcl(nb_topics)[k])
      }
      for (k in 1:true_nb_topics)
        lines(true_documents[,k],type="l",col="black",lty=2)
      dev.off()
    }
  }
  # End of the loop over realizations (loop over j_select)
  
  cat("##############\nEnd of the loop over realizations.\n")
  
  #####################################################
  # Saving the mean estimated alpha over realizations #
  #####################################################
  
  setwd(local_subdirname)
  
  if (Rtopicmodels_VEM)
  {  
    mean_alpha_file = paste0("mean_estimated_alpha_",length_selected_real,"real.txt")
    write(alpha_est_mean,mean_alpha_file)
  }
  
  ##################################
  # Comparing realizations: output #
  ##################################
  
  if (realization_comparison)
  {
    cat("Outputting realization comparison ...\n")
    
    realization_comparison_dirname = paste0(local_subdirname,"Stability_assessment_",MOTU_sample_insert,"_",matching_insert)
    if (!file.exists(realization_comparison_dirname))
      dir.create(realization_comparison_dirname)
    setwd(realization_comparison_dirname)
    
    nES = t(DKL100_allRealPairs)[lower.tri(DKL100_allRealPairs,diag=F)]/KL_allRealPairs_w_rndzations[lower.tri(KL_allRealPairs_w_rndzations,diag=F)]
    nES_byTopic = list()
    for (k in 1:nb_topics)
      nES_byTopic[[k]] = t(DKL100_allRealPairs_byTopic[[k]])[lower.tri(DKL100_allRealPairs_byTopic[[k]],diag=F)]/t(KL_allRealPairs_randomized_byTopic[[k]])[lower.tri(KL_allRealPairs_randomized_byTopic[[k]],diag=F)]
    
    # Comparing realizations averaged over topics:
    save(nES,
         Correlation_allRealPairs,
         KL_allRealPairs,
         llh_differences_allRealPairs,
         p_value_allRealPairs,
         file = paste0("Similarity_allRealPairs_",MOTU_sample_insert,".Rdata"))
    
    save(nES_byTopic,
         Correlation_allRealPairs_byTopic,
         KL_allRealPairs_byTopic,
         llh_differences_allRealPairs,
         p_value_allRealPairs_byTopic,
         #Topic_correspondence_to_best_real,
         file = paste0("Similarity_allRealPairs_byTopic_",MOTU_sample_insert,".Rdata"))
    
    Stability_file = paste0("Stability_",MOTU_sample_insert,"_verbose.txt")
    
    if (exists("data_name"))
      write(paste(data_name,"-",nb_topics,"assemblages -",length_selected_real,"realizations -",MOTU_sample_insert,"\n%%%%%%%%%%%%%\n"),file=Stability_file,append=F)
    else 
      write(paste(nb_topics,"assemblages -",length_selected_real,"realizations -",MOTU_sample_insert,"\n%%%%%%%%%%%%%\n"),file=Stability_file,append=F)
    write(paste("For all assemblages :\n%%%%%%%%%%%%%\n"),file = Stability_file, append=T)
    
    write(paste("Mean correlation =",mean(Correlation_allRealPairs[which(Correlation_allRealPairs!=0)]),"\n"),file=Stability_file,append=T)
    
    write(paste("Mean sKL distance =",mean(KL_allRealPairs[which(!is.na(KL_allRealPairs))]),"\n"),file=Stability_file,append=T)
    
    write(paste("Mean normalized sKL distance difference (Similarity) for",nb_rndzations,"randomizations per assemblage =",mean(nES),"\n"),file=Stability_file,append=T)
    
    # Intercept of nES = f(llh) : 
    fitted_llh = llh_differences_allRealPairs[1,2:length_selected_real]
    fit_llh = lm(nES[1:length(fitted_llh)] ~ fitted_llh)
    write(paste("Similarity = f(llh difference) intercept =",fit_llh$coefficient[1],"\n"),file=Stability_file,append=T)  
    
    # fitted_perplexity = perplexity1[Ordered_realizations$ix[1]]-perplexity1[Ordered_realizations$ix[2:length_selected_real]]
    # fit_perplexity = lm(nES[1:length(fitted_perplexity)] ~ fitted_perplexity)
    # write(paste("Similarity = f(perplexity difference) intercept =",fit_perplexity$coefficient[1],"\n"),file=Stability_file,append=T)  
    
    write(paste("Mean p-value for",nb_rndzations,"randomizations per assemblage =",mean(p_value_allRealPairs[which(!is.na(p_value_allRealPairs))]),"\n"),file=Stability_file,append=T)
    
    stability_data.frame = as.data.frame(matrix(nrow=nb_topics+1,ncol=7,data=0))
    if (exists("data_name"))
      rownames(stability_data.frame)[1] = data_name
    else 
      rownames(stability_data.frame)[1] = ""
    stability_data.frame[1,] = c(diversity$Diversity[nb_topics+1],
                                 mean(Correlation_allRealPairs[which(Correlation_allRealPairs!=0)]),
                                 mean(KL_allRealPairs[which(!is.na(KL_allRealPairs))]),
                                 mean(nES),
                                 fit_llh$coefficient[1],
                                 # fit_perplexity$coefficient[1],
                                 mean(p_value_allRealPairs[which(!is.na(p_value_allRealPairs))]),
                                 ifelse(Rtopicmodels_VEM,alpha_est_mean,alpha))
    
    # Comparing realizations topic by topic:
    for (k in 1:nb_topics)
    {
      #k0 = sort_normal_topic_allreal[[1]]$ix[k]
      k0 = k
      
      p_value = 0
      mean_KL = 0
      mean_nES = 0
      mean_correlation = 0
      for (j_real in 1:(length_selected_real-1))
      {
        # p_value = sum(p_value_allRealPairs_byTopic[[Topic_correspondence_to_best_real[[j_real]][k0]]][j_real,-(1:j_real)])/((length_selected_real-1)*length_selected_real/2) + p_value
        # mean_KL = sum(KL_allRealPairs_byTopic[[Topic_correspondence_to_best_real[[j_real]][k0]]][j_real,-(1:j_real)])/((length_selected_real-1)*length_selected_real/2) + mean_KL
        # mean_nES = sum(DKL100_allRealPairs_byTopic[[Topic_correspondence_to_best_real[[j_real]][k0]]][j_real,-(1:j_real)]/KL_allRealPairs_randomized_byTopic[[Topic_correspondence_to_best_real[[j_real]][k0]]][j_real,-(1:j_real)])/((length_selected_real-1)*length_selected_real/2) + mean_nES
        # mean_correlation = sum(Correlation_allRealPairs_byTopic[[Topic_correspondence_to_best_real[[j_real]][k0]]][j_real,-(1:j_real)])/((length_selected_real-1)*length_selected_real/2) + mean_correlation
        
        p_value = sum(p_value_allRealPairs_byTopic[[k]][j_real,-(1:j_real)])/((length_selected_real-1)*length_selected_real/2) + p_value
        mean_KL = sum(KL_allRealPairs_byTopic[[k]][j_real,-(1:j_real)])/((length_selected_real-1)*length_selected_real/2) + mean_KL
        mean_nES = sum(DKL100_allRealPairs_byTopic[[k]][j_real,-(1:j_real)]/KL_allRealPairs_randomized_byTopic[[k]][j_real,-(1:j_real)])/((length_selected_real-1)*length_selected_real/2) + mean_nES
        mean_correlation = sum(Correlation_allRealPairs_byTopic[[k]][j_real,-(1:j_real)])/((length_selected_real-1)*length_selected_real/2) + mean_correlation
      }
      
      # fit of nES = f(llh) :
      fitted_llh = llh_differences_allRealPairs[1,2:length_selected_real]
      fit_llh = lm(DKL100_allRealPairs_byTopic[[k0]][1,-1]/KL_allRealPairs_randomized_byTopic[[k0]][1,-1] ~ fitted_llh)
      
      # fitted_perplexity = perplexity1[Ordered_realizations$ix[1]]-perplexity1[Ordered_realizations$ix[2:length_selected_real]]
      # fit_perplexity = lm(DKL100_allRealPairs_byTopic[[k0]][1,-1]/KL_allRealPairs_randomized_byTopic[[k0]][1,-1] ~ fitted_perplexity)
      
      write(paste0("For assemblage ",k,":\n%%%%%%%%%%%%%\n"),file=Stability_file,append=T)
      if (exists("data_name"))
        rownames(stability_data.frame)[k+1] = paste("Assemblage",k,data_name)
      else
        rownames(stability_data.frame)[k+1] = paste("Assemblage",k)
      
      write(paste("Mean correlation =",mean_correlation,"\n"),file=Stability_file,append=T)
      write(paste("Mean sKL distance =",mean_KL,"\n"),file=Stability_file,append=T)
      write(paste("Mean normalized sKL difference (Similarity) for",nb_rndzations,"randomizations per assemblage pair =",mean_nES,"\n"),file=Stability_file,append=T)
      write(paste("Similarity = f(llh difference) intercept =",fit_llh$coefficient[1],"\n"),file=Stability_file,append=T)
      # write(paste("Similarity = f(perplexity difference) intercept =",fit_perplexity$coefficient[1],"\n"),file=Stability_file,append=T)
      write(paste("Mean p-value for",nb_rndzations,"randomizations per topic =",p_value,"\n"),file=Stability_file,append=T)
      
      stability_data.frame[k+1,] = c(diversity$Diversity[k],
                                     mean_correlation,
                                     mean_KL,
                                     mean_nES,
                                     fit_llh$coefficient[1],
                                     # fit_perplexity$coefficient[1],
                                     p_value,
                                     NA)
    }
    
    # colnames(stability_data.frame) = c("Diversity","Correlation","SKL","Normalized ES","Normalized ES=f(llh) intercept","Normalized ES=f(perplexity) intercept","p-value","alpha")
    colnames(stability_data.frame) = c("Diversity","Correlation","SKL","Normalized ES","Normalized ES=f(llh) intercept","p-value","alpha")
    saveRDS(stability_data.frame,file=paste0("Stability_",MOTU_sample_insert,".rds"))
    write.csv(stability_data.frame,paste0("Stability_",MOTU_sample_insert,".csv"))
    
    #########
    # Plots #
    #########
    
    pdf(paste0("Similarity_allRealPairs_",MOTU_sample_insert,"_vs_llh-diff.pdf"))  
    par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
    # par(mar = c(bottom, left, top, right))
    par(mar = c(5, 5, 4, 3) + 0.1)
    if (length(which(!is.na(DKL100_allRealPairs))) != length(which(!is.na(llh_differences_allRealPairs))))
      stop("\nError in DKL100_allRealPairs values")
    plot(t(llh_differences_allRealPairs)[lower.tri(llh_differences_allRealPairs,diag=F)],nES,type="p",
         xlab = "Loglikelihood difference", ylab = "Similarity")
    dev.off()
    
    pdf(paste0("Similarity_toBestReal_byTopic_",MOTU_sample_insert,"_vs_llh-diff.pdf"))  
    for (k in 1:nb_topics)
    {
      #k0 = sort_normal_topic_allreal[[1]]$ix[k]
      k0 = k
      par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
      # par(mar = c(bottom, left, top, right))
      par(mar = c(5, 5, 4, 3) + 0.1)
      plot(llh_differences_allRealPairs[1,2:length_selected_real],nES_byTopic[[k0]][1:(length_selected_real-1)],
           type="p", main = paste("Averaged sKL distance\n assemblage", k),
           xlab = "Loglikelihood difference with best realization", ylab = "Similarity to best realization")
    }
    dev.off()
    
    pdf(paste0("Similarity_toBestReal_",MOTU_sample_insert,"_vs_llh-diff.pdf"))  
    par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
    # par(mar = c(bottom, left, top, right))
    par(mar = c(5, 5, 4, 3) + 0.1)
    plot(llh_differences_allRealPairs[1,2:length_selected_real],nES[1:(length_selected_real-1)],type="p",
         xlab = "Loglikelihood difference with best realization", ylab = "Similarity to best realization")
    dev.off()
    
    pdf(paste0("Similarity_toBestReal_",MOTU_sample_insert,"_vs_rank.pdf"))  
    par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
    # par(mar = c(bottom, left, top, right))
    par(mar = c(5, 5, 4, 3) + 0.1)
    plot(seq(2,length_selected_real,1),nES[1:(length_selected_real-1)],type="p",
         xlab = "Realizations sorted by decreasing llh", ylab = "Similarity to best realization")
    dev.off()
    
    pdf(paste0("Correlation_toBestReal_",MOTU_sample_insert,"_vs_llh-diff.pdf"))  
    par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
    plot(llh_differences_allRealPairs[1,2:length_selected_real],Correlation_allRealPairs[1,-1],type="p",
         xlab = "Loglikelihood difference with best realization", ylab = "Correlation coefficient with best realization")
    dev.off()
    
    # Comparing llh between realizations:
    pdf("Realization_llh_comparison.pdf")
    # par(mar = c(bottom, left, top, right))
    par(mar=c(5.1,5.1,4.1,2.1))
    LLH_final_real1_sorted = sort.int(LLH_final_real1,decreasing=T,index.return=T)
    plot(seq(1,nb_real,1),LLH_final_real1_sorted$x,ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
    #title("Loglikelihood for different realizations",cex.main=1.7)
    title(xlab="Ranked realizations",ylab="Loglikelihood value",cex.lab=1.5)
    dev.off()
    
    pdf("Realization_llh_histogram.pdf")
    #par(mar=c(5.1,4.1,4.1,2.1))
    par(mar=c(5.1,5.1,4.1,2.1))
    hist(LLH_final_real1, breaks = 10, freq=T, ann=F, xaxp = c(0, 1, 10))      
    #title("Loglikelihood histogram for different realizations", cex.main=1.7)
    title(ylab="Number of realizations", xlab="Loglikelihood value", cex.lab=1.5)
    dev.off()
    
    if (!grid && !geographic && !testdata)
    {
      ##############################################################
      # Outputting realization comparison for !grid && !geographic #
      ##############################################################
      
      pdf("Samples_composition_allreal_oneplot.pdf")
      par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
      #bottom left top right
      #par(mar=c(5.1,4.1,4.1,2.1)
      par(mar=c(15.1,10.1,4.1,4.1))
      for (k in 1:nb_topics)
      {
        if (k == 1)
        {
          plot(documents_allreal[[1]][,k],type="p",ylim=c(0,1),col=terrain.colors(nb_topics)[k],
               main=paste("Distribution of the ",nb_topics,"\nassemblages",sep=""),yaxt="n",xaxt="n",ann=F)
          axis(2, ylim=c(0,1), col='black')
          axis(1, at=1:19, labels = F)
          text(1:length(colnames_data2m), par("usr")[3], srt = 45, adj = c(1.05,1.75), labels = colnames_data2m, xpd = TRUE, cex=1.1)
        } else
          lines(documents_allreal[[1]][,k],type="p",col=terrain.colors(nb_topics)[k])
        for (j_select in 2:length_selected_real)
        {
          # index0 = which(Correlation_samplewise[[j_select-1]][k,] == max(Correlation_samplewise[[j_select-1]][k,]))
          # lines(documents_allreal[[j_select]][,sort_normal_topic_allreal[[j_select]]$ix[index0]],type="p",col=terrain.colors(nb_topics)[k])
          lines(documents_allreal[[j_select]][,k],type="p",col=terrain.colors(nb_topics)[k])
        }
      }
      title(ylab="Samples composition")
      dev.off()    
    } else if (testdata)
    {
      ##################################################
      # Outputting realization comparison for testdata #
      ##################################################
      
      pdf("Assemblage_distribution_over_samples_allreal.pdf")
      par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
      for (k in 1:nb_topics)
      {
        plot(documents_allreal[[1]][,k],type="l",ylim=c(0,1),
             main=paste("Assemblage #",k,sep=""),
             xlab="Samples",
             ylab="Assemblage proportion")
        for (j_select in 2:length_selected_real)
        {
          # index0 = which(Correlation_samplewise[[j_select-1]][k,] == max(Correlation_samplewise[[j_select-1]][k,]))
          # lines(documents_allreal[[j_select]][,index0],type="l")
          lines(documents_allreal[[j_select]][,k],type="l")
        }
      }
      for (k in 1:true_nb_topics)
        lines(true_documents[,k],type="l",lty=2)
      dev.off()    
      
      pdf("Assemblage_distribution_over_samples_allreal_oneplot.pdf")
      par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
      for (k in 1:nb_topics)
      {
        if (k == 1)
        {
          plot(documents_allreal[[1]][,k],type="l",ylim=c(0,1),col=rainbow_hcl(nb_topics)[k],
               main="Composition of samples",
               xlab="Samples",
               ylab="Assemblage proportions")
        } else
          lines(documents_allreal[[1]][,k],type="l",col=rainbow_hcl(nb_topics)[k])
        for (j_select in 2:length_selected_real)
        {
          # index0 = which(Correlation_samplewise[[j_select-1]][k,] == max(Correlation_samplewise[[j_select-1]][k,]))
          # lines(documents_allreal[[j_select]][,index0],type="l",col=rainbow_hcl(nb_topics)[k])
          lines(documents_allreal[[j_select]][,k],type="l",col=rainbow_hcl(nb_topics)[k])
        }
      }
      for (k in 1:true_nb_topics)
        lines(true_documents[,k],type="l",lty=2)
      dev.off()   
      
      # Computing mean error:
      error_tot = 0
      error_real = vector(length=length_selected_real,mode="numeric")
      documents_mean = matrix(nrow=nb_topics,ncol=nb_doc,data=0)
      for (k in 1:nb_topics)
      {
        index0 = which(Correlation_true_documents[k,] == max(Correlation_true_documents[k,]))
        error_real[1] = mean(abs(documents_allreal[[1]][which(true_documents[,index0]>0),k] 
                                 - true_documents[which(true_documents[,index0]>0),index0]))/nb_topics + error_real[1]
        documents_mean[k,] = documents_allreal[[1]][,k]/length_selected_real 
        for (j_select in 2:length_selected_real)
        {
          # index1 = which(Correlation_samplewise[[j_select-1]][k,] == max(Correlation_samplewise[[j_select-1]][k,]))
          # error_real[j_select] = mean(abs(documents_allreal[[j_select]][which(true_documents[,index0]>0),index1] 
          #                                 - true_documents[which(true_documents[,index0]>0),index0]))/nb_topics + error_real[j_select]
          # documents_mean[k,] = documents_allreal[[j_select]][,index1]/length_selected_real + documents_mean[k,]
          
          error_real[j_select] = mean(abs(documents_allreal[[j_select]][which(true_documents[,index0]>0),k] 
                                          - true_documents[which(true_documents[,index0]>0),index0]))/nb_topics + error_real[j_select]
          documents_mean[k,] = documents_allreal[[j_select]][,k]/length_selected_real + documents_mean[k,]
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
      
      # Plotting the mean error around the mean result over the different runs: 
      pdf("Topic_sample_distribution_error.pdf")
      par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
      for (k in 1:nb_topics)
      {
        if (k == 1)
        {
          plot(documents_mean[k,],type="l",ylim=c(0,1),col=rainbow_hcl(nb_topics)[k],
               main="Composition of samples",
               xlab="Samples",
               ylab="Proportion of the assemblage in the sample")
        } else
          lines(documents_mean[k,],type="l",col=rainbow_hcl(nb_topics)[k])
        lines(documents_mean[k,]+error_tot*documents_mean[k,],type="l",col=rainbow_hcl(nb_topics)[k],lty=2)
        lines(documents_mean[k,]-error_tot*documents_mean[k,],type="l",col=rainbow_hcl(nb_topics)[k],lty=2)
      }
      for (k in 1:true_nb_topics)
        lines(true_documents[,k],type="l",col="black",lty=2)
      dev.off()  
    }
    # Back to the main directory
    setwd(local_subdirname)
  }
  # End of the !mpar condition
} else if (mpar && elbow_aic)
{
  ##########################################
  # Outputting results for mpar            #
  # (comparing several assemblage numbers) #
  ##########################################
  
  # save(nb_topics_range,LLH_final0,LLH_final,perplexity0,perplexity,AIC0,AICc0,AIC,AICc,file="AIC-llh.Rdata")
  
  pdf("Loglikelihood_comparison.pdf")
  #par(mar=c(5.1,4.1,4.1,2.1))
  par(mar=c(5.1,5.1,4.1,2.1))
  errbar(nb_topics_range,LLH_final0[1,],LLH_final0[1,]+LLH_final0[2,],LLH_final0[1,]-LLH_final0[2,],ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
  title(xlab="Number K of assemblages",ylab="Log-likelihood",cex.lab=1.5)
  dev.off()
  
  pdf("Loglikelihood_comparison_allpoints.pdf")
  #par(mar=c(5.1,4.1,4.1,2.1))
  par(mar=c(5.1,5.1,4.1,2.1))
  plot(nb_topics_range,LLH_final[1,],ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
  for (j in 1:nb_real)
    lines(nb_topics_range,LLH_final[j,],lwd=2,type="p",col="black")
  title(xlab="Number K of assemblages",ylab="Log-likelihood",cex.lab=1.5)
  dev.off()
  
  # llh derivative
  ###############
  # derivative = vector(length = length(nb_topics_range)-1, mode="numeric")
  # coor_derivative = vector(length = length(nb_topics_range)-1, mode="numeric")
  # for (i in 2:length(nb_topics_range))
  # {
  #   derivative[i-1] = (LLH_final0[1,i] - LLH_final0[1,i-1])/(nb_topics_range[i] - nb_topics_range[i-1])
  #   coor_derivative[i-1] = (nb_topics_range[i] + nb_topics_range[i-1])/2
  # }
  # pdf("llh_first_derivative.pdf")
  # plot(coor_derivative,derivative,type="o",ylab="First derivative of log-likelihood",xlab="Number K of assemblages",cex.lab=1.3,cex.axis=1.3) 
  # dev.off()
  # 
  # derivative2 = vector(length = length(derivative)-1, mode="numeric")
  # coor_derivative2 = vector(length = length(derivative)-1, mode="numeric")
  # for (i in 2:length(derivative))
  # {
  #   derivative2[i-1] = (derivative[i] - derivative[i-1])/(coor_derivative[i] - coor_derivative[i-1])
  #   coor_derivative2[i-1] = (coor_derivative[i] + coor_derivative[i-1])/2
  # }
  # pdf("llh_second_derivative.pdf")
  # plot(coor_derivative2,derivative2,type="o",ylab="Second derivative of log-likelihood",xlab="Number K of assemblages",cex.lab=1.3,cex.axis=1.3)
  # lines(range(coor_derivative2),c(0,0),lty=2)
  # dev.off()
  # 
  # i = 1
  # while (derivative2[i] < 0 && i < length(derivative2)+1)
  #   i = i+1
  # if (i < length(derivative2)+1)
  # {
  #   optimal_K = coor_derivative2[i]
  # } else
  #   optimal_K = NA
  # write(paste("Llh-based optimal K =",optimal_K),file="optimal_K.txt",append=F)
  # 
  # llh.derivative2.data = data.frame(derivative2 = derivative2, coor_derivative2 = coor_derivative2)
  # saveRDS(llh.derivative2.data, "llh.derivative2.data.rds")
  
  # perplexity derivative
  ###############
  # derivative = vector(length = length(nb_topics_range)-1, mode="numeric")
  # coor_derivative = vector(length = length(nb_topics_range)-1, mode="numeric")
  # for (i in 2:length(nb_topics_range))
  # {
  #   derivative[i-1] = (perplexity0[1,i] - perplexity0[1,i-1])/(nb_topics_range[i] - nb_topics_range[i-1])
  #   coor_derivative[i-1] = (nb_topics_range[i] + nb_topics_range[i-1])/2
  # }
  # pdf("Perplexity_first_derivative.pdf")
  # plot(coor_derivative,derivative,type="o",ylab="First derivative of perplexity",xlab="Number K of assemblages",cex.lab=1.3,cex.axis=1.3) 
  # dev.off()
  
  # derivative2 = vector(length = length(derivative)-1, mode="numeric")
  # coor_derivative2 = vector(length = length(derivative)-1, mode="numeric")
  # for (i in 2:length(derivative))
  # {
  #   derivative2[i-1] = (derivative[i] - derivative[i-1])/(coor_derivative[i] - coor_derivative[i-1])
  #   coor_derivative2[i-1] = (coor_derivative[i] + coor_derivative[i-1])/2
  # }
  # pdf("Perplexity_second_derivative.pdf")
  # plot(coor_derivative2,derivative2,type="o",ylab="Second derivative of perplexity",xlab="Number K of assemblages",cex.lab=1.3,cex.axis=1.3)
  # lines(range(coor_derivative2),c(0,0),lty=2)
  # dev.off()
  
  # i = 1
  # while (derivative2[i] < 0)
  #   i = i+1
  # optimal_K = coor_derivative2[i]
  # write(paste("Perplexity-based optimal K =",optimal_K),file="optimal_K.txt",append=T)
  # 
  # perplexity.derivative2.data = data.frame(derivative2 = derivative2, coor_derivative2 = coor_derivative2)
  # saveRDS(perplexity.derivative2.data, file = "perplexity.derivative2.data.rds")
  
  if (Rtopicmodels_VEM)
  {
    pdf("AIC_comparison.pdf")
    #par(mar=c(5.1,4.1,4.1,2.1))
    par(mar=c(5.1,5.1,4.1,2.1))
    
    errbar(nb_topics_range,AIC0[1,],AIC0[1,]+AIC0[2,],AIC0[1,]-AIC0[2,],ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
    
    title(xlab="Number K of assemblages",ylab="AIC",cex.lab=1.9)
    dev.off()
    
    pdf("AIC_comparison_allpoints.pdf")
    #par(mar=c(5.1,4.1,4.1,2.1))
    par(mar=c(5.1,5.1,4.1,2.1))
    
    plot(nb_topics_range,AIC[1,],ann=F,cex.axis=1.5,lwd=2,type="p",col="black",ylim=range(AIC))
    for (j in 1:nb_real)
      lines(nb_topics_range,AIC[j,],lwd=2,type="p",col="black")
    
    title(xlab="Number K of assemblages",ylab="AIC",cex.lab=1.9)
    dev.off()
    
    pdf("AICc_comparison.pdf")
    #par(mar=c(5.1,4.1,4.1,2.1))
    par(mar=c(5.1,5.1,4.1,2.1))
    
    errbar(nb_topics_range,AICc0[1,],AICc0[1,]+AICc0[2,],AICc0[1,]-AICc0[2,],ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
    
    title(xlab="Number K of assemblages",ylab="AICc",cex.lab=1.9)
    dev.off()
    
    pdf("AICc_comparison_allpoints.pdf")
    #par(mar=c(5.1,4.1,4.1,2.1))
    par(mar=c(5.1,5.1,4.1,2.1))
    
    plot(nb_topics_range,AICc[1,],ann=F,cex.axis=1.5,lwd=2,type="p",col="black",ylim=range(AICc))
    for (j in 1:nb_real)
      lines(nb_topics_range,AICc[j,],lwd=2,type="p",col="black")
    
    title(xlab="Number K of assemblages",ylab="AICc",cex.lab=1.9)
    dev.off()
    
    pdf("Estimated_alpha.pdf")
    #par(mar=c(5.1,4.1,4.1,2.1))
    par(mar=c(5.1,5.1,4.1,2.1))
    
    errbar(nb_topics_range,alpha_est0[1,],alpha_est0[1,]+alpha_est0[2,],alpha_est0[1,]-alpha_est0[2,],ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
    
    title("Estimated alpha value vs number of assemblages",cex.main=1.7)
    title(xlab="Number of assemblages",ylab="Estimated alpha value",cex.lab=1.5)
    dev.off()
  }
} else if (mpar && cross_validation)
{
  if (local_existingresult)
    load("perplexity.Rdata")
  else if (cluster_existingresult)
    load(paste0(cluster_subdirname,"perplexity.Rdata"))
    
  pdf("Perplexity_comparison.pdf")
  #par(mar=c(5.1,4.1,4.1,2.1))
  par(mar=c(5.1,5.1,4.1,2.1))
  errbar(nb_topics_range,colMeans(perplexity_mpar),colMeans(perplexity_mpar)+apply(perplexity_mpar,2,sd),colMeans(perplexity_mpar)-apply(perplexity_mpar,2,sd),ann=F,cex.axis=1.5,lwd=2,type="p",col="black")
  title(xlab="Number K of assemblages",ylab="Perplexity",cex.lab=1.5)
  dev.off()
}
#### end of mpar condition

cat("End of postprocessing.\n")

if (local_computation)
{
  end.time <- Sys.time()
  setwd(local_subdirname)
  time_file = "Computation_time.txt"
  write(end.time-start.time,time_file)
}


