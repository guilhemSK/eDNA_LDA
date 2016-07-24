library(plotrix)
#library(topicmodels,lib.loc="/home/gsommeriaklein/packages_R/")
library(topicmodels,lib.loc="/home/gsommeria/packages_R/")
library(Hmisc)
#library(lattice)

# Exécution lda :
#./lda est 0.1 20 settings.txt /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/H20_GH_sequences_lda_norep.txt random /Users/guilhemsommeria-klein/Desktop/These/Données_H20/algo_result1
#Exécution hdp-faster (default parameter values provided in main.cpp) :
#./hdp --train_data /Users/guilhemsommeria-klein/Desktop/These/Donnees_PetitPlateau/Bacteries_16S/PP_16bact_sequences_lda_norep.txt --directory /Users/guilhemsommeria-klein/Desktop/These/Donnees_PetitPlateau/Bacteries_16S/HDP-faster/hdp-faster_alpha1_gamma1_eta0.01_maxiter100_samplehyper_norep_1/ --max_time -1 --max_iter 100 --save_lag 10 --verbose --sample_hyper
# Exécution hdp :
#./hdp --algorithm train --data /Users/guilhemsommeria-klein/Desktop/These/Donnees_H20/Plantes_GH/H20_GH_sequences_lda_norep.txt --directory /Users/guilhemsommeria-klein/Desktop/These/Donnees_H20/Plantes_GH/hdp_alpha1_gamma1_eta0.01_maxiter100_inittopic100_norep_1/ --max_iter 100 --save_lag 10 --eta 0.01 --split_merge yes/no --sample_hyper yes/no --init_topics 100
#./hdp --algorithm train --data /Users/guilhemsommeria-klein/Desktop/These/Donnees_PetitPlateau/Bacteries_16S/PP_16bact_sequences_lda_norep.txt --directory /Users/guilhemsommeria-klein/Desktop/These/Donnees_PetitPlateau/Bacteries_16S/HDP/hdp_alpha1_gamma1_eta0.01_maxiter100_inittopic3/ --max_iter 100 --save_lag 10 --eta 0.01 --split_merge yes/no --sample_hyper yes/no --init_topics 3

# to know where packages are stored : '.libPaths' in the R console, then access variable .Library

start.time <- Sys.time()

data_h20 = 0
data_pp = 1
# Not implemented:
# data_gs = 0
filled = 1

# Treating the data as occurrence data instead of abundance data
occurrence = 0

# Number of sites an OTU must ocuppy 
# (if nb_occupied_sites_threshold = 1, all OTUs with non-zero abundance are kept)
nb_occupied_sites_threshold = 1
# Removing OTUs with less reads than the number of sites
no_rare = 0

barcode_gh = 0
barcode_ghassigned = 0
barcode_itsfungi = 0
barcode_itsasco = 0
barcode_itsbasidio = 0
barcode_itsglomero= 0
barcode_16sbact = 0
barcode_18s = 0
barcode_18sfungi = 0
barcode_18smetazoa = 0
barcode_18sannelids = 1
barcode_18sarthropods = 0
barcode_18snematodes = 0
barcode_18sprotists = 0
barcode_18splatyhelminthes = 0
barcode_18splants = 0
barcode_16sarch = 0
barcode_itsplant = 0
plantfungi = 0
# Not implemented:
# frogs = 0
testdata = 0
testdata_dir = "Continuous-mixed_samples_nbtopics5_nbmotus1000_randomtopics_sampledreads"

blei = 0
Rtopicmodels_Gibbs = 0
Rtopicmodels_VEM = 1

mpar = 0
nb_topics = 3

#nb_topics_range = c(2,3,4,5,6,7,8)
nb_topics_range = c(2,3,4,5,7,9,11,13,15,17,19)
# nb_topics_range = c(2,3)
# nb_topics_range = c(2,3,4,5,10,15,20,30,40,50,60,70)
#nb_topics_range = c(30,40,50,60,70)
mnb_topics = 1

if (mpar) {
  if (mnb_topics)
    mpar_range = length(nb_topics_range)
  #else if (miter)
  #mpar_range = length(miter_range)
} else mpar_range=1

delta = 0.1
alpha_insert = "alpha0.1"
#alpha_insert = "alpha50_over_nb_topics"

nb_real = 2
# if (mpar && mnb_topics) or if (!mpar && best_keep), we must have best=0 
# best_keep = 1 has no influence if (mpar && mnb_topics)
best = 0
best_keep = 1

# only useful for Gibbs sampling
nb_iter = 2000
llh_keep = 100

# only useful for VEM
#em_tol = 5*(10^-5)
em_tol = 10^-7
var_tol = 10^-8

if (data_h20) {
  data_insert = "Donnees_H20"
  short_data_insert = "H20"
} else if (data_pp)
{
  data_insert = "Donnees_PetitPlateau"
  short_data_insert = "PP"
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
} else if (plantfungi)
{
  barcode_insert = "Plantes_GH-Champignons_ITS"
} else if (testdata)
{
  barcode_insert = paste("Test_data/",testdata_dir,sep="")
} else if (frogs)
{
  barcode_insert = "Anoures"
  short_barcode_insert = "Frogs"
}

if (blei) {
  filename_insert = "Blei_LDA"
} else if (Rtopicmodels_Gibbs) {
  filename_insert = "Rtopicmodels_LDA_Gibbs"
} else if (Rtopicmodels_VEM)
  filename_insert = "Rtopicmodels_LDA_VEM"

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

###############
setwd(paste("/home/gsommeria/work/",data_insert,"/",barcode_insert,"/",sep=""))

if (plantfungi) {
  load("data2m_plantfungi.Rdata")
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
  #load("true_norm_documents.Rdata")
  true_documents = true_norm_documents
} else 
{
  if (filled)
  {
    if (data_pp || data_h20)
      load("data2m_filled.Rdata")
  } else
  {
    data2 = read.table(paste(short_data_insert,"_",short_barcode_insert,"_sequences_counts_norep.txt",sep=""), colClasses="vector", sep=" ")
    data2m = as.matrix(data2[-1,])
    colnames(data2m) = data2[1,] 
  }
}
rownames(data2m) = seq(1,nrow(data2m),1)
data2m = apply(data2m,2,as.numeric)

OTUs_to_be_removed = vector(length=nrow(data2m),mode="logical")
# Prop_OTU_removed = 0
# Prop_reads_removed = 0
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
#   Prop_OTU_removed = length(which(OTUs_to_be_removed))/nrow(data2m)
#   Prop_reads_removed = sum(data2m[which(OTUs_to_be_removed),])/sum(data2m)
  data2m = data2m[-which(OTUs_to_be_removed),]
  rownames(data2m) = seq(1,nrow(data2m),1)
  count2 = rowSums(data2m)
}

if (occurrence)
  data2m[data2m>0] = 1

# Removing empty sites even after filling (this may arise when data2m results from a split of the original dataset)
if (any(colSums(data2m)==0))  
  data2m = data2m[,-which(colSums(data2m)==0)] 

# Start of the loop over parameters
#########################
for (par_index in 1:mpar_range)
{
  if (mpar)
  {
    cat("Parameter value",par_index,"\n")  
  
    if (mnb_topics)
      nb_topics = nb_topics_range[par_index]  
  }
  if (alpha_insert == "alpha50_over_nb_topics")
    alpha = 50/nb_topics
  else if (alpha_insert == "alpha0.1")
    alpha = 0.1

  # beta = #topics x #terms
  # gamma = # documents x #topics
  if (Rtopicmodels_Gibbs)
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
    
    if (par_index == 1)
      Result_mpar = Result
    else if (par_index > 1)
      Result_mpar = c(Result_mpar,Result)
  } else if (Rtopicmodels_VEM)
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
    
    if (par_index == 1)
      Result_mpar = Result
    else if (par_index > 1)
      Result_mpar = c(Result_mpar,Result)
  }
}
 
#########################

if (!mpar)
{
  if (data_h20)
  {
    #dirname = paste("/home/gsommeriaklein/",data_insert,"/plots/",barcode_insert,"/",filename_insert,"/",sep="")
    dirname = paste("/home/gsommeria/work/",data_insert,"/plots/",barcode_insert,"/",filename_insert,"/",sep="")
  } else if (data_pp)
  {
    if (filled || testdata)
    {
      #dirname = paste("/home/gsommeriaklein/",data_insert,"/",barcode_insert,"/",filename_insert,"/",sep="")
      dirname = paste("/home/gsommeria/work/",data_insert,"/",barcode_insert,"/",filename_insert,"/",sep="")
    } else 
      dirname = paste("/home/gsommeria/work/",data_insert,"/",barcode_insert,"/",filename_insert,"/Unfilled/",sep="")
  }
  
  if (!(file.exists(dirname)))
    dir.create(dirname)
  setwd(dirname)  
  
  if (!best && !best_keep)
  {  
    
    if (Rtopicmodels_Gibbs)
    {
      filename = paste(filename_insert,"_alpha",alpha,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,occurrence_insert,remove_single_sites_insert,no_rare_insert,".Rdata",sep="")
      save(Result,file=filename)
    } else if (Rtopicmodels_VEM) {
      filename = paste(filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,occurrence_insert,remove_single_sites_insert,no_rare_insert,".Rdata",sep="")
      save(Result,file=filename)
    }
    cat(filename)

  } else if (best || best_keep)
  {
    
    if (best) 
    {
      if (Rtopicmodels_Gibbs) {
        subdirname = paste(dirname,filename_insert,"_alpha",alpha,"_delta",delta,"_nb_topics",nb_topics,"_nb_iter",nb_iter,"_nb_real",nb_real,"_best",occurrence_insert,remove_single_sites_insert,no_rare_insert,"/",sep="")
      } else if (Rtopicmodels_VEM) {
        subdirname = paste(dirname,filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best",occurrence_insert,remove_single_sites_insert,no_rare_insert,"/",sep="")
      }
      if (!(file.exists(subdirname)))
      {dir.create(subdirname)}
      setwd(subdirname) 
      
      if (Rtopicmodels_Gibbs)
      {
        filename = paste(filename_insert,"_alpha",alpha,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,"_best",occurrence_insert,remove_single_sites_insert,no_rare_insert,".Rdata",sep="")
        save(Result,file=filename)
      } else if (Rtopicmodels_VEM) {
        filename = paste(filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best",occurrence_insert,remove_single_sites_insert,no_rare_insert,".Rdata",sep="")
        save(Result,file=filename)
      }
      cat(filename)
    
      } else if (best_keep)
      {
        if (Rtopicmodels_Gibbs) {
          subdirname = paste(dirname,filename_insert,"_alpha",alpha,"_delta",delta,"_topics",nb_topics,"_nb_iter",nb_iter,"_nb_real",nb_real,"_best_keep",occurrence_insert,remove_single_sites_insert,no_rare_insert,"/",sep="")
        } else if (Rtopicmodels_VEM) {
          subdirname = paste(dirname,filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best_keep",occurrence_insert,remove_single_sites_insert,no_rare_insert,"/",sep="")
        }
        if (!(file.exists(subdirname)))
          dir.create(subdirname)
        setwd(subdirname) 
        
        if (Rtopicmodels_Gibbs)
        {
          filename = paste(filename_insert,"_alpha",alpha,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,"_best_keep",occurrence_insert,remove_single_sites_insert,no_rare_insert,".Rdata",sep="")
          save(Result,file=filename)
        } else if (Rtopicmodels_VEM)
        {
          filename = paste(filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best_keep",occurrence_insert,remove_single_sites_insert,no_rare_insert,".Rdata",sep="")
          save(Result,file=filename)
        }
        cat(filename)
      }
  }
###### end of the !mpar condition
} else if (mpar)
{
  #####
  
  # subdirname = "several_topic_nb/"
  # if (!(file.exists(dirname)))
  #   {dir.create(dirname)}
  # setwd(subdirname)
  
  if (data_h20) {
    #dirname = paste("/home/gsommeriaklein/",data_insert,"/plots/",barcode_insert,"/",filename_insert,"/",sep="")
    dirname = paste("/home/gsommeria/work/",data_insert,"/plots/",barcode_insert,"/",filename_insert,"/",sep="")
  } else if (data_pp)
  {
    if (filled || testdata)
    {
      #dirname = paste("/home/gsommeriaklein/",data_insert,"/",barcode_insert,"/",filename_insert,"/",sep="")
      dirname = paste("/home/gsommeria/work/",data_insert,"/",barcode_insert,"/",filename_insert,"/",sep="")
    } else
      dirname = paste("/home/gsommeria/work/",data_insert,"/",barcode_insert,"/",filename_insert,"/Unfilled/",sep="")
  }
  
  if (!(file.exists(dirname)))
    dir.create(dirname)
  setwd(dirname)  
  
  if (Rtopicmodels_Gibbs)
  {
    subdirname = paste(dirname,filename_insert,"_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_iter",nb_iter,"_nb_real",nb_real,occurrence_insert,remove_single_sites_insert,no_rare_insert,"/",sep="")
  } else if (Rtopicmodels_VEM)
  {
    subdirname = paste(dirname,filename_insert,"_",alpha_insert,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,occurrence_insert,remove_single_sites_insert,no_rare_insert,"/",sep="")
  }
  if (!(file.exists(subdirname)))
    dir.create(subdirname)
  setwd(subdirname)
  
  if (mnb_topics)
  {
    
    if (Rtopicmodels_Gibbs)
    {
      filename = paste(filename_insert,"_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,"_nb_iter",nb_iter,occurrence_insert,remove_single_sites_insert,no_rare_insert,".Rdata",sep="")
      save(Result_mpar,file=filename)
    } else if (Rtopicmodels_VEM)
    {
      filename = paste(filename_insert,"_",alpha_insert,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,occurrence_insert,remove_single_sites_insert,no_rare_insert,".Rdata",sep="")
      save(Result_mpar,file=filename)
    }
    cat(filename)
    
  }
}
    

end.time <- Sys.time()
setwd(subdirname)
time_file = "Computation_time.txt"
write(end.time-start.time,time_file)
    




