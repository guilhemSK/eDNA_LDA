library(plotrix)
library(topicmodels,lib.loc="/home/gsommeriaklein/packages_R/")
library(Hmisc)
#library(lattice)

# Exécution lda :
#./lda est 0.1 20 settings.txt /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/H20_GH_sequences_lda_norep.txt random /Users/guilhemsommeria-klein/Desktop/These/Données_H20/algo_result1
#Exécution hdp-faster :
#./hdp --train_data /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/H20_GH_sequences_lda_norep.txt --directory /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/hdp-faster_alpha1_gamma1_eta0.01_maxiter100_samplehyper_norep_1/ --max_time -1 --max_iter 100 --save_lag 10 --verbose --sample_hyper yes/no
# Exécution hdp :
#./hdp --algorithm train --data /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/H20_GH_sequences_lda_norep.txt --directory /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/hdp_alpha1_gamma1_eta0.01_maxiter100_inittopic100_norep_1/ --max_iter 100 --save_lag 10 --eta 0.01 --split_merge yes/no --sample_hyper yes/no --init_topics 100

# to know where packages are stored : '.libPaths' in the R console, then access variable .Library

# Connexion ssh :
# ssh gsommeriaklein@10.2.3.11
# sshfs gsommeriaklein@10.2.3.11:/home/gsommeriaklein/ /Users/guilhemsommeria-klein/Desktop/Serveur_EDB/ 

local = 1
existingresult = 0
cluster = 0

if (local)
{
  start.time <- Sys.time()
  
  library(topicmodels)
}

data_h20 = 0
data_pp = 0
data_betadiv = 0
data_betadiv_pooled = 1

barcode_gh = 0
barcode_ghassigned = 0
barcode_itsfungi = 0
barcode_16sbact = 1
barcode_18s = 0
barcode_18sfungi = 0
barcode_18smetazoa = 0
barcode_18splants = 0
barcode_16sarch = 0
barcode_itsplant = 0
barcode_16sins = 0
plantfungi = 0
testdata = 0
testdata_dir = "Discrete-mixed1_samples_nbtopics3_nbmotus100_randomtopics"

blei = 0
Rtopicmodels_Gibbs = 0
Rtopicmodels_VEM = 1

mpar = 0
nb_topics = 5

# nb_topics_range = c(2,3)
# nb_topics_range = c(2,3,4,5,10,15,20,30,40,50,60,70)
#nb_topics_range = c(30,40,50,60,70)
mnb_topics = 0

if (mpar) {
  if (mnb_topics)
    mpar_range = length(nb_topics_range)
  #else if (miter)
  #mpar_range = length(miter_range)
} else mpar_range=1

delta = 0.1
#alpha_insert = "alpha50:nb_topics"
alpha_insert = "alpha0.1"

nb_real = 20
# if (mpar && mnb_topics) or if (!mpar && best_keep), we must have best=0 
# best_keep = 1 has no influence if (mpar && mnb_topics)
best = 0
best_keep = 1
# select_real=1 only possible if best_keep=1
select_real = 1
#Selected_real = c(seq(1,10,1),seq(91,100,1))
Selected_real = seq(1,20,1)
#Selected_real = c(1,2)

# only useful for Gibbs sampling
nb_iter = 2000
llh_keep = 100

# only useful for VEM
#em_tol = 5*(10^-5)
em_tol = 10^-7
var_tol = 10^-8

if (data_h20) {
  data_insert = "Données_H20"
  short_data_insert = "H20"
} else if (data_pp)
{
  data_insert = "Données_PetitPlateau"
  short_data_insert = "PP"
} else if (data_betadiv || data_betadiv_pooled)
{
  data_insert = "Données_betadiv"
  short_data_insert = "betadiv"
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
} else if (barcode_16sbact)
{
  barcode_insert = "Bactéries_16S"
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
  barcode_insert = "Métazoaires_18S"
} else if (barcode_18splants) 
{
  barcode_insert = "Plantes_18S"
} else if (barcode_16sarch)
{
  barcode_insert = "Archées_16S"
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
} else if (testdata)
{
  barcode_insert = paste("Test_data/",testdata_dir,sep="")
}

# if (blei) {
#   filename_insert = "Blei_LDA"
# } else if (Rtopicmodels_Gibbs) {
#   filename_insert = "Rtopicmodels_LDA_Gibbs"
# } else if (Rtopicmodels_VEM)
#   filename_insert = "Rtopicmodels_LDA_VEM"

###############
setwd(paste("/home/gsommeriaklein/",data_insert,"/",barcode_insert,"/",sep=""))

data2 <- read.table(paste(short_data_insert,"_",short_barcode_insert,"_sequences_counts_norep.txt",sep=""), sep=" ")
data2m = as.matrix(data2[-1,])
#data2m_num = data.matrix(data2)
data2m_num = matrix(nrow=nrow(data2m),ncol=ncol(data2m))
# count2 contains the total count fo each sequence
# the MOTUs (rows) for which there is zero reads over all samples are removed from data2m and taxo_ref
count2 = vector(length=nrow(data2m),mode="numeric")
for (i in 1:nrow(data2m))
{
  for (j in 1:ncol(data2m))
  {
    data2m_num[i,j] = as.numeric(data2m[i,j])
    count2[i]=data2m_num[i,j]+count2[i]
  }  
}
data2m_num = data2m_num[-which(count2==0),]
data2m = data2m[-which(count2==0),]
count2 = count2[-which(count2==0)]

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
  alpha = 50/nb_topics
  
  getClass(Class="LDA_VEM",package="topicmodels")
  Result_1real = new("LDA_VEM")
  
  write("settings.txt")
  # ...
  
  SEED = vector(length=nb_real,mode="integer")
  for (j in 1:nb_real)
    SEED[j] = as.integer(Sys.time()) - j*10^7
  control_LDA_VEM = list(estimate.alpha=TRUE, alpha=alpha, estimate.beta=TRUE,
                         verbose = 1, prefix = tempfile(), save = 0, keep = 1,
                         seed = SEED, nstart = nb_real, best = best,
                         var = list(iter.max = 500, tol = var_tol),
                         em = list(iter.max = 1000, tol = em_tol),
                         initialize = "random")
  
  command = paste("./lda est 0.1 20 settings.txt /Users/guilhemsommeria-klein/Desktop/These/",data_insert,"Données_betadiv/Plantes_GH/H20_GH_sequences_lda_norep.txt random /Users/guilhemsommeria-klein/Desktop/These/Données_H20/algo_result1")
    paste("cp ",cluster_dirname,filename," ",local_dirname,filename,sep="")
  system(command, intern=TRUE)
  
  
  Result = topicmodels::LDA(x=t(data2m_num),k=nb_topics,method = "VEM",control=control_LDA_VEM,model=NULL)
  
  llh = LLH_final_real1[Ordered_realizations$ix[Selected_real[j_select]]] 
  llh_iterations = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@logLiks
  logbeta = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@beta 
  documents = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@gamma
  nb_terms = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@wordassignments$ncol
  nb_doc = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@wordassignments$nrow
  alpha_est = Result[[Ordered_realizations$ix[Selected_real[j_select]]]]@alpha
  
  if (par_index == 1)
    Result_mpar = Result
  else if (par_index > 1)
    Result_mpar = c(Result_mpar,Result)
  
  if (!best)
  {
    for (j in 1:nb_real)
    { 
      if (max(Result[[j]]@logLiks)>=0)
        cat("Error : positive loglikelihood value\n")
    }
  }
  
}
 
#########################

if (!mpar)
{
  if (data_h20) {
    dirname = paste("/home/gsommeriaklein/",data_insert,"/plots/",barcode_insert,"/",filename_insert,"/",sep="")
  } else if (data_pp)
    dirname = paste("/home/gsommeriaklein/",data_insert,"/",barcode_insert,"/",filename_insert,"/",sep="")
  
  if (!(file.exists(dirname)))
  {dir.create(dirname)}
  setwd(dirname)  
  
  if (!best && !best_keep)
  {  
    
    if (Rtopicmodels_Gibbs)
    {
      filename = paste(filename_insert,"_alpha",alpha,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,".Rdata",sep="")
      save(Result,file=filename)
    } else if (Rtopicmodels_VEM) {
      filename = paste(filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,".Rdata",sep="")
      save(Result,file=filename)
    }
    cat(filename)

  } else if (best || best_keep)
  {
    
    if (best) 
    {
      if (Rtopicmodels_Gibbs) {
        subdirname = paste(dirname,filename_insert,"_alpha",alpha,"_delta",delta,"_nb_topics",nb_topics,"_nb_iter",nb_iter,"_nb_real",nb_real,"_best/",sep="")
      } else if (Rtopicmodels_VEM) {
        subdirname = paste(dirname,filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best/",sep="")
      }
      if (!(file.exists(subdirname)))
      {dir.create(subdirname)}
      setwd(subdirname) 
      
      if (Rtopicmodels_Gibbs)
      {
        filename = paste(filename_insert,"_alpha",alpha,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,"_best.Rdata",sep="")
        save(Result,file=filename)
      } else if (Rtopicmodels_VEM) {
        filename = paste(filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best.Rdata",sep="")
        save(Result,file=filename)
      }
      cat(filename)
    
      } else if (best_keep)
      {
        if (Rtopicmodels_Gibbs) {
          subdirname = paste(dirname,filename_insert,"_alpha",alpha,"_delta",delta,"_topics",nb_topics,"_nb_iter",nb_iter,"_nb_real",nb_real,"_best_keep/",sep="")
        } else if (Rtopicmodels_VEM) {
          subdirname = paste(dirname,filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best_keep/",sep="")
        }
        if (!(file.exists(subdirname)))
        {dir.create(subdirname)}
        setwd(subdirname) 
        
        if (Rtopicmodels_Gibbs)
        {
          filename = paste(filename_insert,"_alpha",alpha,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,"_best_keep.Rdata",sep="")
          save(Result,file=filename)
        } else if (Rtopicmodels_VEM) {
          filename = paste(filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best_keep.Rdata",sep="")
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
    dirname = paste("/home/gsommeriaklein/",data_insert,"/plots/",barcode_insert,"/",filename_insert,"/",sep="")
  } else if (data_pp)
    dirname = paste("/home/gsommeriaklein/",data_insert,"/",barcode_insert,"/",filename_insert,"/",sep="")
  
  if (!(file.exists(dirname)))
  {dir.create(dirname)}
  setwd(dirname)  
  
  if (Rtopicmodels_Gibbs) {
    subdirname = paste(dirname,filename_insert,"_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_iter",nb_iter,"_nb_real",nb_real,"/",sep="")
  } else if (Rtopicmodels_VEM) {
    subdirname = paste(dirname,filename_insert,"_",alpha_insert,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,"/",sep="")
  }
  if (!(file.exists(subdirname)))
  {dir.create(subdirname)}
  setwd(subdirname)
  
  if (mnb_topics)
  {
    
    if (Rtopicmodels_Gibbs) {
      filename = paste(filename_insert,"_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,"_nb_iter",nb_iter,".Rdata",sep="")
      save(Result_mpar,file=filename)
    } else if (Rtopicmodels_VEM)  {
      filename = paste(filename_insert,"_",alpha_insert,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,".Rdata",sep="")
      save(Result_mpar,file=filename)
    }
    cat(filename)
    
  }
}
    
    




