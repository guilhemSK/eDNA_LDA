library(plotrix)
library(topicmodels)
library(Hmisc)
#library(lattice)

# Exécution lda :
#./lda est 0.1 20 settings.txt /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/H20_GH_sequences_lda_norep.txt random /Users/guilhemsommeria-klein/Desktop/These/Données_H20/algo_result1
#Exécution hdp-faster :
#./hdp --train_data /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/H20_GH_sequences_lda_norep.txt --directory /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/hdp-faster_alpha1_gamma1_eta0.01_maxiter100_samplehyper_norep_1/ --max_time -1 --max_iter 100 --save_lag 10 --verbose --sample_hyper yes/no
# Exécution hdp :
#./hdp --algorithm train --data /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/H20_GH_sequences_lda_norep.txt --directory /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/hdp_alpha1_gamma1_eta0.01_maxiter100_inittopic100_norep_1/ --max_iter 100 --save_lag 10 --eta 0.01 --split_merge yes/no --sample_hyper yes/no --init_topics 100

# to know where packages are stored : '.libPaths' in the R console, then access variable .Library

data_h20 = 0
data_pp = 1
barcode_gh = 0
barcode_itsfungi = 1
barcode_16s = 0

blei = 0
Rtopicmodels_Gibbs = 0
Rtopicmodels_VEM = 1

mpar = 0
nb_topics = 5

#nb_topics_range = c(2,3)
nb_topics_range = c(2,3,4,5,10,15,20,30,40,50)
#nb_topics_range = c(30,40,50,60,70)
mnb_topics = 0

if (mpar) {
  if (mnb_topics)
    mpar_range = length(nb_topics_range)
  #else if (miter)
    #mpar_range = length(miter_range)
} else mpar_range=1

delta = 0.1
alpha_insert = "50:nb_topics"

nb_real = 5
# if (mpar && mnb_topics) or ig (!mpar && best_keep), we must have best=0 
# best_keep = 1 has no influence if (mpar && mnb_topics)
best = 0
best_keep = 1

# only useful for Gibbs sampling
nb_iter = 10000
llh_keep = 100

# only useful for VEM
em_tol = 5*(10^-5)
var_tol = 10^-6

if (data_h20) {
  data_insert = "Données_H20"
  short_data_insert = "H20"
} else if (data_pp)
{
  data_insert = "Données_PetitPlateau"
  short_data_insert = "PP"
}
  
if (barcode_gh) {
  barcode_insert = "Plantes_GH"
  short_barcode_insert = "GH"
} else if (barcode_itsfungi)
{
  barcode_insert = "Champignons_ITS"
  short_barcode_insert = "ITSfungi"
} else if (barcode_16s)
{
  barcode_insert = "Bactéries_16S"
  short_barcode_insert = "16Sbact"
}

if (blei) {
  filename_insert = "Blei_LDA"
} else if (Rtopicmodels_Gibbs) {
  filename_insert = "Rtopicmodels_LDA_Gibbs"
} else if (Rtopicmodels_VEM)
  filename_insert = "Rtopicmodels_LDA_VEM"
  
# H20_GH.R (retrieve general information on the original data) :
###############
setwd(paste("/Users/guilhemsommeria-klein/Desktop/These/",data_insert,"/",barcode_insert,"/",sep=""))
taxo_ref <- read.table(paste(short_data_insert,"_",short_barcode_insert,"_taxo_ref_table.txt",sep=""), sep=" ")
# sample counts without replicates and controles - larger number of reads per sample than in _sequences_counts_norep.txt
# sample_counts <- as.vector(t(read.table(paste(short_data_insert,"_",short_barcode_insert,"_sample_counts.txt",sep=""))))

# Originally only computed in "H20_GH.R", now computed here too
# normal_ordered_seq <- as.vector(unlist(read.table("Seq_ordered_according_to_sitenormalized_abund.txt")))
# normal_data2m <- as.matrix(read.table("Site_composition_in_sequences.txt"))
# KL_normal_data2m1 <- as.matrix(read.table("Sequence_composition_in_sites.txt"))

data2 <- read.table(paste(short_data_insert,"_",short_barcode_insert,"_sequences_counts_norep.txt",sep=""), sep=" ")

# if (data_pp)
# {
# #Reminding the spatial positions (column indices) where a blank space needs to be introduced when plotting the spatial distibution of samples
# NEXT_CHAR = function(previous_last_char_index) {
#   if (previous_last_char_index == 29) {
#     next_last_char_index = 1
#   } else {
#     next_last_char_index = previous_last_char_index+1}
#   next_last_char_index
# }
# Sample_name_endings = c("7","8","9","A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")
# # Initialization
# first_miss = 1
# char_chain = as.character(data2[1,1])
# last_char = substr(char_chain,nchar(char_chain),nchar(char_chain))
# last_char_index = which(Sample_name_endings==last_char)
# spatial_index = 1
# if (last_char_index!=1)
# {
#   Missing_positions_indices = 1 
#   first_miss = 0 
# }
# previous_last_char = "7"
# previous_last_char_index = 1 
# # loop over non-empty samples
# for (j in 2:ncol(data2)) 
# {
#   spatial_index = spatial_index+1
#   char_chain = as.character(data2[1,j])
#   last_char = substr(char_chain,nchar(char_chain),nchar(char_chain))
#   last_char_index = which(Sample_name_endings==last_char)
#   if (last_char_index!=NEXT_CHAR(previous_last_char_index))
#   {
#     if (first_miss)
#     {
#       Missing_positions_indices = j
#       first_miss = 0
#     } else Missing_positions_indices = c(Missing_positions_indices,j)
#     iter_last_char_index = NEXT_CHAR(previous_last_char_index)
#     while (last_char_index != NEXT_CHAR(iter_last_char_index))
#     {
#       Missing_positions_indices = c(Missing_positions_indices,j)
#       iter_last_char_index = NEXT_CHAR(iter_last_char_index)
#     }
#   }
#   previous_last_char_index = last_char_index
# }
# }

if (data_pp)
{
  #Reminding the spatial positions (column indices) where a blank space needs to be introduced when plotting the spatial distibution of samples
  NEXT_CHAR = function(previous_last_char_index) {
    if (previous_last_char_index == 29) {
      next_last_char_index = 1
    } else {
      next_last_char_index = previous_last_char_index+1}
    next_last_char_index
  }
  Sample_name_endings = c("7","8","9","A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")
  # Initialization
  Missing_positions_indices = vector(length=29*39,mode="numeric")
  char_chain = as.character(data2[1,1])
  last_char = substr(char_chain,nchar(char_chain),nchar(char_chain))
  last_char_index = which(Sample_name_endings==last_char)
  spatial_index = 1
  if (last_char_index!=1)
    Missing_positions_indices[spatial_index] = 1 
  previous_last_char_index = 1 
  # loop over non-empty samples
  for (j in 2:ncol(data2)) 
  {
    spatial_index = spatial_index+1
    char_chain = as.character(data2[1,j])
    last_char = substr(char_chain,nchar(char_chain),nchar(char_chain))
    last_char_index = which(Sample_name_endings==last_char)
    if (last_char_index!=NEXT_CHAR(previous_last_char_index))
    {
      Missing_positions_indices[spatial_index] = 1
      spatial_index = spatial_index+1
      iter_last_char_index = NEXT_CHAR(previous_last_char_index)
      while (last_char_index != NEXT_CHAR(iter_last_char_index))
      {
        Missing_positions_indices[spatial_index] = 1
        iter_last_char_index = NEXT_CHAR(iter_last_char_index)
        spatial_index = spatial_index+1
      }
    }
    previous_last_char_index = last_char_index
  }
}

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
taxo_ref = taxo_ref[-(which(count2==0)+1),]
for (i in 1:nrow(data2m))
  taxo_ref[i+1,1] = i
rownames(taxo_ref) = seq(1,nrow(taxo_ref),1)
count2 = count2[-which(count2==0)]

# sample_count2 contains the number of reads per sample, in order to normalize the number of reads per sequence
sample_count2 = vector(length=ncol(data2m),mode="numeric")
for (j in 1:ncol(data2m))
{
  for (i in 1:nrow(data2m))
    sample_count2[j] = as.numeric(data2m[i,j])+sample_count2[j]
}
# # computing the variance of the number of reads per sample
# var_sample_count2 = 0
# for (i in 1:ncol(data2m))
#   var_sample_count2 = (sample_count2[i] - sum(sample_count2)/ncol(data2m))^2/ncol(data2m) + var_sample_count2
# std_dev_sample_count2 = sqrt(var_sample_count2)/(sum(sample_count2)/ncol(data2m))

# computing the read proportion per site for each sequence
normal_data2m = matrix(nrow=nrow(data2m),ncol=ncol(data2m),data=0)
# computing the read proportion per sequence for each site
# KL_normal_data2m and smoothed_KL_normal_data2m contains the proportion of the total number of reads belonging to sequence i that is found in sample j
# smoothed_KL2_normal_data2m contains the proportion of the reads in sample j that belong to sequence i
KL_normal_data2m = matrix(nrow=nrow(data2m),ncol=ncol(data2m),data=0)
smoothed_KL_normal_data2m = matrix(nrow=nrow(data2m),ncol=ncol(data2m),data=0)
smoothed_KL2_normal_data2m = matrix(nrow=nrow(data2m),ncol=ncol(data2m),data=0)
# normal_count2 contains the proportion of each sequence average over the sites
normal_count2 = vector(length=nrow(data2m),mode="numeric")
# nrow(data2m) = nb_seq et ncol(data2m) = nb_sites
for (i in 1:nrow(data2m))
{
  norm_factor_KL2 = 0
  for (j in 1:ncol(data2m))
  {
    if (count2[i]!=0)  
    {KL_normal_data2m[i,j] = as.numeric(data2m[i,j])/count2[i]}
    else
    {KL_normal_data2m[i,j] = NA}
    smoothed_KL_normal_data2m[i,j] = (as.numeric(data2m[i,j]) + 1)/(count2[i] + ncol(data2m))
    smoothed_KL2_normal_data2m[i,j] = (as.numeric(data2m[i,j])+1)/(sample_count2[j]+nrow(data2m))
    norm_factor_KL2 = norm_factor_KL2 + smoothed_KL2_normal_data2m[i,j]
    normal_data2m[i,j] = as.numeric(data2m[i,j])/sample_count2[j]
    normal_count2[i] = normal_count2[i] + normal_data2m[i,j]/ncol(data2m)
  }
  smoothed_KL2_normal_data2m[i,] = smoothed_KL2_normal_data2m[i,]/norm_factor_KL2
} 

sorted_normal_abundances2 = sort.int(normal_count2,decreasing=T,index.return=T)
normal_ordered_seq = sorted_normal_abundances2$ix
####################











LLH_final0 = matrix(nrow=2,ncol=mpar_range,data=0)
LLH_final1 = matrix(nrow=nb_real,ncol=mpar_range,data=0)
AIC0 = matrix(nrow=6,ncol=mpar_range,data=0)
AIC1 = matrix(nrow=nb_real,ncol=mpar_range,data=0)
AIC2 = matrix(nrow=nb_real,ncol=mpar_range,data=0)
AIC3 = matrix(nrow=nb_real,ncol=mpar_range,data=0)
if (Rtopicmodels_VEM) 
{
  alpha_est1 = matrix(nrow=nb_real,ncol=mpar_range,data=0)
  alpha_est0 = matrix(nrow=2,ncol=mpar_range,data=0)
}
  
# Start of the loop over parameters
#########################
for (par_index in 1:mpar_range)
{

cat("Parameter value",par_index,"\n")  
  
if (mpar) 
{
  if (mnb_topics)
    nb_topics = nb_topics_range[par_index]  
}
alpha = 50/nb_topics
  
if (blei)
{
  setwd(paste("/Users/guilhemsommeria-klein/Desktop/These/",data_insert,"/",barcode_insert,"/Blei_lda_alphaes",alpha,"_topics",nb_topics_range[i],"_norep_strset_3/"),sep="")
  
  #assignment <- read.table("word-assignments.dat", sep=" ")
  logbeta <- read.table("final.beta", sep=" ")
  logbeta = logbeta[,-c(1,2)]
  gamma <- read.table("final.gamma", sep=" ")
  alpha <- read.table("final.other", sep=" ")
  
  alpha_value = alpha[3,2]
  nb_terms = alpha[2,2]-1
  nb_doc = length(gamma[,1])
  
  documents = as.matrix(gamma-alpha_value)
  documents[documents<0] = 0
}

# beta = #topics x #terms
# gamma = # documents x #topics
else if (Rtopicmodels_Gibbs)
{
  SEED = vector(length=nb_real,mode="integer")
  for (j in 1:nb_real)
    SEED[j] = as.integer(Sys.time()) - j*10^7
  control_LDA_Gibbs = list(alpha=alpha, estimate.beta=TRUE,
                           verbose = 100, prefix = tempfile(), save = 0, keep = llh_keep,
                           seed = SEED, nstart = nb_real, best = best,
                           delta = delta,
                           iter = nb_iter, burnin = 0, thin = nb_iter)
  
  Result = topicmodels::LDA(x=t(data2m_num),k=nb_topics,method = "Gibbs",control=control_LDA_Gibbs,model=NULL)
}

else if (Rtopicmodels_VEM)
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
  
  Result = topicmodels::LDA(x=t(data2m_num),k=nb_topics,method = "VEM",control=control_LDA_VEM,model=NULL)
  if (par_index == 1)
    Result1 = Result
  else if (par_index == 2)
    Result2 = Result
  
  if (!best)
  {
    for (j in 1:nb_real)
    { 
      if (max(Result[[j]]@logLiks)>=0)
        cat("Error : positive loglikelihood value\n")
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
      LLH_final0[1,par_index] = Result[[j]]@fitted[[1]]@loglikelihood/nb_real + LLH_final0[1,par_index]  
      LLH_final1[j,par_index] = Result[[j]]@fitted[[1]]@loglikelihood
      LLH_final0[2,par_index] = Result[[j]]@fitted[[1]]@loglikelihood^2/nb_real + LLH_final0[2,par_index]  
      AIC0[2,par_index] = (2*(nb_topics*(nb_terms+nb_doc) - Result[[j]]@fitted[[1]]@loglikelihood))^2/nb_real + AIC0[2,par_index]
      AIC0[4,par_index] = (2*(nb_words - Result[[j]]@fitted[[1]]@loglikelihood))^2/nb_real + AIC0[4,par_index]
      AIC1[j,par_index] = 2*(nb_topics*(nb_terms+nb_doc) - Result[[j]]@fitted[[1]]@loglikelihood)
      AIC2[j,par_index] = 2*(nb_words - Result[[j]]@fitted[[1]]@loglikelihood)
    } 
    LLH_final0[2,par_index] = LLH_final0[2,par_index] - LLH_final0[1,par_index]^2
    LLH_final0[2,par_index] = sqrt(LLH_final0[2,par_index])
    AIC0[1,par_index] = 2*(nb_topics*(nb_terms+nb_doc) - LLH_final0[1,par_index])
    AIC0[3,par_index] = 2*(nb_words - LLH_final0[1,par_index])
    AIC0[2,par_index] = AIC0[2,par_index] - AIC0[1,par_index]^2
    AIC0[4,par_index] = AIC0[4,par_index] - AIC0[3,par_index]^2
    AIC0[2,par_index] = sqrt(AIC0[2,par_index])
    AIC0[4,par_index] = sqrt(AIC0[4,par_index])
  } else if (Rtopicmodels_VEM)
  {
    nb_terms = Result[[1]]@wordassignments$ncol
    nb_doc = Result[[1]]@wordassignments$nrow
    nb_words = Result[[1]]@n
    for (j in 1:nb_real)
    {
      LLH_final0[1,par_index] = sum(Result[[j]]@loglikelihood)/nb_real + LLH_final0[1,par_index]  
      LLH_final1[j,par_index] = sum(Result[[j]]@loglikelihood)
      LLH_final0[2,par_index] = (sum(Result[[j]]@loglikelihood))^2/nb_real + LLH_final0[2,par_index]  
      AIC0[2,par_index] = (2*(nb_topics*(nb_terms+nb_doc) - sum(Result[[j]]@loglikelihood)))^2/nb_real + AIC0[2,par_index]
      AIC0[4,par_index] = (2*(nb_words - sum(Result[[j]]@loglikelihood)))^2/nb_real + AIC0[4,par_index]
      AIC0[6,par_index] = (2*(nb_topics*nb_terms + 1 - sum(Result[[j]]@loglikelihood)))^2/nb_real + AIC0[6,par_index]
      AIC1[j,par_index] = 2*(nb_topics*(nb_terms+nb_doc) - sum(Result[[j]]@loglikelihood))
      AIC2[j,par_index] = 2*(nb_words - sum(Result[[j]]@loglikelihood))
      AIC3[j,par_index] = 2*(nb_topics*nb_terms + 1 - sum(Result[[j]]@loglikelihood))
      alpha_est1[j,par_index] = Result[[j]]@alpha
      alpha_est0[1,par_index] = Result[[j]]@alpha/nb_real + alpha_est0[1,par_index]
      alpha_est0[2,par_index] = Result[[j]]@alpha^2/nb_real + alpha_est0[2,par_index]
      #cat("\n",alpha_est0[1,i],"\n")
    } 
    LLH_final0[2,par_index] = LLH_final0[2,par_index] - LLH_final0[1,par_index]^2
    LLH_final0[2,par_index] = sqrt(LLH_final0[2,par_index])
    AIC0[1,par_index] = 2*(nb_topics*(nb_terms+nb_doc) - LLH_final0[1,par_index])
    AIC0[3,par_index] = 2*(nb_words - LLH_final0[1,par_index])
    AIC0[5,par_index] = 2*(nb_topics*nb_terms + 1 - LLH_final0[1,par_index])
    AIC0[2,par_index] = AIC0[2,par_index] - AIC0[1,par_index]^2
    AIC0[4,par_index] = AIC0[4,par_index] - AIC0[3,par_index]^2
    AIC0[6,par_index] = AIC0[6,par_index] - AIC0[5,par_index]^2
    AIC0[2,par_index] = sqrt(AIC0[2,par_index])
    AIC0[4,par_index] = sqrt(AIC0[4,par_index])
    AIC0[6,par_index] = sqrt(AIC0[6,par_index])
    alpha_est0[2,par_index] = alpha_est0[2,par_index] - alpha_est0[1,par_index]^2
    alpha_est0[2,par_index] = sqrt(alpha_est0[2,par_index])
  }
} else if (!mpar)
{
  if (best || best_keep)
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
          #LLH_final_real1[j] = Result[[j]]@logLiks[length(Result[[j]]@logLiks)]
        else if (Rtopicmodels_Gibbs)
          LLH_final_real1[j] = Result[[j]]@fitted[[1]]@loglikelihood
        #AIC1[j] = 2*(nb_topics*(nb_terms+nb_doc) - LLH_final_real1[j])
        #AIC2[j] = 2*(nb_words - LLH_final_real1[j])
      }
      for (j in 1:nb_real) 
      {
        Akaike_weights_llh[j] = exp(LLH_final_real1[j] - max(LLH_final_real1))/sum(exp(LLH_final_real1 - max(LLH_final_real1)))
        #Akaike_weights_AIC1[j] = exp(1/2*(min(AIC1) - AIC1[j]))/sum(exp(1/2*(min(AIC1) - AIC1)))
        #Akaike_weights_AIC2[j] = exp(1/2*(min(AIC2) - AIC2[j]))/sum(exp(1/2*(min(AIC2) - AIC2)))
      }
      if (Rtopicmodels_VEM)
      {
        best_real = which(LLH_final_real1 == max(LLH_final_real1))
        llh = LLH_final_real1[best_real] 
        llh_iterations = Result[[best_real]]@logLiks
        logbeta = Result[[best_real]]@beta 
        documents = Result[[best_real]]@gamma
        nb_terms = Result[[best_real]]@wordassignments$ncol
        nb_doc = Result[[best_real]]@wordassignments$nrow
        alpha_est = Result[[best_real]]@alpha
      } else if (Rtopicmodels_Gibbs)
      {
        best_real = which(LLH_final_real1 == max(LLH_final_real1))
        llh = LLH_final_real1[best_real] 
        if (llh_keep)
          llh_iterations = Result[[best_real]]@fitted[[1]]@logLiks 
        logbeta = Result[[best_real]]@fitted[[1]]@beta 
        documents = Result[[best_real]]@fitted[[1]]@gamma
        nb_terms = Result[[best_real]]@fitted[[1]]@wordassignments$ncol
        nb_doc = Result[[best_real]]@fitted[[1]]@wordassignments$nrow
      }
    }
      
    #### August 2014 code #####
    # nrow(data2m) = nb_seq et ncol(data2m) = nb_sites
    # nrow(documents) = nb_doc and ncol(documents) = nb_topics
    
    # propotion of each topic in a document (sums to 1 over topics)
    norm_documents = matrix(nrow=nb_doc,ncol=nb_topics)
    for (i in 1:nb_doc)
    {
      norm_documents[i,] = documents[i,]/sum(documents[i,])
    }
    
    # KL : propotion of each site/document in a topic (sums to 1 over sites/documents)
    # KL2 : propotion of each topic in a site/document (sums to 1 over topics)
    KL_norm_documents = matrix(nrow=nb_doc,ncol=nb_topics)
    smoothed_KL_norm_documents = matrix(nrow=nb_doc,ncol=nb_topics)
    smoothed_KL2_norm_documents = matrix(nrow=nb_doc,ncol=nb_topics)
    for (j in 1:nb_topics)
    {
      KL2_norm_factor = 0
      for (i in 1:nb_doc)
      {
        KL_norm_documents[i,j] = documents[i,j]/sum(documents[,j])
        smoothed_KL_norm_documents[i,j] = (documents[i,j]+1)/(sum(documents[,j])+nb_doc)
        smoothed_KL2_norm_documents[i,j] = (documents[i,j]+1)/(sum(documents[i,])+nb_topics)
        KL2_norm_factor = KL2_norm_factor + smoothed_KL2_norm_documents[i,j]
      }
      smoothed_KL2_norm_documents[,j]=smoothed_KL2_norm_documents[,j]/KL2_norm_factor
    }
    
    dominant_topic=vector(length=nb_doc,mode="integer")
    for (i in 1:nb_doc)
    {
      dominant_topic[i]=which(norm_documents[i,]==max(norm_documents[i,]))
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
      dom_topic=which(norm_documents[i,]==max(norm_documents[i,]))
      for (k in 1:5)
      {
        if (norm_documents[i,dom_topic]>=threshold[k])
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
          if (norm_documents[i,j]>=threshold[k])
            topic_site_nb[j,k]=topic_site_nb[j,k]+1
        }
      }
    }
    
    #tot_reads=sum(sample_counts)
    tot_reads=sum(sample_count2)
    # Computing the proportion of total read number for each topic
    prop_topic = vector(length=nb_topics,mode="numeric")
    for (j in 1:nb_topics)
    { 
      for (i in 1:nb_doc)
      {
        prop_topic[j] =  prop_topic[j] + sample_count2[i]*norm_documents[i,j]/tot_reads 
      }
    }
    sort_prop_topic = sort.int(prop_topic,index.return=T)
    
    # Computing the proportion of site-normalized read number for each topic
    normal_topic = vector(length=nb_topics,mode="numeric")
    for (j in 1:nb_topics)
    { 
      for (i in 1:nb_doc)
      {
        normal_topic[j] =  normal_topic[j] + norm_documents[i,j]/nb_doc 
      }
    }
    sort_normal_topic = sort.int(normal_topic,index.return=T)
    
    topic_compo=matrix(list(),nrow=1,ncol=nb_topics)
    for (k in 1:nb_topics)
    {
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
      topic_compo[[1,k]]=data.frame(proba,mat)
      #     for (i in 1:nb_terms)
      #     {
      #     #topic_compo_index[i]=taxo_ref[i+1,1]
      #     topic_compo[[1,k]][i,2:(length(taxo_ref[1,])+1)]=as.character(taxo_ref[vect$ix[i]+1,])
      #     }
      for (i in 1:nb_terms)
      {
        #topic_compo_index[i]=taxo_ref[i+1,1]
        for (j in 1:length(taxo_ref[1,]))
          topic_compo[[1,k]][i,j+1]=as.character(taxo_ref[vect$ix[i]+1,j])
      }
    }
    # topic_compo[[1,1]][1:20,]
    # composition of first topic by total read proportion ordening :
    # topic_compo[[1,rev(sort_prop_topic$ix)[1]]][1:10,]
    # topic_compo[[1,rev(sort_prop_topic$ix)[2]]][1:10,]
    # topic_compo[[1,rev(sort_prop_topic$ix)[3]]][1:10,]
    # topic_compo[[1,rev(sort_prop_topic$ix)[4]]][1:10,]
    # topic_compo[[1,rev(sort_prop_topic$ix)[5]]][1:10,]
    # topic_compo[[1,rev(sort_prop_topic$ix)[6]]][1:10,]
    # topic_compo[[1,rev(sort_prop_topic$ix)[7]]][1:10,]
    # topic_compo[[1,rev(sort_prop_topic$ix)[8]]][1:10,]
    # topic_compo[[1,rev(sort_prop_topic$ix)[9]]][1:10,]
    # topic_compo[[1,rev(sort_prop_topic$ix)[10]]][1:10,]
    # composition of first topic by site-normalized read proportion ordening :
    # topic_compo[[1,rev(sort_normal_topic$ix)[1]]][1:20,]
    # Ordening of topics by total read proportion :
    # rev(sort_prop_topic$ix)
    # Ordening of topics by site-normalized read proportion :
    # rev(sort_normal_topic$ix)
    
    } else if (!best)
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
          LLH0[1,] = Result[[j]]@fitted[[1]]@logLiks/nb_real + LLH0[1,]
          LLH0[2,] = Result[[j]]@fitted[[1]]@logLiks^2/nb_real + LLH0[2,]
          #llh[j] = Result[[j]]@fitted[[1]]@loglikelihood
        }
        LLH0[2,] = LLH0[2,] - LLH0[1,]^2
        LLH0[2,] = sqrt(LLH0[2,])
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
  
#### end of the !mpar condition
}

#### end of the loop over mpar_range
}
if (Rtopicmodels_VEM)
  if (mpar) {
    alpha_insert = paste("alpha",format(min(alpha_est0[1,]),digits=3),"-",format(max(alpha_est0[1,]),digits=3),sep="")
  } else if (!mpar) {
    if (!best)
      alpha_insert = paste("alpha",format(alpha_est0[1,],digits=3),sep="")
  }












#############################
# ######################### #
# ######## FIGURES ######## #
# ######################### #
#############################

if (!mpar)
{
  if (data_h20) {
    dirname = paste("/Users/guilhemsommeria-klein/Desktop/These/",data_insert,"/plots/",barcode_insert,"/",filename_insert,"/",sep="")
  } else if (data_pp)
    dirname = paste("/Users/guilhemsommeria-klein/Desktop/These/",data_insert,"/",barcode_insert,"/",filename_insert,"/",sep="")
  
  if (!(file.exists(dirname)))
  {dir.create(dirname)}
  setwd(dirname)  
  
  if (!best && !best_keep)
  {  
    
  if (Rtopicmodels_Gibbs)
  {
    filename = paste(filename_insert,"_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,".Rdata",sep="")
    save(Result,LLH0,LLH1,alpha_insert,delta,nb_topics,nb_real,nb_iter,llh_keep,file=filename)
    directory_file = "Directory.txt"
    write(paste(dirname,filename,sep=""),directory_file)
  } else if (Rtopicmodels_VEM) {
    filename = paste(filename_insert,"_",alpha_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,".Rdata",sep="")
    save(Result,LLH0,LLH1,alpha_insert,nb_topics,nb_real,file=filename)
    directory_file = "Directory.txt"
    write(paste(dirname,filename,sep=""),directory_file)
  }
  cat(filename)
  
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
   
  if (best) 
  {
    if (Rtopicmodels_Gibbs) {
      subdirname = paste(dirname,filename_insert,"_alpha",alpha,"_delta",delta,"_nb_topics",nb_topics,"_nb_iter",nb_iter,"_nb_real",nb_real,"_best/",sep="")
    } else if (Rtopicmodels_VEM) {
      subdirname = paste(dirname,filename_insert,"_alpha_est",format(alpha_est,digits=3),"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_best/",sep="")
    }
    if (!(file.exists(subdirname)))
    {dir.create(subdirname)}
    setwd(subdirname) 
    
    if (Rtopicmodels_Gibbs)
    {
      filename = paste(filename_insert,"_alpha",alpha,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,"_best.Rdata",sep="")
      save(Result,alpha,delta,logbeta,documents,nb_doc,nb_terms,nb_topics,nb_real,nb_iter,llh_keep,llh,
           topic_compo,sort_normal_topic,sort_prop_topic,topic_site_nb,dominant_topic_thres,spatial_dominant_topic,smoothed_KL2_norm_documents,norm_documents,file=filename)
      directory_file = "Directory.txt"
      write(paste(subdirname,filename,sep=""),directory_file)
    } else if (Rtopicmodels_VEM) {
      filename = paste(filename_insert,"_alpha_est",format(alpha_est,digits=3),"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_best.Rdata",sep="")
      save(Result,alpha_est,logbeta,documents,nb_doc,nb_terms,nb_topics,nb_real,llh,
           topic_compo,sort_normal_topic,sort_prop_topic,topic_site_nb,dominant_topic_thres,spatial_dominant_topic,smoothed_KL2_norm_documents,norm_documents,file=filename)
      directory_file = "Directory.txt"
      write(paste(subdirname,filename,sep=""),directory_file)
    }
    cat(filename)
  } else if (best_keep)
  {
    if (Rtopicmodels_Gibbs) {
      subdirname = paste(dirname,filename_insert,"_alpha",alpha,"_delta",delta,"_topics",nb_topics,"_nb_iter",nb_iter,"_nb_real",nb_real,"_best_keep/",sep="")
    } else if (Rtopicmodels_VEM) {
      subdirname = paste(dirname,filename_insert,"_alpha_est",format(alpha_est,digits=3),"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_best_keep/",sep="")
    }
    if (!(file.exists(subdirname)))
    {dir.create(subdirname)}
    setwd(subdirname) 
    
    if (Rtopicmodels_Gibbs)
    {
      filename = paste(filename_insert,"_alpha",alpha,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,"_best_keep.Rdata",sep="")
      save(Result,alpha,delta,logbeta,documents,nb_doc,nb_terms,nb_topics,nb_real,nb_iter,llh_keep,llh,AIC1,AIC2,LLH_final_real1,Akaike_weights_llh,
           topic_compo,sort_normal_topic,sort_prop_topic,topic_site_nb,dominant_topic_thres,spatial_dominant_topic,smoothed_KL2_norm_documents,norm_documents,file=filename)
      directory_file = "Directory.txt"
      write(paste(subdirname,filename,sep=""),directory_file)
    } else if (Rtopicmodels_VEM) {
      filename = paste(filename_insert,"_alpha_est",format(alpha_est,digits=3),"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_best_keep.Rdata",sep="")
      save(Result,alpha_est,logbeta,documents,nb_doc,nb_terms,nb_topics,nb_real,llh,AIC1,AIC2,LLH_final_real1,Akaike_weights_llh,
           topic_compo,sort_normal_topic,sort_prop_topic,topic_site_nb,dominant_topic_thres,spatial_dominant_topic,smoothed_KL2_norm_documents,norm_documents,file=filename)
      directory_file = "Directory.txt"
      write(paste(subdirname,filename,sep=""),directory_file)
    }
    cat(filename)
  }
  
# to print the 10 most abundant MOTUs in the most abundant topic, the second most abundant topic ... :
# topic_compo[[1,rev(sort_normal_topic$ix)[1]]][1:10,]
# topic_compo[[1,rev(sort_normal_topic$ix)[2]]][1:10,]

# Saving topic_compo in a txt file does not work for the moment :
# file = "Topic_composition.txt"
# write("test",file)
# for (k in 1:20)
# {
#   for (l in 1:(length(taxo_ref[1,])+1))
#   {
#     list_content = topic_compo[[1,rev(sort_normal_topic$ix)[1]]][k,l]
#     write(list_content,file, sep = "\t", append = TRUE)
#   }
#   write("\n",file, sep = "")
# }
# for (m in 2:nb_topics)
# {
#   write("",append = TRUE)
#   for (k in 1:20)
#   {
#     for (l in 1:(length(taxo_ref[1,])+1))
#     {
#       list_content = topic_compo[[1,rev(sort_normal_topic$ix)[m]]][k,l]
#       write(list_content,file, sep = "\t", append = TRUE)
#     }
#     write("\n",file, sep = "")
#   }
# }

if (best_keep)
{
  #####################
  plotname = paste("Realization_comparison_llh_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,"_em_tol",em_tol,".pdf",sep="")
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
  
  # the probability weights are the same whether they are calculated with the llh or with AIC, since all realizations have the same number of parameters
  #####################
  plotname = paste("Realization_comparison_AIC_weight_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,"_em_tol",em_tol,".pdf",sep="")
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
  plotname = paste("Realization_comparison_log10(AIC_weight)_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,"_em_tol",em_tol,".pdf",sep="")
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
  ######################
}
  
##########################

if (Rtopicmodels_VEM || (Rtopicmodels_Gibbs && llh_keep))
{
plotname = paste("Loglikelihood_convergence_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics,"_nb_real",nb_real,"_nb_iter",nb_iter,"_em_tol",em_tol,".pdf",sep="")
pdf(plotname)
#par(mar=c(5.1,4.1,4.1,2.1))
par(mar=c(5.1,5.1,4.1,2.1))

plot(llh_iterations,ann=F,cex.axis=1.5,lwd=2,type="p",col="black")

legend(x="bottomright",legend=c(filename_insert,alpha_insert,paste("delta =",delta),paste("nb topics =",nb_topics),paste("nb realizations =",nb_real)),col="black",inset=0)
title("Loglikelihood convergence",cex.main=1.7)
title(xlab="Iterations",ylab="Loglikelihood value",cex.lab=1.5)
dev.off()
}

########################################################

col_topic=c("black","blue","red","forestgreen","gold","mediumorchid2","chocolate4","cadetblue4","maroon","tan2","mediumturquoise","lightsalmon2","hotpink","greenyellow","lavender","khaki4")
# legend_topic=c("Topic 1","Topic 2","Topic 3","Topic 4","Topic 5")
# col_vector_rgb=matrix(ncol=3,nrow=nb_topics+1)
# col_vector_rgb[1,]=c(1,0,0)
# col_vector_rgb[2,]=c(0,1,0)
# col_vector_rgb[3,]=c(0,0,1)
# col_vector_rgb[4,]=c(1,0,1)
# col_vector_rgb[5,]=c(1,1,0)
# col_vector_rgb[6,]=c(0,1,1)
# col_vector_rgb[7,]=c(0.5,0,0)
# col_vector_rgb[8,]=c(0,0.5,0)
# col_vector_rgb[9,]=c(0,0,0.5)
# col_vector_rgb[10,]=c(0.5,0,0.5)
# col_vector_rgb[11,]=c(1,1,1)
# last_element=length(col_vector_rgb[,1])

###################################
#
###############################   #
# Topic abundance information #   #
###############################   #
#
###################################
subsubdirname = paste(subdirname,"topics_abundance_info/",sep="")
if (!(file.exists(subsubdirname)))
  {dir.create(subsubdirname)}
setwd(subsubdirname) 

# pdf("Global_topic_proportions.pdf")
# plot(1:nb_topics,prop_topic,ann=FALSE)
# title(xlab="Topic index",ylab="Global topic proportion",main="Global proportion of each topic")
# dev.off()  

pdf("Sorted_global_topic_proportions.pdf")
plot(nb_topics:1,sort_prop_topic$x,ann=FALSE,ylim=c(0,max(prop_topic)))
lines(c(1,nb_topics),c(0,0),type="l",lty=2)
title(ylab="Global topic proportion",main="Global proportion of each topic")
dev.off()  

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

pdf("Sorted_global_topic_log_proportions.pdf")
plot(nb_topics:1,log(sort_prop_topic$x),ann=FALSE)
#lines(c(1,nb_topics),c(1,1),type="l",lty=2)
title(ylab="Global topic proportion (log)",main="Global proportion of each topic (log)")
dev.off()  

######################################################
#
##################################################   #
# Information about topics' sequence composition #   #
##################################################   #
#
######################################################
setwd(subdirname)
subsubdirname = paste(subdirname,"topics_sequence_composition_info/",sep="")
if (!(file.exists(subsubdirname)))
  dir.create(subsubdirname)
setwd(subsubdirname) 

nb_bins = 10
max_topic = vector(length=nb_topics,mode="numeric")
for (j in 1:nb_topics)
  max_topic[j] = max(exp(logbeta)[j,]) 

pdf("topic_dominance_by_sequence_hist.pdf")
hist(max_topic, breaks = nb_bins, freq=T, xlab = "Topic's proportion taken by the most abundant sequence in topic",
     ylab = "Number of topics", main = "Number of topics\n with respect to the topic's proportion\n taken by the most abundant sequence in topic", xaxp = c(0, 1, 10))
#ylab = "Number of sites", main = "Histogramme of the number of sites\n with respect to the proportion\n of the most abundant sequence")
dev.off()

# pdf("topic_dominance_by_sequence_hist_percent.pdf")
# h = hist(max_topic, breaks = nb_bins, plot=F)
# h$density = h$counts/sum(h$counts)
# plot(h, freq=F, xlab = "Proportion of the most abundant sequence",
#      ylab = "Proportion of the total number of topics", 
#      main = "Proportion of the total number of topics\n with respect to the topic's proportion\n
#      taken by the most abundant sequence in topic", xaxp = c(0, 1, 10))
# #ylab = "Number of sites", main = "Proportion of sites\n with respect to the proportion\n of the most abundant sequence")
# dev.off()

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
  # prop_reads_seq_main_topic[j] = sorted_seq$x[1]*prop_topic[sorted_seq$ix[1]]*tot_reads/count_total_seq
}

pdf("topic_exclusivity_hist.pdf")
hist(prop_reads_seq_main_topic, breaks = nb_bins, freq=T, xlab = "Share of the sequence's total number of reads\n in the topic where the sequence is most dominant",
     ylab = "Number of sequences", main = "Number of sequences\n with respect to the share of their total number of reads\n in the topic where they are most dominant", xaxp = c(0, 1, 10))
dev.off()

# pdf("topic_exclusivity_hist_percent.pdf")
# h = hist(prop_reads_seq_main_topic, breaks = nb_bins, plot=F)
# h$density = h$counts/sum(h$counts)
# plot(h, freq=F, xlab = "Share of the sequence's total number of reads\n in the topic where the sequence is most dominant",
#      ylab = "Proportion of sequences", main = "Proportion of the number of sequences\n with respect to the share of their 
#      total number of reads\n in the topic where they are most dominant", xaxp = c(0, 1, 10))
# dev.off()

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

# pdf("Prop_in_main_topic_vs_total_read_share_in_main_topic.pdf")
# par(mar=c(5.1,4.1,4.1,4.1))
# plot(prop_reads_seq_main_topic[-c(242,586,591,702)],prop_max_topic_seq[-c(242,586,591,702)],type="p", yaxt="n",
#      main="Proportion of a sequence in the topic where it is most dominant\n with respect to this topic's share of the sequence's total number of reads",
#      xlab = "Share of the sequence's total number of reads\n in the topic where the sequence is most dominant",
#      ylab = "Proportion of the sequence in the topic where it is most dominant")
# axis(2, ylim=range(prop_max_topic_seq), col='black')
# par(new=T)
# plot(prop_reads_seq_main_topic[-c(242,586,591,702)],log(prop_max_topic_seq[-c(242,586,591,702)])/log(10),type="p",ann=F,yaxt="n",xaxt="n",col="blue")
# axis(4, ylim=range(log(prop_max_topic_seq)/log(10)), col='blue')
# mtext("Proportion of the sequence in the topic where it is most dominant (log_10)",side=4,line=3,col="blue")
# dev.off()
# 
# # figure utilisant le classement "sorted_normal_abundances2$ix" des séquences selon leur abondance site-normalized, obtenu en sortie de H20_GH.R
# pdf(paste("Prop_in_main_topic_vs_total_read_share_in_main_topic_",nb_topics,"firstseq.pdf",sep=""))
# par(mar=c(5.1,4.1,4.1,4.1))
# plot(prop_reads_seq_main_topic[normal_ordered_seq[1:nb_topics]],prop_max_topic_seq[normal_ordered_seq[1:nb_topics]],type="p", yaxt="n",
#      main="Proportion of a sequence in the topic where it is most dominant\n with respect to this topic's share of the sequence's total number of reads",
#      xlab = "Share of the sequence's total number of reads\n in the topic where the sequence is most dominant",
#      ylab = "Proportion of the sequence in the topic where it is most dominant")
# axis(2, ylim=range(prop_max_topic_seq), col='black')
# # par(new=T)
# # plot(prop_reads_seq_main_topic[normal_ordered_seq[1:nb_topics]],log(prop_max_topic_seq[normal_ordered_seq[1:nb_topics]])/log(10),type="p",ann=F,yaxt="n",xaxt="n",col="blue")
# # axis(4, ylim=range(log(prop_max_topic_seq)/log(10)), col='blue')
# # mtext("Proportion of the sequence in the topic where it is most dominant (log_10)",side=4,line=3,col="blue")
# dev.off()
# 
# # figure utilisant le classement "sorted_normal_abundances2$ix" des séquences selon leur abondance site-normalized, obtenu en sortie de H20_GH.R
# pdf(paste("Prop_in_main_topic_vs_total_read_share_in_main_topic_",nb_topics*2,"firstseq.pdf",sep=""))
# par(mar=c(5.1,4.1,4.1,4.1))
# plot(prop_reads_seq_main_topic[normal_ordered_seq[1:nb_topics*2]],prop_max_topic_seq[normal_ordered_seq[1:nb_topics*2]],type="p", yaxt="n",
#      main="Proportion of a sequence in the topic where it is most dominant\n with respect to this topic's share of the sequence's total number of reads",
#      xlab = "Share of the sequence's total number of reads\n in the topic where the sequence is most dominant",
#      ylab = "Proportion of the sequence in the topic where it is most dominant")
# axis(2, ylim=range(prop_max_topic_seq), col='black')
# par(new=T)
# plot(prop_reads_seq_main_topic[normal_ordered_seq[1:nb_topics*2]],log(prop_max_topic_seq[normal_ordered_seq[1:nb_topics*2]])/log(10),type="p",ann=F,yaxt="n",xaxt="n",col="blue")
# axis(4, ylim=range(log(prop_max_topic_seq)/log(10)), col='blue')
# mtext("Proportion of the sequence in the topic where it is most dominant (log_10)",side=4,line=3,col="blue")
# dev.off()

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

pdf("Site-normalized_MOTU_composition_for_each_topic_20firstseq_assignPrecision.pdf")
par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
#par(mar=c(5.1,4.1,4.1,2.1)
#bottom left top right
par(mar=c(14.1,7.1,4.1,4.1))
for (k in 1:nb_topics)
{
  plot(topic_compo[[1,rev(sort_normal_topic$ix)[k]]][1:20,1],type="p",ann=F,yaxt="n",xaxt="n")
  axis(2, ylim=range(topic_compo[[1,rev(sort_normal_topic$ix)[k]]][1:20,1]), col='black')
  axis(1, at=1:20, labels = F)
  labels = vector(length=20,mode="character")
  for (i in 1:20)
    labels[i]=paste(topic_compo[[1,rev(sort_normal_topic$ix)[k]]][i,11],topic_compo[[1,rev(sort_normal_topic$ix)[k]]][i,10],topic_compo[[1,rev(sort_normal_topic$ix)[k]]][i,2],sep=" - ")
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

# KL symmetrised distance between topic's site repartition and sequence site's repartition (result from H20_GH.R)
#############################

KL = matrix(nrow=nb_topics,ncol=nb_terms,data=0)
smoothed_KL = matrix(nrow=nb_topics,ncol=nb_terms,data=0)
smoothed_KL2 = matrix(nrow=nb_topics,ncol=nb_terms,data=0)
# KL = matrix(nrow=nb_topics,ncol=nb_topics,data=0)
# smoothed_KL = matrix(nrow=nb_topics,ncol=nb_topics,data=0)
# smoothed_KL2 = matrix(nrow=nb_topics,ncol=nb_topics,data=0)
for (topic in 1:nb_topics)
{
  # on enlève les quatre derniers séquences qui ont des counts nuls et KL_normal_data2m=NA
  # for (seq in 1:(nb_terms-4))
  for (seq in 1:(10*nb_topics))
  {
    for (site in 1:nb_doc)
    {
      if (KL_norm_documents[site,rev(sort_normal_topic$ix)[topic]]!=0 && !is.na(KL_normal_data2m[normal_ordered_seq[seq],site]))
        #if (KL_norm_documents[site,rev(sort_normal_topic$ix)[topic]]!=0)
        #if (KL_normal_data2m[normal_ordered_seq[seq],site]!=0)
      {
        # KL[topic,seq]=1/(2*log(2))*(norm_documents[site,rev(sort_normal_topic$ix)[topic]]*(log(norm_documents[site,rev(sort_normal_topic$ix)[topic]])-log(normal_data2m[normal_ordered_seq[seq],site]))
        # + normal_data2m[normal_ordered_seq[seq],site]*(log(normal_data2m[normal_ordered_seq[seq],site])-log(norm_documents[site,rev(sort_normal_topic$ix)[topic]]))) 
        KL[topic,seq]=1/(2*log(2))*(KL_norm_documents[site,rev(sort_normal_topic$ix)[topic]]*log(KL_norm_documents[site,rev(sort_normal_topic$ix)[topic]]/KL_normal_data2m[normal_ordered_seq[seq],site])
                                    + KL_normal_data2m[normal_ordered_seq[seq],site]*log(KL_normal_data2m[normal_ordered_seq[seq],site]/KL_norm_documents[site,rev(sort_normal_topic$ix)[topic]])) + KL[topic,seq] 
      }
      smoothed_KL[topic,seq] = (1/(2*log(2))*(smoothed_KL_norm_documents[site,rev(sort_normal_topic$ix)[topic]]*log(smoothed_KL_norm_documents[site,rev(sort_normal_topic$ix)[topic]]/smoothed_KL_normal_data2m[normal_ordered_seq[seq],site])
                                              + smoothed_KL_normal_data2m[normal_ordered_seq[seq],site]*log(smoothed_KL_normal_data2m[normal_ordered_seq[seq],site]/smoothed_KL_norm_documents[site,rev(sort_normal_topic$ix)[topic]])) 
                                + smoothed_KL[topic,seq]) 
      smoothed_KL2[topic,seq] = (1/(2*log(2))*(smoothed_KL2_norm_documents[site,rev(sort_normal_topic$ix)[topic]]*log(smoothed_KL2_norm_documents[site,rev(sort_normal_topic$ix)[topic]]/smoothed_KL2_normal_data2m[normal_ordered_seq[seq],site])
                                               + smoothed_KL2_normal_data2m[normal_ordered_seq[seq],site]*log(smoothed_KL2_normal_data2m[normal_ordered_seq[seq],site]/smoothed_KL2_norm_documents[site,rev(sort_normal_topic$ix)[topic]])) 
                                 + smoothed_KL2[topic,seq]) 
    }
  }
}  

# pdf("KL_distance_between_topic_and_sequence_site_repartitions.pdf")
# color2D.matplot(KL[,1:nb_topics],c(1,0),c(0,0),c(0,1),
#                 extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,
#                 ylab="Topics ranked by site-normalized abundances",
#                 xlab="Sequences ranked by\n site-normalized abundances",
#                 main="KL symmetrized divergence between topics' site distribution\n and sequences' site distribution -\n topics and sequences ranked by site-averaged relative abundances",
#                 #main="KL symmetrized distance between topics' site repartition\n and sequences' site repartition -\n topics and sequences ranked by site-normalized abundances",
#                 do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
#                 border="black",na.color=NA)
# dev.off()
# 
# pdf("KL_distance_between_smoothed_topic_and_sequence_site_repartitions.pdf")
# color2D.matplot(smoothed_KL[,1:nb_topics],c(1,0),c(0,0),c(0,1),
#                 extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,
#                 ylab="Topics ranked by site-normalized abundances",
#                 xlab="Sequences ranked by\n site-normalized abundances",
#                 main="KL symmetrized divergence between smoothed topics' site distribution\n and smoothed sequences' site distribution -\n topics and sequences ranked by site-averaged relative abundances",
#                 do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
#                 border="black",na.color=NA)
# dev.off()

# pdf(paste("KL_distance_between_topic_and_sequence_site_repartitions_",nb_topics,"firstseq.pdf",sep=""))
# color2D.matplot(KL[1:nb_topics,1:nb_topics],c(1,0),c(0,0),c(0,1),
#                 extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,
#                 ylab="Topics ranked by site-normalized abundances",
#                 xlab="Sequences ranked by\n site-normalized abundances",
#                 main="KL symmetrized divergence between topics' site distribution\n and sequences' site distribution -\n topics and sequences ranked by site-averaged relative abundances",
#                 #main="KL symmetrized distance between topics' site repartition\n and sequences' site repartition -\n topics and sequences ranked by site-normalized abundances",
#                 do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
#                 border="black",na.color=NA)
# dev.off()

# pdf("KL_distance_(log)_between_topic_and_sequence_site_repartitions_20firstseq.pdf")
# color2D.matplot(log(KL[1:20,1:20]),c(1,0),c(0,0),c(0,1),
#                 extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,
#                 ylab="Topics ranked by site-normalized abundances",
#                 xlab="Sequences ranked by\n site-normalized abundances",
#                 main="KL symmetrized divergence (log) between topics' site distribution\n and sequences' site distribution -\n topics and sequences ranked by site-averaged relative abundances",
#                 #main="KL symmetrized distance (log) between topics' site repartition\n and sequences' site repartition -\n topics and sequences ranked by site-normalized abundances",
#                 do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
#                 border="black",na.color=NA)
# dev.off()

# KL[which(KL[,1:nb_topics]==0)]=NA
# pdf("KL_distance_(log)_between_topic_and_sequence_site_repartitions.pdf")
# color2D.matplot(log(KL[,1:nb_topics]),c(1,0),c(0,0),c(0,1),
#                 extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,
#                 ylab="Topics ranked by site-normalized abundances",
#                 xlab="Sequences ranked by\n site-normalized abundances",
#                 main="KL symmetrized divergence (log) between topics' site distribution\n and sequences' site distribution -\n topics and sequences ranked by site-averaged relative abundances",
#                 #main="KL symmetrized distance (log) between topics' site repartition\n and sequences' site repartition -\n topics and sequences ranked by site-normalized abundances",
#                 do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
#                 border="black",na.color=NA)
# dev.off()

# pdf("KL2_distance_between_smoothed_topic_and_sequence_site_repartitions.pdf")
# #par(mar=c(5.1,4.1,4.1,2.1))
# par(mar=c(6.1,6.1,5.1,2.1))
# color2D.matplot(smoothed_KL2[,1:nb_topics],c(1,0),c(0,0),c(0,1),
#                 extremes=NA,cellcolors=NA,show.legend=F,nslices=10,
#                 ylab="Classes ranked by\n relative abundances averaged over samples",xlab="",
#                 main="Symmetrized KL divergence\n between smoothed topics' spatial distribution\n and smoothed sequences' spatial distribution",
#                 #main="KL symmetrized distance between smoothed topics' site repartition\n and smoothed sequences' site repartition -\n topics and sequences ranked by site-normalized abundances",
#                 do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2,
#                 border="black",na.color=NA)
# mtext("Sequences ranked by\n relative abundances averaged over samples",side=1,line=4,cex=1.5) 
# dev.off()
# 
# pdf("KL2_distance_between_smoothed_topic_and_sequence_site_repartitions_seq=2nb_topics.pdf")
# color2D.matplot(smoothed_KL2[,1:(2*nb_topics)],c(1,0),c(0,0),c(0,1),
#                 extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,
#                 ylab="Topics ranked by site-normalized abundances",
#                 xlab="Sequences ranked by\n site-normalized abundances",
#                 main="KL symmetrized divergence between smoothed topics' site distribution\n and smoothed sequences' site distribution -\n classes and sequences ranked by site-averaged relative abundances",
#                 #main="KL symmetrized distance between smoothed topics' site repartition\n and smoothed sequences' site repartition -\n topics and sequences ranked by site-normalized abundances",
#                 do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
#                 border="black",na.color=NA)
# dev.off()

pdf("KL2_distance_between_smoothed_topic_and_norm-sorted_sequence_site_repartitions_seq=3nb_topics.pdf")
par(mar=c(6.1,6.1,5.1,2.1))
color2D.matplot(smoothed_KL2[,1:(3*nb_topics)],c(1,0),c(0,0),c(0,1),
                extremes=NA,cellcolors=NA,show.legend=F,nslices=10,
                ylab="Classes ranked by\n relative abundances averaged over samples",xlab="",
                main="Symmetrized KL divergence\n between smoothed classes' spatial distribution\n and smoothed sequences' spatial distribution",
                #main="KL symmetrized distance between smoothed topics' site repartition\n and smoothed sequences' site repartition -\n topics and sequences ranked by site-normalized abundances",
                do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2,
                border="black",na.color=NA)
mtext("Sequences ranked by\n relative abundances averaged over samples",side=1,line=4,cex=1.5) 
dev.off()

# pdf("KL2_distance_between_smoothed_topic_and_prop-sorted_sequence_site_repartitions_seq=3nb_topics.pdf")
# par(mar=c(6.1,6.1,5.1,2.1))
# color2D.matplot(smoothed_KL2[,1:(3*nb_topics)],c(1,0),c(0,0),c(0,1),
#                 extremes=NA,cellcolors=NA,show.legend=F,nslices=10,
#                 ylab="Classes ranked by\n relative abundances averaged over samples",xlab="",
#                 main="Symmetrized KL divergence\n between smoothed classes' spatial distribution\n and smoothed sequences' spatial distribution",
#                 #main="KL symmetrized distance between smoothed topics' site repartition\n and smoothed sequences' site repartition -\n topics and sequences ranked by site-normalized abundances",
#                 do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2,
#                 border="black",na.color=NA)
# mtext("Sequences ranked by\n relative abundances averaged over samples",side=1,line=4,cex=1.5) 
# dev.off()

##################################################
#
##############################################   #
# Information about topics' site repartition #   #
##############################################   #
#
##################################################
setwd(subdirname)
subsubdirname = paste(subdirname,"/topics_site_repartition_info/",sep="")
if (!(file.exists(subsubdirname)))
  dir.create(subsubdirname)
setwd(subsubdirname) 

nb_bins = 10
topic_dominance_by_site = vector(length = nb_doc, mode = "numeric")
for (i in 1:nb_doc)
  topic_dominance_by_site[i] = max(norm_documents[i,])

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

# Figure illisible
# #################################
# pdf("Site_topic_composition.pdf")
# plot(1,norm_documents[1,1],type="n",ann=FALSE,xlim=c(1,nb_doc),ylim=range(norm_documents))
# for (i in 1:nb_topics)
#   lines(1:nb_doc,norm_documents[,i],type="p",col=col_topic[i])
# #lines(range1,c(log(1/8),log(1/8)),type="l",col="black",lty=2)
# #axis(2,at=log(1/8),lab="log(1/8)",col.axis="black",col="black")
# title(xlab="Sample number",ylab="Topic proportion")
# legend(x="topright",legend=legend_topic[1:nb_topics],text.col=col_topic[1:nb_topics],inset=0.1)
# title("Topic composition per site, replicates counted as sites")
# dev.off()
##############################

# The following two figures need manual assignment of colors to topics
#################################
# pdf("Site_dominant_topic_map.pdf")
# # color2D.matplot(spatial_dominant_topic,c(0,1),c(0,1),c(1,0),
# #                   extremes=NA,cellcolors=NA,show.legend=FALSE,nslices=10,xlab="Column",
# #                   ylab="Row",do.hex=FALSE,axes=TRUE,show.values=TRUE,vcol="white",vcex=1,
# #                   border="black",na.color=NA,main="Dominant topic map")
# cellcolors <- matrix(nrow=19,ncol=19)
# # Loop over topics (one picturing the dominant topic for each site)
# for (k in 1:nb_topics)
#   {
#   #cellcolors[spatial_dominant_topic==k] = rgb(k/(nb_topics+1),0,1-k/(nb_topics+1))
#   cellcolors[spatial_dominant_topic==k] = rgb(col_vector_rgb[k,1],col_vector_rgb[k,2],col_vector_rgb[k,3])
# #color.scale(spatial_dominant_topic[spatial_dominant_topic==k],cs1=c(),cs2=c(),cs3=())  
#   }
#   color2D.matplot(spatial_dominant_topic,
#   extremes=NA,cellcolors=cellcolors,show.legend=FALSE,nslices=nb_topics,xlab="Column",
#   ylab="Row",do.hex=FALSE,axes=TRUE,show.values=TRUE,vcol="white",vcex=1,
#   border="black",na.color=NA,main="Dominant topic map")
# dev.off()
# 
# pdf("Site_dominant_topic_threshold_maps.pdf")
# par(mfrow=c(2,2))
# #lattice::levelplot(abund2.pred~x+y, z2, 
# #col.regions=topo.colors(100), aspect = "iso",contour=T,main=binomial)
# # Loop over threshold values (one map for each threshold value)
# for (k in 2:5)
#   {
#   spatial_dominant_topic_thres = matrix(nrow=19,ncol=19)
#   for (j in 1:19)
#     {
#     for (i in 1:19)
#       {
#       spatial_dominant_topic_thres[i,j] = dominant_topic_thres[(j-1)*19+i,k]
#       }
#     }
#   cellcolors <- matrix(nrow=19,ncol=19)
#   # Loop over topics (on each map, the dominant topic of each site is shown, except when its proportion is below the threshold, 
#   # in which case the site's color is white)
#   for (l in 1:nb_topics)
#     {
#     #cellcolors[spatial_dominant_topic==k] = rgb(k/(nb_topics+1),0,1-k/(nb_topics+1))
#     cellcolors[spatial_dominant_topic_thres==l] = rgb(col_vector_rgb[l,1],col_vector_rgb[l,2],col_vector_rgb[l,3])
#     #color.scale(spatial_dominant_topic[spatial_dominant_topic==k],cs1=c(),cs2=c(),cs3=())  
#     }
#   cellcolors[spatial_dominant_topic==0] = rgb(col_vector_rgb[last_element,1],col_vector_rgb[last_element,2],col_vector_rgb[last_element,3])  
#   # colors Red Green Blue
#   color2D.matplot(spatial_dominant_topic_thres,
#                   extremes=NA,cellcolors=cellcolors,show.legend=FALSE,nslices=10,xlab="Column",
#                   ylab="Row",do.hex=FALSE,axes=TRUE,show.values=TRUE,vcol="white",vcex=1,
#                   border="black",na.color=NA,main=paste("Threshold ",threshold[k]))
#   }
# dev.off()

# # cellcolors<-matrix(NA,nrow=10,ncol=10)
# # cellcolors[corr.matrix >= 0]<-
# #   color.scale(corr.matrix[corr.matrix >= 0],
# #               cs1=c(0.7,0),cs2=c(0,0.7),cs3=0)
# # cellcolors[corr.matrix < 0]<-
# #   color.scale(corr.matrix[corr.matrix < 0],
# #               cs1=c(0.7,0),cs2=0,cs3=0.7)
# # color2D.matplot(corr.matrix,cellcolors=cellcolors,
# #                 show.values=TRUE)   
##############################

# Topic composition maps with default topic ordering
#######################################
# pdf("Site_topic_composition_maps.pdf")
# # 2 topics :
# #par(mfrow=c(2,1))
# # 5 topics :
# #par(mfrow=c(3,2))
# # 10 topics
# par(mfrow=c(2,2))
# #lattice::levelplot(abund2.pred~x+y, z2, 
#                #col.regions=topo.colors(100), aspect = "iso",contour=T,main=binomial)
# # Loop over topics (one map per topic, the color stands for the proportion of the topic)
# for (k in 1:nb_topics)
#   {  
#   spatial_topicmix = matrix(nrow=19,ncol=19)
#   for (j in 1:19)
#     {
#     for (i in 1:19)
#       {
#       spatial_topicmix[i,j] = norm_documents[(j-1)*19+i,k]
#       }
#     }
#   # colors Red Green Blue
#   color2D.matplot(spatial_topicmix,c(0,1),c(0,0),c(1,0),
#                 extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,xlab="Column",
#                 ylab="Row",do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
#                 border="black",na.color=NA,main=paste("Topic ",k))
#   }
# dev.off()


############################
pdf("Topic_ordered_by_total_abundance_composition_maps.pdf")
# 2 topics :
#par(mfrow=c(2,1))
# 5 topics :
#par(mfrow=c(3,2))
# 10 topics
par(mfrow=c(2,2))
#lattice::levelplot(abund2.pred~x+y, z2, 
#col.regions=topo.colors(100), aspect = "iso",contour=T,main=binomial)
# Loop over topics (one map per topic, the color stands for the proportion of the topic)
for (k in 1:nb_topics)
{ 
  if (data_h20)
  {
    spatial_topicmix = matrix(nrow=19,ncol=19) 
    for (j in 1:19)
    {
      for (i in 1:19)
      {
        spatial_topicmix[i,j] = norm_documents[(j-1)*19+i,rev(sort_prop_topic$ix)[k]]
      }
    }
  } else if (data_pp)
  { 
    spatial_topicmix = matrix(nrow=29,ncol=39)
    position_shift = 0
    for (j in 1:39)
    {
      for (i in 1:29)
      {
        if (Missing_positions_indices[(j-1)*29+i]==0)
          spatial_topicmix[i,j] = norm_documents[(j-1)*29+i-position_shift,rev(sort_prop_topic$ix)[k]]    
        else if (Missing_positions_indices[(j-1)*29+i]==1)
        {
          spatial_topicmix[i,j] = 0
          position_shift = position_shift+1
        } 
      }
    }
  }
  # colors Red Green Blue
  color2D.matplot(spatial_topicmix,c(0,1),c(0,0),c(1,0),
                  extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,xlab="Column",
                  ylab="Row",do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
                  border="black",na.color=NA,main=paste("Topic ",k," - ",rev(sort_prop_topic$ix)[k],sep=""))
}
dev.off()

# pdf("First_topic_ordered_by_site-normalized_abundance_composition_map.pdf")
# # 2 topics :
# #par(mfrow=c(2,1))
# # 5 topics :
# #par(mfrow=c(3,2))
# # 10 topics
# #par(mfrow=c(2,2))
# #lattice::levelplot(abund2.pred~x+y, z2, 
# #col.regions=topo.colors(100), aspect = "iso",contour=T,main=binomial)
# # Loop over topics (one map per topic, the color stands for the proportion of the topic)
# spatial_topicmix = matrix(nrow=19,ncol=19)
# for (j in 1:19)
# {
#   for (i in 1:19)
#   {
#     spatial_topicmix[i,j] = norm_documents[(j-1)*19+i,rev(sort_normal_topic$ix)[1]]
#   }
# }
# # colors Red Green Blue
# color2D.matplot(spatial_topicmix,c(0,1),c(0,0),c(1,0),
#                 extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,xlab="Column",
#                 ylab="Row",do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2,
#                 border="black",na.color=NA,main="Spatial distribution of the most abundant class")
# dev.off()

pdf("Topic_ordered_by_site-normalized_abundance_composition_maps.pdf")
# 2 topics :
#par(mfrow=c(2,1))
# 5 topics :
#par(mfrow=c(3,2))
# 10 topics
par(mfrow=c(2,2))
#lattice::levelplot(abund2.pred~x+y, z2, 
#col.regions=topo.colors(100), aspect = "iso",contour=T,main=binomial)
# Loop over topics (one map per topic, the color stands for the proportion of the topic)
for (k in 1:nb_topics)
{
  if (data_h20)
  {
    spatial_topicmix = matrix(nrow=19,ncol=19) 
    for (j in 1:19)
    {
      for (i in 1:19)
      {
        spatial_topicmix[i,j] = norm_documents[(j-1)*19+i,rev(sort_normal_topic$ix)[k]]
      }
    }
  } else if (data_pp)
  { 
    spatial_topicmix = matrix(nrow=29,ncol=39)
    position_shift = 0
    for (j in 1:39)
    {
      for (i in 1:29)
      {
        if (Missing_positions_indices[(j-1)*29+i]==0)
          spatial_topicmix[i,j] = norm_documents[(j-1)*29+i-position_shift,rev(sort_normal_topic$ix)[k]]    
        else if (Missing_positions_indices[(j-1)*29+i]==1)
        {
          spatial_topicmix[i,j] = 0
          position_shift = position_shift+1
        } 
      }
    }
  }
  # colors Red Green Blue
  color2D.matplot(spatial_topicmix,c(0,1),c(0,0),c(1,0),
                  extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,xlab="Column",
                  ylab="Row",do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
                  border="black",na.color=NA,main=paste("Topic ",rev(sort_normal_topic$ix)[k]," - #",k,sep=""))
}
dev.off()

pdf("Smoothed_plot-normalized_topic_proportion_ordered_by_site-normalized_abundance_composition_maps.pdf")
# 2 topics :
#par(mfrow=c(2,1))
# 5 topics :
#par(mfrow=c(3,2))
# 10 topics
par(mfrow=c(2,2))
#lattice::levelplot(abund2.pred~x+y, z2, 
#col.regions=topo.colors(100), aspect = "iso",contour=T,main=binomial)
# Loop over topics (one map per topic, the color stands for the proportion of the topic)
for (k in 1:nb_topics)
{  
  if (data_h20)
  {
    spatial_topicmix = matrix(nrow=19,ncol=19) 
    for (j in 1:19)
    {
      for (i in 1:19)
      {
        spatial_topicmix[i,j] = smoothed_KL_norm_documents[(j-1)*19+i,rev(sort_normal_topic$ix)[k]]
      }
    }
  } else if (data_pp)
  { 
    spatial_topicmix = matrix(nrow=29,ncol=39)
    position_shift = 0
    for (j in 1:39)
    {
      for (i in 1:29)
      {
        if (Missing_positions_indices[(j-1)*29+i]==0)
          spatial_topicmix[i,j] = smoothed_KL_norm_documents[(j-1)*29+i-position_shift,rev(sort_normal_topic$ix)[k]]    
        else if (Missing_positions_indices[(j-1)*29+i]==1)
        {
          spatial_topicmix[i,j] = 0
          position_shift = position_shift+1
        }
      }
    }
  }
  # colors Red Green Blue
  color2D.matplot(spatial_topicmix,c(0,1),c(0,0),c(1,0),
                  extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,xlab="Column",
                  ylab="Row",do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
                  border="black",na.color=NA,main=paste("Topic ",rev(sort_normal_topic$ix)[k]," - #",k,sep=""))
}
dev.off()

#   plot(normal_abundances_taxo2[1:20,1],type="p",ann=F,yaxt="n",xaxt="n")
# axis(2, ylim=range(normal_abundances_taxo2[,1]), col='black')
# # for (i in 1:20) 
# #   mtext("aaaaaaaaaa",at=i,side=1,line=1)
# axis(1, at=1:20, labels = F)
# labels = vector(length=20,mode="character")
# for (i in 1:20)
#   labels[i]=paste(normal_abundances_taxo2[i,11],normal_abundances_taxo2[i,10],sep=" - ")
# # text(1:20, par("usr")[3] + 0.25, srt = 45, adj = c(1,1), labels = abundances_taxo2[1:20,10], xpd = TRUE)
# text(1:20, par("usr")[3], srt = 45, adj = c(1.05,1.75), labels = labels, xpd = TRUE, cex=1.1)
# #par(new=T)
# #plot(log(normal_abundances_taxo2[1:20,1])/log(10),type="p",col="blue",ann=F,yaxt="n",xaxt="n")
# #axis(4, ylim=range(log(normal_abundances_taxo2[,1])), col='blue')
# #mtext("Relative abundance\n averaged over samples (log_10)",side=4,line=4,col="blue",cex=1.5) 
# title(ylab="Relative abundance\n averaged over samples")
# title("Relative abundances averaged over samples\n for the 20 most abundant sequences")
# dev.off()

############################

# The three following figures need the results of other runs to be plotted
###############################################

# pdf("Sorted_global_topic_proportions_3runs_lines.pdf")
# plot(nb_topics:1,sort(prop_topic),ann=FALSE,ylim=c(0,max(prop_topic,prop_topic1,prop_topic2)),type="l")
# lines(nb_topics:1,sort(prop_topic1),col="blue",type="l")
# lines(nb_topics:1,sort(prop_topic2),col="forestgreen",type="l")
# lines(c(1,nb_topics),c(0,0),type="l",lty=2)
# title(ylab="Global topic proportion",main="Global proportion of each topic")
# dev.off()  

# pdf("Sorted_global_topic_proportions_4runs_40topics.pdf")
# plot(nb_topics:1,sort(prop_topic),ann=FALSE,ylim=c(0,max(prop_topic,prop_topic1,prop_topic2,prop_topic3)),type="p")
# lines(nb_topics:1,sort(prop_topic1),col="blue",type="p")
# lines(nb_topics:1,sort(prop_topic2),col="forestgreen",type="p")
# lines(nb_topics:1,sort(prop_topic3),col="red",type="p")
# lines(c(1,nb_topics),c(0,0),type="l",lty=2)
# title(ylab="Global topic proportion",main="Global proportion of each topic")
# dev.off()  
# 
# pdf("Sorted_global_topic_log_proportions_4runs_40topics.pdf")
# plot(nb_topics:1,log(sort(prop_topic)),ann=FALSE,ylim=range(log(prop_topic),log(prop_topic1),log(prop_topic2),log(prop_topic3)))
# lines(nb_topics:1,log(sort(prop_topic1)),col="blue",type="p")
# lines(nb_topics:1,log(sort(prop_topic2)),col="forestgreen",type="p")
# lines(nb_topics:1,log(sort(prop_topic3)),col="red",type="p")
# #lines(c(1,nb_topics),c(1,1),type="l",lty=2)
# title(ylab="Global topic proportion (log)",main="Global proportion of each topic (log)")
# dev.off()  

###############################################

pdf("Number_of_dominant_sites_for_each_topic.pdf")
plot(nb_topics:1,topic_site_nb[sort_prop_topic$ix,1],ann=F,type="p",ylim=range(topic_site_nb))
for (k in 2:5)
  lines(nb_topics:1,topic_site_nb[sort_prop_topic$ix,k],type="p",col=col_topic[k])
title(ylab="Number of dominant sites",xlab="Topics ordered according to their global read proportion",main="Number of dominant sites for each topic")
legend(x="topleft",legend=c("Threshold 0.5","Threshold 0.6","Threshold 0.7","Threshold 0.8","Threshold 0.9"),text.col=col_topic[1:5],inset=0.1)
dev.off()

# The following figures need the results of other runs to be plotted
###############################################
# Comparison of the number of dominant sites for each topic and for different runs

# pdf("Number_of_dominant_sites_for_each_topic_readSort_70pc_lines.pdf")
# plot(nb_topics:1,topic_site_nb[sort_prop_topic$ix,3],ann=F,type="n",ylim=range(topic_site_nb[,3],topic_site_nb1[,3],topic_site_nb2[,3],topic_site_nb3[,3]))
# lines(nb_topics:1,topic_site_nb[sort_prop_topic$ix,3],type="l",col=col_topic[1])
# lines(nb_topics:1,topic_site_nb1[sort_prop_topic1$ix,3],type="l",col=col_topic[2])
# lines(nb_topics:1,topic_site_nb2[sort_prop_topic2$ix,3],type="l",col=col_topic[3])
# lines(nb_topics:1,topic_site_nb3[sort_prop_topic3$ix,3],type="l",col=col_topic[4])
# title(ylab="Number of dominant sites",xlab="Topics ordered according to their global read proportion",main="Number of dominant sites for each topic, threshold = 70%")
# dev.off()
# 
# pdf("Number_of_dominant_sites_for_each_topic_readSort_80pc_lines.pdf")
# plot(nb_topics:1,topic_site_nb[sort_prop_topic$ix,4],ann=F,type="n",ylim=range(topic_site_nb[,4],topic_site_nb1[,4],topic_site_nb2[,4],topic_site_nb3[,4]))
# lines(nb_topics:1,topic_site_nb[sort_prop_topic$ix,4],type="l",col=col_topic[1])
# lines(nb_topics:1,topic_site_nb1[sort_prop_topic1$ix,4],type="l",col=col_topic[2])
# lines(nb_topics:1,topic_site_nb2[sort_prop_topic2$ix,4],type="l",col=col_topic[3])
# lines(nb_topics:1,topic_site_nb3[sort_prop_topic3$ix,4],type="l",col=col_topic[4])
# title(ylab="Number of dominant sites",xlab="Topics ordered according to their global read proportion",main="Number of dominant sites for each topic, threshold = 80%")
# dev.off()
# 
# pdf("Number_of_dominant_sites_for_each_topic_readSort_60pc.pdf")
# plot(nb_topics:1,topic_site_nb[sort_prop_topic$ix,2],ann=F,type="n",ylim=range(topic_site_nb[,2],topic_site_nb1[,2],topic_site_nb2[,2],topic_site_nb3[,2]))
# lines(nb_topics:1,topic_site_nb[sort_prop_topic$ix,2],type="p",col=col_topic[1])
# lines(nb_topics:1,topic_site_nb1[sort_prop_topic1$ix,2],type="p",col=col_topic[2])
# lines(nb_topics:1,topic_site_nb2[sort_prop_topic2$ix,2],type="p",col=col_topic[3])
# lines(nb_topics:1,topic_site_nb3[sort_prop_topic3$ix,2],type="p",col=col_topic[4])
# title(ylab="Number of dominant sites",xlab="Topics ordered according to their global read proportion",main="Number of dominant sites for each topic, threshold = 60%")
# dev.off()
#
# pdf("Number_of_dominant_sites_for_each_topic_70pc_lines.pdf")
# plot(nb_topics:1,sort(topic_site_nb[,3]),ann=F,type="n",ylim=range(topic_site_nb[,3],topic_site_nb1[,3],topic_site_nb2[,3],topic_site_nb3[,3]))
# lines(nb_topics:1,sort(topic_site_nb[,3]),type="l",col=col_topic[1])
# lines(nb_topics:1,sort(topic_site_nb1[,3]),type="l",col=col_topic[2])
# lines(nb_topics:1,sort(topic_site_nb2[,3]),type="l",col=col_topic[3])
# lines(nb_topics:1,sort(topic_site_nb3[,3]),type="l",col=col_topic[4])
# title(ylab="Number of dominant sites",main="Number of dominant sites for each topic, threshold = 70%")
# dev.off()
# 
# pdf("Number_of_dominant_sites_for_each_topic_80pc_lines.pdf")
# plot(nb_topics:1,sort(topic_site_nb[,4]),ann=F,type="n",ylim=range(topic_site_nb[,4],topic_site_nb1[,4],topic_site_nb2[,4],topic_site_nb3[,4]))
# lines(nb_topics:1,sort(topic_site_nb[,4]),type="l",col=col_topic[1])
# lines(nb_topics:1,sort(topic_site_nb1[,4]),type="l",col=col_topic[2])
# lines(nb_topics:1,sort(topic_site_nb2[,4]),type="l",col=col_topic[3])
# lines(nb_topics:1,sort(topic_site_nb3[,4]),type="l",col=col_topic[4])
# title(ylab="Number of dominant sites",main="Number of dominant sites for each topic, threshold = 80%")
# dev.off()

# end of (best || best_keep) condition
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
  dirname = paste("/Users/guilhemsommeria-klein/Desktop/These/",data_insert,"/plots/",barcode_insert,"/",filename_insert,"/",sep="")
} else if (data_pp)
  dirname = paste("/Users/guilhemsommeria-klein/Desktop/These/",data_insert,"/",barcode_insert,"/",filename_insert,"/",sep="")
  
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
  save(Result,LLH_final0,AIC0,nb_topics_range,alpha_insert,delta,nb_real,nb_iter,file=filename)
  directory_file = "Directory.txt"
  write(paste(subdirname,filename,sep=""),directory_file)
} else if (Rtopicmodels_VEM)  {
  filename = paste(filename_insert,"_",alpha_insert,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,".Rdata",sep="")
  save(Result,LLH_final0,AIC0,nb_topics_range,alpha_insert,nb_real,alpha_est0,file=filename)
  directory_file = "Directory.txt"
  write(paste(subdirname,filename,sep=""),directory_file)
}
cat(filename)

########
plotname = paste("llh_compare_",alpha_insert,"_delta",delta,"_nb_topics",nb_topics_range[1],"-",nb_topics_range[length(nb_topics_range)],"_nb_real",nb_real,"_nb_iter",nb_iter,".pdf",sep="")
pdf(plotname)
#par(mar=c(5.1,4.1,4.1,2.1))
par(mar=c(5.1,5.1,4.1,2.1))

errbar(nb_topics_range,LLH_final0[1,],LLH_final0[1,]+LLH_final0[2,],LLH_final0[1,]-LLH_final0[2,],ann=F,cex.axis=1.5,lwd=2,type="p",col="black")

legend(x="bottomright",legend=c(filename_insert,alpha_insert,paste("delta =",delta),paste("nb realizations =",nb_real),paste("nb iterations =",nb_iter)),col="black",inset=0)
title("Final log-likelihood value vs number of topics",cex.main=1.7)
title(xlab="Number of topics",ylab="Final log-likelihood value",cex.lab=1.5)
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

legend(x="bottomright",legend=c(filename_insert,alpha_insert,paste("delta =",delta),paste("nb realizations =",nb_real),paste("nb iterations =",nb_iter)),col="black",inset=0)
title("Final log-likelihood value vs number of topics",cex.main=1.7)
title(xlab="Number of topics",ylab="Final log-likelihood value",cex.lab=1.5)
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

plot(nb_topics_range,AIC3[1,],ann=F,cex.axis=1.5,lwd=2,type="p",col="black",ylim=range(AIC1))
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





