# Required to use the function rdirichlet :
library(gtools)

nb_samples=29*39
nb_topics=5
nb_motus=1000
nb_reads=1000

nonmixed = 0
discretemix = 0
continuousmix = 1
randommix = 0
dirichletmix = 0

randomtopics = 0
dirichlettopics = 1

lnor_noise = 1

dirichletalpha = 0.02
sig_noise = 2

if (nonmixed)
{
  sample_compo = "Non-mixed"
} else if (discretemix)
{
  sample_compo = "Discrete-mixed"
  #sample_compo = "Discrete-mixed_equimol"
} else if (continuousmix)
{
  sample_compo = "Continuous-mixed"
  #sample_compo = "Discrete-mixed_equimol"
} else if (randommix)
{
  sample_compo = "Random-mixed"
  #sample_compo = "Discrete-mixed_equimol"
}
  
if (randomtopics)
{
  topic_compo = "randomtopics"
} else if (dirichlettopics)
{
  topic_compo = paste("dirichlettopics",dirichletalpha,sep="")
}  

if (lnor_noise)
{
  lnor_noise_insert = paste("_lnornoise_sig",sig_noise,sep="")
} else if (!lnor_noise)
{
  lnor_noise_insert = ""
}

testdata_dir = paste(sample_compo,"_samples_nbtopics",nb_topics,"_nbmotus",nb_motus,"_",topic_compo,"_",nb_reads,"sampledreads",lnor_noise_insert,"/",sep="")
dirname = paste("/Users/guilhemsommeria-klein/Desktop/These/Donnees_PetitPlateau/Test_data/",testdata_dir,sep="")
if (!(file.exists(dirname)))
  dir.create(dirname)
setwd(dirname) 

Sample_composition = matrix(nrow=nb_topics,ncol=nb_samples,data=0)
if (nonmixed)
{
  # 3 non-mixed communities
  for (j in 1:nb_samples)
  {
    if (j <= 10*39)
      Sample_composition[,j] = c(1,0,0)
    else if (j <= 20*39)
      Sample_composition[,j] = c(0,1,0)
    else 
      Sample_composition[,j] = c(0,0,1)
  }
} else if (discretemix)
{
  # 3 mixed communities
#   for (j in 1:nb_samples)
#     Sample_composition[,j] = c(1/3,1/3,1/3)
    
  for (j in 1:nb_samples)
  {
    if (j <= 10*39)
      Sample_composition[,j] = c(0.8,0.1,0.1)
    else if (j <= 20*39)
      Sample_composition[,j] = c(0.1,0.8,0.1)
    else 
      Sample_composition[,j] = c(0.1,0.1,0.8)
  }
} else if (randommix)
{
  for (j in 1:nb_samples)
  {
    Draw = sort(runif(nb_topics-1),decreasing=F)
    Sample_composition[1,j] = Draw[1]
    for (k in 2:(nb_topics-1))
    {
      Sample_composition[k,j] = Draw[k]-Draw[k-1]
    } 
    Sample_composition[nb_topics,j] = 1-Draw[nb_topics-1]
  }
} else if (dirichletmix)
{
  for (j in 1:nb_samples)
  {
    Sample_composition[,j] = rdirichlet(1,rep(0.5,nb_topics)) 
  }
} else if (continuousmix)
{
  if (nb_topics == 3)
  {
    for (j in 1:nb_samples)
    {
      if (j <= (nb_samples-1)/2)
        Sample_composition[,j] = c((1-cos(2*pi*j/nb_samples))/2,(cos(2*pi*j/nb_samples)+1)/2,0)
      else if (j == (nb_samples+1)/2)
        Sample_composition[,j] = c(1,0,0)
      else if (j > (nb_samples+1)/2)
        Sample_composition[,j] = c((1-cos(2*pi*j/nb_samples))/2,0,(cos(2*pi*j/nb_samples)+1)/2)
    }
  } else if (nb_topics == 5)
  {
    for (j in 1:nb_samples)
    {
      if (j <= (nb_samples-3)/4)
        Sample_composition[,j] = c((1-cos(4*pi*j/nb_samples))/2,(cos(4*pi*j/nb_samples)+1)/2,0,0,0)
      else if (j == ((nb_samples-3)/4+1))
        Sample_composition[,j] = c(1,0,0,0,0)
      else if ((j > ((nb_samples-3)/4+1)) && (j <= (nb_samples-1)/2))
        Sample_composition[,j] = c((1-cos(4*pi*j/nb_samples))/2,0,(cos(4*pi*j/nb_samples)+1)/2,0,0)
      else if (j == ((nb_samples-1)/2+1))
        Sample_composition[,j] = c(0,0,1,0,0)
      else if ((j > ((nb_samples-1)/2+1)) && (j <= (3*(nb_samples-3)/4+2)))
        Sample_composition[,j] = c(0,0,(cos(4*pi*j/nb_samples)+1)/2,(1-cos(4*pi*j/nb_samples))/2,0)
      else if (j == (3*(nb_samples-3)/4+3))
        Sample_composition[,j] = c(0,0,0,1,0)
      else if (j > (3*(nb_samples-3)/4+3))
        Sample_composition[,j] = c(0,0,0,(1-cos(4*pi*j/nb_samples))/2,(cos(4*pi*j/nb_samples)+1)/2)
    }
  }
}
  
true_norm_documents = t(Sample_composition)

Topic_composition = matrix(nrow=nb_motus,ncol=nb_topics,data=0)
if (randomtopics)
{
  for (j in 1:nb_topics)
  {
    Draw = sort(runif(nb_motus-1),decreasing=F)
    Topic_composition[1,j] = Draw[1]
    for (k in 2:(nb_motus-1))
    {
      Topic_composition[k,j] = Draw[k]-Draw[k-1]
    } 
    Topic_composition[nb_motus,j] = 1-Draw[nb_motus-1]
  }
} else if (dirichlettopics)
{
  for (j in 1:nb_topics)
    #Topic_composition[,j] = rmultinom(1, 10000, rdirichlet(1,rep(dirichletalpha,nb_motus)))
    Topic_composition[,j] = rdirichlet(1,rep(dirichletalpha,nb_motus))
}

Read_composition = matrix(nrow=nb_motus,ncol=nb_samples,data=0)
for (i in 1:nb_motus)
{
  for (j in 1:nb_samples)
  {
    for (k in 1:nb_topics)
      Read_composition[i,j] = Topic_composition[i,k]*Sample_composition[k,j] + Read_composition[i,j]
  }
}

if (lnor_noise)
{
  # Adding noise on abundances
  for (j in 1:nb_samples)
  {
    C = rlnorm(nb_motus,0,sig_noise)
    Read_composition[,j] = C*Read_composition[,j]/sum(C*Read_composition[,j])
  }
}

#data2m_num = apply(Read_composition*10000,MARGIN=2,FUN=round)
# Sampling the reads from a multinomial distribution
data2m_num = matrix(nrow=nb_motus,ncol=nb_samples,data=0)
for (j in 1:nb_samples)
  data2m_num[,j] = rmultinom(1, nb_reads, Read_composition[,j])
data2m = apply(data2m_num,MARGIN=2,FUN=as.character)

nb_empty_motus = 0
for (i in 1:nb_motus)
{
  if (sum(data2m_num[i,]) == 0) 
    nb_empty_motus = nb_empty_motus + 1
}
write(nb_empty_motus,"nb_empty_motus.txt")

taxo_ref = matrix(nrow = nb_motus+1, ncol = 11)
taxo_ref = apply(taxo_ref,MARGIN=2,FUN=as.character)
taxo_ref = data.frame(taxo_ref,stringsAsFactors=F)
taxo_ref[1,1] = "Sequence_number"
taxo_ref[2:(nb_motus+1),1] = as.character(seq(1,nb_motus,1))

# vector filled with zeros, since no missing sample
Missing_positions_indices = vector(length=29*39,mode="numeric")

directory_file = "Directory.txt"
write(dirname,directory_file)
save(data2m,file="data2m_testdata.Rdata")
save(taxo_ref,file="taxo_ref_testdata.Rdata")
save(Missing_positions_indices,file="Missing_positions_indices.Rdata")
save(true_norm_documents,file="true_documents.Rdata")


  