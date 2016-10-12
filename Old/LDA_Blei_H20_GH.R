library(plotrix)
library(topicmodels)
#library(lattice)

# Exécution lda :
#./lda est 0.1 20 settings.txt /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/H20_GH_sequences_lda_norep.txt random /Users/guilhemsommeria-klein/Desktop/These/Données_H20/algo_result1
#Exécution hdp-faster :
#./hdp --train_data /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/H20_GH_sequences_lda_norep.txt --directory /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/hdp-faster_alpha1_gamma1_eta0.01_maxiter100_samplehyper_norep_1/ --max_time -1 --max_iter 100 --save_lag 10 --verbose --sample_hyper yes/no
# Exécution hdp :
#./hdp --algorithm train --data /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/H20_GH_sequences_lda_norep.txt --directory /Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/hdp_alpha1_gamma1_eta0.01_maxiter100_inittopic100_norep_1/ --max_iter 100 --save_lag 10 --eta 0.01 --split_merge yes/no --sample_hyper yes/no --init_topics 100

# to know where packages are stored : '.libPaths' in the R console, then access variable .Library

blei = 0
topicmodels = 1

nb_topics = 5

setwd("/Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/")
taxo_ref <- read.table("H20_GH_taxo_ref_table.txt", sep=" ")
# sample counts without replicates and controles
sample_counts <- as.vector(t(read.table("H20_GH_sample_counts.txt")))

normal_ordered_seq <- as.vector(unlist(read.table("Seq_ordered_according_to_sitenormalized_abund.txt")))
# normal_data2m <- as.matrix(read.table("Site_composition_in_sequences.txt"))
# KL_normal_data2m1 <- as.matrix(read.table("Sequence_composition_in_sites.txt"))

# H20_GH.R (retrieve general information on the original data) :
###############
data2 <- read.table("H20_GH_sequences_counts_norep.txt", sep=" ")
data2m = as.matrix(data2[-1,])
# count2 contains the total count fo each sequence
count2 = vector(length=nrow(data2m),mode="numeric")
for (i in 1:nrow(data2m))
{
  for (j in 1:ncol(data2m))
  {
    count2[i]=as.numeric(data2m[i,j])+count2[i]
  }  
}

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
####################

if (blei)
{
  setwd(paste("/Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/Blei_lda_alphaes",alpha,"_topics",nb_topics,"_norep_strset_3/"))

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

if (topicmodels)
{
  control_LDA_Gibbs = list(alpha=50/nb_topics, estimate.beta=TRUE,
                            verbose = 0, prefix = tempfile(), save = 0, keep = 0,
                            seed = as.integer(Sys.time()), nstart = 1, best = TRUE,
                            delta = 0.1,
                            iter = 2000, burnin = 0, thin = 2000)
  
  Result = topicmodels::LDA(x=t(as.integer(data2m)),k=nb_topics,method = "Gibbs",control=control_LDA_Gibbs,model=NULL)
  
}

# propotion of each topic in a document (sums to 1 over topics)
norm_documents = matrix(nrow=nb_doc,ncol=nb_topics)
for (i in 1:nb_doc)
  {
  norm_documents[i,] = documents[i,]/sum(documents[i,])
  }

# propotion of each site/document in a topic (sums to 1 over sites/documents)
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

tot_reads=sum(sample_counts)
# Computing the proportion of total read number for each topic
prop_topic = vector(length=nb_topics,mode="numeric")
for (j in 1:nb_topics)
{ 
  for (i in 1:nb_doc)
  {
    prop_topic[j] =  prop_topic[j] + sample_counts[i]*norm_documents[i,j]/tot_reads 
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
    taxo_ref_names=c(taxo_ref_names,as.character(taxo_ref[1,i]))
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

########################################################
setwd("/Users/guilhemsommeria-klein/Desktop/These/Données_H20/plots/Plantes_GH/Blei_lda_alphaes0.1_topics40_norep_strset_3/")

# file = "Topic_composition.txt"
# write(topic_compo[[1,1]][1:20,],file,ncolumns = length(taxo_ref[1,])+1, sep = "\t")
# for (k in 2:nb_topics)
#   {
#   write("",append = TRUE)
#   write(topic_compo[[1,k]][1:20,],file,ncolumns = length(taxo_ref[1,])+1,append = TRUE, sep = "\t")
#   }

########################################################

col_topic=c("black","blue","red","forestgreen","gold","mediumorchid2","chocolate4","cadetblue4","maroon","tan2","mediumturquoise","lightsalmon2","hotpink","greenyellow","lavender","khaki4")
legend_topic=c("Topic 1","Topic 2","Topic 3","Topic 4","Topic 5")
col_vector_rgb=matrix(ncol=3,nrow=nb_topics+1)
col_vector_rgb[1,]=c(1,0,0)
col_vector_rgb[2,]=c(0,1,0)
col_vector_rgb[3,]=c(0,0,1)
col_vector_rgb[4,]=c(1,0,1)
col_vector_rgb[5,]=c(1,1,0)
col_vector_rgb[6,]=c(0,1,1)
col_vector_rgb[7,]=c(0.5,0,0)
col_vector_rgb[8,]=c(0,0.5,0)
col_vector_rgb[9,]=c(0,0,0.5)
col_vector_rgb[10,]=c(0.5,0,0.5)
col_vector_rgb[11,]=c(1,1,1)
last_element=length(col_vector_rgb[,1])

###################################
                                  #
###############################   #
# Topic abundance information #   #
###############################   #
                                  #
###################################
setwd("topics_abundance_info/")

pdf("Global_topic_proportions.pdf")
plot(1:nb_topics,prop_topic,ann=FALSE)
title(xlab="Topic index",ylab="Global topic proportion",main="Global proportion of each topic")
dev.off()  

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
setwd("../topics_sequence_composition_info/")

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
plot(max_topic[sort_normal_topic$ix],
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

pdf("Prop_in_main_topic_vs_total_read_share_in_main_topic.pdf")
par(mar=c(5.1,4.1,4.1,4.1))
plot(prop_reads_seq_main_topic[-c(242,586,591,702)],prop_max_topic_seq[-c(242,586,591,702)],type="p", yaxt="n",
     main="Proportion of a sequence in the topic where it is most dominant\n with respect to this topic's share of the sequence's total number of reads",
     xlab = "Share of the sequence's total number of reads\n in the topic where the sequence is most dominant",
     ylab = "Proportion of the sequence in the topic where it is most dominant")
axis(2, ylim=range(prop_max_topic_seq), col='black')
par(new=T)
plot(prop_reads_seq_main_topic[-c(242,586,591,702)],log(prop_max_topic_seq[-c(242,586,591,702)])/log(10),type="p",ann=F,yaxt="n",xaxt="n",col="blue")
axis(4, ylim=range(log(prop_max_topic_seq)/log(10)), col='blue')
mtext("Proportion of the sequence in the topic where it is most dominant (log_10)",side=4,line=3,col="blue")
dev.off()

# figure utilisant le classement "sorted_normal_abundances2$ix" des séquences selon leur abondance site-normalized, obtenu en sortie de H20_GH.R
pdf(paste("Prop_in_main_topic_vs_total_read_share_in_main_topic_",nb_topics,"firstseq.pdf",sep=""))
par(mar=c(5.1,4.1,4.1,4.1))
plot(prop_reads_seq_main_topic[normal_ordered_seq[1:nb_topics]],prop_max_topic_seq[normal_ordered_seq[1:nb_topics]],type="p", yaxt="n",
     main="Proportion of a sequence in the topic where it is most dominant\n with respect to this topic's share of the sequence's total number of reads",
     xlab = "Share of the sequence's total number of reads\n in the topic where the sequence is most dominant",
     ylab = "Proportion of the sequence in the topic where it is most dominant")
axis(2, ylim=range(prop_max_topic_seq), col='black')
# par(new=T)
# plot(prop_reads_seq_main_topic[normal_ordered_seq[1:nb_topics]],log(prop_max_topic_seq[normal_ordered_seq[1:nb_topics]])/log(10),type="p",ann=F,yaxt="n",xaxt="n",col="blue")
# axis(4, ylim=range(log(prop_max_topic_seq)/log(10)), col='blue')
# mtext("Proportion of the sequence in the topic where it is most dominant (log_10)",side=4,line=3,col="blue")
dev.off()

# figure utilisant le classement "sorted_normal_abundances2$ix" des séquences selon leur abondance site-normalized, obtenu en sortie de H20_GH.R
pdf(paste("Prop_in_main_topic_vs_total_read_share_in_main_topic_",nb_topics*2,"firstseq.pdf",sep=""))
par(mar=c(5.1,4.1,4.1,4.1))
plot(prop_reads_seq_main_topic[normal_ordered_seq[1:nb_topics*2]],prop_max_topic_seq[normal_ordered_seq[1:nb_topics*2]],type="p", yaxt="n",
     main="Proportion of a sequence in the topic where it is most dominant\n with respect to this topic's share of the sequence's total number of reads",
     xlab = "Share of the sequence's total number of reads\n in the topic where the sequence is most dominant",
     ylab = "Proportion of the sequence in the topic where it is most dominant")
axis(2, ylim=range(prop_max_topic_seq), col='black')
par(new=T)
plot(prop_reads_seq_main_topic[normal_ordered_seq[1:nb_topics*2]],log(prop_max_topic_seq[normal_ordered_seq[1:nb_topics*2]])/log(10),type="p",ann=F,yaxt="n",xaxt="n",col="blue")
axis(4, ylim=range(log(prop_max_topic_seq)/log(10)), col='blue')
mtext("Proportion of the sequence in the topic where it is most dominant (log_10)",side=4,line=3,col="blue")
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

##################################################
                                                 #
##############################################   #
# Information about topics' site repartition #   #
##############################################   #
                                                 #
##################################################
setwd("../topics_site_repartition_info/")

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
  for (seq in 1:(nb_terms-4))
  #for (seq in 1:nb_topics)
    {
    for (site in 1:nb_doc)
      {
      if (KL_norm_documents[site,rev(sort_normal_topic$ix)[topic]]!=0 && KL_normal_data2m[normal_ordered_seq[seq],site]!=0)
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

pdf("KL_distance_between_topic_and_sequence_site_repartitions.pdf")
color2D.matplot(KL[,1:nb_topics],c(1,0),c(0,0),c(0,1),
               extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,
               ylab="Topics ranked by site-normalized abundances",
               xlab="Sequences ranked by\n site-normalized abundances",
               main="KL symmetrized divergence between topics' site distribution\n and sequences' site distribution -\n topics and sequences ranked by site-averaged relative abundances",
               #main="KL symmetrized distance between topics' site repartition\n and sequences' site repartition -\n topics and sequences ranked by site-normalized abundances",
               do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
               border="black",na.color=NA)
dev.off()

pdf("KL_distance_between_smoothed_topic_and_sequence_site_repartitions.pdf")
color2D.matplot(smoothed_KL[,1:nb_topics],c(1,0),c(0,0),c(0,1),
                extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,
                ylab="Topics ranked by site-normalized abundances",
                xlab="Sequences ranked by\n site-normalized abundances",
                main="KL symmetrized divergence between smoothed topics' site distribution\n and smoothed sequences' site distribution -\n topics and sequences ranked by site-averaged relative abundances",
                do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
                border="black",na.color=NA)
dev.off()

pdf("KL2_distance_between_smoothed_topic_and_sequence_site_repartitions.pdf")
#par(mar=c(5.1,4.1,4.1,2.1))
par(mar=c(6.1,6.1,5.1,2.1))
color2D.matplot(smoothed_KL2[,1:nb_topics],c(1,0),c(0,0),c(0,1),
                extremes=NA,cellcolors=NA,show.legend=F,nslices=10,
                ylab="Classes ranked by\n relative abundances averaged over samples",xlab="",
                main="Symmetrized KL divergence\n between smoothed topics' spatial distribution\n and smoothed sequences' spatial distribution",
                #main="KL symmetrized distance between smoothed topics' site repartition\n and smoothed sequences' site repartition -\n topics and sequences ranked by site-normalized abundances",
                do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2,
                border="black",na.color=NA)
mtext("Sequences ranked by\n relative abundances averaged over samples",side=1,line=4,cex=1.5) 
dev.off()

pdf("KL2_distance_between_smoothed_topic_and_sequence_site_repartitions_seq=2nb_topics.pdf")
color2D.matplot(smoothed_KL2[,1:(2*nb_topics)],c(1,0),c(0,0),c(0,1),
                extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,
                ylab="Topics ranked by site-normalized abundances",
                xlab="Sequences ranked by\n site-normalized abundances",
                main="KL symmetrized divergence between smoothed topics' site distribution\n and smoothed sequences' site distribution -\n classes and sequences ranked by site-averaged relative abundances",
                #main="KL symmetrized distance between smoothed topics' site repartition\n and smoothed sequences' site repartition -\n topics and sequences ranked by site-normalized abundances",
                do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
                border="black",na.color=NA)
dev.off()

pdf("KL2_distance_between_smoothed_topic_and_sequence_site_repartitions_seq=3nb_topics.pdf")
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

pdf("KL_distance_between_topic_and_sequence_site_repartitions_15firstseq.pdf")
color2D.matplot(KL[1:15,1:15],c(1,0),c(0,0),c(0,1),
                extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,
                ylab="Topics ranked by site-normalized abundances",
                xlab="Sequences ranked by\n site-normalized abundances",
                main="KL symmetrized divergence between topics' site distribution\n and sequences' site distribution -\n topics and sequences ranked by site-averaged relative abundances",
                #main="KL symmetrized distance between topics' site repartition\n and sequences' site repartition -\n topics and sequences ranked by site-normalized abundances",
                do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
                border="black",na.color=NA)
dev.off()

pdf("KL_distance_(log)_between_topic_and_sequence_site_repartitions_20firstseq.pdf")
color2D.matplot(log(KL[1:20,1:20]),c(1,0),c(0,0),c(0,1),
                extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,
                ylab="Topics ranked by site-normalized abundances",
                xlab="Sequences ranked by\n site-normalized abundances",
                main="KL symmetrized divergence (log) between topics' site distribution\n and sequences' site distribution -\n topics and sequences ranked by site-averaged relative abundances",
                #main="KL symmetrized distance (log) between topics' site repartition\n and sequences' site repartition -\n topics and sequences ranked by site-normalized abundances",
                do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
                border="black",na.color=NA)
dev.off()

KL[which(KL[,1:nb_topics]==0)]=NA
pdf("KL_distance_(log)_between_topic_and_sequence_site_repartitions.pdf")
color2D.matplot(log(KL[,1:nb_topics]),c(1,0),c(0,0),c(0,1),
                extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,
                ylab="Topics ranked by site-normalized abundances",
                xlab="Sequences ranked by\n site-normalized abundances",
                main="KL symmetrized divergence (log) between topics' site distribution\n and sequences' site distribution -\n topics and sequences ranked by site-averaged relative abundances",
                #main="KL symmetrized distance (log) between topics' site repartition\n and sequences' site repartition -\n topics and sequences ranked by site-normalized abundances",
                do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
                border="black",na.color=NA)
dev.off()

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
  spatial_topicmix = matrix(nrow=19,ncol=19)
  for (j in 1:19)
  {
    for (i in 1:19)
    {
      spatial_topicmix[i,j] = norm_documents[(j-1)*19+i,rev(sort_prop_topic$ix)[k]]
    }
  }
  # colors Red Green Blue
  color2D.matplot(spatial_topicmix,c(0,1),c(0,0),c(1,0),
                  extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,xlab="Column",
                  ylab="Row",do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
                  border="black",na.color=NA,main=paste("Topic ",k," - ",rev(sort_prop_topic$ix)[k],sep=""))
}
dev.off()

pdf("First_topic_ordered_by_site-normalized_abundance_composition_map.pdf")
# 2 topics :
#par(mfrow=c(2,1))
# 5 topics :
#par(mfrow=c(3,2))
# 10 topics
#par(mfrow=c(2,2))
#lattice::levelplot(abund2.pred~x+y, z2, 
#col.regions=topo.colors(100), aspect = "iso",contour=T,main=binomial)
# Loop over topics (one map per topic, the color stands for the proportion of the topic)
spatial_topicmix = matrix(nrow=19,ncol=19)
for (j in 1:19)
  {
  for (i in 1:19)
    {
    spatial_topicmix[i,j] = norm_documents[(j-1)*19+i,rev(sort_normal_topic$ix)[1]]
    }
  }
  # colors Red Green Blue
  color2D.matplot(spatial_topicmix,c(0,1),c(0,0),c(1,0),
                  extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,xlab="Column",
                  ylab="Row",do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2,
                  border="black",na.color=NA,main="Spatial distribution of the most abundant class")
dev.off()

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
  spatial_topicmix = matrix(nrow=19,ncol=19)
  for (j in 1:19)
  {
    for (i in 1:19)
    {
      spatial_topicmix[i,j] = norm_documents[(j-1)*19+i,rev(sort_normal_topic$ix)[k]]
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
  spatial_topicmix = matrix(nrow=19,ncol=19)
  for (j in 1:19)
  {
    for (i in 1:19)
    {
      spatial_topicmix[i,j] = smoothed_KL_norm_documents[(j-1)*19+i,rev(sort_normal_topic$ix)[k]]
    }
  }
  # colors Red Green Blue
  color2D.matplot(spatial_topicmix,c(0,1),c(0,0),c(1,0),
                  extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,xlab="Column",
                  ylab="Row",do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
                  border="black",na.color=NA,main=paste("Topic ",rev(sort_normal_topic$ix)[k]," - #",k,sep=""))
}
dev.off()
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

