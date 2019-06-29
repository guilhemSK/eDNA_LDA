

LDA_best_real_fun <- function(nb_rndzations_best_real,nrow_data2m,KL_documents,KL_topic_compo,assemblage_names_vect,documents,topic_compo,sum_data2m){
  
KL_topic_compo_randomized = vector(length=nb_rndzations_best_real,mode="list")
# Performing permutations
for (rndzation in 1:nb_rndzations_best_real)
{
  KL_topic_compo_randomized[[rndzation]] = matrix(nrow=nrow_data2m,ncol=nb_topics,data=0)
  for (k in 1:nb_topics)
    KL_topic_compo_randomized[[rndzation]][,k] = KL_topic_compo[sample(seq(1,nrow_data2m,1),nrow_data2m),k]
}

KL_topic_comparison_samplewise = matrix(nrow=nb_topics,ncol=nb_topics,data=NA)
KL_topic_comparison_MOTUwise = matrix(nrow=nb_topics,ncol=nb_topics,data=NA)
nES_topic_comparison_MOTUwise = matrix(nrow=nb_topics,ncol=nb_topics,data=NA)
progress_index = 0
for (i_topic in 1:(nb_topics-1))
{
  for (j_topic in (i_topic+1):nb_topics)
  {
    cat("Assemblage comparison within the best realization:",format(progress_index/(nb_topics*(nb_topics-1)/2)*100,digits=2),"%\n")
    progress_index = progress_index+1
    
    #k = sort_normal_topic$ix[i_topic]
    #k1 = sort_normal_topic$ix[j_topic]
    k = i_topic
    k1 = j_topic
    KL_topic_comparison_samplewise[i_topic,j_topic] = 1/2*(KL.plugin(KL_documents[,k],KL_documents[,k1]) + KL.plugin(KL_documents[,k1],KL_documents[,k]))
    KL_topic_comparison_samplewise[j_topic,i_topic] = KL_topic_comparison_samplewise[i_topic,j_topic]
    KL_topic_comparison_MOTUwise[i_topic,j_topic] = 1/2*(KL.plugin(KL_topic_compo[,k],KL_topic_compo[,k1]) + KL.plugin(KL_topic_compo[,k1],KL_topic_compo[,k]))
    KL_topic_comparison_MOTUwise[j_topic,i_topic] = KL_topic_comparison_MOTUwise[i_topic,j_topic]
    
    KL_topic_comparison_MOTUwise_randomized = vector(length=nb_rndzations_best_real,mode="numeric")
    for (rndzation in 1:nb_rndzations_best_real)
      KL_topic_comparison_MOTUwise_randomized[rndzation] = 1/2*(KL.plugin(KL_topic_compo[,k],KL_topic_compo_randomized[[rndzation]][,k1]) + KL.plugin(KL_topic_compo_randomized[[rndzation]][,k1],KL_topic_compo[,k]))
    
    nES_topic_comparison_MOTUwise[i_topic,j_topic] = (mean(KL_topic_comparison_MOTUwise_randomized) - KL_topic_comparison_MOTUwise[i_topic,j_topic])/mean(KL_topic_comparison_MOTUwise_randomized)
    nES_topic_comparison_MOTUwise[j_topic,i_topic] = nES_topic_comparison_MOTUwise[i_topic,j_topic]
  }
}

colnames(KL_topic_comparison_samplewise) = assemblage_names_vect
colnames(KL_topic_comparison_MOTUwise) = assemblage_names_vect
colnames(nES_topic_comparison_MOTUwise) = assemblage_names_vect
rownames(KL_topic_comparison_samplewise) = assemblage_names_vect
rownames(KL_topic_comparison_MOTUwise) = assemblage_names_vect
rownames(nES_topic_comparison_MOTUwise) = assemblage_names_vect

#documents1 = documents[,sort_normal_topic$ix]
#topic_compo1 = topic_compo[,sort_normal_topic$ix]
documents1 = documents
topic_compo1 = topic_compo
colnames(documents1) = assemblage_names_vect
colnames(topic_compo1) = assemblage_names_vect
Corr_topic_comparison_samplewise = cor(documents1)
Corr_topic_comparison_MOTUwise = cor(topic_compo1)

Hellinger_topic_comparison_MOTUwise = 1/sqrt(2)*dist(t(sqrt(topic_compo1)), method = "euclidean", diag = FALSE, upper = FALSE)
#KL_topic_compo1 = KL_topic_compo[,sort_normal_topic$ix]
KL_topic_compo1 = KL_topic_compo
topic_compo1_bin = topic_compo1
topic_compo1_bin[topic_compo1 < 1/sum_data2m] = 0
Jaccard_topic_comparison_MOTUwise = vegan::designdist(t(topic_compo1_bin), "(b+c)/(a+b+c)", abcd=TRUE)
Beta.sim_topic_comparison_MOTUwise = vegan::designdist(t(topic_compo1_bin), "pmin(b,c)/(pmin(b,c)+a)", abcd=TRUE)

return(list(Hellinger_topic_comparison_MOTUwise = Hellinger_topic_comparison_MOTUwise,
       Jaccard_topic_comparison_MOTUwise = Jaccard_topic_comparison_MOTUwise,
       Beta.sim_topic_comparison_MOTUwise = Beta.sim_topic_comparison_MOTUwise,
       Corr_topic_comparison_samplewise = Corr_topic_comparison_samplewise,
       Corr_topic_comparison_MOTUwise = Corr_topic_comparison_MOTUwise,
       KL_topic_comparison_samplewise = KL_topic_comparison_samplewise,
       KL_topic_comparison_MOTUwise = KL_topic_comparison_MOTUwise,
       nES_topic_comparison_MOTUwise = nES_topic_comparison_MOTUwise))
}