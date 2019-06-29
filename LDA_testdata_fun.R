
LDA_testdata_fun <- function(nb_rndzations_true_documents,samplewise,MOTUwise,
                             ncol_data2m,nrow_data2m,sum_data2m,nb_topics,
                             KL_documents,documents,true_documents,KL_topic_compo){

# true_KL_topic_compo is missing.  
  
# Computing spatial randomizations for the first realization so as to compare the LDA-retrieved topics with the true topics 
KL_documents_randomized = vector(length=nb_rndzations_true_documents,mode="list")
for (rndzation in 1:nb_rndzations_true_documents)
{
  if (samplewise)
  {
    KL_documents_randomized[[rndzation]] = matrix(nrow=ncol_data2m,ncol=nb_topics,data=0)
    # Performing "unidimensional spatial randomizations" (i.e., circular permutations)
    for (k in 1:nb_topics)
    {
      a = ncol_data2m
      while (a==ncol_data2m)
        a=sample(1:ncol_data2m,1)
      KL_documents_randomized[[rndzation]][,k] = KL_documents[c((a+1):ncol_data2m,1:a),k] 
    }
  } else if (MOTUwise)
  {
    KL_topic_compo_randomized = vector(length=nb_rndzations_true_documents,mode="list")
    # Performing permutations
    for (rndzation in 1:nb_rndzations_true_documents)
      KL_topic_compo_randomized[[rndzation]] = KL_topic_compo[sample(seq(1,nrow_data2m,1),nrow_data2m),]
  }
}

Correlation_true_documents = cor(documents,true_documents) 

# Computing a one-to-one topic correspondence between the topics in the first realization
# and the true documents
True_topic_correspondence = vector(length=nb_topics,mode="numeric")
i = 1
nb_topics_assigned = 0
k_vect = vector(length=nb_topics,mode="numeric")
k1_vect = vector(length=nb_topics,mode="numeric")
sorted_correlation_true_documents = sort.int(Correlation_true_documents,index.return=T)
while (nb_topics_assigned < nb_topics)
{
  k = sorted_correlation_true_documents$ix[i] %% nb_topics 
  k1 = sorted_correlation_true_documents$ix[i] %/% nb_topics
  if (k==0)
    k = nb_topics
  else 
    k1 = k1 + 1
  if (!any(k_vect == k) && !any(k1_vect == k1))
  {
    True_topic_correspondence[k] = k1
    nb_topics_assigned = nb_topics_assigned + 1
    k_vect[nb_topics_assigned] = k
    k1_vect[nb_topics_assigned] = k1
  }
  i = i+1
}

true_KL_documents = sweep(true_documents,MARGIN=2,colSums(true_documents),`/`)
true_KL_documents[true_KL_documents < 1/sum_data2m] = 1/sum_data2m

Mean_KL_topic_comparison_true_documents_randomized = vector(length=nb_topics,mode="numeric")
p_value_true_documents = vector(length = nb_topics, mode = "numeric")
KL_topic_comparison_true_documents = vector(length = nb_topics, mode="numeric")
DKL100_true_documents = vector(length = nb_topics, mode="numeric")
nES_true_documents = vector(length = nb_topics, mode="numeric")
for (k in 1:nb_topics)
{
  if (samplewise)
    KL_topic_comparison_true_documents[k] = 1/2*(KL.plugin(KL_documents[,k],true_KL_documents[,True_topic_correspondence[k]]) +
                                                   KL.plugin(true_KL_documents[,True_topic_correspondence[k]],KL_documents[,k]))
  else if (MOTUwise)
    KL_topic_comparison_true_documents[k] = 1/2*(KL.plugin(KL_topic_compo[,k],true_KL_topic_compo[,True_topic_correspondence[k]]) +
                                                   KL.plugin(true_KL_topic_compo[,True_topic_correspondence[k]],KL_topic_compo[,k]))
  
  nb_non_significant_rndzations = 0
  for (rndzation in 1:nb_rndzations_true_documents) 
  {
    if (samplewise)
      KL_topic_comparison_true_documents_randomized = 1/2*(KL.plugin(KL_documents_randomized[[rndzation]][,k],true_KL_documents[,True_topic_correspondence[k]]) +
                                                             KL.plugin(true_KL_documents[,True_topic_correspondence[k]],KL_documents_randomized[[rndzation]][,k]))
    else if (MOTUwise)
      KL_topic_comparison_true_documents_randomized = 1/2*(KL.plugin(KL_topic_compo_randomized[[rndzation]][,k],true_KL_topic_compo[,True_topic_correspondence[k]]) +
                                                             KL.plugin(true_KL_topic_compo[,True_topic_correspondence[k]],KL_topic_compo_randomized[[rndzation]][,k]))
    
    if (KL_topic_comparison_true_documents[k] > KL_topic_comparison_true_documents_randomized)
      nb_non_significant_rndzations = nb_non_significant_rndzations + 1
    Mean_KL_topic_comparison_true_documents_randomized[k] = 1/nb_rndzations_true_documents*KL_topic_comparison_true_documents_randomized + Mean_KL_topic_comparison_true_documents_randomized[k]
  }
  p_value_true_documents[k] = nb_non_significant_rndzations/nb_rndzations_true_documents
  DKL100_true_documents[k] = (Mean_KL_topic_comparison_true_documents_randomized[k] - KL_topic_comparison_true_documents[k])
  nES_true_documents[k] = DKL100_true_documents[k]/Mean_KL_topic_comparison_true_documents_randomized[k]
}

return(list(Correlation_true_documents = Correlation_true_documents,
            Mean_KL_topic_comparison_true_documents_randomized = Mean_KL_topic_comparison_true_documents_randomized,
            p_value_true_documents = p_value_true_documents,
            KL_topic_comparison_true_documents = KL_topic_comparison_true_documents,
            DKL100_true_documents = DKL100_true_documents,
            nES_true_documents = nES_true_documents))
}