
LDA_testdata_fun <- function(nb_rndzations_true_documents,samplewise,MOTUwise,ncol_data2m,nrow_data2m,nb_topics,KL_documents,documents,true_documents,KL_topic_compo){

# !!! true_KL_topic_compo is missing !!!  
  
# Computing randomizations in a spacialized way for the first realization for testdata=1 so as to compare the LDA results with the true topics 
KL_documents_randomized = vector(length=nb_rndzations_true_documents,mode="list")
for (rndzation in 1:nb_rndzations_true_documents)
{
  if (samplewise)
  {
    KL_documents_randomized[[rndzation]] = matrix(nrow=ncol_data2m,ncol=nb_topics,data=0)
    # Performing spatial randomizations
    for (k in 1:nb_topics)
    {
      Spatial_KL_documents = matrix(KL_documents[,k],ncol=39)
      a=29
      b=39
      while ((a==29) && (b==39))
      {
        a=sample(1:29,1)
        b=sample(1:39,1)
      }  
      if ((a!=29) && (b!=39))
        Spatial_KL_documents_randomized = Spatial_KL_documents[c((a+1):29,1:a),c((b+1):39,1:b)] 
      else if (b==39)
        Spatial_KL_documents_randomized = Spatial_KL_documents[c((a+1):29,1:a),]
      else if (a==29)
        Spatial_KL_documents_randomized = Spatial_KL_documents[,c((b+1):39,1:b)] 
      KL_documents_randomized[[rndzation]][,k] = as.vector(Spatial_KL_documents_randomized)
    }
  } else if (MOTUwise)
  {
    KL_topic_compo_randomized = vector(length=nb_rndzations_true_documents,mode="list")
    # Performing permutations (non-spatial randomizations)
    for (rndzation in 1:nb_rndzations_true_documents)
      KL_topic_compo_randomized[[rndzation]] = KL_topic_compo[sample(seq(1,nrow_data2m,1),nrow_data2m),]
  }
}

#         # Computing randomizations for the first realization for testdata=1 so as to compare the LDA results with the true topics 
#         KL_documents_jselect_randomized = vector(length=nb_rndzations_true_documents,mode="list")
#         for (rndzation in 1:nb_rndzations_true_documents)
#         {
#           KL_documents_jselect_randomized[[rndzation]] = matrix(nrow=ncol(data2m),ncol=nb_topics,data=0)
#           # Performing spatial randomizations
#           for (k in 1:nb_topics)
#           {
#             a=1131
#             while (a==1131)
#               a=sample(1:1131,1)
#             KL_documents_jselect_randomized[[rndzation]][,k] = KL_documents_allreal[[j_select]][c((a+1):1131,1:a),k]
#           }
#         }

Correlation_true_documents = cor(documents,true_documents) 

# Computing a oneto-one topic correspondence between the first realization
# and the true documents
#Topic_correspondence_MOTUwise = vector(length=nb_topics,mode="numeric")
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

#Best_true_topic_comparison_samplewise = apply(Correlation_true_documents,1,max)

true_KL_documents = sweep(true_documents,MARGIN=2,colSums(true_documents),`/`)
true_KL_documents[true_KL_documents < 1/sum(data2m)] = 1/sum(data2m)

#True_topic_correspondence_samplewise = vector(length=nb_topics,mode="numeric")
Mean_KL_topic_comparison_true_documents_randomized = vector(length=nb_topics,mode="numeric")
p_value_true_documents = vector(length = nb_topics, mode = "numeric")
KL_topic_comparison_true_documents = vector(length = nb_topics, mode="numeric")
DKL100_true_documents = vector(length = nb_topics, mode="numeric")
nES_true_documents = vector(length = nb_topics, mode="numeric")
for (k in 1:nb_topics)
{
  #True_topic_correspondence_samplewise[k] = which(Best_true_topic_comparison_samplewise[k] == Correlation_true_documents[k,])
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

# end function
}

