# Variables values for testing:
# locally_based = 1
# genotoul_cluster_based = 0
# EDB_cluster_based = 0
# j_select = 2
# j_real = 1
# nb_topics = 3
# samplewise = 1
# MOTUwise = 0
# nb_rndzations = 10
# documents_allreal = list()
# for (i in 1:5) 
#   documents_allreal[[i]] = replicate(nb_topics,abs(rnorm(1131)))
# KL_documents_allreal = list()
# for (i in 1:5) 
#   KL_documents_allreal[[i]] = replicate(nb_topics,abs(rnorm(1131)))
# data2m = replicate(1131,abs(rnorm(100)))
# testdata = 0
# filled = 1
# filled_with_gaps =0
# bij = 1

LDA_realization_comparison_fun <- function(j_select,j_real,KL_documents_jselect_randomized,KL_topic_compo_jselect_randomized,nb_topics,samplewise,MOTUwise,
                                           documents_allreal,KL_documents_allreal,topic_compo_allreal,KL_topic_compo_allreal,
                                           testdata,filled,filled_with_gaps,nb_rndzations,ncol_data2m,nrow_data2m,Missing_positions_indices,bij){
  
  #KL_MOTUwise_topic_comparison = matrix(nrow=nb_topics,ncol=nb_topics,data=0)
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
  
  # Randomizations:
  if (j_real == 1) 
  {
    if (!testdata)
    {
      if (samplewise)
      {
        if (filled && !filled_with_gaps)
        {
          #KL_documents_jselect_randomized = vector(length=nb_rndzations,mode="list")
          for (rndzation in 1:nb_rndzations)
          {
            KL_documents_jselect_randomized[[rndzation]] = matrix(nrow=ncol_data2m,ncol=nb_topics,data=0)
            # Performing spatial randomizations
            for (k in 1:nb_topics)
            {
              Spatial_KL_documents = matrix(KL_documents_allreal[[j_select]][,k],ncol=39)
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
              KL_documents_jselect_randomized[[rndzation]][,k] = as.vector(Spatial_KL_documents_randomized)
            }
          }
        } else if (!filled || filled_with_gaps)
        {
          Spatial_KL_documents = vector(length=nb_topics,mode="list")
          for (k in 1:nb_topics)
          {
            Spatial_KL_documents[[k]] = matrix(nrow=29,ncol=39,data=0)
            position_shift = 0
            for (j in 1:39)
            {
              for (i in 1:29)
              {
                if (Missing_positions_indices[(j-1)*29+i]==0)
                  Spatial_KL_documents[[k]][i,j] = KL_documents_allreal[[j_select]][(j-1)*29+i-position_shift,k]    
                else if (Missing_positions_indices[(j-1)*29+i]==1)
                {
                  Spatial_KL_documents[[k]][i,j] = NA
                  position_shift = position_shift+1
                }
              }
            }
          }
          #KL_documents_jselect_randomized = vector(length=nb_rndzations,mode="list")
          for (rndzation in 1:nb_rndzations)
          {
            KL_documents_jselect_randomized[[rndzation]] = matrix(nrow=ncol(data2m),ncol=nb_topics,data=0)
            # Performing spatial randomizations
            for (k in 1:nb_topics)
            {
              a=29
              b=39
              while ((a==29) && (b==39))
              {
                a=sample(1:29,1)
                b=sample(1:39,1)
              }  
              if ((a!=29) && (b!=39))
                Spatial_KL_documents_randomized = Spatial_KL_documents[[k]][c((a+1):29,1:a),c((b+1):39,1:b)] 
              else if (b==39)
                Spatial_KL_documents_randomized = Spatial_KL_documents[[k]][c((a+1):29,1:a),]
              else if (a==29)
                Spatial_KL_documents_randomized = Spatial_KL_documents[[k]][,c((b+1):39,1:b)] 
              KL_documents_jselect_randomized[[rndzation]][,k] = as.vector(Spatial_KL_documents_randomized[!is.na(Spatial_KL_documents_randomized)])
            }
          }
        }
      } else if (MOTUwise)
      {
        #KL_topic_compo_jselect_randomized = vector(length=nb_rndzations,mode="list")
        # Performing permutations (non-spatial randomizations)
        for (rndzation in 1:nb_rndzations)
          KL_topic_compo_jselect_randomized[[rndzation]] = KL_topic_compo_allreal[[j_select]][sample(seq(1,nrow_data2m,1),nrow_data2m),]
      }
    } else if (testdata)
    {
      if (samplewise)
      {
        #KL_documents_jselect_randomized = vector(length=nb_rndzations,mode="list")
        for (rndzation in 1:nb_rndzations)
        {
          KL_documents_jselect_randomized[[rndzation]] = matrix(nrow=ncol_data2m,ncol=nb_topics,data=0)
          # Performing spatial randomizations
          for (k in 1:nb_topics)
          {
            a=1131
            while (a==1131)
              a=sample(1:1131,1)
            KL_documents_jselect_randomized[[rndzation]][,k] = KL_documents_allreal[[j_select]][c((a+1):1131,1:a),k]
          }
        }
      } else if (MOTUwise)
      {
        #KL_topic_compo_jselect_randomized = vector(length=nb_rndzations,mode="list")
        
        # Performing permutations (non-spatial randomizations)
        for (rndzation in 1:nb_rndzations)
          KL_topic_compo_jselect_randomized[[rndzation]] = KL_topic_compo_allreal[[j_select]][sample(seq(1,nrow_data2m,1),nrow_data2m),]
      }
    } 
  }
  
  #Best_topic_comparison_MOTUwise = apply(KL_MOTUwise_topic_comparison,1,min)
  #Best_topic_comparison_samplewise = apply(KL_samplewise_topic_comparison,1,min)
  #KL_documents_allreal[[j_real]]
  nb_non_significant_rndzations = 0
  Var_KL_topic_comparison_randomized = vector(length=nb_topics,mode="numeric")
  Mean_KL_topic_comparison_randomized = vector(length=nb_topics,mode="numeric")
  p_value = vector(length=nb_topics,mode="numeric")
  
  if (bij)
  {
    # Computing a oneto-one topic correspondence between j_real and j_select (no topic of j_real can be associated with two different topics in j_select)
    #Topic_correspondence_MOTUwise = vector(length=nb_topics,mode="numeric")
    Topic_correspondence = vector(length=nb_topics,mode="numeric")
    i = 1
    nb_topics_assigned = 0
    k_vect = vector(length=nb_topics,mode="numeric")
    k1_vect = vector(length=nb_topics,mode="numeric")
    sorted_KL_topic_comparison = sort.int(KL_topic_comparison,index.return=T)
    while (nb_topics_assigned < nb_topics)
    {
      k = sorted_KL_topic_comparison$ix[i] %% nb_topics 
      k1 = sorted_KL_topic_comparison$ix[i] %/% nb_topics
      if (k==0)
        k = nb_topics
      else 
        k1 = k1 + 1
      if (!any(k_vect == k) && !any(k1_vect == k1))
      {
        Topic_correspondence[k] = k1
        nb_topics_assigned = nb_topics_assigned+1
        k_vect[nb_topics_assigned] = k
        k1_vect[nb_topics_assigned] = k1
      }
      i = i+1
    }
  } else if (!bij)
  {
    Topic_correspondence = matrix(nrow=nb_topics,ncol=2,data=NA)
    sorted_KL_topic_comparison = sort.int(KL_topic_comparison,index.return=T)
    for (i in 1:nb_topics)
    {
      k = sorted_KL_topic_comparison$ix[i] %% nb_topics 
      k1 = sorted_KL_topic_comparison$ix[i] %/% nb_topics
      if (k==0)
        k = nb_topics
      else 
        k1 = k1 + 1
      Topic_correspondence[i,1] = k
      Topic_correspondence[i,2] = k1
    }
  }
  
  if (bij)
  {
    for (k in 1:nb_topics)
    {          
      for (rndzation in 1:nb_rndzations) 
      {
        if (samplewise)
          KL_topic_comparison_randomized = 1/2*(KL.plugin(KL_documents_allreal[[j_real]][,k],KL_documents_jselect_randomized[[rndzation]][,Topic_correspondence[k]]) +
                                                  KL.plugin(KL_documents_jselect_randomized[[rndzation]][,Topic_correspondence[k]],KL_documents_allreal[[j_real]][,k]))
        else if (MOTUwise)
          KL_topic_comparison_randomized = 1/2*(KL.plugin(KL_topic_compo_allreal[[j_real]][,k],KL_topic_compo_jselect_randomized[[rndzation]][,Topic_correspondence[k]]) +
                                                  KL.plugin(KL_topic_compo_jselect_randomized[[rndzation]][,Topic_correspondence[k]],KL_topic_compo_allreal[[j_real]][,k]))
        if (KL_topic_comparison[k,Topic_correspondence[k]] > KL_topic_comparison_randomized)
          nb_non_significant_rndzations = nb_non_significant_rndzations + 1
        Mean_KL_topic_comparison_randomized[k] = 1/nb_rndzations*KL_topic_comparison_randomized +
          Mean_KL_topic_comparison_randomized[k]
        Var_KL_topic_comparison_randomized[k] = KL_topic_comparison_randomized^2/nb_rndzations + 
          Var_KL_topic_comparison_randomized[k]
      }
      p_value[k] = nb_non_significant_rndzations/nb_rndzations
    }
  } else if (!bij)
  {
    for (k in 1:nb_topics)
    {
      k0 = Topic_correspondence[i,1]
      k1 = Topic_correspondence[i,2]
      for (rndzation in 1:nb_rndzations) 
      {
        KL_topic_comparison_randomized = 1/2*(KL.plugin(KL_documents_allreal[[j_real]][,k0],KL_documents_jselect_randomized[[rndzation]][,k1]) +
                                                KL.plugin(KL_documents_jselect_randomized[[rndzation]][,k1],KL_documents_allreal[[j_real]][,k0]))
        if (KL_topic_comparison[k0,k1] > KL_topic_comparison_randomized)
          nb_non_significant_rndzations = nb_non_significant_rndzations + 1
        Mean_KL_topic_comparison_randomized[k] = 1/nb_rndzations*KL_topic_comparison_randomized +
          Mean_KL_topic_comparison_randomized[k]
        Var_KL_topic_comparison_randomized[k] = KL_topic_comparison_randomized^2/nb_rndzations + 
          Var_KL_topic_comparison_randomized[k]
      }
      p_value[k] = nb_non_significant_rndzations/nb_rndzations
    }
  }
  
  #if (length(which(Best_topic_comparison_MOTUwise != Best_topic_comparison_samplewise)) > 0)
  #cat("Unconsistent topic correspondence between real.",j_select,"and",j_real,"\n")        
  
  return(list(Topic_correspondence = Topic_correspondence,
              KL_topic_comparison = KL_topic_comparison,
              Mean_KL_topic_comparison_randomized = Mean_KL_topic_comparison_randomized,
              p_value = p_value,
              Var_KL_topic_comparison_randomized = Var_KL_topic_comparison_randomized,
              KL_documents_jselect_randomized = KL_documents_jselect_randomized,
              KL_topic_compo_jselect_randomized = KL_topic_compo_jselect_randomized))
  
  # end function  
}


