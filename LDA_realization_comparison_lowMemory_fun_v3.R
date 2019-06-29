# Variables values for testing:
# j_select = 2
# j_real = 1
# nb_topics = 3
# samplewise = 1
# MOTUwise = 0
# nb_rndzations = 10
# KL_documents_allreal = list()
# for (i in 1:5) 
#   KL_documents_allreal[[i]] = replicate(nb_topics,abs(rnorm(ncol_grid*nrow_grid)))
# data2m = replicate(ncol_grid*nrow_grid,abs(rnorm(100)))
# testdata = 0
# missing = 0

LDA_realization_comparison_lowMemory_fun <- function(j_select,j_real,nb_topics,samplewise,MOTUwise,
                                                     KL_documents_allreal,KL_topic_compo_allreal,byRow,
                                                     testdata,missing,nb_rndzations,nrow_data2m,ncol_data2m,Missing_positions_indices,
                                                     KL_topic_comparison,Topic_correspondence){
  
  Mean_KL_topic_comparison_randomized = vector(length=nb_topics,mode="numeric")
  p_value = vector(length=nb_topics,mode="numeric")
  nb_non_significant_rndzations = vector(length=nb_topics,mode="numeric")
  
  if (grid && !geographic && samplewise)
    Spatial_KL_documents = vector(length=nb_topics,mode="list")
  for (rndzation in 1:nb_rndzations) 
  {
    for (k in 1:nb_topics)
    {
      k1 = Topic_correspondence[k]
      if (samplewise)
      {
        if (grid && !geographic)
        {
          if (rndzation == 1)
          {
            if (!missing) 
            {
              if (byRow)
                Spatial_KL_documents[[k]] = t(matrix(KL_documents_allreal[[j_select]][,k1],ncol=nrow_grid))
              else
                Spatial_KL_documents[[k]] = matrix(KL_documents_allreal[[j_select]][,k1],ncol=ncol_grid)
            } else if (missing)
            {
              Spatial_KL_documents[[k]] = matrix(nrow=nrow_grid, ncol=ncol_grid, data=0)
              position_shift = 0
              if (byRow)
              {
                for (i in 1:nrow_grid)
                {
                  for (j in 1:ncol_grid)
                  {
                    missing_index = (i-1)*ncol_grid+j
                    if (Missing_positions_indices[missing_index]==0)
                    {
                      Spatial_KL_documents[[k]][i,j] = KL_documents_allreal[[j_select]][missing_index-position_shift,k]
                    } else if (Missing_positions_indices[missing_index]==1)
                    {
                      Spatial_KL_documents[[k]][i,j] = NA
                      position_shift = position_shift+1
                    }
                  }
                }
              } else
              {
                for (j in 1:ncol_grid)
                {
                  for (i in 1:nrow_grid)
                  {
                    missing_index = missing_index = (j-1)*nrow_grid+i
                    if (Missing_positions_indices[missing_index]==0)
                    {
                      Spatial_KL_documents[[k]][i,j] = KL_documents_allreal[[j_select]][missing_index-position_shift,k]
                    } else if (Missing_positions_indices[missing_index]==1)
                    {
                      Spatial_KL_documents[[k]][i,j] = NA
                      position_shift = position_shift+1
                    }
                  }
                }
              }
            }
          }
          # Performing spatial randomizations over the spatial distribution of topics
          a = nrow_grid
          b = ncol_grid
          while (a == nrow_grid && b == ncol_grid)
          {
            a = sample(1:nrow_grid,1)
            b = sample(1:ncol_grid,1)
          }  
          if (a != nrow_grid && b != ncol_grid)
            Spatial_KL_documents_randomized = Spatial_KL_documents[[k]][c((a+1):nrow_grid,1:a),c((b+1):ncol_grid,1:b)] 
          else if (b == ncol_grid)
            Spatial_KL_documents_randomized = Spatial_KL_documents[[k]][c((a+1):nrow_grid,1:a),]
          else if (a == nrow_grid)
            Spatial_KL_documents_randomized = Spatial_KL_documents[[k]][,c((b+1):ncol_grid,1:b)] 
          
          if (byRow)
            KL_documents_jselect_randomized = as.vector(t(Spatial_KL_documents_randomized))[!is.na(as.vector(t(Spatial_KL_documents_randomized)))]
          else  
            KL_documents_jselect_randomized = Spatial_KL_documents_randomized[!is.na(Spatial_KL_documents_randomized)]
        } else if (testdata)
        {
          # Performing spatial randomizations over the spatial distribution of topics
          a = ncol_data2m
          while (a == ncol_data2m)
            a = sample(1:ncol_data2m,1)
          KL_documents_jselect_randomized = KL_documents_allreal[[j_select]][c((a+1):ncol_data2m,1:a),k1]
        } else
        {
          # Non-spatial permutations:
          KL_documents_jselect_randomized = KL_documents_allreal[[j_select]][sample(seq(1,ncol_data2m,1),ncol_data2m),k1]
        }
        
        KL_topic_comparison_randomized = 1/2*(KL.plugin(KL_documents_allreal[[j_real]][,k],KL_documents_jselect_randomized) +
                                                KL.plugin(KL_documents_jselect_randomized,KL_documents_allreal[[j_real]][,k]))
        if (KL_topic_comparison[k,k1] > KL_topic_comparison_randomized)
          nb_non_significant_rndzations[k] = nb_non_significant_rndzations[k] + 1
      } else if (MOTUwise)
      {
        # Performing randomizations over the taxonomic composition of topics
        KL_topic_compo_jselect_randomized = KL_topic_compo_allreal[[j_select]][sample(seq(1,nrow_data2m,1),nrow_data2m),k1]
        
        KL_topic_comparison_randomized = 1/2*(KL.plugin(KL_topic_compo_allreal[[j_real]][,k],KL_topic_compo_jselect_randomized) +
                                                KL.plugin(KL_topic_compo_jselect_randomized,KL_topic_compo_allreal[[j_real]][,k]))
        if (KL_topic_comparison[k,k1] > KL_topic_comparison_randomized)
          nb_non_significant_rndzations[k] = nb_non_significant_rndzations[k] + 1
      }
      
      Mean_KL_topic_comparison_randomized[k] = 1/nb_rndzations*KL_topic_comparison_randomized +
        Mean_KL_topic_comparison_randomized[k]
    }
  }
  p_value = nb_non_significant_rndzations/nb_rndzations   
  
  return(list(Mean_KL_topic_comparison_randomized = Mean_KL_topic_comparison_randomized,
              p_value = p_value))
}


