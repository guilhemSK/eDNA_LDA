# Variables values for testing:
# j_select = 2
# j_real = 1
# nb_topics = 3
# samplewise = 1
# MOTUwise = 0
# nb_rndzations = 10
# documents_allreal = list()
# for (i in 1:5) 
#   documents_allreal[[i]] = replicate(nb_topics,abs(rnorm(nrow_grid*ncol_grid)))
# KL_documents_allreal = list()
# for (i in 1:5) 
#   KL_documents_allreal[[i]] = replicate(nb_topics,abs(rnorm(nrow_grid*ncol_grid)))
# data2m = replicate(nrow_grid*ncol_grid,abs(rnorm(100)))
# testdata = 0
# missing = 0

LDA_realization_comparison_fun <- function(j_select,j_real,KL_documents_jselect_randomized,KL_topic_compo_jselect_randomized,nb_topics,samplewise,MOTUwise,
                                           KL_documents_allreal,KL_topic_compo_allreal,nrow_grid,ncol_grid,byRow,
                                           testdata,missing,nb_rndzations,ncol_data2m,nrow_data2m,Missing_positions_indices,
                                           KL_topic_comparison,Topic_correspondence){
  
  # Computing randomizations:
  if (j_real == 1) 
  {
    if (testdata)
    {
      if (samplewise)
      {
        for (rndzation in 1:nb_rndzations)
        {
          KL_documents_jselect_randomized[[rndzation]] = matrix(nrow=ncol_data2m,ncol=nb_topics,data=0)
          # Performing spatial randomizations over the spatial distribution of topics
          for (k in 1:nb_topics)
          {
            a=ncol_data2m
            while (a==ncol_data2m)
              a=sample(1:ncol_data2m,1)
            KL_documents_jselect_randomized[[rndzation]][,k] = KL_documents_allreal[[j_select]][c((a+1):ncol_data2m,1:a),k]
          }
        }
      } else if (MOTUwise)
      {
        # Performing randomizations over the taxonomic composition of topics
        for (rndzation in 1:nb_rndzations)
          KL_topic_compo_jselect_randomized[[rndzation]] = KL_topic_compo_allreal[[j_select]][sample(seq(1,nrow_data2m,1),nrow_data2m),]
      }
    } else 
    {
      if (samplewise)
      {
        if (grid && !geographic)
        {  
          Spatial_KL_documents = vector(length=nb_topics,mode="list")
          for (k in 1:nb_topics)
          {
            if (!missing) 
            {
              if (byRow)
                Spatial_KL_documents[[k]] = t(matrix(KL_documents_allreal[[j_select]][,k],ncol=nrow_grid))
              else
                Spatial_KL_documents[[k]] = matrix(KL_documents_allreal[[j_select]][,k],ncol=ncol_grid)
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
                    missing_index = (j-1)*nrow_grid+i
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
          for (rndzation in 1:nb_rndzations)
          {
            KL_documents_jselect_randomized[[rndzation]] = matrix(nrow=ncol_data2m,ncol=nb_topics,data=0)
            # Performing spatial randomizations over the spatial distribution of topics
            for (k in 1:nb_topics)
            {
              a = nrow_grid
              b = ncol_grid
              while (a == nrow_grid && b == ncol_grid)
              {
                a = sample(1:nrow_grid,1)
                b = sample(1:ncol_grid,1)
              }  
              if (a != nrow_grid && b != ncol_grid)
                Spatial_KL_documents_randomized = Spatial_KL_documents[[k]][c((a+1):nrow_grid,1:a),c((b+1):ncol_grid,1:b)] 
              else if (b==ncol_grid)
                Spatial_KL_documents_randomized = Spatial_KL_documents[[k]][c((a+1):nrow_grid,1:a),]
              else if (a==nrow_grid)
                Spatial_KL_documents_randomized = Spatial_KL_documents[[k]][,c((b+1):ncol_grid,1:b)]
              
              if (byRow)
                KL_documents_jselect_randomized[[rndzation]][,k] = as.vector(t(Spatial_KL_documents_randomized))[!is.na(as.vector(t(Spatial_KL_documents_randomized)))]
              else  
                KL_documents_jselect_randomized[[rndzation]][,k] = Spatial_KL_documents_randomized[!is.na(Spatial_KL_documents_randomized)]
            }
          }
        } else
        {
          # Non-spatial permutations:
          for (rndzation in 1:nb_rndzations)
            KL_documents_jselect_randomized[[rndzation]] = KL_documents_allreal[[j_select]][sample(seq(1,ncol_data2m,1),ncol_data2m),]
        }
      } else if (MOTUwise)
      {
        # Performing randomizations over the taxonomic composition of topics
        for (rndzation in 1:nb_rndzations)
          KL_topic_compo_jselect_randomized[[rndzation]] = KL_topic_compo_allreal[[j_select]][sample(seq(1,nrow_data2m,1),nrow_data2m),]
      }
    }
  }
  
  nb_non_significant_rndzations = 0
  Mean_KL_topic_comparison_randomized = vector(length=nb_topics,mode="numeric")
  p_value = vector(length=nb_topics,mode="numeric")
  
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
    }
    p_value[k] = nb_non_significant_rndzations/nb_rndzations
  }
  
  return(list(Mean_KL_topic_comparison_randomized = Mean_KL_topic_comparison_randomized,
              p_value = p_value,
              KL_documents_jselect_randomized = KL_documents_jselect_randomized,
              KL_topic_compo_jselect_randomized = KL_topic_compo_jselect_randomized))
}


