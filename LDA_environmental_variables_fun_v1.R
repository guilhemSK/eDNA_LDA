# Variables values for testing :
# grid = 1
# # Random matrix with 1131 columns and 100 rows:
# data2m = replicate(1131,abs(rnorm(100)))
# pca_abiotic = 0
# nb_topics = 3
# nb_abiotic_rndzations = 5
# documents = replicate(nb_topics,abs(rnorm(1131)))
# KL_documents = documents

LDA_environmental_variables_fun <- function(grid,missing,Missing_positions_indices,byRow,
                                        ncol_data2m,pca_abiotic,nb_topics,nb_abiotic_rndzations,
                                        documents,data_abiotic){
  
  ncol0 = ncol(data_abiotic)
  
  if (missing)
    data_abiotic = data_abiotic[-which(Missing_positions_indices==1),]
  
  if (pca_abiotic)
  {
    # PCA on standardized abiotic data, ie centered and divided by std. dev. for each column 
    abiotic_PCA = dudi.pca(as.data.frame(data_abiotic), row.w = rep(1, nrow(data_abiotic))/nrow(data_abiotic), 
                           col.w = rep(1, ncol(data_abiotic)), center = TRUE, scale = TRUE, 
                           scannf = F, nf = ncol(data_abiotic)) 
  }
  
  Correlation_abiotic = list()
  Correlation_abiotic[[1]] = matrix(nrow = nb_topics, ncol = ncol0)
  Correlation_abiotic[[2]] = matrix(nrow = nb_topics, ncol = ncol0)
  if (pca_abiotic)
  {
    Correlation_abiotic[[3]] = matrix(nrow = nb_topics, ncol = ncol0)
    Correlation_abiotic[[4]] = matrix(nrow = nb_topics, ncol = ncol0)
  }
  
  Mean_cor_abiotic_comparison_randomized = matrix(nrow=nb_topics,ncol=ncol0,data=0)
  Sd_cor_abiotic_comparison_randomized = matrix(nrow=nb_topics,ncol=ncol0,data=0)
  p_value_abiotic = matrix(nrow=nb_topics,ncol=ncol0,data=0)
  for (j_abiotic in 1:ncol0)
  {
    if ((j_abiotic == 1) && grid && !geographic)
    {
      if (!missing)
      { 
        documents_randomized = vector(length=nb_abiotic_rndzations,mode="list")
        KL_documents_randomized = vector(length=nb_abiotic_rndzations,mode="list")
        for (rndzation in 1:nb_abiotic_rndzations)
        {
          documents_randomized[[rndzation]] = matrix(nrow=ncol_data2m,ncol=nb_topics,data=0)
          KL_documents_randomized[[rndzation]] = matrix(nrow=ncol_data2m,ncol=nb_topics,data=0)
          # Performing spatial randomizations
          for (k in 1:nb_topics)
          {
            if (byRow)
              Spatial_documents = t(matrix(documents[,k],ncol=nrow_grid))
            else
              Spatial_documents = matrix(documents[,k],ncol=ncol_grid)
            
            a=nrow_grid
            b=ncol_grid
            while ((a==nrow_grid) && (b==ncol_grid))
            {
              a=sample(1:nrow_grid,1)
              b=sample(1:ncol_grid,1)
            }  
            if ((a!=nrow_grid) && (b!=ncol_grid))
              Spatial_documents_randomized = Spatial_documents[c((a+1):nrow_grid,1:a),c((b+1):ncol_grid,1:b)] 
            else if (b==ncol_grid)
              Spatial_documents_randomized = Spatial_documents[c((a+1):nrow_grid,1:a),]
            else if (a==nrow_grid)
              Spatial_documents_randomized = Spatial_documents[,c((b+1):ncol_grid,1:b)] 
            
            if (byRow)
              documents_randomized[[rndzation]][,k] = as.vector(t(Spatial_documents_randomized))
            else
              documents_randomized[[rndzation]][,k] = as.vector(Spatial_documents_randomized)
            KL_documents_randomized[[rndzation]][,k] = documents_randomized[[rndzation]][,k]/sum(documents_randomized[[rndzation]][,k])
          }
        }
      } else if (missing)
      { 
        Spatial_documents = vector(length=nb_topics,mode="list")
        for (k in 1:nb_topics)
        {
          Spatial_documents[[k]] = matrix(nrow=nrow_grid, ncol=ncol_grid, data=0)
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
                  Spatial_documents[[k]][i,j] = documents[missing_index-position_shift,k]
                } else if (Missing_positions_indices[missing_index]==1)
                {
                  Spatial_documents[[k]][i,j] = NA
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
                  Spatial_documents[[k]][i,j] = documents[missing_index-position_shift,k]
                } else if (Missing_positions_indices[missing_index]==1)
                {
                  Spatial_documents[[k]][i,j] = NA
                  position_shift = position_shift+1
                }
              }
            }
          }
        }
        documents_randomized = vector(length=nb_abiotic_rndzations,mode="list")
        KL_documents_randomized = vector(length=nb_abiotic_rndzations,mode="list")
        for (rndzation in 1:nb_abiotic_rndzations)
        {
          documents_randomized[[rndzation]] = matrix(nrow=ncol(data2m),ncol=nb_topics,data=0)
          KL_documents_randomized[[rndzation]] = matrix(nrow=ncol(data2m),ncol=nb_topics,data=0)
          # Performing spatial randomizations
          for (k in 1:nb_topics)
          {
            a=nrow_grid
            b=ncol_grid
            while ((a==nrow_grid) && (b==ncol_grid))
            {
              a=sample(1:nrow_grid,1)
              b=sample(1:ncol_grid,1)
            }  
            if ((a!=nrow_grid) && (b!=ncol_grid))
              Spatial_documents_randomized = Spatial_documents[[k]][c((a+1):nrow_grid,1:a),c((b+1):ncol_grid,1:b)] 
            else if (b==ncol_grid)
              Spatial_documents_randomized = Spatial_documents[[k]][c((a+1):nrow_grid,1:a),]
            else if (a==nrow_grid)
              Spatial_documents_randomized = Spatial_documents[[k]][,c((b+1):ncol_grid,1:b)] 
            
            if (byRow)
              documents_randomized[[rndzation]][,k] = as.vector(t(Spatial_documents_randomized))[!is.na(as.vector(t(Spatial_documents_randomized)))]
            else
              documents_randomized[[rndzation]][,k] = Spatial_documents_randomized[!is.na(Spatial_documents_randomized)]
            KL_documents_randomized[[rndzation]][,k] = documents_randomized[[rndzation]][,k]/sum(documents_randomized[[rndzation]][,k])
          }
        }
      }
    }
    
    for (k in 1:nb_topics)
    {
      Cor_test = cor.test(documents[,k],data_abiotic[,j_abiotic])
      Correlation_abiotic[[1]][k,j_abiotic] = Cor_test$estimate
      Correlation_abiotic[[2]][k,j_abiotic] = Cor_test$p.value
      
      if (pca_abiotic)
      {
        Cor_test = cor.test(documents[,k],abiotic_PCA$li[,j_abiotic])
        Correlation_abiotic[[3]][k,j_abiotic] = Cor_test$estimate
        Correlation_abiotic[[4]][k,j_abiotic] = Cor_test$p.value
      }
      
      # Randomizations not implemented for grid && !geographic
      if (grid && !geographic)
      {
        nb_non_significant_rndzations = 0
        #nb_non_significant_rndzations_PCA = 0
        cor_abiotic_comparison_randomized = vector(length=nb_abiotic_rndzations,mode="numeric")
        for (rndzation in 1:nb_abiotic_rndzations) 
        {
          cor_abiotic_comparison_randomized[rndzation] = cor(documents_randomized[[rndzation]][,k],data_abiotic[,j_abiotic])
          
          if (Correlation_abiotic[[1]][k,j_abiotic] < cor_abiotic_comparison_randomized[rndzation])
            nb_non_significant_rndzations = nb_non_significant_rndzations + 1
        }
        Mean_cor_abiotic_comparison_randomized[k,j_abiotic] = mean(cor_abiotic_comparison_randomized)
        Sd_cor_abiotic_comparison_randomized[k,j_abiotic] = sd(cor_abiotic_comparison_randomized)
        
        if (Correlation_abiotic[[1]][k,j_abiotic] > 0)
          p_value_abiotic[k,j_abiotic] = nb_non_significant_rndzations/nb_abiotic_rndzations
        else if (Correlation_abiotic[[1]][k,j_abiotic] < 0)
          p_value_abiotic[k,j_abiotic] = 1 - nb_non_significant_rndzations/nb_abiotic_rndzations 
      }
    }
  }
  # End of the loop over j_abiotic
  
    return(list(Correlation_abiotic=Correlation_abiotic,
                Mean_cor_abiotic_comparison_randomized=Mean_cor_abiotic_comparison_randomized,
                Sd_cor_abiotic_comparison_randomized=Sd_cor_abiotic_comparison_randomized,
                p_value_abiotic=p_value_abiotic,
                colnames_abiotic=colnames_abiotic))
}
