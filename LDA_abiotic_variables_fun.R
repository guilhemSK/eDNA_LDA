# Variables values for testing :
# locally_based = 1
# genotoul_cluster_based = 0
# EDB_cluster_based = 0
# data_pp = 1
# testdata = 0
# data_gs = 0
# case = "lidar"
# local_prefix = "/Users/guilhemsommeria-klein/Desktop/These/"
# data_insert = "Donnees_PetitPlateau"
# filled = 1
# filled_with_gaps =0
# # Random matrix with 1131 columns and 100 rows:
# data2m = replicate(1131,abs(rnorm(100)))
# pca_abiotic = 0
# nb_topics = 3
# nb_abiotic_rndzations = 5
# documents = replicate(nb_topics,abs(rnorm(1131)))
# KL_documents = documents
# occurrence = 1

LDA_abiotic_variables_fun <- function(data_pp,testdata,data_gs,case,local_prefix,data_insert,coordGS,filled,filled_with_gaps,Missing_positions_indices,
                                      ncol_data2m,sum_data2m,pca_abiotic,nb_topics,nb_abiotic_rndzations,documents,KL_documents,occurrence,case_list){
  
  if (case == "chemi")
  {
    setwd(paste(local_prefix,data_insert,"/Chemistery/",sep=""))
    # produces a dataframe with columns as vectors instead of factors (otherwise one cannot direclty apply "as.numeric" to the values)
    data_abiotic = read.table("chemistry_pred_CRIJ.txt",sep=" ",colClasses="vector")
    data_abiotic = apply(data_abiotic,2,as.numeric)
    colnames_abiotic = colnames(data_abiotic)
  } else if (case == "lidar")
  {
    # Lidar data
    setwd(paste(local_prefix,data_insert,"/Lidar/",sep=""))
    data_abiotic_data.frame = read.table("Lidar_locsites_E20_origin_mean_10_extended.csv",colClasses="vector",sep=";")
    data_abiotic = matrix(nrow=length(data_abiotic_data.frame$V5[-1]),ncol=11)
    colnames(data_abiotic) = c("Topography","Canopy height","Light","Wetness","Slope","Canopy height difference 2008-2012","Tree density","Number of tree deaths 2008-2012","Loss of above-ground biomass 2008-2012","Slope standard deviation","Slope (log)")
    colnames_abiotic = colnames(data_abiotic)
    
    # Removing from data_lidar: canop08, which is already contained in dif_can, tree_density2008, AGBloss_mort2008, ndeath2008 
    data_abiotic[,1] = as.vector(data_abiotic_data.frame$V5[-1])
    data_abiotic[,2] = as.vector(data_abiotic_data.frame$V4[-1])
    data_abiotic[,3] = as.vector(data_abiotic_data.frame$V6[-1])
    data_abiotic[,4] = as.vector(data_abiotic_data.frame$V7[-1])
    data_abiotic[,5] = as.vector(data_abiotic_data.frame$V8[-1])
    data_abiotic[,6] = as.vector(data_abiotic_data.frame$V10[-1])
    data_abiotic[,7] = as.vector(data_abiotic_data.frame$V11[-1])
    data_abiotic[,8] = as.vector(data_abiotic_data.frame$V12[-1])
    data_abiotic[,9] = as.vector(data_abiotic_data.frame$V13[-1])
    data_abiotic[,10] = as.vector(data_abiotic_data.frame$V17[-1])
    # log of slope
    data_abiotic[,11] = data_abiotic[,5]
    
    data_abiotic = apply(data_abiotic,2,as.numeric)
    data_abiotic[,10] = sqrt(data_abiotic[,10])
    data_abiotic[,11] = log(data_abiotic[,11])
  } else if (case == "climate")
  {
    # Lidar data
    setwd(paste(local_prefix,data_insert,"/Climate/",sep=""))
    #data_abiotic = GDALinfo("Isothermality.grd")
    isothermality = as.data.frame(raster("Isothermality.grd"), row.names=NULL, optional=T, xy=TRUE, na.rm=T, long=FALSE)
    data_abiotic = matrix(nrow = nrow(isothermality), ncol = 5, data=0)
    data_abiotic[,1] = isothermality$bio3
    coordGS1 = data.frame(x=isothermality$x*10,y=isothermality$y*10)
    data_abiotic[,2] = as.data.frame(raster("Max_Temperature_of_Warmest_Month.grd"), row.names=NULL, optional=T, xy=TRUE, na.rm=T, long=FALSE)$bio5
    data_abiotic[,3] = as.data.frame(raster("Min_Temperature_of_Coldest_Month.grd"), row.names=NULL, optional=T, xy=TRUE, na.rm=T, long=FALSE)$bio6
    data_abiotic[,4] = as.data.frame(raster("Precipitation_of_Driest_Quarter.grd"), row.names=NULL, optional=T, xy=TRUE, na.rm=T, long=FALSE)$bio17
    data_abiotic[,5] = as.data.frame(raster("Precipitation_of_Wettest_Quarter.grd"), row.names=NULL, optional=T, xy=TRUE, na.rm=T, long=FALSE)$bio16
    data_abiotic = apply(data_abiotic,2,as.numeric)          
    
    colnames_abiotic = c("Isothermality","Max. temperature of warmest month","Min. temperature of coldest month","Precipitation of driest quarter","Precipitation of wettest quarter")
    colnames(data_abiotic) = colnames_abiotic
    
    #           pdf("Climate_rasters.pdf")
    #           plot(raster("Isothermality.grd"))
    #           title("Isothermality")
    #           plot(raster("Max_Temperature_of_Warmest_Month.grd"))
    #           title("Max. temperature of warmest month")
    #           plot(raster("Min_Temperature_of_Coldest_Month.grd"))
    #           title("Min. temperature of coldest month")
    #           plot(raster("Precipitation_of_Driest_Quarter.grd"))
    #           title("Precipitation of driest quarter")
    #           plot(raster("Precipitation_of_Wettest_Quarter.grd"))
    #           title("Precipitation of wettest quarter")
    #           dev.off()
    
    # The following replaces "data_abiotic = data_abiotic[-which(Missing_positions_indices==1),]" for data_gs:
    sites_to_remove = rep(1,nrow(coordGS1)) 
    for (i in 1:nrow(coordGS1))
    {
      for (j in 1:nrow(coordGS))
      {
        if ((coordGS1[i,1] == coordGS[j,1]) && (coordGS1[i,2] == coordGS[j,2]))
          sites_to_remove[i] = 0
      }
    }
    data_abiotic = data_abiotic[-which(sites_to_remove==1),]
    
    if (nrow(data_abiotic)!=nrow(coordGS))   
      stop("Error: Abiotic and biotic data do not match")      
  }
  
  #         if (case == "chemi")
  #           browser()
  
  ncol0 = ncol(data_abiotic)
  
  if (data_pp)
  {
    # Switching from indexing row by row (first the 39 elements of the first row, then the 39 elements of the second row, etc.)
    # to indexing column by column (first the 29 elements of the first column, then the 29 elements of the second column, etc.)
    data_abiotic_wrongindexing = data_abiotic
    for (i in 1:29)
    {
      for (j in 1:39)
      {
        data_abiotic[(j-1)*29+i,] = data_abiotic_wrongindexing[(i-1)*39+j,]   
      }
    }
  }
  
  #         if (case == "chemi")
  #           browser()
  
  if ((!filled || filled_with_gaps) && !data_gs)
  {
    data_abiotic = data_abiotic[-which(Missing_positions_indices==1),]
  }
  
  if (nb_topics == 3 && data_pp)
  {
    data_abiotic1 = matrix(nrow = ncol_data2m, ncol = ncol0+3)
    data_abiotic1[,1:ncol0] = data_abiotic
    
    setwd(paste0(local_prefix,data_insert))
    if (!occurrence)
      documents_bacteria = readRDS("documents_bacteria_PP_3topics_best_real.rds")
    else if (occurrence)
      documents_bacteria = readRDS("documents_bacteria_PP_3topics_occ_best_real.rds")
    if (!filled || filled_with_gaps) 
      documents_bacteria = documents_bacteria[-which(Missing_positions_indices==1),] 
    if (case == "lidar") {
      setwd(paste0(local_prefix,data_insert,"/Lidar/"))
    } else if (case == "chemi")
      setwd(paste(local_prefix,data_insert,"/Chemistery/",sep=""))
    data_abiotic1[,(ncol0+1):(ncol0+3)] = documents_bacteria
    colnames(data_abiotic1) = c(colnames(data_abiotic),"Terra firme","Hydromorphic","Exposed rock")
    colnames_abiotic = colnames(data_abiotic1)
    ncol1 = ncol(data_abiotic1)
    data_abiotic = data_abiotic1
    
    if (case == case_list[1])
    {
      KL_documents_bacteria = documents_bacteria
      #           if (occurrence)
      #             KL_documents_bacteria[which(documents_bacteria < 1/2008486)] = 1/2008486
      #           else if (!occurrence)
      #             KL_documents_bacteria[which(documents_bacteria < 1/21850334)] = 1/21850334
      # Bacterial data are applied the same abundance threshold factor as the barcode of interest, so as to have the same minimal value in both distributions
      KL_documents_bacteria[which(documents_bacteria < 1/sum_data2m)] = 1/sum_data2m
      KL_documents_bacteria = sweep(KL_documents_bacteria,MARGIN=2,colSums(KL_documents_bacteria),`/`)
    }
  }
  
  #         if (case == "chemi")
  #           browser()
  
  if (pca_abiotic)
  {
    # PCA on standardized ("normed") data, ie centered and divided by std. dev. for each column 
    abiotic_PCA = dudi.pca(as.data.frame(data_abiotic), row.w = rep(1, nrow(data_abiotic))/nrow(data_abiotic), 
                           col.w = rep(1, ncol(data_abiotic)), center = TRUE, scale = TRUE, 
                           scannf = F, nf = ncol(data_abiotic)) 
  }
  
  Correlation_abiotic = list()
  if (nb_topics == 3 && data_pp)
  {
    Correlation_abiotic[[1]] = matrix(nrow = nb_topics, ncol = ncol1)
    Correlation_abiotic[[2]] = matrix(nrow = nb_topics, ncol = ncol1)
    
    if (case == case_list[1])
    {
      SKL_bacteria = list()
      SKL_bacteria[[1]] = matrix(nrow = nb_topics, ncol = nb_topics)
      SKL_bacteria[[2]] = matrix(nrow = nb_topics, ncol = nb_topics)
    }
  } else
  {
    Correlation_abiotic[[1]] = matrix(nrow = nb_topics, ncol = ncol0)
    Correlation_abiotic[[2]] = matrix(nrow = nb_topics, ncol = ncol0)
  }
  if (pca_abiotic)
  {
    Correlation_abiotic[[3]] = matrix(nrow = nb_topics, ncol = ncol0)
    Correlation_abiotic[[4]] = matrix(nrow = nb_topics, ncol = ncol0)
  }
  
  if (nb_topics == 3 && data_pp)
  {
    Mean_cor_abiotic_comparison_randomized = matrix(nrow=nb_topics,ncol=ncol1,data=0)
    Sd_cor_abiotic_comparison_randomized = matrix(nrow=nb_topics,ncol=ncol1,data=0)
    p_value_abiotic = matrix(nrow=nb_topics,ncol=ncol1,data=0)
    j_abiotic_range = ncol1
    
    if (case == case_list[1])
    {
      Mean_SKL_bacteria_comparison_randomized = matrix(nrow=nb_topics,ncol=nb_topics,data=0)
      Sd_SKL_bacteria_comparison_randomized = matrix(nrow=nb_topics,ncol=nb_topics,data=0)
      p_value_SKL_bacteria = matrix(nrow=nb_topics,ncol=nb_topics,data=0)
    }
  } else
  {
    Mean_cor_abiotic_comparison_randomized = matrix(nrow=nb_topics,ncol=ncol0,data=0)
    Sd_cor_abiotic_comparison_randomized = matrix(nrow=nb_topics,ncol=ncol0,data=0)
    p_value_abiotic = matrix(nrow=nb_topics,ncol=ncol0,data=0)
    j_abiotic_range = ncol0
  }
  
  #         if (case == "chemi")
  #           browser()
  
  for (j_abiotic in 1:j_abiotic_range)
  {
    if ((j_abiotic == 1) && data_pp && (case == case_list[1]))
    {
      if (filled && !filled_with_gaps)
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
            Spatial_documents = matrix(documents[,k],ncol=39)
            a=29
            b=39
            while ((a==29) && (b==39))
            {
              a=sample(1:29,1)
              b=sample(1:39,1)
            }  
            if ((a!=29) && (b!=39))
              Spatial_documents_randomized = Spatial_documents[c((a+1):29,1:a),c((b+1):39,1:b)] 
            else if (b==39)
              Spatial_documents_randomized = Spatial_documents[c((a+1):29,1:a),]
            else if (a==29)
              Spatial_documents_randomized = Spatial_documents[,c((b+1):39,1:b)] 
            documents_randomized[[rndzation]][,k] = as.vector(Spatial_documents_randomized)
            KL_documents_randomized[[rndzation]][,k] = documents_randomized[[rndzation]][,k]/sum(documents_randomized[[rndzation]][,k])
          }
        }
      } else if (!filled || filled_with_gaps)
      { 
        Spatial_documents = vector(length=nb_topics,mode="list")
        for (k in 1:nb_topics)
        {
          Spatial_documents[[k]] = matrix(nrow=29,ncol=39,data=0)
          position_shift = 0
          for (j in 1:39)
          {
            for (i in 1:29)
            {
              if (Missing_positions_indices[(j-1)*29+i]==0)
                Spatial_documents[[k]][i,j] = documents[(j-1)*29+i-position_shift,k]    
              else if (Missing_positions_indices[(j-1)*29+i]==1)
              {
                Spatial_documents[[k]][i,j] = NA
                position_shift = position_shift+1
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
            a=29
            b=39
            while ((a==29) && (b==39))
            {
              a=sample(1:29,1)
              b=sample(1:39,1)
            }  
            if ((a!=29) && (b!=39))
              Spatial_documents_randomized = Spatial_documents[[k]][c((a+1):29,1:a),c((b+1):39,1:b)] 
            else if (b==39)
              Spatial_documents_randomized = Spatial_documents[[k]][c((a+1):29,1:a),]
            else if (a==29)
              Spatial_documents_randomized = Spatial_documents[[k]][,c((b+1):39,1:b)] 
            documents_randomized[[rndzation]][,k] = as.vector(Spatial_documents_randomized[!is.na(Spatial_documents_randomized)])
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
      
      if (pca_abiotic && (j_abiotic < ncol0+1))
      {
        Cor_test = cor.test(documents[,k],abiotic_PCA$li[,j_abiotic])
        Correlation_abiotic[[3]][k,j_abiotic] = Cor_test$estimate
        Correlation_abiotic[[4]][k,j_abiotic] = Cor_test$p.value
      }
      
      # Randomizations not implemented yet for data_gs
      if (data_pp)
      {
        nb_non_significant_rndzations = 0
        #nb_non_significant_rndzations_PCA = 0
        cor_abiotic_comparison_randomized = vector(length=nb_abiotic_rndzations,mode="numeric")
        for (rndzation in 1:nb_abiotic_rndzations) 
        {
          cor_abiotic_comparison_randomized[rndzation] = cor(documents_randomized[[rndzation]][,k],data_abiotic[,j_abiotic])
          
          #               cor_abiotic_comparison_randomized_PCA = cor.test(documents_randomized[[rndzation]][,rev(sort_normal_topic$ix)[k]],lidar_PCA$li[,j])
          #               if (Correlation_lidar[[3]][k,j] < cor_abiotic_comparison_randomized_PCA)
          #                 nb_non_significant_rndzations_PCA = nb_non_significant_rndzations_PCA + 1
          
          
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
    # End loop over j_abiotic  :
  }
  
  if (data_pp && (nb_topics == 3) && (case == case_list[1]))
  {
    for (j_abiotic in 1:nb_topics)
    {
      for (k in 1:nb_topics)
      {
        SKL_bacteria[[1]][k,j_abiotic] = 1/2*(KL.plugin(KL_documents[,k],KL_documents_bacteria[,j_abiotic]) + KL.plugin(KL_documents_bacteria[,j_abiotic],KL_documents[,k]))
        
        nb_non_significant_rndzations_SKL = 0
        SKL_bacteria_comparison_randomized = vector(length=nb_abiotic_rndzations,mode="numeric")
        for (rndzation in 1:nb_abiotic_rndzations) 
        {
          SKL_bacteria_comparison_randomized[rndzation] = 1/2*(KL.plugin(KL_documents_randomized[[rndzation]][,k],KL_documents_bacteria[,j_abiotic]) + KL.plugin(KL_documents_bacteria[,j_abiotic],KL_documents_randomized[[rndzation]][,k]))
          
          #               cor_abiotic_comparison_randomized_PCA = cor.test(documents_randomized[[rndzation]][,rev(sort_normal_topic$ix)[k]],lidar_PCA$li[,j])
          #               if (Correlation_lidar[[3]][k,j] < cor_abiotic_comparison_randomized_PCA)
          #                 nb_non_significant_rndzations_PCA = nb_non_significant_rndzations_PCA + 1
          
          if (SKL_bacteria[[1]][k,j_abiotic] > SKL_bacteria_comparison_randomized[rndzation])
            nb_non_significant_rndzations_SKL = nb_non_significant_rndzations_SKL + 1
        }
        Mean_SKL_bacteria_comparison_randomized[k,j_abiotic] = mean(SKL_bacteria_comparison_randomized)
        Sd_SKL_bacteria_comparison_randomized[k,j_abiotic] = sd(SKL_bacteria_comparison_randomized)
        SKL_bacteria[[2]][k,j_abiotic] = (Mean_SKL_bacteria_comparison_randomized[k,j_abiotic] - SKL_bacteria[[1]][k,j_abiotic])/Mean_SKL_bacteria_comparison_randomized[k,j_abiotic]
        
        if (SKL_bacteria[[2]][k,j_abiotic] > 0)
          p_value_SKL_bacteria[k,j_abiotic] = nb_non_significant_rndzations_SKL/nb_abiotic_rndzations
        else if (SKL_bacteria[[2]][k,j_abiotic] < 0)
          p_value_SKL_bacteria[k,j_abiotic] = 1 - nb_non_significant_rndzations_SKL/nb_abiotic_rndzations 
      }
    }
  }
  
  if (data_pp && (nb_topics == 3))
    return(list(Correlation_abiotic=Correlation_abiotic,ncol0=ncol0,Mean_cor_abiotic_comparison_randomized=Mean_cor_abiotic_comparison_randomized,
                Sd_cor_abiotic_comparison_randomized=Sd_cor_abiotic_comparison_randomized,p_value_abiotic=p_value_abiotic,colnames_abiotic=colnames_abiotic,
                SKL_bacteria=SKL_bacteria,p_value_SKL_bacteria=p_value_SKL_bacteria))
  else
    return(list(Correlation_abiotic=Correlation_abiotic,ncol0=ncol0,Mean_cor_abiotic_comparison_randomized=Mean_cor_abiotic_comparison_randomized,
           Sd_cor_abiotic_comparison_randomized=Sd_cor_abiotic_comparison_randomized,p_value_abiotic=p_value_abiotic,colnames_abiotic=colnames_abiotic))
  
  #End function
}





