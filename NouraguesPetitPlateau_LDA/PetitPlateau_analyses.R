
library(ggplot2)
library(gridExtra)
library(vegan)
library(gstat)
library(ape)

test_data = 0
test_error = 0

if (test_data)
{
  taxo_vect = c("Test_data/Continuous-mixed_samples_nbtopics5_nbmotus1000_randomtopics_sampledreads","Test_data/Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics0.02_1000sampledreads")
  taxo_names = c("Uniformly-distributed assemblages","Dirichlet-distributed assemblages")
  diversity = c(1000,1000)
} else
{
  taxo_vect = c("Bacteries_16S","Protistes_18S","Champignons_18S","Arthropodes_18S","Nematodes_18S","Platyhelminthes_18S","Annelides_18S")
  taxo_names = c("Bacteria 16S","Protists 18S","Fungi 18S","Arthropods 18S","Nematodes 18S","Platyhelminthes 18S","Annelids 18S")
  diversity = c(20165,1650,1142,1882,378,126,51)
  optimalK = c(5,2,4,2,2,2,2)
}
# mean_sim = matrix(nrow = length(taxo_insert_vect), ncol = 3, dimnames = list(taxo_insert_vect,inference_method), data=0)
# sim_intercept = matrix(nrow = length(taxo_insert_vect), ncol = 3, dimnames = list(taxo_insert_vect,inference_method), data=0)

nrow_grid = 39
ncol_grid = 29

figure_folder = "/Users/guilhemsommeria-klein/Desktop/Manuscrits/ThesisPaper2_Topic_Modelling_PP/Figures/Submission2_analyses/"

optimalK_comput = 1
if (test_data)
{
  topic_range = "2-8"
} else 
  topic_range = "2-15"
# topic_range = "2-15"
# topic_range = "2-30"
optimalK_mpar_perplexity = 1
optimalK_mpar_aic = 1
spline_df = 5

variableK = 0
nb_topics = 3

environmental_data = 0
abiotic_pca = 1
spatial_randomizations = 0
nb_abiotic_rndzations = 10000
maps = 0
test_plots = 0

spatial_structure = 0

assemblage_diversity = 0

stability = 0

if (test_data)
{
  variableK_insert = ""
} else
{
  if (variableK)
  {
    variableK_insert = "_optimalK"
  } else
    variableK_insert = "_K=3"
}
if (optimalK_comput)
{
  if (optimalK_mpar_perplexity)
  {
    K.perplexity.folder_name = c(paste0("/Rtopicmodels_LDA_VEM_alpha0.1_nb_topics",topic_range,"_post-predictive-cross-valid_fold_size10_em_tol1e-06_var_tol1e-08_occurrence/"),
                      paste0("/Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics",topic_range,"_post-predictive-cross-valid_nb_iter1000_thin25_burnin2000_fold_size10_occurrence/"),
                      paste0("/Rtopicmodels_LDA_VEM_alpha0.1_nb_topics",topic_range,"_post-predictive-cross-valid_fold_size10_em_tol1e-06_var_tol1e-08/"),
                      paste0("/Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics",topic_range,"_post-predictive-cross-valid_nb_iter1000_thin25_burnin2000_fold_size10/"))
    
    optimalK_min.crossValid = list()
    # optimalK_crossValid_corCumul = list()
    plot.spline.perplexity = list()
    plot.spline.perplexity.mean.median = list()
    # plot.spline.perplexity.derivative = list()
    # plot.spline.perplexity.derivative2 = list()
    # plot.spline.perplexity.corrCumul = list()
  }
  if (optimalK_mpar_aic)
  {
    K.aic.folder_name = c(paste0("/Rtopicmodels_LDA_VEM_alpha0.1_nb_topics",topic_range,"_elbow-AIC_nb_real",ifelse(test_data,50,10),"_em_tol1e-06_var_tol1e-08_occurrence/"),
                          paste0("/Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics",topic_range,"_elbow-AIC_nb_iter1000_thin25_burnin2000_nb_real50_occurrence/"),
                      paste0("/Rtopicmodels_LDA_VEM_alpha0.1_nb_topics",topic_range,"_elbow-AIC_nb_real",ifelse(test_data,50,10),"_em_tol1e-06_var_tol1e-08/"),
                      paste0("/Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics",topic_range,"_elbow-AIC_nb_iter1000_thin25_burnin2000_nb_real50/"))
    
    optimalK_aic = list()
    # optimalK_crossValid_corCumul = list()
    plot.spline.aic = list()
  }
  if (optimalK_mpar_perplexity && optimalK_mpar_aic)
  {
    plot.spline.both = list()
  }
}
if (environmental_data || spatial_structure)
{
  if (spatial_structure)
  {
    I_linear.observed_w.mean = I_inverse.observed_w.mean = I_square.observed_w.mean = I_cubic.observed_w.mean = list()
    I_linear.p.value_w.mean = I_inverse.p.value_w.mean = I_square.p.value_w.mean = I_cubic.p.value_w.mean = list()
    
    # Factor 10 is to obtain the distances in meters, since there are 10 m between the nearest sampling points:
    coord = expand.grid(x = 10*seq(1,ncol_grid,1), y = 10*seq(1,nrow_grid,1))
    spatial_distances = as.matrix(dist(coord, method = "euclidean"))
  }
  if (assemblage_diversity)
  {
    nb.OTUS = list()
  }
  if (environmental_data)
  {
    env_regression = list()
    env_selected = list()
    env_regression.axes = list()
    env_regression.topics = list()
    env_selected.topics = list()
    env_regression.topics.axes.adjR2 = list()
    env_regression.topics.axes.sign = list()
    env_regression.topics.axes.pval = list()
    
    lidar_data_path = "/Users/guilhemsommeria-klein/Desktop/These/Donnees_PetitPlateau/Lidar/"
    lidar_data_file = "Lidar_locsites_E20_origin_mean_10_extended.csv"
    
    # Lidar data
    lidar_data = read.table(paste0(lidar_data_path,lidar_data_file), sep=";", colClasses="vector", row.names = 1, header = T)
    # Removing from data_lidar: canop08, which is already contained in dif_can, tree_density2008, AGBloss_mort2008, ndeath2008: 
    lidar_data = apply(lidar_data[,!colnames(lidar_data) %in% c("ym","xm","dif_can","canop08","tree_density2008","AGBloss_mort2008","ndeath2008")],2,as.numeric)
    # Taking the square root of the slope variance:
    lidar_data[,9] = sqrt(lidar_data[,9])
    # data_abiotic[,5] = log(data_abiotic[,5])
    # colnames(lidar_data) = c("Canopy.height","Topography","Light","Wetness","Slope","Canopy.height.difference.2008-2012","Tree.density","Number.of.tree.deaths.2008-2012","Loss.of.above-ground.biomass.2008-2012","Slope.standard.deviation")
    colnames(lidar_data) = c("Canopy.height","Topography","Light","Wetness","Slope","Tree.density","Number.of.tree.deaths.2008.2012","Loss.of.above.ground.biomass.2008.2012","Slope.standard.deviation")
    # Reordering the stations to match documents:
    for (j in 1:ncol(lidar_data))
      lidar_data[,j]  = as.vector(t(matrix(data=lidar_data[,j],nrow=nrow_grid)))
    
    if (abiotic_pca)
    {
      pca_insert = ".pca"
      library(ade4)
      lidar_pca = dudi.pca(as.data.frame(lidar_data), row.w = rep(1, nrow(lidar_data))/nrow(lidar_data), 
                           col.w = rep(1, ncol(lidar_data)), center = TRUE, scale = TRUE, 
                           scannf = F, nf = ncol(lidar_data))
      lidar_data = lidar_pca$li[,1:5]
    } else
      pca_insert = ""
  }
}
if (stability)
{
  mean_sim = list()
  mean_sim.k = list()
  sim_intercept = list()
  sim_intercept.k = list()
  alpha_all_real = list()
  stability.plot.llh.diff = list()
  stability.plot.rank = list()
}

for (case in c(1,3))
# for (case in 1:4)
{
  if (case %in% c(1,3))
  {
    if (optimalK_comput)
      inference_insert = "_VEM_EM.tol1e-06Var.tol1e-08"
    else
      inference_insert = "_VEM_EM.tol1e-07Var.tol1e-08"
    if (case == 1)
      occ_insert = "_occurrence"
    else if (case == 3)
      occ_insert = ""
  } else if (case %in% c(2,4))
  {
    inference_insert = "_Gibbs_iter1000thin25burnin2000"
    if (case == 2)
      occ_insert = "_occurrence"
    else if (case == 4)
      occ_insert = ""
  }
  
  if (test_data)
  {
    nb_real = 50
  } else
  {
    nb_real = c(10,50,10,50)[case]
  }
  
  if (optimalK_comput)
  {
    if (optimalK_mpar_perplexity)
    {
      optimalK_min.crossValid[[case]] = matrix(nrow = length(taxo_vect), ncol = 2, dimnames = list(taxo_names,c("Mean","Median")), data = 0)
      # optimalK_crossValid_corCumul[[case]] = matrix(nrow = length(taxo_vect), ncol = 3, data = 0)
      plot.spline.perplexity[[case]] = list()
      # plot.spline.perplexity.mean.median[[case]] = list()
      # plot.spline.perplexity.derivative[[case]] = list()
      # plot.spline.perplexity.derivative2[[case]] = list()
      # plot.spline.perplexity.corrCumul[[case]] = list()
    }
    if (optimalK_mpar_aic)
    {
      optimalK_aic[[case]] = matrix(nrow = length(taxo_vect), ncol = 2, dimnames = list(taxo_names,c("Mean","Median")), data = 0)
      plot.spline.aic[[case]] = list()
    }
    if (optimalK_mpar_aic && optimalK_mpar_perplexity)
    {
      plot.spline.both[[case]] = list()
    }
  }
  if (environmental_data)
  {
    # env_regression[[case]] = matrix(nrow = 7, ncol = 5, data = NA, dimnames = list(taxo_names,c("adj.R.squared","R.squared","p.value","spatial.p.value.adj","spatial.p.value")))
    env_regression[[case]] = matrix(nrow = 7, ncol = 3, data = NA, dimnames = list(taxo_names,c("adj.R.squared","R.squared","p.value")))
    env_selected[[case]] = env_regression.axes[[case]] = env_regression.topics[[case]] = 
      env_selected.topics[[case]] = env_regression.topics.axes.adjR2[[case]] = env_regression.topics.axes.pval[[case]] = env_regression.topics.axes.sign[[case]] = vector(length = 7, mode = "list")
    names(env_selected[[case]])  = names(env_regression.axes[[case]]) = names(env_regression.topics[[case]]) = names(env_selected.topics[[case]]) = taxo_names
  }
  if (spatial_structure)
  {
    I_linear.observed_w.mean[[case]] = I_inverse.observed_w.mean[[case]] = I_square.observed_w.mean[[case]] = I_cubic.observed_w.mean[[case]] = vector(length = 7, mode = "list")  
    names(I_linear.observed_w.mean[[case]]) = names(I_inverse.observed_w.mean[[case]]) = names(I_square.observed_w.mean[[case]]) = names(I_cubic.observed_w.mean[[case]]) = taxo_names
    I_linear.p.value_w.mean[[case]] = I_inverse.p.value_w.mean[[case]] = I_square.p.value_w.mean[[case]] = I_cubic.p.value_w.mean[[case]] = vector(length = 7, mode = "list")  
    names(I_linear.p.value_w.mean[[case]]) = names(I_inverse.p.value_w.mean[[case]]) = names(I_square.p.value_w.mean[[case]]) = names(I_cubic.p.value_w.mean[[case]]) = taxo_names
  }
  if (assemblage_diversity)
  {
    nb.OTUS[[case]] = vector(length = length(taxo_vect), mode = "list")
    names(nb.OTUS[[case]]) = taxo_names
  }
  if (stability)
  {
    mean_sim[[case]] = vector(length = length(taxo_vect), mode = "numeric")
    names(mean_sim[[case]]) = taxo_names
    mean_sim.k[[case]] = vector(length = length(taxo_vect), mode = "list")
    names(mean_sim.k[[case]]) = taxo_names
    sim_intercept[[case]] = vector(length = length(taxo_vect), mode = "numeric")
    names(sim_intercept[[case]]) = taxo_names
    sim_intercept.k[[case]] = vector(length = length(taxo_vect), mode = "list")
    names(sim_intercept.k[[case]]) = taxo_names
    if (case %in% c(1,3))
      alpha_all_real[[case]] = matrix(nrow = 100, ncol = length(taxo_vect), dimnames = list(1:100,taxo_names), data = 0)
    stability.plot.llh.diff[[case]] = list()
    stability.plot.rank[[case]] = list()
  }
  i_taxon = 0
  for (taxon in taxo_vect)
  {
    i_taxon = i_taxon+1
    cat(taxon,"\n")
    data.folder_name = paste0("/Users/guilhemsommeria-klein/Desktop/Serveur_Bioclust.data/Donnees_PetitPlateau/",taxon)
    
    if (optimalK_comput)
    {
      if (optimalK_mpar_aic)
      {
        K.folder_name_case = K.aic.folder_name[case]
        data.folder_name_case = paste0(data.folder_name,K.folder_name_case)
        aic_file = paste0(data.folder_name_case,"AIC-llh.Rdata")
        if (file.exists(aic_file))
        {
          load(aic_file)
          
          if (case %in% c(1,2))
          {
            nb_topics_range = 2:8
            AIC = AIC[,1:7]
          }
          
          spline_mean = smooth.spline(nb_topics_range,colMeans(AIC),df=spline_df)
          spline_median = smooth.spline(nb_topics_range,apply(AIC,2,median),df=spline_df)
          
          index = sort.int(spline_mean$y,decreasing = F,index.return = T)$ix[1]
          optimalK_aic[[case]][i_taxon,1] = nb_topics_range[index]
          
          index = sort.int(spline_median$y,decreasing = F,index.return = T)$ix[1]
          optimalK_aic[[case]][i_taxon,2] = nb_topics_range[index]
          
          # llh = LLH_final0[1,]
          # sd_llh = LLH_final0[2,]
          
          # llh_allpoints = as.vector(LLH_final)
          nb_topics_range_allpoints = unlist(lapply(nb_topics_range,rep,nb_real))
          
          plot.spline.aic[[case]][[i_taxon]] = ggplot(data=data.frame(K=nb_topics_range_allpoints,llh=as.vector(AIC))) +
            geom_boxplot(aes(x=as.factor(K),y=llh), size = 0.4, outlier.size = 0.4, outlier.color = "darkgrey") +
            geom_line(data = data.frame(K = seq_along(nb_topics_range), llh = spline_median$y), aes(x = K, y = llh), size = 0.4) +
            geom_vline(xintercept = which(nb_topics_range == optimalK_aic[[case]][i_taxon,2]), linetype = "dashed", size = 0.4) +
            scale_x_discrete(breaks = seq(from = nb_topics_range[1], to = nb_topics_range[length(nb_topics_range)], by = 3)) +
            scale_y_continuous(labels = scales::scientific) +
            theme_bw() +
            theme(axis.title=element_text(size=11),
                  axis.text.x=element_text(size=11),
                  axis.text.y=element_text(size=8),
                  plot.title=element_text(hjust=0, size=12),
                  plot.margin=unit(c(1,2,0.1,0.1),"mm")) +
            labs(x="Number K of assemblages", y="AIC")
          
          
          # With bacteria:
          if (optimalK_mpar_perplexity)
            plot.spline.both[[case]][[2*i_taxon-1]] = plot.spline.aic[[case]][[i_taxon]] + ggtitle(paste(LETTERS[2*i_taxon-1],"-",taxo_names[i_taxon]))
          # No bacteria:
          # if (optimalK_mpar_perplexity && i_taxon>1)
          #   plot.spline.both[[case]][[2*i_taxon-1]] = plot.spline.aic[[case]][[i_taxon]] + ggtitle(paste(LETTERS[2*(i_taxon-1)-1],"-",taxo_names[i_taxon]))
          # Bacteria only for AIC:
          # if (optimalK_mpar_perplexity)
          # {
          #   if (i_taxon == 1)
          #     plot.spline.both[[case]][[1]] = plot.spline.aic[[case]][[i_taxon]] + ggtitle(paste(LETTERS[1],"-",taxo_names[i_taxon]))
          #   else
          #     plot.spline.both[[case]][[2*(i_taxon-1)]] = plot.spline.aic[[case]][[i_taxon]] + ggtitle(paste(LETTERS[2*(i_taxon-1)],"-",taxo_names[i_taxon]))
          # }
          
          plot.spline.aic[[case]][[i_taxon]] = plot.spline.aic[[case]][[i_taxon]] + ggtitle(paste(LETTERS[i_taxon],"-",taxo_names[i_taxon]))
        } else 
        {
          cat("Missing result\n")
          optimalK_aic[[case]][i_taxon,] = c(NA,NA)
          plot.spline.aic[[case]][[i_taxon]] = NA
        }
      }
      
      if (optimalK_mpar_perplexity)
      {
        K.folder_name_case = K.perplexity.folder_name[case]
        data.folder_name_case = paste0(data.folder_name,K.folder_name_case)
        perplexity_file = paste0(data.folder_name_case,"perplexity.Rdata")
        if (file.exists(perplexity_file))
        {
          load(perplexity_file)
          if (nrow(perplexity_mpar) > 1)
          {
            if (case %in% c(1,2))
            {
              if (i_taxon == 1 && case == 1)
              {
                perplexity_mpar0 = perplexity_mpar[,1:3]
                
                K.folder_name_case = paste0("/Rtopicmodels_LDA_VEM_alpha0.1_nb_topics5-15_post-predictive-cross-valid_fold_size100_em_tol1e-06_var_tol1e-08_occurrence/") 
                data.folder_name_case = paste0(data.folder_name,K.folder_name_case)
                perplexity_file = paste0(data.folder_name_case,"perplexity.Rdata")
                load(perplexity_file)
                perplexity_mpar = cbind(perplexity_mpar0,rbind(perplexity_mpar[,1:2],matrix(nrow = nrow(perplexity_mpar0) - nrow(perplexity_mpar), ncol = 2, data = NA)))
                nb_topics_range = 2:6
              } else
              {
                nb_topics_range = 2:8
                perplexity_mpar = perplexity_mpar[,1:7]
              }
            }
            
            spline_mean = smooth.spline(nb_topics_range,colMeans(perplexity_mpar,na.rm=T),df=spline_df)
            spline_median = smooth.spline(nb_topics_range,apply(perplexity_mpar,2,median,na.rm=T),df=spline_df)
            
            index = sort.int(spline_mean$y,decreasing = F,index.return = T)$ix[1]
            optimalK_min.crossValid[[case]][i_taxon,1] = nb_topics_range[index]
            
            index = sort.int(spline_median$y,decreasing = F,index.return = T)$ix[1]
            optimalK_min.crossValid[[case]][i_taxon,2] = nb_topics_range[index]
            
            plot.spline.perplexity[[case]][[i_taxon]] = ggplot(data = data.frame(K = unlist(lapply(nb_topics_range,rep,nrow(perplexity_mpar))), llh = as.vector(perplexity_mpar)))
            if (case == 1)
            {
              plot.spline.perplexity[[case]][[i_taxon]] = plot.spline.perplexity[[case]][[i_taxon]] +
                ylim(min(perplexity_mpar),if (i_taxon == 1) 8500 else if (i_taxon == 2) 1000 else if (i_taxon == 4) 1200 else if (i_taxon == 5) 500 else if (i_taxon == 6) 60 else if (i_taxon == 7) 35 else max(perplexity_mpar))
            } else if (case == 3)
              plot.spline.perplexity[[case]][[i_taxon]] = plot.spline.perplexity[[case]][[i_taxon]] +
                ylim(min(perplexity_mpar),if (i_taxon == 2) 250 else if (i_taxon == 3) 75 else if (i_taxon == 4) 325 else if (i_taxon == 5) 75 else if (i_taxon == 6) 22 else if (i_taxon == 7) 7 else max(perplexity_mpar))
            plot.spline.perplexity[[case]][[i_taxon]] = plot.spline.perplexity[[case]][[i_taxon]] +
              geom_boxplot(aes(x=as.factor(K),y=llh), size = 0.4, outlier.size = 0.4, outlier.color = "darkgrey") +
              geom_line(data = data.frame(K = seq_along(nb_topics_range), llh = spline_median$y), aes(x = K, y = llh), size = 0.4) +
              # geom_line(data = data.frame(K = seq_along(nb_topics_range), llh = apply(perplexity_mpar,2,median,na.rm=T)), aes(x = K, y = llh), size = 0.4, linetype = "dashed") +
              geom_vline(xintercept = which(nb_topics_range == optimalK_min.crossValid[[case]][i_taxon,2]), linetype = "dashed", size = 0.4) +
              # geom_vline(xintercept = which(nb_topics_range == optimalK_min.crossValid[[case]][i_taxon,1]), linetype = "dashed", size = 0.4, col = "blue") +
              scale_x_discrete(breaks = seq(from = nb_topics_range[1], to = nb_topics_range[length(nb_topics_range)], by = 3)) +
              theme_bw() +
              theme(axis.title=element_text(size=11),
                    axis.text=element_text(size=11),
                    plot.title=element_text(hjust=0, size=12),
                    plot.margin=unit(c(1,1.1,0.1,0.1),"mm")) +
              labs(x="Number K of assemblages", y="Cross-validation perplexity")
            
              # With bacteria:
              if (optimalK_mpar_aic)
                plot.spline.both[[case]][[2*i_taxon]] = plot.spline.perplexity[[case]][[i_taxon]] + ggtitle(paste(LETTERS[2*i_taxon],"-",taxo_names[i_taxon]))
              # No bacteria:
              # if (optimalK_mpar_aic && i_taxon>1)
               # plot.spline.both[[case]][[2*i_taxon]] = plot.spline.perplexity[[case]][[i_taxon]] + ggtitle(paste(LETTERS[2*(i_taxon-1)],"-",taxo_names[i_taxon]))
              # Bacteria only for AIC:
              # if (optimalK_mpar_aic && i_taxon>1)
              #   plot.spline.both[[case]][[2*i_taxon - 1]] = plot.spline.perplexity[[case]][[i_taxon]] + ggtitle(paste(LETTERS[2*i_taxon - 1],"-",taxo_names[i_taxon]))
            
            plot.spline.perplexity[[case]][[i_taxon]] = plot.spline.perplexity[[case]][[i_taxon]] + ggtitle(paste(LETTERS[i_taxon],"-",taxo_names[i_taxon]))
            
            # plot.spline.perplexity.mean.median[[case]][[i_taxon]] = ggplot(data = data.frame(K = unlist(lapply(nb_topics_range,rep,nrow(perplexity_mpar))), llh = as.vector(perplexity_mpar))) +
            #   ylim(range(c(spline_mean$y,spline_median$y))) +
            #   geom_boxplot(aes(x=as.factor(K),y=llh), size = 0.4, outlier.size = 0.4, outlier.color = "darkgrey") +
            #   # geom_boxplot(aes(x=as.factor(K),y=llh), size = 0.4, outlier.size = 0.4, outlier.shape = NA) +
            #   # geom_point(aes(x=K,y=llh), size = 0.4) +
            #   geom_line(data = data.frame(K = seq_along(nb_topics_range), llh = spline_mean$y), aes(x = K, y = llh), size = 0.4) +
            #   geom_line(data = data.frame(K = seq_along(nb_topics_range), llh = spline_median$y), aes(x = K, y = llh), size = 0.4, col = "blue") +
            #   geom_line(data = data.frame(K = seq_along(nb_topics_range), llh = apply(perplexity_mpar,2,median)), aes(x = K, y = llh), size = 0.4, col = "blue", linetype = "dashed") +
            #   geom_vline(xintercept = which(nb_topics_range == optimalK_min.crossValid[[case]][i_taxon,1]), linetype = "dotted", size = 0.4) +
            #   geom_vline(xintercept = which(nb_topics_range == optimalK_min.crossValid[[case]][i_taxon,2]), linetype = "dotted", size = 0.4, col = "blue") +
            #   scale_x_discrete(breaks = seq(from = nb_topics_range[1], to = nb_topics_range[length(nb_topics_range)], by = 3)) +
            #   theme_bw() +
            #   theme(axis.title=element_text(size=9),
            #         axis.text=element_text(size=9),
            #         plot.title=element_text(hjust=0, size=10),
            #         plot.margin=unit(c(1,1,1,0.5),"mm")) +
            #   ggtitle(paste(taxon,"-",diversity[i_taxon],"OTUs")) +
            #   labs(x="Number K of assemblages", y=paste0(spline_df,"-df-splined 10-sample-fold perplexity"))
            
            # plot.spline.perplexity.mean.median[[case]][[i_taxon]] = ggplot(data=data.frame(x=nb_topics_range,y1=spline_mean$y,y2=spline_median$y,y3=colMeans(perplexity_mpar),y4=apply(perplexity_mpar,2,median))) +
            #   geom_line(aes(x,y1)) +
            #   geom_line(aes(x,y3),linetype="dashed") +
            #   geom_vline(xintercept = optimalK_min.crossValid[[case]][i_taxon,1], linetype = "dotted", size = 0.4) +
            #   scale_x_discrete(breaks = seq(from = nb_topics_range[1], to = nb_topics_range[length(nb_topics_range)], by = 3)) +
            #   ylim(range(spline_mean$y)) +
            #   theme_bw() +
            #   theme(axis.title=element_text(size=11),
            #         axis.text=element_text(size=11),
            #         plot.title=element_text(hjust=0, size=14),
            #         plot.margin=unit(c(0.1,0.1,0.1,0.1),"mm")) +
            #   ggtitle(paste(LETTERS[i_taxon],"-",taxo_names[i_taxon])) +
            #   labs(x="Number K of assemblages", y="Cross-validation perplexity")
            
            # plot.spline.perplexity.corrCumul[[case]][[i_taxon]] = ggplot(data = data.frame(K = unlist(lapply(nb_topics_range,rep,nrow(perplexity_mpar))), llh = as.vector(perplexity_mpar))) +
            #   geom_boxplot(aes(x=as.factor(K),y=llh), size = 0.4, outlier.size = 0.4, outlier.color = "darkgrey") +
            #   geom_line(data = data.frame(K = seq_along(nb_topics_range), llh = spline_mean$y), aes(x = K, y = llh), size = 0.4) +
            #   geom_vline(xintercept = which(nb_topics_range == optimalK_min.crossValid[[case]][i_taxon,1]), linetype = "dotted", size = 0.4) +
            #   geom_vline(xintercept = which(nb_topics_range == optimalK_crossValid_corCumul[[case]][i_taxon,1]), linetype = "longdash", size = 0.4) +
            #   geom_vline(xintercept = which(nb_topics_range == optimalK_crossValid_corCumul[[case]][i_taxon,2]), linetype = "longdash", size = 0.4, col = "blue") +
            #   geom_vline(xintercept = which(nb_topics_range == optimalK_crossValid_corCumul[[case]][i_taxon,3]), linetype = "longdash", size = 0.4, col = "green") +
            #   scale_x_discrete(breaks = seq(from = nb_topics_range[1], to = nb_topics_range[length(nb_topics_range)], by = 3)) +
            #   theme_bw() +
            #   theme(axis.title=element_text(size=9),
            #         axis.text=element_text(size=9),
            #         plot.title=element_text(hjust=0, size=10),
            #         plot.margin=unit(c(1,1,1,0.5),"mm")) +
            #   ggtitle(paste(taxon,"-",diversity[i_taxon],"OTUs")) +
            #   labs(x="Number K of assemblages", y=paste0(spline_df,"-df-splined 10-sample-fold perplexity"))
            # 
            # plot.spline.perplexity.derivative[[case]][[i_taxon]] = ggplot(data=data.frame(x=coor_derivative,y1=derivative[,1],y2=derivative[,2])) +
            #   geom_line(aes(x,y1)) +
            #   geom_line(aes(x,y2),col="blue") +
            #   geom_vline(xintercept = optimalK_min.crossValid[[case]][i_taxon,1], linetype = "dotted", size = 0.4) +
            #   geom_vline(xintercept = optimalK_min.crossValid[[case]][i_taxon,2], linetype = "dotted", size = 0.4, col = "blue") +
            #   theme_bw() +
            #   ggtitle(paste(taxon,"-",diversity[i_taxon],"OTUs")) +
            #   theme(axis.title=element_text(size=9),
            #         plot.title=element_text(hjust=0, size=15),
            #         plot.margin=unit(c(1,1,1,0.5),"mm")) +
            #   labs(x="Number K of assemblages", y=paste0("1st derivative of ",spline_df,"-df-splined 10-sample-fold perplexity"))
            # 
            # plot.spline.perplexity.derivative2[[case]][[i_taxon]] = ggplot(data=data.frame(x=coor_derivative2,y1=derivative2[,1],y2=derivative2[,2])) +
            #   geom_line(aes(x,y1)) +
            #   geom_line(aes(x,y2),col="blue") +
            #   theme_bw() +
            #   ggtitle(paste(taxon,"-",diversity[i_taxon],"OTUs")) +
            #   theme(axis.title=element_text(size=9),
            #         plot.title=element_text(hjust=0, size=15),
            #         plot.margin=unit(c(1,1,1,0.5),"mm")) +
            #   labs(x="Number K of assemblages", y=paste0("2nd derivative of ",spline_df,"-df-splined 10-sample-fold perplexity"))
          } else
          {
            optimalK_min.crossValid[[case]][i_taxon,] = c(NA,NA)
            # optimalK_crossValid_corCumul[[case]][i_taxon,] = c(NA,NA,NA)
            # plot.spline.perplexity.corrCumul[[case]][[i_taxon]] = NA
            plot.spline.perplexity[[case]][[i_taxon]] = NA
            plot.spline.perplexity.mean.median[[case]][[i_taxon]] = NA
            # plot.spline.perplexity.derivative[[case]][[i_taxon]] = NA
            # plot.spline.perplexity.derivative2[[case]][[i_taxon]] = NA
          }
        } else 
        {
          cat("Missing result\n")
          optimalK_min.crossValid[[case]][i_taxon,] = c(NA,NA)
          # optimalK_crossValid_corCumul[[case]][i_taxon,] = c(NA,NA,NA)
          # plot.spline.perplexity.corrCumul[[case]][[i_taxon]] = NA
          plot.spline.perplexity[[case]][[i_taxon]] = NA
          # plot.spline.perplexity.mean.median[[case]][[i_taxon]] = NA
          # plot.spline.perplexity.derivative[[case]][[i_taxon]] = NA
          # plot.spline.perplexity.derivative2[[case]][[i_taxon]] = NA
        }
      }
    }
    
    if (environmental_data || spatial_structure || assemblage_diversity)
    {
      if (variableK)
      {
        nb_topics = optimalK[i_taxon]
      }
      
      K.folder_name = c(paste0("/Rtopicmodels_LDA_VEM_nb_topics",nb_topics,"_nb_real100_em_tol1e-07_var_tol1e-08_best_keep_occurrence/"),
                        paste0("/Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics",nb_topics,"_nb_iter1000_nb_real100_meanPosteriorDistributedLlh_thin25_burnin2000_occurrence/"),
                        paste0("/Rtopicmodels_LDA_VEM_nb_topics",nb_topics,"_nb_real100_em_tol1e-07_var_tol1e-08_best_keep/"),
                        paste0("/Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics",nb_topics,"_nb_iter1000_nb_real100_meanPosteriorDistributedLlh_thin25_burnin2000/"))
      Result.file_name = c(paste0("Rtopicmodels_LDA_VEM_nb_topics",nb_topics,"_nb_real100_em_tol1e-07_var_tol1e-08_best_keep_occurrence.Rdata"),
                           paste0("Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics",nb_topics,"_nb_iter1000_nb_real100_meanPosteriorDistributedLlh_thin25_burnin2000_occurrence.Rdata"),
                           paste0("Rtopicmodels_LDA_VEM_nb_topics",nb_topics,"_nb_real100_em_tol1e-07_var_tol1e-08_best_keep.Rdata"),
                           paste0("Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics",nb_topics,"_nb_iter1000_nb_real100_meanPosteriorDistributedLlh_thin25_burnin2000.Rdata"))
      
      K.folder_name_case = K.folder_name[case]
      data.folder_name_case = paste0(data.folder_name,K.folder_name_case)
      # if (case %in% c(2,4))
      #   documents_file_name = paste0(data.folder_name_case,"1st_closestToMean_realization/Spatial_topicmix_kriged.rds")
      # else if (case %in% c(1,3))
      #   documents_file_name = paste0(data.folder_name_case,"1st_best_realization/Spatial_topicmix_kriged.rds")
      result.file = paste0(data.folder_name_case,Result.file_name[case])
      if (file.exists(result.file))
      {
        # spatial_topicmix_kriged = readRDS(documents_file_name)
        # # selecting the z.pred columns in all topics:
        # documents = unlist(lapply(spatial_topicmix_kriged,function(g) g$z.pred))
        # # setting one topic per column
        # documents = matrix(documents,ncol=nb_topics,dimnames=list(rownames(spatial_topicmix_kriged[[1]]),paste0("assemblage",1:nb_topics)))
        
        load(result.file)
        Ordered_realizations = readRDS(paste0(data.folder_name_case,"Ordered_realizations.rds"))
        # One topic per columns; sums to one over each row/site:
        documents = Result[[Ordered_realizations$ix[ifelse(case == 1 && taxon == "Platyhelminthes_18S",8,1)]]]@gamma
        sort_normal_topic = sort.int(apply(documents,2,mean), index.return=T, decreasing = T)
        documents = documents[,sort_normal_topic$ix]
        
        # Removing missing stations from lidar_data:
        load(paste0(data.folder_name,"/data2m_filled.Rdata"))
        Missing_positions_indices = vector(length = ncol(data2m), mode = "numeric")
        missing = F
        cat(taxo_names[i_taxon],ncol(data2m),"\n")
        if (length(which(colSums(data2m)==0)) != 0)
        {
          missing = T
          Missing_positions_indices[colSums(data2m)==0] = 1
        }
        if (case %in% c(1,2))
          data2m[data2m>0] = 1
        sum_data2m = sum(data2m)
        remove(data2m)
        
        if (assemblage_diversity)
        {
          logbeta = Result[[Ordered_realizations$ix[ifelse(case == 1 && taxon == "Platyhelminthes_18S",8,1)]]]@beta 
          # topic_compo : proportion of each MOTU in each topic (sums to 1 over MOTUs for each topic)
          topic_compo = exp(t(logbeta))
          topic_compo = topic_compo[,sort_normal_topic$ix]
          nb.OTUS[[case]][[i_taxon]] = vector(length = nb_topics, mode = "numeric")
          names(nb.OTUS[[case]][[i_taxon]]) = paste("Assemblage",1:nb_topics)
          for (k in 1:nb_topics)
            nb.OTUS[[case]][[i_taxon]][k] = length(which(topic_compo[,k] > 1/sum_data2m))
        }
        
        if (spatial_structure)
        {
          selected_spatial_distances = spatial_distances[Missing_positions_indices==0,Missing_positions_indices==0]
          diag(selected_spatial_distances) = NA
          max_dist = max(selected_spatial_distances,na.rm=T)
          min_dist = min(selected_spatial_distances,na.rm=T)
          weight_matrix_linear = (1 - selected_spatial_distances/max_dist)/(1 - min_dist/max_dist)
          diag(weight_matrix_linear) = 0
          weight_matrix_inverse = (max_dist/selected_spatial_distances - 1)/(max_dist/min_dist - 1)
          diag(weight_matrix_inverse) = 0
          weight_matrix_square = ((max_dist/selected_spatial_distances)^2 - 1)/((max_dist/min_dist)^2 - 1)
          diag(weight_matrix_square) = 0
          weight_matrix_cubic = ((max_dist/selected_spatial_distances)^3 - 1)/((max_dist/min_dist)^3 - 1)
          diag(weight_matrix_cubic) = 0
          
          I_linear = matrix(nrow =  ncol(documents), ncol = 3, dimnames = list(paste("Assemblage",1:nb_topics),c("I.observed","I.expected","I.p.value")), data = 0)
          I_inverse = matrix(nrow =  ncol(documents), ncol = 3, dimnames = list(paste("Assemblage",1:nb_topics),c("I.observed","I.expected","I.p.value")), data = 0)
          I_square = matrix(nrow =  ncol(documents), ncol = 3, dimnames = list(paste("Assemblage",1:nb_topics),c("I.observed","I.expected","I.p.value")), data = 0)
          I_cubic = matrix(nrow =  ncol(documents), ncol = 3, dimnames = list(paste("Assemblage",1:nb_topics),c("I.observed","I.expected","I.p.value")), data = 0)
          for (j in 1:ncol(documents))
          {
            I = ape::Moran.I(documents[,j],weight_matrix_linear)
            I_linear[j,1] = I$observed
            I_linear[j,2] = I$expected
            I_linear[j,3] = I$p.value
            
            I = ape::Moran.I(documents[,j],weight_matrix_inverse)
            I_inverse[j,1] = I$observed
            I_inverse[j,2] = I$expected
            I_inverse[j,3] = I$p.value
            
            I = ape::Moran.I(documents[,j],weight_matrix_square)
            I_square[j,1] = I$observed
            I_square[j,2] = I$expected
            I_square[j,3] = I$p.value
            
            I = ape::Moran.I(documents[,j],weight_matrix_cubic)
            I_cubic[j,1] = I$observed
            I_cubic[j,2] = I$expected
            I_cubic[j,3] = I$p.value
          }
          
          I_linear.observed_w.mean[[case]][[i_taxon]] = c(I_linear[,1],weighted.mean(I_linear[,1], w = colSums(documents)))
          I_inverse.observed_w.mean[[case]][[i_taxon]] = c(I_inverse[,1],weighted.mean(I_inverse[,1], w = colSums(documents)))
          I_square.observed_w.mean[[case]][[i_taxon]] = c(I_square[,1],weighted.mean(I_square[,1], w = colSums(documents)))
          I_cubic.observed_w.mean[[case]][[i_taxon]] = c(I_cubic[,1],weighted.mean(I_cubic[,1], w = colSums(documents)))
          I_linear.p.value_w.mean[[case]][[i_taxon]] = c(I_linear[,3],weighted.mean(I_linear[,3], w = colSums(documents)))
          I_inverse.p.value_w.mean[[case]][[i_taxon]] = c(I_inverse[,3],weighted.mean(I_inverse[,3], w = colSums(documents)))
          I_square.p.value_w.mean[[case]][[i_taxon]] = c(I_square[,3],weighted.mean(I_square[,3], w = colSums(documents)))
          I_cubic.p.value_w.mean[[case]][[i_taxon]] = c(I_cubic[,3],weighted.mean(I_cubic[,3], w = colSums(documents)))
          names(I_linear.observed_w.mean[[case]][[i_taxon]]) = names(I_inverse.observed_w.mean[[case]][[i_taxon]]) = 
            names(I_square.observed_w.mean[[case]][[i_taxon]]) = names(I_cubic.observed_w.mean[[case]][[i_taxon]]) =
            names(I_linear.p.value_w.mean[[case]][[i_taxon]]) = names(I_inverse.p.value_w.mean[[case]][[i_taxon]]) = 
            names(I_square.p.value_w.mean[[case]][[i_taxon]]) = names(I_cubic.p.value_w.mean[[case]][[i_taxon]]) = c(paste("Assemblage",1:nb_topics),"Weighted mean")
          # I_linear_mean[i_taxon,case] = mean(I_linear_taxon[[case]][,1])
          # I_inverse_mean[i_taxon,case] = mean(I_inverse_taxon[[case]][,1])
          # I_square_mean[i_taxon,case] = mean(I_square_taxon[[case]][,1])
        }
        
        if (environmental_data)
        {
          env_regression.axes[[case]][[i_taxon]] = matrix(nrow = ncol(lidar_data), ncol = 3, data = NA,
                                                          dimnames = list(colnames(lidar_data),
                                                                          # c("adj.R.squared","R.squared","p.value","spatial.p.value.adj","spatial.p.value")))
                                                                          c("adj.R.squared","R.squared","p.value")))
          env_regression.topics[[case]][[i_taxon]] = matrix(nrow = nb_topics, ncol = 3, data = NA,
                                                            dimnames = list(paste("Assemblage",1:nb_topics),
                                                                            c("adj.R.squared","R.squared","p.value")))
          env_selected.topics[[case]][[i_taxon]] = vector(length = nb_topics, mode = "list")
          names(env_selected.topics[[case]][[i_taxon]]) = paste("Assemblage",1:nb_topics)
          env_regression.topics.axes.adjR2[[case]][[i_taxon]] = matrix(nrow = nb_topics, ncol = ncol(lidar_data), data = NA,
                                                                 dimnames = list(paste("Assemblage",1:nb_topics),colnames(lidar_data)))
          env_regression.topics.axes.pval[[case]][[i_taxon]] = matrix(nrow = nb_topics, ncol = ncol(lidar_data), data = NA,
                                                                 dimnames = list(paste("Assemblage",1:nb_topics),colnames(lidar_data)))
          env_regression.topics.axes.sign[[case]][[i_taxon]] = matrix(nrow = nb_topics, ncol = ncol(lidar_data), data = NA,
                                                                      dimnames = list(paste("Assemblage",1:nb_topics),colnames(lidar_data)))
          
          RDA0 = rda(documents ~ 1, data = as.data.frame(scale(lidar_data[Missing_positions_indices==0,])), na.action = na.exclude)
          RDA_global_test = rda(documents ~ scale(lidar_data[Missing_positions_indices==0,]), na.action = na.exclude)
          # env_selected[[i_taxon]][[case]] = rep(F,ncol(env_descriptor))
          # Variable selection is only performed if the global RDA is significant; otherwise only the individually significant axes are kept
          if (length(anova(RDA_global_test)$'Pr(>F)'[1]) > 0 && !is.na(anova(RDA_global_test)$'Pr(>F)'[1]) && anova(RDA_global_test)$'Pr(>F)'[1] < 0.05)
          {
            RDA_selected = ordiR2step(RDA0, scope = as.formula(paste("documents ~", paste(colnames(lidar_data), collapse = " + "))), 
                                      # direction = "forward", trace = F)
                                      direction = "both", trace = F)
            # The axes are selected if they are either individually significant or selected by the variable selection procedure
            if (length(rownames(RDA_selected$CCA$biplot)) > 0 && !is.na(rownames(RDA_selected$CCA$biplot)))
            {
              env_selected[[case]][[i_taxon]] = colnames(lidar_data) %in% rownames(RDA_selected$CCA$biplot)
              env_regression[[case]][i_taxon,1] = as.numeric(unlist(RsquareAdj(RDA_selected)[2]))
              env_regression[[case]][i_taxon,2] = sum(RDA_selected$CCA$eig)/RDA_selected$tot.chi
              env_regression[[case]][i_taxon,3] = anova(RDA_selected)$'Pr(>F)'[1]
            } else
            {
              env_selected[[case]][[i_taxon]] = rep(F,ncol(lidar_data))
              env_regression[[case]][i_taxon,1:3] = rep(NA,3)
            }
            names(env_selected[[case]][[i_taxon]]) = colnames(lidar_data)
          } else
          {
            env_selected[[case]][[i_taxon]] = NA
            env_regression[[case]][i_taxon,1:3] = rep(NA,3)
          }
          # {
          #   significant_global_model[[case]][i_taxon] = F
          # }
          
          for (j in 1:ncol(documents))
          {
            RDA0 = rda(documents[,j] ~ 1, data = as.data.frame(scale(lidar_data[Missing_positions_indices==0,])), na.action = na.exclude)
            RDA_global_test = rda(documents[,j] ~ scale(lidar_data[Missing_positions_indices==0,]), na.action = na.exclude)
            # env_selected[[i_taxon]][[case]] = rep(F,ncol(env_descriptor))
            # Variable selection is only performed if the global RDA is significant; otherwise only the individually significant axes are kept
            if (length(anova(RDA_global_test)$'Pr(>F)'[1]) > 0 && !is.na(anova(RDA_global_test)$'Pr(>F)'[1]) && anova(RDA_global_test)$'Pr(>F)'[1] < 0.05)
            {
              RDA_selected = ordiR2step(RDA0, scope = as.formula(paste("documents[,j] ~", paste(colnames(lidar_data), collapse = " + "))), 
                                        # direction = "forward", trace = F)
                                        direction = "both", trace = F)
              # The axes are selected if they are either individually significant or selected by the variable selection procedure
              if (length(rownames(RDA_selected$CCA$biplot)) > 0 && !is.na(rownames(RDA_selected$CCA$biplot)))
              {
                env_selected.topics[[case]][[i_taxon]][[j]] = colnames(lidar_data) %in% rownames(RDA_selected$CCA$biplot)
                env_regression.topics[[case]][[i_taxon]][j,1] = as.numeric(unlist(RsquareAdj(RDA_selected)[2]))
                env_regression.topics[[case]][[i_taxon]][j,2] = sum(RDA_selected$CCA$eig)/RDA_selected$tot.chi
                env_regression.topics[[case]][[i_taxon]][j,3] = anova(RDA_selected)$'Pr(>F)'[1]
              } else
              {
                env_selected.topics[[case]][[i_taxon]][[j]] = rep(F,ncol(lidar_data))
                env_regression.topics[[case]][[i_taxon]][j,1:3] = rep(NA,3)
              }
              names(env_selected.topics[[case]][[i_taxon]][[j]]) = colnames(lidar_data)
            } else
            {
              env_selected.topics[[case]][[i_taxon]][[j]] = NA
              env_regression.topics[[case]][[i_taxon]][j,1:3] = rep(NA,3)
            }
            # {
            #   significant_global_model[[case]][i_taxon] = F
            # }
          }
          
          for (j_lidar in 1:ncol(lidar_data))
          {
            RDA = rda(scale(documents) ~ scale(lidar_data[Missing_positions_indices==0,j_lidar]), na.action = na.exclude)
            env_regression.axes[[case]][[i_taxon]][j_lidar,1] = as.numeric(unlist(RsquareAdj(RDA)[2]))
            env_regression.axes[[case]][[i_taxon]][j_lidar,2] = sum(RDA$CCA$eig)/RDA$tot.chi
            env_regression.axes[[case]][[i_taxon]][j_lidar,3] = anova(RDA)$'Pr(>F)'[1]
            for (j_doc in 1:ncol(documents))
            {
              RDA = rda(scale(documents[,j_doc]) ~ scale(lidar_data[Missing_positions_indices==0,j_lidar]), na.action = na.exclude)
              env_regression.topics.axes.adjR2[[case]][[i_taxon]][j_doc,j_lidar] = as.numeric(unlist(RsquareAdj(RDA)[2]))
              env_regression.topics.axes.sign[[case]][[i_taxon]][j_doc,j_lidar] = sign(cor(scale(documents[,j_doc]),scale(lidar_data[Missing_positions_indices==0,j_lidar])))
              env_regression.topics.axes.pval[[case]][[i_taxon]][j_doc,j_lidar] = anova(RDA)$'Pr(>F)'[1]
            }
          }
          
          if (!abiotic_pca)
          {
            env_regression.topics.axes.pval[[case]][[i_taxon]] = env_regression.topics.axes.pval[[case]][[i_taxon]][,c("Topography","Wetness","Slope","Slope.standard.deviation","Canopy.height")]
            env_regression.topics.axes.adjR2[[case]][[i_taxon]] = env_regression.topics.axes.adjR2[[case]][[i_taxon]][,c("Topography","Wetness","Slope","Slope.standard.deviation","Canopy.height")]
            env_regression.topics.axes.sign[[case]][[i_taxon]] = env_regression.topics.axes.sign[[case]][[i_taxon]][,c("Topography","Wetness","Slope","Slope.standard.deviation","Canopy.height")]
          }
          
          if (spatial_randomizations)
          {
            if (!missing)
            {
              # Spatial randomizations of documents:
              documents_randomized = vector(length=nb_abiotic_rndzations,mode="list")
              for (rndzation in 1:nb_abiotic_rndzations)
              {
                documents_randomized[[rndzation]] = matrix(nrow=nrow(documents),ncol=nb_topics,data=0)
                # Performing spatial randomizations
                for (k in 1:nb_topics)
                {
                  Spatial_documents = t(matrix(documents[,k],ncol=nrow_grid))
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
                  
                  documents_randomized[[rndzation]][,k] = as.vector(t(Spatial_documents_randomized))
                }
              }
            } else
            {
              ######
              Spatial_documents = vector(length=nb_topics,mode="list")
              for (k in 1:nb_topics)
              {
                dim1 = nrow_grid
                dim2 = ncol_grid
                Spatial_documents[[k]] = matrix(nrow=dim2,ncol=dim1,data=0)
                position_shift = 0
                for (j in 1:dim1)
                {
                  for (i in 1:dim2)
                  {
                    if (Missing_positions_indices[(j-1)*dim2+i]==0)
                      Spatial_documents[[k]][i,j] = documents[(j-1)*dim2+i-position_shift,k]    
                    else if (Missing_positions_indices[(j-1)*dim2+i]==1)
                    {
                      Spatial_documents[[k]][i,j] = NA
                      position_shift = position_shift+1
                    }
                  }
                }
                Spatial_documents[[k]] = t(Spatial_documents[[k]])
              }
              documents_randomized = vector(length=nb_abiotic_rndzations,mode="list")
              for (rndzation in 1:nb_abiotic_rndzations)
              {
                documents_randomized[[rndzation]] = matrix(nrow=nrow(documents),ncol=nb_topics,data=0)
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
                  
                  documents_randomized[[rndzation]][,k] = as.vector(t(Spatial_documents_randomized))[!is.na(as.vector(t(Spatial_documents_randomized)))]
                }
              }
              ######
            }
            
            nb_non_significant_rndzations = 0
            nb_non_significant_rndzations_adj = 0
            cor_abiotic_comparison_randomized = vector(length=nb_abiotic_rndzations,mode="numeric")
            for (rndzation in 1:nb_abiotic_rndzations) 
            {
              RDA = rda(scale(documents_randomized[[rndzation]]) ~ scale(lidar_data[Missing_positions_indices==0,]), na.action = na.exclude)
              if (env_regression[[case]][i_taxon,1] < as.numeric(unlist(RsquareAdj(RDA)[2])))
                nb_non_significant_rndzations_adj = nb_non_significant_rndzations_adj + 1
              if (env_regression[[case]][i_taxon,2] < sum(RDA$CCA$eig)/RDA$tot.chi)
                nb_non_significant_rndzations = nb_non_significant_rndzations + 1
            }
            if (env_regression[[case]][i_taxon,1] > 0)
              env_regression[[case]][i_taxon,4] = nb_non_significant_rndzations/nb_abiotic_rndzations
            else
              env_regression[[case]][i_taxon,4] = 1 - nb_non_significant_rndzations/nb_abiotic_rndzations 
            if (env_regression[[case]][i_taxon,2] > 0)
              env_regression[[case]][i_taxon,5] = nb_non_significant_rndzations/nb_abiotic_rndzations
            else
              env_regression[[case]][i_taxon,5] = 1 - nb_non_significant_rndzations/nb_abiotic_rndzations 
            
            for (j in 1:ncol(lidar_data))
            {
              for (rndzation in 1:nb_abiotic_rndzations) 
              {
                RDA = rda(scale(documents_randomized[[rndzation]]) ~ scale(lidar_data[Missing_positions_indices==0,j]), na.action = na.exclude)
                if (env_regression.detail[[case]][[i_taxon]][j,1] < as.numeric(unlist(RsquareAdj(RDA)[2])))
                  nb_non_significant_rndzations_adj = nb_non_significant_rndzations_adj + 1
                if (env_regression.detail[[case]][[i_taxon]][j,2] < sum(RDA$CCA$eig)/RDA$tot.chi)
                  nb_non_significant_rndzations = nb_non_significant_rndzations + 1
              }
              
              if (env_regression.detail[[case]][[i_taxon]][j,1] > 0)
                env_regression.detail[[case]][[i_taxon]][j,4] = nb_non_significant_rndzations/nb_abiotic_rndzations
              else
                env_regression.detail[[case]][[i_taxon]][j,4] = 1 - nb_non_significant_rndzations/nb_abiotic_rndzations 
              if (env_regression.detail[[case]][[i_taxon]][j,2] > 0)
                env_regression.detail[[case]][[i_taxon]][j,5] = nb_non_significant_rndzations/nb_abiotic_rndzations
              else
                env_regression.detail[[case]][[i_taxon]][j,5] = 1 - nb_non_significant_rndzations/nb_abiotic_rndzations 
            }
          }
        }
      } else
      {
        if (environmental_data)
        {
          cat("Missing result\n")
          env_regression[[case]][i_taxon,] = rep(NA,3)
          env_selected[[case]][[i_taxon]] = NA
          env_regression.axes[[case]][[i_taxon]] = NA
          env_regression.topics[[case]][[i_taxon]] = NA
          env_selected.topics[[case]][[i_taxon]] = NA
          env_regression.topics.axes.adjR2[[case]][[i_taxon]] = NA
          env_regression.topics.axes.pval[[case]][[i_taxon]] = NA
          env_regression.topics.axes.sign[[case]][[i_taxon]] = NA
        }
        if (spatial_structure)
        {
          I_linear.observed_w.mean[[case]][[i_taxon]] = I_inverse.observed_w.mean[[case]][[i_taxon]] = I_square.observed_w.mean[[case]][[i_taxon]] = I_cubic.observed_w.mean[[case]][[i_taxon]] = NA
          I_linear.p.value_w.mean[[case]][[i_taxon]] = I_inverse.p.value_w.mean[[case]][[i_taxon]] = I_square.p.value_w.mean[[case]][[i_taxon]] = I_cubic.p.value_w.mean[[case]][[i_taxon]] = NA
        }
        if (assemblage_diversity)
          nb.OTUS[[case]][[i_taxon]] = NA
      }
    }
    
    if (stability)
    {
      if (variableK)
      {
        nb_topics = optimalK[i_taxon]
      }
      
      K.folder_name = c(paste0("/Rtopicmodels_LDA_VEM_nb_topics",nb_topics,"_nb_real100_em_tol1e-07_var_tol1e-08_best_keep_occurrence/"),
                        paste0("/Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics",nb_topics,"_nb_iter1000_nb_real100_meanPosteriorDistributedLlh_thin25_burnin2000_occurrence/"),
                        paste0("/Rtopicmodels_LDA_VEM_nb_topics",nb_topics,"_nb_real100_em_tol1e-07_var_tol1e-08_best_keep/"),
                        paste0("/Rtopicmodels_LDA_Gibbs_alpha0.1_delta0.1_nb_topics",nb_topics,"_nb_iter1000_nb_real100_meanPosteriorDistributedLlh_thin25_burnin2000/"))
      
      K.folder_name_case = K.folder_name[case]
      data.folder_name_case = paste0(data.folder_name,K.folder_name_case)
      
      # stability_file = paste0(data.folder_name_case,"Stability_assessment_samplewise_maxmatching/Stability_samplewise.rds")
      stability_file = paste0(data.folder_name_case,"Stability_assessment_samplewise_maxmatching/Similarity_allRealPairs_samplewise.Rdata")
      stability_file_byTopic = paste0(data.folder_name_case,"Stability_assessment_samplewise_maxmatching/Similarity_allRealPairs_byTopic_samplewise.Rdata")
      if (file.exists(stability_file) && file.exists(stability_file_byTopic))
      {
        # stability_data.frame = readRDS(stability_file)
        Ordered_realizations = readRDS(paste0(data.folder_name_case,"Ordered_realizations.rds"))
        if (case == 1 && taxon == "Bacteries_16S")
          included_real = 1:44
        else if (case == 1 && taxon == "Champignons_18S")
          included_real = 1:86
        else if (case == 1 && taxon == "Platyhelminthes_18S")
          included_real = 8:100
        else
          included_real = 1:100
        
        load(stability_file)
        DKL100_allRealPairs = DKL100_allRealPairs[included_real,included_real]
        KL_allRealPairs_w_rndzations = KL_allRealPairs_w_rndzations[included_real,included_real]
        nES = t(DKL100_allRealPairs)[lower.tri(DKL100_allRealPairs,diag=F)]/KL_allRealPairs_w_rndzations[lower.tri(KL_allRealPairs_w_rndzations,diag=F)]
        # mean_sim[[case]][i_taxon] = stability_data.frame$`Normalized ES`[1]
        mean_sim[[case]][i_taxon] = mean(nES)
        # sim_intercept[[case]][i_taxon] = stability_data.frame$`Normalized ES=f(llh) intercept`[1]
        sim_intercept[[case]][i_taxon] = lm(nES[1:(length(included_real)-1)] ~ abs(Ordered_realizations$x[included_real[1]]-Ordered_realizations$x[included_real[-1]]))$coefficient[1]
        
        load(stability_file_byTopic)
        mean_sim.k[[case]][[i_taxon]] = vector(length = nb_topics, mode = "numeric")
        sim_intercept.k[[case]][[i_taxon]] = vector(length = nb_topics, mode = "numeric")
        names(mean_sim.k[[case]][[i_taxon]]) = paste("Assemblage",1:nb_topics)
        names(sim_intercept.k[[case]][[i_taxon]]) = paste("Assemblage",1:nb_topics)
        for (k in 1:nb_topics)
        {
            DKL100_allRealPairs_byTopic[[k]] = DKL100_allRealPairs_byTopic[[k]][included_real,included_real]
            KL_allRealPairs_randomized_byTopic[[k]] = KL_allRealPairs_randomized_byTopic[[k]][included_real,included_real]
            nES_topic = t(DKL100_allRealPairs_byTopic[[k]])[lower.tri(DKL100_allRealPairs_byTopic[[k]],diag=F)]/
                  t(KL_allRealPairs_randomized_byTopic[[k]])[lower.tri(KL_allRealPairs_randomized_byTopic[[k]],diag=F)]
            mean_sim.k[[case]][[i_taxon]][k] = mean(nES_topic)
            sim_intercept.k[[case]][[i_taxon]][k] = lm(nES_topic[1:(length(included_real)-1)] ~ abs(Ordered_realizations$x[included_real[1]]-Ordered_realizations$x[included_real[-1]]))$coefficient[1]
        }
        
        stability.plot.llh.diff[[case]][[i_taxon]] = ggplot(data = data.frame(x = abs(Ordered_realizations$x[included_real[1]]-Ordered_realizations$x[included_real[-1]]),
                                                                              y = nES[1:(length(included_real)-1)])) +
          ylim(0,1) +
          geom_smooth(aes(x,y),method='lm', size = 0.4, col="black", linetype = "dashed", se = F) +
          geom_point(aes(x,y), size = 0.4) +
          theme_bw() +
          ggtitle(paste(LETTERS[i_taxon],"-",taxo_names[i_taxon])) +
          theme(axis.title=element_text(size=10),
                plot.title=element_text(hjust=0, size=13),
                plot.margin=unit(c(1,1,1,0.5),"mm")) +
          labs(x=if (case %in% c(2,4)) "Llh difference to mean-llh real." else "Llh difference to best real.",
               y=if (case %in% c(2,4)) "Similarity to mean-llh real." else "Similarity to best real.")
        
        stability.plot.rank[[case]][[i_taxon]] = ggplot(data = data.frame(x = included_real[-1], y = nES[1:(length(included_real)-1)])) +
          ylim(0,1) +
          geom_smooth(aes(x,y),method='lm', size = 0.4, col="black", linetype = "dashed", se = F) +
          geom_point(aes(x,y), size = 0.4) +
          theme_bw() +
          ggtitle(paste(LETTERS[i_taxon],"-",taxo_names[i_taxon])) +
          theme(axis.title=element_text(size=10),
                plot.title=element_text(hjust=0, size=13),
                plot.margin=unit(c(1,1,1,0.5),"mm")) +
          labs(x=if (case %in% c(2,4)) "Real. sorted by decreasing\n difference to mean-llh" else "Real. sorted by decreasing llh", 
               y=if (case %in% c(2,4)) "Similarity to mean-llh real." else "Similarity to best real.")
      } else
      {
        cat("No stability file found:",taxon,"\n")
        mean_sim[[case]][i_taxon] = NA
        mean_sim.k[[case]][[i_taxon]] = NA
        sim_intercept[[case]][i_taxon] = NA
        sim_intercept.k[[case]][[i_taxon]] = NA
        stability.plot.llh.diff[[case]][[i_taxon]] = NA
        stability.plot.rank[[case]][[i_taxon]] = NA
      }
      if (case %in% c(1,3))
      {
        for (r in 1:100) 
        {
          if (r == 1)
            alpha = read.table(file = paste0(data.folder_name_case,"1st_best_realization/estimated_alpha.txt"))
          else 
            alpha = read.table(file = paste0(data.folder_name_case,if (r == 2) "2nd" else if (r == 3) "3rd" else paste0(r,"th"),"_best_realization_maxmatching/estimated_alpha.txt"))
          alpha_all_real[[case]][r,i_taxon] = unlist(alpha)
        }
      }
    }
  }
  
  if (optimalK_comput)
  {
    if (optimalK_mpar_perplexity && !optimalK_mpar_aic)
    {
      ggsave(filename = paste0(figure_folder,spline_df,"dfSplined10sampleFoldPerplexity_errbar_vs_K",if (test_data) "_testData" else "","_",topic_range,"t",inference_insert,occ_insert,".pdf"), 
             do.call("arrangeGrob", c(plot.spline.perplexity[[case]][!is.na(plot.spline.perplexity[[case]])], nrow = ifelse(test_data,1,3), ncol = ifelse(test_data,2,3))),
             height = ifelse(test_data,5/2,5/2*3), width = ifelse(test_data,5,5/2*3))
      
      # ggsave(filename = paste0(figure_folder,spline_df,"dfSplined10sampleFoldPerplexity_vs_K",if (test_data) "_testData" else "","_",topic_range,"t",inference_insert,occ_insert,".pdf"), 
      #        do.call("arrangeGrob", c(plot.spline.perplexity.mean.median[[case]][!is.na(plot.spline.perplexity.mean.median[[case]])], nrow = ifelse(test_data,1,3), ncol = ifelse(test_data,2,3))),
      #        height = ifelse(test_data,5/2,5/2*3), width = ifelse(test_data,5,5/2*3))
      
      # spl = split(plot.spline.perplexity.corrCumul[[case]][!is.na(plot.spline.perplexity.corrCumul[[case]])],
      #             (seq_along(plot.spline.perplexity.corrCumul[[case]][!is.na(plot.spline.perplexity.corrCumul[[case]])])-1) %/% 20)
      # ppl = lapply(spl, function(g) marrangeGrob(grobs = g, nrow = 2, ncol = 4, layout_matrix = matrix(data=1:7,nrow=2,byrow=T)))
      # pdf(paste0(figure_folder,spline_df,"dfSplined10sampleFoldPerplexity_errbar_vs_K_corrCumul_",topic_range,"t",case_insert,".pdf"),height = 5, width = 10)
      # print(ppl)
      # dev.off()
      # 
      # spl = split(plot.spline.perplexity.derivative[[case]][!is.na(plot.spline.perplexity.derivative[[case]])],
      #             (seq_along(plot.spline.perplexity.derivative[[case]][!is.na(plot.spline.perplexity.derivative[[case]])])-1) %/% 20)
      # ppl = lapply(spl, function(g) marrangeGrob(grobs = g, nrow = 2, ncol = 4, layout_matrix = matrix(data=1:7,nrow=2,byrow=T)))
      # pdf(paste0(figure_folder,spline_df,"dfSplined10sampleFoldPerplexity_derivative_vs_K_",topic_range,"t",case_insert,".pdf"),height = 5, width = 10)
      # print(ppl)
      # dev.off()
      # 
      # spl = split(plot.spline.perplexity.derivative2[[case]][!is.na(plot.spline.perplexity.derivative2[[case]])],
      #             (seq_along(plot.spline.perplexity.derivative2[[case]][!is.na(plot.spline.perplexity.derivative2[[case]])])-1) %/% 20)
      # ppl = lapply(spl, function(g) marrangeGrob(grobs = g, nrow = 2, ncol = 4, layout_matrix = matrix(data=1:7,nrow=2,byrow=T)))
      # pdf(paste0(figure_folder,spline_df,"dfSplined10sampleFoldPerplexity_derivative2_vs_K_",topic_range,"t",case_insert,"pdf"),height = 5, width = 10)
      # print(ppl)
      # dev.off()
    }
    
    if (optimalK_mpar_aic && !optimalK_mpar_perplexity)
    {
      ggsave(filename = paste0(figure_folder,spline_df,"dfSplined",nb_real,"realAIC_errbar_vs_K",if (test_data) "_testData" else "","_",topic_range,"t",inference_insert,occ_insert,".pdf"), 
             do.call("arrangeGrob", c(plot.spline.aic[[case]][!is.na(plot.spline.aic[[case]])], nrow = ifelse(test_data,1,3), ncol = ifelse(test_data,2,3))),
             height = ifelse(test_data,5/2,5/2*3), width = ifelse(test_data,5,5/2*3))
      
      # ggsave(filename = paste0(figure_folder,spline_df,"dfSplined50realAIC_vs_K_testData_",topic_range,"t",inference_insert,occ_insert,".pdf"), 
      #        do.call("arrangeGrob", c(plot.spline.perplexity.mean.median[[case]][!is.na(plot.spline.perplexity.mean.median[[case]])], nrow = 2, ncol = 4)),
      #        height = 5, width = 10)
    }
    
    if (optimalK_mpar_perplexity && optimalK_mpar_aic)
    {
      #With bacteria:
      ggsave(filename = paste0(figure_folder,spline_df,"dfSplined",nb_real,"realAIC-10sampleFoldPerplexity_errbar_vs_K",if (test_data) "_testData" else "","_",topic_range,"t",inference_insert,occ_insert,".pdf"),
             do.call("arrangeGrob", c(plot.spline.both[[case]][!is.na(plot.spline.both[[case]])], nrow = ifelse(test_data,1,4), ncol = ifelse(test_data,2,4))),
             height = ifelse(test_data,5/2,5/2*4), width = ifelse(test_data,5,5/2*4))
      
      #No bacteria:
      # ggsave(filename = paste0(figure_folder,spline_df,"dfSplined",nb_real,"realAIC-10sampleFoldPerplexity_errbar_vs_K",if (test_data) "_testData" else "_noBacteria","_",topic_range,"t",inference_insert,occ_insert,".pdf"),
      #        do.call("arrangeGrob", c(plot.spline.both[[case]][!is.na(plot.spline.both[[case]])][-c(1,2)], nrow = ifelse(test_data,1,3), ncol = ifelse(test_data,2,4))),
      #        height = ifelse(test_data,5/2,5/2*3), width = ifelse(test_data,5,5/2*4))
      
      #Bacteria only for AIC:
      # ggsave(filename = paste0(figure_folder,spline_df,"dfSplined",nb_real,"realAIC-10sampleFoldPerplexity_errbar_vs_K",if (test_data) "_testData" else "_AIC.Bacteria.only","_",topic_range,"t",inference_insert,occ_insert,".pdf"),
      #        do.call("arrangeGrob", c(plot.spline.both[[case]][!is.na(plot.spline.both[[case]])], nrow = ifelse(test_data,1,4), ncol = ifelse(test_data,2,4))),
      #        height = ifelse(test_data,5/2,5/2*4), width = ifelse(test_data,5,5/2*4))
    }
  }
  
  if (environmental_data)
  {
    saveRDS(list(env_regression[[case]],env_selected[[case]]),file=paste0(figure_folder,"Lidar",pca_insert,"_RDA",variableK_insert,"_all.variables",inference_insert,occ_insert,".rds"))
    saveRDS(env_regression.axes[[case]],file=paste0(figure_folder,"Lidar",pca_insert,"_RDA",variableK_insert,"_axis.by.axis",inference_insert,occ_insert,".rds"))
    saveRDS(list(env_regression.topics[[case]],env_selected.topics[[case]]),file=paste0(figure_folder,"Lidar",pca_insert,"_RDA",variableK_insert,"_topic.by.topic",inference_insert,occ_insert,".rds"))
  }
  
  if (stability)
  {
    ggsave(filename = paste0(figure_folder,"/Stability_taxoGroups_llh.diff_",if (case %in% c(2,4)) "Gibbs" else "VEM","100r",if (test_data) "_testData" else "",variableK_insert,inference_insert,occ_insert,".pdf"), 
           do.call("arrangeGrob", c(stability.plot.llh.diff[[case]][!is.na(stability.plot.llh.diff[[case]])], nrow = ifelse(test_data,1,3), ncol = ifelse(test_data,2,3))),
           height = ifelse(test_data,5/2,5/2*3), width = ifelse(test_data,5,5/2*3))
    
    ggsave(filename = paste0(figure_folder,"/Stability_taxoGroups_rank_",if (case %in% c(2,4)) "Gibbs" else "VEM","100r",if (test_data) "_testData" else "",variableK_insert,inference_insert,occ_insert,".pdf"), 
           do.call("arrangeGrob", c(stability.plot.rank[[case]][!is.na(stability.plot.rank[[case]])], nrow = ifelse(test_data,1,3), ncol = ifelse(test_data,2,3))),
           height = ifelse(test_data,5/2,5/2*3), width = ifelse(test_data,5,5/2*3))
  }
}

if (test_data && test_error)
{
  library(maxmatching)
  library(ggplot2)
  source("/Users/guilhemsommeria-klein/Desktop/Code/R/LDA_project/LDA_topic_correspondence_fun.R")
  
  data.folder_name = "/Users/guilhemsommeria-klein/Desktop/Serveur_Bioclust.data/Donnees_PetitPlateau/Test_data/"
  figure_folder = "/Users/guilhemsommeria-klein/Desktop/Manuscrits/ThesisPaper2_Topic_Modelling_PP/Figures/Submission2_analyses/"
  
  test_taxo_vect = c("Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics0.001_1000sampledreads",
                    "Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics0.005_1000sampledreads",
                    "Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics0.01_1000sampledreads",
                    "Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics0.02_1000sampledreads",
                    "Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics0.05_1000sampledreads",
                    "Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics0.1_1000sampledreads",
                    "Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics0.5_1000sampledreads",
                    "Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics1_1000sampledreads")
  
  Rsquared = list()
  mean_nb_OTUs = vector(length = length(test_taxo_vect), mode = "numeric")
  for (case in 1:2)
  {
    Rsquared[[case]] = vector(length = length(test_taxo_vect), mode = "numeric")
    names(Rsquared[[case]]) = c("0.001","0.005","0.01", "0.02", "0.05","0.1","0.5","1")
    if (case == 1)
      result.folder = "/Rtopicmodels_LDA_VEM_nb_topics5_nb_real100_em_tol1e-06_var_tol1e-08_best_keep_occurrence"
    else 
      result.folder = "/Rtopicmodels_LDA_VEM_nb_topics5_nb_real100_em_tol1e-06_var_tol1e-08_best_keep"
    for (taxon in test_taxo_vect)
    {
      taxon.folder_name = paste0(data.folder_name,taxon)
      load(paste0(taxon.folder_name,"/true_documents.Rdata"))
      
      load(paste0(taxon.folder_name,result.folder,result.folder,".Rdata"))
      Ordered_realizations = readRDS(paste0(taxon.folder_name,result.folder,"/Ordered_realizations.rds"))
      documents = Result[[Ordered_realizations$ix[1]]]@gamma
      
      cor_documents = 1-cor(documents,true_norm_documents)
      Topic_correspondence = LDA_topic_correspondence_fun(cor_documents,maxmatching=1,greedymatching=0)
      Rsquared[[case]][which(test_taxo_vect == taxon)] = 1 - sum((documents-true_norm_documents[,Topic_correspondence])^2)/
                                                              sum(sweep(true_norm_documents[,Topic_correspondence],2,colMeans(true_norm_documents[,Topic_correspondence]),'-')^2)
      
      if (case == 1)
      {
        load(paste0(taxon.folder_name,"/data2m_testdata.Rdata"))
        nb_OTUs = vector(length = ncol(data2m), mode = "numeric")
        for (i in 1:ncol(data2m))
          nb_OTUs[i] = length(which(data2m[,i]>0))
        mean_nb_OTUs[which(test_taxo_vect == taxon)] = mean(nb_OTUs)
      }
      
      # logbeta = Result[[Ordered_realizations$ix[1]]]@beta 
      # # topic_compo : proportion of each MOTU in each topic (sums to 1 over MOTUs for each topic)
      # topic_compo = exp(t(logbeta))
      # topic_compo = topic_compo[,sort_normal_topic$ix]
      # nb.OTUS[[case]][[i_taxon]] = vector(length = nb_topics, mode = "numeric")
      # names(nb.OTUS[[case]][[i_taxon]]) = paste("Assemblage",1:nb_topics)
      # for (k in 1:nb_topics)
      #   nb.OTUS[[case]][[i_taxon]][k] = length(which(topic_compo[,k] > 1/sum_data2m))
    }
  }
  
  tmp.plot = ggplot(data = data.frame(y = c(Rsquared[[2]],Rsquared[[1]]), x = rep(mean_nb_OTUs,2),
                                      shape = c(rep("Read count",length(test_taxo_vect)),rep("Occurrence",length(test_taxo_vect))))) +
    geom_point(aes(x,y,shape=shape), col = "black", size = 2) +
    geom_line(data = data.frame(y = Rsquared[[1]], x = mean_nb_OTUs), aes(x,y), linetype = "dashed", size = 0.5) +
    geom_line(data = data.frame(y = Rsquared[[2]], x = mean_nb_OTUs), aes(x,y), linetype = "dashed", size = 0.5) +
    # geom_point(aes(x=mean.sample.size,y=occ.Rsquared,fill=occ.fill), col = "black", size = 2, shape = 15) +
    # geom_line(aes(x=mean.sample.size,y=occ.Rsquared), linetype = "dashed", size = 0.5) +
    theme_bw() +
    # ggtitle(paste(LETTERS[i_taxon],"-",taxo_names[i_taxon])) +
    theme(axis.title = element_text(size=20),
          axis.text = element_text(size=15),
          plot.title=element_text(hjust=0, size=13),
          legend.position = c(0.8, 0.2),
          # legend.position = "right",
          legend.text=element_text(size=15), 
          legend.title=element_blank(),
          legend.background = element_rect(colour = "white"),
          legend.box.background = element_rect(colour = "black"),
          plot.margin=unit(c(1,1,1,0.5),"mm")) +
    labs(x="Mean nb. of distinct OTUs per sample\n (out of 1,000 simulated OTUs)", 
    # labs(x="Mean nb. of OTUs per sample\n wrt. the total nb. of simulated OTUs", 
         y="Model's R squared")
  
  pdf(paste0(figure_folder,"Test_data_Rsquared_ab-occ.pdf"))
  print(tmp.plot)
  dev.off()
}

if (environmental_data && maps)
{
  if (abiotic_pca)
  {
    pdf(paste0(figure_folder,"Explained_variance_lidar_PCA_",nrow(lidar_pca$co),"variables.pdf"))
    par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
    plot(lidar_pca$eig/sum(lidar_pca$eig),ylab="Proportion of variance",xlab="PCA axis")
    abline(h = 1/nrow(lidar_pca$co), lty = 2)
    dev.off()
    
    pdf(paste0(figure_folder,"Explained_variance_lidar_PCA_",nrow(lidar_pca$co),"variables_barplot.pdf"))
    par(cex.lab=1.5,cex.main=1.7,cex.axis=2,lwd=2)
    barplot(lidar_pca$eig,ylab="",xlab="",col="black")
    abline(h = 1, lty = 2)
    dev.off()
    
    pdf(paste0(figure_folder,"lidar_PCA_axes_",nrow(lidar_pca$co),"variables_noweights.pdf"))
    #bottom left top right
    #par(mar=c(5.1,4.1,4.1,2.1)
    par(mar=c(7.1,4.1,4.1,2.1))
    par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
    for (j in 1:ncol(lidar_pca$li))
    { 
      variable_sort = sort.int(lidar_pca$c1[,j]^2,decreasing=T,index.return=T)
      # variable_weights contains the squared normed loadings of the current axis on all input variables sorted by decreasing order,
      # the signed normed loadings of the current axis in the same order, and the names of all input variables: 
      variable_weights = data.frame(variable_sort$x,lidar_pca$c1[variable_sort$ix,j],rownames(lidar_pca$co)[variable_sort$ix])
      variable_color = vector(length=nrow(lidar_pca$co),mode="character")
      variable_color[which(variable_weights[,2]>0)]="black"
      variable_color[which(variable_weights[,2]<0)]="red"
      
      # The "eigenvalue" of each axis is given as the actual eigenvalue divided by the trace, times the number of input variables
      x = plot(variable_weights[,1],type="p",ann=T,xaxt="n",ylab = "Environmental variable proportion in axis", xlab="",col=variable_color,
               main = paste("Axis #",j," - Eigenvalue = ",format(lidar_pca$eig[j]/sum(lidar_pca$eig)*nrow(lidar_pca$co),digits=4)))
      axis(1, at=1:nrow(lidar_pca$co), labels = F)
      labels = variable_weights[,3]
      # text(1:ncol(abiotic_data), par("usr")[3], adj = 0, srt = 45, pos=1, labels = labels, xpd = T, cex=1.2)
      text(1:nrow(lidar_pca$co), par("usr")[3], srt = 45, pos=1, offset=3, labels = labels, xpd = T, cex=1)
      #title=paste("Axis #",element_index," ",format(chemi_PCA$eig[element_index]/ncol(data_chemi)*100,digits=2),"% variance",sep="")
    }
    dev.off()
  }
  
  color.pal = colorRampPalette(c("#00007F", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),space = "Lab")
  grain = 10
  lidar_data1 = lidar_data
  # for (j in 1:ncol(lidar_data))
  #   lidar_data1[,j]  = as.vector(t(matrix(data=lidar_data[,j],nrow=nrow_grid)))
  coord = expand.grid(x = grain*seq(1,ncol_grid,1), y = grain*seq(1,nrow_grid,1))
  # coord = expand.grid(y = grain*seq(1,nrow_grid,1), x = grain*seq(1,ncol_grid,1))[,c(2,1)]
  coord0 = coord
  newdata = expand.grid(x = seq(min(coord$x) - grain, max(coord$x) + grain, 1), y = seq(min(coord$y) - grain, max(coord$y) + grain, 1))
  # newdata = expand.grid(y = seq(min(coord$y) - grain, max(coord$y) + grain, 1), x = seq(min(coord$x) - grain, max(coord$x) + grain, 1))[,c(2,1)]
  tmp.plot = list()
  spatial_topicmix_kriged = list()
  for (j in 1:ncol(lidar_data))
  {
    # Building a kriged spatial map of the grid, topic by topic
    tmp = data.frame(x=coord$x,y=coord$y,z=lidar_data1[,j])
    mod = gstat(id = "z", formula = z~1, locations = ~x+y, data=tmp, model=vgm(10, "Exp", 20))
    spatial_topicmix_kriged[[j]] = predict(mod, newdata, na.action=na.omit)
    
    tmp.plot[[j]] = ggplot(data=spatial_topicmix_kriged[[j]]) +
      geom_raster(aes(x,y,fill=z.pred)) +
      # qplot(x, y, data=spatial_topicmix_kriged, geom="raster", fill=z.pred) +
      scale_fill_gradientn(colours=color.pal(7)) +
      coord_equal() + theme_minimal() +
      labs(fill = colnames(lidar_data)) + 
      ggtitle(if (abiotic_pca) paste("Axis #",j," - Eigenvalue = ",format(lidar_pca$eig[j]/sum(lidar_pca$eig)*nrow(lidar_pca$co),digits=4)) else colnames(lidar_data)[j]) +
      theme(legend.position="bottom", legend.text=element_text(size=7), 
            legend.title=element_text(size=8), axis.title=element_blank(), 
            axis.text = element_blank(),
            plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
      guides(fill = guide_colorbar(barwidth = 8, barheight = 0.4, title.position="bottom")) +
      scale_y_continuous(limits = c(min(coord$y) - floor(grain/2), max(coord$y) + floor(grain/2)), expand = c(0,0)) +
      scale_x_continuous(limits = c(min(coord$x) - floor(grain/2), max(coord$x) + floor(grain/2)), expand = c(0,0)) + 
      geom_point(data = data.frame(coord0,z.pred=rep(0,nrow(coord0))), aes(x,y), color="black", size=1.5, alpha=0.3)
  }
  pdf(paste0(figure_folder,"lidar",pca_insert,"_spatial_maps.pdf"))
  for (j in 1:ncol(lidar_data))
    print(tmp.plot[[j]])
  dev.off()
}

if (maps && test_plots)
{
  ###############
  # Topic by topic:
  color.pal = colorRampPalette(c("#00007F", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),space = "Lab")
  grain = 2
  nb_topics = 3
  data.folder_name = paste0("/Users/guilhemsommeria-klein/Desktop/Serveur_Bioclust.data/Donnees_PetitPlateau/Platyhelminthes_18S")
  data.folder_name1 = paste0("/Rtopicmodels_LDA_VEM_nb_topics",nb_topics,"_nb_real100_em_tol1e-07_var_tol1e-08_best_keep/")
  
  # data.file_name = paste0("Rtopicmodels_LDA_VEM_nb_topics",nb_topics,"_nb_real100_em_tol1e-07_var_tol1e-08_best_keep.Rdata")
  # load(paste0(data.folder_name,data.folder_name1,data.file_name))
  # Ordered_realizations = readRDS(paste0(data.folder_name_case,"Ordered_realizations.rds"))
  # # One topic per columns; sums to one over each row/site:
  # documents = Result[[Ordered_realizations$ix[1]]]@gamma
  
  spatial_topicmix_kriged = readRDS(paste0(data.folder_name,data.folder_name1,"1st_best_realization/Spatial_topicmix_kriged.rds"))
  
  coord = expand.grid(x = grain*seq(1,ncol_grid,1), y = grain*seq(1,nrow_grid,1))
  # coord = expand.grid(y = grain*seq(1,nrow_grid,1), x = grain*seq(1,ncol_grid,1))[,c(2,1)]
  coord0 = coord
  newdata = expand.grid(x = seq(min(coord$x) - grain, max(coord$x) + grain, 1), y = seq(min(coord$y) - grain, max(coord$y) + grain, 1))
  # newdata = expand.grid(y = seq(min(coord$y) - grain, max(coord$y) + grain, 1), x = seq(min(coord$x) - grain, max(coord$x) + grain, 1))[,c(2,1)]
  
  tmp.plot = list()
  for (j in 1:ncol(documents))
  {
    # Building a kriged spatial map of the grid, topic by topic
    # tmp = data.frame(x=coord$x,y=coord$y,z=documents[,j])
    # mod = gstat(id = "z", formula = z~1, locations = ~x+y, data=tmp, model=vgm(10, "Exp", 20))
    # spatial_topicmix_kriged = predict(mod, newdata, na.action=na.omit)

    tmp.plot[[j]] = ggplot(data=spatial_topicmix_kriged[[j]]) +
      geom_raster(aes(x,y,fill=z.pred)) +
      # qplot(x, y, data=spatial_topicmix_kriged, geom="raster", fill=z.pred) +
      scale_fill_gradientn(colours=color.pal(7)) +
      coord_equal() + theme_minimal() +
      labs(fill = colnames(documents)) +
      ggtitle(colnames(documents)[j]) +
      theme(legend.position="bottom", legend.text=element_text(size=7),
            legend.title=element_text(size=8), axis.title=element_blank(),
            axis.text = element_blank(),
            plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
      guides(fill = guide_colorbar(barwidth = 8, barheight = 0.4, title.position="bottom")) +
      scale_y_continuous(limits = c(min(coord$y) - floor(grain/2), max(coord$y) + floor(grain/2)), expand = c(0,0)) +
      scale_x_continuous(limits = c(min(coord$x) - floor(grain/2), max(coord$x) + floor(grain/2)), expand = c(0,0)) +
      geom_point(data = data.frame(coord0,z.pred=rep(0,nrow(coord0))), aes(x,y), color="black", size=0.2, alpha=0.3)
  }
  pdf(paste0(figure_folder,"Platy_3t100r_ab_map.pdf"))
  for (j in 1:ncol(documents))
    print(tmp.plot[[j]])
  dev.off()

  #######################
  # 15t prevalence plots:
  pdf(paste0(figure_folder,"/15t_prevalence_barplot_allGroups.pdf"),height=7*0.9*3/2)
  par(mfrow = c(3,3))
  #bottom left top right
  par(mar=c(1.5,4.5,4.5,2.4))
  ii_taxon = 0
  for (taxon in taxo_vect)
  {
    ii_taxon = ii_taxon+1
    data.folder_name = paste0("/Users/guilhemsommeria-klein/Desktop/Serveur_Bioclust.data/Donnees_PetitPlateau/",taxon,"/Rtopicmodels_LDA_VEM_nb_topics15_nb_real10_em_tol1e-07_var_tol1e-08_best_keep/")
    result.file = paste0(data.folder_name,"Rtopicmodels_LDA_VEM_nb_topics15_nb_real10_em_tol1e-07_var_tol1e-08_best_keep.Rdata")
    load(result.file)
    Ordered_realizations = readRDS(paste0(data.folder_name,"Ordered_realizations.rds"))
    documents = Result[[Ordered_realizations$ix[1]]]@gamma
    par(cex.lab=1.7,cex.main=1.7,cex.axis=1.4,lwd=2)
    barplot(sort(apply(documents,2,mean),decreasing=T),space=0.5,ylab=ifelse(ii_taxon == 1,"Assemblage prevalence",""),xlab=ifelse(ii_taxon == 1,"Assemblages",""),border=NA,col="grey")
    title(paste(LETTERS[which(taxo_vect == taxon)],"-",taxo_names[which(taxo_vect == taxon)]))
    # abline(h = 1, lty = 2)
  }
  dev.off()
  
  ###############
  # One plot:
  grain = 10
  nb_topics = 3
  data.folder_name = paste0("/Users/guilhemsommeria-klein/Desktop/Serveur_Bioclust.data/Donnees_PetitPlateau/Champignons_18S")
  data.folder_name1 = paste0("/Rtopicmodels_LDA_VEM_nb_topics",nb_topics,"_nb_real100_em_tol1e-07_var_tol1e-08_best_keep_occurrence/")
  
  # data.file_name = paste0("Rtopicmodels_LDA_VEM_nb_topics",nb_topics,"_nb_real100_em_tol1e-07_var_tol1e-08_best_keep.Rdata")
  # load(paste0(data.folder_name,data.folder_name1,data.file_name))
  # Ordered_realizations = readRDS(paste0(data.folder_name_case,"Ordered_realizations.rds"))
  # # One topic per columns; sums to one over each row/site:
  # documents = Result[[Ordered_realizations$ix[1]]]@gamma
  
  spatial_topicmix_kriged = readRDS(paste0(data.folder_name,data.folder_name1,"1st_best_realization/Spatial_topicmix_kriged.rds"))
  # spatial_topicmix_kriged = readRDS(paste0(data.folder_name,data.folder_name1,"8th_best_realization_maxmatching/Spatial_topicmix_kriged.rds"))
  
  load(paste0(data.folder_name,"/Missing_positions_indices.Rdata"))
  load(paste0(data.folder_name,"/data2m_filled.Rdata"))
  colSums_data2m = colSums(data2m)
  spatial_colSums_data2m = matrix(nrow=nrow_grid, ncol=ncol_grid, data=0)
  position_shift = 0
  for (i in 1:nrow_grid)
  {
    for (j in 1:ncol_grid)
    {
      missing_index = (i-1)*ncol_grid+j
      if (Missing_positions_indices[missing_index]==0)
      {
        spatial_colSums_data2m[i,j] = colSums_data2m[missing_index-position_shift]    
      } else if (Missing_positions_indices[missing_index]==1)
      {
        spatial_colSums_data2m[i,j] = NA
        position_shift = position_shift+1
      }
    }
  }
  spatial_colSums_data2m = as.vector(t(spatial_colSums_data2m))
  Missing_positions_indices0 = Missing_positions_indices
  Missing_positions_indices0[which(spatial_colSums_data2m==0)] = 1
  coord = expand.grid(x = grain*seq(1,ncol_grid,1), y = grain*seq(1,nrow_grid,1))
  coord0 = coord[-which(Missing_positions_indices0==1),]
  
  # Removing interpolated relative abundances falling outside of [0,1]
  # The sum of all assemblages may then not be 1 for all pixels
  for (k in 1:nb_topics)
  {
    spatial_topicmix_kriged[[k]]$z.pred[spatial_topicmix_kriged[[k]]$z.pred<0] = 0
    spatial_topicmix_kriged[[k]]$z.pred[spatial_topicmix_kriged[[k]]$z.pred>1] = 1
  }
  # Storing all assemblages into a single dataframe:
  spatial_topicmix_kriged_all_topics = data.frame(x=spatial_topicmix_kriged[[1]]$x,y=spatial_topicmix_kriged[[1]]$y)
  for (k in 1:nb_topics)
  {
    spatial_topicmix_kriged_all_topics = cbind(spatial_topicmix_kriged_all_topics,spatial_topicmix_kriged[[k]]$z.pred)
    colnames(spatial_topicmix_kriged_all_topics)[2+k] = paste0("z.pred",k)
  }
  # Defining the color in each location based on the assemblage composition:
  spatial_topicmix_kriged_all_topics_colors = data.frame(spatial_topicmix_kriged_all_topics[,1:2],red=0,green=0,blue=0)
  coord_matrix = t(as.matrix(spatial_topicmix_kriged_all_topics[,3:(2+nb_topics)]))
  # col_matrix = matrix(data=c(c(0,0,255),c(0,255,0),c(255,0,0)),ncol=3)
  col_matrix = matrix(data=c(c(0,0,255),c(255,0,0),c(0,255,0)),ncol=3)
  rownames(col_matrix) = c("red","green","blue")
  spatial_topicmix_kriged_all_topics_colors[,3:5] = t(col_matrix %*% coord_matrix)
    
  # Plotting each topic as a distinct color on a single map
  tmp.plot = ggplot(data = spatial_topicmix_kriged_all_topics_colors) +
    geom_raster(data = spatial_topicmix_kriged_all_topics_colors, aes(x, y), 
                fill=rgb(spatial_topicmix_kriged_all_topics_colors[,3:5]/max(255,max(spatial_topicmix_kriged_all_topics_colors[,3:5]))), inherit.aes = F) +
    #geom_raster() +
    coord_equal() +
    theme_minimal() +
    #labs(fill=paste0("Assemblage ",k)) + ggtitle(letters[k]) +
    theme(axis.title=element_blank(), axis.text = element_blank(),
          #                 legend.position="bottom", legend.text=element_text(size=13), 
          #                 legend.title=element_text(size=20),
          # margin: top right bottom left
          plot.title=element_text(hjust=0,size=11), plot.margin=unit(c(0.1,0.1,0.1,0.1),"mm")) +
    #plot.title=element_text(hjust=0,size=20), plot.margin=unit(c(2,-10,0,-10),"mm")) +
    #plot.title=element_text(hjust=0,size=20)) +
    #guides(fill = guide_colorbar(barwidth = 8, barheight = 0.5, title.position="bottom")) +
    scale_y_continuous(limits=c(5,395), expand = c(0,0)) +
    scale_x_continuous(limits=c(5,295), expand = c(0,0)) + 
    geom_point(data = data.frame(coord0,z.pred=rep(0,nrow(coord0))), aes(x,y), color="black", size=1.5, alpha=0.3, inherit.aes = F)
  
  pdf(paste0(figure_folder,"Fungi_3t100r_occ_map.pdf"))
  print(tmp.plot)
  dev.off()
}

