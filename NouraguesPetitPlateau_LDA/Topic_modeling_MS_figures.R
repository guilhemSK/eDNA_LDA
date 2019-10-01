corrplot = 0
kriged_topics = 1
lidar_maps = 0
stability = 0
size_diversity = 0
taxo_tree = 0
testdata = 0
rank_abundance = 0
mock_stability = 0

figure_folder = "/Users/guilhemsommeria-klein/Desktop/Manuscrits/ThesisPaper2_Topic_Modelling_PP/Figures/Submission2_analyses"

if (corrplot)
{
  library(corrplot)
  library(reshape2)
  library(ggplot2)
  
  ##########
  tableS1 = read.table(paste0(figure_folder,"/tableS1_corrplot.csv"), sep=";", colClasses="vector", row.names = 1, header = T)
  tableS1 = apply(tableS1,2,as.numeric)
  rownames(tableS1) = c("Wetness",  "Slope", "Slope std. dev.", "Canopy height", "Light", "Tree density", "Number of\n tree deaths", "Loss of\n above-ground biomass")
  colnames(tableS1) = c("Topography", "Wetness",  "Slope", "Slope std. dev.", "Canopy height", "Light", "Tree density", "Number of\n tree deaths")
  tableS1[is.na(tableS1)] = 0
  
  col2 = colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                                      "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                                      "#4393C3", "#2166AC", "#053061"))
  
  pdf(paste0(figure_folder,"/tableS1_corrplot.pdf"), height = 8, width = 11)
  corrplot(tableS1,
           method = "color",
           col = rev(col2(101)),
           cl.lim=c(-max(abs(tableS1)),max(abs(tableS1))), 
           is.corr=F,
           # Hor. labels, black labels, no grid:
           tl.cex = 1.4, tl.srt = 30, tl.col = "black", tl.offset = 1, addgrid.col = "grey",
           cl.cex = 1.3, cl.ratio = 0.2,
           mar = c(0.1, 0.1, 0.1, 0.1))
  dev.off()
  #########
  
  tableS2 = read.table(paste0(figure_folder,"/tableS2_corrplot.csv"), sep=";", colClasses="vector", row.names = 1, header = T)
  rownames_tableS2 = rownames(tableS2)
  tableS2 = apply(tableS2,2,as.numeric)
  rownames(tableS2) = rownames_tableS2
  tableS2[is.na(tableS2)] = 0
  tableS2.color = tableS2
  tableS2.color[1,] = 0
  
  plotS2.heatmap = ggplot(data = melt(tableS2.color[nrow(tableS2.color):1,])) +
    geom_tile(aes(x = Var2, y = Var1, fill = value), color = "grey") +
    geom_text(data = melt(tableS2[nrow(tableS2):1,]),
              aes(x = Var2, y = Var1, label = value)) +
    scale_x_discrete(position = "top", labels = c("Axis 1","Axis 2","Axis 3","Axis 4","Axis 5")) +
    scale_y_discrete(labels = rev(c("Eigenvalue",
                                    "Topography",
                                    "Wetness",
                                    "Slope",
                                    "Slope std. dev.",
                                    "Canopy height",
                                    "Light",
                                    "Tree density",
                                    "Number of tree deaths",
                                    "Loss of above-ground biomass"))) +
    # scale_fill_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = 0, 
    # scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B" , midpoint = 0, 
    # limits = c(-max(abs(table2)),max(abs(table2))), 
    # na.value = "white" , name = expression("Adj. R"^2*" x sign(cor.)")) +
    # scale_fill_gradientn(colours = rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
    #                                      "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
    #                                      "#4393C3", "#2166AC", "#053061")),
    scale_fill_gradientn(colours = rev(c("#B2182B", "#D6604D", "#F4A582",
                                         "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                         "#4393C3", "#2166AC")),
                         limits = c(-max(abs(tableS2.color)),max(abs(tableS2.color))),
                         #values = rescale(1:11,c(-max(abs(table2)),max(abs(table2)))), 
                         name = expression("Squared loading\n x sign(loading)")) +
    guides(fill = guide_colorbar(barwidth = 0.8, barheight = 12, title.position="top")) +
    #coord_equal() +
    theme_void() +
    theme(axis.title=element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0, size=12),
          axis.text.y = element_text(size = 12, hjust = 0),
          legend.text=element_text(size=10),
          legend.title=element_text(size=12),
          plot.margin=unit(c(0.1,0.1,0.1,0.1),"mm"))
  pdf(paste0(figure_folder,"/tableS2_heatmap.pdf"), height = 5.5)
  print(plotS2.heatmap)
  dev.off()
  ##########
  
  table2 = read.table(paste0(figure_folder,"/table2_corrplot.csv"), sep=";", colClasses="vector", row.names = 1, header = T)
  rownames_table2 = rownames(table2)
  table2.color = apply(table2,2,as.numeric)
  table2 = apply(table2,2,as.character)
  rownames(table2.color) = rownames_table2
  rownames(table2) = rownames_table2
  table2.color[is.na(table2.color)] = 0
  
  # col1 = colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
  #                           "cyan", "#007FFF", "blue","#00007F"))
  # col2 = colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
  #                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
  #                           "#4393C3", "#2166AC", "#053061"))
  # col3 = colorRampPalette(c("#67001F","white", "#053061"))
  # 
  # pdf(paste0(figure_folder,"/table2_corrplot.pdf"))
  # corrplot(table2,
  #          col = c("red","white","blue"),
  #          # Hor. labels, black labels, no grid:
  #          tl.srt = 30, tl.col = "black", addgrid.col = NA, tl.offset = 2,
  #          # not a corr. matrix, range:
  #          is.corr=F, cl.lim = range(table2), cl.cex = 0.7, cl.ratio = 0.2,
  #          # shape the matrix into a square:
  #          win.asp = ncol(table2)/nrow(table2))
  # dev.off()
  
  plot2.heatmap = ggplot(data = melt(table2.color[nrow(table2.color):1,])) +
    geom_tile(aes(x = Var2, y = Var1, fill = value), color = "grey") +
    geom_text(data = melt(table2[nrow(table2):1,]),
              aes(x = Var2, y = Var1, label = value), size = 2.7) +
    scale_x_discrete(position = "top", labels = c("All selected axes","Axis 1","Axis 2","Axis 3","Axis 4","Axis 5")) +
    scale_y_discrete(labels = rev(c(expression(paste("Bacteria 16S - ",italic("all assemb."))),
                                expression(italic("      Blue assemb. - Terra firme")),
                                expression(italic("      Green assemb. - Hydromorphic")),
                                expression(italic("      Red assemb. - Exposed rock")),
                                expression(paste("Protists 18S - ",italic("all assemb."))),
                                expression(italic("      Blue assemb. - Terra firme")),
                                expression(italic("      Green assemb. - Hydromorphic")),
                                expression(italic("      Red assemb. - Exposed rock")),
                                expression(paste("Fungi 18S - ",italic("all assemb."))),
                                expression(italic("      Blue assemb. - Terra firme")),
                                expression(italic("      Green assemb. - Hydromorphic")),
                                expression(italic("      Red assemb. - Exposed rock")),
                                expression(paste("Arthropods 18S - ",italic("all assemb."))),
                                expression(italic("      Blue assemb.")),
                                expression(italic("      Green assemb.")),
                                expression(italic("      Red assemb. - Exposed rock")),
                                expression(paste("Nematodes 18S - ",italic("all assemb."))),
                                expression(italic("      Blue assemb.")),
                                expression(italic("      Green assemb.")),
                                expression(italic("      Red assemb.")),
                                expression(paste("Platyhelm. 18S - ",italic("all assemb."))),
                                expression(italic("      Blue assemb.")),
                                expression(italic("      Green assemb.")),
                                expression(italic("      Red assemb.")),
                                expression(paste("Annelids 18S - ",italic("all assemb."))),
                                expression(italic("      Blue assemb.")),
                                expression(italic("      Green assemb.")),
                                expression(italic("      Red assemb."))))) +
    # scale_fill_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = 0, 
    # scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B" , midpoint = 0, 
                         # limits = c(-max(abs(table2)),max(abs(table2))), 
                         # na.value = "white" , name = expression("Adj. R"^2*" x sign(cor.)")) +
    # scale_fill_gradientn(colours = rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
    #                                  "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
    #                                  "#4393C3", "#2166AC", "#053061")),
    scale_fill_gradientn(colours = rev(c("#B2182B", "#D6604D", "#F4A582",
                                         "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                         "#4393C3", "#2166AC")),
                         limits = c(-max(abs(table2.color)),max(abs(table2.color))),
                         #values = rescale(1:11,c(-max(abs(table2)),max(abs(table2)))), 
                         name = expression("Adj. R"^2*" x sign(cor.)")) +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 12, title.position="top")) +
    #coord_equal() +
    theme_void() +
    theme(axis.title=element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 0.05, hjust = 0, size=10),
          axis.text.y = element_text(size = 10, hjust = 0),
          legend.text=element_text(size=10),
          legend.title=element_text(size=10),
          plot.margin=unit(c(0.1,0.1,0.1,0.1),"mm"))
  pdf(paste0(figure_folder,"/table2_heatmap.pdf"))
  print(plot2.heatmap)
  dev.off()
  ###############
  
  tableS3 = read.table(paste0(figure_folder,"/tableS3_corrplot.csv"), sep=";", colClasses="vector", header = T)
  tableS3.color = apply(tableS3[,-1],2,as.numeric)
  tableS3 = apply(tableS3[,-1],2,as.character)
  rownames(tableS3.color) = rownames(tableS3) = c(letters,"A","B")
  tableS3.color[is.na(tableS3.color)] = 0
  
  plot.S3.heatmap = ggplot(data = melt(tableS3.color[nrow(tableS3.color):1,])) +
    geom_tile(aes(x = Var2, y = Var1, fill = value), color = "grey") +
    geom_text(data = melt(tableS3[nrow(tableS3):1,]),
              aes(x = Var2, y = Var1, label = value), size = 2.7) +
    scale_x_discrete(position = "top", labels = c("All selected axes","Axis 1","Axis 2","Axis 3","Axis 4","Axis 5")) +
    scale_y_discrete(labels = rev(c(expression(paste("Bacteria 16S - ",italic("all assemb."))),
                                    expression(italic("      Blue assemb. - Terra firme")),
                                    expression(italic("      Green assemb. - Hydromorphic")),
                                    expression(italic("      Red assemb. - Exposed rock")),
                                    expression(paste("Protists 18S - ",italic("all assemb."))),
                                    expression(italic("      Blue assemb. - Terra firme")),
                                    expression(italic("      Green assemb.")),
                                    expression(italic("      Red assemb. - Exposed rock")),
                                    expression(paste("Fungi 18S - ",italic("all assemb."))),
                                    expression(italic("      Blue assemb.")),
                                    expression(italic("      Green assemb.")),
                                    expression(italic("      Red assemb.")),
                                    expression(paste("Arthropods 18S - ",italic("all assemb."))),
                                    expression(italic("      Blue assemb.")),
                                    expression(italic("      Green assemb.")),
                                    expression(italic("      Red assemb.")),
                                    expression(paste("Nematodes 18S - ",italic("all assemb."))),
                                    expression(italic("      Blue assemb.")),
                                    expression(italic("      Green assemb.")),
                                    expression(italic("      Red assemb.")),
                                    expression(paste("Platyhelm. 18S - ",italic("all assemb."))),
                                    expression(italic("      Blue assemb.")),
                                    expression(italic("      Green assemb.")),
                                    expression(italic("      Red assemb.")),
                                    expression(paste("Annelids 18S - ",italic("all assemb."))),
                                    expression(italic("      Blue assemb.")),
                                    expression(italic("      Green assemb.")),
                                    expression(italic("      Red assemb."))))) +
    # scale_fill_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = 0, 
    #                      limits = c(-max(abs(tableS3)),max(abs(tableS3))), 
    #                      na.value = "white" , name = expression("Adj. R"^2*" x sign(cor.)")) +
    # scale_fill_gradientn(colours = rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
    #                                      "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
    #                                      "#4393C3", "#2166AC", "#053061")),
    scale_fill_gradientn(colours = rev(c("#B2182B", "#D6604D", "#F4A582",
                                         "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                         "#4393C3", "#2166AC")),
                         limits = c(-max(abs(tableS3.color)),max(abs(tableS3.color))),
                         #values = rescale(1:11,c(-max(abs(table2)),max(abs(table2)))), 
                         name = expression("Adj. R"^2*" x sign(cor.)")) +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 12, title.position="top")) +
    #coord_equal() +
    theme_void() +
    theme(axis.title=element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 0.05, hjust = 0, size=10),
          axis.text.y = element_text(size = 10, hjust = 0),
          legend.text=element_text(size=10),
          legend.title=element_text(size=10),
          plot.margin=unit(c(0.1,0.1,0.1,0.1),"mm"))
  pdf(paste0(figure_folder,"/tableS3_heatmap.pdf"))
  print(plot.S3.heatmap)
  dev.off()
  ###############
}
if (kriged_topics)
{
  library(raster)
  library(ggplot2)
  library(gridExtra)
  library(scales)
  library(gstat)
  
  horizontal = 0
  vertical = 1
  
  one_plot = 1
  
  discrete = 0
  
  raw_lidar = 0
  pca_axes = 1
  
  legend_bar = 1
  
  barcode_insert_list = rep(c("Bacteries_16S","Protistes_18S","Champignons_18S"),2)
  barcode_labels_list = rep(c("Bacteria","Protists","Fungi"),2)
  
  nb_occupied_sites_threshold_vect = rep(1,6)
  no_rare_vect = rep(0,6)
  
  # Number of sites an OTU must ocuppy to be kept
  # (if nb_occupied_sites_threshold = 1, all OTUs with non-zero abundance are kept)
  
  # Occ. vs. ab. K=3:
  ###################
  nb_topics_vect = rep(3,6)
  inverted_red_green = c(F,T,F,F,F,T)
  occurrence_vect = c(0,0,0,1,1,1)
  color.pal = colorRampPalette(c("#00007F", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), space = "Lab")
  
  # nb_topics_vect = c(5,2,4)
  
  # Optimal K vs. K = 3 occ.:
  ######################
  # nb_topics_vect = c(3,3,3,5,2,4)
  # inverted_red_green = c(F,F,T,F,F,F)
  # occurrence_vect = c(1,1,1,1,1,1)
  # color.pal = colorRampPalette(c("blue","cyan","green","#FF7F00","red"))
  
  nb_real = 100
  em_tol = 10^-7
  var_tol = 10^-8
  grain = 10
  nrow_grid = 39
  ncol_grid = 29
  
  if (vertical)
  {
    horizVert_insert = "_vertical"
  } else if (horizontal)
    horizVert_insert = "_horizontal"
  
  if (discrete)
  {
    discrete_insert = "_discrete"
  } else
    discrete_insert = ""
  
  if (raw_lidar)
  {
    pca_insert = "_lidarVariables"
  } else if (pca_axes)
  {
    pca_insert = "_environmentalAxes"
  }
  
  if (legend_bar)
    legend_bar_insert = "_withLegend"
  else
    legend_bar_insert = ""
  
  data_insert = "Donnees_PetitPlateau"  
  
  ########
  if (raw_lidar)
  {
    setwd(paste0("/Users/guilhemsommeria-klein/Desktop/These/",data_insert,"/Lidar/"))
    r_topo = raster("r_topol_0.asc")
    r_topo_transposed_flipped = flip(t(r_topo),'y')
    r_slope = raster("r_slopel_0.asc")
    r_slope_transposed_flipped = flip(t(r_slope),'y')
    r_wetness = raster("r_wetnessl_0.asc")
    r_wetness_transposed_flipped = flip(t(r_wetness),'y')
    
    #TopoPP = as.data.frame(r_topo, row.names=NULL, optional=FALSE, xy=TRUE, na.rm=FALSE, long=FALSE)
    TopoPP = as.data.frame(r_topo_transposed_flipped, row.names=NULL, optional=FALSE, xy=TRUE, na.rm=FALSE, long=FALSE)
    SlopePP = as.data.frame(r_slope_transposed_flipped, row.names=NULL, optional=FALSE, xy=TRUE, na.rm=FALSE, long=FALSE)
    SlopePPmod = SlopePP 
    SlopePPmod$layer[SlopePP$layer>0.50] = 0.50
    #     TopoPP_inverse = TopoPP
    #     TopoPP_inverse$layer = max(TopoPP_inverse$layer)-TopoPP_inverse$layer
    WetnessPP = as.data.frame(r_wetness_transposed_flipped, row.names=NULL, optional=FALSE, xy=TRUE, na.rm=FALSE, long=FALSE)
    data_lidar_3topics = list()
    data_lidar_3topics[[1]] = TopoPP
    #data_lidar_3topics[[2]] = TopoPP_inverse
    data_lidar_3topics[[2]] = WetnessPP
    #data_lidar_3topics[[3]] = SlopePP
    data_lidar_3topics[[3]] = SlopePPmod
    
    #Assemblage_labels_list = c("Terra firme","Hydromorphic","Exposed rock")
    
    #map_labels = c("Topography","Inverse topography","Slope")
    map_labels = c("Topography (m a.s.l.)","Topographic Wetness Index","Slope")
  } else if (pca_axes)
  {
    lidar_data_path = "/Users/guilhemsommeria-klein/Desktop/These/Donnees_PetitPlateau/Lidar/"
    lidar_data_file = "Lidar_locsites_E20_origin_mean_10_extended.csv"
    
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
    
    library(ade4)
    lidar_pca = dudi.pca(as.data.frame(lidar_data), row.w = rep(1, nrow(lidar_data))/nrow(lidar_data), 
                         col.w = rep(1, ncol(lidar_data)), center = TRUE, scale = TRUE, 
                         scannf = F, nf = ncol(lidar_data))
    lidar_data = lidar_pca$li[,1:5]
    
    map_labels = paste("Env. axis",1:5)
  }
  ###########
  
  coord = expand.grid(x = grain*seq(1,ncol_grid,1), y = grain*seq(1,nrow_grid,1))
  
  tmp.plot = list()
  lidar.plot = list()
  #   chemi.plot = list()
  if (pca_axes)
  {
    newdata = expand.grid(x = seq(min(coord$x) - grain, max(coord$x) + grain, 1), y = seq(min(coord$y) - grain, max(coord$y) + grain, 1))
    #spatial_lidarAxes_kriged = list()
  }
  barcode_index = 0
  for (barcode_insert in barcode_insert_list)
  {
    barcode_index = barcode_index+1
    
    nb_topics = nb_topics_vect[barcode_index]
    
    if (occurrence_vect[barcode_index])
      occurrence_insert = "_occurrence"
    else
      occurrence_insert = ""
    
    if (nb_occupied_sites_threshold_vect[barcode_index] > 1)
    {
      remove_single_sites_insert = paste0("_woOTUsWithLessThan",nb_occupied_sites_threshold_vect[barcode_index],"Sites")
    } else 
    {
      remove_single_sites_insert = ""
    }
    
    if (no_rare_vect[barcode_index])
    {
      no_rare_insert = "_noRareOTU"
    } else 
    {
      no_rare_insert = ""
    }
    
    filename_insert = "Rtopicmodels_LDA_VEM"
    # local_prefix = "/Users/guilhemsommeria-klein/Desktop/These/" 
    local_prefix = "/Users/guilhemsommeria-klein/Desktop/Serveur_Bioclust.data/" 
    local_dirname = paste0(local_prefix,data_insert,"/",barcode_insert,"/")
    local_subdirname = paste0(local_dirname,filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best_keep",occurrence_insert,remove_single_sites_insert,no_rare_insert)
    #     if (occurrence_vect[barcode_index] && (nb_occupied_sites_threshold_vect[barcode_index] == 1) && !no_rare_vect[barcode_index] && (barcode_insert == "Platyhelminthes_18S"))
    #       subsubdirname = paste(local_subdirname,"ordered_realization_number_5/",sep="")
    #     else 
    # subsubdirname = paste(local_subdirname,"/ordered_realization_number_1",sep="") 
    subsubdirname = paste0(local_subdirname,"/1st_best_realization") 
    # subsubdirname = paste(local_subdirname,"/ordered_realization_number_1",sep="")  
    # subsubsubdirname = paste(subsubdirname,"/topics_site_repartition_info/",sep="")  
    
    load(paste0(local_dirname,"Missing_positions_indices.Rdata"))
    load(paste0(local_dirname,"data2m_filled.Rdata"))
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
    coord0 = coord[-which(Missing_positions_indices0==1),]
    
    # setwd(subsubsubdirname)
    spatial_topicmix_kriged = readRDS(paste0(subsubdirname,"/Spatial_topicmix_kriged.rds"))
    
    if (vertical)
      tmp.plot[[barcode_index]] = list()
    
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
    if (nb_topics == 3)
    {
      if (inverted_red_green[barcode_index])
      {
        col_matrix = matrix(data=c(c(0,0,255),c(255,0,0),c(0,255,0)),ncol=3)
      } else
        col_matrix = matrix(data=c(c(0,0,255),c(0,255,0),c(255,0,0)),ncol=3)
    } else
    {
      col_matrix = col2rgb(color.pal(nb_topics))
    }
    rownames(col_matrix) = c("red","green","blue")
    spatial_topicmix_kriged_all_topics_colors[,3:5] = t(col_matrix %*% coord_matrix)

    ##################
    
    # Plotting each topic as a distinct color on a single map
    current_plot = ggplot(data = spatial_topicmix_kriged_all_topics_colors) +
      geom_raster(data = spatial_topicmix_kriged_all_topics_colors, aes(x, y), 
                  fill=rgb(spatial_topicmix_kriged_all_topics_colors[,3:5]/max(255,max(spatial_topicmix_kriged_all_topics_colors[,3:5]))), inherit.aes = F) +
      #geom_raster() +
      coord_equal() +
      theme_minimal() +
      labs(fill=barcode_labels_list[barcode_index]) + 
      ggtitle(LETTERS[barcode_index]) +
      #labs(fill=paste0("Assemblage ",k)) + ggtitle(letters[k]) +
      theme(axis.title=element_blank(),
            axis.text = element_blank(),
            #                 legend.position="bottom", legend.text=element_text(size=13), 
            #                 legend.title=element_text(size=20),
            # margin: top right bottom left
            plot.title=element_text(hjust=0,size=14),
            plot.margin=unit(c(1.1,0.1,0.1,0.1),"mm")) +
      #plot.title=element_text(hjust=0,size=20), plot.margin=unit(c(2,-10,0,-10),"mm")) +
      #plot.title=element_text(hjust=0,size=20)) +
      #guides(fill = guide_colorbar(barwidth = 8, barheight = 0.5, title.position="bottom")) +
      scale_y_continuous(limits=c(5,395), expand = c(0,0)) +
      scale_x_continuous(limits=c(5,295), expand = c(0,0)) + 
      geom_point(data = data.frame(coord0,z.pred=rep(0,nrow(coord0))), aes(x,y), color="black", size=0.6, alpha=0.3, inherit.aes = F)
    
    if (horizontal)
    {
      current_plot = current_plot + ggtitle(letters[(k-1)*(length(barcode_insert_list)+1) + barcode_index])
      tmp.one.plot[[barcode_index]] = current_plot
    }
    else if (vertical)
    {
      #current_plot = current_plot + ggtitle(letters[(barcode_index-1)*nb_topics+k])
      tmp.plot[[barcode_index]] = current_plot
    }
    
    if (vertical && barcode_index<4)
    {
      if (raw_lidar)
      {
        lidar.plot[[barcode_index]] = qplot(x, y, data = data_lidar_3topics[[barcode_index]], geom="raster", fill=layer) +
          #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", function(x) round(x,0))) + 
          #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", scientific_format())) +
          scale_fill_gradientn(colours=color.pal(7)) +
          coord_equal() +
          labs(fill=map_labels[barcode_index]) + 
          #ggtitle(letters[k+length(barcode_insert_list)*nb_topics]) +
          ggtitle(LETTERS[length(barcode_insert_list)+barcode_index]) +
          scale_y_continuous(limits=c(5,395), expand = c(0,0)) +
          scale_x_continuous(limits=c(5,295), expand = c(0,0)) + 
          theme_minimal() +
          theme(legend.position="bottom", 
                legend.text=element_text(size=11), 
                legend.title=element_text(size=14), 
                axis.title=element_blank(), 
                #           theme(legend.position="bottom", legend.text=element_text(size=11), 
                #                 legend.title=element_text(size=12), axis.title=element_blank(), 
                axis.text = element_blank(),
                #                 plot.title=element_text(hjust=0,size=16), plot.margin=unit(c(0,1,-2,2),"mm")) +
                plot.title=element_text(hjust=0,size=14), 
                plot.margin=unit(c(0.1,0.1,0.1,0.1),"mm")) +
          #                 plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
          guides(fill = guide_colorbar(barwidth = 8, barheight = 0.5, title.position="bottom"))
        # size = 14 and 15 for vertical 5 topics
      } else if (pca_axes)
      {
        # Building a kriged spatial map of the grid, topic by topic
        # tmp = data.frame(x=coord$x,y=coord$y,z=lidar_data[,barcode_index])
        # mod = gstat(id = "z", formula = z~1, locations = ~x+y, data=tmp, model=vgm(10, "Exp", 20))
        # spatial_lidarAxes_kriged[[barcode_index]] = predict(mod, newdata, na.action=na.omit)
        
        plotted_spatial_lidarAxes_kriged = spatial_lidarAxes_kriged[[barcode_index]]
        plotted_spatial_lidarAxes_kriged$z.pred[plotted_spatial_lidarAxes_kriged$z.pred > 2.5] = 2.5
        plotted_spatial_lidarAxes_kriged$z.pred[plotted_spatial_lidarAxes_kriged$z.pred < -2.5] = -2.5
        
        lidar.plot[[barcode_index]] = ggplot(data=plotted_spatial_lidarAxes_kriged) +
          geom_raster(aes(x,y,fill=z.pred)) +
          # qplot(x, y, data=spatial_topicmix_kriged, geom="raster", fill=z.pred) +
          #scale_fill_gradientn(colours=color.pal(7), guide=FALSE) +
          scale_fill_gradientn(colours=color.pal(7), guide = ifelse(legend_bar,T,F)) +
          coord_equal() +
          labs(fill = map_labels[barcode_index]) + 
          ggtitle(LETTERS[length(barcode_insert_list)+barcode_index]) +
          theme_minimal() +
          theme(axis.title=element_blank(), 
                axis.text = element_blank(),
                legend.position="bottom",
                legend.text=element_text(size=11),
                legend.title=element_text(size=14),
                plot.title=element_text(hjust=0,size=14), 
                plot.margin=unit(c(1.1,0.1,0.1,0.1),"mm")) +
          # guides(fill = guide_colorbar(barwidth = 8, barheight = 0.5, title.position="bottom")) +
          scale_y_continuous(limits = c(min(coord$y) - floor(grain/2), max(coord$y) + floor(grain/2)), expand = c(0,0)) +
          scale_x_continuous(limits = c(min(coord$x) - floor(grain/2), max(coord$x) + floor(grain/2)), expand = c(0,0)) +
          geom_point(data = data.frame(coord,z.pred=rep(0,nrow(coord))), aes(x,y), color="black", size=0.6, alpha=0.3)
        if (legend_bar)
          lidar.plot[[barcode_index]] = lidar.plot[[barcode_index]] + guides(fill = guide_colorbar(barwidth = 8, barheight = 0.4, title.position="top"))
      }
    }
    # End of the barcode_insert loop:
  }
  
  figure_folder = "/Users/guilhemsommeria-klein/Desktop/Manuscrits/ThesisPaper2_Topic_Modelling_PP/Figures/Submission2_analyses"
  if (vertical)
  {
    if (one_plot)
    {
      # Figure 4:
      #Old version:
      #grid.arrange(grobs = c(tmp.plot, lidar.plot),ncol=3,heights = c(10*4/3*0.72,10*4/3*0.72,10*4/3), width = 10)
      # ggsave(filename = paste0(figure_folder,"/Kriged_topics_with_Lidar_bact-prot-fungi18S_occ-ab",horizVert_insert,"_oneplot",discrete_insert,".pdf"),
      #        do.call("arrangeGrob",c(tmp.plot,lidar.plot,ncol=3,nrow=3)),width=5/2*3,height=5/2*3*4/3)
      
      # Occ. vs ab. figure, new version:
      plot = grid.arrange(grobs = c(tmp.plot,lidar.plot), heights = c(1,1,ifelse(raw_lidar || pca_axes && legend_bar,1.26,1)), nrow = 3, layout_matrix = matrix(1:9,nrow=3,byrow=T))
      ggsave(filename = paste0(figure_folder,"/Kriged_topics_with_Lidar_bact-prot-fungi18S_occ-ab",horizVert_insert,"_oneplot",pca_insert,discrete_insert,legend_bar_insert,".pdf"),
             plot=plot,width=ifelse(raw_lidar || pca_axes && legend_bar,7.5,8.5),height=11)
      
      # Optimal K figure:
      # plot = grid.arrange(grobs = tmp.plot, nrow = 2, layout_matrix = matrix(1:6,nrow=2,byrow=T))
      # ggsave(filename = paste0(figure_folder,"/Kriged_topics_optimalK_bact-prot-fungi18S_occ-ab",horizVert_insert,"_oneplot",pca_insert,discrete_insert,".pdf"),
      #        plot=plot,width=7.5,height=6.748)
      
      # Fig S3:
      #grid.arrange(grobs = tmp.plot, ncol=3, width = 10, height = 10/3*4*4/3, layout_matrix = rbind(c(1,2,3),c(4,5,6),c(7,8,NA),c(9,10,NA)))
      #grid.arrange(grobs = as.list(tmp.plot), ncol=4)
      #     ggsave(filename = paste0("Kriged_topics_arth-nema-plath-anne_occ-ab",horizVert_insert,"_oneplot",discrete_insert,".pdf"),
      #            do.call("arrangeGrob",c(as.list(tmp.plot),ncol=4)),width=4)
      #dev.off()
    } else
    {
      ggsave(filename = paste0(figure_folder,"/Kriged_topics_bacteria.pdf"), do.call("arrangeGrob", 
                                                              c(tmp.plot[[1]], ncol=nb_topics)),
             height = 10/3*4/3*(length(barcode_insert_list)), width = 10)
      #     ggsave(filename = "Kriged_topics_bacteria_labels.pdf", do.call("arrangeGrob", 
      #             c(tmp.plot[[1]], ncol=nb_topics)),
      #            height = 10/3*4/3*length(barcode_insert_list), width = 10)
      #     ggsave(filename = "Kriged_topics_arthropods.pdf", do.call("arrangeGrob", 
      #           c(tmp.plot[[1]], ncol=nb_topics)),
      #            height = 10/3*4/3*length(barcode_insert_list), width = 10)
    }
  } else if (horizontal)
    #     ggsave(filename = paste0("Kriged_topics_with_Lidar_",length(barcode_insert_list),"barcodes_",horizVert_insert,"_SlopeMod.pdf"), do.call("arrangeGrob", 
    #                                                               c(tmp.plot[[1]], tmp.plot[[2]], tmp.plot[[3]], lidar.plot, nrow=nb_topics)),
    #            height = 10*4/3, width = 10*(length(barcode_insert_list)+1)/3)
    ggsave(filename = paste0(figure_folder,"/Kriged_topics_with_Lidar_",length(barcode_insert_list),"barcodes",occurrence_insert,remove_single_sites_insert,no_rare_insert,horizVert_insert,".pdf"),
           do.call("arrangeGrob",c(tmp.plot[[1]], tmp.plot[[2]], tmp.plot[[3]], nrow=nb_topics)),
           height = 10*4/3, width = 10*length(barcode_insert_list)/3)
  
}
####################################################################################################
####################################################################################################
if (lidar_maps)
{
  library(ggplot2)
  library(raster)
  library(gridExtra)
  library(scales)
  library(gstat)
  
  lidar_data_path = "/Users/guilhemsommeria-klein/Desktop/These/Donnees_PetitPlateau/Lidar/"
  
  #########
  # Raw Lidara data (not avergaed per cell - only 6 are available):
 
  # r_topo = raster(paste0(lidar_data_path,"r_topol_0.asc"))
  # r_topo_transposed_flipped = flip(t(r_topo),'y')
  # r_slope = raster(paste0(lidar_data_path,"r_slopel_0.asc"))
  # r_slope_transposed_flipped = flip(t(r_slope),'y')
  # r_wetness = raster(paste0(lidar_data_path,"r_wetnessl_0.asc"))
  # r_wetness_transposed_flipped = flip(t(r_wetness),'y')
  # 
  # #TopoPP = as.data.frame(r_topo, row.names=NULL, optional=FALSE, xy=TRUE, na.rm=FALSE, long=FALSE)
  # TopoPP = as.data.frame(r_topo_transposed_flipped, row.names=NULL, optional=FALSE, xy=TRUE, na.rm=FALSE, long=FALSE)
  # SlopePP = as.data.frame(r_slope_transposed_flipped, row.names=NULL, optional=FALSE, xy=TRUE, na.rm=FALSE, long=FALSE)
  # SlopePPmod = SlopePP
  # SlopePPmod$layer[SlopePP$layer>0.50] = 0.50
  # #     TopoPP_inverse = TopoPP
  # #     TopoPP_inverse$layer = max(TopoPP_inverse$layer)-TopoPP_inverse$layer
  # WetnessPP = as.data.frame(r_wetness_transposed_flipped, row.names=NULL, optional=FALSE, xy=TRUE, na.rm=FALSE, long=FALSE)
  # data_lidar_3topics = list()
  # data_lidar_3topics[[1]] = TopoPP
  # #data_lidar_3topics[[2]] = TopoPP_inverse
  # data_lidar_3topics[[2]] = WetnessPP
  # #data_lidar_3topics[[3]] = SlopePP
  # data_lidar_3topics[[3]] = SlopePPmod
  
  ##########
  # Lidar data averaged per cell (those used to compute the correlations):
  
  nrow_grid = 39
  ncol_grid = 29
  
  lidar_data_file = "Lidar_locsites_E20_origin_mean_10_extended.csv"
  # Lidar data
  lidar_data = read.table(paste0(lidar_data_path,lidar_data_file), sep=";", colClasses="vector", row.names = 1, header = T)
  # Removing from data_lidar: canop08, which is already contained in dif_can, tree_density2008, AGBloss_mort2008, ndeath2008: 
  lidar_data = apply(lidar_data[,!colnames(lidar_data) %in% c("ym","xm","dif_can","canop08","tree_density2008","AGBloss_mort2008","ndeath2008")],2,as.numeric)
  # Taking the square root of the slope variance:
  lidar_data[,9] = sqrt(lidar_data[,9])
  # data_abiotic[,5] = log(data_abiotic[,5])
  # colnames(lidar_data) = c("Canopy.height","Topography","Light","Wetness","Slope","Canopy.height.difference.2008-2012","Tree.density","Number.of.tree.deaths.2008-2012","Loss.of.above-ground.biomass.2008-2012","Slope.standard.deviation")
  colnames(lidar_data) = c("Canopy height (m)","Topography (m a.s.l.)","Light","Topographic Wetness Index","Slope","Tree density",
                           "Nb. of recent tree deaths","Recent loss of a.g.b.","Slope std. dev.")
  # Reordering the stations to match documents:
  #lidar_data.frame = list()
  for (j in 1:ncol(lidar_data))
    lidar_data[,j]  = as.vector(t(matrix(data=lidar_data[,j],nrow=nrow_grid)))
    #lidar_data.frame[[j]]  = t(matrix(data=lidar_data[,j],nrow=nrow_grid))
  #########
  
  color.pal = colorRampPalette(c("#00007F", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),space = "Lab")
  grain = 10
  # for (j in 1:ncol(lidar_data))
  #   lidar_data1[,j]  = as.vector(t(matrix(data=lidar_data[,j],nrow=nrow_grid)))
  coord = expand.grid(x = grain*seq(1,ncol_grid,1), y = grain*seq(1,nrow_grid,1))
  # coord = expand.grid(y = grain*seq(1,nrow_grid,1), x = grain*seq(1,ncol_grid,1))[,c(2,1)]
  newdata = expand.grid(x = seq(min(coord$x) - grain, max(coord$x) + grain, 1), y = seq(min(coord$y) - grain, max(coord$y) + grain, 1))
  # newdata = expand.grid(y = seq(min(coord$y) - grain, max(coord$y) + grain, 1), x = seq(min(coord$x) - grain, max(coord$x) + grain, 1))[,c(2,1)]
  lidar.plot = list()
  #spatial_topicmix_kriged = list()
  for (j in 1:ncol(lidar_data))
  {
    # Building a kriged spatial map of the grid, topic by topic
    # tmp = data.frame(x=coord$x,y=coord$y,z=lidar_data[,j])
    # mod = gstat(id = "z", formula = z~1, locations = ~x+y, data=tmp, model=vgm(10, "Exp", 20))
    # spatial_topicmix_kriged[[j]] = predict(mod, newdata, na.action=na.omit)
    
    lidar.plot[[j]] = ggplot(data=spatial_topicmix_kriged[[j]]) +
      geom_raster(aes(x,y,fill=z.pred)) +
      # qplot(x, y, data=spatial_topicmix_kriged, geom="raster", fill=z.pred) +
      scale_fill_gradientn(colours=color.pal(7)) +
      coord_equal() + 
      theme_minimal() +
      labs(fill = colnames(lidar_data)[j]) + 
      ggtitle(LETTERS[j]) +
      #ggtitle(if (abiotic_pca) paste("Axis #",j," - Eigenvalue = ",format(lidar_pca$eig[j]/sum(lidar_pca$eig)*nrow(lidar_pca$co),digits=4)) else colnames(lidar_data)[j]) +
      theme(legend.position="bottom", 
            legend.text=element_text(size=11),
            legend.title=element_text(size=14),
            axis.title=element_blank(), 
            axis.text = element_blank(),
            plot.title=element_text(hjust=0, size = 14), 
            plot.margin=unit(c(0.1,0.1,0.1,0.1),"mm")) +
            #plot.margin=unit(c(0,1,-2,2),"mm")) +
      guides(fill = guide_colorbar(barwidth = 8, barheight = 0.4, title.position="bottom")) +
      scale_y_continuous(limits = c(min(coord$y) - floor(grain/2), max(coord$y) + floor(grain/2)), expand = c(0,0)) +
      scale_x_continuous(limits = c(min(coord$x) - floor(grain/2), max(coord$x) + floor(grain/2)), expand = c(0,0)) + 
      geom_point(data = data.frame(coord,z.pred=rep(0,nrow(coord))), aes(x,y), color="black", size=0.6, alpha=0.3)
  }
  
  plot = grid.arrange(grobs = lidar.plot, nrow = 3, layout_matrix = matrix(1:9,nrow=3,byrow=T))
  ggsave(filename = paste0(figure_folder,"/Lidar_spatial_maps_figS2.pdf"),
         plot=plot,width=7.5,height=12)
}
####################################################################################################
####################################################################################################
if (stability)
{
  library(ggplot2)
  library(gridExtra)
  
  #   barcode_insert_list = c("Bacteries_16S","Archees_16S","Protistes_18S","Champignons_ITS","Plantes_GH","Arthropodes_18S")
  #   barcode_labels_list = c("Bacteria 16S","Archaea 16S","Protists 18S","Fungi ITS","Plants trnL","Arthropods 18S")
  
  # All barcodes :
  #   barcode_insert_list = c("Bacteries_16S","Archees_16S","Protistes_18S","Champignons_ITS","Glomeromycetes_ITS","Basidiomycetes_ITS","Ascomycetes_ITS","Plantes_GH","Arthropodes_18S","Nematodes_18S","PLatyhelminthes_18S","Annelides_18S")
  #   barcode_labels_list = c("Bacteria 16S","Archaea 16S","Protists 18S","Fungi ITS","Glomeromycota ITS","Basidiomycota ITS","Ascomycota ITS","Plants trnL","Arthropods 18S","Nematodes 18S","Platyhelminthes 18S","Annelids 18S")
  
  #     barcode_insert_list = c("Bacteries_16S","Protistes_18S","Champignons_ITS")
  #     barcode_labels_list = c("Bacteria","Protists","Fungi")
  
  # barcode_insert_list = c("Bacteries_16S","Protistes_18S","Champignons_ITS","Arthropodes_18S","Nematodes_18S","PLatyhelminthes_18S","Annelides_18S")
  # barcode_labels_list = c("Bacteria","Protists","Fungi","Arthropods","Nematodes","Platyhelminthes","Annelids")
  
  barcode_insert_list = c("Bacteries_16S","Protistes_18S","Arthropodes_18S","Nematodes_18S","PLatyhelminthes_18S","Annelides_18S")
  barcode_labels_list = c("Bacteria","Protists","Arthropods","Nematodes","Platyhelminthes","Annelids")
  
  #   barcode_insert_list = c("Champignons_ITS","Bacteries_16S","Arthropodes_18S")
  #   barcode_labels_list = c("Fungi","Bacteria","Arthropods") 
  
  # barcode_insert_list = c("Champignons_ITS","Bacteries_16S")
  # barcode_labels_list = c("Fungi","Bacteria")  
  
  #     barcode_insert_list = c("Bacteries_16S")
  #     barcode_labels_list = c("Bacteria")
  
  nb_topics = 3
  nb_real = 100
  em_tol = 10^-7
  var_tol = 10^-8
  
  occurrence = 1
  
  bij = 1
  samplewise = 1
  MOTUwise = 0
  
  if (bij) {
    bij_insert = "bijective_correspondence"
  } else if (!bij)
    bij_insert = "non-bijective_correspondence"
  
  if (samplewise) {
    MOTU_sample_insert = "samplewise"
  } else if (MOTUwise)
    MOTU_sample_insert = "MOTUwise"
  
  if (occurrence)
    occurrence_insert = "_occurrence"
  else
    occurrence_insert = ""
  
  #   pdf("SKL_vs_rank.pdf",width=12.5/2,height=12.5/4)
  #   par(mfrow=c(2,3))
  #   
  #   par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
  #   # par(mar = c(bottom, left, top, right))
  #   par(mar = c(4, 4, 3, 1) + 0.1)
  
  plot.skl = list()
  plot.skl.llh = list()
  plot.dskl.llh = list()
  plot.dskln.llh = list()
  plot.pval.llh = list()
  plot.ses.llh = list()
  
  plot.skl.all = list()
  plot.pval.all = list()
  
  barcode_index = 0
  for (barcode_insert in barcode_insert_list)
  {
    barcode_index = barcode_index+1
    
    if (occurrence)
    {
      if (barcode_insert == "Bacteries_16S")
      {
        occurrence_insert = "_occurrence_54r"
        length_selected_real = 54
      }
      else
      {
        occurrence_insert = "_occurrence"
        length_selected_real = 100
      }
    } else
    {
      occurrence_insert = ""
      length_selected_real = 100
    }
    
    data_insert = "Donnees_PetitPlateau"
    filename_insert = "Rtopicmodels_LDA_VEM"
    local_prefix = "/Users/guilhemsommeria-klein/Desktop/These/" 
    local_dirname = paste(local_prefix,data_insert,"/",barcode_insert,"/",filename_insert,"/",sep="")
    local_subdirname = paste(local_dirname,filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best_keep",occurrence_insert,"/",sep="")
    realization_comparison_dirname = paste0(local_subdirname,"Realization_comparison_",MOTU_sample_insert,"_",bij_insert)
    #     subsubdirname = paste(local_subdirname,"ordered_realization_number_1/",sep="")  
    #     subsubsubdirname = paste(subsubdirname,"/topics_site_repartition_info/",sep="")  
    
    setwd(realization_comparison_dirname)
    #load("sort_normal_topic.Rdata")
    #filename1 = paste(filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best_keep.Rdata",sep="")
    #filename = "Topic_KL-correlation_matrices_between_all_real.Rdata"
    filename = paste0("SKL-correlation_allRealPairs_",MOTU_sample_insert,"_",bij_insert,".Rdata")
    load(filename)
    #Ordered_realizations = readRDS("Ordered_realizations.rds")
    data.skl = data.frame(x=seq(2,length_selected_real,1),y=KL_allRealPairs_w_rndzations[1,-1])
    #     if ((barcode_insert == "Bacteries_16S") && occurrence)
    #       data.skl.llh = data.frame(x=llh_differences_allRealPairs[1,c(-1,-ncol(llh_differences_allRealPairs))],y=KL_allRealPairs_w_rndzations[1,c(-1,-ncol(KL_allRealPairs_w_rndzations))])
    #     else
    data.skl.llh = data.frame(x=llh_differences_allRealPairs[1,-1],y=KL_allRealPairs_w_rndzations[1,-1])
    data.dskl.llh = data.frame(x=llh_differences_allRealPairs[1,-1],y=DKL100_allRealPairs[1,-1])
    data.dskln.llh = data.frame(x=llh_differences_allRealPairs[1,-1],y=DKL100_allRealPairs[1,-1]/KL_allRealPairs_w_rndzations[-1,1])
    data.ses.llh = data.frame(x=llh_differences_allRealPairs[1,-1],y=SES_allRealPairs[1,-1])
    #data.skl.llh = data.frame(x=Ordered_realizations$x[1]-Ordered_realizations$x[2:length_selected_real],y=KL_samplewise_allRealPairs_w_rndzations[1,-1])
    data.pval.llh = data.frame(x=llh_differences_allRealPairs[1,-1],y=p_value_allRealPairs[1,-1])
    
    data.skl.all = data.frame(x=llh_differences_allRealPairs[upper.tri(llh_differences_allRealPairs,diag=F)],y=KL_allRealPairs_w_rndzations[upper.tri(KL_allRealPairs_w_rndzations,diag=F)])
    data.pval.all = data.frame(x=llh_differences_allRealPairs[upper.tri(llh_differences_allRealPairs,diag=F)],y=p_value_allRealPairs[upper.tri(p_value_allRealPairs,diag=F)])
    
    plot.skl[[barcode_index]] = qplot(x, y, data = data.skl, geom="point") +
      #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", function(x) round(x,0))) + 
      #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", scientific_format())) +
      theme_bw() + ggtitle(paste(letters[barcode_index],"-",barcode_labels_list[barcode_index])) +
      #ylim(0,max(data.skl$y)) 
      # for bij samplewise:  
      #ylim(0,4.5)
      # for MOTUwise:
      #ylim(0,8)
      #for !bij samplewise:
      ylim(0,3.5)
    #     scale_y_continuous(limits=c(5,395), expand = c(0,0)) +
    #     scale_x_continuous(limits=c(5,295), expand = c(0,0)) + 
    #     theme(legend.position="bottom", legend.text=element_text(size=11), 
    #         legend.title=element_text(size=12), axis.title=element_blank(), 
    #         axis.text = element_blank(),
    #         plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm"))
    plot.skl.llh[[barcode_index]] = qplot(x, y, data = data.skl.llh, geom="point") +
      theme_bw() + ggtitle(paste(letters[barcode_index],"-",barcode_labels_list[barcode_index])) +
      #ylim(0,max(data.skl$y)) 
      # for bij samplewise:  
      ylim(0,4.5) + xlim(0,35000) 
    # for MOTUwise:
    #ylim(0,8)
    #for !bij samplewise: 
    #ylim(0,3.5)
    #for bij samplewise occurrence:
    #ylim(0,2.8) #+ xlim(0,3200)
    plot.dskl.llh[[barcode_index]] = qplot(x, y, data = data.dskl.llh, geom="point") +
      theme_bw() + ggtitle(paste(letters[barcode_index],"-",barcode_labels_list[barcode_index])) +
      ylim(min(0,min(data.dskl.llh$y)),max(data.dskl.llh$y))    
    plot.dskln.llh[[barcode_index]] = qplot(x, y, data = data.dskln.llh, geom="point") +
      theme_bw() + ggtitle(paste(letters[barcode_index],"-",barcode_labels_list[barcode_index])) +
      #ylim(min(0,min(data.dskln.llh$y)),max(data.dskln.llh$y))
      ylim(min(0,min(data.dskln.llh$y)),1)
    plot.pval.llh[[barcode_index]] = qplot(x, y, data = data.pval.llh, geom="point") +
      theme_bw() + ggtitle(paste(letters[barcode_index],"-",barcode_labels_list[barcode_index])) +
      ylim(0,max(data.pval.llh$y)) 
    plot.ses.llh[[barcode_index]] = qplot(x, y, data = data.ses.llh, geom="point") +
      theme_bw() + ggtitle(paste(letters[barcode_index],"-",barcode_labels_list[barcode_index])) +
      ylim(0,max(data.ses.llh$y)) 
    
    plot.skl.all[[barcode_index]] = qplot(x, y, data = data.skl.all, geom="point") +
      theme_bw() + ggtitle(paste(letters[barcode_index],"-",barcode_labels_list[barcode_index])) +
      ylim(0,max(data.skl.all$y))
    plot.pval.all[[barcode_index]] = qplot(x, y, data = data.pval.all, geom="point") +
      theme_bw() + ggtitle(paste(letters[barcode_index],"-",barcode_labels_list[barcode_index])) +
      ylim(0,max(data.pval.all$y))
    
    plot.skl[[barcode_index]] = plot.skl[[barcode_index]] + theme(axis.title=element_text(size=11), 
                                                                  plot.title=element_text(hjust=0, size=12),
                                                                  plot.margin=unit(c(1,1,1,0),"mm"))
    #plot.margin=unit(c(1,1.5,2,0),"mm"))
    plot.skl.llh[[barcode_index]] = plot.skl.llh[[barcode_index]] + theme(axis.title=element_text(size=11), 
                                                                          plot.title=element_text(hjust=0, size=12),
                                                                          plot.margin=unit(c(1,1,1,0),"mm")) 
    #if (!(barcode_insert == "Bacteries_16S") || !occurrence)
    plot.skl.llh[[barcode_index]] = plot.skl.llh[[barcode_index]] + geom_smooth(method='lm',se=F,linetype="dashed",size=0.5)
    
    plot.dskl.llh[[barcode_index]] = plot.dskl.llh[[barcode_index]] + theme(axis.title=element_text(size=11), 
                                                                            plot.title=element_text(hjust=0, size=12),
                                                                            plot.margin=unit(c(1,1,1,0),"mm"))
    plot.dskln.llh[[barcode_index]] = plot.dskln.llh[[barcode_index]] + theme(axis.title=element_text(size=18), 
                                                                              plot.title=element_text(hjust=0, size=18),
                                                                              plot.margin=unit(c(1,1,1,0),"mm"),
                                                                              axis.text=element_text(size=14),
                                                                              axis.line=element_line(size=0.8)) +
      geom_smooth(method='lm',se=F,linetype="dashed",size=0.8)
    plot.pval.llh[[barcode_index]] = plot.pval.llh[[barcode_index]] + theme(axis.title=element_text(size=11), 
                                                                            plot.title=element_text(hjust=0, size=12),
                                                                            plot.margin=unit(c(1,1,1,0),"mm"))
    plot.ses.llh[[barcode_index]] = plot.ses.llh[[barcode_index]] + theme(axis.title=element_text(size=11), 
                                                                          plot.title=element_text(hjust=0, size=12),
                                                                          plot.margin=unit(c(1,1,1,0),"mm"))
    
    plot.skl.all[[barcode_index]] = plot.skl.all[[barcode_index]] + theme(axis.title=element_text(size=11), 
                                                                          plot.title=element_text(hjust=0, size=12),
                                                                          plot.margin=unit(c(1,1,1,0),"mm"))
    plot.pval.all[[barcode_index]] = plot.pval.all[[barcode_index]] + theme(axis.title=element_text(size=11), 
                                                                            plot.title=element_text(hjust=0, size=12),
                                                                            plot.margin=unit(c(1,1,1,0),"mm"))
    
    #if (barcode_index == 1)
    if (barcode_index == 1 || barcode_index == 3 || barcode_index == 5)
    {
      plot.skl[[barcode_index]] = plot.skl[[barcode_index]] + labs(x="Realizations\n sorted by decreasing likelihood", y="SKL to best realization")
      plot.skl.llh[[barcode_index]] = plot.skl.llh[[barcode_index]] + labs(x="Llh difference\n with best realization", y="SKL to best realization")
      #dskl.plot.llh[[barcode_index]] = dskl.plot.llh[[barcode_index]] + labs(x="Loglikelihood difference\n with best realization", y=expression(paste("SKL-<SKL>",[random.]*" to best realization"))) bquote(.(labNames[1]) ~ x^2)
      plot.dskl.llh[[barcode_index]] = plot.dskl.llh[[barcode_index]] + labs(x="Llh difference\n with best realization", y=bquote(.("SKL-<SKL>")[random.]*" to best realization"))
      #plot.dskln.llh[[barcode_index]] = plot.dskln.llh[[barcode_index]] + labs(x="Llh difference\n with best realization", y=bquote(.("(SKL-<SKL>")[rnd]*")/<SKL>"[rnd]))
      plot.dskln.llh[[barcode_index]] = plot.dskln.llh[[barcode_index]] + labs(x="Llh difference\n with best realization\n", y="Spatial similarity\n to best realization")
      plot.pval.llh[[barcode_index]] = plot.pval.llh[[barcode_index]] + labs(x="Llh difference\n with best realization", y="p-value")
      plot.ses.llh[[barcode_index]] = plot.ses.llh[[barcode_index]] + labs(x="Llh difference\n with best realization", y="SES")
      
      plot.skl.all[[barcode_index]] = plot.skl.all[[barcode_index]] + labs(x="Llh difference", y="SKL")
      plot.pval.all[[barcode_index]] = plot.pval.all[[barcode_index]] + labs(x="Llh difference", y="p-value")
    }
    else
    {
      plot.skl[[barcode_index]] = plot.skl[[barcode_index]] + labs(x="\n", y="")
      plot.skl.llh[[barcode_index]] = plot.skl.llh[[barcode_index]] + labs(x="\n", y="")
      plot.dskl.llh[[barcode_index]] = plot.dskl.llh[[barcode_index]] + labs(x="\n", y="")
      #plot.dskln.llh[[barcode_index]] = plot.dskln.llh[[barcode_index]] + labs(x="\n", y="")
      plot.dskln.llh[[barcode_index]] = plot.dskln.llh[[barcode_index]] + labs(x="Llh difference\n with best realization\n", y="")
      plot.pval.llh[[barcode_index]] = plot.pval.llh[[barcode_index]] + labs(x="\n", y="")
      plot.ses.llh[[barcode_index]] = plot.ses.llh[[barcode_index]] + labs(x="\n", y="")
      
      plot.skl.all[[barcode_index]] = plot.skl.all[[barcode_index]] + labs(x="\n", y="")
      plot.pval.all[[barcode_index]] = plot.pval.all[[barcode_index]] + labs(x="\n", y="")
    }
    #guides(fill = guide_colorbar(barwidth = 8, barheight = 0.4, title.position="bottom"))
    
    #     if (barcode_index == 1)
    #     {
    #       plot(seq(2,length_selected_real,1),KL_samplewise_allRealPairs_w_rndzations[1,-1],type="p",
    #         xlab = "Realizations\n sorted by decreasing llh", ylab = "SKL to best realization", main = barcode_labels_list[barcode_index],
    #         ylim = c(0,max(KL_samplewise_allRealPairs_w_rndzations[1,-1])))
    #       mtext(text="a",side=3,line=1.3,cex=1,at=xmin-(xmax-xmin)*0.1,font=2)
    #     }
    #     else 
    #       plot(seq(2,length_selected_real,1),KL_samplewise_allRealPairs_w_rndzations[1,-1],type="p", 
    #            main = barcode_labels_list[barcode_index], xlab="", ylab="",
    #            ylim = c(0,max(KL_samplewise_allRealPairs_w_rndzations[1,-1])))
  }
  
  setwd("/Users/guilhemsommeria-klein/Desktop/These/Manuscrits/ThesisPaper2_Topic_Modelling_PP/Figures/")
  #   ggsave(filename = paste0("SKL_vs_rank",occurrence_insert,"_",MOTU_sample_insert,"_",bij_insert,".pdf"), do.call("arrangeGrob", 
  #                                                                                              c(plot.skl, nrow=2)),
  #          height = 10*2/3, width = 10)
  #   ggsave(filename = paste0("SKL_vs_llh-diff_bact-prot-fungi",occurrence_insert,"_",MOTU_sample_insert,"_",bij_insert,"_llh0-3200.pdf"), do.call("arrangeGrob", 
  #                                                                                                  c(plot.skl.llh, nrow=1)),
  #          height = 12/3, width = 12)
  #   ggsave(filename = paste0("SKL_vs_llh-diff_bact-prot-fungi",occurrence_insert,"_",MOTU_sample_insert,"_",bij_insert,"_llh0-35000.pdf"), do.call("arrangeGrob", 
  #                                                                                                                                     c(plot.skl.llh, nrow=1)),
  #        height = 12/3, width = 12)
  #   ggsave(filename = paste0("SKL_vs_llh-diff_bact-prot-fungi",occurrence_insert,"_",MOTU_sample_insert,"_",bij_insert,".pdf"), do.call("arrangeGrob", 
  #                                                                                                                                     c(plot.skl.llh, nrow=1)),
  #        height = 12/3, width = 12)
  #   ggsave(filename = paste0("DSKL_vs_llh-diff",occurrence_insert,"_",MOTU_sample_insert,"_",bij_insert,".pdf"), do.call("arrangeGrob", 
  #                                                                                                   c(plot.dskl.llh, nrow=2)),
  #          height = 10*2/3, width = 10)
  #   ggsave(filename = paste0("p-value_vs_llh-diff",occurrence_insert,"_",MOTU_sample_insert,"_",bij_insert,".pdf"), do.call("arrangeGrob", 
  #                                                                                                      c(plot.pval.llh, nrow=2)),
  #          height = 10*2/3, width = 10)
  #   ggsave(filename = paste0("SES_vs_llh-diff",occurrence_insert,"_",MOTU_sample_insert,"_",bij_insert,".pdf"), do.call("arrangeGrob", 
  #                                                                                                                         c(plot.ses.llh, nrow=2)),
  #          height = 10*2/3, width = 10)
  # ggsave(filename = paste0("SES_vs_llh-diff_allBarcodes",occurrence_insert,"_",MOTU_sample_insert,"_",bij_insert,".pdf"), do.call("arrangeGrob", 
  #                                                                                                                     c(plot.ses.llh, nrow=2)),
  #        height = 15/3, width = 15)
  # ggsave(filename = paste0("NormES_vs_llh-diff_allBarcodes_withLM",occurrence_insert,"_",MOTU_sample_insert,"_",bij_insert,".pdf"), do.call("arrangeGrob", 
  #                                                                                                                                 c(plot.dskln.llh, nrow=2)),
  #        height = 15/3, width = 15)
  
  # ggsave(filename = paste0("NormES_vs_llh-diff_-restrictedBact-prot-fungi-arth-nem-platy-anne_withLM",occurrence_insert,"_",MOTU_sample_insert,"_",bij_insert,".pdf"), do.call("arrangeGrob", 
  #                                                                                                                                                         c(plot.dskln.llh, nrow=4)),
  #        height = 10/2*4, width = 10)
  ggsave(filename = paste0("NormES_vs_llh-diff_-restrictedBact-prot-arth-nem-platy-anne_withLM",occurrence_insert,"_",MOTU_sample_insert,"_",bij_insert,".pdf"), do.call("arrangeGrob", 
                                                                                                                                                                         c(plot.dskln.llh, nrow=3)),
         height = 10/2*3, width = 10)
  
  #   ggsave(filename = paste0("SKL_vs_llh-diff_allReal",occurrence_insert,"_",MOTU_sample_insert,"_",bij_insert,".pdf"), do.call("arrangeGrob", 
  #                                                                                                                     c(plot.skl.all, nrow=2)),
  #        height = 10*2/3, width = 10)
  #   ggsave(filename = paste0("p-value_vs_llh-diff_allReal",occurrence_insert,"_",MOTU_sample_insert,"_",bij_insert,".pdf"), do.call("arrangeGrob", 
  #                                                                                                                     c(plot.pval.all, nrow=2)),
  #        height = 10*2/3, width = 10)
  #dev.off()
  #   
  #   pdf("Topic_KL_difference_samplewise_between_real_averaged_w_randomizations_vs_llh_difference.pdf")  
  #   par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
  #   # par(mar = c(bottom, left, top, right))
  #   par(mar = c(5, 5, 4, 3) + 0.1)
  #   plot(Ordered_realizations$x[1]-Ordered_realizations$x[2:length_selected_real],DKL100_samplewise_allRealPairs[1,-1],type="p", main = "Averaged samplewise KL distance diffrence between\n the component communities in different realizations",
  #        xlab = "Loglikelihood difference with best realization", ylab = "Difference in KL distance to best realization samplewise")
  #   dev.off()
  
}
####################################################################################################
####################################################################################################
if (size_diversity)
{
  library(plotrix)
  library(ggplot2)
  library(gridExtra)
  
  setwd("/Users/guilhemsommeria-klein/Desktop/These/Manuscrits/ThesisPaper2_Topic_Modelling_PP/Figures/")
  # SKL_bij_samplewise = c(0.27,0.58,1.9,1.5,1.1,2.0,3.0,3.1,3.0)
  # ES_bij_samplewise = c(2.0,1.5,3.1,1.6,0.86,1.6,2.5,1.4,0.95)
  # SES_bij_samplewise = c(1.3,1.3,0.81,0.47,0.93,0.95,0.75,0.86,1.5)
  # pval_bij_samplewise = c(0,0.029,0.069,0.17,0.082,0.13,0.042,0.091,0.092)
  # 
  # SKL_bij_samplewise1 = c(0.27,0.58,1.9,0.025,1.1,2.0,3.0,3.1,3.0)
  # ES_bij_samplewise1 = c(2.0,1.5,3.1,2.7,0.86,1.6,2.5,1.4,0.95)
  # SES_bij_samplewise1 = c(1.3,1.3,0.81,3.5,0.93,0.95,0.75,0.86,1.5)
  # pval_bij_samplewise1 = c(0,0.029,0.069,0,0.082,0.13,0.042,0.091,0.092)
  
  #Occ.
  #   nES_bij_samplewise = c(0.32,0.68,0.41,0.88,0.33,0.52,0.62,0.70,0.10)
  #   nES_bij_samplewise1 = c(0.85,0.68,0.41,0.88,0.33,0.52,0.62,0.70,0.10)
  # nES_bij_samplewise = c(0.32,0.68,0.41,0.33,0.52,0.62,0.70)
  # nES_bij_samplewise1 = c(0.85,0.68,0.41,0.33,0.52,0.62,0.70)
  # Without Fungi ITS:
  nES_bij_samplewise = c(0.32,0.68,0.41,0.33,0.52,0.62)
  nES_bij_samplewise1 = c(0.85,0.68,0.41,0.33,0.52,0.62)
  
  # SKL_nobij_samplewise = c(0.18,0.48,1.2,0.44,0.80,1.3,2.5,2.4,2.8)
  # ES_nobij_samplewise = c(1.6,1.5,2.5,1.8,0.75,1.7,2.3,1.1,0.97)
  # SES_nobij_samplewise = c(1.2,1.2,0.68,1.4,1.4,1.2,0.71,1.0,1.7)
  # pval_nobij_samplewise = c(0,0.032,0.19,0.041,0.072,0.074,0.15,0.12,0.065)
  # 
  # SKL_nobij_samplewise1 = c(0.18,0.48,1.2,0.025,0.80,1.3,2.5,2.4,2.8)
  # ES_nobij_samplewise1 = c(1.6,1.5,2.5,2.7,0.75,1.7,2.3,1.1,0.97)
  # SES_nobij_samplewise1 = c(1.2,1.2,0.68,3.5,1.4,1.2,0.71,1.0,1.7)
  # pval_nobij_samplewise1 = c(0,0.032,0.19,0,0.072,0.074,0.15,0.12,0.065)
  # 
  # SKL_bij_MOTUwise = c(0.11,0.56,2.6,2.2,2.4,2.8,3.6,2.2,5.4)
  # ES_bij_MOTUwise = c(7.0,6.3,7.8,9.8,5.2,5.6,6.7,4.3,5.9)
  # SES_bij_MOTUwise = c(9.5,4.4,2.7,5.1,4.7,3.2,3.4,9.6,7.3)
  # pval_bij_MOTUwise = c(0,0,0.045,0.032,0.0091,0.028,0.0071,0,0.00056)
  # 
  # SKL_bij_MOTUwise1 = c(0.11,0.56,2.6,0.031,2.4,2.8,3.6,2.2,5.4)
  # ES_bij_MOTUwise1 = c(7.0,6.3,7.8,12,5.2,5.6,6.7,4.3,5.9)
  # SES_bij_MOTUwise1 = c(9.5,4.4,2.7,7,4.7,3.2,3.4,9.6,7.3)
  # pval_bij_MOTUwise1 = c(0,0,0.045,0,0.0091,0.028,0.0071,0,0.00056)
  
  #Occ.
  #   nES_bij_MOTUwise = c(0.86,0.95,0.83,0.98,0.88,0.86,0.91,0.90,0.85)
  #   nES_bij_MOTUwise1 = c(0.95,0.95,0.83,0.98,0.88,0.86,0.91,0.90,0.85)
  # 
  # nES_bij_MOTUwise = c(0.86,0.95,0.83,0.88,0.86,0.91,0.90)
  # nES_bij_MOTUwise1 = c(0.95,0.95,0.83,0.88,0.86,0.91,0.90)
  # Without Fungi ITS:
  nES_bij_MOTUwise = c(0.86,0.95,0.83,0.88,0.86,0.91)
  nES_bij_MOTUwise1 = c(0.95,0.95,0.83,0.88,0.86,0.91)
  
  #taxo_names = c("Bacteria","Protists","Annelids","Archaea","Nematodes","Platyhelminthes","Arthropods","Fungi","Plants")
  #taxo_names = c("Bacteria","Protists","Annelids","Archaea","Nematodes","Platyhelminthes","Arthropods","Fungi")
  #taxo_names = c("Bacteria","Protists","Annelids","Nematodes","Platyhelminthes","Arthropods","Fungi")
  taxo_names = c("Bacteria","Protists","Annelids","Nematodes","Platyhelminthes","Arthropods")
  #   size = log10(c(mean(c(0.5*10^-6,5*10^-6,0.5*10^-6,10^-6,10^-6,0.5*10^-6,5*10^-6,0.5*10^-6,5*10^-6,0.5*10^-6)),
  #                  100*10^-6,
  #                  20e-3,
  #                  0.5*10^-6,
  #                  100*10^-6,
  #                  20*10^-3,
  #                  10^-3,
  #                  mean(c(100*10^-6,10*10^-6,200*10^-6)),
  #                  0.5*10^2))
  # size = log10(c(mean(c(0.5*10^-6,5*10^-6,0.5*10^-6,10^-6,10^-6,0.5*10^-6,5*10^-6,0.5*10^-6,5*10^-6,0.5*10^-6)),
  #                100*10^-6,
  #                20e-3,
  #                100*10^-6,
  #                20*10^-3,
  #                10^-3,
  #                mean(c(100*10^-6,10*10^-6,200*10^-6))))
  # Without fungi ITS:
  size = log10(c(mean(c(0.5*10^-6,5*10^-6,0.5*10^-6,10^-6,10^-6,0.5*10^-6,5*10^-6,0.5*10^-6,5*10^-6,0.5*10^-6)),
                 100*10^-6,
                 20e-3,
                 100*10^-6,
                 20*10^-3,
                 10^-3))
  
  #diversity0 = log10(c(20162,1648,51,4101,378,126,1881,9855,1360))
  #diversity = log10(c(20162,1648,51,378,126,1881,9855))
  diversity = log10(c(20162,1648,51,378,126,1881))
  
  data.tmp = list()
  #   data.tmp[[1]] = data.frame(x=size,y=SKL_bij_samplewise)
  #   data.tmp[[2]] = data.frame(x=size,y=ES_bij_samplewise)
  #   data.tmp[[3]] = data.frame(x=size,y=SES_bij_samplewise)
  #   data.tmp[[4]] = data.frame(x=size,y=pval_bij_samplewise)
  #   data.tmp[[5]] = data.frame(x=size,y=SKL_bij_MOTUwise)
  #   data.tmp[[6]] = data.frame(x=size,y=ES_bij_MOTUwise)
  #   data.tmp[[7]] = data.frame(x=size,y=SES_bij_MOTUwise)
  #   data.tmp[[8]] = data.frame(x=size,y=pval_bij_MOTUwise)
  #   data.tmp[[9]] = data.frame(x=size,y=SKL_nobij_samplewise)
  #   data.tmp[[10]] = data.frame(x=size,y=ES_nobij_samplewise)
  #   data.tmp[[11]] = data.frame(x=size,y=SES_nobij_samplewise)
  #   data.tmp[[12]] = data.frame(x=size,y=pval_nobij_samplewise)
  data.tmp[[1]] = data.frame(x=size,y=nES_bij_samplewise)
  data.tmp[[2]] = data.frame(x=size,y=nES_bij_MOTUwise)
  
  data.tmp1 = list()
  #   data.tmp1[[1]] = data.frame(x=size,y=SKL_bij_samplewise1)
  #   data.tmp1[[2]] = data.frame(x=size,y=ES_bij_samplewise1)
  #   data.tmp1[[3]] = data.frame(x=size,y=SES_bij_samplewise1)
  #   data.tmp1[[4]] = data.frame(x=size,y=pval_bij_samplewise1)
  #   data.tmp1[[5]] = data.frame(x=size,y=SKL_bij_MOTUwise1)
  #   data.tmp1[[6]] = data.frame(x=size,y=ES_bij_MOTUwise1)
  #   data.tmp1[[7]] = data.frame(x=size,y=SES_bij_MOTUwise1)
  #   data.tmp1[[8]] = data.frame(x=size,y=pval_bij_MOTUwise1)
  #   data.tmp1[[9]] = data.frame(x=size,y=SKL_nobij_samplewise1)
  #   data.tmp1[[10]] = data.frame(x=size,y=ES_nobij_samplewise1)
  #   data.tmp1[[11]] = data.frame(x=size,y=SES_nobij_samplewise1)
  #   data.tmp1[[12]] = data.frame(x=size,y=pval_nobij_samplewise1)
  #   data.tmp1[[1]] = data.frame(x=size[-length(size)],y=nES_bij_samplewise1[-length(nES_bij_samplewise1)])
  #   data.tmp1[[2]] = data.frame(x=size[-length(size)],y=nES_bij_MOTUwise1[-length(nES_bij_MOTUwise1)])
  data.tmp1[[1]] = data.frame(x=size,y=nES_bij_samplewise1)
  data.tmp1[[2]] = data.frame(x=size,y=nES_bij_MOTUwise1)
  
  #   plot_labels = c(rep("One-to-one correpondence -\n spatial distribution",4),rep("One-to-one correspondence -\n MOTU composition",4),rep("K best correspondances -\n spatial distribution",4),
  #                   "Spatial similarity","Taxonomic similarity")
  #y_labels = rep(c("SKL","ES","SES","p-value"),3)
  #y_labels = c(rep(c("SKL","ES","SES","Stability p-value"),3),"Spatial stability","Taxonomic stability")
  y_labels = c("Spatial stability","Taxonomic stability")
  
  plot.size = list()
  plot.size1 = list()
  plot.diversity = list()
  plot.diversity1 = list()
  rsquared = list()
  # i == 1: spatial stability
  # i == 2: taxonomic stability
  # rsquared[[i]] : size-stability ; size-stability with restricted bact ; diversity-stability ; diversity-stability with restricted bact
  for (i in 1:2)
  {
    #diversity = diversity0[-length(diversity0)]
    
    rsquared[[i]] = vector(length=4,mode="numeric")
    
    plot.size[[i]] = qplot(x, y, data = data.tmp[[i]], geom="point",label=taxo_names) +
      #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", function(x) round(x,0))) + 
      #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", scientific_format())) +
      theme_bw() + #ggtitle(paste(letters[i],"-",plot_labels[i])) +
      #ylim(0,max(data.skl$y)) 
      #     scale_y_continuous(limits=c(5,395), expand = c(0,0)) +
      #     scale_x_continuous(limits=c(5,295), expand = c(0,0)) + 
      #     theme(legend.position="bottom", legend.text=element_text(size=11), 
      #         legend.title=element_text(size=12), axis.title=element_blank(), 
      #         axis.text = element_blank(),
      #         plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm"))
      theme(axis.title=element_text(size=8), 
            plot.title=element_text(hjust=0, size=8),
            plot.margin=unit(c(1,1,1,0),"mm")) +
      labs(x="Size (log10)", y=y_labels[i]) +
      xlim(min(data.tmp[[i]]$x),3) +
      #       geom_text(mapping = NULL, data = data.frame(x=data.tmp[[i]]$x,y=data.tmp[[i]]$y,label=taxo_names), x=x, y=y, label=label, stat = "identity", position = "identity", 
      #               parse = FALSE, check_overlap = T, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE)
      geom_text(mapping = NULL, stat = "identity",size=2, hjust = 0, nudge_x = 0.2,
                parse = FALSE, check_overlap = F, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) +
      geom_smooth(method='lm')
    
    rsquared[[i]][1] = cor(data.tmp[[i]][,2],data.tmp[[i]][,1])^2
    
    plot.size1[[i]] = qplot(x, y, data = data.tmp1[[i]], geom="point", label=taxo_names) +
      #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", function(x) round(x,0))) + 
      #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", scientific_format())) +
      theme_bw() +
      #ggtitle(paste(letters[i],"-",plot_labels[i])) +
      #ylim(0,max(data.skl$y)) 
      #     scale_y_continuous(limits=c(5,395), expand = c(0,0)) +
      #     scale_x_continuous(limits=c(5,295), expand = c(0,0)) + 
      #     theme(legend.position="bottom", legend.text=element_text(size=11), 
      #         legend.title=element_text(size=12), axis.title=element_blank(), 
      #         axis.text = element_blank(),
      #         plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm"))
      theme(axis.title=element_text(size=9), 
            #plot.title=element_text(hjust=0, size=9),
            plot.margin=unit(c(1,1,1,0),"mm")) +
      labs(x="Organism size (log10)", y=y_labels[i]) +
      xlim(min(data.tmp1[[i]]$x),3) +
      #xlim(min(data.tmp1[[i]]$x),max(data.tmp1[[i]]$x)) +
      geom_smooth(method='lm')
    if (i == 1)
    {
      plot.size1[[i]] = plot.size1[[i]] + geom_text(mapping = NULL, stat = "identity",size=2, 
                                                    #hjust=0, nudge_x = 0.2,
                                                    #nudge_y = c(0,-0.01,0,0.02,0,0,0,0.02), hjust=0, nudge_x = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2),
                                                    nudge_y = c(0,-0.01,0,0,0,0,0.02), hjust=0, nudge_x = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2),
                                                    parse = FALSE, check_overlap = F, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) 
    } else if (i == 2)
    {
      plot.size1[[i]] = plot.size1[[i]] + geom_text(mapping = NULL, stat = "identity",size=2, 
                                                    #hjust=0, nudge_x = 0.2,
                                                    #nudge_y = c(0.01,0,0,0.01,0,0,0,0), hjust=0, nudge_x = c(0,0.2,0.2,0,0.2,0.2,0.2,0.2),
                                                    nudge_y = c(0.01,0,0,0,0,0,0), hjust=0, nudge_x = c(0,0.2,0.2,0.2,0.2,0.2,0.2),
                                                    parse = FALSE, check_overlap = F, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) 
    }
    #c("Bacteria","Protists","Annelids","Archaea","Nematodes","Platyhelminthes","Arthropods","Fungi","Plants")
    
    rsquared[[i]][2] = cor(data.tmp1[[i]][,2],data.tmp1[[i]][,1])^2
    
    plot.diversity[[i]] = qplot(x = diversity, y, data=data.tmp[[i]], label=taxo_names, geom="point") +
      #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", function(x) round(x,0))) + 
      #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", scientific_format())) +
      theme_bw() + #ggtitle(paste(letters[i],"-",plot_labels[i])) +
      #ylim(0,max(data.skl$y)) 
      #     scale_y_continuous(limits=c(5,395), expand = c(0,0)) +
      #     scale_x_continuous(limits=c(5,295), expand = c(0,0)) + 
      #     theme(legend.position="bottom", legend.text=element_text(size=11), 
      #         legend.title=element_text(size=12), axis.title=element_blank(), 
      #         axis.text = element_blank(),
      #         plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm"))
      theme(axis.title=element_text(size=9), 
            plot.title=element_text(hjust=0, size=9),
            plot.margin=unit(c(1,1,1,0),"mm")) +
      labs(x="Number of OTUs (log10)", y=y_labels[i]) +
      xlim(min(diversity),5) +
      geom_text(mapping = NULL, stat = "identity",size=2, hjust = 0, nudge_x = 0.1,
                parse = FALSE, check_overlap = F, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) +
      geom_smooth(method='lm')
    
    rsquared[[i]][3] = cor(data.tmp[[i]][,2],diversity)^2
    
    plot.diversity1[[i]] = qplot(x = diversity, y, data=data.tmp1[[i]], label=taxo_names, geom="point") +
      #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", function(x) round(x,0))) + 
      #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", scientific_format())) +
      theme_bw() +
      #ggtitle(paste(letters[i],"-",plot_labels[i])) +
      ggtitle(paste(letters[i])) +
      #ylim(0,max(data.skl$y)) 
      #     scale_y_continuous(limits=c(5,395), expand = c(0,0)) +
      #     scale_x_continuous(limits=c(5,295), expand = c(0,0)) + 
      #     theme(legend.position="bottom", legend.text=element_text(size=11), 
      #         legend.title=element_text(size=12), axis.title=element_blank(), 
      #         axis.text = element_blank(),
      #         plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm"))
      theme(axis.title=element_text(size=9), 
            plot.title=element_text(hjust=0, size=9),
            plot.margin=unit(c(1,1,1,0),"mm")) +
      labs(x="Number of OTUs (log10)", y=y_labels[i]) +
      xlim(min(diversity),5) +
      geom_smooth(method='lm') 
    
    if (i == 1)
    {
      plot.diversity1[[i]] = plot.diversity1[[i]] + geom_text(mapping = NULL, stat = "identity",size=2, 
                                                              #hjust=0, nudge_x = 0.2,
                                                              #nudge_y = c(0,0,0,0,0,0,0,0), hjust=0, nudge_x = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1),
                                                              #nudge_y = c(0,0,0,0,0,0,0), hjust=0, nudge_x = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1),
                                                              nudge_y = c(0,0,0,0,0,0), hjust=0, nudge_x = c(0.1,0.1,0.1,0.1,0.1,0.1),
                                                              parse = FALSE, check_overlap = F, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) 
    } else if (i == 2)
    {
      plot.diversity1[[i]] = plot.diversity1[[i]] + geom_text(mapping = NULL, stat = "identity",size=2, 
                                                              #hjust=0, nudge_x = 0.2,
                                                              #nudge_y = c(0,0,0,0,0,0,0), hjust=0, nudge_x = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1),
                                                              nudge_y = c(0,0,0,0,0,0), hjust=0, nudge_x = c(0.1,0.1,0.1,0.1,0.1,0.1),
                                                              parse = FALSE, check_overlap = F, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) 
    }
    
    #c("Bacteria","Protists","Annelids","Archaea","Nematodes","Platyhelminthes","Arthropods","Fungi","Plants")
    
    rsquared[[i]][4] = cor(data.tmp1[[i]][,2],diversity)^2
  }
  
  #   ggsave(filename = "Stability_vs_size.pdf", do.call("arrangeGrob", c(plot.size, nrow=3)),
  #         height = 10*3/4, width = 10)
  #   
  #   ggsave(filename = "Stability_vs_diversity.pdf", do.call("arrangeGrob", c(plot.diversity, nrow=3)),
  #         height = 10*3/4, width = 10)
  # 
  #   ggsave(filename = "Stability_vs_diversity_restricted_archaea.pdf", do.call("arrangeGrob", c(plot.diversity1, nrow=3)),
  #         height = 10*3/4, width = 10)
  # 
  #   ggsave(filename = "Stability_vs_size_restricted_archaea.pdf", do.call("arrangeGrob", c(plot.size1, nrow=3)),
  #          height = 10*3/4, width = 10)
  
  # Short version without "Best K values", ES and SES:
  #   ii = 1
  #   new_plot_labels = c(rep("Spatial distribution",2),rep("MOTU composition",2))
  #   for (i in c(1,4,5,8))
  #   {
  #     #     plot.diversity1[[i]] = plot.diversity1[[i]] + ggtitle(paste(letters[ii],"-",new_plot_labels[ii]))
  #     #     plot.size1[[i]] = plot.size1[[i]] + ggtitle(paste(letters[ii],"-",new_plot_labels[ii]))
  #     ii = ii+1
  #   }
  
  #   ggsave(filename = "Stability_vs_diversity_SKL-pval_restricted_archaea.pdf", do.call("arrangeGrob", c(plot.diversity1[c(1,4,5,8)], nrow=2)),
  #        height = 10/2, width = 10/2)
  
  #   ggsave(filename = "Stability_vs_size_SKL-pval_restricted_archaea_spatial-only.pdf", do.call("arrangeGrob", c(plot.size1[c(1,4)], nrow=1)),
  #        height = 10/4, width = 10/2)
  
  #   ggsave(filename = "Stability_vs_size_pval_restricted_archaea_spatial-only.pdf", do.call("arrangeGrob", c(plot.size1[4], nrow=1)),
  #          height = 10/4, width = 10/4)
  #   
  #   ggsave(filename = "Stability_vs_diversity_pval_restricted_archaea_spatial-only.pdf", do.call("arrangeGrob", c(plot.diversity1[4], nrow=1)),
  #          height = 10/4, width = 10/4)
  
  #   ggsave(filename = "Stability_vs_size_nES_restricted_bacteria_spatial.pdf", do.call("arrangeGrob", c(plot.size1[13], nrow=1)),
  #          height = 10/4, width = 2*10/4)
  #   ggsave(filename = "Stability_vs_size_nES_restricted_bacteria_spatial.pdf", do.call("arrangeGrob", c(plot.size1[13], nrow=1)),
  #        height = 10/4, width = 10/4)
  #   rsquared[[13]][2]
  #   ggsave(filename = "Stability_vs_size_nES_restricted_bacteria_taxo.pdf", do.call("arrangeGrob", c(plot.size1[14], nrow=1)),
  #        height = 10/4, width = 10/4)
  #   rsquared[[14]][2]
  
  ggsave(filename = "Stability_vs_diversity_nES_restricted_bacteria_spatial-taxo_noplants_noarchaea1_noFungiITS.pdf", do.call("arrangeGrob", c(plot.diversity1[1], plot.diversity1[2], nrow=1)),
         height = 10/4, width = 2*10/4)
  #   ggsave(filename = "Stability_vs_diversity_nES_restricted_bacteria_spatial.pdf", do.call("arrangeGrob", c(plot.diversity1[13], nrow=1)),
  #        height = 10/4, width = 10/4)
  #   rsquared[[13]][4]  
  #   ggsave(filename = "Stability_vs_diversity_nES_restricted_bacteria_taxo.pdf", do.call("arrangeGrob", c(plot.diversity1[14], nrow=1)),
  #        height = 10/4, width = 10/4)
  #   rsquared[[14]][4]
  
  #   ggsave(filename = "Stability_vs_size_SKL-pval_restricted_archaea.pdf", do.call("arrangeGrob", c(plot.size1[c(1,4,5,8)], nrow=2)),
  #        height = 10/2, width = 10/2)
  
  #   pdf("SKL_vs_size.pdf",width=10,height=10/4)
  #   par(mfrow=c(1,4))
  #   
  #   plot(size,SKL_bij_samplewise,ylab="SKL",xlab="Size (log10)",main="One-to-one correspondence",type="p",pch=16)
  #   thigmophobe.labels(size, SKL_bij_samplewise, labels = taxo_names, cex=0.8, offset=0.5)
  #   plot(size,ES_bij_samplewise,ylab="ES",xlab="Size (log10)",main="One-to-one correspondence",type="p",pch=16)
  #   thigmophobe.labels(size, ES_bij_samplewise, labels = taxo_names, cex=0.8, offset=0.5)
  #   plot(size,SES_bij_samplewise,ylab="SES",xlab="Size (log10)",main="One-to-one correspondence",type="p",pch=16)
  #   thigmophobe.labels(size, SES_bij_samplewise, labels = taxo_names, cex=0.8, offset=0.5)
  #   plot(size,pval_bij_samplewise,ylab="pval",xlab="Size (log10)",main="One-to-one correspondence",type="p",pch=16)
  #   thigmophobe.labels(size, pval_bij_samplewise, labels = taxo_names, cex=0.8, offset=0.5)
  #   
  # #   plot(size,SKL_nobij_samplewise,ylab="SKL",xlab="Size (log10)",main="K best correspondences",type="p",pch=16)
  # #   thigmophobe.labels(size, SKL_nobij_samplewise, labels = taxo_names, cex=0.8, offset=0.5)
  #   #plot(size,SKL_bij_MOTUwise)
  #   dev.off()
  
  
  
}
####################################################################################################
####################################################################################################
if (taxo_tree)
{ 
  library(ggplot2)
  library(ggdendro)
  library(cluster)
  #   barcode_insert_list = c("Champignons_ITS","Bacteries_16S","Arthropodes_18S")
  #   barcode_labels_list = c("Fungi","Bacteria","Arthropods") 
  
  #   barcode_insert_list = c("Bacteries_16S","Protistes_18S","Champignons_ITS")
  #   barcode_labels_list = c("Bacteria","Protists","Fungi")  
  
  #   barcode_insert_list = c("Bacteries_16S","Protistes_18S","Champignons_ITS","Arthropodes_18S","Nematodes_18S","Platyhelminthes_18S","Annelides_18S")
  #   barcode_labels_list = c("Bacteria","Protists","Fungi","Arthropods","Nematodes","Platyhelminthes","Annelids") 
  
  barcode_insert_list = c("Bacteries_16S","Protistes_18S","Champignons_ITS","Arthropodes_18S","Nematodes_18S")
  barcode_labels_list = c("Bacteria","Protists","Fungi","Arthropods","Nematodes") 
  
  nb_topics = 3
  nb_real = 100
  em_tol = 10^-7
  var_tol = 10^-8
  
  occurrence = 1
  
  if (occurrence)
    occurrence_insert = "_occurrence"
  else
    occurrence_insert = ""
  
  taxotree.plot = list()
  barcode_index = 0
  upgma_Hellinger = list()   
  
  for (barcode_insert in barcode_insert_list)
  {    
    barcode_index = barcode_index+1  
    
    data_insert = "Donnees_PetitPlateau"
    filename_insert = "Rtopicmodels_LDA_VEM"
    local_prefix = "/Users/guilhemsommeria-klein/Desktop/These/" 
    local_dirname = paste(local_prefix,data_insert,"/",barcode_insert,"/",filename_insert,"/",sep="")
    local_subdirname = paste(local_dirname,filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best_keep",occurrence_insert,"/",sep="")
    setwd(paste0(local_subdirname,"/Comparison_in_best_real"))
    
    load("Topic_comparison_inBestReal.Rdata")  
    #   save(KL_topic_comparison_samplewise, KL_topic_comparison_MOTUwise, nES_topic_comparison_MOTUwise, 
    #        Corr_topic_comparison_samplewise, Corr_topic_comparison_MOTUwise,
    #        Hellinger_topic_comparison_MOTUwise, Jaccard_topic_comparison_MOTUwise, 
    #        Beta.sim_topic_comparison_MOTUwise, diversity, file = "Topic_comparison_inBestReal.Rdata")
    
    #     documents1 = documents[,rev(sort_normal_topic$ix)]
    #     topic_compo1 = topic_compo[,rev(sort_normal_topic$ix)]
    #     colnames(documents1) = assemblage_names_vect
    #     colnames(topic_compo1) = assemblage_names_vect
    #     Corr_topic_comparison_samplewise = cor(documents1)
    #     Corr_topic_comparison_MOTUwise = cor(topic_compo1)
    
    #     Hellinger_topic_comparison_MOTUwise = 1/sqrt(2)*dist(t(sqrt(topic_compo1)), method = "euclidean", diag = FALSE, upper = FALSE)
    #     KL_topic_compo1 = KL_topic_compo[,rev(sort_normal_topic$ix)]
    #     #Hellinger_topic_comparison_MOTUwise1 = 1/sqrt(2)*dist(t(sqrt(KL_topic_compo1)), method = "euclidean", diag = FALSE, upper = FALSE)
    #     topic_compo1_bin = topic_compo1
    #     topic_compo1_bin[topic_compo1 < 1/sum(data2m)] = 0
    #     #Hellinger_topic_comparison_MOTUwise2 = 1/sqrt(2)*dist(t(sqrt(topic_compo1_bin)), method = "euclidean", diag = FALSE, upper = FALSE)
    #     Jaccard_topic_comparison_MOTUwise = vegan::designdist(t(topic_compo1_bin), "(b+c)/(a+b+c)", abcd=TRUE)
    #     Beta.sim_topic_comparison_MOTUwise = vegan::designdist(t(topic_compo1_bin), "pmin(b,c)/(pmin(b,c)+a)", abcd=TRUE)
    
    # Code Antoine :
    #       upgma_Hellinger2 <- agnes(Hellinger_topic_comparison_MOTUwise2, diss =T, method = "average", keep.diss =F, keep.data =F)
    #       upgma_Hellinger1 <- agnes(Hellinger_topic_comparison_MOTUwise1, diss =T, method = "average", keep.diss =F, keep.data =F)
    # method = "average" : UPGMA
    upgma_Hellinger[[barcode_index]] = agnes(Hellinger_topic_comparison_MOTUwise, diss =T, method = "average", keep.diss =F, keep.data =F)
    #     upgma_Jaccard[[barcode_index]] = agnes(Jaccard_topic_comparison_MOTUwise, diss =T, method = "average", keep.diss =F, keep.data =F)
    #     upgma_Beta.sim[[barcode_index]] = agnes(Beta.sim_topic_comparison_MOTUwise, diss =T, method = "average", keep.diss =F, keep.data =F)
    
    #         upgma_data.frame = dendro_data(upgma_Hellinger[[barcode_index]])
    #     upgma_data.frame = dendro_data(upgma_Jaccard)
    #     upgma_data.frame = dendro_data(upgma_Beta.sim)
    
    taxotree.plot[[barcode_index]] = ggdendrogram(upgma_Hellinger[[barcode_index]], rotate = FALSE) +
      ggtitle(paste(letters[barcode_index],"-",barcode_labels_list[[barcode_index]])) +
      theme(plot.margin=unit(c(2,2,0,15),"mm"), plot.title = element_text(hjust=0,size=12), axis.text=element_text(size=10))
    
    #         taxotree.plot[[barcode_index]] = ggplot() +
    #         geom_segment(data = upgma_data.frame$segments, 
    #                    aes(x = x, y = y, xend = xend, yend = yend)) + 
    #         geom_text(data = upgma_data.frame$labels, 
    #                 aes(x = x, y = y, label = label), size = 4, vjust = 1) +
    #     #     geom_text(data = upgma_data.frame$leaf_labels, 
    #     #             aes(x = x, y = y, label = label), size = 3, vjust = 1) +
    #         theme_dendro() +
    #     # "top", "right", "bottom", "left"
    #     theme(plot.margin=unit(c(0,5,0,5),"mm"))
    
    #     geom_raster(data = spatial_topicmix_kriged_all_topics_colors, aes(x, y), 
    #               fill=rgb(spatial_topicmix_kriged_all_topics_colors[,3:5]/max(255,max(spatial_topicmix_kriged_all_topics_colors[,3:5]))), inherit.aes = F) +
    #geom_raster() +
    #     coord_equal() + theme_minimal() +
    #     labs(fill=barcode_labels_list[barcode_index]) + ggtitle(letters[barcode_index]) +
    #     #labs(fill=paste0("Assemblage ",k)) + ggtitle(letters[k]) +
    #     theme(legend.position="bottom", legend.text=element_text(size=13), 
    #         legend.title=element_text(size=16), axis.title=element_blank(), 
    #         axis.text = element_blank(),
    #         #               plot.title=element_text(hjust=0,size=16), plot.margin=unit(c(0,1,-2,2),"mm")) +
    #         plot.title=element_text(hjust=0,size=10), plot.margin=unit(c(2,-10,0,-10),"mm")) +
    #     guides(fill = guide_colorbar(barwidth = 8, barheight = 0.5, title.position="bottom")) +
    #     scale_y_continuous(limits=c(5,395), expand = c(0,0)) +
    #     scale_x_continuous(limits=c(5,295), expand = c(0,0)) + 
    #     geom_point(data = data.frame(coordPP,z.pred=rep(0,nrow(coordPP))), aes(x,y), color="black", size=0.1, alpha=0.3, inherit.aes = F)
  }
  
  setwd("/Users/guilhemsommeria-klein/Desktop/These/Manuscrits/Paper2/Figures/")
  
  #pdf("UPGMA_Helinger_bact-prot-fungi.pdf")
  #grid.arrange(grobs = taxotree.plot, ncol=3, height = 1, width = 3)
  #grid.arrange(grobs = taxotree.plot, ncol=3, height = 1, width = 3)
  # par(mfrow = c(1,3))
  # for (barcode_index in 1:3)
  #   plot(upgma_Hellinger[[barcode_index]], which.plots=2, ann=F)
  #dev.off()
  
  #   ggsave(filename = "UPGMA_Helinger_bact-prot-fungi-arth-nema-platy-anne.pdf",
  #         do.call("arrangeGrob",c(taxotree.plot,ncol=3)),width=9,height=8)
  
  ggsave(filename = "UPGMA_Helinger_bact-prot-fungi-arth-nema.pdf",
         do.call("arrangeGrob",c(taxotree.plot,ncol=3)),width=9,height=8/3*2)
  
  #     pdf("Hierarchical_clustering_of_topics_Hellinger.pdf")
  #     plot(upgma_Hellinger, which.plots=2, ann=F)
  #     title("Average Hellinger")
  #     dev.off()
  #     
  #     pdf("Hierarchical_clustering_of_topics_Jaccard.pdf")
  #     plot(upgma_Jaccard, which.plots=2, ann=F)
  #     title("Average Jaccard")
  #     dev.off()
  #     
  #     pdf("Hierarchical_clustering_of_topics_beta.sim.pdf")
  #     plot(upgma_Beta.sim, which.plots=2, ann=F)
  #     title("Average beta.sim")
  #     dev.off()
  
}
if (testdata)
{
  library(colorspace)
  library(ggplot2)
  library(gridExtra)
  library(scales)
  
  stability_plots = 0
  data_insert = "Donnees_PetitPlateau"
  filename_insert1 = "Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics0.02_1000sampledreads"
  # filename_insert1 = "Continuous-mixed_samples_nbtopics5_nbmotus1000_randomtopics_sampledreads"
  # filename_insert2 = "Rtopicmodels_LDA_VEM"
  filename_insert2 = ""
  # local_prefix = "/Users/guilhemsommeria-klein/Desktop/Serveur_Bioclust.data/" 
  local_prefix = "/Users/guilhemsommeria-klein/Desktop/These/" 
  barcode_insert = "Test_data"
  local_dirname = paste0(local_prefix,data_insert,"/",barcode_insert,"/",filename_insert1,"/")
  
  true_nb_topics = 5
  nb_topics = 5
  #   nb_real = 100
  #   em_tol = 10^-7
  #   var_tol = 10^-8
  nb_topics_range = 2:8
  
  # occurrence_vect = c(0,1,0,1,0,1)
  occurrence_vect = c(0,0,0,1,1,1)
  # mpar_vect = c(0,0,0,0,1,1)
  mpar_vect = c(0,1,1,0,1,1)
  plot.testdata = list()
  
  # setwd("/Users/guilhemsommeria-klein/Desktop/Manuscrits/ThesisPaper2_Topic_Modelling_PP/Figures/Submission2_analyses")
  # pdf("Fig_testdata_dirichlet_perplexity-AIC.pdf",height=9)
  #pdf("FigS1_testdata_AIC.pdf",height=9)
  # pdf("FigS2_testdata_AIC.pdf",height=6)
  # par(mfrow=c(3,2))
  # par(cex.lab=1,cex.main=1,cex.axis=1,lwd=1,font.main=1)
  # par(mar = c(bottom, left, top, right))
  #par(mar = c(5, 6, 1.5, 2) + 0.1)
  for (i in 1:6)
  {  
    occurrence = occurrence_vect[i]
    mpar = mpar_vect[i]  
    
    if (occurrence)
      occurrence_insert = "_occurrence"
    else
      occurrence_insert = ""
    
    setwd(local_dirname)
    load("true_documents.Rdata")
    true_documents = true_norm_documents
    
    if (i==1 || i==4)
    {
      local_subdirname = paste0(local_dirname,"Rtopicmodels_LDA_VEM/Rtopicmodels_LDA_VEM_nb_topics5_nb_real100_em_tol1e-07_var_tol1e-08_best_keep",occurrence_insert,"/")
      # if (i==1)
      #   par(mar = c(3, 4, 2.5, 2) + 0.1)
      # else if (i==2)
      #   par(mar = c(3, 2, 2.5, 2) + 0.1)
      load(paste0(local_subdirname,"Rtopicmodels_LDA_VEM_nb_topics5_nb_real100_em_tol1e-07_var_tol1e-08_best_keep",occurrence_insert,".Rdata"))
      Ordered_realizations = readRDS(paste0(local_subdirname,"Ordered_realizations.rds"))
      documents = Result[[Ordered_realizations$ix[1]]]@gamma
      load(paste0(local_subdirname,"sort_normal_topic.Rdata"))
      if (i==1)
      {
        # title = "Abundance"
        topic_order = rev(sort_normal_topic$ix)
      } else if (i==2)
      {
        # title = "Occurrence"
        topic_order = rev(sort_normal_topic$ix)
        topic_order[4] = rev(sort_normal_topic$ix)[5]
        topic_order[5] = rev(sort_normal_topic$ix)[4]
      }
      plot.testdata[[i]] = ggplot()
      for (k in 1:nb_topics)
      {
        if (k==1)
        {
          plot.testdata[[i]] = plot.testdata[[i]] + 
            geom_line(data = data.frame(y=documents[,topic_order[k]],x=seq_along(documents[,topic_order[k]])),aes(x,y),col=rainbow_hcl(nb_topics)[k]) +
            theme_bw() +
            theme(axis.title=element_text(size=11),
                  axis.text=element_text(size=11),
                  plot.title=element_text(hjust=0, size=14),
                  plot.margin=unit(c(0.1,1,0.1,1),"mm")) +
            ggtitle(LETTERS[i]) +
            # labs(x="Samples", y=ifelse(i==1,"Assemblage proportions\n in samples",""))
            labs(x="Samples", y=if (i == 1) "Assemb. proportions in samples\n - read-count data" else if (i == 4) "Assemb. proportions in samples\n - occurrence-transformed data")
        } else
          plot.testdata[[i]] = plot.testdata[[i]] + 
            geom_line(data = data.frame(y=documents[,topic_order[k]],x=seq_along(documents[,topic_order[k]])),aes(x,y),col=rainbow_hcl(nb_topics)[k])
      }
      for (k in 1:true_nb_topics)
        plot.testdata[[i]] = plot.testdata[[i]] + 
        geom_line(data = data.frame(y=true_documents[,k],x=seq_along(true_documents[,k])),aes(x,y),col="black",linetype="dashed")
      
      # for (k in 1:nb_topics)
      # {
      #   if (k==1)
      #   {
      #     plot(documents[,topic_order[k]],type="l",col=rainbow_hcl(nb_topics)[k],ann=F,tck=-0.02)
      #     if (i==1)
      #       title(ylab="Assemblage proportions\n in samples",line=2.5)
      #     title(paste(letters[i],"-",title),line=1,adj=0)
      #     title(xlab="Samples",line=2)
      #     #           axis(side = 1, padj = -1, tck=-0.02)
      #     #           axis(side = 2, padj = -1, tck=-0.02)
      #   } else
      #     lines(documents[,topic_order[k]],type="l",col=rainbow_hcl(nb_topics)[k])
      # }
      # for (k in 1:true_nb_topics)
      #   lines(true_documents[,k],type="l",col="black",lty=2)
    } else if (stability_plots == 1 && (i==3 || i==4))
    {
      local_subdirname = paste0(local_dirname,"Rtopicmodels_LDA_VEM/Rtopicmodels_LDA_VEM_nb_topics5_nb_real100_em_tol1e-07_var_tol1e-08_best_keep",occurrence_insert,"/")
      if (i==3)
        par(mar = c(5, 6, 1.5, 2) + 0.1)
      else if (i==4)
        par(mar = c(5, 2, 1.5, 6) + 0.1)
      setwd(local_subdirname)
      Ordered_realizations = readRDS("Ordered_realizations.rds")
      setwd(paste0(local_subdirname,"Realization_comparison_samplewise_bijective_correspondence/"))
      load("SKL-correlation_allRealPairs_samplewise_bijective_correspondence.Rdata")
      nES = t(DKL100_allRealPairs)[lower.tri(DKL100_allRealPairs,diag=F)]/KL_allRealPairs_w_rndzations[lower.tri(KL_allRealPairs_w_rndzations,diag=F)]  
      # par(mar = c(bottom, left, top, right))
      #       par(mar = c(4, 6, 4, 3) + 0.1)
      plot(Ordered_realizations$x[1]-Ordered_realizations$x[2:100],nES[1:(100-1)],type="p",pch=20,xlab="",ylab="")
      if (i==3)
        title(ylab = "Spatial similarity\n to best realization")
      if (i == 3)
      {
        #title = "Abundance"
        llh_dif = Ordered_realizations$x[1]-Ordered_realizations$x[2:100]
        fit_ab = lm(nES[1:(100-1)] ~ llh_dif)
        lines(llh_dif,llh_dif*fit_ab$coefficient[2]+fit_ab$coefficient[1]*rep(1,99),type="l",col="blue",lty=2)
      } else if (i == 4)
      {
        #title = "Occurrence"
        llh_dif = Ordered_realizations$x[1]-Ordered_realizations$x[2:100]
        fit_oc = lm(nES[1:(100-1)] ~ llh_dif)
        lines(llh_dif,llh_dif*fit_oc$coefficient[2]+fit_oc$coefficient[1]*rep(1,99),type="l",col="blue",lty=2)
      }
      title(xlab="Llh difference\n with best realization",line=4)
      title(letters[i-2],line=1,adj=0)
    } else if (stability_plots == 0 && (i==2 || i==5))
    {
      local_subdirname = paste0(local_dirname,"Rtopicmodels_LDA_VEM_alpha0.1_nb_topics2-8_post-predictive-cross-valid_fold_size10_em_tol1e-06_var_tol1e-08",occurrence_insert,"/")
      # if (i==3)
      #   par(mar = c(4, 6, 1.5, 2) + 0.1)
      # else if (i==4)
      #   par(mar = c(4, 2, 1.5, 6) + 0.1)
      load(paste0(local_subdirname,"perplexity.Rdata"))
      
      plot.testdata[[i]] = ggplot(data = data.frame(K = unlist(lapply(nb_topics_range,rep,nrow(perplexity_mpar))), llh = as.vector(perplexity_mpar))) +
        geom_boxplot(aes(x=as.factor(K),y=llh), size = 0.4, outlier.size = 0.4, outlier.color = "darkgrey") +
        # scale_x_discrete(breaks = seq(from = nb_topics_range[1], to = nb_topics_range[length(nb_topics_range)], by = 3)) +
        theme_bw() +
        theme(axis.title=element_text(size=11),
              axis.text=element_text(size=11),
              plot.title=element_text(hjust=0, size=14),
              plot.margin=unit(c(0.1,1,0.1,1),"mm")) +
        ggtitle(LETTERS[i]) +
        # labs(x="Number K of assemblages", y=ifelse(i==3,"Cross-validation perplexity",""))
        labs(x="Number K of assemblages", y="Cross-validation perplexity")
      
      # plot(nb_topics_range,AIC3[1,],ann=F,type="p",pch=20,col="black",ylim=range(AIC3))
      # for (j in 1:nb_real)
      #   lines(nb_topics_range,AIC3[j,],pch=20,col="black",lty=0,type="p")
      # 
      # title(letters[i-2],line=1,adj=0)
      # title(xlab="Number K of assemblages")
      # if (i==5)
      #   title(ylab="Cross-validation perplexity")
    } else if (i==3 || i==6)
    {
      local_subdirname = paste0(local_dirname,"Rtopicmodels_LDA_VEM_alpha0.1_nb_topics2-8_elbow-AIC_nb_real50_em_tol1e-06_var_tol1e-08",occurrence_insert,"/")
      # if (i==5)
      #   par(mar = c(4, 6, 1.5, 2) + 0.1)
      # else if (i==6)
      #   par(mar = c(4, 2, 1.5, 6) + 0.1)
      load(paste0(local_subdirname,"AIC-llh.Rdata"))
      # nb_real = 50
      # AIC3 = matrix(nrow=nb_real,ncol=7,data=0)
      # #AICc = matrix(nrow=nb_real,ncol=7,data=0)
      # for (par_index in 1:7)
      # {
      #   nb_topics = nb_topics_range[par_index] 
      #   for (j in 1:nb_real)
      #   {
      #     if (j==1)
      #     {
      #       Result = Result_mpar[[(par_index-1)*nb_real+j]]
      #       nb_terms = Result@wordassignments$ncol
      #       nb_doc = Result@wordassignments$nrow
      #       #       nb_words = Result[[1]]@n
      #       AIC3[j,par_index] = 2*(nb_topics*(nb_terms-1) + 1 - sum(Result@loglikelihood))
      #       #AICc[j,par_index] = 2*((nb_topics*(nb_terms-1) + 1)*(1+1/nb_doc) - sum(Result@loglikelihood)) 
      #     }
      #     else if (j>1)
      #     {
      #       Result = c(Result,Result_mpar[[(par_index-1)*nb_real+j]])
      #       AIC3[j,par_index] = 2*(nb_topics*(nb_terms-1) + 1 - sum(Result[[j]]@loglikelihood))
      #       #AICc[j,par_index] = 2*((nb_topics*(nb_terms-1) + 1)*(1+1/nb_doc) - sum(Result[[j]]@loglikelihood)) 
      #     }
      #   }
      # }
      #par(mar=c(5.1,4.1,4.1,2.1))
      #       par(mar=c(5.1,5.1,4.1,2.1))
      
      # llh = LLH_final0[1,]
      # sd_llh = LLH_final0[2,]
      
      # llh_allpoints = as.vector(LLH_final)
      nb_real = 50
      nb_topics_range_allpoints = unlist(lapply(nb_topics_range,rep,nb_real))
      
      plot.testdata[[i]] = ggplot(data=data.frame(K=nb_topics_range_allpoints,llh=as.vector(AIC))) +
        geom_boxplot(aes(x=as.factor(K),y=llh), size = 0.4, outlier.size = 0.4, outlier.color = "darkgrey") +
        #scale_x_discrete(breaks = seq(from = nb_topics_range[1], to = nb_topics_range[length(nb_topics_range)], by = 3)) +
        scale_y_continuous(labels = scales::scientific) +
        theme_bw() +
        # ggtitle(paste(taxo_names[i_taxon],"-",as.vector(diversity)[i_taxon],"OTUs")) +
        ggtitle(LETTERS[i]) +
        theme(axis.title=element_text(size=11),
              axis.text.x=element_text(size=11),
              axis.text.y=element_text(size=8),
              plot.title=element_text(hjust=0, size=14),
              plot.margin=unit(c(0.1,1,0.1,1),"mm")) +
        labs(x="Number K of assemblages", y="AIC")
        # labs(x="Number K of assemblages", y=ifelse(i==5,"AIC",""))
      
      # plot(nb_topics_range,AIC3[1,],ann=F,type="p",pch=20,col="black",ylim=range(AIC3))
      # for (j in 1:nb_real)
      #   lines(nb_topics_range,AIC3[j,],pch=20,col="black",lty=0,type="p")
      # 
      # title(letters[i-2],line=1,adj=0)
      # title(xlab="Number K of assemblages")
      # if (i==5)
      #   title(ylab="AIC")
    }
  }

  #####################
  library(maxmatching)
  library(vegan)
  source("/Users/guilhemsommeria-klein/Desktop/Code/R/LDA_project/LDA_topic_correspondence_fun.R")
  
  # data.folder_name = "/Users/guilhemsommeria-klein/Desktop/Serveur_Bioclust.data/Donnees_PetitPlateau/Test_data/"
  data.folder_name = "/Users/guilhemsommeria-klein/Desktop/These/Donnees_PetitPlateau/Test_data/"
  
  test_taxo_vect = c("Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics0.001_1000sampledreads",
                     "Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics0.005_1000sampledreads",
                     "Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics0.01_1000sampledreads",
                     "Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics0.02_1000sampledreads",
                     "Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics0.05_1000sampledreads",
                     "Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics0.1_1000sampledreads",
                     "Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics0.5_1000sampledreads",
                     "Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics1_1000sampledreads",
                     "Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics5_1000sampledreads")
                     # "Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics10_1000sampledreads")
  
  Rsquared = list()
  mean_nb_OTUs = vector(length = length(test_taxo_vect), mode = "numeric")
  tot_nb_OTUs = vector(length = length(test_taxo_vect), mode = "numeric")
  Sorensen_samples = vector(length = length(test_taxo_vect), mode = "numeric")
  test_taxo_names = c("0.001","0.005","0.01","0.02","0.05","0.1","0.5","1","5")
  names(mean_nb_OTUs) = names(tot_nb_OTUs) = names(Sorensen_samples) = test_taxo_names
  mean_rank_abundance = matrix(nrow = 1000, ncol = length(test_taxo_vect), dimnames = list(NULL,test_taxo_names), data = 0)
  mean_S100 = I100 = Sorensen_assemblages = nb_OTUs_per_topic = matrix(nrow = length(test_taxo_vect), ncol = 2, dimnames = list(test_taxo_names,c("Occurrence","Read-count")), data = 0)
  for (case in 1:2)
  {
    Rsquared[[case]] = vector(length = length(test_taxo_vect), mode = "numeric")
    names(Rsquared[[case]]) = test_taxo_names
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
      
      stability_samplewise = readRDS(paste0(taxon.folder_name,result.folder,"/Stability_assessment_samplewise_maxmatching/Stability_samplewise.rds"))
      mean_S100[which(test_taxo_vect == taxon),case] = stability_samplewise[1,4]
      I100[which(test_taxo_vect == taxon),case] = stability_samplewise[1,5]
      
      # logbeta = Result[[Ordered_realizations$ix[1]]]@beta 
      # # topic_compo : proportion of each MOTU in each topic (sums to 1 over MOTUs for each topic)
      # topic_compo = exp(t(logbeta))
      # topic_compo = topic_compo[,sort_normal_topic$ix]
      # nb.OTUS[[case]][[i_taxon]] = vector(length = nb_topics, mode = "numeric")
      # names(nb.OTUS[[case]][[i_taxon]]) = paste("Assemblage",1:nb_topics)
      # for (k in 1:nb_topics)
      #   nb.OTUS[[case]][[i_taxon]][k] = length(which(topic_compo[,k] > 1/sum_data2m))
      
      if (case == 1)
      {
        load(paste0(taxon.folder_name,"/data2m_testdata.Rdata"))
        nb_OTUs = vector(length = ncol(data2m), mode = "numeric")
        for (i in 1:ncol(data2m))
          nb_OTUs[i] = length(which(data2m[,i]>0))
        mean_nb_OTUs[which(test_taxo_vect == taxon)] = mean(nb_OTUs)
        tot_nb_OTUs[which(test_taxo_vect == taxon)] = length(which(rowSums(data2m) != 0))

        taxon.folder_name1 = paste0(data.folder_name,taxon,"1")
        load(paste0(taxon.folder_name1,"/true_topic_compo.Rdata"))
        load(paste0(taxon.folder_name1,"/data2m_testdata.Rdata"))
        data2m0 = data2m
        data2m[data2m>0] = 1
        effective_OTUs_ab = true_topic_compo > 1/sum(data2m0)
        effective_OTUs_oc = true_topic_compo > 1/sum(data2m)
        nb_OTUs_ab = nb_OTUs_oc = vector(length = 5, mode = "numeric")
        for (j in 1:5)
        {
          nb_OTUs_ab[j] = length(which(effective_OTUs_ab[,j]))
          nb_OTUs_oc[j] = length(which(effective_OTUs_oc[,j]))
        }
        
        Sorensen_assemblages[which(test_taxo_vect == taxon),] = c(mean(1 - vegan::designdist(t(effective_OTUs_oc), "(b+c)/(2*a+b+c)", abcd=TRUE)),
                                                                  mean(1 - vegan::designdist(t(effective_OTUs_ab), "(b+c)/(2*a+b+c)", abcd=TRUE)))
        Sorensen_samples[which(test_taxo_vect == taxon)] =  mean(1 - vegan::designdist(t(data2m), "(b+c)/(2*a+b+c)", abcd=TRUE))
        
        nb_OTUs_per_topic[which(test_taxo_vect == taxon),2] = mean(nb_OTUs_ab)
        nb_OTUs_per_topic[which(test_taxo_vect == taxon),1] = mean(nb_OTUs_oc)
        mean_rank_abundance[,which(test_taxo_vect == taxon)] = sort(rowMeans(apply(true_topic_compo,2,sort,decreasing=T)),decreasing=T)
        mean_rank_abundance[mean_rank_abundance[,which(test_taxo_vect == taxon)] < 1/sum(data2m0),which(test_taxo_vect == taxon)] = NA
      }
    }
  }
  
  ###############
  plot.testdata[[7]] = ggplot(data = data.frame(y = c(Rsquared[[2]],Rsquared[[1]]), x = c(nb_OTUs_per_topic[,2],nb_OTUs_per_topic[,1]),
                                      shape = c(rep("Read count",length(test_taxo_vect)),rep("Occurrence",length(test_taxo_vect))))) +
    geom_point(aes(x,y,shape=shape), col = rep(c("red",rep("black",2),"blue",rep("black",4),"forestgreen"),2), size = 2) +
    geom_line(data = data.frame(y = Rsquared[[1]][1:(length(test_taxo_vect)-1)], x = nb_OTUs_per_topic[1:(length(test_taxo_vect)-1),1]), aes(x,y), linetype = "dashed", size = 0.5) +
    geom_line(data = data.frame(y = Rsquared[[2]][1:(length(test_taxo_vect)-1)], x = nb_OTUs_per_topic[1:(length(test_taxo_vect)-1),2]), aes(x,y), linetype = "dashed", size = 0.5) +
    # geom_vline(xintercept = nb_OTUs_per_topic_oc[4], linetype = "longdash", size = 0.25) +
    # scale_shape_manual(values=c(16, 15)) +
    # geom_point(aes(x=mean.sample.size,y=occ.Rsquared,fill=occ.fill), col = "black", size = 2, shape = 15) +
    # geom_line(aes(x=mean.sample.size,y=occ.Rsquared), linetype = "dashed", size = 0.5) +
    theme_bw() +
    ggtitle(LETTERS[7]) +
    #theme(axis.title.x = element_text(hjust = 0.8, size=11),
    theme(axis.title.x = element_text(size=11),
          axis.title.y = element_text(size=11),
          axis.text = element_text(size=11),
          plot.title=element_text(hjust=0, size=14),
          legend.position = c(0.6, 0.2),
          # legend.position = "right",
          legend.text=element_text(size=11), 
          legend.title=element_blank(),
          legend.background = element_rect(colour = "white"),
          legend.box.background = element_rect(colour = "black"),
          plot.margin=unit(c(0.1,2,1,1),"mm")) +
    labs(x="Mean nb. of OTUs per assemblage\n out of 1,000 simulated OTUs",
    # labs(x="Mean nb. of OTUs per sample,\n out of 1,000 simulated OTUs", 
    # labs(x="Mean nb. of OTUs per sample\n wrt. the total nb. of simulated OTUs", 
         y=bquote("Model's "*R^2))
  
  plot.testdata[[8]] = ggplot(data = data.frame(y = c(mean_S100[,1],mean_S100[,2]),
  # plot.testdata[[8]] = ggplot(data = data.frame(y = c(I100[,1],I100[,2]), 
                                                x = c(nb_OTUs_per_topic[,1],nb_OTUs_per_topic[,2]),
                                                shape = c(rep("Occurrence",length(test_taxo_vect)),rep("Read count",length(test_taxo_vect))))) +
    geom_point(aes(x,y,shape=shape), col = rep(c("red",rep("black",2),"blue",rep("black",4),"forestgreen"),2), size = 2) +
    geom_line(data = data.frame(y = mean_S100[1:(length(test_taxo_vect)-1),1], x = nb_OTUs_per_topic[1:(length(test_taxo_vect)-1),1]), aes(x,y), linetype = "dashed", size = 0.5) +
    geom_line(data = data.frame(y = mean_S100[1:(length(test_taxo_vect)-1),2], x = nb_OTUs_per_topic[1:(length(test_taxo_vect)-1),2]), aes(x,y), linetype = "dashed", size = 0.5) +
    theme_bw() +
    ggtitle(LETTERS[8]) +
    #theme(axis.title.x = element_text(hjust = 0.8, size=11),
    theme(axis.title.x = element_text(size=11),
          axis.title.y = element_text(size=11),
          axis.text = element_text(size=11),
          plot.title=element_text(hjust=0, size=14),
          legend.position = c(0.6, 0.2),
          # legend.position = "right",
          legend.text=element_text(size=11),
          legend.title=element_blank(),
          legend.background = element_rect(colour = "white"),
          legend.box.background = element_rect(colour = "black"),
          plot.margin=unit(c(0.1,2,1,1),"mm")) +
    labs(y=bquote("<"*S*">"[100]),
    # labs(y=bquote(I[100]),
         x="Mean nb. of OTUs per assemblage\n out of 1,000 simulated OTUs")
  
  plot.testdata[[9]] = ggplot(data = data.frame(mean_rank_abundance)) +
    theme_bw() +
    ggtitle(LETTERS[9]) +
    scale_y_log10() +
    theme(axis.title = element_text(size=11),
          axis.text.x = element_text(size=11),
          axis.text.y = element_text(size=8),
          plot.title=element_text(hjust=0, size=14),
          # legend.position = c(0.6, 0.2),
          # legend.position = "right",
          # legend.text=element_text(size=11), 
          # legend.title=element_blank(),
          # legend.background = element_rect(colour = "white"),
          # legend.box.background = element_rect(colour = "black"),
          plot.margin=unit(c(0.1,2,1,1),"mm")) +
    labs(y="OTU proportion in assemblage",x="OTU abundance rank in assemblage")
  for (j in 1:ncol(mean_rank_abundance))
    plot.testdata[[9]] = plot.testdata[[9]] + 
      geom_line(data = data.frame(x = (1:1000)[!is.na(mean_rank_abundance[,j])],
                                  y = mean_rank_abundance[!is.na(mean_rank_abundance[,j]),j]), 
                aes(x,y), size = 0.25, col = c("red",rep("black",2),"blue",rep("black",4),"forestgreen")[j])
  
  # bquote('Assimilation ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')
  
  # spl = split(plot.posteriorLlh.errbar.prevalence[c(T,selected_groups) & !is.na(plot.posteriorLlh.errbar.prevalence)], 
  #             (seq_along(plot.posteriorLlh.errbar.prevalence[c(T,selected_groups) & !is.na(plot.posteriorLlh.errbar.prevalence)])-1) %/% 20)
  ppl = marrangeGrob(grobs = plot.testdata, nrow = 3, ncol = 3, layout_matrix = matrix(data=1:9,nrow=3,byrow=T), top = NULL)
  figure_folder = "/Users/guilhemsommeria-klein/Desktop/Manuscrits/ThesisPaper2_Topic_Modelling_PP/Figures/Submission2_analyses/"
  # pdf(paste0(figure_folder,"Fig_testdata_dirichlet_perplexity-AIC.pdf"),height=9)
  pdf(paste0(figure_folder,"Fig_testdata_random_perplexity-AIC_withParamComparison1.pdf"),height=9,width=9)
  print(ppl)
  dev.off()
  
}
if (rank_abundance)
{
  #library("utils")
  
  barcode_insert_list = c("Bacteries_16S","Protistes_18S","Champignons_ITS")
  barcode_labels_list = c("Bacteria","Protists","Fungi")
  short_barcode_insert_list = c("16Sbact","","ITSfungi")
  
  data_insert = "Donnees_PetitPlateau"
  short_data_insert = "PP"
  
  local_prefix = "/Users/guilhemsommeria-klein/Desktop/These/" 
  
  Ranked_abundances = list()
  barcode_index = 0
  for (barcode_insert in barcode_insert_list)
  {    
    barcode_index = barcode_index+1  
    
    short_barcode_insert = short_barcode_insert_list[barcode_index]
    
    local_dirname = paste0(local_prefix,data_insert,"/",barcode_insert,"/")
    setwd(local_dirname)
    
    if (barcode_insert == "Bacteries_16S" || barcode_insert == "Champignons_ITS")
    {
      data2 = read.table(paste(short_data_insert,"_",short_barcode_insert,"_sequences_counts_norep.txt",sep=""), colClasses="vector", sep=" ")
      data2m = apply(data2[-1,],2,as.numeric)
    } else if (barcode_insert == "Protistes_18S")
    {
      load("data2m_euka_protists.Rdata")
      data2m = data2m_euka_protists
    }
    
    Ranked_abundances[[barcode_index]] = sort.int(rowSums(data2m),decreasing=T)
  }
  
  setwd("/Users/guilhemsommeria-klein/Desktop/These/Manuscrits/Paper2/Figures/")
  pdf("Rank-abundance_bact-prot-fungi.pdf",width=9,height=3)
  par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
  #par(mar=c(5.1,4.1,4.1,2.1)
  #bottom, left, top and right
  par(mar=c(2.1,3.1,2.1,1.5))
  par(mfrow=c(1,3))
  for (i in 1:barcode_index)
  {
    plot(seq(1,length(Ranked_abundances[[i]])),log10(Ranked_abundances[[i]]),type="l",ann=F)
    #plot(abundances_taxo2[,1],type="p",ann=F,yaxt="n")
    #axis(2, ylim=range(abundances_taxo2[,1]), col='black')
    #par(new=T)
    #plot(log(abundances_taxo2[,1])/log(10),type="p",col="blue",ann=F,yaxt="n",xaxt="n")
    #axis(4, ylim=range(log(abundances_taxo2[,1])), col='blue')
    #mtext("Read number (log_10)",side=4,line=3,col="blue",cex=1.5) 
    #title(xlab="Sequences ranked by abundances",ylab="Read number")
    #legend(x="topright",legend=legend_topic[1:nb_topics],text.col=col_topic[1:nb_topics],inset=0.1)
    #title("Total number of reads per sequence")
  }
  dev.off()
}
if (mock_stability)
{ 
  interval = seq(1,100)
  llh.diff = vector(length=100,mode="numeric")
  llh.diff[1] = 1
  for (i in 1:99)
    llh.diff[i+1] = interval[i]+llh.diff[i]
  llh.diff2 = seq(1,5000,50)
  #mock.data.convergence = data.frame(x=llh.diff,y=(-llh.diff*0.5/llh.diff[length(llh.diff)]+1)*rnorm(100,1,llh.diff/20000))
  #mock.data.noConvergence = data.frame(x=llh.diff2,y=rnorm(100,0.5,0.15))
  load("/Users/guilhemsommeria-klein/Desktop/Manuscrits/ThesisPaper2_Topic_Modelling_PP/Figures/mock.data_fig2.Rdata")
  
  plot.mock.data = list()
  
  plot.mock.data[[1]] = qplot(x, y, data = mock.data.convergence, geom="point") +
    theme_bw() + ggtitle(LETTERS[1]) +
    #ylim(min(0,min(data.dskln.llh$y)),max(data.dskln.llh$y))
    ylim(0,1) +
    theme(axis.title=element_text(size=22), 
          plot.title=element_text(hjust=0, size=22),
          # t, r, b, l
          plot.margin=unit(c(1,7,1,1),"mm"),
          axis.text=element_text(size=18),
          axis.line=element_line(size=0.8)) +
    geom_smooth(method='lm',se=F,linetype="dashed",size=0.8,col="black") +
    labs(x="\nLlh difference\n with best realization", y="Similarity S\n to best realization\n")
  
  plot.mock.data[[2]] = qplot(x, y, data = mock.data.noConvergence, geom="point") +
    theme_bw() + ggtitle(LETTERS[2]) +
    #ylim(min(0,min(data.dskln.llh$y)),max(data.dskln.llh$y))
    ylim(0,1) +
    theme(axis.title=element_text(size=22), 
          plot.title=element_text(hjust=0, size=22),
          plot.margin=unit(c(1,7,1,1),"mm"),
          axis.text=element_text(size=18),
          axis.line=element_line(size=0.8)) +
    geom_smooth(method='lm',se=F,linetype="dashed",size=0.8,col="black") +
    labs(x="\nLog-likelihood difference\n with best realization", y="")
  
  #save(mock.data.convergence,mock.data.noConvergence,file="mock.data_fig2.Rdata")
  
  ggsave(filename = "/Users/guilhemsommeria-klein/Desktop/Manuscrits/ThesisPaper2_Topic_Modelling_PP/Figures/Submission2_analyses/Figure2_mock_data4.pdf",
         do.call("arrangeGrob", c(plot.mock.data, nrow=1)),
         height = 12/2, width = 12)
}

