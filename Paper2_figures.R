kriged_topics = 1
stability = 0
size = 0
taxo_tree = 0

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
  
  color.pal = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  #barcode_insert_list = c("Bacteries_16S","Archees_16S","Protistes_18S","Arthropodes_18S","Nematodes_18S","Annelides_18S","Champignons_ITS","Plantes_GH")
  #   barcode_insert_list = c("Bacteries_16S","Archees_16S","Protistes_18S","Champignons_ITS","Plantes_GH")
  #   barcode_labels_list = c("Bacteria","Archaea","Protists","Fungi","Plants")
  
  #   barcode_insert_list = c("Arthropodes_18S","Nematodes_18S","Annelides_18S","Platyhelminthes_18S")
  #   barcode_labels_list = c("Arthropods","Nematodes","Annelids","Platyhelminthes")
  
  #     barcode_insert_list = "Arthropodes_18S"
  #     barcode_labels_list = "Arthropods"
  
#       barcode_insert_list = "Bacteries_16S"
#       barcode_labels_list = "Bacteria"
  
#   barcode_insert_list = c("Bacteries_16S","Archees_16S","Protistes_18S","Champignons_ITS","Plantes_GH","Arthropodes_18S","Nematodes_18S","PLatyhelminthes_18S")
#   barcode_labels_list = c("Bacteria 16S","Archaea 16S","Protists 18S","Fungi ITS","Plants trnL","Arthropods 18S","Nematodes 18S","Platyhelminthes 18S")
  
  # All barcodes :
#   barcode_insert_list = c("Bacteries_16S","Archees_16S","Protistes_18S","Champignons_ITS","Glomeromycetes_ITS","Basidiomycetes_ITS","Ascomycetes_ITS","Plantes_GH","Arthropodes_18S","Nematodes_18S","PLatyhelminthes_18S","Annelides_18S")
#   barcode_labels_list = c("Bacteria 16S","Archaea 16S","Protists 18S","Fungi ITS","Glomeromycota ITS","Basidiomycota ITS","Ascomycota ITS","Plants trnL","Arthropods 18S","Nematodes 18S","Platyhelminthes 18S","Annelids 18S")
  
  # All barcodes except Nema :
#   barcode_insert_list = c("Bacteries_16S","Archees_16S","Protistes_18S","Champignons_ITS","Glomeromycetes_ITS","Basidiomycetes_ITS","Ascomycetes_ITS","Plantes_GH","Arthropodes_18S","PLatyhelminthes_18S","Annelides_18S")
#   barcode_labels_list = c("Bacteria 16S","Archaea 16S","Protists 18S","Fungi ITS","Glomeromycota ITS","Basidiomycota ITS","Ascomycota ITS","Plants trnL","Arthropods 18S","Platyhelminthes 18S","Annelids 18S")

#   barcode_insert_list = c("Bacteries_16S","Protistes_18S","Champignons_ITS")
#   barcode_labels_list = c("Bacteria","Protists","Fungi")

#   barcode_insert_list = rep(c("Bacteries_16S","Protistes_18S","Champignons_ITS"),3)
#   barcode_labels_list = rep(c("Bacteria","Protists","Fungi"),3)

barcode_insert_list = c(rep(c("Archees_16S","Arthropodes_18S","Nematodes_18S"),2),rep(c("Platyhelminthes_18S","Annelides_18S"),2))
barcode_labels_list = c(rep(c("Archaea","Arthropods","Nematodes"),2),rep(c("Platyhelminthes","Annelids"),2))
  
  #   barcode_insert_list = c("Champignons_ITS","Plantes_GH")
  #   barcode_labels_list = c("Fungi","Plants")
  
  #   barcode_insert_list = c("Archees_16S","Protistes_18S")
  #   barcode_labels_list = c("Archaea","Protists")
  
  #   barcode_insert_list = c("Champignons_ITS","Plantes_GH","Arthropodes_18S")
  #   barcode_labels_list = c("Fungi ITS","Plants trnL","Arthropods 18S")
  
  #occurrence_vect = c(0,0,0,1,1,1)
  #occurrence_vect = c(1,1,1,1,1,1,1,1,1)
  occurrence_vect = c(1,1,1,0,0,0,1,1,0,0)
#    occurrence_vect = 1
  # Number of sites an OTU must ocuppy to be kept
  # (if nb_occupied_sites_threshold = 1, all OTUs with non-zero abundance are kept)
  #nb_occupied_sites_threshold_vect = c(1,1,1,2,2,2,1,1,1)
  nb_occupied_sites_threshold_vect = rep(1,10)
  # Removing OTUs with less reads than the number of sites
  #no_rare_vect = c(0,0,0,0,0,0,1,1,1)
  no_rare_vect = rep(0,10)

  nb_topics = 3
  nb_real = 100
  em_tol = 10^-7
  var_tol = 10^-8
  
  if (vertical)
    horizVert_insert = "_vertical"
  else if (horizontal)
    horizVert_insert = "_horizontal"

  if (discrete)
    discrete_insert = "_discrete"
  else
    discrete_insert = ""
  
########
setwd(paste("/Users/guilhemsommeria-klein/Desktop/These/",data_insert,"/Lidar/",sep=""))
r_topo <- raster("r_topol_0.asc")
r_topo_transposed_flipped = flip(t(r_topo),'y')
r_slope <- raster("r_slopel_0.asc")
r_slope_transposed_flipped = flip(t(r_slope),'y')
r_wetness <- raster("r_wetnessl_0.asc")
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

Assemblage_labels_list = c("Terra firme","Hydromorphic","Exposed rock")

#map_labels = c("Topography","Inverse topography","Slope")
map_labels = c("Topography (m a.s.l.)","Wetness","Slope")

##########
# if (occurrence && (any(barcode_insert_list == "Plantes_GH")))
# {
#   setwd(paste("/Users/guilhemsommeria-klein/Desktop/These/",data_insert,"/Chemistery/",sep=""))
#   # produces a dataframe with columns as vectors instead of factors (otherwise one cannot direclty apply "as.numeric" to the values)
#   data_abiotic = read.table("chemistry_pred_CRIJ.txt",sep=" ",colClasses="vector")
#   data_abiotic = apply(data_abiotic,2,as.numeric)
#   
#   abiotic_PCA = dudi.pca(as.data.frame(data_abiotic), row.w = rep(1, nrow(data_abiotic))/nrow(data_abiotic), 
#                          col.w = rep(1, ncol(data_abiotic)), center = TRUE, scale = TRUE, 
#                          scannf = F, nf = ncol(data_abiotic)) 
#   
#   newdata = expand.grid(x=seq(0,300,2), y=seq(0,400,2))
#   coordPP = expand.grid(x=seq(10,290,10), y=seq(10,390,10))
#   spatial_chemi_kriged = list()
#   for (k in 1:nb_topics)
#   {
#     data_chemi = data.frame(x=coordPP$x,y=coordPP$y,z=abiotic_PCA$li[,k])
#     mod = gstat(id = "z", formula = z~1, locations = ~x+y, data=data_chemi, model=vgm(10, "Exp", 20))
#     spatial_chemi_kriged[[k]] = predict(mod, newdata, na.action=na.omit)
#   }
# }

###########

  tmp.plot = list()
  lidar.plot = list()
#   chemi.plot = list()
  barcode_index = 0
  for (barcode_insert in barcode_insert_list)
  {
    barcode_index = barcode_index+1
    
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
    
    data_insert = "Donnees_PetitPlateau"
    filename_insert = "Rtopicmodels_LDA_VEM"
    local_prefix = "/Users/guilhemsommeria-klein/Desktop/These/" 
    local_dirname = paste(local_prefix,data_insert,"/",barcode_insert,"/",filename_insert,"/",sep="")
    local_subdirname = paste(local_dirname,filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best_keep",occurrence_insert,remove_single_sites_insert,no_rare_insert,"/",sep="")
    if (occurrence_vect[barcode_index] && (nb_occupied_sites_threshold_vect[barcode_index] == 1) && !no_rare_vect[barcode_index] && (barcode_insert == "Platyhelminthes_18S"))
      subsubdirname = paste(local_subdirname,"ordered_realization_number_5/",sep="")
    else 
      subsubdirname = paste(local_subdirname,"ordered_realization_number_1/",sep="")  
    subsubsubdirname = paste(subsubdirname,"/topics_site_repartition_info/",sep="")  
    
    #newdata = expand.grid(x=seq(0,300,2), y=seq(0,400,2))
    coordPP = expand.grid(x=seq(10,290,10), y=seq(10,390,10))
    # coordPP = coordPP[-which(Missing_positions_indices==1),-which(Missing_positions_indices==1)]
    
    # loading lda results so as to acces sort_normal_topic
    setwd(local_subdirname)
    load("sort_normal_topic.Rdata")
    #   filename = paste(filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best_keep.Rdata",sep="")
    #   load(filename)
    if (!occurrence_vect[barcode_index])
    {
      if (barcode_insert == "Archees_16S") 
        topic_order = rev(sort_normal_topic$ix)[c(1,3,2)]
      else if (barcode_insert == "Protistes_18S")
        topic_order = rev(sort_normal_topic$ix)[c(1,3,2)]
      else if (barcode_insert == "Plantes_GH")
        topic_order = c(3,2,1)
#       else if (barcode_insert == "Nematodes_18S")
#         topic_order = rev(sort_normal_topic$ix)[c(2,3,1)]
      else if (barcode_insert == "Annelides_18S")
        topic_order = rev(sort_normal_topic$ix)[c(1,3,2)]
      else 
        topic_order = rev(sort_normal_topic$ix)
    } else if (occurrence_vect[barcode_index])
    {
      if (barcode_insert == "Platyhelminthes_18S") 
        topic_order = rev(sort_normal_topic$ix)[c(1,3,2)]
      else if (barcode_insert == "Annelides_18S")
        topic_order = rev(sort_normal_topic$ix)[c(1,3,2)]
      else
        topic_order = rev(sort_normal_topic$ix)
    }
    
    setwd(subsubsubdirname)
    spatial_topicmix_kriged = readRDS("Spatial_topicmix_kriged.rds")
    
    if (vertical)
      tmp.plot[[barcode_index]] = list()
    
    # Plotting all topics on the same plot
    if (one_plot)
    {
      color.pal = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
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
        #k0 = rev(sort_normal_topic$ix)[k]
        k0 = topic_order[k]
        spatial_topicmix_kriged_all_topics = cbind(spatial_topicmix_kriged_all_topics,spatial_topicmix_kriged[[k0]]$z.pred)
        colnames(spatial_topicmix_kriged_all_topics)[2+k] = paste0("z.pred",k)
      }
      
      # Assigning a single assemblage (the most abudant) to each location:
      spatial_topicmix_kriged_all_topics_discrete = spatial_topicmix_kriged_all_topics
      spatial_topicmix_kriged_all_topics_discrete[,3:(2+nb_topics)] = 0
      #test = apply(spatial_topicmix_kriged_all_topics[,3:(2+nb_topics)],1,sort.int,index.return = T,decreasing=T)
      dominant_topic_proportion = apply(spatial_topicmix_kriged_all_topics[,3:(2+nb_topics)],1,max)
      dominant_topic_index = apply((spatial_topicmix_kriged_all_topics[,3:(2+nb_topics)] == dominant_topic_proportion),1,which)
      for (i in 1:nrow(spatial_topicmix_kriged_all_topics))
        spatial_topicmix_kriged_all_topics_discrete[i,2+dominant_topic_index[[i]]] = 1
      
      #       spatial_topicmix_kriged_all_topics_gradient = data.frame(spatial_topicmix_kriged_all_topics[,1:2],z.pred.grad=vector(length=length(spatial_topicmix_kriged_all_topics[,1]),mode="numeric"))
      #       for (k in 1:nb_topics)
      #         spatial_topicmix_kriged_all_topics_gradient$z.pred.grad = c(0,abs(diff(spatial_topicmix_kriged_all_topics[,2+k],lag=1))) +
      #                                                                   spatial_topicmix_kriged_all_topics_gradient$z.pred.grad
      #       spatial_topicmix_kriged_all_topics_gradient$z.pred.grad[spatial_topicmix_kriged_all_topics_gradient$x == -720] = 0
      
      # Defining the color in each location based on the assemblage composition:
      spatial_topicmix_kriged_all_topics_colors = data.frame(spatial_topicmix_kriged_all_topics[,1:2],red=0,green=0,blue=0)
      coord_matrix = t(as.matrix(spatial_topicmix_kriged_all_topics[,3:(2+nb_topics)]))
      spatial_topicmix_kriged_all_topics_discrete_colors = data.frame(spatial_topicmix_kriged_all_topics[,1:2],red=0,green=0,blue=0)
      coord_matrix_discrete = t(as.matrix(spatial_topicmix_kriged_all_topics_discrete[,3:(2+nb_topics)]))
      #       if (nb_topics == 8)
      #       {
      #         col_matrix = col2rgb(color.pal(nb_topics)[c(1,5,2,6,3,7,4,8)])
      #       } else 
      if (nb_topics == 3)
      {
        #col_matrix = col2rgb(color.pal(7)[c(2,4,6)])
        col_matrix = matrix(data=c(c(0,0,255),c(0,255,0),c(255,0,0)),ncol=3)
        rownames(col_matrix) = c("red","green","blue")
      } else
        col_matrix = col2rgb(color.pal(nb_topics))
      spatial_topicmix_kriged_all_topics_colors[,3:5] = t(col_matrix %*% coord_matrix)
      spatial_topicmix_kriged_all_topics_discrete_colors[,3:5] = t(col_matrix %*% coord_matrix_discrete)
      
      if (!discrete)
      {
      # Plotting each topic as a distinct color on a single map
      current_plot = ggplot(data = spatial_topicmix_kriged_all_topics_colors) +
        geom_raster(data = spatial_topicmix_kriged_all_topics_colors, aes(x, y), 
                    fill=rgb(spatial_topicmix_kriged_all_topics_colors[,3:5]/max(255,max(spatial_topicmix_kriged_all_topics_colors[,3:5]))), inherit.aes = F) +
        #geom_raster() +
        coord_equal() + theme_minimal() +
        labs(fill=barcode_labels_list[barcode_index]) + ggtitle(letters[barcode_index]) +
        #labs(fill=paste0("Assemblage ",k)) + ggtitle(letters[k]) +
        theme(legend.position="bottom", legend.text=element_text(size=13), 
              legend.title=element_text(size=16), axis.title=element_blank(), 
              axis.text = element_blank(),
#               plot.title=element_text(hjust=0,size=16), plot.margin=unit(c(0,1,-2,2),"mm")) +
              plot.title=element_text(hjust=0,size=10), plot.margin=unit(c(2,-10,0,-10),"mm")) +
        guides(fill = guide_colorbar(barwidth = 8, barheight = 0.5, title.position="bottom")) +
        scale_y_continuous(limits=c(5,395), expand = c(0,0)) +
        scale_x_continuous(limits=c(5,295), expand = c(0,0)) + 
        geom_point(data = data.frame(coordPP,z.pred=rep(0,nrow(coordPP))), aes(x,y), color="black", size=0.1, alpha=0.3, inherit.aes = F)
      
      } else if (discrete)
      {
        # Plotting only the dominant topic as a distinct color in each pixel
        current_plot = ggplot(data = spatial_topicmix_kriged_all_topics_discrete_colors) +
          geom_raster(data = spatial_topicmix_kriged_all_topics_discrete_colors, aes(x, y),
                      fill=rgb(spatial_topicmix_kriged_all_topics_discrete_colors[,3:5]/max(255,max(spatial_topicmix_kriged_all_topics_discrete_colors[,3:5]))), inherit.aes = F) +
          #geom_raster() +
          coord_equal() + theme_minimal() +
          labs(fill="") + ggtitle(letters[barcode_index]) +
          theme(legend.position="bottom", legend.text=element_text(size=13), 
                legend.title=element_text(size=16), axis.title=element_blank(), 
                #           theme(legend.position="bottom", legend.text=element_text(size=11), 
                #                 legend.title=element_text(size=12), axis.title=element_blank(), 
                axis.text = element_blank(),
                plot.title=element_text(hjust=0,size=16), plot.margin=unit(c(0,1,-2,2),"mm")) +
          #                 plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
          guides(fill = guide_colorbar(barwidth = 8, barheight = 0.5, title.position="bottom")) +
          scale_y_continuous(limits=c(5,395), expand = c(0,0)) +
          scale_x_continuous(limits=c(5,295), expand = c(0,0)) + 
          geom_point(data = data.frame(coordPP,z.pred=rep(0,nrow(coordPP))), aes(x,y), color="black", size=0.2, alpha=0.3, inherit.aes = F)
      }
      
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
        lidar.plot[[barcode_index]] = qplot(x, y, data = data_lidar_3topics[[barcode_index]], geom="raster", fill=layer) +
          #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", function(x) round(x,0))) + 
          #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", scientific_format())) +
          scale_fill_gradientn(colours=color.pal(7)) +
          coord_equal() +
          labs(fill=map_labels[barcode_index]) +  theme_minimal() +
          #ggtitle(letters[k+length(barcode_insert_list)*nb_topics]) +
          ggtitle(letters[length(barcode_insert_list)+barcode_index]) +
          scale_y_continuous(limits=c(5,395), expand = c(0,0)) +
          scale_x_continuous(limits=c(5,295), expand = c(0,0)) + 
          theme(legend.position="bottom", legend.text=element_text(size=8), 
                legend.title=element_text(size=10), axis.title=element_blank(), 
                #           theme(legend.position="bottom", legend.text=element_text(size=11), 
                #                 legend.title=element_text(size=12), axis.title=element_blank(), 
                axis.text = element_blank(),
#                 plot.title=element_text(hjust=0,size=16), plot.margin=unit(c(0,1,-2,2),"mm")) +
                plot.title=element_text(hjust=0,size=10), plot.margin=unit(c(2,-10,0,-10),"mm")) +
          #                 plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
          guides(fill = guide_colorbar(barwidth = 6, barheight = 0.3, title.position="bottom"))
        # size = 14 and 15 for vertical 5 topics
      }
    } else
    #####################################################
    {
      
      for (k in 1:nb_topics)
      {
        k0 = topic_order[k]
        
        if (horizontal)
        {
          if (barcode_index == 1)
            tmp.plot[[k]] = list()
        }
        
        current_plot = qplot(x, y, data=spatial_topicmix_kriged[[k0]], geom="raster", fill=z.pred) +
          #tmp.plot[[k]] = qplot(x, y, spatial_topicmix_Blaise, geom="raster", fill=z.pred) +
          #tmp.plot[[k]] = qplot(x, y, data=spatial_topicmix_kriged, geom="raster", fill=z.pred) +
          #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", function(x) round(x,0))) + 
          #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", scientific_format())) +
          scale_fill_gradientn(colours=color.pal(7)) +
          coord_equal() +
          #labs(fill=paste("Assemblage",k,barcode_labels_list[barcode_index])) +
          labs(fill=paste("Assemblage",k)) +
          #labs(fill=paste(Assemblage_labels_list[k],"-",barcode_labels_list[barcode_index])) +
          #labs(fill=Assemblage_labels_list[k]) +
          theme_minimal() + 
          scale_y_continuous(limits=c(5,395), expand = c(0,0)) +
          scale_x_continuous(limits=c(5,295), expand = c(0,0)) + 
          theme(legend.position="bottom", legend.text=element_text(size=13), 
                legend.title=element_text(size=25), axis.title=element_blank(), 
                #           theme(legend.position="bottom", legend.text=element_text(size=11), 
                #                 legend.title=element_text(size=12), axis.title=element_blank(), 
                axis.text = element_blank(),
                #               plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
                #               plot.title=element_text(hjust=0,size=16), plot.margin=unit(c(0,1,-2,2),"mm")) +
                plot.title=element_text(hjust=0,size=16), plot.margin=unit(c(0,2,-2,2),"mm")) +
          guides(fill = guide_colorbar(barwidth = 10, barheight = 0.5, title.position="bottom")) +
          geom_point(data = data.frame(coordPP,z.pred=rep(0,nrow(coordPP))), aes(x,y), color="black", size=0.2, alpha=0.3)
        
        if (horizontal)
        {
          current_plot = current_plot + ggtitle(letters[(k-1)*(length(barcode_insert_list)+1) + barcode_index])
          tmp.plot[[k]][[barcode_index]] = current_plot
        }
        else if (vertical)
        {
          #current_plot = current_plot + ggtitle(letters[(barcode_index-1)*nb_topics+k])
          tmp.plot[[barcode_index]][[k]] = current_plot
        }
        
        #sauver l'objet en pdf
        #ggsave(filename = "Topic_ordered_by_site-normalized_abundance_composition_maps_kriged.pdf", tmp.plot)
        
        #Tu peux mettre plusieurs cartes (plusieurs tmp.plot) dans une liste et arranger la sortie sur plusieurs panels de la manière suivante (ici mettons que tu as 14 cartes et que tu veux les montrer sur 7 colonnes (et donc 2 lignes):
        #         ggsave(filename = "Figchim.pdf", do.call("arrangeGrob", c(list.tmp.plot, ncol=7)), height = 8, width = 10)
        
        if ((barcode_index == length(barcode_insert_list)) && horizontal)
        {
          tmp.plot[[k]][[barcode_index+1]] = qplot(x, y, data = data_lidar_3topics[[k]], geom="raster", fill=layer) +
            #lidar.plot[[k]] = qplot(x, y, data = data_lidar_3topics[[k]], geom="raster", fill=layer) +
            #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", function(x) round(x,0))) + 
            #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", scientific_format())) +
            scale_fill_gradientn(colours=color.pal(7)) +
            coord_equal() +
            labs(fill=map_labels[k]) +  theme_minimal() + ggtitle(letters[(k-1)*(length(barcode_insert_list)+1) + barcode_index + 1]) +
            scale_y_continuous(limits=c(5,395), expand = c(0,0)) +
            scale_x_continuous(limits=c(5,295), expand = c(0,0)) + 
            theme(legend.position="bottom", legend.text=element_text(size=13), 
                  legend.title=element_text(size=16), axis.title=element_blank(), 
                  #           theme(legend.position="bottom", legend.text=element_text(size=11), 
                  #                 legend.title=element_text(size=12), axis.title=element_blank(), 
                  axis.text = element_blank(),
                  plot.title=element_text(hjust=0,size=16), plot.margin=unit(c(0,2,-2,2),"mm")) +
            #                 plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
            guides(fill = guide_colorbar(barwidth = 10, barheight = 0.5, title.position="bottom"))
          # size = 14 and 15 for vertical 5 topics
        } else if ((barcode_index == 1) && vertical)
        {          
          lidar.plot[[k]] = qplot(x, y, data = data_lidar_3topics[[k]], geom="raster", fill=layer) +
            #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", function(x) round(x,0))) + 
            #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", scientific_format())) +
            scale_fill_gradientn(colours=color.pal(7)) +
            labs(fill=map_labels[k]) +  theme_minimal() +
            #ggtitle(letters[k+length(barcode_insert_list)*nb_topics]) +
            scale_y_continuous(limits=c(5,395), expand = c(0,0)) +
            scale_x_continuous(limits=c(5,295), expand = c(0,0)) + 
            theme(legend.position="bottom", legend.text=element_text(size=16), 
                  legend.title=element_text(size=20), axis.title=element_blank(), 
                  #           theme(legend.position="bottom", legend.text=element_text(size=11), 
                  #                 legend.title=element_text(size=12), axis.title=element_blank(), 
                  axis.text = element_blank(),
                  plot.title=element_text(hjust=0,size=16), plot.margin=unit(c(0,1,-2,2),"mm")) +
            #                 plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
            guides(fill = guide_colorbar(barwidth = 8, barheight = 0.5, title.position="bottom"))
          # size = 14 and 15 for vertical 5 topics
        }
      }
    }
    # End of the barcode_insert loop:
  }
  
  setwd("/Users/guilhemsommeria-klein/Desktop/These/Manuscrits/Paper2/Figures/")
  if (vertical)
  {
    if (one_plot)
    {
#       ggsave(filename = paste0("Kriged_topics_with_Lidar_bact-prot-fungi_ab-",occurrence_insert,remove_single_sites_insert,no_rare_insert,horizVert_insert,"_oneplot",discrete_insert,".pdf"),
#              do.call("arrangeGrob", 
#                                  c(tmp.plot, lidar.plot, ncol=3)),
# #                                  c(tmp.plot, ncol=3)),
# #                                    c(tmp.plot, ncol=1)),
# #              height = 10/3*4/3*(2+1), width = 10)
#               heights = c(rep(10/3*4/3*0.8,6),rep(10/3*4/3,6)))
# #              height = 10/3*4/3, width = 10/3)
      
      #pdf(paste0("Kriged_topics_with_Lidar_bact-prot-fungi_ab-",occurrence_insert,remove_single_sites_insert,no_rare_insert,horizVert_insert,"_oneplot",discrete_insert,".pdf"))
      #pdf(paste0("Kriged_topics_bact-prot-fungi_occ_noSingleSite_noRare",horizVert_insert,"_oneplot",discrete_insert,".pdf"))
      pdf(paste0("Kriged_topics_arch-arth-nema-plath-anne_occ-ab",horizVert_insert,"_oneplot",discrete_insert,".pdf"))
      #grid.arrange(grobs = c(tmp.plot, lidar.plot),ncol=3,heights = unit(c(10*4/3*0.7,10*4/3*0.7,10*4/3),"cm"), widths=unit(rep(10,3),"cm"))
      # Figure 4 :
      #grid.arrange(grobs = c(tmp.plot, lidar.plot),ncol=3,heights = c(10*4/3*0.72,10*4/3*0.72,10*4/3), width = 10)
      grid.arrange(grobs = tmp.plot, ncol=3, width = 10, height = 10/3*4*4/3, layout_matrix = rbind(c(1,2,3),c(4,5,6),c(7,8,NA),c(9,10,NA)))
      dev.off()
    } else
    {
    #     ggsave(filename = "Kriged_topics_with_Lidar_5barcodes_vertical.pdf", do.call("arrangeGrob", 
    #                 c(tmp.plot[[1]], tmp.plot[[2]], tmp.plot[[3]], tmp.plot[[4]], tmp.plot[[5]], lidar.plot, ncol=nb_topics)),
    #                 height = 10/3*4/3*(length(barcode_insert_list)+1), width = 10)
    #     ggsave(filename = "Kriged_topics_with_Lidar.pdf", do.call("arrangeGrob", 
    #                 c(tmp.plot[[1]], tmp.plot[[2]], tmp.plot[[3]], lidar.plot, ncol=nb_topics)),
    #                 height = 10/3*4/3*(length(barcode_insert_list)+1), width = 10)
    #     ggsave(filename = "Kriged_topics_archaea-protists.pdf", do.call("arrangeGrob", 
    #             c(tmp.plot[[1]], tmp.plot[[2]], ncol=nb_topics)),
    #            height = 10/3*4/3*length(barcode_insert_list), width = 10)
#       ggsave(filename = "Kriged_topics_with_Lidar_bacteria.pdf", do.call("arrangeGrob", 
#                                                                        c(tmp.plot[[1]], lidar.plot, ncol=nb_topics)),
#            height = 10/3*4/3*(length(barcode_insert_list)+1), width = 10)
      ggsave(filename = "Kriged_topics_bacteria.pdf", do.call("arrangeGrob", 
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
    ggsave(filename = paste0("Kriged_topics_with_Lidar_",length(barcode_insert_list),"barcodes",occurrence_insert,remove_single_sites_insert,no_rare_insert,horizVert_insert,".pdf"), do.call("arrangeGrob", 
                                                                                                                                                    c(tmp.plot[[1]], tmp.plot[[2]], tmp.plot[[3]], nrow=nb_topics)),
           height = 10*4/3, width = 10*length(barcode_insert_list)/3)
  
}
####################################################################################################
####################################################################################################
if (stability)
{
  
#   barcode_insert_list = c("Bacteries_16S","Archees_16S","Protistes_18S","Champignons_ITS","Plantes_GH","Arthropodes_18S")
#   barcode_labels_list = c("Bacteria 16S","Archaea 16S","Protists 18S","Fungi ITS","Plants trnL","Arthropods 18S")
  
  # All barcodes :
#   barcode_insert_list = c("Bacteries_16S","Archees_16S","Protistes_18S","Champignons_ITS","Glomeromycetes_ITS","Basidiomycetes_ITS","Ascomycetes_ITS","Plantes_GH","Arthropodes_18S","Nematodes_18S","PLatyhelminthes_18S","Annelides_18S")
#   barcode_labels_list = c("Bacteria 16S","Archaea 16S","Protists 18S","Fungi ITS","Glomeromycota ITS","Basidiomycota ITS","Ascomycota ITS","Plants trnL","Arthropods 18S","Nematodes 18S","Platyhelminthes 18S","Annelids 18S")
  
#     barcode_insert_list = c("Bacteries_16S","Protistes_18S","Champignons_ITS")
#     barcode_labels_list = c("Bacteria","Protists","Fungi")
 
barcode_insert_list = c("Champignons_ITS","Bacteries_16S","Arthropodes_18S")
barcode_labels_list = c("Fungi","Bacteria","Arthropods") 

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
  
    if (barcode_index == 1)
    {
      plot.skl[[barcode_index]] = plot.skl[[barcode_index]] + labs(x="Realizations\n sorted by decreasing likelihood", y="SKL to best realization")
      plot.skl.llh[[barcode_index]] = plot.skl.llh[[barcode_index]] + labs(x="Llh difference\n with best realization", y="SKL to best realization")
      #dskl.plot.llh[[barcode_index]] = dskl.plot.llh[[barcode_index]] + labs(x="Loglikelihood difference\n with best realization", y=expression(paste("SKL-<SKL>",[random.]*" to best realization"))) bquote(.(labNames[1]) ~ x^2)
      plot.dskl.llh[[barcode_index]] = plot.dskl.llh[[barcode_index]] + labs(x="Llh difference\n with best realization", y=bquote(.("SKL-<SKL>")[random.]*" to best realization"))
      #plot.dskln.llh[[barcode_index]] = plot.dskln.llh[[barcode_index]] + labs(x="Llh difference\n with best realization", y=bquote(.("(SKL-<SKL>")[rnd]*")/<SKL>"[rnd]))
      plot.dskln.llh[[barcode_index]] = plot.dskln.llh[[barcode_index]] + labs(x="\nLlh difference\n with best realization", y="Spatial similarity\n to best realization\n")
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
      plot.dskln.llh[[barcode_index]] = plot.dskln.llh[[barcode_index]] + labs(x="\nLlh difference\n with best realization", y="")
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
  
  setwd("/Users/guilhemsommeria-klein/Desktop/These/Manuscrits/Paper2/Figures/")
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
ggsave(filename = paste0("NormES_vs_llh-diff_fungi-restrictedBact-arth_withLM",occurrence_insert,"_",MOTU_sample_insert,"_",bij_insert,".pdf"), do.call("arrangeGrob", 
                                                                                                                                c(plot.dskln.llh, nrow=1)),
       height = 15/3, width = 15)

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
if (size)
{
  library(plotrix)
  
  setwd("/Users/guilhemsommeria-klein/Desktop/These/Manuscrits/Paper2/Figures/")
  SKL_bij_samplewise = c(0.27,0.58,1.9,1.5,1.1,2.0,3.0,3.1,3.0)
  ES_bij_samplewise = c(2.0,1.5,3.1,1.6,0.86,1.6,2.5,1.4,0.95)
  SES_bij_samplewise = c(1.3,1.3,0.81,0.47,0.93,0.95,0.75,0.86,1.5)
  pval_bij_samplewise = c(0,0.029,0.069,0.17,0.082,0.13,0.042,0.091,0.092)
  
  SKL_bij_samplewise1 = c(0.27,0.58,1.9,0.025,1.1,2.0,3.0,3.1,3.0)
  ES_bij_samplewise1 = c(2.0,1.5,3.1,2.7,0.86,1.6,2.5,1.4,0.95)
  SES_bij_samplewise1 = c(1.3,1.3,0.81,3.5,0.93,0.95,0.75,0.86,1.5)
  pval_bij_samplewise1 = c(0,0.029,0.069,0,0.082,0.13,0.042,0.091,0.092)
  
  #Occ.
  nES_bij_samplewise = c(0.32,0.68,0.41,0.88,0.33,0.52,0.62,0.70,0.10)
  nES_bij_samplewise1 = c(0.85,0.68,0.41,0.88,0.33,0.52,0.62,0.70,0.10)
  
  SKL_nobij_samplewise = c(0.18,0.48,1.2,0.44,0.80,1.3,2.5,2.4,2.8)
  ES_nobij_samplewise = c(1.6,1.5,2.5,1.8,0.75,1.7,2.3,1.1,0.97)
  SES_nobij_samplewise = c(1.2,1.2,0.68,1.4,1.4,1.2,0.71,1.0,1.7)
  pval_nobij_samplewise = c(0,0.032,0.19,0.041,0.072,0.074,0.15,0.12,0.065)
  
  SKL_nobij_samplewise1 = c(0.18,0.48,1.2,0.025,0.80,1.3,2.5,2.4,2.8)
  ES_nobij_samplewise1 = c(1.6,1.5,2.5,2.7,0.75,1.7,2.3,1.1,0.97)
  SES_nobij_samplewise1 = c(1.2,1.2,0.68,3.5,1.4,1.2,0.71,1.0,1.7)
  pval_nobij_samplewise1 = c(0,0.032,0.19,0,0.072,0.074,0.15,0.12,0.065)
  
  SKL_bij_MOTUwise = c(0.11,0.56,2.6,2.2,2.4,2.8,3.6,2.2,5.4)
  ES_bij_MOTUwise = c(7.0,6.3,7.8,9.8,5.2,5.6,6.7,4.3,5.9)
  SES_bij_MOTUwise = c(9.5,4.4,2.7,5.1,4.7,3.2,3.4,9.6,7.3)
  pval_bij_MOTUwise = c(0,0,0.045,0.032,0.0091,0.028,0.0071,0,0.00056)
  
  SKL_bij_MOTUwise1 = c(0.11,0.56,2.6,0.031,2.4,2.8,3.6,2.2,5.4)
  ES_bij_MOTUwise1 = c(7.0,6.3,7.8,12,5.2,5.6,6.7,4.3,5.9)
  SES_bij_MOTUwise1 = c(9.5,4.4,2.7,7,4.7,3.2,3.4,9.6,7.3)
  pval_bij_MOTUwise1 = c(0,0,0.045,0,0.0091,0.028,0.0071,0,0.00056)
  
  #Occ.
  nES_bij_MOTUwise = c(0.86,0.95,0.83,0.98,0.88,0.86,0.91,0.90,0.85)
  nES_bij_MOTUwise1 = c(0.95,0.95,0.83,0.98,0.88,0.86,0.91,0.90,0.85)
  
  #taxo_names = c("Bacteria","Protists","Annelids","Archaea","Nematodes","Platyhelminthes","Arthropods","Fungi","Plants")
  taxo_names = c("Bacteria","Protists","Annelids","Archaea","Nematodes","Platyhelminthes","Arthropods","Fungi")
  size = log10(c(mean(c(0.5*10^-6,5*10^-6,0.5*10^-6,10^-6,10^-6,0.5*10^-6,5*10^-6,0.5*10^-6,5*10^-6,0.5*10^-6)),
                 100*10^-6,
                 20e-3,
                 0.5*10^-6,
                 100*10^-6,
                 20*10^-3,
                 10^-3,
                 mean(c(100*10^-6,10*10^-6,200*10^-6)),
                 0.5*10^2))
  
  diversity0 = log10(c(20162,1648,51,4101,378,126,1881,9855,1360))
  
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
  data.tmp1[[1]] = data.frame(x=size[-length(size)],y=nES_bij_samplewise1[-length(nES_bij_samplewise1)])
  data.tmp1[[2]] = data.frame(x=size[-length(size)],y=nES_bij_MOTUwise1[-length(nES_bij_MOTUwise1)])
  
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
  for (i in 1:2)
  {
    diversity = diversity0[-length(diversity0)]
    
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
                nudge_y = c(0,-0.01,0,0.02,0,0,0,0.02), hjust=0, nudge_x = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2),
                parse = FALSE, check_overlap = F, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) 
    } else if (i == 2)
    {
      plot.size1[[i]] = plot.size1[[i]] + geom_text(mapping = NULL, stat = "identity",size=2, 
                #hjust=0, nudge_x = 0.2,
                nudge_y = c(0.01,0,0,0.01,0,0,0,0), hjust=0, nudge_x = c(0,0.2,0.2,0,0.2,0.2,0.2,0.2),
                parse = FALSE, check_overlap = F, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) 
    }
    #c("Bacteria","Protists","Annelids","Archaea","Nematodes","Platyhelminthes","Arthropods","Fungi","Plants")
    
    rsquared[[i]][2] = cor(data.tmp1[[i]][,2],data.tmp1[[i]][,1])^2
    
    plot.diversity[[i]] = qplot(x = diversity0, y, data=data.tmp[[i]], label=taxo_names, geom="point") +
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
      labs(x="Number of MOTUs (log10)", y=y_labels[i]) +
      xlim(min(diversity0),5) +
      geom_text(mapping = NULL, stat = "identity",size=2, hjust = 0, nudge_x = 0.1,
                parse = FALSE, check_overlap = F, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) +
      geom_smooth(method='lm')
    
    rsquared[[i]][3] = cor(data.tmp[[i]][,2],diversity0)^2
    
    plot.diversity1[[i]] = qplot(x = diversity, y, data=data.tmp1[[i]], label=taxo_names, geom="point") +
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
      labs(x="Number of OTUs (log10)", y=y_labels[i]) +
      xlim(min(diversity),5) +
      geom_smooth(method='lm') 
      
    
    if (i == 1)
    {
      plot.diversity1[[i]] = plot.diversity1[[i]] + geom_text(mapping = NULL, stat = "identity",size=2, 
                                                    #hjust=0, nudge_x = 0.2,
                                                    nudge_y = c(0,0,0,0,0,0,0,0), hjust=0, nudge_x = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1),
                                                    parse = FALSE, check_overlap = F, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) 
    } else if (i == 2)
    {
      plot.diversity1[[i]] = plot.diversity1[[i]] + geom_text(mapping = NULL, stat = "identity",size=2, 
                                                    #hjust=0, nudge_x = 0.2,
                                                    nudge_y = c(0,0,0,0,0,0,0,0), hjust=0, nudge_x = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1),
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
  
  ggsave(filename = "Stability_vs_diversity_nES_restricted_bacteria_spatial-taxo_noplants.pdf", do.call("arrangeGrob", c(plot.diversity1[1], plot.diversity1[2], nrow=1)),
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
#   barcode_insert_list = c("Champignons_ITS","Bacteries_16S","Arthropodes_18S")
#   barcode_labels_list = c("Fungi","Bacteria","Arthropods") 
  
  barcode_insert_list = c("Bacteries_16S","Protistes_18S","Champignons_ITS")
  barcode_labels_list = c("Bacteria","Protists","Fungi")  
  
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
    setwd(local_subdirname)
    
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
    
#     upgma_data.frame = dendro_data(upgma_Hellinger)
# #     upgma_data.frame = dendro_data(upgma_Jaccard)
# #     upgma_data.frame = dendro_data(upgma_Beta.sim)
# 
#     taxotree.plot[[barcode_index]] = ggplot() +
#     geom_segment(data = upgma_data.frame$segments, 
#                aes(x = x, y = y, xend = xend, yend = yend)) + 
#     geom_text(data = upgma_data.frame$labels, 
#             aes(x = x, y = y, label = label), size = 3, vjust = 1) +
# #     geom_text(data = upgma_data.frame$leaf_labels, 
# #             aes(x = x, y = y, label = label), size = 3, vjust = 1) +
#     theme_dendro() 
    #theme(plot.margin=unit(c(0,2,0,2),"mm"))

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

pdf("UPGMA_Helinger_bact-prot-fungi.pdf")
#grid.arrange(grobs = taxotree.plot, ncol=3, height = 1, width = 3)
grid.arrange(grobs = taxotree.plot, ncol=3, height = 1, width = 3)
# par(mfrow = c(1,3))
# for (barcode_index in 1:3)
#   plot(upgma_Hellinger[[barcode_index]], which.plots=2, ann=F)
dev.off()

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




