# Variables values for testing :
# locally_based = 1
# genotoul_cluster_based = 0
# EDB_cluster_based = 0
# Selected_real = 1
# j_select = 1
# kriged_real = 1
# kriging = 1
# nb_topics = 3
# data_pp = 1
# data_h20 = 0
# data_gs = 0
# filled = 1
# filled_with_gaps = 0
# documents = replicate(nb_topics,abs(rnorm(1131)))
# sort_normal_topic = sort.int(c(1,2,3),index.return=T)
# nb_real = 100
# em_tol = 10^-7
# var_tol = 10^-8
# #
# local_prefix = "/Users/guilhemsommeria-klein/Desktop/These/"
# data_insert = "Donnees_PetitPlateau"
# barcode_insert = "Protistes_18S"
# filename_insert = "Rtopicmodels_LDA_VEM"
# occurrence_insert = ""
# remove_single_sites_insert = ""
# no_rare_insert = ""
# occurrence_insert = "_occurrence"
# local_dirname = paste(local_prefix,data_insert,"/",barcode_insert,"/",filename_insert,"/",sep="")
# local_subdirname = paste0(local_dirname,filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best_keep",occurrence_insert,remove_single_sites_insert,no_rare_insert,"/")
# subsubdirname = paste0(local_subdirname,"ordered_realization_number_",Selected_real[j_select],"/")
# subsubsubdirname = paste0(subsubdirname,"/topics_site_repartition_info/")
# #

LDA_spatial_maps_fun <- function(j_select,kriged_real,kriging,nb_topics,data_pp,data_h20,
                                                data_gs,filled,filled_with_gaps,Missing_positions_indices,
                                                documents,coordGS,local_prefix,data_insert,
                                                sort_normal_topic,subsubsubdirname){
  
  spatial_topicmix = vector(length=nb_topics,mode="list")
  if (kriging && any(kriged_real == j_select))
  {
    if (!file.exists(paste0(subsubsubdirname,"Spatial_topicmix_kriged.rds")))
    {
      if (data_pp)
      {
        newdata = expand.grid(x=seq(0,300,2), y=seq(0,400,2))
        #newdata_Blaise = expand.grid(x=seq(0,400,2), y=rev(seq(0,300,2)))
      } else if (data_h20)
      {
        newdata = expand.grid(x=seq(0,100,1), y=seq(0,100,1))
      } else if (data_gs)
      {
        newdata = expand.grid(x=seq(-720,-470,1), y=seq(-110,90,1))
      }
      spatial_topicmix_kriged = vector(length=nb_topics,mode="list")
    } else if (file.exists(paste0(subsubsubdirname,"Spatial_topicmix_kriged.rds")))
    {
      spatial_topicmix_kriged = readRDS(paste0(subsubsubdirname,"Spatial_topicmix_kriged.rds"))
      #spatial_topicmix_kriged_Blaise = readRDS("Spatial_topicmix_kriged_Blaise.rds")
    }
    if (data_pp)
    {
      coordPP = expand.grid(x=seq(10,290,10), y=seq(10,390,10))
      #coordPP_Blaise = expand.grid(x=seq(10,390,10), y=rev(seq(10,290,10)))
      if (!filled || filled_with_gaps)
        coordPP = coordPP[-which(Missing_positions_indices==1),-which(Missing_positions_indices==1)]
    } else if (data_h20)
    {
      coordH20 = expand.grid(x=seq(1,19,1)*5, y=seq(1,19,1)*5)
    } 
  }
  
  # Loop over topics (one map per topic, the color stands for the proportion of the topic)
  for (k in 1:nb_topics)
  {  
    if (data_h20)
    {
      #             spatial_topicmix = matrix(nrow=19,ncol=19) 
      spatial_topicmix[[k]] = matrix(documents[,k],ncol=19)
      #             for (j in 1:19)
      #             {
      #               for (i in 1:19)
      #               {
      #                 spatial_topicmix[i,j] = documents[(j-1)*19+i,rev(sort_normal_topic$ix)[k]]
      #               }
      #             }
    } else if (data_pp)
    {             
      #             spatial_topicmix = matrix(nrow=29,ncol=39)
      if (filled && !filled_with_gaps)
      {
        spatial_topicmix[[k]] = matrix(documents[,k],ncol=39)
        #               for (j in 1:39)
        #               {
        #                 for (i in 1:29)
        #                 {
        #                   spatial_topicmix[i,j] = documents[(j-1)*29+i,rev(sort_normal_topic$ix)[k]]
        #                 }
        #               }
      } else
      {
        spatial_topicmix[[k]] = matrix(nrow=29,ncol=39,data=0)
        position_shift = 0
        for (j in 1:39)
        {
          for (i in 1:29)
          {
            if (Missing_positions_indices[(j-1)*29+i]==0)
              spatial_topicmix[[k]][i,j] = documents[(j-1)*29+i-position_shift,k]    
            else if (Missing_positions_indices[(j-1)*29+i]==1)
            {
              spatial_topicmix[[k]][i,j] = NA
              position_shift = position_shift+1
            }
          }
        }
      }
    }
    if (kriging && any(kriged_real == j_select))
    {
      #         ycoor = rep(1,nrow(spatial_topicmix))
      #         for (i in 2:ncol(spatial_topicmix))
      #           ycoor = c(ycoor,rep(i,nrow(spatial_topicmix)))
      #         coordPP = data.frame(x=rep(seq(1,nrow(spatial_topicmix),1),ncol(spatial_topicmix))*10,y=ycoor*10)
      if (!file.exists(paste0(subsubsubdirname,"Spatial_topicmix_kriged.rds")))
      {
        if (data_pp)
          tmp = data.frame(x=coordPP$x,y=coordPP$y,z=documents[,k])
        else if (data_h20)
          tmp = data.frame(x=coordH20$x,y=coordH20$y,z=documents[,k])
        else if (data_gs)
          tmp = data.frame(x=coordGS$x,y=coordGS$y,z=documents[,k])
        mod = gstat(id = "z", formula = z~1, locations = ~x+y, data=tmp, model=vgm(10, "Exp", 20))
        spatial_topicmix_kriged[[k]] = predict(mod, newdata, na.action=na.omit)
      }
      #spatial_topicmix_kriged[[k]] = kriging(rep(seq(1,nrow(spatial_topicmix),1),ncol(spatial_topicmix)),ycoor,
      #as.vector(spatial_topicmix), model = "spherical", lags = 10, pixels = 500, polygons = NULL)
    }
  } 
  
  
  if (kriging && any(kriged_real == j_select))
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
      k0 = rev(sort_normal_topic$ix)[k]
      spatial_topicmix_kriged_all_topics = cbind(spatial_topicmix_kriged_all_topics,spatial_topicmix_kriged[[k0]]$z.pred)
      colnames(spatial_topicmix_kriged_all_topics)[2+k] = paste0("z.pred",k)
    }
    
    # Assigning a single assemblage (the most abudant) to each location:
    spatial_topicmix_kriged_all_topics_discrete = spatial_topicmix_kriged_all_topics
    spatial_topicmix_kriged_all_topics_discrete[,3:(2+nb_topics)] = 0
    #test = apply(spatial_topicmix_kriged_all_topics[,3:(2+nb_topics)],1,sort.int,index.return = T,decreasing=T)
    dominant_topic_proportion = apply(spatial_topicmix_kriged_all_topics[,3:(2+nb_topics)],1,max)
    dominant_topic_index = apply((spatial_topicmix_kriged_all_topics[,3:(2+nb_topics)] == dominant_topic_proportion),1,which)
    # dominant_topic_index is a list if two topics have equal proportions in one site, otherwise it is a vector  
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
    
    tmp.plot = list()  
    for (k in 1:nb_topics)
    {
      
      k0 = rev(sort_normal_topic$ix)[k]
      # image(spatial_topicmix_kriged[[k]], col = terrain.colors(100))
      #         axis(2, at=c(0,5,10,15,20,25,30), labels=c("","50","100","150","200","250","300"))
      #         axis(1, at=c(0,5,10,15,20,25,30,35,40), labels=c("0","","100","","200","","300","","400"))
      
      # Shifting from Lucie's indexing to Blaise's indexing
      #         spatial_topicmix_Blaise = spatial_topicmix_kriged[[k0]]
      #         spatial_topicmix_Blaise$y = spatial_topicmix_kriged[[k0]]$x
      #         spatial_topicmix_Blaise$x = spatial_topicmix_kriged[[k0]]$y
      #         spatial_topicmix_Blaise_wrongindexing = spatial_topicmix_Blaise
      #         coordPP_Blaise = coordPP
      #         coordPP_Blaise$x = coordPP$y
      #         coordPP_Blaise$y = coordPP$x
      #         coordPP_Blaise_wrongindexing = coordPP_Blaise
      #         for (i in 1:201)
      #         {
      #           for (j in 1:151)
      #           {
      #             spatial_topicmix_Blaise[(151-j)*201+i,] = spatial_topicmix_Blaise_wrongindexing[(i-1)*151+j,]
      #           }
      #         }
      #         for (i in 1:39)
      #         {
      #           for (j in 1:29)
      #             coordPP_Blaise[(29-j)*39+i,] = coordPP_Blaise_wrongindexing[(i-1)*29+j,]
      #         }
      
      if (data_gs)
      {
        #plotting each topic separately as a color gradient 
        df.plot = data.frame(x=spatial_topicmix_kriged[[k0]]$x/10,y=spatial_topicmix_kriged[[k0]]$y/10,z.pred=spatial_topicmix_kriged[[k0]]$z.pred)
        df.plot.raster = rasterFromXYZ(df.plot)
        df.plot.extruded = extract(df.plot.raster,land,df=T,cellnumbers=T)
        df.plot.extruded.data.frame = data.frame(x=spatial_topicmix_kriged_all_topics_colors$x[df.plot.extruded$cell],
                                                 y=spatial_topicmix_kriged_all_topics_colors$y[df.plot.extruded$cell],df.plot.extruded[,3])
        df.plot.extruded.raster = rasterFromXYZ(df.plot.extruded.data.frame)
        df.plot.extruded.raster = flip(df.plot.extruded.raster,'y')
        df.plot.extruded.data.frame = as.data.frame(df.plot.extruded.raster, row.names=NULL, optional=T, xy=TRUE, na.rm=T, long=FALSE)
        df.plot = data.frame(x=df.plot.extruded.data.frame$x,y=df.plot.extruded.data.frame$y,z.pred=df.plot.extruded.data.frame[,3])
        
        tmp.plot[[k]] = ggplot(data = bordersGgplot) +
          #theme(legend.position = "none") +
          geom_raster(data = df.plot, aes(x, y, fill=z.pred), inherit.aes = F) +
          scale_fill_gradientn(colours=color.pal(7)) +
          #geom_polygon(data = oceanGgplotNew, aes(x=long, y=lat, group=group, fill=piece)) +
          geom_path(data = bordersGgplot, aes(x=long*10, y=lat*10, group=group), 
                    color = "black", size=0.1, inherit.aes = F) +
          geom_path(data = riversGgplot, aes(x=long*10, y=lat*10, group=group), 
                    color = "blue", size=0.4, inherit.aes = F, alpha=1) +
          scale_y_continuous(limits=c(-110,90), expand = c(0,0)) +
          scale_x_continuous(limits=c(-720,-470), expand = c(0,0)) +
          geom_point(data = data.frame(coordGS,z.pred=rep(0,nrow(coordGS))), 
                     aes(x,y), color="black", size=0.2, alpha=0.5, inherit.aes = F) +   
          coord_equal() +
          labs(fill=paste0("Assemblage ",k)) +  theme_minimal() + ggtitle(letters[k]) +
          #           scale_x_continuous(limits=c(5,395), expand = c(0,0)) +
          #           scale_y_continuous(limits=c(5,295), expand = c(0,0)) + 
          theme(legend.position="bottom", legend.text=element_text(size=7), 
                legend.title=element_text(size=8), axis.title=element_blank(), 
                axis.text = element_blank(), panel.background = element_blank(),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
          guides(fill = guide_colorbar(barwidth = 8, barheight = 0.4, title.position="bottom"))
        
        if (k==1)
        {
          # Plotting each topic as a distinct color on a single map
          df.plot = data.frame(x=spatial_topicmix_kriged_all_topics_colors$x/10,y=spatial_topicmix_kriged_all_topics_colors$y/10,spatial_topicmix_kriged_all_topics_colors[,3:5])
          df.plot.raster = rasterFromXYZ(df.plot)
          df.plot.extruded = extract(df.plot.raster,land,df=T,cellnumbers=T)
          df.plot.extruded.data.frame = data.frame(x=spatial_topicmix_kriged_all_topics_colors$x[df.plot.extruded$cell],
                                                   y=spatial_topicmix_kriged_all_topics_colors$y[df.plot.extruded$cell],df.plot.extruded[,3:5])
          df.plot.extruded.raster = rasterFromXYZ(df.plot.extruded.data.frame)
          df.plot.extruded.raster = flip(df.plot.extruded.raster,'y')
          df.plot.extruded.data.frame = as.data.frame(df.plot.extruded.raster, row.names=NULL, optional=T, xy=TRUE, na.rm=T, long=FALSE)
          df.plot = data.frame(x=df.plot.extruded.data.frame$x,y=df.plot.extruded.data.frame$y,df.plot.extruded.data.frame[,3:5])
          
          #               land.plot = ggplot(data = landGgplot) +
          #                 geom_path(data = landGgplot, aes(x=long*10, y=lat*10, group=group), 
          #                          color = "black", size=0.1) 
          #               ggsave(filename = "Test2.pdf", land.plot, width = 10)
          
          tmp.one.plot = ggplot(data = bordersGgplot) +
            #geom_raster(data = spatial_topicmix_kriged_all_topics_colors, aes(x, y), fill=rgb(spatial_topicmix_kriged_all_topics_colors[,3:5]/max(255,max(spatial_topicmix_kriged_all_topics_colors[,3:5]))), inherit.aes = F) +
            geom_raster(data = df.plot, aes(x, y), fill=rgb(df.plot[,3:5]/max(255,max(df.plot[,3:5]))), inherit.aes = F) +
            #geom_raster(data = test_crop_data.frame, aes(x, y), fill=rgb(test_crop_data.frame[,3:5]/max(255,max(test_crop_data.frame[,3:5]))), inherit.aes = F) +
            #geom_raster(data = spatial_topicmix_kriged_all_topics_gradient, aes(x, y, fill=spatial_topicmix_kriged_all_topics[,3:(2+nb_topics)]), inherit.aes = F) +
            #geom_raster(data = spatial_topicmix_kriged_all_topics_gradient, aes(x, y, fill=z.pred.grad), inherit.aes = F) +
            #scale_fill_gradientn(colours=color.pal(7)) +
            #                 geom_path(data = landGgplot, aes(x=long*10, y=lat*10, group=group), 
            #                           color = "black", size=0.1, inherit.aes = F) +
            geom_path(data = bordersGgplot, aes(x=long*10, y=lat*10, group=group), 
                      color = "black", size=0.1, inherit.aes = F) +
            geom_path(data = riversGgplot, aes(x=long*10, y=lat*10, group=group), 
                      color = "blue", size=0.4, inherit.aes = F, alpha=1) +
            scale_y_continuous(limits=c(-110,90), expand = c(0,0)) +
            scale_x_continuous(limits=c(-720,-470), expand = c(0,0)) +
            geom_point(data = data.frame(coordGS,z.pred=rep(0,nrow(coordGS))), 
                       aes(x,y), color="black", size=0.2, alpha=0.5, inherit.aes = F) +   
            coord_equal() +
            theme_minimal() +
            # labs(fill=paste0("Assemblage ",k)) + ggtitle(letters[k]) +
            #           scale_x_continuous(limits=c(5,395), expand = c(0,0)) +
            #           scale_y_continuous(limits=c(5,295), expand = c(0,0)) + 
            theme(legend.position="bottom", legend.text=element_text(size=7), 
                  legend.title=element_text(size=8), axis.title=element_blank(), 
                  axis.text = element_blank(), panel.background = element_blank(),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
            guides(fill = guide_colorbar(barwidth = 8, barheight = 0.4, title.position="bottom"))
          
          #ggsave(filename = "Test1.pdf", tmp.one.plot, width = 10) 
          
          ##########################################
          # Plotting only the dominant topic as a distinct color in each pixel
          
          df.plot = data.frame(x=spatial_topicmix_kriged_all_topics_colors$x/10,y=spatial_topicmix_kriged_all_topics_colors$y/10,spatial_topicmix_kriged_all_topics_discrete_colors[,3:5])
          df.plot.raster = rasterFromXYZ(df.plot)
          df.plot.extruded = extract(df.plot.raster,land,df=T,cellnumbers=T)
          df.plot.extruded.data.frame = data.frame(x=spatial_topicmix_kriged_all_topics_colors$x[df.plot.extruded$cell],
                                                   y=spatial_topicmix_kriged_all_topics_colors$y[df.plot.extruded$cell],df.plot.extruded[,3:5])
          df.plot.extruded.raster = rasterFromXYZ(df.plot.extruded.data.frame)
          df.plot.extruded.raster = flip(df.plot.extruded.raster,'y')
          df.plot.extruded.data.frame = as.data.frame(df.plot.extruded.raster, row.names=NULL, optional=T, xy=TRUE, na.rm=T, long=FALSE)
          df.plot = data.frame(x=df.plot.extruded.data.frame$x,y=df.plot.extruded.data.frame$y,df.plot.extruded.data.frame[,3:5])
          
          tmp.one.plot.discrete = ggplot(data = bordersGgplot) +
            geom_raster(data = df.plot, aes(x, y), 
                        fill=rgb(df.plot[,3:5]/max(255,max(df.plot[,3:5]))), inherit.aes = F) +
            #geom_raster(data = spatial_topicmix_kriged_all_topics_gradient, aes(x, y, fill=spatial_topicmix_kriged_all_topics[,3:(2+nb_topics)]), inherit.aes = F) +
            #geom_raster(data = spatial_topicmix_kriged_all_topics_gradient, aes(x, y, fill=z.pred.grad), inherit.aes = F) +
            #scale_fill_gradientn(colours=color.pal(7)) +
            geom_path(data = bordersGgplot, aes(x=long*10, y=lat*10, group=group), 
                      color = "black", size=0.1, inherit.aes = F) +
            geom_path(data = riversGgplot, aes(x=long*10, y=lat*10, group=group), 
                      color = "blue", size=0.4, inherit.aes = F, alpha=1) +
            scale_y_continuous(limits=c(-110,90), expand = c(0,0)) +
            scale_x_continuous(limits=c(-720,-470), expand = c(0,0)) +
            geom_point(data = data.frame(coordGS,z.pred=rep(0,nrow(coordGS))), 
                       aes(x,y), color="black", size=0.2, alpha=0.5, inherit.aes = F) +   
            coord_equal() +
            theme_minimal() +
            # labs(fill=paste0("Assemblage ",k)) + ggtitle(letters[k]) +
            #           scale_x_continuous(limits=c(5,395), expand = c(0,0)) +
            #           scale_y_continuous(limits=c(5,295), expand = c(0,0)) + 
            theme(legend.position="bottom", legend.text=element_text(size=7), 
                  legend.title=element_text(size=8), axis.title=element_blank(), 
                  axis.text = element_blank(), panel.background = element_blank(),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
            guides(fill = guide_colorbar(barwidth = 8, barheight = 0.4, title.position="bottom"))
          
          #ggsave(filename = "Test2.pdf", tmp.one.plot.discrete, width = 10) 
        }  
      } else 
      {
        #plotting each topic separately as a color gradient
        tmp.plot[[k]] = qplot(x, y, data=spatial_topicmix_kriged[[k0]], geom="raster", fill=z.pred) +
          #tmp.plot[[k]] = qplot(x, y, spatial_topicmix_Blaise, geom="raster", fill=z.pred) +
          #tmp.plot[[k]] = qplot(x, y, data=spatial_topicmix_kriged, geom="raster", fill=z.pred) +
          #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", function(x) round(x,0))) + 
          #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", scientific_format())) +
          scale_fill_gradientn(colours=color.pal(7)) +
          coord_equal() + theme_minimal() +
          labs(fill=paste0("Assemblage ",k)) + ggtitle(letters[k]) +
          theme(legend.position="bottom", legend.text=element_text(size=7), 
                legend.title=element_text(size=8), axis.title=element_blank(), 
                axis.text = element_blank(),
                plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
          guides(fill = guide_colorbar(barwidth = 8, barheight = 0.4, title.position="bottom"))
        
        if (data_pp)
        {
          tmp.plot[[k]] = tmp.plot[[k]] + 
            scale_y_continuous(limits=c(5,395), expand = c(0,0)) +
            scale_x_continuous(limits=c(5,295), expand = c(0,0)) + 
            geom_point(data = data.frame(coordPP,z.pred=rep(0,nrow(coordPP))), aes(x,y), color="black", size=0.2, alpha=0.3)
        } else if (data_h20)
        {
          tmp.plot[[k]] = tmp.plot[[k]] + 
            scale_y_continuous(limits=c(5,95), expand = c(0,0)) +
            scale_x_continuous(limits=c(5,95), expand = c(0,0)) + 
            geom_point(data = data.frame(coordH20,z.pred=rep(0,nrow(coordH20))), aes(x,y), color="black", size=0.2, alpha=0.3)
        }
        
        if (k==1)
        {
          # Plotting each topic as a distinct color on a single map
          tmp.one.plot = ggplot(data = spatial_topicmix_kriged_all_topics_colors) +
            geom_raster(data = spatial_topicmix_kriged_all_topics_colors, aes(x, y), 
                        fill=rgb(spatial_topicmix_kriged_all_topics_colors[,3:5]/max(255,max(spatial_topicmix_kriged_all_topics_colors[,3:5]))), inherit.aes = F) +
            #geom_raster() +
            coord_equal() + theme_minimal() +
            #labs(fill=paste0("Assemblage ",k)) + ggtitle(letters[k]) +
            theme(legend.position="bottom", legend.text=element_text(size=7), 
                  legend.title=element_text(size=8), axis.title=element_blank(), 
                  axis.text = element_blank(),
                  plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
            guides(fill = guide_colorbar(barwidth = 8, barheight = 0.4, title.position="bottom"))
          
          # Plotting only the dominant topic as a distinct color in each pixel
          tmp.one.plot.discrete = ggplot(data = spatial_topicmix_kriged_all_topics_discrete_colors) +
            geom_raster(data = spatial_topicmix_kriged_all_topics_discrete_colors, aes(x, y),
                        fill=rgb(spatial_topicmix_kriged_all_topics_discrete_colors[,3:5]/max(255,max(spatial_topicmix_kriged_all_topics_discrete_colors[,3:5]))), inherit.aes = F) +
            #geom_raster() +
            coord_equal() + theme_minimal() +
            #labs(fill=paste0("Assemblage ",k)) + ggtitle(letters[k]) +
            theme(legend.position="bottom", legend.text=element_text(size=7), 
                  legend.title=element_text(size=8), axis.title=element_blank(), 
                  axis.text = element_blank(),
                  plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
            guides(fill = guide_colorbar(barwidth = 8, barheight = 0.4, title.position="bottom"))
          if (data_pp)
          {
            tmp.one.plot = tmp.one.plot + 
              scale_y_continuous(limits=c(5,395), expand = c(0,0)) +
              scale_x_continuous(limits=c(5,295), expand = c(0,0)) + 
              geom_point(data = data.frame(coordPP,z.pred=rep(0,nrow(coordPP))), aes(x,y), color="black", size=0.2, alpha=0.3, inherit.aes = F)
            
            tmp.one.plot.discrete = tmp.one.plot.discrete + 
              scale_y_continuous(limits=c(5,395), expand = c(0,0)) +
              scale_x_continuous(limits=c(5,295), expand = c(0,0)) + 
              geom_point(data = data.frame(coordPP,z.pred=rep(0,nrow(coordPP))), aes(x,y), color="black", size=0.2, alpha=0.3, inherit.aes = F)
          } else if (data_h20)
          {
            tmp.one.plot = tmp.one.plot + 
              scale_y_continuous(limits=c(5,95), expand = c(0,0)) +
              scale_x_continuous(limits=c(5,95), expand = c(0,0)) + 
              geom_point(data = data.frame(coordH20,z.pred=rep(0,nrow(coordH20))), aes(x,y), color="black", size=0.2, alpha=0.3)
            
            tmp.one.plot.discrete = tmp.one.plot.discrete + 
              scale_y_continuous(limits=c(5,95), expand = c(0,0)) +
              scale_x_continuous(limits=c(5,95), expand = c(0,0)) + 
              geom_point(data = data.frame(coordH20,z.pred=rep(0,nrow(coordH20))), aes(x,y), color="black", size=0.2, alpha=0.3)
          }
        }
      }
    } 
  }
  
  if (kriging && (j_select == kriged_real[1]) && (nb_topics == 3) && data_pp)
  { 
    lidar_prefix = paste0(local_prefix,data_insert,"/Lidar/")
    
    ########
    library(raster)
    
    r_topo <- raster(paste0(lidar_prefix,"r_topol_0.asc"))
    r_topo_transposed_flipped = flip(t(r_topo),'y')
    r_slope <- raster(paste0(lidar_prefix,"r_slopel_0.asc"))
    r_slope_transposed_flipped = flip(t(r_slope),'y')
    
    #TopoPP = as.data.frame(r_topo, row.names=NULL, optional=FALSE, xy=TRUE, na.rm=FALSE, long=FALSE)
    TopoPP = as.data.frame(r_topo_transposed_flipped, row.names=NULL, optional=FALSE, xy=TRUE, na.rm=FALSE, long=FALSE)
    SlopePP = as.data.frame(r_slope_transposed_flipped, row.names=NULL, optional=FALSE, xy=TRUE, na.rm=FALSE, long=FALSE)
    TopoPP_inverse = TopoPP
    TopoPP_inverse$layer = max(TopoPP_inverse$layer)-TopoPP_inverse$layer
    data_lidar_3topics = list()
    data_lidar_3topics[[1]] = TopoPP
    data_lidar_3topics[[2]] = TopoPP_inverse
    data_lidar_3topics[[3]] = SlopePP
    
    #         old_data_lidar_3topics = data_lidar_3topics
    #         newdata_lidar = expand.grid(x=TopoPP[1:420,1], y=rev(TopoPP[420*(1:320),2]))
    #         for (i in 1:nrow(newdata_lidar))
    #         {
    #           found = 0
    #           j = 1
    #           while (!found)
    #           {
    #             if (TopoPP[j,1]!=newdata_lidar[i,1] || TopoPP[j,2]!=newdata_lidar[i,2])
    #               j = j+1
    #             else
    #               found = 1
    #           }
    #           data_lidar_3topics[i,] = old_data_lidar_3topics[j,]   
    #         }
    
    map_labels = c("Topography","Inverse topography","Slope")
    lidar.plot = list()
    for (k in 1:3)
    {
      lidar.plot[[k]] = qplot(x, y, data = data_lidar_3topics[[k]], geom="raster", fill=layer) +
        #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", function(x) round(x,0))) + 
        #scale_fill_gradientn(colours=color.pal(7), labels=trans_format("identity", scientific_format())) +
        scale_fill_gradientn(colours=color.pal(7)) +
        coord_equal() + theme_minimal() +
        labs(fill=map_labels[k]) +  ggtitle(letters[k+3]) +
        #                   scale_y_continuous(limits=c(5,395), expand = c(0,0)) +
        #                   scale_x_continuous(limits=c(5,295), expand = c(0,0)) +
        scale_y_continuous(limits=c(5,395), expand = c(0,0)) +
        scale_x_continuous(limits=c(5,295), expand = c(0,0)) + 
        theme(legend.position="bottom", legend.text=element_text(size=7), 
              legend.title=element_text(size=8), axis.title=element_blank(), 
              axis.text = element_blank(),
              plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
        guides(fill = guide_colorbar(barwidth = 8, barheight = 0.4, title.position="bottom"))
      
      #                 if (data_pp)
      #                   tmp.plot[[k]] = tmp.plot[[k]] + geom_point(data = data.frame(coordPP,z.pred=rep(0,nrow(coordPP))), aes(x,y), color="black", size=0.2, alpha=0.3)
      #                 else if (data_h20)
      #                   tmp.plot[[k]] = tmp.plot[[k]] + geom_point(data = data.frame(coordH20,z.pred=rep(0,nrow(coordH20))), aes(x,y), color="black", size=0.2, alpha=0.3)
    }
    
    # setwd("/Users/guilhemsommeria-klein/Desktop/")
    # ggsave(filename = "Test.pdf", do.call("arrangeGrob", c(lidar.plot, ncol=3)), width = 10, height = 10/3*4/3)
  }
  
  return(list(lidar.plot = lidar.plot,
              tmp.plot = tmp.plot,
              tmp.one.plot = tmp.one.plot,
              tmp.one.plot.discrete = tmp.one.plot.discrete,
              spatial_topicmix = spatial_topicmix,
              spatial_topicmix_kriged = ifelse(!file.exists(paste0(subsubsubdirname,"Spatial_topicmix_kriged.rds")),spatial_topicmix_kriged,NA)))
  # spatial_topicmix_kriged is returned only if it had not been saved previously
  
}
