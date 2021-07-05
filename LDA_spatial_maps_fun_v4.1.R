# Variables values for testing :
# j_select = 1
# plotted_real = 1
# nb_topics = 10
# grid = 0
# geographic = 1
# documents = replicate(nb_topics,abs(rnorm(1131)))
# sort_normal_topic = sort.int(c(1,2,3),index.return=T)
# nb_real = 100
# em_tol = 10^-7
# var_tol = 10^-8
# Selected_real = 1:100
# #
# local_prefix = "/Users/guilhemsommeria-klein/Desktop/Post-doc_2017-2018/"
# filename_insert = "Rtopicmodels_LDA_VEM"
# data_insert = "Donnees_Tara"
# barcode_insert = "18S_V9_TARA_CompleteSizeRange_byStationByDepth_AllTaxa_noArcticNoBiomark_noLagoon"
# occurrence_insert = ""
# remove_single_sites_insert = ""
# no_rare_insert = ""
# occurrence_insert = "_occurrence"
# local_dirname = paste(local_prefix,data_insert,"/",barcode_insert,"/",filename_insert,"/",sep="")
# local_subdirname = paste0(local_dirname,filename_insert,"_nb_topics",nb_topics,"_nb_real",nb_real,"_em_tol",em_tol,"_var_tol",var_tol,"_best_keep",occurrence_insert,remove_single_sites_insert,no_rare_insert,"/")
# outputdirname = paste0(local_subdirname,"ordered_realization_number_",Selected_real[j_select],"/")
# setwd(outputdirname)

LDA_spatial_maps_fun <- function(grid,geographic,oneplot,
                                 documents,coord,coord0,
                                 backgroundGgplot,riversGgplot,land,bat,
                                 outputdirname,
                                 UnlimitedPointsPerStation,
                                 SurDCMperStation,
                                 legend_labels,
                                 col_range){
  
  nb_topics = ncol(documents)
  
  # if (!file.exists(paste0(outputdirname,"Spatial_topicmix_kriged.rds")))
  # {
    # newdata contains the coordinates of all pixels in the kriged map, with grain pixels between two neighbouring sampling sites
    if (grid)
    {
      if (byRow)
        newdata = expand.grid(x = seq(min(coord$x) - grain, max(coord$x) + grain, 1), y = seq(min(coord$y) - grain, max(coord$y) + grain, 1))
      else
        newdata = expand.grid(y = seq(min(coord$y) - grain, max(coord$y) + grain, 1), x = seq(min(coord$x) - grain, max(coord$x) + grain, 1))[,c(2,1)]
    }
    
    spatial_topicmix_kriged = vector(length=nb_topics,mode="list")
    # Loop over topics
    for (k in 1:nb_topics)
    {   
      if (grid)
      {
        # Building a kriged spatial map of the grid, topic by topic
        tmp = data.frame(x=coord$x,y=coord$y,z=documents[,k])
        mod = gstat(id = "z", formula = z~1, locations = ~x+y, data=tmp, model=vgm(10, "Exp", 20))
        spatial_topicmix_kriged[[k]] = predict(mod, newdata, na.action=na.omit)
        # Removing interpolated relative abundances falling outside of [0,1]
        # Indeed, the sum of all assemblages may then not be 1 for all pixels
        spatial_topicmix_kriged[[k]]$z.pred[spatial_topicmix_kriged[[k]]$z.pred<0] = 0
        spatial_topicmix_kriged[[k]]$z.pred[spatial_topicmix_kriged[[k]]$z.pred>1] = 1
      } else if (!grid && geographic)
      {
        # For !grid && geographic, spatial_topicmix_kriged is actually not kriged, but the name of the variable is kept for consistency with the remaining of the pipeline
        spatial_topicmix_kriged[[k]] = data.frame(x=coord$x,y=coord$y,z.pred=documents[,k])
        rownames(spatial_topicmix_kriged[[k]]) = rownames(coord)
      }
    } 
  # } else 
  #   spatial_topicmix_kriged = readRDS(paste0(outputdirname,"Spatial_topicmix_kriged.rds"))
  
  # Storing all assemblages into a single dataframe:
  spatial_topicmix_kriged_all_topics = data.frame(x = spatial_topicmix_kriged[[1]]$x, y = spatial_topicmix_kriged[[1]]$y)
  rownames(spatial_topicmix_kriged_all_topics) = rownames(spatial_topicmix_kriged[[1]])
  for (k in 1:nb_topics)
  {
    spatial_topicmix_kriged[[k]]$z.pred[spatial_topicmix_kriged[[k]]$z.pred < 0.001] = 0
    spatial_topicmix_kriged_all_topics = cbind(spatial_topicmix_kriged_all_topics,spatial_topicmix_kriged[[k]]$z.pred)
    colnames(spatial_topicmix_kriged_all_topics)[2+k] = paste0("z.pred",k)
  }
  
  # Defining the color palette:
  # color.pal = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  color.pal = colorRampPalette(c("#00007F", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),space = "Lab")
  
  if (oneplot && grid)
  {
    # Discrete spatial distribution: assigning a single assemblage (the most abundant) to each location
    spatial_topicmix_kriged_all_topics_discrete = spatial_topicmix_kriged_all_topics
    spatial_topicmix_kriged_all_topics_discrete[,3:(2+nb_topics)] = 0
    dominant_topic_proportion = apply(spatial_topicmix_kriged_all_topics[,3:(2+nb_topics)],1,max)
    dominant_topic_index = apply((spatial_topicmix_kriged_all_topics[,3:(2+nb_topics)] == dominant_topic_proportion),1,which)
    # dominant_topic_index is a list if two topics have equal proportions in one site, otherwise it is a vector  
    for (i in 1:nrow(spatial_topicmix_kriged_all_topics))
      spatial_topicmix_kriged_all_topics_discrete[i,2+dominant_topic_index[[i]]] = 1
    if (nb_topics == 3)
    {
      # Using a red-green-blue color code fro 3 topics, as for the Petit Plateau data
      col_matrix = matrix(data=c(c(0,0,255),c(0,255,0),c(255,0,0)),ncol=3)
      rownames(col_matrix) = c("red","green","blue")
    } else
    {
      col_matrix = col2rgb(color.pal(nb_topics))
    }
    # Defining the color of each pixel based on the assemblage composition:
    spatial_topicmix_kriged_all_topics_colors = data.frame(spatial_topicmix_kriged_all_topics[,1:2],red=0,green=0,blue=0)
    coord_matrix = t(as.matrix(spatial_topicmix_kriged_all_topics[,3:(2+nb_topics)]))
    spatial_topicmix_kriged_all_topics_discrete_colors = data.frame(spatial_topicmix_kriged_all_topics[,1:2],red=0,green=0,blue=0)
    coord_matrix_discrete = t(as.matrix(spatial_topicmix_kriged_all_topics_discrete[,3:(2+nb_topics)]))
    spatial_topicmix_kriged_all_topics_colors[,3:5] = t(col_matrix %*% coord_matrix)
    spatial_topicmix_kriged_all_topics_discrete_colors[,3:5] = t(col_matrix %*% coord_matrix_discrete)
  }
  
  if (geographic && !grid)
  {
    if (SurDCMperStation)
    {
      stations_depths = as.data.frame(matrix(nrow=nrow(coord),ncol=2))
      colnames(stations_depths) = c("Station","Depth")
      for (station_depth_index in 1:nrow(coord))
      {
        stations_depths[station_depth_index,] = c(strsplit(rownames(coord)[station_depth_index],split=" ")[[1]][1],
                                                  strsplit(rownames(coord)[station_depth_index],split=" ")[[1]][2])
      }
      stations_names = levels(as.factor(stations_depths$Station))
    } else if (UnlimitedPointsPerStation)
    {
      sample_ref = read.table(paste0(local_prefix,data_insert,"/TV9_corr_st_lat_long_nocross.txt"), colClasses="vector", sep=",")
      # stations_names stores stations ordered as in sample_ref
      stations_names = levels(as.factor(sample_ref$V11[-1]))
      #stations_names = stations_names[-c(which(stations_names=="TARA_017"),which(stations_names=="TARA_062"))]
      #stations_names = stations_names[-which(stations_names=="TARA_017")]
    }
    
    if (SurDCMperStation)
    {
      blues = c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
      greys = c(grey(0.6), grey(0.93), grey(0.99))
      
      pdf("Spatial_distribution_maps_SUR_onebyone.pdf")
      for (k in 1:nb_topics)
      {
        # Plotting Surface:
        plot(bat, image = T, land = T, lty = 0, bty = "n", xaxt = "n", yaxt = "n", ann = F, bpal = list(c(0, max(bat), greys[1]), c(min(bat), 0, blues))) #plot map without isobaths
        plot(bat, lwd = 0.8, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline
        points(cbind(spatial_topicmix_kriged_all_topics$x[stations_depths[,2] == "SUR"],spatial_topicmix_kriged_all_topics$y[stations_depths[,2] == "SUR"]),
               pch = 21, cex=1.2, bg = rgb(colorRamp(c("darkblue","firebrick2"),space = "Lab")(spatial_topicmix_kriged_all_topics[stations_depths[,2] == "SUR",2+k]),maxColorValue = 255))
      }
      dev.off()
      # pdf("AllTaxa_GibbsAlpha0.1Best100r10t_DCMplots_topicByTopic_marmap.pdf")
      pdf("Spatial_distribution_maps_DCM_onebyone.pdf")
      for (k in 1:nb_topics)
      {
        # Plotting DCM:
        plot(bat, image = T, land = T, lty = 0, bty = "n", xaxt = "n", yaxt = "n", ann = F, bpal = list(c(0, max(bat), greys[1]), c(min(bat), 0, blues))) #plot map without isobaths
        plot(bat, lwd = 0.8, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline
        points(cbind(spatial_topicmix_kriged_all_topics$x[stations_depths[,2] == "DCM"],spatial_topicmix_kriged_all_topics$y[stations_depths[,2] == "DCM"]),
               pch = 21, cex=1.2, bg = rgb(colorRamp(c("darkblue","firebrick2"),space = "Lab")(spatial_topicmix_kriged_all_topics[stations_depths[,2] == "DCM",2+k]),maxColorValue = 255))
      }
      dev.off()
      
      if (oneplot)
      {
        present_topics = apply(spatial_topicmix_kriged_all_topics[,3:(2+nb_topics)],1,function(g) which(g>0))
        pure_stations = unlist(lapply(present_topics,function(g) length(g))) == 1
        
        pdf("Spatial_distribution_maps_SUR.DCM_oneplot.pdf")
        # Plotting Surface:
        plot(bat, image = T, land = T, lty = 0, bty = "n", xaxt = "n", yaxt = "n", ann = F, bpal = list(c(0, max(bat), greys[1]), c(min(bat), 0, blues))) #plot map without isobaths
        plot(bat, lwd = 0.8, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline
        if (length(which(stations_depths[,2] == "SUR" & pure_stations)) > 0)
          points(cbind(spatial_topicmix_kriged_all_topics$x[stations_depths[,2] == "SUR" & pure_stations],spatial_topicmix_kriged_all_topics$y[stations_depths[,2] == "SUR" & pure_stations]),
               pch = 21, cex=1.2, bg = color.pal(nb_topics)[unlist(present_topics[stations_depths[,2] == "SUR" & pure_stations])])
        if (length(which(stations_depths[,2] == "SUR" & !pure_stations)) > 0)
        {
          pie.colors.matrix = data.frame(matrix(nrow = length(which(stations_depths[,2] == "SUR" & !pure_stations)), ncol=nb_topics, data = NA))
          pie.slices.matrix = data.frame(matrix(nrow = length(which(stations_depths[,2] == "SUR" & !pure_stations)), ncol=nb_topics, data = NA))
          ii_pie = 0
          for (i_pie in which(stations_depths[,2] == "SUR" & !pure_stations))
          {
            ii_pie = ii_pie+1
            station_topics = which(spatial_topicmix_kriged_all_topics[i_pie,3:(2+nb_topics)] > 0)
            pie.colors.matrix[ii_pie,] = c(color.pal(nb_topics)[station_topics],rep(NA,nb_topics-length(station_topics)))
            pie.slices.matrix[ii_pie,] = c(spatial_topicmix_kriged_all_topics[i_pie,2+station_topics],rep(NA,nb_topics-length(station_topics)))
          }
          space.pies(spatial_topicmix_kriged_all_topics$x[stations_depths[,2] == "SUR" & !pure_stations], spatial_topicmix_kriged_all_topics$y[stations_depths[,2] == "SUR" & !pure_stations],
                     pie.slices = pie.slices.matrix, pie.colors = pie.colors.matrix, 
                     pie.radius=3, pie.space=0.01,
                     link=TRUE, seg.lwd=0.5, seg.col=1, seg.lty=1, coord=NULL)
        }
        # Plotting DCM:
        plot(bat, image = T, land = T, lty = 0, bty = "n", xaxt = "n", yaxt = "n", ann = F, bpal = list(c(0, max(bat), greys[1]), c(min(bat), 0, blues))) #plot map without isobaths
        plot(bat, lwd = 0.8, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline
        if (length(which(stations_depths[,2] == "DCM" & pure_stations)) > 0)
          points(cbind(spatial_topicmix_kriged_all_topics$x[stations_depths[,2] == "DCM" & pure_stations],spatial_topicmix_kriged_all_topics$y[stations_depths[,2] == "DCM" & pure_stations]),
               pch = 21, cex=1.2, bg = color.pal(nb_topics)[unlist(present_topics[stations_depths[,2] == "DCM" & pure_stations])])
        if (length(which(stations_depths[,2] == "DCM" & !pure_stations)) > 0)
        {
          pie.colors.matrix = data.frame(matrix(nrow = length(which(stations_depths[,2] == "DCM" & !pure_stations)), ncol=nb_topics, data = NA))
          pie.slices.matrix = data.frame(matrix(nrow = length(which(stations_depths[,2] == "DCM" & !pure_stations)), ncol=nb_topics, data = NA))
          ii_pie = 0
          for (i_pie in which(stations_depths[,2] == "DCM" & !pure_stations))
          {
            ii_pie = ii_pie+1
            station_topics = which(spatial_topicmix_kriged_all_topics[i_pie,3:(2+nb_topics)] > 0)
            pie.colors.matrix[ii_pie,] = c(color.pal(nb_topics)[station_topics],rep(NA,nb_topics-length(station_topics)))
            pie.slices.matrix[ii_pie,] = c(spatial_topicmix_kriged_all_topics[i_pie,2+station_topics],rep(NA,nb_topics-length(station_topics)))
          }
          space.pies(spatial_topicmix_kriged_all_topics$x[stations_depths[,2] == "DCM" & !pure_stations], spatial_topicmix_kriged_all_topics$y[stations_depths[,2] == "DCM" & !pure_stations],
                     pie.slices = pie.slices.matrix, pie.colors = pie.colors.matrix, 
                     pie.radius=3, pie.space=0.01,
                     link=TRUE, seg.lwd=0.5, seg.col=1, seg.lty=1, coord=NULL)
        }
        dev.off()
        
        tmp.one.plot.discrete = NA
        tmp.one.plot = NA
      }
    }
  }
  
  tmp.plot = list()
  for (k in 1:nb_topics)
  {
    if (geographic && !grid) 
    {
      if (!UnlimitedPointsPerStation && !SurDCMperStation)
      {
        # Plotting each topic separately as a color gradient from blue (absent) to red (100% abundance)
        if (!is.na(backgroundGgplot))
        {
          tmp.plot[[k]] = ggplot(data = backgroundGgplot) +
                          geom_path(aes(x=long, y=lat, group=group), color = "black", size = 0.3, inherit.aes=T)
        } else
          tmp.plot[[k]] = ggplot(data = spatial_topicmix_kriged[[k]])
        
        tmp.plot[[k]] =  tmp.plot[[k]] +
          geom_point(data = spatial_topicmix_kriged[[k]], aes(x,y,colour=z.pred), size=2, alpha=1, inherit.aes = F) +
          scale_colour_gradientn(colours = c("darkblue","firebrick2"), space = "Lab") +
          coord_equal() +
          # scale_y_continuous(limits=range(coord$y), expand = c(0,0)) +
          # scale_x_continuous(limits=range(coord$x), expand = c(0,0)) +
          labs(fill=legend_labels[k]) +  theme_minimal() + ggtitle(letters[k]) +
          theme(legend.position="bottom", legend.text=element_text(size=7), 
                legend.title=element_text(size=8), axis.title=element_blank(), 
                axis.text = element_blank(), panel.background = element_blank(),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
          guides(fill = guide_colorbar(barwidth = 8, barheight = 0.4, title.position="bottom"))
        
        if (is.list(riversGgplot))
          tmp.plot[[k]] = tmp.plot[[k]]+
            geom_path(data = riversGgplot, aes(x=long*grain, y=lat*grain, group=group), 
                      color = "blue", size=0.4, inherit.aes = F, alpha=1)
      } else if (UnlimitedPointsPerStation)
      {
        tmp.plot[[k]] = ggplot(data = backgroundGgplot, aes(x=long, y=lat, group=group)) + geom_path(color = "black", size=0.1, inherit.aes = T)
        z.pred.k = spatial_topicmix_kriged[[k]]$z.pred
        
        for (station in stations_names)
        {
            #station_index = station_index+1
            # Looking for the samples corresponding to the current station in sample_ref:
            samples_indices = which(sample_ref$V11[-1]==station)
            # Storing the names of the samples corresponding to the current station, so as to look for them in data2m:
            samples_names = sample_ref$V1[samples_indices+1]
            # Storing the indices of those samples in data2m
            samples_indices_data2m = vector(length=length(samples_indices),mode="numeric")
            sample_index = 0
            # Loop over all samples cooresponding to the current station
            for (sample in samples_names)
            {
              sample_index = sample_index+1
              samples_indices_data2m[sample_index] = which(colnames_data2m==sample)
            }

          pie = ggplot(data.frame(x=z.pred.k[samples_indices_data2m],y=rep(1,length(samples_indices_data2m))), aes(x=1, y, fill=x))
          
          if (col_range == "natural")
            range_int = range(z.pred.k[samples_indices_data2m])
          else if (col_range == "standard")
            range_int = c(0,1)
          pie = pie + 
            #             scale_fill_brewer(palette="color.pal.heat") +
            #             scale_fill_brewer(palette="heat.colors") +
            scale_fill_gradient(limits = range_int, low="darkblue", high="firebrick2", na.value="transparent") +
            #scale_fill_gradientn(colours=color.pal.heat(2)) +
            #             scale_fill_manual(values=z.pred.k[samples_indices_data2m]) +
            geom_bar(stat="identity", width=1) + 
            coord_polar(theta="y",start=3*pi/2) + theme_tree() + 
            xlab(NULL) + ylab(NULL) + 
            theme_transparent()
          #           tmp.one.plot.piechart %<>% subview(tmp.one.plot.piechart,pie,coordTara$Long[samples_indices_data2m[1]],coordTara$Lat[samples_indices_data2m[1]],width=0.1, height=0.1)
          tmp.plot[[k]] = geom_subview(tmp.plot[[k]],pie,coord$x[samples_indices_data2m[1]],coord$y[samples_indices_data2m[1]],width=0.09, height=0.09)        
        }
        
        tmp.plot[[k]] = tmp.plot[[k]] +
          coord_equal() +
          #labs(fill=paste0("Assemblage ",k)) +  
          theme_minimal() + ggtitle(legend_labels[k]) +
          #           scale_x_continuous(limits=c(5,395), expand = c(0,0)) +
          #           scale_y_continuous(limits=c(5,295), expand = c(0,0)) + 
          theme(legend.position="bottom", legend.text=element_text(size=7), 
                legend.title=element_text(size=8), axis.title=element_blank(), 
                axis.text = element_blank(), panel.background = element_blank(),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
          guides(fill = guide_colorbar(barwidth = 8, barheight = 0.4, title.position="bottom"))
      }
    } else if (grid && geographic)
    {
      #plotting each topic separately as a color gradient 
      if (isS4(land))
      {
        df.plot = data.frame(x=spatial_topicmix_kriged[[k]]$x/grain,y=spatial_topicmix_kriged[[k]]$y/grain,z.pred=spatial_topicmix_kriged[[k]]$z.pred)
        df.plot.raster = rasterFromXYZ(df.plot)
        df.plot.extruded = extract(df.plot.raster,land,df=T,cellnumbers=T)
        df.plot.extruded.data.frame = data.frame(x=spatial_topicmix_kriged_all_topics_colors$x[df.plot.extruded$cell],
                                                 y=spatial_topicmix_kriged_all_topics_colors$y[df.plot.extruded$cell],df.plot.extruded[,3])
        df.plot.extruded.raster = rasterFromXYZ(df.plot.extruded.data.frame)
        df.plot.extruded.raster = flip(df.plot.extruded.raster,'y')
        df.plot.extruded.data.frame = as.data.frame(df.plot.extruded.raster, row.names=NULL, optional=T, xy=TRUE, na.rm=T, long=FALSE)
        df.plot = data.frame(x=df.plot.extruded.data.frame$x,y=df.plot.extruded.data.frame$y,z.pred=df.plot.extruded.data.frame[,3])
      } else
        df.plot = data.frame(x=spatial_topicmix_kriged[[k]]$x,y=spatial_topicmix_kriged[[k]]$y,z.pred=spatial_topicmix_kriged[[k]]$z.pred)
      
      tmp.plot[[k]] = ggplot(data = backgroundGgplot) +
        geom_raster(data = df.plot, aes(x, y, fill=z.pred), inherit.aes = F) +
        scale_fill_gradientn(colours=color.pal(7)) +
        geom_path(data = backgroundGgplot, aes(x=long*grain, y=lat*grain, group=group), 
                  color = "black", size=0.1, inherit.aes = F) +
        scale_y_continuous(limits=range(coord$y), expand = c(0,0)) +
        scale_x_continuous(limits=range(coord$x), expand = c(0,0)) +
        geom_point(data = data.frame(coord,z.pred=rep(0,nrow(coord))), 
                   aes(x,y), color="black", size=0.2, alpha=0.5, inherit.aes = F) +   
        coord_equal() +
        labs(fill=legend_labels[k]) +  theme_minimal() + ggtitle(letters[k]) +
        theme(legend.position="bottom", legend.text=element_text(size=7), 
              legend.title=element_text(size=8), axis.title=element_blank(), 
              axis.text = element_blank(), panel.background = element_blank(),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
        guides(fill = guide_colorbar(barwidth = 8, barheight = 0.4, title.position="bottom"))
      
      if (is.list(riversGgplot))
        tmp.plot[[k]] = tmp.plot[[k]]+
        geom_path(data = riversGgplot, aes(x=long*grain, y=lat*grain, group=group), 
                  color = "blue", size=0.4, inherit.aes = F, alpha=1)
      
      if (oneplot && k==1)
      {
        # Plotting each topic as a distinct color on a single map
        if (isS4(land))
        {
          df.plot = data.frame(x=spatial_topicmix_kriged_all_topics_colors$x/grain,y=spatial_topicmix_kriged_all_topics_colors$y/grain,spatial_topicmix_kriged_all_topics_colors[,3:5])
          df.plot.raster = rasterFromXYZ(df.plot)
          df.plot.extruded = extract(df.plot.raster,land,df=T,cellnumbers=T)
          df.plot.extruded.data.frame = data.frame(x=spatial_topicmix_kriged_all_topics_colors$x[df.plot.extruded$cell],
                                                   y=spatial_topicmix_kriged_all_topics_colors$y[df.plot.extruded$cell],df.plot.extruded[,3:5])
          df.plot.extruded.raster = rasterFromXYZ(df.plot.extruded.data.frame)
          df.plot.extruded.raster = flip(df.plot.extruded.raster,'y')
          df.plot.extruded.data.frame = as.data.frame(df.plot.extruded.raster, row.names=NULL, optional=T, xy=TRUE, na.rm=T, long=FALSE)
          df.plot = data.frame(x=df.plot.extruded.data.frame$x,y=df.plot.extruded.data.frame$y,df.plot.extruded.data.frame[,3:5])
        } else
          df.plot = data.frame(x=spatial_topicmix_kriged_all_topics_colors$x,y=spatial_topicmix_kriged_all_topics_colors$y,spatial_topicmix_kriged_all_topics_colors[,3:5])
        
        tmp.one.plot = ggplot(data = backgroundGgplot) +
          geom_raster(data = df.plot, aes(x, y), fill=rgb(df.plot[,3:5]/max(255,max(df.plot[,3:5]))), inherit.aes = F) +
          geom_path(data = backgroundGgplot, aes(x=long*grain, y=lat*grain, group=group), 
                    color = "black", size=0.1, inherit.aes = F) +
          scale_y_continuous(limits=range(coord$y), expand = c(0,0)) +
          scale_x_continuous(limits=range(coord$x), expand = c(0,0)) +
          geom_point(data = data.frame(coord,z.pred=rep(0,nrow(coord))), 
                     aes(x,y), color="black", size=0.2, alpha=0.5, inherit.aes = F) +   
          coord_equal() +
          theme_minimal() +
          theme(legend.position="bottom", legend.text=element_text(size=7), 
                legend.title=element_text(size=8), axis.title=element_blank(), 
                axis.text = element_blank(), panel.background = element_blank(),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
          guides(fill = guide_colorbar(barwidth = 8, barheight = 0.4, title.position="bottom"))
        
        if (is.list(riversGgplot))
          tmp.one.plot = tmp.one.plot +
          geom_path(data = riversGgplot, aes(x=long*grain, y=lat*grain, group=group), 
                    color = "blue", size=0.4, inherit.aes = F, alpha=1)
        
        # Plotting only the dominant topic as a distinct color in each pixel
        if (isS4(land))
        {
          df.plot = data.frame(x=spatial_topicmix_kriged_all_topics_colors$x/grain,y=spatial_topicmix_kriged_all_topics_colors$y/grain,spatial_topicmix_kriged_all_topics_discrete_colors[,3:5])
          df.plot.raster = rasterFromXYZ(df.plot)
          df.plot.extruded = extract(df.plot.raster,land,df=T,cellnumbers=T)
          df.plot.extruded.data.frame = data.frame(x=spatial_topicmix_kriged_all_topics_colors$x[df.plot.extruded$cell],
                                                   y=spatial_topicmix_kriged_all_topics_colors$y[df.plot.extruded$cell],
                                                   df.plot.extruded[,3:5])
          df.plot.extruded.raster = rasterFromXYZ(df.plot.extruded.data.frame)
          df.plot.extruded.raster = flip(df.plot.extruded.raster,'y')
          df.plot.extruded.data.frame = as.data.frame(df.plot.extruded.raster, row.names=NULL, optional=T, xy=TRUE, na.rm=T, long=FALSE)
          df.plot = data.frame(x=df.plot.extruded.data.frame$x, y=df.plot.extruded.data.frame$y, df.plot.extruded.data.frame[,3:5])
        } else
          df.plot = data.frame(x=spatial_topicmix_kriged_all_topics_colors$x,y=spatial_topicmix_kriged_all_topics_colors$y,spatial_topicmix_kriged_all_topics_discrete_colors[,3:5])
        
        tmp.one.plot.discrete = ggplot(data = backgroundGgplot) +
          geom_raster(data = df.plot, aes(x, y), 
                      fill=rgb(df.plot[,3:5]/max(255,max(df.plot[,3:5]))), inherit.aes = F) +
          geom_path(data = backgroundGgplot, aes(x=long*grain, y=lat*grain, group=group), 
                    color = "black", size=0.1, inherit.aes = F) +
          scale_y_continuous(limits=range(coord$y), expand = c(0,0)) +
          scale_x_continuous(limits=range(coord$x), expand = c(0,0)) +
          geom_point(data = data.frame(coord,z.pred=rep(0,nrow(coord))), 
                     aes(x,y), color="black", size=0.2, alpha=0.5, inherit.aes = F) +   
          coord_equal() +
          theme_minimal() +
          theme(legend.position="bottom", legend.text=element_text(size=7), 
                legend.title=element_text(size=8), axis.title=element_blank(), 
                axis.text = element_blank(), panel.background = element_blank(),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
          guides(fill = guide_colorbar(barwidth = 8, barheight = 0.4, title.position="bottom"))
        
        if (is.list(riversGgplot))
          tmp.one.plot.discrete = tmp.one.plot.discrete +
          geom_path(data = riversGgplot, aes(x=long*grain, y=lat*grain, group=group), 
                    color = "blue", size=0.4, inherit.aes = F, alpha=1)
      }  
    } else if (grid && !geographic)
    {
      # coord0$y = coord0$y+2/(max(coord0$y)-1)*(1-coord0$y/grain)+1
      # coord0$x = coord0$x+2/(max(coord0$x)-1)*(1-coord0$x/grain)+1
      
      # Plotting each topic separately as a color gradient
      tmp.plot[[k]] = qplot(x, y, data=spatial_topicmix_kriged[[k]], geom="raster", fill=z.pred) +
        scale_fill_gradientn(colours=color.pal(7)) +
        coord_equal() + theme_minimal() +
        labs(fill=legend_labels) + ggtitle(letters[k]) +
        theme(legend.position="bottom", legend.text=element_text(size=7), 
              legend.title=element_text(size=8), axis.title=element_blank(), 
              axis.text = element_blank(),
              plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
        guides(fill = guide_colorbar(barwidth = 8, barheight = 0.4, title.position="bottom")) +
        scale_y_continuous(limits = c(min(coord$y) - floor(grain/2), max(coord$y) + floor(grain/2)), expand = c(0,0)) +
        scale_x_continuous(limits = c(min(coord$x) - floor(grain/2), max(coord$x) + floor(grain/2)), expand = c(0,0)) + 
        geom_point(data = data.frame(coord0+1,z.pred=rep(0,nrow(coord0))), aes(x,y), color="black", size=0.5, alpha=0.3)
      
      if (oneplot && k==1)
      {
        # Plotting each topic as a distinct color on a single map
        tmp.one.plot = ggplot(data = spatial_topicmix_kriged_all_topics_colors) +
          geom_raster(data = spatial_topicmix_kriged_all_topics_colors, aes(x, y), 
                      fill=rgb(spatial_topicmix_kriged_all_topics_colors[,3:5]/max(255,max(spatial_topicmix_kriged_all_topics_colors[,3:5]))), inherit.aes = F) +
          coord_equal() + theme_minimal() +
          theme(legend.position="bottom", legend.text=element_text(size=7), 
                legend.title=element_text(size=8), axis.title=element_blank(), 
                axis.text = element_blank(),
                plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
          guides(fill = guide_colorbar(barwidth = 8, barheight = 0.4, title.position="bottom")) +
          scale_y_continuous(limits = c(min(coord$y) - floor(grain/2), max(coord$y) + floor(grain/2)), expand = c(0,0)) +
          scale_x_continuous(limits = c(min(coord$x) - floor(grain/2), max(coord$x) + floor(grain/2)), expand = c(0,0)) + 
          geom_point(data = data.frame(coord0+1,z.pred=rep(0,nrow(coord0))), aes(x,y), color="black", size=1.5, alpha=0.3, inherit.aes = F)
        
        # Plotting only the dominant topic as a distinct color in each pixel
        tmp.one.plot.discrete = ggplot(data = spatial_topicmix_kriged_all_topics_discrete_colors) +
          geom_raster(data = spatial_topicmix_kriged_all_topics_discrete_colors, aes(x, y),
                      fill=rgb(spatial_topicmix_kriged_all_topics_discrete_colors[,3:5]/max(255,max(spatial_topicmix_kriged_all_topics_discrete_colors[,3:5]))), inherit.aes = F) +
          coord_equal() + theme_minimal() +
          theme(legend.position="bottom", legend.text=element_text(size=7), 
                legend.title=element_text(size=8), axis.title=element_blank(), 
                axis.text = element_blank(),
                plot.title=element_text(hjust=0), plot.margin=unit(c(0,1,-2,2),"mm")) +
          guides(fill = guide_colorbar(barwidth = 8, barheight = 0.4, title.position="bottom")) +
          scale_y_continuous(limits = c(min(coord$y) - floor(grain/2), max(coord$y) + floor(grain/2)), expand = c(0,0)) +
          scale_x_continuous(limits = c(min(coord$x) - floor(grain/2), max(coord$x) + floor(grain/2)), expand = c(0,0)) + 
          geom_point(data = data.frame(coord0,z.pred=rep(0,nrow(coord0))), aes(x,y), color="black", size=1.5, alpha=0.3, inherit.aes = F)
      }
    }
  } 
  
  if (oneplot)
    return(list(tmp.plot = tmp.plot,
              tmp.one.plot = tmp.one.plot,
              tmp.one.plot.discrete = tmp.one.plot.discrete,
              spatial_topicmix_kriged = spatial_topicmix_kriged))
  else 
    return(list(tmp.plot = tmp.plot,
                spatial_topicmix_kriged = spatial_topicmix_kriged))
}
