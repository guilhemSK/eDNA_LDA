library(plotrix)

# to execute : Rscript progname.R in shell

nb_topics=40
cat(nb_topics,"topics\n")
k_max = 4
cat(k_max,"runs\n\n")

setwd(paste("/Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/Blei_lda_alphaes0.1_topics",nb_topics,"_norep_strset_1/",sep=""))
logbeta1 <- read.table("final.beta", sep=" ")
logbeta1 = logbeta1[,-c(1,2)]

nb_terms=length(logbeta1[1,])

gamma <- read.table("final.gamma", sep=" ")
alpha <- read.table("final.other", sep=" ")

alpha_value = alpha[3,2]
#nb_terms = alpha[2,2]-1
nb_doc = length(gamma[,1])

documents = as.matrix(gamma-alpha_value)
documents[documents<0] = 0
# propotion of each topic in a document (sums to 1 over topics)
norm_documents = matrix(nrow=nb_doc,ncol=nb_topics)
for (i in 1:nb_doc)
{
  norm_documents[i,] = documents[i,]/sum(documents[i,])
}

# Computing the proportion of site-normalized read number for each topic
normal_topic = vector(length=nb_topics,mode="numeric")
for (j in 1:nb_topics)
{ 
  for (i in 1:nb_doc)
  {
    normal_topic[j] =  normal_topic[j] + norm_documents[i,j]/nb_doc 
  }
}
sort_normal_topic = sort.int(normal_topic,index.return=T)

#############
averaged_KL1 = matrix(nrow=nb_topics,ncol=nb_topics,data=0)
#averaged_KL2 = matrix(nrow=nb_topics,ncol=nb_topics,data=0)
for (k in 2:k_max)
  {
  cat("\n###########\nrun ",k-1,sep="")
  setwd(paste("/Users/guilhemsommeria-klein/Desktop/These/Données_H20/Plantes_GH/Blei_lda_alphaes0.1_topics",nb_topics,"_norep_strset_",k,"/",sep=""))
  logbeta <- read.table("final.beta", sep=" ")
  logbeta = logbeta[,-c(1,2)]
  KL = matrix(nrow=nb_topics,ncol=nb_topics,data=0)
  new_j1 = vector(length=nb_topics,mode="numeric")
  #new_j2 = vector(length=nb_topics,mode="numeric")
  for (i in 1:nb_topics)
    {
    cat("\n\n#i=",i,"\n\n",sep="")
    for (j in 1:nb_topics)
      {
      cat("j=",j,"\n",sep="")
      for (seq in 1:nb_terms)
        {
        #cat("seq=",seq,"\n",sep="")
        KL[i,j] = 1/(2*log(2))*(exp(logbeta1)[rev(sort_normal_topic$ix)[i],seq]*(logbeta1[rev(sort_normal_topic$ix)[i],seq]-logbeta[j,seq]) 
                                + exp(logbeta)[j,seq]*(logbeta[j,seq]-logbeta1[rev(sort_normal_topic$ix)[i],seq])) + KL[i,j]
        }
      }
    # First sorting : topics of first run are ranked by site-norm abundances, 
    # then a topic of the second run is assigned to each topic of the first run by following the ranking,
    # so that no 2nd run's topic is assigned twice to the same 1st run's topic
    sorting=sort.int(KL[i,],index.return=T)
    go_on = 1
    l = 1
    while (go_on)
      {
      j = 1
      break_var = 0
      while (j <= (i-1) && (!break_var))
        { 
        if (sorting$ix[l]==new_j1[j])
          break_var = 1
        else j = j+1
        }
      if (break_var)
        l = l+1
      else go_on = 0  
      }
    new_j1[i] = sorting$ix[l]
    }
  # Second sorting : all KL values are sorted, and second run's and first run's topics are paired following that ranking
#   sorting=sort.int(KL,index.return=T)
#   m = 1  
#   while (length(which(new_j2==0))!=0)
#     {
#     i = sorting$ix[m]%%nb_topics
#     j = (sorting$ix[m]-i)/nb_topics  
#     if (new_j2[i]==0)
#       {
#       l = 1
#       break_var = 0
#       while (l<=nb_topics && (!break_var))
#         { 
#         if (j==new_j2[l])
#           break_var = 1
#         else l = l+1
#         }
#       if (!break_var)
#         new_j2[i] = j
#       }
#     m = m+1
#     }
  
  for (i in 1:nb_topics)
    {
    for (j in 1:nb_topics)
      {
      averaged_KL1[i,j] = KL[i,new_j1[j]] + averaged_KL1[i,j]
      # averaged_KL2[i,j] = KL[i,new_j2[j]] + averaged_KL2[i,j]
      }
    }
  setwd(paste("/Users/guilhemsommeria-klein/Desktop/These/Données_H20/plots/Plantes_GH/run_to_run_Blei_lda_alphaes0.1_topics",nb_topics,"_norep_strset/",sep=""))
  
  file = paste("averaged_KL_divergence_between_",k,"runs_1stSorting.txt",sep="")
  write(averaged_KL1[1,]/(k-1),file,ncolumns=nb_topics)
  for (i in 2:nb_topics)
    write(averaged_KL1[i,]/(k-1),file,ncolumns=nb_topics,append = TRUE)
  
#   file = paste("averaged_KL_divergence_between_",k,"runs_2ndSorting.txt",sep="")
#   write(averaged_KL2[1,]/(k-1),file,ncolumns=nb_topics)
#   for (i in 2:nb_topics)
#     write(averaged_KL2[i,]/(k-1),file,ncolumns=nb_topics,append = TRUE)
#   
  pdf(paste("averaged_KL_divergence_between_",k,"runs_1stSorting.pdf",sep=""))
  par(mar=c(5.1,6.1,5.1,2.1))
  color2D.matplot(averaged_KL1/(k-1),c(1,0),c(0,0),c(0,1),
                  extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,
                  ylab="Classes ranked by\n relative abundances averaged over samples, first run",
                  xlab="Classes ranked by best match, following runs",
                  main=paste("Symmetrized KL divergence\n between samples' distributions\n over sequences for",k,"runs"),
                  do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
                  border="black",na.color=NA)
  dev.off()
  
#   pdf(paste("averaged_KL_divergence_between_",k,"runs_2ndSorting.pdf",sep=""))
#   color2D.matplot(averaged_KL2/(k-1),c(1,0),c(0,0),c(0,1),
#                   extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,
#                   ylab="Topics ranked by relative abundances averaged over samples for reference run",
#                   xlab="Topics of following runs",
#                   main=paste("Averaged KL symmetrized divergence\n between topics' distributions over sequences\n for",k,"runs"),
#                   do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
#                   border="black",na.color=NA)
#   dev.off()
  
  pdf(paste("averaged_KL_divergence_log_between_",k,"runs_1stSorting.pdf",sep=""))
  par(mar=c(5.1,6.1,5.1,2.1))
  color2D.matplot(log(averaged_KL1/(k-1))/log(10),c(1,0),c(0,0),c(0,1),
                  extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,
                  ylab="Classes ranked by\n relative abundances averaged over samples, first run",
                  xlab="Classes ranked by best match, following runs",
                  main=paste("Symmetrized KL divergence (log)\n between samples' distributions\n over sequences for",k,"runs"),
                  do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2,
                  border="black",na.color=NA)
  dev.off()
  
#   pdf(paste("averaged_KL_divergence_log_between_",k,"runs_2ndSorting.pdf",sep=""))
#   color2D.matplot(log(averaged_KL2/(k-1))/log(10),c(1,0),c(0,0),c(0,1),
#                   extremes=NA,cellcolors=NA,show.legend=TRUE,nslices=10,
#                   ylab="Topics ranked by relative abundances averaged over samples for reference run",
#                   xlab="Topics of following runs",
#                   main=paste("Averaged KL symmetrized divergence (log)\n between topics' distributions over sequences\n for",k,"runs"),
#                   do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,
#                   border="black",na.color=NA)
#   dev.off()
  }

# nb_topics = 40
# setwd(paste("/Users/guilhemsommeria-klein/Desktop/These/Données_H20/plots/Plantes_GH/run_to_run_Blei_lda_alphaes0.1_topics",nb_topics,"_norep_strset/",sep=""))
# 
# ##############################################
# 
# setwd("/Users/guilhemsommeria-klein/Desktop/These/Données_H20/plots/Plantes_GH/run_to_run_Blei_lda_alphaes0.1_topics40_norep_strset/")
# averaged_KL = read.table("averaged_KL_divergence_between_2runs_1stSorting.txt")
# averaged_KL = as.matrix(averaged_KL)
# 
# pdf("averaged_KL_divergence_(log)_between_2runs_1stSorting.pdf")
# par(mar=c(5.1,6.1,5.1,2.1)) 
# color2D.matplot(log(averaged_KL)/log(10),c(1,0),c(0,0),c(0,1),
#                 extremes=NA,cellcolors=NA,show.legend=F,nslices=10,
#                 ylab="Classes ranked by\n relative abundances averaged over samples, first run",
#                 xlab="Classes ranked by best match, second run",
#                 main="Symmetrized KL divergence (log)\n between samples' distributions\n over sequences for 2 runs",
#                 do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol="white",vcex=1,cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2,
#                 border="black",na.color=NA)
# dev.off()




  
  