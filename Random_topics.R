# Required to use the function rdirichlet :
library(gtools)

nb_samples = 1131
nb_topics=3
nb_motus=10000
nb_real = 100
dirich_alpha = 0.05

Hellinger = vector(mode="list",length=nb_real-1)
for (r in 1:nb_real)
{
  Sample_composition = matrix(nrow=nb_samples,ncol=nb_topics,data=0)
  # Random samples: 
#   for (j in 1:nb_samples)
#   {
#     Draw = sort(runif(nb_topics-1),decreasing=F)
#     Sample_composition[j,1] = Draw[1]
#     for (k in 2:(nb_topics-1))
#     {
#       Sample_composition[j,k] = Draw[k]-Draw[k-1]
#     } 
#     Sample_composition[j,nb_topics] = 1-Draw[nb_topics-1]
#   }
  # Dirichlet samples:
  for (j in 1:nb_samples)
  {
    Sample_composition[j,] = rdirichlet(1,rep(dirich_alpha,nb_topics)) 
  }
  KL_Sample_composition = matrix(nrow=nb_samples,ncol=nb_topics)
  for (j in 1:nb_samples)
  {
    for (k in 1:nb_topics)
      KL_Sample_composition[j,k] = Sample_composition[j,k]/sum(Sample_composition[,k])
  }
  if (r==1)
  {
    KL_Sample_composition1 = KL_Sample_composition
  } else 
  {
    Similarity = 1 - 1/sqrt(2)*dist(t(sqrt(cbind(KL_Sample_composition1,KL_Sample_composition))),
                                    method = "euclidean", diag = FALSE, upper = FALSE)
    Hellinger[[r-1]] = as.matrix(Similarity)[(nb_topics+1):(2*nb_topics),1:nb_topics]
  }
}

setwd("/Users/guilhemsommeria-klein/Desktop/These/Donnees_PetitPlateau/Test_data")

# Plotting one of the topic in one of the realization, to see how it looks like
#library(plotrix)
#pdf("Spatial distribution.pdf")
# spatial_topicmix = matrix(nrow=29,ncol=39)
# for (j in 1:39)
# {
#   for (i in 1:29)
#   {
#     spatial_topicmix[i,j] = Sample_composition[(j-1)*29+i,1]
#   }
# }
# par(cex.lab=1.5,cex.main=1.7,cex.axis=2)
# color2D.matplot(spatial_topicmix,c(226/255,2/255),c(226/255,63/255),c(226/255,165/255),
#                 extremes=NA,cellcolors=NA,show.legend=F,nslices=10,xlab="",xrange=c(0,1),
#                 ylab="",do.hex=FALSE,axes=FALSE,show.values=FALSE,vcol="white",vcex=1,
#                 border=NA,na.color="white")
# axis(2, at=c(0,5,10,15,20,25,30), labels=c("","50","100","150","200","250","300"))
# axis(1, at=c(0,5,10,15,20,25,30,35,40), labels=c("0","","100","","200","","300","","400"))
#dev.off()

Topic_comparison = matrix(nrow=nb_real-1,ncol=nb_topics)
for (r in 2:nb_real)
  Topic_comparison[r-1,] = apply(Hellinger[[r-1]],1,max)
#save(Topic_comparison,file="Topic_comparison.Rdata")

#mean(Topic_comparison[,1])
mean(apply(Topic_comparison,1,mean))

pdf("Topic_Hellinger_similarity_between_random_real_3t10000mot100r_dirichletsamples0.05.pdf")
par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
# par(mar = c(bottom, left, top, right))
par(mar = c(5, 5, 4, 3) + 0.1)
for (k in 1:nb_topics)
{
  if (k == 1)
    tag = "st"
  else if (k == 2)
    tag = "nd"
  else if (k == 3)
    tag = "rd"
  else tag = "th"
  plot(seq(2,nb_real,1),Topic_comparison[,k],ylim=c(0,1),type="p",main = paste("Hellinger similarity between the\n ",k,tag," component community in different realizations",sep=""),
       xlab = "Realizations ranked by decreasing likelihood", ylab = "Hellinger similarity")
}
dev.off()

pdf("Topic_Hellinger_similarity_between_real_averaged_3t10000mot100r_dirichletsamples0.05.pdf")  
par(cex.lab=1.5,cex.main=1.7,cex.axis=1.5,lwd=2)
# par(mar = c(bottom, left, top, right))
par(mar = c(5, 5, 4, 3) + 0.1)
plot(seq(2,nb_real,1),apply(Topic_comparison,1,mean),ylim=c(0,1),type="p", main = "Averaged Hellinger similarity between\n the component communities in different realizations",
     xlab = "Realizations ranked by decreasing likelihood", ylab = "Hellinger similarity")
dev.off()

