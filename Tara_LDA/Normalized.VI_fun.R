normalized_VI_fun = function(documents1,documents2,stations.means="different",SUR.DCM.all,stations_depths,coord)
{
  if (stations.means=="different")
  {
    Pk_1 = colMeans(documents1)
    Pk1_2 = colMeans(documents2)
  } 
  
  if (SUR.DCM.all %in% c("SUR","DCM"))
  {
    documents1 = documents1[rownames(documents1) %in% rownames(documents2),][stations_depths[rownames(coord) %in% rownames(documents1) & rownames(coord) %in% rownames(documents2),2] == SUR.DCM.all,]
    documents2 = documents2[rownames(documents2) %in% rownames(documents1),][stations_depths[rownames(coord) %in% rownames(documents1) & rownames(coord) %in% rownames(documents2),2] == SUR.DCM.all,]
  } else if (SUR.DCM.all == "all")
  {
    documents1 = documents1[rownames(documents1) %in% rownames(documents2),]
    documents2 = documents2[rownames(documents2) %in% rownames(documents1),]
  }
  
  if (nrow(documents1)>0 && nrow(documents2)>0)
  {
    if (stations.means=="same")
    {
      Pk_1 = colMeans(documents1)
      Pk1_2 = colMeans(documents2)
    }
    
    #########
    # Pkk1 = matrix(nrow = ncol(documents1), ncol = ncol(documents2), data = NA)
    # Pkk1_over_PkPk1 = matrix(nrow = ncol(documents1), ncol = ncol(documents2), data = NA)
    # for (k in 1:ncol(documents1))
    # {
    #   for (k1 in 1:ncol(documents2))
    #   {
    #     Pkk1[k,k1] = mean(documents1[,k]*documents2[,k1])
    #     Pkk1_over_PkPk1[k,k1] = Pkk1[k,k1]/(Pk_1[k]*Pk1_2[k1])
    #   }
    # }
    # I = sum(Pkk1*log(Pkk1_over_PkPk1))
    ########
    sum_Pkk1_times_log_Pkk1_over_PkPk1 = 0
    for (k in 1:ncol(documents1))
    {
      for (k1 in 1:ncol(documents2))
      {
        Pkk1 = mean(documents1[,k]*documents2[,k1])
        if (Pkk1>0)
          sum_Pkk1_times_log_Pkk1_over_PkPk1 = sum_Pkk1_times_log_Pkk1_over_PkPk1 + Pkk1*(log(Pkk1) - log(Pk_1[k]) - log(Pk1_2[k1]))
      }
    }
    # Computing the Variation of Information using probability values 
    H_1 = -sum(Pk_1*log(Pk_1))
    H_2 = -sum(Pk1_2*log(Pk1_2))
    I = sum_Pkk1_times_log_Pkk1_over_PkPk1
    Normalized_VI = (H_1 + H_2 - 2*I)/(H_1 + H_2 - I)
  } else
    Normalized_VI = NA
  
  return(Normalized_VI)
}