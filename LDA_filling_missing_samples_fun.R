
LDA_filling_missing_samples_fun <- function(Missing_positions_indices,data2m,nrow_grid,ncol_grid,byRow){
   
spatial_data2m = matrix(nrow=nrow(data2m),ncol=nrow_grid*ncol_grid,data=0)
spatial_data2m[,which(Missing_positions_indices==0)] = data2m

# If samples are listed by row in data2m and Missing_positions_indices, they need to be reordered since LDA_filling_missing_samples_fun assumes samples to be listed by column: 
if (byRow)
{ 
  spatial_data2m = spatial_data2m[,as.vector(t(matrix(1:(nrow_grid*ncol_grid),ncol=nrow_grid)))]
  Missing_positions_indices = Missing_positions_indices[as.vector(t(matrix(1:(nrow_grid*ncol_grid),ncol=nrow_grid)))]
}

# Looping over all sites of the nrow_grid x ncol_grid sampling grid, and filling the missing samples based on the samples at neighbouring sites
# Grid edges are treated separately
data2m_filled = spatial_data2m
for (i in 1:nrow_grid)
{
  for (j in 1:ncol_grid)
  {
    if (Missing_positions_indices[(j-1)*nrow_grid+i] == 1)
    {
      # upper left corner
      if (i==1 && j==1)
      {
        # Looking for the non-empty neighbouring sites
        nonempty_neighb_indices = which(c(Missing_positions_indices[(j-1)*nrow_grid+i+1],Missing_positions_indices[j*nrow_grid+i],Missing_positions_indices[j*nrow_grid+i+1])==0)
        nb_nonempty_neighbs = length(nonempty_neighb_indices)
        neighbours = cbind(spatial_data2m[,(j-1)*nrow_grid+i+1],spatial_data2m[,j*nrow_grid+i],spatial_data2m[,j*nrow_grid+i+1])
        nb_reads = sum(neighbours)/nb_nonempty_neighbs
        # motu_compo contains the MOTU proportions from which the new read abundances are sampled
        # motu_compo corresponds to the mean taxonomic composition of the non-empty neighbouring samples
        motu_compo = vector(length=nrow(data2m),mode="numeric")
        for (k in nonempty_neighb_indices)
        {
          if (sum(neighbours[,k]) != 0)
            motu_compo = motu_compo + neighbours[,k]/sum(neighbours[,k])/nb_nonempty_neighbs
        }
        data2m_filled[,(j-1)*nrow_grid+i] = rmultinom(1, nb_reads, motu_compo) 
      }
      # upper right corner
      else if (i==1 && j==ncol_grid)
      {
        nonempty_neighb_indices = which(c(Missing_positions_indices[(j-2)*nrow_grid+i],Missing_positions_indices[(j-2)*nrow_grid+i+1],Missing_positions_indices[(j-1)*nrow_grid+i+1])==0)
        nb_nonempty_neighbs = length(nonempty_neighb_indices)
        neighbours = cbind(spatial_data2m[,(j-2)*nrow_grid+i],spatial_data2m[,(j-2)*nrow_grid+i+1],spatial_data2m[,(j-1)*nrow_grid+i+1])
        nb_reads = sum(neighbours)/nb_nonempty_neighbs
        motu_compo = vector(length=nrow(data2m),mode="numeric")
        for (k in nonempty_neighb_indices)
        {
          if (sum(neighbours[,k]) != 0)
            motu_compo = motu_compo + neighbours[,k]/sum(neighbours[,k])/nb_nonempty_neighbs
        }
        data2m_filled[,(j-1)*nrow_grid+i] = rmultinom(1, nb_reads, motu_compo)
      }
      # lower left corner
      else if (i==nrow_grid && j==1)
      {
        nonempty_neighb_indices = which(c(Missing_positions_indices[(j-1)*nrow_grid+i-1],Missing_positions_indices[j*nrow_grid+i-1],Missing_positions_indices[j*nrow_grid+i])==0)
        nb_nonempty_neighbs = length(nonempty_neighb_indices)
        neighbours = cbind(spatial_data2m[,(j-1)*nrow_grid+i-1],spatial_data2m[,j*nrow_grid+i-1],spatial_data2m[,j*nrow_grid+i])
        nb_reads = sum(neighbours)/nb_nonempty_neighbs
        motu_compo = vector(length=nrow(data2m),mode="numeric")
        for (k in nonempty_neighb_indices)
        {
          if (sum(neighbours[,k]) != 0)
            motu_compo = motu_compo + neighbours[,k]/sum(neighbours[,k])/nb_nonempty_neighbs
        }
        data2m_filled[,(j-1)*nrow_grid+i] = rmultinom(1, nb_reads, motu_compo)
      }
      # lower right corner
      else if (i==nrow_grid && j==ncol_grid)
      {
        nonempty_neighb_indices = which(c(Missing_positions_indices[(j-2)*nrow_grid+i],Missing_positions_indices[(j-2)*nrow_grid+i-1],Missing_positions_indices[(j-1)*nrow_grid+i-1])==0)
        nb_nonempty_neighbs = length(nonempty_neighb_indices)
        neighbours = cbind(spatial_data2m[,(j-2)*nrow_grid+i],spatial_data2m[,(j-2)*nrow_grid+i-1],spatial_data2m[,(j-1)*nrow_grid+i-1])
        nb_reads = sum(neighbours)/nb_nonempty_neighbs
        motu_compo = vector(length=nrow(data2m),mode="numeric")
        for (k in nonempty_neighb_indices)
        {
          if (sum(neighbours[,k]) != 0)
            motu_compo = motu_compo + neighbours[,k]/sum(neighbours[,k])/nb_nonempty_neighbs
        }
        data2m_filled[,(j-1)*nrow_grid+i] = rmultinom(1, nb_reads, motu_compo)
      }
      # first row
      else if (i==1 && j!=1 && j!=ncol_grid)
      {
        nonempty_neighb_indices = which(c(Missing_positions_indices[(j-2)*nrow_grid+i],Missing_positions_indices[(j-2)*nrow_grid+i+1],
                                          Missing_positions_indices[(j-1)*nrow_grid+i+1],Missing_positions_indices[j*nrow_grid+i+1],Missing_positions_indices[j*nrow_grid+i])==0)
        nb_nonempty_neighbs = length(nonempty_neighb_indices)
        neighbours = cbind(spatial_data2m[,(j-2)*nrow_grid+i],spatial_data2m[,(j-2)*nrow_grid+i+1],
                           spatial_data2m[,(j-1)*nrow_grid+i+1],spatial_data2m[,j*nrow_grid+i+1],spatial_data2m[,j*nrow_grid+i])
        nb_reads = sum(neighbours)/nb_nonempty_neighbs
        motu_compo = vector(length=nrow(data2m),mode="numeric")
        for (k in nonempty_neighb_indices)
        {
          if (sum(neighbours[,k]) != 0)
            motu_compo = motu_compo + neighbours[,k]/sum(neighbours[,k])/nb_nonempty_neighbs
        }
        data2m_filled[,(j-1)*nrow_grid+i] = rmultinom(1, nb_reads, motu_compo)
      }
      # last row
      else if (i==nrow_grid && j!=1 && j!=ncol_grid)
      {
        nonempty_neighb_indices = which(c(Missing_positions_indices[(j-2)*nrow_grid+i],Missing_positions_indices[(j-2)*nrow_grid+i-1],
                                          Missing_positions_indices[(j-1)*nrow_grid+i-1],Missing_positions_indices[j*nrow_grid+i-1],Missing_positions_indices[j*nrow_grid+i])==0)
        nb_nonempty_neighbs = length(nonempty_neighb_indices)
        neighbours = cbind(spatial_data2m[,(j-2)*nrow_grid+i],spatial_data2m[,(j-2)*nrow_grid+i-1],
                           spatial_data2m[,(j-1)*nrow_grid+i-1],spatial_data2m[,j*nrow_grid+i-1],spatial_data2m[,j*nrow_grid+i])
        nb_reads = sum(neighbours)/nb_nonempty_neighbs
        motu_compo = vector(length=nrow(data2m),mode="numeric")
        for (k in nonempty_neighb_indices)
        {
          if (sum(neighbours[,k]) != 0)
            motu_compo = motu_compo + neighbours[,k]/sum(neighbours[,k])/nb_nonempty_neighbs
        }
        data2m_filled[,(j-1)*nrow_grid+i] = rmultinom(1, nb_reads, motu_compo)
      }
      # first column
      else if (j==1 && i!=1 && i!=nrow_grid)
      {
        nonempty_neighb_indices = which(c(Missing_positions_indices[(j-1)*nrow_grid+i-1],Missing_positions_indices[j*nrow_grid+i-1],
                                          Missing_positions_indices[j*nrow_grid+i],Missing_positions_indices[j*nrow_grid+i+1],Missing_positions_indices[(j-1)*nrow_grid+i+1])==0)
        nb_nonempty_neighbs = length(nonempty_neighb_indices)
        neighbours = cbind(spatial_data2m[,(j-1)*nrow_grid+i-1],spatial_data2m[,j*nrow_grid+i-1],
                           spatial_data2m[,j*nrow_grid+i],spatial_data2m[,j*nrow_grid+i+1],spatial_data2m[,(j-1)*nrow_grid+i+1])
        nb_reads = sum(neighbours)/nb_nonempty_neighbs
        motu_compo = vector(length=nrow(data2m),mode="numeric")
        for (k in nonempty_neighb_indices)
        {
          if (sum(neighbours[,k]) != 0)
            motu_compo = motu_compo + neighbours[,k]/sum(neighbours[,k])/nb_nonempty_neighbs
        }
        data2m_filled[,(j-1)*nrow_grid+i] = rmultinom(1, nb_reads, motu_compo)
      }
      # last column
      else if (j==ncol_grid && i!=1 && i!=nrow_grid)
      {
        nonempty_neighb_indices = which(c(Missing_positions_indices[(j-1)*nrow_grid+i-1],Missing_positions_indices[(j-2)*nrow_grid+i-1],
                                          Missing_positions_indices[(j-2)*nrow_grid+i],Missing_positions_indices[(j-2)*nrow_grid+i+1],Missing_positions_indices[(j-1)*nrow_grid+i+1])==0)
        nb_nonempty_neighbs = length(nonempty_neighb_indices)
        neighbours = cbind(spatial_data2m[,(j-1)*nrow_grid+i-1],spatial_data2m[,(j-2)*nrow_grid+i-1],
                           spatial_data2m[,(j-2)*nrow_grid+i],spatial_data2m[,(j-2)*nrow_grid+i+1],spatial_data2m[,(j-1)*nrow_grid+i+1])
        nb_reads = sum(neighbours)/nb_nonempty_neighbs
        motu_compo = vector(length=nrow(data2m),mode="numeric")
        for (k in nonempty_neighb_indices)
        {
          if (sum(neighbours[,k]) != 0)
            motu_compo = motu_compo + neighbours[,k]/sum(neighbours[,k])/nb_nonempty_neighbs
        }
        data2m_filled[,(j-1)*nrow_grid+i] = rmultinom(1, nb_reads, motu_compo)
      }
      # regular entries
      else 
      {
        nonempty_neighb_indices = which(c(Missing_positions_indices[(j-1)*nrow_grid+i-1],Missing_positions_indices[j*nrow_grid+i-1],
                                          Missing_positions_indices[j*nrow_grid+i],Missing_positions_indices[j*nrow_grid+i+1],
                                          Missing_positions_indices[(j-1)*nrow_grid+i+1],Missing_positions_indices[(j-2)*nrow_grid+i+1],
                                          Missing_positions_indices[(j-2)*nrow_grid+i],Missing_positions_indices[(j-2)*nrow_grid+i-1])==0)
        nb_nonempty_neighbs = length(nonempty_neighb_indices)
        neighbours = cbind(spatial_data2m[,(j-1)*nrow_grid+i-1],spatial_data2m[,j*nrow_grid+i-1],
                           spatial_data2m[,j*nrow_grid+i],spatial_data2m[,j*nrow_grid+i+1],
                           spatial_data2m[,(j-1)*nrow_grid+i+1],spatial_data2m[,(j-2)*nrow_grid+i+1],
                           spatial_data2m[,(j-2)*nrow_grid+i],spatial_data2m[,(j-2)*nrow_grid+i-1])
        nb_reads = sum(neighbours)/nb_nonempty_neighbs
        motu_compo = vector(length=nrow(data2m),mode="numeric")
        for (k in nonempty_neighb_indices)
        {
          if (sum(neighbours[,k]) != 0)
            motu_compo = motu_compo + neighbours[,k]/sum(neighbours[,k])/nb_nonempty_neighbs
        }
        data2m_filled[,(j-1)*nrow_grid+i] = rmultinom(1, nb_reads, motu_compo)
      }
    }
  }
} 

# Ensuring that samples are listed in their original order in data2m:
if (byRow)
  data2m = data2m_filled[,as.vector(t(matrix(1:(nrow_grid*ncol_grid),ncol=ncol_grid)))]
else
  data2m = data2m_filled

return(data2m)

}
