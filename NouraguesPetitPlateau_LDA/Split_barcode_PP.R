PP_euka = 0
PP_fungi = 1

path.to.data.folder = ""

################################################################################
# Extracting plants, fungi and metazoa from 18S barcode for Petit Plateau data #
################################################################################
if (PP_euka)
{
setwd(paste0(path.to.data.folder,"/Eukaryotes_18S/"))

taxo_ref = read.table("PP_18Seuka_taxo_ref_table_assignScore.txt", sep=" ", colClasses="vector")
data2m = read.table("PP_18Seuka_sequences_counts_norep.txt", sep=" ", colClasses="numeric", header = T)
taxo_ref_protists = read.table("protists1.csv", sep=";", colClasses="vector")
full_data_euka = read.table("PP_18Seuka_sequences_sep.txt", sep=";", colClasses="vector")

i_fungi=2
i_metazoa=2
i_plants=2
i_annelids=2
i_arthropods=2
i_nematodes=2
i_protists=2
i_platyhelminthes=2
taxo_ref_euka_fungi = taxo_ref
taxo_ref_euka_metazoa = taxo_ref
taxo_ref_euka_plants = taxo_ref
taxo_ref_euka_annelids = taxo_ref
taxo_ref_euka_arthropods = taxo_ref
taxo_ref_euka_nematodes = taxo_ref
taxo_ref_euka = taxo_ref
taxo_ref_euka_protists = taxo_ref
taxo_ref_euka_platyhelminthes = taxo_ref
data2m_euka_fungi = data2m
data2m_euka_metazoa = data2m
data2m_euka_plants = data2m
data2m_euka_annelids = data2m
data2m_euka_arthropods = data2m
data2m_euka_nematodes = data2m
data2m_euka = data2m
data2m_euka_protists = data2m
data2m_euka_platyhelminthes = data2m
for (i in 2:length(taxo_ref[,1]))
{
  OTU_assigned = 0
  if (!is.na(taxo_ref[i,2]))
  {
    if (taxo_ref[i,2]=="Fungi")
    {
      taxo_ref_euka_fungi[i_fungi,]=taxo_ref[i,]
      taxo_ref_euka_fungi[i_fungi,1]=i_fungi-1
      data2m_euka_fungi[i_fungi-1,]=data2m[i-1,]
      i_fungi=i_fungi+1
      OTU_assigned = 1
    } else if (taxo_ref[i,2]=="Metazoa")
    {
      taxo_ref_euka_metazoa[i_metazoa,]=taxo_ref[i,]
      taxo_ref_euka_metazoa[i_metazoa,1]=i_metazoa-1
      data2m_euka_metazoa[i_metazoa-1,]=data2m[i-1,]
      i_metazoa=i_metazoa+1
      OTU_assigned = 1
      if (!is.na(taxo_ref[i,3]))
      {
        if (taxo_ref[i,3]=="Annelida")
        {
          taxo_ref_euka_annelids[i_annelids,]=taxo_ref[i,]
          taxo_ref_euka_annelids[i_annelids,1]=i_annelids-1
          data2m_euka_annelids[i_annelids-1,]=data2m[i-1,]
          i_annelids=i_annelids+1
        } else if (taxo_ref[i,3]=="Arthropoda")
        {
          taxo_ref_euka_arthropods[i_arthropods,]=taxo_ref[i,]
          taxo_ref_euka_arthropods[i_arthropods,1]=i_arthropods-1
          data2m_euka_arthropods[i_arthropods-1,]=data2m[i-1,]
          i_arthropods=i_arthropods+1
        } else if (taxo_ref[i,3]=="Nematoda")
        {
          taxo_ref_euka_nematodes[i_nematodes,]=taxo_ref[i,]
          taxo_ref_euka_nematodes[i_nematodes,1]=i_nematodes-1
          data2m_euka_nematodes[i_nematodes-1,]=data2m[i-1,]
          i_nematodes=i_nematodes+1
        } else if (taxo_ref[i,3]=="Platyhelminthes")
        {
          taxo_ref_euka_platyhelminthes[i_platyhelminthes,]=taxo_ref[i,]
          taxo_ref_euka_platyhelminthes[i_platyhelminthes,1]=i_platyhelminthes-1
          data2m_euka_platyhelminthes[i_platyhelminthes-1,]=data2m[i-1,]
          i_platyhelminthes=i_platyhelminthes+1
        }
      }
    } else if (taxo_ref[i,2]=="Viridiplantae")
    {
      taxo_ref_euka_plants[i_plants,]=taxo_ref[i,]
      taxo_ref_euka_plants[i_plants,1]=i_plants-1
      data2m_euka_plants[i_plants-1,]=data2m[i-1,]
      i_plants=i_plants+1
      OTU_assigned = 1
    }  
  }
  if (!OTU_assigned)
  {
    prot = 2
    while ((prot < (length(taxo_ref_protists[,1])+1)) && !OTU_assigned) 
    {
      if (taxo_ref_protists[prot,1]==full_data_euka[i,1])
      {
        taxo_ref_euka_protists[i_protists,]=taxo_ref[i,]
        taxo_ref_euka_protists[i_protists,1]=i_protists-1
        data2m_euka_protists[i_protists-1,]=data2m[i-1,]
        i_protists = i_protists+1
        OTU_assigned = 1
      }
      prot = prot+1
    }
  }
}
# taxo_ref : data frame with column names in the first line
taxo_ref_euka_fungi = taxo_ref_euka_fungi[1:(i_fungi-1),]
taxo_ref_euka_metazoa = taxo_ref_euka_metazoa[1:(i_metazoa-1),]
taxo_ref_euka_plants = taxo_ref_euka_plants[1:(i_plants-1),]
taxo_ref_euka_annelids = taxo_ref_euka_annelids[1:(i_annelids-1),]
taxo_ref_euka_arthropods = taxo_ref_euka_arthropods[1:(i_arthropods-1),]
taxo_ref_euka_nematodes = taxo_ref_euka_nematodes[1:(i_nematodes-1),]
taxo_ref_euka_protists = taxo_ref_euka_protists[1:(i_protists-1),]
taxo_ref_euka_platyhelminthes = taxo_ref_euka_platyhelminthes[1:(i_platyhelminthes-1),]
# data2m : matrix without first line for the colomn names
data2m_euka_fungi = data2m_euka_fungi[1:(i_fungi-2),]
data2m_euka_metazoa = data2m_euka_metazoa[1:(i_metazoa-2),]
data2m_euka_plants = data2m_euka_plants[1:(i_plants-2),]
data2m_euka_annelids = data2m_euka_annelids[1:(i_annelids-2),]
data2m_euka_arthropods = data2m_euka_arthropods[1:(i_arthropods-2),]
data2m_euka_nematodes = data2m_euka_nematodes[1:(i_nematodes-2),]
data2m_euka_protists = data2m_euka_protists[1:(i_protists-2),]
data2m_euka_platyhelminthes = data2m_euka_platyhelminthes[1:(i_platyhelminthes-2),]

###################
setwd(paste0(path.to.data.folder,"/Champignons_18S/"))
taxo_ref = taxo_ref_euka_fungi
data2m = data2m_euka_fungi
save(taxo_ref, file="taxo_ref.Rdata")
save(data2m, file="data2m.Rdata")

setwd(paste0(path.to.data.folder,"/Metazoaires_18S/"))
taxo_ref = taxo_ref_euka_metazoa
data2m = data2m_euka_metazoa
save(taxo_ref, file="taxo_ref.Rdata")
save(data2m, file="data2m.Rdata")

setwd(paste0(path.to.data.folder,"/Arthropodes_18S/"))
taxo_ref = taxo_ref_euka_arthropods
data2m = data2m_euka_arthropods
save(taxo_ref, file="taxo_ref.Rdata")
save(data2m, file="data2m.Rdata")

setwd(paste0(path.to.data.folder,"/Nematodes_18S/"))
taxo_ref = taxo_ref_euka_nematodes
data2m = data2m_euka_nematodes
save(taxo_ref, file="taxo_ref.Rdata")
save(data2m, file="data2m.Rdata")

setwd(paste0(path.to.data.folder,"/Annelides_18S/"))
taxo_ref = taxo_ref_euka_annelids
data2m = data2m_euka_annelids
save(taxo_ref, file="taxo_ref.Rdata")
save(data2m, file="data2m.Rdata")

setwd(paste0(path.to.data.folder,"/Plantes_18S/"))
taxo_ref = taxo_ref_euka_plants
data2m = data2m_euka_plants
save(taxo_ref, file="taxo_ref.Rdata")
save(data2m, file="data2m.Rdata")

setwd(paste0(path.to.data.folder,"/Protistes_18S/"))
taxo_ref = taxo_ref_euka_protists
data2m = data2m_euka_protists
save(taxo_ref, file="taxo_ref.Rdata")
save(data2m, file="data2m.Rdata")

setwd(paste0(path.to.data.folder,"/Platyhelminthes_18S/"))
taxo_ref = taxo_ref_euka_platyhelminthes
data2m = data2m_euka_platyhelminthes
save(taxo_ref, file="taxo_ref.Rdata")
save(data2m, file="data2m.Rdata")

} else if (PP_fungi)
{
  setwd(paste0(path.to.data.folder,"/Champignons_ITS/"))
  
  taxo_ref = read.table("PP_ITSfungi_taxo_ref_table_assignScore.txt", sep=" ", colClasses="vector")
  data2m = read.table("PP_ITSfungi_sequences_counts_norep.txt", sep=" ", colClasses="numeric", header = T)
  
  i_asco = 2
  i_basidio = 2
  i_glomero = 2
  taxo_ref_asco = taxo_ref
  taxo_ref_basidio = taxo_ref
  taxo_ref_glomero = taxo_ref
  data2m_asco = data2m
  data2m_basidio = data2m
  data2m_glomero = data2m
  for (i in 2:length(taxo_ref[,1]))
  {
    if (!is.na(taxo_ref[i,3]))
    {
      if (taxo_ref[i,3]=="Ascomycota")
      {
        taxo_ref_asco[i_asco,]=taxo_ref[i,]
        taxo_ref_asco[i_asco,1]=i_asco-1
        data2m_asco[i_asco-1,]=data2m[i-1,]
        i_asco=i_asco+1
      } else if (taxo_ref[i,3]=="Basidiomycota")
      {
        taxo_ref_basidio[i_basidio,]=taxo_ref[i,]
        taxo_ref_basidio[i_basidio,1]=i_basidio-1
        data2m_basidio[i_basidio-1,]=data2m[i-1,]
        i_basidio=i_basidio+1
      } else if (taxo_ref[i,3]=="Glomeromycota")
      {
        taxo_ref_glomero[i_glomero,]=taxo_ref[i,]
        taxo_ref_glomero[i_glomero,1]=i_glomero-1
        data2m_glomero[i_glomero-1,]=data2m[i-1,]
        i_glomero=i_glomero+1
      }  
    }
  }
  # taxo_ref : data frame with column names in the first line
  taxo_ref_asco = taxo_ref_asco[1:(i_asco-1),]
  taxo_ref_basidio = taxo_ref_basidio[1:(i_basidio-1),]
  taxo_ref_glomero = taxo_ref_glomero[1:(i_glomero-1),]
  # data2m : matrix without first line for the colomn names
  data2m_asco = data2m_asco[1:(i_asco-2),]
  data2m_basidio = data2m_basidio[1:(i_basidio-2),]
  data2m_glomero = data2m_glomero[1:(i_glomero-2),]
  
  ###################
  setwd(paste0(path.to.data.folder,"/Ascomycota_ITS/"))
  taxo_ref = taxo_ref_asco
  data2m = data2m_asco
  save(taxo_ref,file="taxo_ref.Rdata")
  save(data2m,file="data2m.Rdata")
  
  setwd(paste0(path.to.data.folder,"/Basidiomycota_ITS/"))
  taxo_ref = taxo_ref_basidio
  data2m = data2m_basidio
  save(taxo_ref,file="taxo_ref.Rdata")
  save(data2m,file="data2m.Rdata")
  
  setwd("/Users/guilhemsommeria-klein/Desktop/These/Donnees_PetitPlateau/Glomeromycetes_ITS/")
  taxo_ref = taxo_ref_glomero
  data2m = data2m_glomero
  save(taxo_ref,file="taxo_ref.Rdata")
  save(data2m,file="data2m.Rdata")
}





