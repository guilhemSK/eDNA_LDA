
bact = 1
prot = 0

occ = 1
abund = 0

nb_topics = 3

# family_level produces tables and phylum_level figures
# phylum_levels also details classes in proteobacteria
family_level = 0
phylum_level = 1

titles = c("Terra firme","Hydromorphic","Exposed rock")
short_titles = c("T.F.","H.","E.R.")

figure_folder = "/Users/guilhemsommeria-klein/Desktop/Manuscrits/ThesisPaper2_Topic_Modelling_PP/Figures/Submission2_analyses"

# Only for protists abundances:
# titles = c("Terra firme","Exposed rock","Hydromorphic")
##############

if (bact)
{
  setwd("/Users/guilhemsommeria-klein/Desktop/These/Donnees_PetitPlateau/Bacteries_16S")
  new_table = read.table("pp_bact---ssu---otus.csv", colClasses="vector", sep="\t", fill=T)
  old_table = read.table("metabar_PP_BA.uniq.trim.clust.taxo_agg.tabRF.txt", colClasses="vector", sep="\t")
  
  taxo_ref = read.table("PP_16Sbact_taxo_ref_table_assignScore.txt", colClasses="vector", sep=" ")
  
  if (occ)
  {
    #Bacteria occurrences:
    load("/Users/guilhemsommeria-klein/Desktop/These/Donnees_PetitPlateau/Bacteries_16S/Rtopicmodels_LDA_VEM_nb_topics3_nb_real100_em_tol1e-07_var_tol1e-08_best_keep_occurrence/ordered_realization_number_1/topics_sequence_composition_info/Taxonomic_composition.Rdata")
  } else if (abund)
  {
    #Bacteria abundances:
    load("/Users/guilhemsommeria-klein/Desktop/These/Donnees_PetitPlateau/Bacteries_16S/Rtopicmodels_LDA_VEM_nb_topics3_nb_real100_em_tol1e-07_var_tol1e-08_best_keep/ordered_realization_number_1/topics_sequence_composition_info/Taxonomic_composition.Rdata")
  }
  
  # test = lapply(as.character(new_table1$silva),strsplit,split="|",fixed=T)
  # bacteria_assign1 = vector(length = length(new_table1$silva), mode = "character")
  # bacteria_assign2 = vector(length = length(new_table1$silva), mode = "character")
  # for (i in 1:length(new_table1$silva))
  # {
  #   bacteria_assign1[i] = unlist(strsplit(unlist(test[[i]])[length(unlist(test[[i]]))],";",fixed=T))[2]
  #   bacteria_assign2[i] = unlist(strsplit(unlist(test[[i]])[length(unlist(test[[i]]))],";",fixed=T))[2]
  # }
  # all(bacteria_assign1=="Bacteria")
  # (sum(table(bacteria_assign2)) == length(bacteria_assign))
  
} else if (prot)
{
  # setwd("/Users/guilhemsommeria-klein/Desktop/These/Donnees_PetitPlateau/Protistes_18S")
  # load("taxo_ref_euka_protists.Rdata")
  
  setwd("/Users/guilhemsommeria-klein/Desktop/These/Donnees_PetitPlateau/Eukaryotes_18S")
  new_table = read.table("pp_euk---ssu---otus.csv", colClasses="vector", sep="\t", fill=T)
  old_table = read.table("metabar_PP_18s.uniq.trim.clust.taxo_agg.tabRF.txt", colClasses="vector", sep="\t")
  
  taxo_ref = read.table("PP_18Seuka_taxo_ref_table_assignScore.txt", sep=" ", colClasses="vector")
  taxo_ref_protists = read.table("protists1.csv", sep=";", colClasses="vector")
  full_data_euka = read.table("PP_18Seuka_sequences_sep.txt", sep=";", colClasses="vector")
  
  i_protists=2
  taxo_ref_euka_protists = taxo_ref
  for (i in 2:length(taxo_ref[,1]))
  {
    OTU_assigned = 0
    if (!(!is.na(taxo_ref[i,2]) && (taxo_ref[i,2]=="Fungi" || taxo_ref[i,2]=="Metazoa" || taxo_ref[i,2]=="Viridiplantae")))
    {
      prot = 2
      while ((prot < (length(taxo_ref_protists[,1])+1)) && !OTU_assigned) 
      {
        if (taxo_ref_protists[prot,1]==full_data_euka[i,1])
        {
          taxo_ref_euka_protists[i_protists,]=taxo_ref[i,]
          #taxo_ref_euka_protists[i_protists,1]=i_protists-1
          i_protists = i_protists+1
          OTU_assigned = 1
        }
        prot = prot+1
      }
    }
  }
  # taxo_ref : data frame with column names in the first line
  taxo_ref_euka_protists = taxo_ref_euka_protists[1:(i_protists-1),]
  taxo_ref = taxo_ref_euka_protists
  
  setwd("/Users/guilhemsommeria-klein/Desktop/These/Donnees_PetitPlateau/Protistes_18S")
  
  if (occ)
  {
    load("/Users/guilhemsommeria-klein/Desktop/These/Donnees_PetitPlateau/Protistes_18S/Rtopicmodels_LDA_VEM/Rtopicmodels_LDA_VEM_nb_topics3_nb_real100_em_tol1e-07_var_tol1e-08_best_keep_occurrence/ordered_realization_number_1/topics_sequence_composition_info/Taxonomic_composition.Rdata")
  } else if (abund)
  {
    load("/Users/guilhemsommeria-klein/Desktop/These/Donnees_PetitPlateau/Protistes_18S/Rtopicmodels_LDA_VEM/Rtopicmodels_LDA_VEM_nb_topics3_nb_real100_em_tol1e-07_var_tol1e-08_best_keep/ordered_realization_number_1/topics_sequence_composition_info/Taxonomic_composition.Rdata")
  }
}

colnames(taxo_ref) = taxo_ref[1,]
taxo_ref = taxo_ref[-1,]

load("data2m_filled.Rdata")
data2m = apply(data2m,2,as.numeric)

if (occ && any(data2m != 1 & data2m != 0))
  data2m[data2m > 0] = 1

OTUs_to_be_removed = vector(length=nrow(data2m),mode="logical")
Prop_OTU_removed = 0
Prop_reads_removed = 0
# This step slows down the computation for large datasets:
for (OTU in 1:nrow(data2m))
  # Removing all OTUs which are present in less than nb_occupied_sites_threshold, including OTUs with 0 total abundance counts
  OTUs_to_be_removed[OTU] = (length(which(data2m[OTU,]>0)) < 1)

if (any(OTUs_to_be_removed))
{
  Prop_OTU_removed = length(which(OTUs_to_be_removed))/nrow(data2m)
  Prop_reads_removed = sum(data2m[which(OTUs_to_be_removed),])/sum(data2m)
  data2m = data2m[-which(OTUs_to_be_removed),]
  taxo_ref = taxo_ref[-which(OTUs_to_be_removed),]
  rownames(data2m) = seq(1,nrow(data2m),1)
  rownames(taxo_ref) = seq(1,nrow(data2m),1)
  #taxo_ref[1:nrow(data2m),1] = seq(1,nrow(data2m),1)
}

###############

new_table = new_table[-1,]
colnames(old_table) = old_table[1,]
old_table = old_table[-1,]

new_table1 = data.frame(sequence_nb = vector(length = nrow(new_table), mode = "numeric"), silva = new_table$V9, ncbi = new_table$V8)
for (i in 1:nrow(new_table))
{
  new_table1[i,1] = as.numeric(unlist(strsplit(new_table[i,3],split="_",fixed=T))[4])
}

old_table1 = data.frame(sequence_nb = vector(length = nrow(old_table), mode = "numeric"), sequence_order = 1:nrow(old_table), previous_phylum = old_table$phylum_name_ok, previous_name = old_table$scientific_name_ok)
for (i in 1:nrow(old_table))
{
  old_table1[i,1] = as.numeric(unlist(strsplit(old_table[i,1],split="_",fixed=T))[4])
} 

###############

# output = list()
# for (k in 1:3)
# {
#   output[[k]] = list()
#   for (i in 1:5)
#   {
#     index = as.numeric(taxo_ref$Sequence_number[as.numeric(Ordered_topic_compo_taxo[[k]]$Sequence_number[i])])
#     output[[k]][[i]] = list(Ordered_topic_compo_taxo[[k]]$proba[i],new_table1[which(new_table1[,1] == old_table1[index,1]),], 
#                   old_table1[index,3:4])
#   }
# }

if (prot)
{
  taxon_insert = "prot"
} else if (bact)
  taxon_insert = "bact"

if (occ)
{
  abund_insert = "occ"
} else if (abund)
  abund_insert = "abund"

# output1 contains the taxonomic composition of each topic (k = 1:3) and of the whole dataset (k = 4), with in each case the proportion, sequence and taxonomic assignment of every OTU, the OTUs being sorted by proportion 
# For the whole dataset, proportions are computed based on occurrences or sequence reads depending on the option (occ or abund) 
output1 = list()
for (k in 1:3)
{
  cum_prop = 0
  for (i in 1:length(Ordered_topic_compo_taxo[[k]]$Sequence_number))
  {
    taxo_ref_index = as.numeric(Ordered_topic_compo_taxo[[k]]$Sequence_number[i])
    old_table_index = as.numeric(taxo_ref$Sequence_number[taxo_ref_index])
    new_table_index = which(new_table1[,1] == old_table1[old_table_index,1])
    cum_prop = Ordered_topic_compo_taxo[[k]]$proba[i] + cum_prop
    if (length(new_table_index) != 0)
    {
      if (i == 1)
      {
        output1[[k]] = data.frame(proportion = Ordered_topic_compo_taxo[[k]]$proba[i], cumulative_proportion = cum_prop, new_table1[new_table_index,], 
                                  old_table1[old_table_index,3:4])
      } else 
        output1[[k]] = rbind(output1[[k]], data.frame(proportion = Ordered_topic_compo_taxo[[k]]$proba[i], cumulative_proportion = cum_prop, new_table1[new_table_index,], 
                                                      old_table1[old_table_index,3:4]))
    } else
    {
      if (i == 1)
      {
        output1[[k]] = data.frame(proportion = Ordered_topic_compo_taxo[[k]]$proba[i], cumulative_proportion = cum_prop, sequence_nb = old_table1[old_table_index,1], silva = NA, ncbi = NA, 
                                  old_table1[old_table_index,3:4])
      } else 
        output1[[k]] = rbind(output1[[k]],data.frame(proportion = Ordered_topic_compo_taxo[[k]]$proba[i], cumulative_proportion = cum_prop, sequence_nb = old_table1[old_table_index,1], silva = NA, ncbi = NA, 
                                                     old_table1[old_table_index,3:4]))
    }
  }
  # Replacing ";" with " - " in the output:
  output1[[k]]$silva = gsub(";"," - ",output1[[k]]$silva)
  output1[[k]]$ncbi = gsub(";"," - ",output1[[k]]$ncbi)
  
  if (titles[k] == "Terra firme")
    topic_insert = "terraFirme"
  else if (titles[k] == "Hydromorphic")
    topic_insert = "hydromorphic"
  else if (titles[k] == "Exposed rock")
    topic_insert = "exposedRock"
  
  #write.table(output1[[k]], file = paste0("Topic_compo_",taxon_insert,"_",abund_insert,"_",topic_insert,".csv"), quote = F, sep = ";", col.names = T, row.names = F)
}

cum_read_count = 0
OTU_read_counts = rowSums(data2m)
tot_read_count = sum(OTU_read_counts)
ordered_taxo_ref = taxo_ref[sort.int(OTU_read_counts,decreasing=T,index.return = T)$ix,]
ordered_read_counts = sort(OTU_read_counts,decreasing =T)
for (i in 1:nrow(ordered_taxo_ref))
{
  old_table_index = as.numeric(ordered_taxo_ref$Sequence_number[i])
  new_table_index = which(new_table1[,1] == old_table1[old_table_index,1])
  cum_read_count = ordered_read_counts[i] + cum_read_count
  if (length(new_table_index) != 0)
  {
    if (i == 1)
    {
      output1[[4]] = data.frame(proportion = ordered_read_counts[i]/tot_read_count, cumulative_proportion = cum_read_count/tot_read_count, new_table1[new_table_index,], 
                                old_table1[old_table_index,3:4])
    } else 
      output1[[4]] = rbind(output1[[4]], data.frame(proportion = ordered_read_counts[i]/tot_read_count, cumulative_proportion = cum_read_count/tot_read_count, new_table1[new_table_index,], 
                                                    old_table1[old_table_index,3:4]))
  } else
  {
    if (i == 1)
    {
      output1[[4]] = data.frame(proportion = ordered_read_counts[i]/tot_read_count, cumulative_proportion = cum_read_count/tot_read_count, sequence_nb = old_table1[old_table_index,1], silva = NA, ncbi = NA, 
                                old_table1[old_table_index,3:4])
    } else 
      output1[[4]] = rbind(output1[[4]],data.frame(proportion = ordered_read_counts[i]/tot_read_count, cumulative_proportion = cum_read_count/tot_read_count, sequence_nb = old_table1[old_table_index,1], silva = NA, ncbi = NA, 
                                                   old_table1[old_table_index,3:4]))
  }
  if (nrow(output1[[4]]) != i)
    stop("error")
  
  # Replacing ";" with " - " in the output:
  output1[[4]]$silva = gsub(";"," - ",output1[[4]]$silva)
  output1[[4]]$ncbi = gsub(";"," - ",output1[[4]]$ncbi)
}

# pdf(paste0("Topic_RAD_",taxon_insert,"_",abund_insert,".pdf"),width=7*3/2,height=7/2)
# par(mfrow = c(1,3))
# for (k in 1:3)
#   plot(1:nrow(data2m),c(log10(output1[[k]]$proportion), rep(NA,nrow(data2m) - length(output1[[k]]$proportion))),ylab="log10(proportion)",xlab="Rank",main = titles[k])
# dev.off()

#save(output1,file = "output1_bact_occ.Rdata")

# write.table(data.frame(slope = signi_test_distance[[i_barcode]][,1], intercept = signi_test_distance[[i_barcode]][,2] + signi_test_distance[[i_barcode]][,1]*min(log10(distances1[[i_barcode]])),
#                        mean = signi_test_distance[[i_barcode]][,3], mantel.r = signi_test_distance[[i_barcode]][,4], mantel.p.value = signi_test_distance[[i_barcode]][,5]),
#             file = "Simple_Mantel_distance.csv",
#             quote = F, sep = ";",
#             col.names = T, row.names = paste0(barcode_insert_vect[i_barcode],".",distances_vect))

# silva_phylum contains the phylum assignment of each OTU in each topic and in the whole dataset
# Occurrence topic comparaison:
########################################
silva_phylum = list()
silva_phylum_endemicProp = list()
silva_phylum_ubiquistProp = list()
silva_phylum_semiEndemic1Prop = list()
silva_phylum_semiEndemic2Prop = list()
silva_phylum_prop = list()
silva_phylum_prop_endemicProp = list()
silva_phylum_prop_ubiquistProp = list()
silva_phylum_prop_semiEndemic1Prop = list()
silva_phylum_prop_semiEndemic2Prop = list()
silva_phylum_mat = matrix(nrow = max(unlist(lapply(output1,nrow))), ncol = 3, data = "NA")
silva_phylum_endemicProp_allPhyla = vector(length = 3, mode = "numeric")
silva_phylum_ubiquistProp_allPhyla = vector(length = 3, mode = "numeric")
silva_phylum_prop_endemicProp_allPhyla = vector(length = 3, mode = "numeric")
silva_phylum_prop_ubiquistProp_allPhyla = vector(length = 3, mode = "numeric")
for (k in 1:4)
{
  # silva_assign contains the full silva assignment of each OTU for the current topic
  silva_assign = lapply(as.character(output1[[k]]$silva),strsplit,split="|",fixed=T)
  # silva_phylum_k and silva_phylum_mat register the phylum assignment of each OTU for the current topic
  silva_phylum_k = vector(length = nrow(output1[[k]]), mode = "character")
  if (phylum_level)
    silva_subphylum_k = vector(length = nrow(output1[[k]]), mode = "character")
  if (k < 4)
  {
    endemic_k = vector(length = nrow(output1[[k]]), mode = "logical")
    ubiquist_k = vector(length = nrow(output1[[k]]), mode = "logical")
    semiEndemic1_k = vector(length = nrow(output1[[k]]), mode = "logical")
    semiEndemic2_k = vector(length = nrow(output1[[k]]), mode = "logical")
  } else if (k == 4)
  {
    TF = vector(length = nrow(output1[[k]]), mode = "logical")
    H = vector(length = nrow(output1[[k]]), mode = "logical")
    ER = vector(length = nrow(output1[[k]]), mode = "logical")
    TFH = vector(length = nrow(output1[[k]]), mode = "logical")
    HER = vector(length = nrow(output1[[k]]), mode = "logical")
    TFER = vector(length = nrow(output1[[k]]), mode = "logical")
    TFHER = vector(length = nrow(output1[[k]]), mode = "logical")
  }
  for (i in 1:nrow(output1[[k]]))
  {
    if (family_level)
    {
      silva_phylum_k[i] = paste(unlist(strsplit(unlist(silva_assign[[i]])[length(unlist(silva_assign[[i]]))]," - ",fixed=T))[2:5],collapse = " ; ")
    } else if (phylum_level)
    {
      silva_phylum_k[i] = unlist(strsplit(unlist(silva_assign[[i]])[length(unlist(silva_assign[[i]]))]," - ",fixed=T))[2]
      silva_subphylum_k[i] = unlist(strsplit(unlist(silva_assign[[i]])[length(unlist(silva_assign[[i]]))]," - ",fixed=T))[3]
    }
    
    # if (k < 4)
    #   silva_phylum_mat[i,k] = silva_phylum_k[i]  
    if (k < 4)
    {
      # OTU i in assemblage k is endemic if it is not found in the other assemblages:
      endemic_k[i] = !any(output1[[(1:3)[-k][1]]]$sequence_nb == output1[[k]]$sequence_nb[i]) && !any(output1[[(1:3)[-k][2]]]$sequence_nb == output1[[k]]$sequence_nb[i])
      # OTU i in assemblage k is ubiquist if it found in all the other assemblages:
      ubiquist_k[i] = any(output1[[(1:3)[-k][1]]]$sequence_nb == output1[[k]]$sequence_nb[i]) && any(output1[[(1:3)[-k][2]]]$sequence_nb == output1[[k]]$sequence_nb[i])
      # OTU i in assemblage k is semi_endemic1 if it is found in assemblage (1:3)[-k][1] but not in assemblage (1:3)[-k][2], and the conversely for semi_endemic2
      semiEndemic1_k[i] = any(output1[[(1:3)[-k][1]]]$sequence_nb == output1[[k]]$sequence_nb[i]) && !any(output1[[(1:3)[-k][2]]]$sequence_nb == output1[[k]]$sequence_nb[i])
      semiEndemic2_k[i] = !any(output1[[(1:3)[-k][1]]]$sequence_nb == output1[[k]]$sequence_nb[i]) && any(output1[[(1:3)[-k][2]]]$sequence_nb == output1[[k]]$sequence_nb[i])
    } else if (k == 4)
    {
      TF_i = any(output1[[1]]$sequence_nb == output1[[k]]$sequence_nb[i])
      H_i = any(output1[[2]]$sequence_nb == output1[[k]]$sequence_nb[i])
      ER_i = any(output1[[3]]$sequence_nb == output1[[k]]$sequence_nb[i])
      TFHER[i] = TF_i && H_i && ER_i
      TF[i] = TF_i && !H_i && !ER_i
      H[i] = !TF_i && H_i && !ER_i
      ER[i] = !TF_i && !H_i && ER_i
      TFH[i] = TF_i && H_i && !ER_i
      TFER[i] = TF_i && !H_i && ER_i
      HER[i] = !TF_i && H_i && ER_i
    }
  }
  # phylum_levels also details classes in proteobacteria
  if (phylum_level)
    silva_phylum_k[silva_phylum_k == "Proteobacteria" & !(is.na(silva_phylum_k))] = silva_subphylum_k[silva_phylum_k == "Proteobacteria" & !(is.na(silva_subphylum_k))]
  silva_phylum[[k]] = table(silva_phylum_k)
  silva_phylum_prop[[k]] = vector(length = length(silva_phylum[[k]]), mode = "numeric")
  if (k < 4)
  {
    silva_phylum_endemicProp[[k]] = vector(length = length(silva_phylum[[k]]), mode = "numeric")
    silva_phylum_ubiquistProp[[k]] = vector(length = length(silva_phylum[[k]]), mode = "numeric")
    silva_phylum_semiEndemic1Prop[[k]] = vector(length = length(silva_phylum[[k]]), mode = "numeric")
    silva_phylum_semiEndemic2Prop[[k]] = vector(length = length(silva_phylum[[k]]), mode = "numeric")
    silva_phylum_prop_endemicProp[[k]] = vector(length = length(silva_phylum[[k]]), mode = "numeric")
    silva_phylum_prop_ubiquistProp[[k]] = vector(length = length(silva_phylum[[k]]), mode = "numeric")
    silva_phylum_prop_semiEndemic1Prop[[k]] = vector(length = length(silva_phylum[[k]]), mode = "numeric")
    silva_phylum_prop_semiEndemic2Prop[[k]] = vector(length = length(silva_phylum[[k]]), mode = "numeric")
  } else if (k == 4)
  {
    silva_phylum_TFendemicProp = vector(length = length(silva_phylum[[k]]), mode = "numeric")
    silva_phylum_HendemicProp = vector(length = length(silva_phylum[[k]]), mode = "numeric")
    silva_phylum_ERendemicProp = vector(length = length(silva_phylum[[k]]), mode = "numeric")
    silva_phylum_TFHERProp = vector(length = length(silva_phylum[[k]]), mode = "numeric")
    silva_phylum_prop_TFendemicProp = vector(length = length(silva_phylum[[k]]), mode = "numeric")
    silva_phylum_prop_HendemicProp = vector(length = length(silva_phylum[[k]]), mode = "numeric")
    silva_phylum_prop_ERendemicProp = vector(length = length(silva_phylum[[k]]), mode = "numeric")
    silva_phylum_prop_TFHERProp = vector(length = length(silva_phylum[[k]]), mode = "numeric")
    silva_phylum_prop_TFHendemicProp = vector(length = length(silva_phylum[[k]]), mode = "numeric")
    silva_phylum_prop_HERendemicProp = vector(length = length(silva_phylum[[k]]), mode = "numeric")
    silva_phylum_prop_TFERendemicProp = vector(length = length(silva_phylum[[k]]), mode = "numeric")
  }
  for (p in 1:length(silva_phylum[[k]]))
  {
    silva_phylum_prop[[k]][p] = sum(output1[[k]]$proportion[which(silva_phylum_k == names(silva_phylum[[k]])[p])]) 
    # silva_phylum_endemicProp[[k]][p] + silva_phylum_ubiquistProp[[k]][p] + silva_phylum_semiEndemic1Prop[[k]][p] + silva_phylum_semiEndemic2Prop[[k]][p] = 1 for each k and p
    # silva_phylum_prop_endemicProp[[k]][p] + silva_phylum_prop_ubiquistProp[[k]][p] + silva_phylum_prop_semiEndemic1Prop[[k]][p] + silva_phylum_prop_semiEndemic2Prop[[k]][p] = 1 for each k and p
    if (k < 4)
    {
      silva_phylum_endemicProp[[k]][p] = length(which(endemic_k[silva_phylum_k == names(silva_phylum[[k]])[p]]))/silva_phylum[[k]][p]
      silva_phylum_ubiquistProp[[k]][p] = length(which(ubiquist_k[silva_phylum_k == names(silva_phylum[[k]])[p]]))/silva_phylum[[k]][p]
      silva_phylum_semiEndemic1Prop[[k]][p] = length(which(semiEndemic1_k[silva_phylum_k == names(silva_phylum[[k]])[p]]))/silva_phylum[[k]][p]
      silva_phylum_semiEndemic2Prop[[k]][p] = length(which(semiEndemic2_k[silva_phylum_k == names(silva_phylum[[k]])[p]]))/silva_phylum[[k]][p]
      # which() is required here because silva_phylum_k contains some NA
      silva_phylum_prop_endemicProp[[k]][p] = sum(output1[[k]]$proportion[which(silva_phylum_k == names(silva_phylum[[k]])[p] & endemic_k)])/silva_phylum_prop[[k]][p]  
      silva_phylum_prop_ubiquistProp[[k]][p] = sum(output1[[k]]$proportion[which(silva_phylum_k == names(silva_phylum[[k]])[p] & ubiquist_k)])/silva_phylum_prop[[k]][p]
      silva_phylum_prop_semiEndemic1Prop[[k]][p] = sum(output1[[k]]$proportion[which(silva_phylum_k == names(silva_phylum[[k]])[p] & semiEndemic1_k)])/silva_phylum_prop[[k]][p]  
      silva_phylum_prop_semiEndemic2Prop[[k]][p] = sum(output1[[k]]$proportion[which(silva_phylum_k == names(silva_phylum[[k]])[p] & semiEndemic2_k)])/silva_phylum_prop[[k]][p]
      # if (is.na(silva_phylum_prop_endemicProp[[k]][p]) || is.na(silva_phylum_prop_ubiquistProp[[k]][p]) || is.na(silva_phylum_prop_semiEndemic1Prop[[k]][p]) || is.na(silva_phylum_prop_semiEndemic1Prop[[k]][p]))
      #   stop("NA value")
    } else if (k == 4)
    {
      silva_phylum_TFendemicProp[p] = length(which(TF[silva_phylum_k == names(silva_phylum[[k]])[p]]))/silva_phylum[[k]][p]
      silva_phylum_HendemicProp[p] = length(which(H[silva_phylum_k == names(silva_phylum[[k]])[p]]))/silva_phylum[[k]][p]
      silva_phylum_ERendemicProp[p] = length(which(ER[silva_phylum_k == names(silva_phylum[[k]])[p]]))/silva_phylum[[k]][p]
      silva_phylum_TFHERProp[p] = length(which(TFHER[silva_phylum_k == names(silva_phylum[[k]])[p]]))/silva_phylum[[k]][p]
      silva_phylum_prop_TFendemicProp[p] = sum(output1[[k]]$proportion[which(silva_phylum_k == names(silva_phylum[[k]])[p] & TF)])/silva_phylum_prop[[k]][p] 
      silva_phylum_prop_HendemicProp[p] = sum(output1[[k]]$proportion[which(silva_phylum_k == names(silva_phylum[[k]])[p] & H)])/silva_phylum_prop[[k]][p] 
      silva_phylum_prop_ERendemicProp[p] = sum(output1[[k]]$proportion[which(silva_phylum_k == names(silva_phylum[[k]])[p] & ER)])/silva_phylum_prop[[k]][p]
      silva_phylum_prop_TFHendemicProp[p] = sum(output1[[k]]$proportion[which(silva_phylum_k == names(silva_phylum[[k]])[p] & TFH)])/silva_phylum_prop[[k]][p] 
      silva_phylum_prop_HERendemicProp[p] = sum(output1[[k]]$proportion[which(silva_phylum_k == names(silva_phylum[[k]])[p] & HER)])/silva_phylum_prop[[k]][p] 
      silva_phylum_prop_TFERendemicProp[p] = sum(output1[[k]]$proportion[which(silva_phylum_k == names(silva_phylum[[k]])[p] & TFER)])/silva_phylum_prop[[k]][p] 
      silva_phylum_prop_TFHERProp[p] = sum(output1[[k]]$proportion[which(silva_phylum_k == names(silva_phylum[[k]])[p] & TFHER)])/silva_phylum_prop[[k]][p] 
    }
  }
  # Not used in the plots:
  if (k < 4)
  {
    silva_phylum_endemicProp_allPhyla[k] = length(which(endemic_k))/sum(silva_phylum[[k]])
    silva_phylum_ubiquistProp_allPhyla[k] = length(which(ubiquist_k))/sum(silva_phylum[[k]])
    silva_phylum_prop_endemicProp_allPhyla[k] = sum(output1[[k]]$proportion[endemic_k])/sum(output1[[k]]$proportion)
    silva_phylum_prop_ubiquistProp_allPhyla[k] = sum(output1[[k]]$proportion[ubiquist_k])/sum(output1[[k]]$proportion)
  }
}
#silva_phylum_table = table(silva_phylum_mat)[-which(names(table(silva_phylum_mat)) == "NA")]
#silva_phylum_table = table(silva_phylum_mat)[which(table(silva_phylum_mat)>threshold & names(table(silva_phylum_mat)) != "NA")]
#silva_phylum_table = table(silva_phylum_mat)[which(names(table(silva_phylum_mat)) != "NA")]
#silva_phylum_table = silva_phylum[[4]][which(names(silva_phylum[[4]]) != "NA")]

silva_phylum_mat1 = matrix(nrow = length(silva_phylum[[4]]), ncol = 4, data = 0, dimnames = list(names(silva_phylum[[4]])))
silva_phylum_prop_mat1 = matrix(nrow = length(silva_phylum[[4]]), ncol = 4, data = 0, dimnames = list(names(silva_phylum[[4]])))
silva_phylum_endemicProp_mat1 = matrix(nrow = length(silva_phylum[[4]]) + 1, ncol = 3, data = 0, dimnames = list(c(names(silva_phylum[[4]]),"All phyla")))
silva_phylum_ubiquistProp_mat1 = matrix(nrow = length(silva_phylum[[4]]) + 1, ncol = 3, data = 0, dimnames = list(c(names(silva_phylum[[4]]),"All phyla")))
silva_phylum_prop_endemicProp_mat1 = matrix(nrow = length(silva_phylum[[4]]) + 1, ncol = 3, data = 0, dimnames = list(c(names(silva_phylum[[4]]),"All phyla")))
silva_phylum_prop_ubiquistProp_mat1 = matrix(nrow = length(silva_phylum[[4]]) + 1, ncol = 3, data = 0, dimnames = list(c(names(silva_phylum[[4]]),"All phyla")))
silva_phylum_semiEndemic1Prop_mat1 = matrix(nrow = length(silva_phylum[[4]]) + 1, ncol = 3, data = 0, dimnames = list(c(names(silva_phylum[[4]]),"All phyla")))
silva_phylum_semiEndemic2Prop_mat1 = matrix(nrow = length(silva_phylum[[4]]) + 1, ncol = 3, data = 0, dimnames = list(c(names(silva_phylum[[4]]),"All phyla")))
silva_phylum_prop_semiEndemic1Prop_mat1 = matrix(nrow = length(silva_phylum[[4]]) + 1, ncol = 3, data = 0, dimnames = list(c(names(silva_phylum[[4]]),"All phyla")))
silva_phylum_prop_semiEndemic2Prop_mat1 = matrix(nrow = length(silva_phylum[[4]]) + 1, ncol = 3, data = 0, dimnames = list(c(names(silva_phylum[[4]]),"All phyla")))
for (p in 1:length(silva_phylum[[4]]))
{
  for (k in 1:4)
  {
    index = which(names(silva_phylum[[k]]) == names(silva_phylum[[4]])[p])
    if (length(index) != 0)
    {
      silva_phylum_mat1[p,k] =  silva_phylum[[k]][index]
      silva_phylum_prop_mat1[p,k] =  silva_phylum_prop[[k]][index]
      if (k < 4)
      {
        silva_phylum_endemicProp_mat1[p,k] = silva_phylum_endemicProp[[k]][index]
        silva_phylum_ubiquistProp_mat1[p,k] = silva_phylum_ubiquistProp[[k]][index]
        silva_phylum_prop_endemicProp_mat1[p,k] = silva_phylum_prop_endemicProp[[k]][index]
        silva_phylum_prop_ubiquistProp_mat1[p,k] = silva_phylum_prop_ubiquistProp[[k]][index]
        silva_phylum_semiEndemic1Prop_mat1[p,k] = silva_phylum_semiEndemic1Prop[[k]][index]
        silva_phylum_semiEndemic2Prop_mat1[p,k] = silva_phylum_semiEndemic2Prop[[k]][index]
        silva_phylum_prop_semiEndemic1Prop_mat1[p,k] = silva_phylum_prop_semiEndemic1Prop[[k]][index]
        silva_phylum_prop_semiEndemic2Prop_mat1[p,k] = silva_phylum_prop_semiEndemic2Prop[[k]][index]
      } 
    }
  }
}
# Not used in the plots:
silva_phylum_endemicProp_mat1[p+1,] = silva_phylum_endemicProp_allPhyla
silva_phylum_ubiquistProp_mat1[p+1,] = silva_phylum_ubiquistProp_allPhyla
silva_phylum_prop_endemicProp_mat1[p+1,] = silva_phylum_prop_endemicProp_allPhyla
silva_phylum_prop_ubiquistProp_mat1[p+1,] = silva_phylum_prop_ubiquistProp_allPhyla

# Computing the overall proportion of ubiquists and endemics in the dataset:
# ubiquist = vector(length = nrow(output1[[4]]), mode = "logical")
# endemic = vector(length = nrow(output1[[4]]), mode = "logical")
# for (i in 1:nrow(output1[[4]]))
# {
#   ubiquist[i] = any(output1[[1]]$sequence_nb == output1[[4]]$sequence_nb[i]) && any(output1[[2]]$sequence_nb == output1[[4]]$sequence_nb[i]) && any(output1[[3]]$sequence_nb == output1[[4]]$sequence_nb[i])
#   endemic[i] = length(which(c(output1[[4]]$sequence_nb[i] %in% output1[[1]]$sequence_nb,output1[[4]]$sequence_nb[i] %in% output1[[2]]$sequence_nb,output1[[4]]$sequence_nb[i] %in% output1[[3]]$sequence_nb))) == 1
# }
  
# silva_phylum_mat1_filtered = silva_phylum_mat1[which(rowSums(silva_phylum_mat1) > richness_threshold),]
# silva_phylum_prop_mat1_filtered = silva_phylum_prop_mat1[which(rowSums(silva_phylum_prop_mat1) > prop_threshold),]

#sort(rowSums(silva_phylum_prop_mat1),decreasing =T)

labels = c(titles,"Whole dataset")

# Computing for each clade (down to family level) in each assemblage the proportion, number of OTUs, relative difference to mean proportion and relative difference to mean proportion of OTU number
if (family_level && occ && bact)
{
  rel_diff_to_mean_prop = list()
  rel_diff_to_mean_prop_richness = list()
  rel_diff_to_mean_richness = list()
  for (k in 1:3)
  {
    names(silva_phylum_prop[[k]]) = names(silva_phylum[[k]])
    rel_diff_to_mean_prop[[k]] = (silva_phylum_prop[[k]]-rowMeans(silva_phylum_prop_mat1[names(silva_phylum[[k]]),1:3]))/rowMeans(silva_phylum_prop_mat1[names(silva_phylum[[k]]),1:3])*100
    rel_diff_to_mean_prop_richness[[k]] = (silva_phylum[[k]]/sum(silva_phylum[[k]])-rowMeans(sweep(silva_phylum_mat1,2,colSums(silva_phylum_mat1),"/")[names(silva_phylum[[k]]),1:3]))/rowMeans(sweep(silva_phylum_mat1,2,colSums(silva_phylum_mat1),"/")[names(silva_phylum[[k]]),1:3])*100
    rel_diff_to_mean_richness[[k]] = (silva_phylum[[k]]-rowMeans(silva_phylum_mat1[names(silva_phylum[[k]]),1:3]))/rowMeans(silva_phylum_mat1[names(silva_phylum[[k]]),1:3])*100
    assemblage_compo = paste(";",silva_phylum_prop[[k]]*100,";",silva_phylum[[k]],
                             ";",rel_diff_to_mean_prop[[k]],
                             ";",rel_diff_to_mean_prop_richness[[k]],
                             ";",rel_diff_to_mean_richness[[k]])
    names(assemblage_compo) = names(silva_phylum[[k]])
    write.table(paste0("\n",labels[k]),"Assemblage_compo.csv",quote = F, col.names = F, append = ifelse(k == 1, F, T))
    write.table(assemblage_compo,"Assemblage_compo.csv",quote = F, col.names = F, append = T)
    
    # Selecting families above 50 OTUs in diversity with an above-50% increase in proportion of occurrences, and families above 100 in diversity with an above-15% increase in proportion of occurrences
    # selected_families = rel_diff_to_mean_prop[[k]]>15 & silva_phylum[[k]]>100 | rel_diff_to_mean_prop[[k]]>30 & silva_phylum[[k]]>50
    selected_families = rel_diff_to_mean_prop[[k]] > 15 & silva_phylum_prop[[k]]*100 > 1 | rel_diff_to_mean_prop[[k]] > 30 & silva_phylum_prop[[k]]*100 > 0.5
    # selected_families = rel_diff_to_mean_prop_richness[[k]]>15 & silva_phylum[[k]]>100 | rel_diff_to_mean_prop_richness[[k]]>30 & silva_phylum[[k]]>50
    assemblage_compo_selected = paste(";",silva_phylum_prop[[k]][selected_families]*100,
                                      ";",silva_phylum[[k]][selected_families],
                                      ";",rel_diff_to_mean_prop[[k]][selected_families],
                                      ";",rel_diff_to_mean_prop_richness[[k]][selected_families],
                                      ";",rel_diff_to_mean_richness[[k]][selected_families])
    names(assemblage_compo_selected) = names(silva_phylum[[k]][selected_families])
    write.table(paste0("\n",labels[k]),"Assemblage_compo_selected.csv",quote = F, col.names = F, append = ifelse(k == 1, F, T))
    write.table(assemblage_compo_selected,"Assemblage_compo_selected.csv",quote = F, col.names = F, append = T)
  }
} else if (phylum_level && bact)
{
  
  # Barplot taxo composition 
  ##########################
  
  if (abund)
  {
    richness_threshold = 50
  } else if (occ)
    richness_threshold = 100
  read_threshold = 1
  selected_groups_indices = which(silva_phylum_mat1[,4] > richness_threshold & silva_phylum_prop_mat1[,4]*tot_read_count > read_threshold)
  silva_phylum_mat1_filtered = silva_phylum_mat1[selected_groups_indices,]
  
  pdf(paste0("Assemblage_comparison_phylum_proteo_",taxon_insert,"_",abund_insert,"_",richness_threshold,"OTU",read_threshold,"readsThres_propRichness.pdf"))
  par(mar=c(7.1,4.1,4.1,2.1))
  for (k in 1:4)
  {
    x = barplot(c(silva_phylum_mat1_filtered[,k],sum(silva_phylum_mat1[,k])-sum(silva_phylum_mat1_filtered[,k]))/sum(silva_phylum_mat1[,k]), space = 1, xaxt="n", 
                ylim = c(0,max(sweep(silva_phylum_mat1_filtered,2,colSums(silva_phylum_mat1),"/"))))
    labs = c(rownames(silva_phylum_mat1_filtered),"Others")
    #text(cex=0.8, x=x+1, y=-70, labels = labs, xpd=TRUE, srt=45, pos=2)
    text(cex=1.2, x=x+0.5, y=-0.01, labels = labs, xpd=TRUE, srt=45, pos=2)
    title(ylab = "Proportion of OTU richness", main = labels[k])
  }
  dev.off()
  
  richness_threshold = 1
  if (abund)
  {
    read_threshold = 100000
  } else if (occ)
    read_threshold = 5000
  selected_groups_indices = which(silva_phylum_mat1[,4] > richness_threshold & silva_phylum_prop_mat1[,4]*tot_read_count > read_threshold)
  silva_phylum_prop_mat1_filtered = silva_phylum_prop_mat1[selected_groups_indices,]
  
  pdf(paste0("Assemblage_comparison_phylum_proteo_",taxon_insert,"_",abund_insert,"_",richness_threshold,"OTU",read_threshold,"readsThres_prop.pdf"))
  par(mar=c(7.1,4.1,4.1,2.1))
  for (k in 1:4)
  {
    x = barplot(c(silva_phylum_prop_mat1_filtered[,k],1-sum(silva_phylum_prop_mat1_filtered[,k])), space = 1, xaxt="n",
                ylim = c(0,max(silva_phylum_prop_mat1_filtered)))
    labs = c(rownames(silva_phylum_prop_mat1_filtered),"Others")
    #text(cex=0.8, x=x+1, y=-70, labels = labs, xpd=TRUE, srt=45, pos=2)
    text(cex=1.2, x=x+0.5, y=-0.01, labels = labs, xpd=TRUE, srt=45, pos=2)
    if (abund)
      title(ylab = "Proportion of read count", main = labels[k])
    else if (occ)
      title(ylab = "Proportion of occurrences", main = labels[k])
  }
  dev.off()
  
  # Barplot taxo composition oneplot
  ##################################
  
  library(colorspace)
  
  if (abund)
  {
    richness_threshold = 50
  } else if (occ)
    richness_threshold = 100
  read_threshold = 1
  selected_groups_indices = which(silva_phylum_mat1[,4] > richness_threshold & silva_phylum_prop_mat1[,4]*tot_read_count > read_threshold)
  silva_phylum_mat1_filtered = silva_phylum_mat1[selected_groups_indices,]
  
  pdf(paste0("Assemblage_comparison_phylum_",taxon_insert,"_",abund_insert,"_",richness_threshold,"OTU",read_threshold,"readsThres_propRichness_oneplot.pdf"))
  par(mar=c(7.1,4.1,4.1,10.1))
  #barplot(VarPart_data.frame_reordered,col=terrain.colors(3),ann=F,names.arg=colnames(VarPart_data.frame_reordered),las=3,legend.text = T, args.legend = list(bty = "n"))
  x = barplot(sweep(rbind(silva_phylum_mat1_filtered, Others = colSums(silva_phylum_mat1)-colSums(silva_phylum_mat1_filtered)),2,colSums(silva_phylum_mat1),"/"),
              col=c(rainbow_hcl(nrow(silva_phylum_mat1))[sample(1:nrow(silva_phylum_mat1),nrow(silva_phylum_mat1))],"blue"),
              legend.text = T, xaxt="n", space = 0.5, 
              args.legend = list(bty = "n", x= 8 + 1, y=1))
  #colours = terrain.colors(nrow(silva_phylum_mat1))[sample(1:nrow(silva_phylum_mat1),nrow(silva_phylum_mat1))]
  #x = barplot(sweep(silva_phylum_mat1,2,colSums(silva_phylum_mat1),"/"),col=colours,
  #xaxt="n", space = 0.5)
  #legend("topright",rownames(silva_phylum_mat1),col=colours)
  title(ylab="Proportion of OTU richness",cex.lab=1.3)
  labs = labels
  text(cex=1, x=x+0.3, y=-0.02, labels = labs, xpd=TRUE, srt=45, pos=2)
  dev.off()
  
  richness_threshold = 1
  if (abund)
  {
    read_threshold = 100000
  } else if (occ)
    read_threshold = 5000
  selected_groups_indices = which(silva_phylum_mat1[,4] > richness_threshold & silva_phylum_prop_mat1[,4]*tot_read_count > read_threshold)
  silva_phylum_prop_mat1_filtered = silva_phylum_prop_mat1[selected_groups_indices,]
  
  pdf(paste0("Assemblage_comparison_phylum_",taxon_insert,"_",abund_insert,"_",richness_threshold,"OTU",read_threshold,"readsThres_prop_oneplot.pdf"))
  par(mar=c(7.1,4.1,4.1,10.1))
  #barplot(VarPart_data.frame_reordered,col=terrain.colors(3),ann=F,names.arg=colnames(VarPart_data.frame_reordered),las=3,legend.text = T, args.legend = list(bty = "n"))
  x = barplot(rbind(silva_phylum_prop_mat1_filtered,Others = colSums(silva_phylum_prop_mat1)-colSums(silva_phylum_prop_mat1_filtered)),
              col=c(rainbow_hcl(nrow(silva_phylum_prop_mat1_filtered))[sample(1:nrow(silva_phylum_prop_mat1_filtered),nrow(silva_phylum_prop_mat1_filtered))],"blue"),
              legend.text = T, xaxt="n", space = 0.5, 
              args.legend = list(bty = "n", x= 8 + 1, y=1))
  #colours = terrain.colors(nrow(silva_phylum_mat1))[sample(1:nrow(silva_phylum_mat1),nrow(silva_phylum_mat1))]
  #x = barplot(sweep(silva_phylum_mat1,2,colSums(silva_phylum_mat1),"/"),col=colours,
  #xaxt="n", space = 0.5)
  #legend("topright",rownames(silva_phylum_mat1),col=colours)
  if (abund)
  {
    title(ylab="Proportion of read count",cex.lab=1.3)
  } else if (occ)
    title(ylab="Proportion of occurrences",cex.lab=1.3)
  labs = labels
  text(cex=1, x=x+0.3, y=-0.02, labels = labs, xpd=TRUE, srt=45, pos=2)
  dev.off()
  
  # Relative difference to mean proportion
  ########################################
  
  # richness_threshold = 10
  richness_threshold = 300
  if (abund)
  {
    read_threshold = 100
  } else if (occ)
  {
    # read_threshold = 10
    read_threshold = 1
  }
  selected_groups_indices = which(silva_phylum_mat1[,4] > richness_threshold & silva_phylum_prop_mat1[,4]*tot_read_count > read_threshold)
  silva_phylum_mat1_filtered = silva_phylum_mat1[selected_groups_indices,]
  silva_phylum_prop_mat1_filtered = silva_phylum_prop_mat1[selected_groups_indices,]
  
  pdf(paste0("Assemblage_comparison_phylum_",taxon_insert,"_",abund_insert,"_",richness_threshold,"OTU",read_threshold,"readsThres_relDiffToMeanPropRichness.pdf"))
  par(mar=c(7.1,4.1,4.1,2.1))
  silva_phylum_mat1_prop = sweep(silva_phylum_mat1_filtered,2,colSums(silva_phylum_mat1),"/")
  diffToMean = sweep(sweep(silva_phylum_mat1_prop,1,rowMeans(silva_phylum_mat1_prop),"-"),1,rowMeans(silva_phylum_mat1_prop),"/")
  for (k in 1:3)
  {
    colours = rep("blue",nrow(diffToMean))
    #diffToMean = sweep(silva_phylum_mat1,1,rowMeans(silva_phylum_mat1),"-")/silva_phylum_mat1[,k]
    colours[which(diffToMean[,k] < 0)] = "red"
    x = barplot(abs(diffToMean[,k])*100, space = 1, xaxt="n", col = colours, ylim = c(0,max(abs(diffToMean)*100)))
    labs = rownames(silva_phylum_mat1_filtered)
    #text(cex=0.8, x=x+1, y=-70, labels = labs, xpd=TRUE, srt=45, pos=2)
    text(cex=1.2, x=x+1, y=-3, labels = labs, xpd=TRUE, srt=45, pos=2)
    title(ylab = "Relative difference to the mean proportion of OTU richness (%)", main = labels[k])
  }
  dev.off()
  
  pdf(paste0("Assemblage_comparison_phylum_",taxon_insert,"_",abund_insert,"_",richness_threshold,"OTU",read_threshold,"readsThres_relDiffToMeanRichness.pdf"))
  par(mar=c(7.1,4.1,4.1,2.1))
  #silva_phylum_mat1_prop = sweep(silva_phylum_mat1,2,colSums(silva_phylum_mat1),"/")
  diffToMean = sweep(sweep(silva_phylum_mat1_filtered,1,rowMeans(silva_phylum_mat1_filtered),"-"),1,rowMeans(silva_phylum_mat1_filtered),"/")
  for (k in 1:3)
  {
    colours = rep("blue",nrow(diffToMean))
    #diffToMean = sweep(silva_phylum_mat1,1,rowMeans(silva_phylum_mat1),"-")/silva_phylum_mat1[,k]
    colours[which(diffToMean[,k] > 0)] = "red"
    x = barplot(abs(diffToMean[,k])*100, space = 1, xaxt="n", col = colours, ylim = c(0,max(abs(diffToMean)*100)))
    labs = rownames(silva_phylum_mat1_filtered)
    #text(cex=0.8, x=x+1, y=-70, labels = labs, xpd=TRUE, srt=45, pos=2)
    text(cex=1.2, x=x+1, y=-3, labels = labs, xpd=TRUE, srt=45, pos=2)
    title(ylab = "Relative difference to mean OTU richness (%)", main = labels[k])
  }
  dev.off()
  
  pdf(paste0("Assemblage_comparison_phylum_proteo_",taxon_insert,"_",abund_insert,"_",richness_threshold,"OTU",read_threshold,"readsThres_relDiffToMeanProp.pdf"))
  par(mar=c(11.1,7.1,4.1,2.1))
  par(cex.axis = 1.5)
  diffToMean = sweep(sweep(silva_phylum_prop_mat1_filtered,1,rowMeans(silva_phylum_prop_mat1_filtered),"-"),1,rowMeans(silva_phylum_prop_mat1_filtered),"/")
  for (k in 1:3)
  {
    colours = rep("blue",nrow(diffToMean))
    #diffToMean = sweep(silva_phylum_mat1,1,rowMeans(silva_phylum_mat1),"-")/silva_phylum_mat1[,k]
    colours[which(diffToMean[,k] < 0)] = "red"
    x = barplot(abs(diffToMean[,k])*100, space = 1, xaxt="n", col = colours, ylim = c(0,max(abs(diffToMean)*100)))
    labs = rownames(silva_phylum_prop_mat1_filtered)
    #text(cex=0.8, x=x+1, y=-70, labels = labs, xpd=TRUE, srt=45, pos=2)
    text(cex=1.7, x=x+0.5, y=-3, labels = labs, xpd=TRUE, srt=45, pos=2)
    title(main = paste(LETTERS[k],"-", labels[k]), cex.main = 1.7, adj = 0)
    if (k == 1)
    {
      legend(x = "topright", legend = c("Enriched wrt average","Depleted wrt average"), fill = c("blue","red"), bty = "n", cex = 1.5)
      title(ylab = "Relative difference to\n average proportion (%)", cex.lab = 1.7)
    }
  }
  dev.off()
  
  if (abund)
  {
    richness_threshold = 50
  } else if (occ)
    richness_threshold = 100
  read_threshold = 1
  selected_groups_indices = which(silva_phylum_mat1[,4] > richness_threshold & silva_phylum_prop_mat1[,4]*tot_read_count > read_threshold)
  silva_phylum_mat1_filtered = silva_phylum_mat1[selected_groups_indices,]
  
  pdf(paste0("Assemblage_comparison_phylum_",taxon_insert,"_",abund_insert,"_",richness_threshold,"OTU",read_threshold,"readsThres_diffToMeanPropRichness.pdf"))
  par(mar=c(7.1,4.1,4.1,2.1))
  silva_phylum_mat1_prop = sweep(silva_phylum_mat1_filtered,2,colSums(silva_phylum_mat1),"/")
  diffToMean = sweep(silva_phylum_mat1_prop,1,rowMeans(silva_phylum_mat1_prop),"-")
  for (k in 1:3)
  {
    colours = rep("blue",nrow(diffToMean))
    #diffToMean = sweep(silva_phylum_mat1,1,rowMeans(silva_phylum_mat1),"-")/silva_phylum_mat1[,k]
    colours[which(diffToMean[,k] < 0)] = "red"
    x = barplot(abs(diffToMean[,k]), space = 1, xaxt="n", col = colours, ylim = c(0,max(abs(diffToMean))))
    labs = rownames(silva_phylum_mat1_filtered)
    #text(cex=0.8, x=x+1, y=-70, labels = labs, xpd=TRUE, srt=45, pos=2)
    text(cex=1.2, x=x+1, y=-0.001, labels = labs, xpd=TRUE, srt=45, pos=2)
    title(ylab = "Difference to the mean proportion of OTU richness", main = labels[k])
  }
  dev.off()
  
  # richness_threshold = 1
  richness_threshold = 300
  if (abund)
  {
    read_threshold = 100000
  } else if (occ)
  {
    # read_threshold = 5000
    read_threshold = 1
  }
  selected_groups_indices = which(silva_phylum_mat1[,4] > richness_threshold & silva_phylum_prop_mat1[,4]*tot_read_count > read_threshold)
  silva_phylum_prop_mat1_filtered = silva_phylum_prop_mat1[selected_groups_indices,]
  
  pdf(paste0("Assemblage_comparison_phylum_proteo_",taxon_insert,"_",abund_insert,"_",richness_threshold,"OTU",read_threshold,"readsThres_diffToMeanProp2.pdf"))
  par(mar=c(11.1,5.1,2.1,2.1))
  par(cex.axis = 1.5)
  # diffToMean = sweep(silva_phylum_prop_mat1_filtered,1,rowMeans(silva_phylum_prop_mat1_filtered),"-")
  diffToMean = matrix(nrow = nrow(silva_phylum_prop_mat1_filtered), ncol = 3, data = 0)
  for (k in 1:3)
    diffToMean[,k] = (silva_phylum_prop_mat1_filtered[,k] - rowMeans(silva_phylum_prop_mat1_filtered[,(1:3)[-k]]))*100
  for (k in 1:3)
  {
    colours = rep("blue",length(diffToMean))
    colours[which(diffToMean[,k] < 0)] = "red"
    x = barplot(abs(diffToMean[,k]), space = 1, xaxt="n", col = colours, ylim = c(0,max(abs(diffToMean))))
    labs = rownames(silva_phylum_prop_mat1_filtered)
    #text(cex=0.8, x=x+1, y=-70, labels = labs, xpd=TRUE, srt=45, pos=2)
    text(cex=1.7, x=x+0.6, y=-0.003*100, labels = labs, xpd=TRUE, srt=45, pos=2)
    title(main = paste(LETTERS[k],"-", labels[k]), cex.main = 1.7, adj = 0)
    if (k == 1)
    {
      legend(x = 4, y = 9, legend = c("Enriched wrt other assemblages","Depleted wrt other assemblages"), fill = c("blue","red"), bty = "n", cex = 1.5)
      # title(ylab = "Difference to average proportion", cex.lab = 1.7)
      title(ylab = "Variation in proportion", cex.lab = 1.7)
    }
  }
  dev.off()
  
  # Endemism
  ###########
  
  richness_threshold = 10
  read_threshold = 100
  selected_groups_indices = which(silva_phylum_mat1[,4] > richness_threshold & silva_phylum_prop_mat1[,4]*tot_read_count > read_threshold)
  silva_phylum_endemicProp_mat1_filtered = silva_phylum_endemicProp_mat1[c(selected_groups_indices,nrow(silva_phylum_endemicProp_mat1)),]
  
  pdf(paste0("Assemblage_comparison_phylum_",taxon_insert,"_",abund_insert,"_",richness_threshold,"OTU",read_threshold,"readsThres_endemism.pdf"))
  par(mar=c(7.1,4.1,4.1,2.1))
  for (k in 1:3)
  {
    x = barplot(silva_phylum_endemicProp_mat1_filtered[,k], space = 1, xaxt="n", ylim = c(0,max(silva_phylum_endemicProp_mat1_filtered)), 
                col = c(rep("grey",nrow(silva_phylum_endemicProp_mat1_filtered)-1),"blue"))
    labs = rownames(silva_phylum_endemicProp_mat1_filtered)
    #text(cex=0.8, x=x+1, y=-70, labels = labs, xpd=TRUE, srt=45, pos=2)
    text(cex=1.2, x=x+1, y=-0.02, labels = labs, xpd=TRUE, srt=45, pos=2)
    title(ylab = "Proportion of endemic OTUs", main = labels[k])
  }
  dev.off()
  
  # OTU richness-based proportions, with intersetions across assemblages only for the assemblage-specific plots, with one colour per assemblage
  if (abund)
  {
    richness_threshold = 50
  } else if (occ)
    richness_threshold = 300
  read_threshold = 1
  selected_groups_indices = which(silva_phylum_mat1[,4] > richness_threshold & silva_phylum_prop_mat1[,4]*tot_read_count > read_threshold)
  silva_phylum_mat1_filtered = silva_phylum_mat1[selected_groups_indices,]
  silva_phylum_endemicProp_mat1_filtered = silva_phylum_endemicProp_mat1[selected_groups_indices,]
  # Computing the proportion of endemics OTUs in the non selected phyla, removing those absent in each topic:
  silva_phylum_endemicProp_mat1_others = colSums(silva_phylum_endemicProp_mat1[-selected_groups_indices,])/colSums(apply(silva_phylum_mat1[-selected_groups_indices,-4],2,function(g) g !=0))
  silva_phylum_ubiquistProp_mat1_filtered = silva_phylum_ubiquistProp_mat1[selected_groups_indices,]
  silva_phylum_ubiquistProp_mat1_others = colSums(silva_phylum_ubiquistProp_mat1[-selected_groups_indices,])/colSums(apply(silva_phylum_mat1[-selected_groups_indices,-4],2,function(g) g !=0))
  silva_phylum_semiEndemic1Prop_mat1_filtered = silva_phylum_semiEndemic1Prop_mat1[selected_groups_indices,]
  silva_phylum_semiEndemic1Prop_mat1_others = colSums(silva_phylum_semiEndemic1Prop_mat1[-selected_groups_indices,])/colSums(apply(silva_phylum_mat1[-selected_groups_indices,-4],2,function(g) g !=0))
  silva_phylum_semiEndemic2Prop_mat1_filtered = silva_phylum_semiEndemic2Prop_mat1[selected_groups_indices,]
  silva_phylum_semiEndemic2Prop_mat1_others = colSums(silva_phylum_semiEndemic2Prop_mat1[-selected_groups_indices,])/colSums(apply(silva_phylum_mat1[-selected_groups_indices,-4],2,function(g) g !=0))
  #
  silva_phylum_TFendemicProp_filtered = silva_phylum_TFendemicProp[selected_groups_indices]
  silva_phylum_HendemicProp_filtered = silva_phylum_HendemicProp[selected_groups_indices]
  silva_phylum_ERendemicProp_filtered = silva_phylum_ERendemicProp[selected_groups_indices]
  silva_phylum_TFHERProp_filtered = silva_phylum_TFHERProp[selected_groups_indices]
  silva_phylum_others = c(sum(silva_phylum_TFendemicProp[-selected_groups_indices]),sum(silva_phylum_HendemicProp[-selected_groups_indices]),
                          sum(silva_phylum_ERendemicProp[-selected_groups_indices]),sum(silva_phylum_TFHERProp[-selected_groups_indices]))/nrow(silva_phylum_mat1)
  
  #colours = c("chocolate1","dodgerblue1","darkgoldenrod1")
  colours = c("blue","chartreuse3","red")
  pdf(paste0(figure_folder,"/Assemblage_comparison_phylum_proteo_",taxon_insert,"_",abund_insert,"_",richness_threshold,"OTU",read_threshold,"readsThres_propRichness_withEndemicsSemiendemicsUbiquists2.pdf"))
  par(mar=c(11.1,5.1,4.1,2.1))
  par(cex.axis = 1.5)
  for (k in 1:4)
  {
    heights = c(silva_phylum_mat1_filtered[,k],sum(silva_phylum_mat1[,k])-sum(silva_phylum_mat1_filtered[,k]))/sum(silva_phylum_mat1[,k])*100
    if (k < 4)
      x = barplot(t(matrix(c(heights[-length(heights)]*silva_phylum_endemicProp_mat1_filtered[,k],heights[length(heights)]*silva_phylum_endemicProp_mat1_others[k],
                             # heights[-length(heights)]*(1-silva_phylum_endemicProp_mat1_filtered[,k]-silva_phylum_ubiquistProp_mat1_filtered[,k]),heights[length(heights)]*(1-silva_phylum_endemicProp_mat1_others[,k]-silva_phylum_ubiquistProp_mat1_others[,k]),
                             heights[-length(heights)]*silva_phylum_semiEndemic1Prop_mat1_filtered[,k],heights[length(heights)]*silva_phylum_semiEndemic1Prop_mat1_others[k],
                             heights[-length(heights)]*silva_phylum_semiEndemic2Prop_mat1_filtered[,k],heights[length(heights)]*silva_phylum_semiEndemic2Prop_mat1_others[k],
                             heights[-length(heights)]*silva_phylum_ubiquistProp_mat1_filtered[,k],heights[length(heights)]*silva_phylum_ubiquistProp_mat1_others[k]),ncol=4)),
                  space = 1, xaxt="n", col = t(matrix(c(rep("chartreuse3",length(heights)),rep(colours[(1:3)[-k][1]],length(heights)),rep(colours[(1:3)[-k][2]],length(heights)),rep("black",length(heights))),ncol=4)),
                  # space = 1, xaxt="n", col = t(matrix(c(rep(colours[k],length(heights)),rep("grey",length(heights)),rep("black",length(heights))),ncol=3)),
                  ylim = c(0,max(sweep(silva_phylum_mat1_filtered,2,colSums(silva_phylum_mat1),"/"))*100))
    else 
      x = barplot(t(matrix(c(heights[-length(heights)]*silva_phylum_TFendemicProp_filtered,heights[length(heights)]*silva_phylum_others[1],
                             heights[-length(heights)]*silva_phylum_HendemicProp_filtered,heights[length(heights)]*silva_phylum_others[2],
                             heights[-length(heights)]*silva_phylum_ERendemicProp_filtered,heights[length(heights)]*silva_phylum_others[3],
                             heights[-length(heights)]*(1-silva_phylum_ERendemicProp_filtered-silva_phylum_TFendemicProp_filtered-silva_phylum_HendemicProp_filtered-silva_phylum_TFHERProp_filtered),
                             heights[length(heights)]*(1 - sum(silva_phylum_others)),
                             heights[-length(heights)]*silva_phylum_TFHERProp_filtered,heights[length(heights)]*silva_phylum_others[4])*sum(silva_phylum[[4]])/100,ncol=5)),
                  space = 1, xaxt="n", col = t(matrix(c(rep(colours[1],length(heights)),rep(colours[2],length(heights)),rep(colours[3],length(heights)),rep("grey",length(heights)),rep("black",length(heights))),ncol=5)),
                  # ylim = c(0,max(sweep(silva_phylum_mat1_filtered,2,colSums(silva_phylum_mat1),"/"))*100))
                  ylim = c(0,max(sweep(silva_phylum_mat1_filtered,2,colSums(silva_phylum_mat1),"/"))*sum(silva_phylum[[4]])))
      # x = barplot(heights,
      #             space = 1, xaxt="n", col = rep("grey",length(heights)),
      #             ylim = c(0,max(sweep(silva_phylum_mat1_filtered,2,colSums(silva_phylum_mat1),"/"))*100))
    labs = c(rownames(silva_phylum_mat1_filtered),"Others")
    
    #text(cex=0.8, x=x+1, y=-70, labels = labs, xpd=TRUE, srt=45, pos=2)
    if (k == 1)
      title(ylab = "Proportion of OTU richness (%)", cex.lab = 1.7)
    if (k < 4)
    {
      title(main = paste(LETTERS[k],"-", labels[k]), cex.main = 1.7, adj = 0)
      text(cex=1.7, x=x+0.5, y=-0.01*100, labels = labs, xpd=TRUE, srt=45, pos=2)
      # legend(x = "topright", legend = c("Ubiquitous OTUs","Endemic OTUs","Other OTUs"), fill = c("blue","green","grey"), bty = "n", cex = 1.5)
      legend(y = 21, x = 6, legend = rev(c("Assemblage-specific OTUs",paste("OTUs also found in",short_titles[(1:3)[-k][1]]),paste("OTUs also found in",short_titles[(1:3)[-k][2]]),
                                        "OTUs found in all 3 assemblages")),
             fill = rev(c("chartreuse3",colours[(1:3)[-k][1]],colours[(1:3)[-k][2]],"black")), bty = "n", cex = 1.3)
    } else if (k == 4)
    {
      title(main = paste(LETTERS[1],"-", labels[k]), cex.main = 1.7, adj = 0)
      text(cex=1.7, x=x+0.5, y=-0.01*sum(silva_phylum[[4]]), labels = labs, xpd=TRUE, srt=45, pos=2)
      legend(y = 21/100*sum(silva_phylum[[4]]), x = 6, legend = rev(c("T.F.-specific OTUs","H.-specific OTUs","E.R.-specific OTUs","OTUs found in 2 assemblages","OTUs found in all 3 assemblages")),
             fill = rev(c(colours,"grey","black")), bty = "n", cex = 1.3)
      title(ylab = "OTU richness", cex.lab = 1.7)
    }
  }
  dev.off()
  
  # Occurrence-based proportions, with intersetions across assemblages only for the assemblage-specific plots, with one colour per assemblage
  # richness_threshold = 1
  # if (abund)
  # {
  #   read_threshold = 100000
  # } else if (occ)
  #   read_threshold = 5000
  if (abund)
  {
    richness_threshold = 50
  } else if (occ)
    richness_threshold = 300
  read_threshold = 1
  selected_groups_indices = which(silva_phylum_mat1[,4] > richness_threshold & silva_phylum_prop_mat1[,4]*tot_read_count > read_threshold)
  silva_phylum_prop_mat1_filtered = silva_phylum_prop_mat1[selected_groups_indices,]
  silva_phylum_prop_endemicProp_mat1_filtered = silva_phylum_prop_endemicProp_mat1[selected_groups_indices,]
  silva_phylum_prop_ubiquistProp_mat1_filtered = silva_phylum_prop_ubiquistProp_mat1[selected_groups_indices,]
  silva_phylum_prop_endemicProp_mat1_others = colSums(silva_phylum_prop_endemicProp_mat1[-selected_groups_indices,])/colSums(apply(silva_phylum_prop_mat1[-selected_groups_indices,-4],2,function(g) g !=0))
  silva_phylum_prop_ubiquistProp_mat1_filtered = silva_phylum_prop_ubiquistProp_mat1[selected_groups_indices,]
  silva_phylum_prop_ubiquistProp_mat1_others = colSums(silva_phylum_prop_ubiquistProp_mat1[-selected_groups_indices,])/colSums(apply(silva_phylum_prop_mat1[-selected_groups_indices,-4],2,function(g) g !=0))
  silva_phylum_prop_semiEndemic1Prop_mat1_filtered = silva_phylum_prop_semiEndemic1Prop_mat1[selected_groups_indices,]
  silva_phylum_prop_semiEndemic1Prop_mat1_others = colSums(silva_phylum_prop_semiEndemic1Prop_mat1[-selected_groups_indices,])/colSums(apply(silva_phylum_prop_mat1[-selected_groups_indices,-4],2,function(g) g !=0))
  silva_phylum_prop_semiEndemic2Prop_mat1_filtered = silva_phylum_prop_semiEndemic2Prop_mat1[selected_groups_indices,]
  silva_phylum_prop_semiEndemic2Prop_mat1_others = colSums(silva_phylum_prop_semiEndemic2Prop_mat1[-selected_groups_indices,])/colSums(apply(silva_phylum_prop_mat1[-selected_groups_indices,-4],2,function(g) g !=0))
  #
  silva_phylum_prop_TFendemicProp_filtered = silva_phylum_prop_TFendemicProp[selected_groups_indices]
  silva_phylum_prop_HendemicProp_filtered = silva_phylum_prop_HendemicProp[selected_groups_indices]
  silva_phylum_prop_ERendemicProp_filtered = silva_phylum_prop_ERendemicProp[selected_groups_indices]
  silva_phylum_prop_TFHERProp_filtered = silva_phylum_prop_TFHERProp[selected_groups_indices]
  silva_phylum_prop_others = c(sum(silva_phylum_prop_TFendemicProp[-selected_groups_indices]),sum(silva_phylum_prop_HendemicProp[-selected_groups_indices]),
                          sum(silva_phylum_prop_ERendemicProp[-selected_groups_indices]),sum(silva_phylum_prop_TFHERProp[-selected_groups_indices]))/nrow(silva_phylum_prop_mat1)
  
  # colours = c("chocolate1","dodgerblue1","darkgoldenrod1")
  colours = c("blue","chartreuse3","red")
  pdf(paste0(figure_folder,"/Assemblage_comparison_phylum_proteo_",taxon_insert,"_",abund_insert,"_",richness_threshold,"OTU",read_threshold,"readsThres_prop_withEndemicsSemiendemicsUbiquists3.pdf"))
  par(mar=c(11.1,5.1,4.1,2.1))
  par(cex.axis = 1.5)
  for (k in 1:4)
  {
    heights = c(silva_phylum_prop_mat1_filtered[,k],1-sum(silva_phylum_prop_mat1_filtered[,k]))*100
    if (k < 4)
      x = barplot(t(matrix(c(heights[-length(heights)]*silva_phylum_prop_endemicProp_mat1_filtered[,k],heights[length(heights)]*silva_phylum_prop_endemicProp_mat1_others[k],
                             heights[-length(heights)]*(1-silva_phylum_prop_endemicProp_mat1_filtered[,k]-silva_phylum_prop_ubiquistProp_mat1_filtered[,k]),heights[length(heights)]*(1-silva_phylum_prop_endemicProp_mat1_others[k]-silva_phylum_prop_ubiquistProp_mat1_others[k]),
                             # heights[-length(heights)]*silva_phylum_prop_semiEndemic1Prop_mat1_filtered[,k],heights[length(heights)]*silva_phylum_prop_semiEndemic1Prop_mat1_others[k],
                             # heights[-length(heights)]*silva_phylum_prop_semiEndemic2Prop_mat1_filtered[,k],heights[length(heights)]*silva_phylum_prop_semiEndemic2Prop_mat1_others[k],
                             heights[-length(heights)]*silva_phylum_prop_ubiquistProp_mat1_filtered[,k],heights[length(heights)]*silva_phylum_prop_ubiquistProp_mat1_others[k]),ncol=3)),
                  # space = 1, xaxt="n", col = t(matrix(c(rep("green",length(heights)),rep(colours[(1:3)[-k][1]],length(heights)),rep(colours[(1:3)[-k][2]],length(heights)),rep("black",length(heights))),ncol=4)),
                  space = 1, xaxt="n", col = t(matrix(c(rep(colours[k],length(heights)),rep("grey",length(heights)),rep("black",length(heights))),ncol=3)),
                  ylim = c(0,max(silva_phylum_prop_mat1_filtered))*100)
    else 
      x = barplot(t(matrix(c(heights[-length(heights)]*silva_phylum_prop_TFendemicProp_filtered,heights[length(heights)]*silva_phylum_prop_others[1],
                             heights[-length(heights)]*silva_phylum_prop_HendemicProp_filtered,heights[length(heights)]*silva_phylum_prop_others[2],
                             heights[-length(heights)]*silva_phylum_prop_ERendemicProp_filtered,heights[length(heights)]*silva_phylum_prop_others[3],
                             heights[-length(heights)]*(1 - silva_phylum_prop_ERendemicProp_filtered-silva_phylum_prop_TFendemicProp_filtered-silva_phylum_prop_HendemicProp_filtered-silva_phylum_prop_TFHERProp_filtered),
                             heights[length(heights)]*(1 - sum(silva_phylum_prop_others)),
                             heights[-length(heights)]*silva_phylum_prop_TFHERProp_filtered,heights[length(heights)]*silva_phylum_prop_others[4]),ncol=5)),
                  space = 1, xaxt="n", col = t(matrix(c(rep(colours[1],length(heights)),rep(colours[2],length(heights)),rep(colours[3],length(heights)),rep("grey",length(heights)),rep("black",length(heights))),ncol=5)),
                  ylim = c(0,max(silva_phylum_prop_mat1_filtered))*100)
    
    labs = c(rownames(silva_phylum_prop_mat1_filtered),"Others")
    # text(cex=1.2, x=x+0.5, y=-0.01, labels = labs, xpd=TRUE, srt=45, pos=2)
    text(cex=1.7, x=x+0.5, y=-0.01*100, labels = labs, xpd=TRUE, srt=45, pos=2)
    title(main = paste(LETTERS[k+1],"-", labels[k]), cex.main = 1.7, adj = 0)
    # if (k == 1)
    # {
    if (abund)
    {
      title(ylab = "Proportion of read count (%)", cex.lab = 1.7)
    } else if (occ)
      title(ylab = "Proportion (%)", cex.lab = 1.7)
    # }
    # if (k < 4)
    # {
    #   # legend(x = "topright", legend = c("Ubiquitous OTUs","Endemic OTUs","Other OTUs"), fill = c("blue","green","grey"), bty = "n", cex = 1.5)
    #   legend(y = 25, x = 6, legend = rev(c("Assemblage-specific OTUs",paste("OTUs also found in",short_titles[(1:3)[-k][1]]),paste("OTUs also found in",short_titles[(1:3)[-k][2]]),
    #                                        "OTUs found in all 3 assemblages")),
    #          fill = rev(c("green",colours[(1:3)[-k][1]],colours[(1:3)[-k][2]],"black")), bty = "n", cex = 1.5)
    # } else if (k == 4)
    # {
    #   legend(y = 25, x = 6, legend = rev(c("T.F.-specific OTUs","H.-specific OTUs","E.R.-specific OTUs","OTUs found in 2 assemblages","OTUs found in all 3 assemblages")),
    #          fill = rev(c(colours,"grey","black")), bty = "n", cex = 1.5)
    # }
  }
  dev.off()
  
  # Occurrence-based proportions, with intersetions across assemblages on all plots, with one colour per intersection
  ################
  if (abund)
  {
    richness_threshold = 50
  } else if (occ)
    richness_threshold = 300
  read_threshold = 1
  selected_groups_indices = which(silva_phylum_mat1[,4] > richness_threshold & silva_phylum_prop_mat1[,4]*tot_read_count > read_threshold)
  silva_phylum_prop_mat1_filtered = silva_phylum_prop_mat1[selected_groups_indices,]
  silva_phylum_prop_endemicProp_mat1_filtered = silva_phylum_prop_endemicProp_mat1[selected_groups_indices,]
  silva_phylum_prop_ubiquistProp_mat1_filtered = silva_phylum_prop_ubiquistProp_mat1[selected_groups_indices,]
  silva_phylum_prop_endemicProp_mat1_others = colSums(silva_phylum_prop_endemicProp_mat1[-selected_groups_indices,])/colSums(apply(silva_phylum_prop_mat1[-selected_groups_indices,-4],2,function(g) g !=0))
  silva_phylum_prop_ubiquistProp_mat1_filtered = silva_phylum_prop_ubiquistProp_mat1[selected_groups_indices,]
  silva_phylum_prop_ubiquistProp_mat1_others = colSums(silva_phylum_prop_ubiquistProp_mat1[-selected_groups_indices,])/colSums(apply(silva_phylum_prop_mat1[-selected_groups_indices,-4],2,function(g) g !=0))
  silva_phylum_prop_semiEndemic1Prop_mat1_filtered = silva_phylum_prop_semiEndemic1Prop_mat1[selected_groups_indices,]
  silva_phylum_prop_semiEndemic1Prop_mat1_others = colSums(silva_phylum_prop_semiEndemic1Prop_mat1[-selected_groups_indices,])/colSums(apply(silva_phylum_prop_mat1[-selected_groups_indices,-4],2,function(g) g !=0))
  silva_phylum_prop_semiEndemic2Prop_mat1_filtered = silva_phylum_prop_semiEndemic2Prop_mat1[selected_groups_indices,]
  silva_phylum_prop_semiEndemic2Prop_mat1_others = colSums(silva_phylum_prop_semiEndemic2Prop_mat1[-selected_groups_indices,])/colSums(apply(silva_phylum_prop_mat1[-selected_groups_indices,-4],2,function(g) g !=0))
  #
  silva_phylum_prop_TFendemicProp_filtered = silva_phylum_prop_TFendemicProp[selected_groups_indices]
  silva_phylum_prop_HendemicProp_filtered = silva_phylum_prop_HendemicProp[selected_groups_indices]
  silva_phylum_prop_ERendemicProp_filtered = silva_phylum_prop_ERendemicProp[selected_groups_indices]
  silva_phylum_prop_TFERendemicProp_filtered = silva_phylum_prop_TFERendemicProp[selected_groups_indices]
  silva_phylum_prop_HERendemicProp_filtered = silva_phylum_prop_HERendemicProp[selected_groups_indices]
  silva_phylum_prop_TFHendemicProp_filtered = silva_phylum_prop_TFHendemicProp[selected_groups_indices]
  silva_phylum_prop_TFHERProp_filtered = silva_phylum_prop_TFHERProp[selected_groups_indices]
  silva_phylum_prop_others = c(sum(sum(silva_phylum_prop_TFendemicProp[-selected_groups_indices]),sum(silva_phylum_prop_HendemicProp[-selected_groups_indices]),sum(silva_phylum_prop_ERendemicProp[-selected_groups_indices])),
                               sum(silva_phylum_prop_TFERendemicProp[-selected_groups_indices]),sum(silva_phylum_prop_HERendemicProp[-selected_groups_indices]),
                               sum(silva_phylum_prop_TFHendemicProp[-selected_groups_indices]),sum(silva_phylum_prop_TFHERProp[-selected_groups_indices]))/nrow(silva_phylum_prop_mat1)
  
  # TFH, HER, TFER
  # colours = c("dodgerblue1","darkgoldenrod1","chocolate1")
  colours = c("blue","chartreuse3","red")
  pdf(paste0(figure_folder,"/Assemblage_comparison_phylum_proteo_",taxon_insert,"_",abund_insert,"_",richness_threshold,"OTU",read_threshold,"readsThres_prop_withEndemicsSemiendemicsUbiquists2.pdf"))
  par(mar=c(11.1,5.1,4.1,2.1))
  par(cex.axis = 1.5)
  for (k in 1:4)
  {
    heights = c(silva_phylum_prop_mat1_filtered[,k],1-sum(silva_phylum_prop_mat1_filtered[,k]))*100
    if (k < 4)
      x = barplot(t(matrix(c(heights[-length(heights)]*silva_phylum_prop_endemicProp_mat1_filtered[,k],heights[length(heights)]*silva_phylum_prop_endemicProp_mat1_others[k],
                             #heights[-length(heights)]*(1-silva_phylum_endemicProp_mat1_filtered[,k]-silva_phylum_ubiquistProp_mat1_filtered[,k]),heights[length(heights)],
                             heights[-length(heights)]*silva_phylum_prop_semiEndemic1Prop_mat1_filtered[,k],heights[length(heights)]*silva_phylum_prop_semiEndemic1Prop_mat1_others[k],
                             heights[-length(heights)]*silva_phylum_prop_semiEndemic2Prop_mat1_filtered[,k],heights[length(heights)]*silva_phylum_prop_semiEndemic2Prop_mat1_others[k],
                             heights[-length(heights)]*silva_phylum_prop_ubiquistProp_mat1_filtered[,k],heights[length(heights)]*silva_phylum_prop_ubiquistProp_mat1_others[k]),ncol=4)),
                  space = 1, xaxt="n", col = t(matrix(c(rep("green",length(heights)),rep(colours[c(1,1,3)[k]],length(heights)),rep(colours[c(3,2,2)[k]],length(heights)),rep("black",length(heights))),ncol=4)),
                  ylim = c(0,max(silva_phylum_prop_mat1_filtered))*100)
    else 
      x = barplot(t(matrix(c(heights[-length(heights)]*(silva_phylum_prop_ERendemicProp_filtered + silva_phylum_prop_TFendemicProp_filtered + silva_phylum_prop_HendemicProp_filtered),heights[length(heights)]*silva_phylum_prop_others[1],
                             heights[-length(heights)]*silva_phylum_prop_TFHendemicProp_filtered,heights[length(heights)]*silva_phylum_prop_others[2],
                             heights[-length(heights)]*silva_phylum_prop_HERendemicProp_filtered,heights[length(heights)]*silva_phylum_prop_others[3],
                             heights[-length(heights)]*silva_phylum_prop_TFERendemicProp_filtered,heights[length(heights)]*silva_phylum_prop_others[4],
                             heights[-length(heights)]*silva_phylum_prop_TFHERProp_filtered,heights[length(heights)]*silva_phylum_prop_others[5]),ncol=5)),
                  space = 1, xaxt="n", col = t(matrix(c(rep("chartreuse3",length(heights)),rep(colours[1],length(heights)),rep(colours[2],length(heights)),rep(colours[3],length(heights)),rep("black",length(heights))),ncol=5)),
                  ylim = c(0,max(silva_phylum_prop_mat1_filtered))*100)
    
    labs = c(rownames(silva_phylum_prop_mat1_filtered),"Others")
    # text(cex=1.2, x=x+0.5, y=-0.01, labels = labs, xpd=TRUE, srt=45, pos=2)
    text(cex=1.7, x=x+0.5, y=-0.01*100, labels = labs, xpd=TRUE, srt=45, pos=2)
    title(main = paste(LETTERS[k],"-", labels[k]), cex.main = 1.7, adj = 0)
    if (k == 1)
    {
      if (abund)
        title(ylab = "Proportion of read count (%)", cex.lab = 1.7)
      else if (occ)
        title(ylab = "Proportion of occurrences (%)", cex.lab = 1.7)
    }
    if (k < 4)
    {
      # legend(x = "topright", legend = c("Ubiquitous OTUs","Endemic OTUs","Other OTUs"), fill = c("blue","green","grey"), bty = "n", cex = 1.5)
      legend(y = 25, x = 6, legend = rev(c("Assemblage-specific OTUs",paste("OTUs found in",short_titles[k],"and",short_titles[(1:3)[-k][1]]),paste("OTUs found in",short_titles[k],"and",short_titles[(1:3)[-k][2]]),
                                           "OTUs found in all 3 assemblages")),
             fill = rev(c("chartreuse3",colours[c(1,1,3)[k]],colours[c(3,2,2)[k]],"black")), bty = "n", cex = 1.5)
    } else if (k == 4)
    {
      legend(y = 25, x = 6, legend = rev(c("Assemblage-specific OTUs",paste("OTUs found in",short_titles[1],"and",short_titles[2]),paste("OTUs found in",short_titles[2],"and",short_titles[3]),
                                           paste("OTUs found in",short_titles[1],"and",short_titles[3]),"OTUs found in all 3 assemblages")),
             fill = rev(c("chartreuse3",colours,"black")), bty = "n", cex = 1.5)
    }
  }
  dev.off()
  
  # pdf("Functional_assignations.pdf")
  # par(mar=c(11.1,4.1,4.1,2.1))
  # x = barplot(table(data_Federico_selected[[3]]$Function), space = 1, xaxt="n")
  # labs = names(table(data_Federico_selected[[3]]$Function))
  # #text(cex=0.8, x=x+1, y=-70, labels = labs, xpd=TRUE, srt=45, pos=2)
  # text(cex=1.2, x=x+0.5, y=-1000, labels = labs, xpd=TRUE, srt=45, pos=2)
  # title(ylab = "Number of OTUs")
  # dev.off()
  
}


