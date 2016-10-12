#####################
#####################
## PARAMETER PANEL ##
#####################
#####################

# Where the code is run from:
locally_based = 1
genotoul_cluster_based = 0
EDB_cluster_based = 0

# Local LDA calculation:
local = 0
# Existing LDA decomposition stored locally:
existingresult = 1
# Existing LDA decomposition stored on the Genotoul cluster, only available for locally_based = 1:
cluster = 0

if (local)
{
  start.time <- Sys.time()
  
  library(topicmodels)
}

data_h20 = 0
data_pp = 1
# data_gs = 1 requires filled = 1:
data_gs = 0
filled = 1
data_betadiv = 0
data_betadiv_pooled = 0

# Treating the data as occurrence data instead of abundance data
occurrence = 1

# Number of sites an OTU must ocuppy to be kept
# (if nb_occupied_sites_threshold = 1, all OTUs with non-zero abundance are kept)
nb_occupied_sites_threshold = 1
# Removing OTUs with less reads than the number of sites
no_rare = 0

barcode_gh = 0
barcode_ghassigned = 0
barcode_itsfungi = 0
barcode_itsasco = 0
barcode_itsbasidio = 0
barcode_itsglomero = 0
bacteriafungiarchaea = 0
barcode_16sbact = 1
barcode_18s = 0
barcode_18sfungi = 0
barcode_18smetazoa = 0
barcode_18sannelids = 0
barcode_18sarthropods = 0
barcode_18snematodes = 0
barcode_18sprotists = 0
barcode_18splatyhelminthes = 0
barcode_18splants = 0
barcode_16sarch = 0
barcode_itsplant = 0
barcode_16sins = 0
plantfungi = 0
#frogs
frogs_iucn_all = 0
frogs_iucn_forest = 0
frogs_iucn_open = 0
frogs_poly1_all = 0
frogs_poly1_forest = 0
frogs_poly1_open = 0
frogs_poly2_all = 0
frogs_poly2_forest = 0
frogs_poly2_open = 0   

# testdata = 1 requires data_pp = 1 and filled = 1
testdata = 0
testdata_dir = "Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics0.02_1000sampledreads"
#testdata_dir = "Continuous-mixed_samples_nbtopics5_nbmotus1000_dirichlettopics0.02_1000sampledreads_lnornoise_sig0.5"
#testdata_dir = "Continuous-mixed_samples_nbtopics5_nbmotus1000_randomtopics_sampledreads"
# testdata_dir = "Continuous-mixed_samples_nbtopics5_nbmotus1000_randomtopics_sampledreads_lnornoise_sig1"
nb_rndzations_true_documents = 10000
# The baseline proportion of the topics in true_KL_documents in samples where they are actuallt absent
# (added due to KL.pugin() taking only non-zero values as input)
#true_KL_documents_baseline = 10^-6
#true_KL_documents_baseline = 1/1131000

blei = 0
Rtopicmodels_Gibbs = 0
Rtopicmodels_VEM = 1
hdp_faster = 0

mpar = 0
nb_topics = 3
# For testdata only :
true_nb_topics = 5

#nb_topics_range = c(2,3)
nb_topics_range = c(2,3,4,5,6,7)
#nb_topics_range = c(2,3,4,5,6,7,8,9,10,11,13,15)
#nb_topics_range = c(2,3,4,5,7,9,11,13,15,17,19)
#nb_topics_range = c(2,3,4,5,10,15,20,30,40,50)
#nb_topics_range = c(2,3,4,5,7,9,11,13,15,17,19,21,23,25,27,29,31)
#nb_topics_range = c(2,3,4,5,6,7,8,9,10,11,12,13,15,17,19)
#nb_topics_range = c(30,40,50,60,70)
# mnb_topics = 1 only possible if mpar = 1
mnb_topics = 1

if (mpar) {
  if (mnb_topics)
    mpar_range = length(nb_topics_range)
  #else if (miter)
  #mpar_range = length(miter_range)
} else mpar_range=1

delta = 0.1
#alpha_insert = "alpha50:nb_topics"
#alpha_insert = "alpha50_over_nb_topics"
alpha_insert = "alpha0.1"

nb_real = 100
# if (mpar && mnb_topics) or if (!mpar && best_keep), we must have best=0 
# best_keep = 1 has no influence if (mpar && mnb_topics)
best = 0
best_keep = 1
# select_real=1 only possible if best_keep=1
select_real = 1
# Selected_real must start at 1
#Selected_real = c(seq(1,10,1),seq(91,100,1))
#Selected_real = seq(1,100,1)
#Selected_real = seq(1,5,1)
Selected_real = 1

# only useful for Gibbs sampling
nb_iter = 2000
llh_keep = 100

# only useful for VEM
#em_tol = 5*(10^-5)
em_tol = 10^-7
var_tol = 10^-8

# Calculation associated to topic distribution and associated plots
dominance_calculation = 0
# Calculating the proportion of the total read number taken by each topic (useful for ordering the topics as an alternative to normalized abundance)
topic_read_proportion_comput = 0

# Comparing topics to chemistery and lidar data (for data_pp only):
abiotic_variables = 0
# Number of randomizations for lidar data, required for abiotic_variables = 1
# Randomizations not implemented yet for data_gs
nb_abiotic_rndzations = 100000
# Perform a PCA over the abiotic variables before computing the correlation to the assemblages (only useful for abiotic_variables = 1) :
pca_abiotic = 0

# Comparing realizations to each other
realization_comparison = 0
# Comparing the realizations MOTUwise or samplewise (with spatial randomizations)
# For testdata = 1, MOTUwise = 1 would require loading true_topic_compo (not done yet)
samplewise = 1
MOTUwise = 0
# How the topic correspondence between realizations is computed: in a bijective way if bij=1 or taking into account the 3 lowest SKL if bij=0 
# (only useful if realization_comparison=1)
bij = 1
# Number of randomizations in computing the distance between realizations (only useful if realization_comparison=1)
nb_rndzations = 1000

# Comparing topics within the best realization:
best_real_comparison = 1
# Number of randomization in comparing the topics to each other
nb_rndzations_best_real = 10000

# kriging for the spatial distribution of topics in the best realization
kriging = 0
# required if kriging = 1, only the kriged_real are kriged (all kriged realizations must be included in Selected_real)
#kriged_real = seq(1,5,1)
kriged_real = c(54,55)

#########################
#########################
## END PARAMETER PANEL ##
#########################
#########################