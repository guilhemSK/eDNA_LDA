
R code and input and output data associated with the article 'Global drivers of eukaryotic plankton biogeography in the sunlit ocean', Sommeria-Klein et al. Science 2021.

Formatted.input.data_and_output.biogeographies contains, bundled together in a zip file, the input data (properly formatted for the 'LDA_main' R code used to generate the results) and the output biogeographies (as data frames). The file tree structure is as outputted by the code: there is one folder per DNA marker and major plankton group, which contains:
- the input community matrix ('data2m.Rdata'),
- the coordinates of the corresponding Tara stations ('coord.Rdata'), 
- taxonomic information for the corresponding OTUs ('taxo_ref.Rdata'), 
and within a subfolder:
- the spatial distribution of the assemblages outputted by Latent Dirichlet Allocation (i.e., the group's biogeography, 'Spatial_topicmix_kriged.rds', which does not involve kriging despite the file name), 
- and the OTU composition of the assemblages ('assemblage_composition.rds'). 

In the case of the 18S-V9 marker the folder '18S_V9_TARA_CompleteSizeRange_byStationByDepth_AllTaxa_noLagoon' contains the same files but for the whole dataset, as used in Figure 1 of the article. Finally, the 'Abiotic_data' folder contains the World Ocean Atlas 2013 data ('woa13_env_tara_all.csv') and the simulated minimum transport times by currents (at 5 m depth, 'minAijji_tarrive_min_surface_10.csv', and at 75 m depth, 'tarrive_min_75m_10.csv') used in the article.  

In addition, the file "Biogeographies18SV9_assemblageByAssemblage_10majorGroups.pdf" shows assemblage by assemblage representations of the biogeographies illustrated in Figures 2 and S5 for the sake of clarity.

'Multigroup_Tara_analyses.R' contains all the scripts used to generate the results presented in the article and its Supplementary Materials based on the Latent Dirichlet Allocation model outputs mentioned above. The corresponding results are provided in the 'Saved_results' folder. 'Tara_figures.R' contains all the scripts used to generate the manuscript's main and supplementary figures from the results in 'Saved_results'. 'Normalized.VI_fun.R' is a function called by 'Multigroup_Tara_analyses.R'.

Suitable absolute paths to data, code and figure folders need to be defined at the beginning of 'Multigroup_Tara_analyses.R' (lines 2-14) and 'Tara_figures.R' (lines 69-77). In particular, 'code_folder' should point to the folder to which the eDNA_LDA repository has been cloned, and 'data_folder' (as well as all 'data_folder_workspacex' variables) should point to the folder where the content of 'Formatted.input.data_and_output.biogeographies' has been downloaded and expanded.

