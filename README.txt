R code associated with the submitted manuscript 'Latent Dirichlet Allocation of a DNA-based census reveals the spatial and taxonomic structure of soil biodiversity in a tropical forest'.

The eDNA_LDA repository contains the main code used to obtain the results presented in the paper. The code consists of a main part LDA_main_v5.R, of a parameter panel LDA_parameters.R, and of several functions called from the main code. Documentation and parameter settings are found in the LDA_parameters.R file. This code is aimed at being used for the analysis of environmental DNA data beyond the Petit Plateau dataset explored in the manuscript. Accordingly, it is provided with a range of options.

The Test_data_protists_PP folder contains the properly formatted protist data of the Petit Plateau dataset. To test the code on these data with default settings, run the following command from the eDNA_LDA directory:
> Rscript --vanilla LDA_main_v5.R LDA_parameters.R .
