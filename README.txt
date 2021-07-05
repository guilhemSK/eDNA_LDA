The eDNA_LDA repository contains all the code used for the two following papers:
- "Latent Dirichlet Allocation reveals spatial and taxonomic structure in a DNA-based census of soil biodiversity from a tropical forest", Sommeria-Klein et al. MER 2020
- "Global drivers of eukaryotic plankton biogeography in the sunlit ocean", Sommeria-Klein et al. bioRxiv 2020
More generally, this code is aimed at being used for the analysis of environmental DNA data using Latent Dirichlet Allocation. Accordingly, it is provided with a range of options.

The repository consists of a main script LDA_main_v6.2.R, of a parameter panel LDA_parameters.R, and of several functions called from the main script. Documentation and parameter settings are found in the LDA_parameters.R file. The folders NouraguesPetitPlateau_LDA and Tara_LDA contain code more specifically related to the corresponding papers.

The Test_data_protists_PP folder contains the properly formatted protist data of the Petit Plateau dataset. To test the code on these data with default settings, run the following command from the eDNA_LDA directory:
> Rscript --vanilla LDA_main_v6.2.R ./Test_data_protists_PP LDA_parameters.R

