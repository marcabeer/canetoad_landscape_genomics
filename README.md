# canetoad_landscape_genomics
This repository contains scripts and data underlying population genomics analyses for a ddRADseq dataset of cane toads across their invasive Australian range

The two major subdirectories contain the following:
- data
  - ct.vcf.gz (a gzipped VCF file of filtered genotypes for 932 cane toads and 5723 SNPs)
  - metadata_individual (a tab-delimited file giving geographic coordinates of individuals)

- analyses
	- colonization year (analysis of occurrence records to reconstruct the cane toad's spread across Australia)
	- env_data (processing of environmental data for use in genetic-environment association tests)
	- fig1_sampling_map (scripts and data underlying the creation of Fig.1 in the main text of this project's peer-reviewed publication)
	- gea_tests (analysis of signatures of selection in the form of genetic-environment associations)
	- genetic_diversity (analysis of various metrics concerning genetic variation)
	- pop_struct (analysis of population genetic differentiation among geographic sampling units)