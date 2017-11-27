
## *Clostridium difficile* differentially alters the structure and metabolism of distinct cecal microbiomes to promote sustained colonization during infection

### Abstract

*Clostridium difficile* has become the most common single cause of hospital-acquired infection over the last decade in the United States. Susceptibility is primarily associated with previous exposure to antibiotics, which compromise the structure and function of the gut bacterial community. Furthermore, specific classes correlate more strongly with recurrent or persistent *C. difficile* infection. We used a murine model of infection to explore the effect of distinct antibiotic classes on sustained *C. difficile* colonization, as well as the impact of infection on community-level gene expression and metabolism 18 hours post-infection. Utilizing untargeted metabolomic analysis, we found that *C. difficile* infection has larger impacts on the metabolic activity of the microbiota across cefoperazone and streptomycin-treated mice, which become persistently colonized. Using metagenome-enabled metatranscriptomics, we observed that the infected communities were enriched in pathways associated with amino acid metabolism and particularly in non-dominant species relative to mock-infected controls. Conversely, in clindamycin pretreatment where *C. difficile* is cleared within 8 days, the effect of infection on the microbiota was only detectable in changes to the community structure but not in metabolic activity or gene expression. Our results suggest that *C. difficile* is able to restructure the nutrient-niche landscape in certain gut environments in order to promote persistent infection.

### Summary

Colonization resistance to the nosocomial pathogen *Clostridium difficile* is primarily driven by the gut microbiota. When the intact community of bacteria in the gastrointestinal tract is disrupted by treatments like antibiotics for previous infections, *C. difficile* is subsequently permitted to colonize and cause disease. Further complicating matters, the microbiota is also involved in clearing the infection as the community recovers from perturbation. As distinct antibiotics are associated with different risk levels for CDI, we utilized a mouse model of infection with 3 separate antibiotic pretreatment regimes to generate alternative gut microbiomes that each allowed for *C. difficile* colonization but vary in clearance rate. To assess community-level dynamics, we implemented a multi-omic approach by integrating both metatranscriptomics and untargeted metabolomics that revealed infection significantly shifted many aspects of the gut ecosystem. Additionally, the degree to which this change occurred inversely correlated with clearance during the first six days of infection. Following targeted analysis of Stickland fermentation byproducts, we found that *C. difficile* may differentially modify the gut environment to promote persistence based on availability of a preferred nutrient niche. Our results improve understanding of the ecology associated with *C. difficile* infection and provides groundwork for identification of context-specific probiotic therapies.


### Overview

	project
	|- README          # the top level description of content
	|- LICENSE         # the license for this project
	|
	|- Jenior_Metatransciptomics_PLOSPathogens_2017.Rmd		# executable Rmarkdown for this study
	|- Jenior_Metatransciptomics_PLOSPathogens_2017.md		# Markdown (GitHub) version of the *.Rmd file
	|- Jenior_Metatransciptomics_PLOSPathogens_2017.docx		# Word document of study manuscript
	|- manuscript_format.docx				# Word document associated with Rmd for proper formatting
	|- references.bib					# BibTeX formatted references
	|- plos-pathogens.csl						# csl file to format references
	|
	|- doc/		# Metagenome & metatranscriptome bioinformatic processing pipelines
	|
	|- data/          	# raw and primary data
	| |- /16S_analysis			# 
	| |- /kegg			# subset of KEGG reference files used
	| |- /metabolome			# results from untargeted LCMS analysis
	| |- /metagenomes			# metagenomic contigs
	| |- /read_mapping			# normalized metagenomic and mettranscriptomic read abundances
	| |- /old_versions			# legacy versions of wetlab analyses
	| |- /assembly_stats			# quality metrics for metagenomic contigs
	| |- /cfu_time			# cfu counts over time for each antibiotic
	| |- /metadata			# experimental metadata from all animal experiments
	| |- /taxonomy_color			# colors used for taxonomic groups for OTU bar plots
	| +- /wetlab_assays		# cfu counts and toxin titers for all pretreatment groups
	|
	|- code/	# all programmatic code (python, R, pbs, mothur_batch)
	|
	|- results/		# all output from workflows and analyses
	| |- figures/		# manuscript figures
	| |- /table_1			# antibiotic regime information
	| +- supplement/	# supplementary tables and figures

