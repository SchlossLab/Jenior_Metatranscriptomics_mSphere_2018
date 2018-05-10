
## *Clostridium difficile* differentially alters the structure and metabolism of distinct cecal microbiomes to promote sustained colonization during infection

### Abstract

Susceptibility to *Clostridium difficile* infection is primarily associated with previous exposure to antibiotics, which compromise the structure and function of the gut bacterial community. Specific antibiotic classes correlate more strongly with recurrent or persistent *C. difficile* infection. As such, we utilized a mouse model of infection to explore the effect of distinct antibiotic classes on the impact that infection has on community-level transcription and metabolic signatures shortly following pathogen colonization and how those changes may associate with persistence of *C. difficile*. Untargeted metabolomic analysis revealed that *C. difficile* infection had significantly larger impacts on the metabolic environment across cefoperazone and streptomycin-pretreated mice, which become persistently colonized compared to clindamycin-pretreated mice where infection quickly became undetectable. Through metagenome-enabled metatranscriptomics we observed that transcripts for genes associated with carbon and energy acquisition were greatly reduced in infected animals, suggesting those niches were instead occupied by *C. difficile*. Furthermore, the largest changes in transcription were seen in the least abundant species indicating that *C. difficile* may "attack the loser" in gut environments where sustained infection occurs more readily. Overall, our results suggest that *C. difficile* is able to restructure the nutrient-niche landscape in the gut to promote persistent infection.

#### Importance

*Clostridium difficile* has become the most common single cause of hospital-acquired infection over the last decade in the United States. Colonization resistance to the nosocomial pathogen is primarily provided by the gut microbiota, which is also involved in clearing the infection as the community recovers from perturbation. As distinct antibiotics are associated with different risk levels for CDI, we utilized a mouse model of infection with 3 separate antibiotic pretreatment regimes to generate alternative gut microbiomes that each allowed for *C. difficile* colonization but varied in clearance rate. To assess community-level dynamics, we implemented an integrative multi-omic approach that revealed that infection significantly changed many aspects of the gut community. The degree to which the community changed was inversely correlated with clearance during the first six days of infection, suggesting that *C. difficile* differentially modifies the gut environment to promote persistence. This is the first time metagenome-enabled metatranscriptomics have been employed to study the behavior of a host-associated microbiota in response to an infection. Our results allow for a previously unseen understanding of the ecology associated with *C. difficile* infection and provides groundwork for identification of context-specific probiotic therapies.


### Overview

	project
	|- README          # the top level description of content
	|- LICENSE         # the license for this project
	|
	|- Jenior_Metatransciptomics_mSphere_2018.Rmd		# executable Rmarkdown for this study
	|- Jenior_Metatransciptomics_mSphere_2018.md		# Markdown (GitHub) version of the *.Rmd file
	|- Jenior_Metatransciptomics_mSphere_2018.docx		# Word document of study manuscript
	|- manuscript_format.docx				# Word document associated with Rmd for proper formatting
	|- references.bib					# BibTeX formatted references
	|- msphere.csl						# csl file to format references
	|
	|- doc/		# Metagenome & metatranscriptome bioinformatic processing pipelines
	|
	|- data/          		# raw and primary data
	| |- /16S_analysis		# 16S results for all pretreatments
	| |- /kegg			# subset of KEGG reference files used
	| |- /metabolome		# results from untargeted LCMS analysis
	| |- /metagenomes		# metagenomic contigs
	| |- /read_mapping		# normalized metagenomic and mettranscriptomic read abundances
	| |- /old_versions		# legacy versions of wetlab analyses
	| |- /assembly_stats		# quality metrics for metagenomic contigs
	| |- /cfu_time			# cfu counts over time for each antibiotic
	| |- /metadata			# experimental metadata from all animal experiments
	| |- /taxonomy_color		# colors used for taxonomic groups for OTU bar plots
	| +- /wetlab_assays		# cfu counts and toxin titers for all pretreatment groups
	|
	|- code/	# all programmatic code (python, R, pbs, mothur_batch)
	|
	|- results/		# all output from workflows and analyses
	| |- figures/		# manuscript figures
	| +- supplement/	# supplementary tables and figures

