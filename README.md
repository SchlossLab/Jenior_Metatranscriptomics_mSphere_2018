
## *Clostridium difficile* differentially alters the structure and metabolism of distinct cecal microbiomes to promote sustained colonization during infection


*Clostridium difficile* has become the most common single cause of hospital-acquired infection over the last decade in the United States. Susceptibility is primarily associated with previous exposure to antibiotics, which compromise the structure and function of the gut bacterial community. Furthermore, specific classes correlate more strongly with recurrent or persistent *C. difficile* infection. We used a murine model of infection to explore the effect of distinct antibiotic classes on sustained *C. difficile* colonization, as well as the impact of infection on community-level gene expression and metabolism 18 hours post-infection. Utilizing untargeted metabolomic analysis, we found that *C. difficile* infection has larger impacts on the metabolic activity of the microbiota across cefoperazone and streptomycin-treated mice, which become persistently colonized. Using metagenome-enabled metatranscriptomics, we observed that the infected communities were enriched in pathways associated with amino acid metabolism and particularly in non-dominant species relative to mock-infected controls. Conversely, in clindamycin pretreatment where *C. difficile* is cleared within 8 days, the effect of infection on the microbiota was only detectable in changes to the community structure but not in metabolic activity or gene expression. Our results suggest that *C. difficile* is able to restructure the nutrient-niche landscape in certain gut environments in order to promote persistent infection.


### Overview

	project
	|- README          # the top level description of content
	|- LICENSE         # the license for this project
	|
	|- Jenior_Metatransciptomics_eLife_2017.Rmd		# executable Rmarkdown for this study
	|- Jenior_Metatransciptomics_eLife_2017.md		# Markdown (GitHub) version of the *.Rmd file
	|- Jenior_Metatransciptomics_eLife_2017.docx		# Word document of study manuscript
	|- manuscript_format.docx				# Word document associated with Rmd for proper formatting
	|- references.bib					# BibTeX formatted references
	|- eLife.csl						# csl file to format references
	|
	|- doc/		# Metagenome & metatranscriptome bioinformatic processing pipelines
	|
	|- data/          # raw and primary data
	| |- references/  # reference files used in analysis
	| |- raw/         # raw data, unaltered
	| +- process/     # cleaned data, unaltered once created
	|
	|- code/	# all programmatic code (python, R, pbs, mothur_batch)
	|
	|- results/			# all output from workflows and analyses
	| |- figures/		# manuscript figures
	| +- supplement/	# supplementary tables and figures

