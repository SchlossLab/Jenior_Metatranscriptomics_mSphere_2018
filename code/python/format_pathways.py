#!/usr/bin/env python
'''USAGE: 
python format_taxonomy.py input_pathways.tsv output_pathways.tsv
'''
# Pools pathway annotations for metagenome-normalized metatranscriptomic mappings of mock and c diff infected

import sys
import re


with open(sys.argv[1],'r') as infile:

	pathway_630_dictionary = {}
	pathway_mock_dictionary = {}
	header = 0

	for line in infile:

		if header == 0:
			header += 1
			continue

		line = line.split()

		pathway = line[-2]
		pathway = pathway.split('|')
		pathways = []
		for index1 in pathway:
			index1 = index1.split(';')[-1]
			pathways.append(index1)

		for index2 in pathways:
			if not index2 in pathway_630_dictionary.keys():
				pathway_630_dictionary[index2] = float(line[0])
			else:
				pathway_630_dictionary[index2] += float(line[0])

			if not index2 in pathway_mock_dictionary.keys():
				pathway_mock_dictionary[index2] = float(line[1])
			else:
				pathway_mock_dictionary[index2] += float(line[1])


with open(sys.argv[2],'w') as outfile:

	outfile.write('norm_reads_630\tnorm_reads_mock\tpathway\n')

	all_keys = set(pathway_630_dictionary.keys() + pathway_mock_dictionary.keys())
	omit_paths = set(['Nervous_system', 'Cardiovascular_diseases', 'Endocrine_and_metabolic_diseases', 'Immune_diseases', 'Environmental_adaptation',
		'Substance_dependence', 'Excretory_system', 'Circulatory_system', 'Infectious_diseases:_Parasitic', 'Endocrine_system', 'Immune_system', 'Neurodegenerative_diseases',
		'Infectious_diseases:_Viral', 'Infectious_diseases:_Bacterial', 'Cancers:_Specific_types', 'Digestive_system', 'Cancers:_Overview','Development', 'Sensory_system'])
	all_keys = list(all_keys - omit_paths)

	for index3 in all_keys:

		try:
			reads_630 = pathway_630_dictionary[index3]
		except keyError:
			reads_630 = 0
		try:
			reads_mock = pathway_mock_dictionary[index3]
		except keyError:
			reads_mock = 0

		outline = '\t'.join([str(reads_630), str(reads_mock), index3 + '\n'])
		outfile.write(outline)