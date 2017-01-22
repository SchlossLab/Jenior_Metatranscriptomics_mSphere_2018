#!/usr/bin/env python
'''USAGE: python python competition.py interaction.files --p n.s.
Calculates putative community-level competition between species from aggregated bigSMALL analysis
'''
# work on acronym...


# 
# 
# 
# 
# 
# 
# 
# 
# 


# Initialize all modules, functions, and compound dictionary
import sys
import math
import os
import pickle
import argparse
from operator import itemgetter


#---------------------------------------------------------------------------------------#


# Find the directory where interaction.py is being run from
script_path = str(os.path.dirname(os.path.realpath(__file__)))

# Create a string with the name of the initial directory you start running the script from
starting_directory = str(os.getcwd())


#---------------------------------------------------------------------------------------#


# Define functions

# Function to read in all species cominations from interaction file
def read_files(file_names):

	names = []
	omit = []

	for line in file_names:

		species = line.strip()
		if species in omit:
			continue

		if os.path.exists(species) == False:
			print('WARNING: ' + species + ' does not exist. Omitting combinations.')
			omit.append(species)
			continue

		names.append(species)
	
	return names


# Function to get the organism name from the bigSMALL parameters file
def read_parameters(parameters):

	for line in parameters:

		line = line.split()

		if line[0] == 'KO expression file:':
			file_name = str(line[1])

		elif line[0] == 'Graph name:':
			graph_name = str(line[1])

	if graph_name != 'organism':
		name = graph_name
	else:
		name = file_name

	return name


# Create dictionary for full compound names
def read_names(importance_scores):

	name_dictionary = {}

	for line in compound_scores:
		
		line = line.split()
		if line[0] == 'Compound_code':
			continue
		
		compound_name = str(line[1])
		name_dictionary[compound_code] = compound_name
		
	return name_dictionary


# Reads importance files, applying p-value filter, generates a dictionary for compound names, and returns a list of compound scores
def read_scores(importance_scores, p_cutoff):

	name_dictionary = {}
	score_list = []

	for line in compound_scores:
		
		line = line.split()
		if line[0] == 'Compound_code':
			continue
		
		compound_code = str(line[0])

		score = float(line[2])
		if score < 0:
			continue

		p_value = str(line[3]).strip('<')
		if p_value == 'n.s.': 
			p_value = 1.0
		else:
			p_value = float(p_value)
		if p_value > p_cutoff:
			continue
		else:
			score_list.append([compound_code, score])
		
	return score_list


# Sort, Rank, and Normalize importance scores
def convert_scores(scores, factor):

	sorted_scores = sorted(scores, key=itemgetter(1), reverse=True)
	important_total = len(sorted_scores)
	new_score_dictionary = {}

	for index in sorted_scores:
		ranked_score = important_total
		important_total = important_total - 1
		new_score = ranked_score * factor
		new_score_dictionary[index[0]] = new_score

	return new_score_dictionary


# Function for calculating edges of metabolic competition
def competition(score_dict_1, score_dict_2):
	
	all_compounds = list(set(score_dict_1.keys() + score_dict_2.keys()))
	
	competition_dictionary = {}

	for index in all_compounds:

		try:
			score_1 = float(score_dict_1[index])
		except keyError:
			continue
		try:
			score_2 = float(score_dict_2[index])
		except keyError:
			continue

		score_1 = score_dict_1[index][1]
		score_2 = score_dict_2[index][1]

		competition_score = (score_1 / score_2) + (score_2 / score_1)
		competition_score = 100 - competition_score

		competition_dictionary[index] = competition_score

	return competition_dictionary


# Calculate cumulative impostance of each compound across the groups of organisms tested
def community_competition(community, normalized):

	all_compounds = list(set(community.keys() + normalized.keys()))

	for compound in all_compounds:
		try:
			community[compound] = community[compound] + normalized[compound]
		except keyError:
			community[compound] = normalized[compound]

	return community


# Function to write data to output file
def write_output(header, competiton, names, file_name):

	with open(file_name, 'w') as outfile: 
		
		outfile.write(header)
			
		for compound_code in competiton.keys():
			competition_score = competiton[compound_code]
			compound_name = names[compound_code]
			entry = '\t'.join([compound_code, compound_name, str(competition_score) + '\n'])
			outfile.write(entry)	
	

#---------------------------------------------------------------------------------------#


# Set up arguments
parser = argparse.ArgumentParser(description='Calculate metabolic pair-wise interactions of species from the output of bigSMALL.')
parser.add_argument('input_file')
parser.add_argument('--p', default='n.s.', help='Minimum p-value for metabolites to be considered in calculations')

args = parser.parse_args()
interactions = str(args.input_file)
p_value = float(args.p_value)

if interactions ==  'none': sys.exit('WARNING: Missing input file, quitting')
if p_value != 'n.s.' and p_value < 0.0: sys.exit('WARNING: invalid p-value cutoff, quitting')


#---------------------------------------------------------------------------------------#


# Create a new directory for community interaction files
directory = str(os.getcwd()) + '/community.files'
if not os.path.exists(directory):	
	os.makedirs(directory)
os.chdir(directory)

#---------------------------------------------------------------------------------------#


# Retrieve and read in the necessary files
print('\nReading in interaction file.\n')
species_list = read_files(interactions)

species = {}
name_dictionary = {}
for index in species_list:
	os.chdir(species_1)

	importances = open('importances.tsv','r')
	organism_name = read_parameters(open('parameters.txt','r'))
	scores = read_scores(importances, p_value)
	species = species + scores
	compound_names = read_names(importances)
	name_dictionary = name_dictionary + compound_names

	os.chdir(starting_directory)


current = 0
for x in range(0, len(species)):

	y = x + 1
	if y == len(species): continue

	for y in range(0, len(species))
		current += 1
		print('Calculating interaction ' + str(current))

		# Retrieve individual species information
		scores_1 = species[x]
		scores_2 = species[y]

		# Normalize scores to the same length ##### NEED TO FIX TO ONLY NORMALIZE ONE
		if len(scores_1.keys()) > len(scores_2.keys()):
			reduction_factor = len(scores_2.keys()) / len(scores_1.keys())
			norm_scores_1 = convert_scores(scores_1, reduction_factor)
			norm_scores_2 = scores_2
		elif len(scores_1.keys()) < len(scores_2.keys())
			reduction_factor = len(scores_1.keys()) / len(scores_2.keys())
			norm_scores_2 = convert_scores(scores_2, reduction_factor)
			norm_scores_1 = scores_1
		else:
			norm_scores_1 = scores_1
			norm_scores_2 = scores_2

		# Calculate putative metabolic interactions and parse the output
		crosstalk = competiton(norm_scores_1, norm_scores_2)
		
		# Write output tables and summary to files
		head = 'Compound_code\tCompound_name\tCompetition_score\n'
		file_name = name_1 + '.vs.' + name_2 + '.competition.tsv'
		write_output(head, crosstalk, file_name)

		# Add interaction to community summary


		###### need and if statement deciding if an organism has already been added

		if current == 1:
			community_scores = crosstalk
		else:
			community_scores = community_competition(community_scores, norm_scores_1)
			community_scores = community_competition(community_scores, norm_scores_2)


# Write cumulative scores to a file
out_file_name = 'Community.' + str(p_value) + '.tsv'
head = 'Compound_code\tCompound_name\tCommunity_Metabolite_Score\tp_value\n'
write_output(head, crosstalk, out_file_name)


# Navigate back to starting directory
os.chdir(starting_directory)	

print('Done.\n')
