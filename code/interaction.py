#!/usr/bin/env python
'''USAGE: python python interaction.py interaction.files --p n.s.
Calculates putative community-level and pair-wise metabolic interations between species from aggregated bigSMALL analysis
'''


# work on acronym...
#comparison of bigsmall out for ecological interaction inference....



# Initialize all modules, functions, and compound dictionary
import sys
import numpy
import os
import argparse

#---------------------------------------------------------------------------------------#

# Find the directory where interaction.py is being run from
script_path = str(os.path.dirname(os.path.realpath(__file__)))

# Create a string with the name of the initial directory you start running the script from
starting_directory = str(os.getcwd())

#---------------------------------------------------------------------------------------#

# Set up arguments
parser = argparse.ArgumentParser(description='Calculate metabolic pair-wise interactions of species from the output of bigSMALL.')
parser.add_argument('input_file')
parser.add_argument('--p_max', default='n.s.', help='Minimum p-value for metabolites to be considered in calculations')

args = parser.parse_args()
interactions = args.input_file
p_value = args.p_max

if os.stat(interactions).st_size == 0 : sys.exit('WARNING: Input file empty, quitting')
if p_value != 'n.s.' and p_value < 0.0: sys.exit('WARNING: Invalid p-value cutoff, quitting')
print('\n')

#---------------------------------------------------------------------------------------#

# Define functions

# Function to read in all species cominations from interaction file
def read_files(files):

	with open(files, 'r'):

		interactions = []
		community = []

		for line in files:
			line = line = split()

			species1 = line[0]
			species2 = line[1]

			if os.path.exists(species1) == False:
				print('WARNING: ' + species1 + ' does not exist. Omitting combination.')
				continue
			elif os.path.exists(species2) == False:
				print('WARNING: ' + species2 + ' does not exist. Omitting combination.')
				continue

			interactions.append([species1, species2])

	return list(set(interactions))


# Reads importance files, applying p-value filter, normalizes score to reads, and generates a dictionary for compound names and compound scores
def read_scores(importance_scores, p_cutoff):

	score_dictionary = {}
	name_dictionary = {}

	for line in importance_scores:
		
		line = line.split()
		if line[0] == 'Compound_code': continue
		
		compound_code = str(line[0])
		compound_name = str(line[1])
		
		score = float(line[2])
		if score < 0.0:
			score = -(2 ^ -score)
		else:
			score = 2 ^ score
		
		p_value = str(line[3]).lstrip('<')
		if p_value == 'n.s.': p_value = 1
		p_value = float(p_value)
		if p_value > p_cutoff: continue

		score_dictionary[compound_code] = score
		name_dictionary[compound_code] = compound_name	

	final_score_dictionary = {}
	score_sum = sum([abs(x) for x in score_dictionary.values()])
	for index in score_dictionary.keys()
		final_score = score_dictionary[index] / score_sum
		final_score_dictionary[index] = [name_dictionary[index], final_score]

	return final_score_dictionary


# Function for calculating edges of metabolic competition
def single_interaction(score_dict_1, score_dict_2):
	
	all_compounds = list(set(score_dict_1.keys() + score_dict_2.keys()))
	
	interaction_dictionary = {}
	for index in all_compounds:

		# Consider only those compounds in the network of both organisms
		try:
			score_1 = float(score_dict_1[index][1])
			name = str(score_dict_1[index][0])
		except keyError:
			continue
		try:
			score_2 = float(score_dict_2[index][1])
		except keyError:
			continue

		# Determine type of interaction
		if score_1 > 0.0 and score_2 < 0.0:
			relationship = 'synergy'
		elif score_1 < 0.0 and score_2 > 0.0:
			relationship = 'synergy'
		elif score_1 > 0.0 and score_2 > 0.0:
			relationship = 'competition'
		else:
			relationship = 'none'

		# Report size of difference (arbitrary)
		magnitude = abs(score_1 - score_2)
		if relationship == 'synergy':
			if magnitude >= 100:
				strength = 'strong'
			elif magnitude >= 50:
				strength = 'moderate'
			else:
				strength = 'weak'
		elif relationship == 'competition':
			if magnitude <= 50:
				strength = 'strong'
			elif magnitude <= 100:
				strength = 'moderate'
			else:
				strength = 'weak'
		else:
			strength = 'none'

		interaction_dictionary[index] = [name, score_1, score_2, magnitude, relationship, strength]

	return interaction_dictionary


# Calculate cumulative importance of each compound across the groups of models tested
def community_interaction(community_dict, member_dict):

	for compound in member_dict.keys():
		name = member_dict[compound][0]
		score = member_dict[compound][1]
	
		try:
			cumulative_score = community_dict[compound][2] + score
			community_dict[compound] = [name, cumulative_score]
		except keyError:
			community_dict[compound] = [name, score]

	return community_dict


# Calculates the percentile for the given type of score
def calc_percentile(dictionary, type_index):

	all_values = []
	for index in dictionary.keys():
		all_values.append(dictionary[index][type_index])

	per_90 = [numpy.percentile(current_list, 90), numpy.percentile(current_list, 10)]
	per_80 = [numpy.percentile(current_list, 80), numpy.percentile(current_list, 20)]
	per_70 = [numpy.percentile(current_list, 70), numpy.percentile(current_list, 30)]
	per_60 = [numpy.percentile(current_list, 60), numpy.percentile(current_list, 40)]

	for index in dictionary.keys():
		if dictionary[index][type_index] < per_90[0] or dictionary[index][type_index] > per_90[1]:
			dictionary[index] = dictionary[index].append('90th')
			continue
		elif dictionary[index][type_index] < per_80[0] or dictionary[index][type_index] > per_80[1]:
			dictionary[index] = dictionary[index].append('80th')
			continue
		elif dictionary[index][type_index] < per_70[0] or dictionary[index][type_index] > per_70[1]:
			dictionary[index] = dictionary[index].append('70th')
			continue
		elif dictionary[index][type_index] < per_60[0] or dictionary[index][type_index] > per_60[1]:
			dictionary[index] = dictionary[index].append('60th')
			continue
		else:
			dictionary[index] = dictionary[index].append('none')

	return dictionary


# Function to write data to output file
def write_output(header, dictionary, file_name, type_output):

	with open(file_name, 'w') as outfile:
		outfile.write(header)
		
		if type_output == 'single':

			for index in dictionary.keys():
			
				name = dictionary[index][0]
				score_1 = dictionary[index][1]
				score_2 = dictionary[index][2]
				magnitude = dictionary[index][3]
				strength = dictionary[index][4]
				percentile = dictionary[index][5]

				# Transform scores back to log2
				if score_1 == 0.0:
					score_1 = 0.0
				elif score_1 < 0.0:
					score_1 = numpy.log2(abs(score_1)) * -1
				else:
					score_1 = numpy.log2((score_1))
				if score_2 == 0.0:
					score_2 = 0.0
				elif score_2 < 0.0:
					score_2 = numpy.log2(abs(score_2)) * -1
				else:
					score_2 = numpy.log2((score_2))

				entry = '\t'.join([index, name, str(score_1), str(score_2), str(magnitude), strength, percentile]) + '\n'
				outfile.write(entry)

		elif type_output == 'community':

			for index in dictionary.keys():
			
				name = dictionary[index][0]
				score = dictionary[index][1]
				magnitude = dictionary[index][2]
				percentile = dictionary[index][3]

				# Transform scores back to log2
				if score == 0.0:
					score = 0.0
				elif score < 0.0:
					score = numpy.log2(abs(score)) * -1
				else:
					score = numpy.log2((score))

				entry = '\t'.join([index, name, str(score), percentile]) + '\n'
				outfile.write(entry)


#---------------------------------------------------------------------------------------#


# Create a new directory for community interaction files
directory = str(os.getcwd()) + '/community.files'
if not os.path.exists(directory):	
	os.makedirs(directory)
os.chdir(directory)

#---------------------------------------------------------------------------------------#

# Worflow

# Retrieve and read in the necessary files
interactions_list = read_files(interactions)
community_dictionary = {}
community = []
current = 0
for index in interactions_list:

	os.chdir(index[0])
	scores_1 = read_scores(open('importances.tsv','r'), p_value)
	os.chdir(starting_directory)
	os.chdir(index[1])
	scores_2 = read_scores(open('importances.tsv','r'), p_value)
	os.chdir(starting_directory)

	current += 1
	print('Calculating interaction ' + str(current) + ' of ' + str(len(interactions_list)) + '.')
	interaction = single_interaction(scores_1, scores_2)
	interaction = calc_percentile(interaction, 3)

	if not str(index[0]) in community:
		community_dictionary = community_interaction(community_dictionary, scores_1)
		community.append(str(index[0]))
	if not str(index[1]) in community:
		community_dictionary = community_interaction(community_dictionary, scores_2)
		community.append(str(index[1]))

	header = 'Compound_code\tCompound_name\tScore_1\tScore_2\tMagnitude\tPercentile\n'
	file_name = str(index[0]) + '.and.' + str(index[1]) + '.interaction.txt'
	write_output(header, interaction, file_name, 'single')

# Write cumulative scores to a file
community = calc_percentile(community, 1)
header = 'Compound_code\tCompound_name\tCommunity_Metabolite_Score\tPercentile\n'
file_name = 'community_importance.tsv'
write_output(header, community, file_name, 'community')
print('Done\n')