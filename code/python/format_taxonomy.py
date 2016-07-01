#!/usr/bin/env python
'''USAGE: 
python format_taxonomy.py final.cons.taxonomy --label 1-3-5
'''
# This script reformats the taxonomy files output by mothur to be more readable (helps with plot labels)

import sys
import re
import argparse

parser = argparse.ArgumentParser(description='Reformats taxonomy files to indicate desired clasification level neatly.')
parser.add_argument('infile')
parser.add_argument('--label', default='0', help='5=phylum, 4=class, 3=order, 2=family, 1=genus, 0=last classified')

args = parser.parse_args()
labels = str(args.label).split('-')

tax_dict={}
tax_dict['5'] = 'phylum'
tax_dict['4'] = 'class'
tax_dict['3'] = 'order'
tax_dict['2'] = 'family'
tax_dict['1'] = 'genus'
tax_dict['0'] = 'last'

print('\nOutput file names:')

for index in labels:

	try:
		suffix = '.' + tax_dict[str(index)] + '.format.taxonomy'
	except KeyError:
		sys.exit('Error:  Invalid taxonomy label provided.')
	outfile_name = str(args.infile).replace('.taxonomy', suffix)
	outfile = open(outfile_name,'w')

	if index == '0':

		with open(args.infile,'r') as infile:

			header = 0
			taxon_lst = []

			for line in infile:
		
				if header == 0:
					outfile.write(line)
					header += 1
					continue
		
				line = line.split()
				
				otu = line[0].lstrip('Otu')
				otu = str(otu).lstrip('0')
				otu = ' (OTU' + str(otu) + ')'
		
				temp_taxon = line[2].split(';')
		
				taxon = [re.sub('\(.*?\)', '', x) for x in temp_taxon]
		
				if taxon[1] != 'unclassified':
			
					if taxon[2] != 'unclassified':
					
						if taxon[3] != 'unclassified':
					
							if taxon[4] != 'unclassified':
						
								if taxon[5] != 'unclassified':
									best_resolution = taxon[5] + otu
								else:
									best_resolution = taxon[4] + otu
							else:
								best_resolution = taxon[3] + otu
						else:
							best_resolution = taxon[2] + otu
					else:
						best_resolution = taxon[1] + otu
			
					final_taxon = best_resolution
			
				else:
					final_taxon = 'Bacteria_unclassified' + otu


				final_taxon = final_taxon.replace('/', '_')
				entry = '\t'.join([line[0], line[1], final_taxon + '\n'])
				outfile.write(entry)
				
	else:
	
		with open(args.infile,'r') as infile:

			header = 0
			taxon_lst = []

			for line in infile:
		
				if header == 0:
					outfile.write(line)
					header += 1
					continue
		
				line = line.split()
		
				temp_taxon = line[2].split(';')
		
				taxon = [re.sub('\(.*?\)', '', x) for x in temp_taxon]
		
				if index == '5':
					final_taxon = taxon[1]
					
				elif index == '4':
					final_taxon = '_'.join([taxon[1], taxon[2]])
					
				elif index == '3':
					final_taxon = '_'.join([taxon[1], taxon[3]])
					
				elif index == '2':
					final_taxon = '_'.join([taxon[1], taxon[4]])
					
				elif index == '1':
					final_taxon = '_'.join([taxon[1], taxon[5]])
					
				if taxon[1] == 'unclassified':
					final_taxon = 'Bacteria_unclassified'
	
				if final_taxon in taxon_lst:
					taxon_count = taxon_lst.count(final_taxon) + 1
					taxon_lst.append(final_taxon)
					final_taxon = ''.join([final_taxon, '_', str(taxon_count)])
				else:
					taxon_lst.append(final_taxon)
				
				final_taxon = final_taxon.replace('/', '')
				entry = '\t'.join([line[0], line[1], final_taxon + '\n'])
				outfile.write(entry)
	
	outfile.close()
	print(outfile_name)


print('\n')

