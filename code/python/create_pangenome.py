#!/bin/python

# USAGE:  create_pangenome.py  contigs.files


import sys


def read_assembly(contig_file, output_file):

	file_name = str(contig_file)
	contig_count = 1
	
	with open(contig_file, 'r') as contigs:
	
		for line in contigs:
		
			if line == '\n': continue
			
			elif line.startswith('>'):
				line = line.strip()
				line = line + '|' + file_name + ':' + str(contig_count) + '\n'
				output_file.write(line)
				contig_count += 1
				
			else:
				output_file.write(line) 

	
with open(sys.argv[1], 'r') as assemblies:

	assembly_list = []
	
	for line in assemblies:
		assembly_list.append(line.strip())
		
		
pangenome = open('all_contigs.fasta', 'w')
for index in assembly_list:
	read_assembly(index, pangenome)