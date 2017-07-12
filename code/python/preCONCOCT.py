#!/bin/python
# preCONCOCT.py raw.fasta formatted.fasta
# Formats metagenomic contigs for input to CONCOCT according to the author's recommendations

import sys

output_fasta = open(sys.argv[2], 'w')

with open(sys.argv[1], 'r') as input_fasta:

	for line in input_fasta:

		if line == '\n':
			continue

		elif line[0] == '>':
			current_name = line.replace(' ', '_')
			continue

		else:
			if len(line) < 1000:
				continue
			elif len(line) <= 10000:
				output_fasta.write(current_name)
				output_fasta.write(line)
				output_fasta.write('\n')
				continue
			
			else:
				line = line.strip()
				line = [line[x:x+10000] for x in range(0, len(line), 10000)]
				current_name = current_name.strip() + '_'

				for y in range(0, len(line)):
					if len(line[y]) < 1000:
						continue
					else:
						output_fasta.write(current_name + str(y) + '\n')
						output_fasta.write(line[y] + '\n\n')


output_fasta.close()
