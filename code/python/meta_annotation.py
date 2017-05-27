#!/usr/bin/python
'''USAGE: meta_annotation.py input_annotation output_annotation
Reformats gene annotations used in BLAST database
'''
import sys

with open(sys.argv[1],'r') as raw_annotation:
	output_annotation = open(sys.argv[2], 'w')
	output_annotation.write('entry\torganism\tgene\tdescription\tko\n')

        for line in raw_annotation:
        	if line == '\n': continue
        	line = line.strip()
        	entry = line.split('|')

        	if entry[-1].split('_')[0] == 'unknown':
        		gene = entry[0]
                        org = 'NA'
        		description = entry[-1]
        		ko = 'NA\n'

        	elif len(entry) == 4:
        		gene = entry[0]
                        org = gene.split(':')[0]
        		description = entry[1]
        		ko = entry[2] + '\n'

        	else:
        		gene = entry[0]
                        org = gene.split(':')[0]
        		description = entry[2]
        		ko = entry[3] + '\n'
        	
        	output = '\t'.join([line, org, gene, description, ko])
        	output_annotation.write(output)

output_annotation.close()
