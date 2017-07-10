#!/bin/python
# slice_fasta.py raw.fasta 1000000000

import sys

output_fasta = str(sys.argv[1]).rstrip('fasta') + 'pick.fasta'
output_fasta = open(output_fasta, 'w')
max_seqs = int(sys.argv[2])
current_seqs = 0

with open(sys.argv[1], 'r') as input_fasta:

	for line in input_fasta:

		if line == '\n':
			continue

		elif line[0] == '>' and current_seqs < max_seqs:
			output_fasta.write(line)
			current_seqs += 1
			continue

		elif line[0] != '>' and current_seqs <= max_seqs:
			output_fasta.write(line)
			continue

		else:
			break


output_fasta.close()
