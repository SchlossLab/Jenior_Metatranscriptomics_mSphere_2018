#!/usr/bin/env python 

# USAGE: QC_mapping.py input_file
# Finds the percentage of reads that mapped to the reference

import sys

with open(sys.argv[1], 'r') as infile:
	
	mapped = 0 

	for line in infile:
		
		line = line.split()
		
		if line[0] == '*':
			unmapped = float(line[3])
			continue
		else:
			mapped += float(line[3])

total = mapped + unmapped
percent = (mapped / total) * 100
percent = round(percent, 2)

print('\nTotal reads: ' + str(int(total)))	
print('Mapped reads: ' + str(int(mapped)))	
print('Unmapped reads: ' + str(int(unmapped)))	
print('\nPercent mapped of total reads: ' + str(percent) + '%\n')

