#!/usr/bin/python
'''USAGE: rmdup_SQ.py input.sam output.same
removes duplicate @SQ entries from SAM files
'''
import sys
import os

# Read SAM file header looking for duplicates
with open(sys.argv[1], 'r') as input_sam:
	output_sam = open(sys.argv[2], 'w')
	headers = set()

	for line in input_sam:

		if not line.split()[0] == '@SQ':
			output_sam.write(line)
			continue
		else:
			header = line.split()[1]
			if not header in headers:
				output_sam.write(line)
				headers.add(header)
				continue
			else:
				print('Duplicated header: ' + header)
				continue

output_sam.close()
