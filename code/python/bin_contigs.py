#!/usr/bin/python
'''USAGE: bin_contigs.py concoct_clusters contig_fasta
Bins contigs based on output of CONCOCT clustering into separate fasta files
'''
import sys
import os


# Create a new output directory
output_directory = str(os.getcwd()) + '/binned_contigs/'
if not os.path.exists(output_directory):
	os.makedirs(output_directory)


# Create dictionary of contigs and corresponding clusters
with open(sys.argv[2], 'r') as clusters:
	cluster_dict = {}
	unique_clusters = set()
	for line in clusters:
		contig = line.strip().split(',')[0]
		clusters = set(line.strip().split(',')[1:])
		unique_clusters |= clusters
		clusters = [int(x) for x in list(clusters)]
		cluster_dict[contig] = clusters


# Read through contig file and bin them to correct cluster fasta 
contig_fasta = open(sys.argv[1],'r')
total = 0
count = 0
omit = 0

for line in contig_fasta:

	if line == '\n': continue

	elif str(line)[0] == '>':
		total += 1
		count += 1
		if count == 1000:
			print(str(total))
			count = 0
		contig_name = str(line).lstrip('>').split()[0]
		try:
			current_clusters = cluster_dict[contig_name]
		except KeyError:
			omit += 1
			output_fasta_name = output_directory + 'unclustered.contigs.fasta'
			with open(output_fasta_name, 'a+') as current_output:
				current_output.write('>' + contig_name + '\n')
			continue

		for cluster in current_clusters:
			output_fasta_name = output_directory + 'cluster_' + str(cluster) + '.contigs.fasta'
			with open(output_fasta_name, 'a+') as current_output:
				current_output.write('>' + contig_name + '\n')
		continue

	else:
		for cluster in current_clusters:
			output_fasta_name = output_directory + 'cluster_' + str(cluster) + '.contigs.fasta'
			with open(output_fasta_name, 'a+') as current_output:
				current_output.write(line + '\n')

contig_fasta.close()


# Give the user some feedback in a log file
log_name = output_directory + 'log.txt'
with open(log_name, 'w') as log_file:
	log_file.write('Contig clusters: ' + str(len(list(unique_clusters))) + '\n')
	log_file.write('Successfully clustered contigs: ' + str(total - omit) + '\n')
	log_file.write('Unclustered contigs: ' + str(omit))

print('Done')