#!/usr/bin/python
'''USAGE: get_common_taxa.py fasta_file
Finds most common KEGG taxon label from annotated fasta
'''
import sys
from operator import itemgetter

with open(sys.argv[1],'r') as input_fasta:

        taxa_dictionary = {}
        total_genes = 0

        for line in input_fasta:

                if str(line)[0] == '>':
                        total_genes += 1

                        taxa = str(line).strip('>').split(':')[0]

                        if not taxa in taxa_dictionary.keys():
                                taxa_dictionary[taxa] = 1
                        else:
                                taxa_dictionary[taxa] = taxa_dictionary[taxa] + 1

                else:
                        continue

taxa_list = []
for taxa in taxa_dictionary.keys():
        percentage = round(((float(taxa_dictionary[taxa]) / float(total_genes)) * 100), 1)
        taxa_list.append([taxa, taxa_dictionary[taxa], percentage])

taxa_list = sorted(taxa_list, key=itemgetter(1), reverse=True)[0:5]

print('Fasta file: ' + str(sys.argv[1]))
print('Total genes: ' + str(total_genes))
print('1st most frequent:\t' + taxa_list[0][0] + '\t' + str(taxa_list[0][1]) + '\t' + str(taxa_list[0][2]) + '%')
print('2nd most frequent:\t' + taxa_list[1][0] + '\t' + str(taxa_list[1][1]) + '\t' + str(taxa_list[1][2]) + '%')
print('3rd most frequent:\t' + taxa_list[2][0] + '\t' + str(taxa_list[2][1]) + '\t' + str(taxa_list[2][2]) + '%')
print('4th most frequent:\t' + taxa_list[3][0] + '\t' + str(taxa_list[3][1]) + '\t' + str(taxa_list[3][2]) + '%')
print('5th most frequent:\t' + taxa_list[4][0] + '\t' + str(taxa_list[4][1]) + '\t' + str(taxa_list[4][2]) + '%\n')
