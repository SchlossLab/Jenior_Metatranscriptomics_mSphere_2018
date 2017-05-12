#!/usr/bin/python
'''USAGE: get_ko.py idxstats_file
Retreives KOs and their mapped reads from my idxstats files
'''
import sys

with open(sys.argv[1],'r') as idxstats:

        ko_dictionary = {}
        for line in idxstats:
                line = line.strip().split('|')
                if len(line) < 4: continue

                reads = float(line[-1].split()[1])

                if line[2][0] == 'K':
                        if not line[2] in ko_dictionary.keys():
                                ko_dictionary[line[2]] = reads
                        else:
                                ko_dictionary[line[2]] = ko_dictionary[line[2]] + reads
                        continue
                elif line[3][0] == 'K':
                        if not line[3] in ko_dictionary.keys():
                                ko_dictionary[line[3]] = reads
                        else:
                                ko_dictionary[line[3]] = ko_dictionary[line[3]] + reads
                        continue
                else:
                        continue

output_name = str(sys.argv[1]).rstrip('txt') + 'ko.txt'
with open(output_name, 'w') as output_file:
        total_kos = 0
        total_reads = 0
        for index in ko_dictionary.keys():
                total_kos += 1
                total_reads += ko_dictionary[index]
                entry = '\t'.join([index, str(ko_dictionary[index])+'\n'])
                output_file.write(entry)

print('Total KOs: ' + str(total_kos))
print('Total normalized reads: ' + str(total_reads))