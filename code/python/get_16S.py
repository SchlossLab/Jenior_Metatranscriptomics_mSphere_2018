#!/usr/bin/python
'''USAGE: get_16S.py fasta_file
Find genes annotated as 16S rRNA genes in a .fasta and writes them to a new file
'''
import sys

with open(sys.argv[1],'r') as input_fasta:

        rRNA_16S_genes = 0
        include = 0
        outfile_name = str(sys.argv[1]).rstrip('fastn') + 'rRNA.fasta'
        output_fasta = open(outfile_name, 'w')

        for line in input_fasta:

                if str(line)[0] == '>':

                        include = 0
                        name = str(line).strip('>').strip()

                        if 'rps' in name or 'rrs' in name or 'rpl' in name or 'rpm' in name:
                                rRNA_16S_genes += 1
                                include = 1
                                output_fasta.write('>' + name + '\n')
                                continue
                        else:
                                continue

                elif include == 1:
                        output_fasta.write(line)
                        continue

output_fasta.close()

