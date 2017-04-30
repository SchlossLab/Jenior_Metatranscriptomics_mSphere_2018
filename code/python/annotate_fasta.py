#!/usr/bin/python
'''USAGE: annotate_fasta.py BLASToutput Fasta
This script annotates a fasta file with sequence names determined from a BLAST of that file.
'''
import sys


# Write dictionary from BLAST output
with open(sys.argv[1],'r') as blast_output:
        blast_dictionary = {}
        gene_count = 0
        for line in blast_output:
                gene_count += 1
                current = str(gene_count).zfill(9)
                line = line.split()
                query = line[0]
                target = line[1]
                evalue = line[10]
                entry = target + '|evalue_' + str(evalue) + '|' + current
                blast_dictionary[query] = entry


# Parse the fasta file and print translated names and unchanged sequence data
with open(sys.argv[2],'r') as input_fasta:
        seq_count = 0
        annotated = 0
        key_errors = 0
        hypothetical = 0
        unknown = 0
        uncharacterized = 0
        outfile_name = str(sys.argv[2]).rstrip('fastn') + 'annotated.fasta'
        output_fasta = open(outfile_name, 'w')

        for line in input_fasta:

                if str(line)[0] == '>':
                        seq_count += 1
                        entry = str(line).strip('>').strip()
                        entry = entry.rstrip('\n')

                        try:
                                blast_hit = blast_dictionary[entry]
                                annotated += 1
                                if 'hypothetical' in blast_hit: hypothetical += 1
                                if 'unknown' in blast_hit: unknown += 1
                                if 'uncharacterized' in blast_hit: uncharacterized += 1
                        except KeyError:
                                key_errors += 1
                                blast_hit = entry + '|unknown:' + str(key_errors)

                        final_entry = '>' + blast_hit + '\n'
                        output_fasta.write(final_entry)
                        continue

                else:
                        output_fasta.write(line + '\n')

output_fasta.close()


# Write some stats to a file about the success of the annotation
logfile_name = str(sys.argv[2]).rstrip('fastn') + 'annotation_logfile.txt'
with open(logfile_name,'w') as logfile:
        logfile.write(str(seq_count) + ' total sequences.\n')
        logfile.write(str(key_errors) + ' sequences could not be annotated.\n')
        logfile.write(str(seq_count - key_errors) + ' sequences recieved some annotation.\n')
        logfile.write(str(hypothetical) + ' sequences annotated as hypothetical.\n')
        logfile.write(str(unknown) + ' sequences annotated as unknown function.\n')
        logfile.write(str(uncharacterized) + ' sequences annotated as uncharacterized.\n')
        logfile.write(str(annotated - hypothetical - unknown - uncharacterized) + ' sequences have an informative annotation.\n')

