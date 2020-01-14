#!/usr/bin/env python3      
'''
Iz FASTA datoteke poravnanja genierira originalne slijedove.
Primjer pokretanja:

python3 genes_from_alignment.py HIV1_REF_2010_genome_DNA.fasta -no_gaps/-gaps

'''
import os.path
import sys
__author__ = "Stjepan Jakovac"
__credits__ = ["Jurij Kos", "Frane Polić"]
__license__ = "MIT"


ALGMNT_START_MARKER = '>'
algmnt_fn = sys.argv[1]

remove_gaps = True if sys.argv[2] == '-no_gaps' else False
base_dir = "../data/"
name = "g{}.txt" if remove_gaps else "a{}.txt"

with open(algmnt_fn, 'r') as algmnt:
	
	out_file = None
	count = -1
	for line in algmnt.readlines():
		if line.startswith(ALGMNT_START_MARKER):
			count += 1
			if out_file is not None:
				out_file.close()
			out_file = open(os.path.expanduser(os.path.join(base_dir ,name.format(count))), 'w')
			continue
		
		line = line.strip()
		if remove_gaps:
			line = line.replace('-', '')
		
		out_file.write(line)
		
		

