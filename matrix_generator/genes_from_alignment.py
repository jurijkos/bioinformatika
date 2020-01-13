#!/usr/bin/env python3
import sys

'''
Iz FASTA datoteke poravnanja genierira originalne slijedove.
Primjer pokretanja:

python3 genes_from_alignment.py HIV1_ALL_2010_genome_DNA.fasta gH1.txt gH2.txt -no_gaps/-gaps

'''
ALGMNT_START_MARKER = '>'

algmnt_fn = sys.argv[1]
g_1_fn = sys.argv[2]
g_2_fn = sys.argv[3]

remove_gaps = True if sys.argv[4] == '-no_gaps' else False

with open(algmnt_fn, 'r') as algmnt, \
		open(g_1_fn, 'w') as g_1, \
		open(g_2_fn, 'w') as g_2:
	
	curr = None
	for line in algmnt.readlines():
		if line.startswith(ALGMNT_START_MARKER):
			curr = g_1 if curr is None else g_2
			continue
		
		line = line.strip()
		if remove_gaps:
			line = line.replace('-', '')
		
		curr.write(line)
		
		

