#!/usr/bin/env python3
'''
Primjer pokretanja:
-sam gradi poravnanja od parova gena
python3 matrix_generator -g gene_pairs.txt 10

-poravnanja predana u datoteci
python3 matrix_generator -a algnmt_pairs.txt
'''
import sys
import numpy as np
from Bio import Align
__author__ = "Stjepan Jakovac"
__credits__ = ["Jurij Kos", "Frane Polić"]
__license__ = "MIT"


#emiss_i = {'-':0, 'A':1, 'B':2, 'C':3}  #mapira baze na indexe za matrice
emiss_i = {'-':0, 'A':1, 'C':2, 'G':3, 'T':4}
emiss_matrix = np.zeros((len(emiss_i),len(emiss_i)))
trans_matrix = np.zeros((5,5)) # 5 jer start, gap_x, gap_y, miss/match, end

def generate_empty_matrices():
	return np.zeros((len(emiss_i),len(emiss_i))), np.zeros((5,5))  #emiss, trans
	
def cleaned_fasta(file_name, chars_to_remove = ('x','N', 'Y')):
	'''
	Iz datoteke ucita genom, izbaci iz nj. s pojave znakov iz tupla chars_to_remove.
	Vrati string reprezentaciju ociscenog genoma.
	'''
	genome = None
	with open("../data/"+file_name, 'r') as source:
		genome = source.readline()
	replace_table = str.maketrans({c:None for c in chars_to_remove})
	
	genome = genome.translate(replace_table)
	return genome.strip()

def populate_from_pair(g_1, g_2, trans_m, emiss_m, N = 10):
	'''
	Popunjava trans i emis matricu iz najvise N optimalnih poravnanja koja se dobivaju
	od genoma g1 i g2. Tablicu popunjavaju pojavama(count) te se ne radi pretovrba u vjerojatnosu
	matricu
	'''
	aligner = Align.PairwiseAligner()
	alignments = aligner.align(g_1, g_2)
	
	for i,alignment in enumerate(alignments):
		#print(alignment)
		if N is not None and i>=N:
			break
		populate_from_aligment(alignment, trans_m, emiss_m)
	
	
def populate_from_aligment(alignment, trans_m, emiss_m, emiss_i = emiss_i):
	'''
	Iz jednog poravnanja popunjuje trans i emiss matricu
	'''
	s_algmnt = alignment.format()
	algmnt_len = len(s_algmnt)//3-1
	prev_state = 0 #start | gap_a/emit_y | gap_b/emit_x | m (match/miss) | end
	state = None
	for i in range(algmnt_len):
		c1 = s_algmnt[i]
		c2 = s_algmnt[i+(algmnt_len+1)*2]
		
			
		if c1 == '-':
			state = 1
		elif c2 == '-': 
			state = 2
		elif c2 != '\n': 
			state = 3
		trans_m[prev_state][state] += 1
		prev_state = state
		
		i1 = emiss_i[c1]; i2 = emiss_i[c2]
		emiss_m[i1][i2] += 1			
	trans_m[prev_state][4] += 1
	
def populate_from_raw_aligment(raw_1, raw_2, trans_m, emiss_m, emiss_i = emiss_i):
	'''
	Iz dva pojedinacna string poravnanja popunjava matrice transmisje i emisije
	'''
	prev_state = 0 #start | gap_a | gap_b | m (match/miss) | end
	state = None
	for c1, c2 in zip(raw_1, raw_2):
		if c1 == '-' and c2 == '-':
			continue
			
		if c1 == '-':
			state = 1
		elif c2 == '-': 
			state = 2
		else: 
			state = 3
		trans_m[prev_state][state] += 1
		prev_state = state
		
		i1 = emiss_i[c1]; i2 = emiss_i[c2]
		emiss_m[i1][i2] += 1			
	trans_m[prev_state][4] += 1
	
def proba_from_counts(trans_m, emiss_m):
	'''
	Matrice u kojima je biljezen broj pojava transformia u vjerojatnosne matrice.
	(in-place zasada)
	'''
	row_sums = trans_m.sum(axis=1)
	row_sums[4] = 1 #da zadnji ne bude ZeroDiv (zadnji su P(#nesto|end_state))
	trans_m /= row_sums[:,None]
	
	gap_a_sum = emiss_matrix[0,:].sum()
	gap_b_sum = emiss_matrix[:,0].sum()
	mm_sum = emiss_matrix.sum() - gap_a_sum - gap_b_sum
	emiss_matrix[0,:] /= gap_a_sum
	emiss_matrix[:,0] /= gap_b_sum
	emiss_matrix[1:,1:] /= mm_sum


mode = 0 if sys.argv[1] == '-g' else 1
pairs = []
config_file_name = sys.argv[2]
with open(config_file_name, 'r') as config_file:
	for line in config_file.readlines():
		pairs.append(line.split())
		
		
#iz datoteke zadane preko komandne linije ucitava parove iz kojih ce se graditi matrice
if mode == 0:
	N = int(sys.argv[3])
	#iterira po parovima i popunjava tablice
	for fn_1, fn_2 in pairs:
		g_1 = cleaned_fasta(fn_1)
		g_2 = cleaned_fasta(fn_2)
		populate_from_pair(g_1, g_2, trans_matrix, emiss_matrix, N = N)

#popunjava tablice iz gotovog poravnanja
else:
	for fn_1, fn_2 in pairs:
		g_1 = cleaned_fasta(fn_1)
		g_2 = cleaned_fasta(fn_2)
		populate_from_raw_aligment(g_1, g_2, trans_matrix, emiss_matrix)

'''
print(trans_matrix)
print()
print(emiss_matrix)
'''
proba_from_counts(trans_matrix, emiss_matrix)

formating = {	
				'prefix' : '+', \
				'separator' : ' ', \
				'formatter' : {'float' : (lambda x : ' {:.6f}'.format(x))}
			}
			
tm_file = open("../data/tran_mat.txt", 'w')
tm_file.write(np.array2string(trans_matrix, **formating).replace('[', '').replace(']', ''))
tm_file.close()

em_file = open("../data/emis_mat.txt", 'w')
em_file.write(np.array2string(emiss_matrix, **formating).replace('[', '').replace(']', ''))
em_file.close()



