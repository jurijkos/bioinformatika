# Pairwise alignment using Hidden Markov Models(Group /) - Bionformatics

Project for Bioinformatics course at FER, University of Zagreb.

## Dependencies
- **cmake**
```
sudo apt-get install cmake
```
- **python3** with **biopython** (https://github.com/biopython/biopython) and **numpy**
```
pip3 install biopython
sudo pip3 install numpy
```

##Instalation and usage
Clone git repository:
```
git clone https://github.com/jurijkos/bioinformatika
```

Generating `Makefile` and building from source:
```
cd bioinformatika
bash build.sh
make
```

Program outputs the alignment and running time. 
Running HMM aligner:
```
./hmm-pairwise-alignment data/g1.txt data/g2.txt
```

Using flags `-t` and `-e` it's possible to specify path to transition and emission matrices represented
as text-file.
Example:
```
./hmm-pairwise-alignment -t data/tran_mat.txt -e data/emis_mat.txt data/g1.txt data/g2.txt
```



## Optional initialization of parameters

### From clean DNA sequences
To initialize parameters from clean sequences use `matrix_generator.py`with `-g` flag. Script
reads configuration file which specifies sequence pairs to align with N.W. and learn parameters from.
Additional paramater **N** defines who many aligments from each pair wil be used.
Example (has to be run from scripts/` directory):
```
./matrix_generator.py -g ../data/HIV_gene_pairs.txt 1
```
Scripts generates files `data/tran_matrix.txt` and `data/tran_matrix.txt`.

###  From FASTA alignment file
Building transition and emission matrices from FASTA alignment file is done using
**genes_from_alignment.py** along `matrix_generator.py` script. `genes_from_alignment.py`
using flag `-no-gaps` generates raw alignment files from FASTA aligment file provided as argument.
Calling `matrix_generator.py` with flag `-a` and alignment pairs file produces parameter matrices
as textual files. FASTA alignment available [here](https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html).

Example (has to be run from `scripts/` directory`):
```
./genes_from_alignment.py ..data/HIV1_REF_2010_genome_DNA.fasta -gaps
./matrix_generator.py -a ../data/HIV_algnm_pairs.txt
```

### Creating clean sequences from FASTA alignment file
Generating clean DNA sequences from FASTA alignment file using `/genes_from_alignment.py` script
with **-no-gaps** option. Script leavs only A, C, G and T character and creates output file in `data/` folder
from each alignment in original file. 
Example (has to be run from `scripts/` directory):
```
./genes_from_alignment.py ..data/HIV1_REF_2010_genome_DNA.fasta -no-gaps
```

### Coverting raw files to FASTA pair
Generating one file for pair from two separate files:
```
./genes_to_pair.sh ../data/g1.txt ../data/g2.txt g3.txt
```

License
---------
MIT License
---------
2020. [Stjepan Jakovac](https://github.com/sJakovac), [Jurij Kos](https://github.com/jurijkos), [Frane Polić](https://github.com/fPolic)
