# A simple script to calculate residue contact matrices in multiple threads

## dependencies
 * executable to find contacts: https://github.com/LilySnow/PDB_related/blob/master/contact-chainID_allAtoms
 * pytorch 2.1 : https://pytorch.org/
 * python 3.10 : https://www.python.org/
 * biopython 1.75 : https://biopython.org/

## running
./count_all.py <path to contact-chainID_allAtoms> <list of pdb files> <number of simultaneous workers> <output matrix file>
