#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

#==============================================================
# Main program
#==============================================================
def read_fasta(amplicon_file, minseqlen):
	"""prend deux arguments correspondant au fichier fasta et à la longueur minimale des séquences 
	et retourne un générateur de séquences de longueur l >= minseqlen: yield sequence
	"""
	with open(amplicon_file, 'r') as file:
		seq = []
		for line in file:
			if line[0] == ">" and bool(seq):
				if ((len(seq)-1)*80) + len(seq[-1]) >= minseqlen:
					yield "".join(seq)
				seq = []
			elif line[0] == ">":
				continue
			else:	
				seq.append(line.strip())
		else:
			yield "".join(seq)

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
	"""Prend trois arguments correspondant au fichier fasta,  la longueur minimale des séquences 
	et leur comptage minimum. Elle fait appel au générateur fourni par read_fasta et retourne 
	un générateur des séquences uniques ayant une occurrence O>=mincount ainsi que leur occurrence. 
	Les séquences seront retournées par ordre décroissant d’occurrence: yield [sequence, count]
	"""
	# occurence = {}
	# for seq in read_fasta(amplicon_file,minseqlen):
	# 	if seq not in occurence:
	# 		occurence[seq] = [0,seq]
	# 	occurence[seq][0] += 1
	# 	if occurence[seq][0] >= mincount:
	# 		yield occurence[seq]
	occurence = {}
	for seq in read_fasta(amplicon_file,minseqlen):
		if seq not in occurence:
			occurence[seq] = [0,seq]
		occurence[seq][0] += 1
	occu_sorted = sorted(list(occurence.values()),reverse = True)
	for couple in occu_sorted:
		if couple[0] >= mincount:
			yield [couple[1],couple[0]]

def get_chunks(sequence, chunk_size):
	"""prend une séquence et un longueur de segment l: chunk_size et retourne cette une liste 
	de sous-séquence de taille l non chevauchant. A minima 4 segments doivent être obtenus par séquence
	"""
	nchunk = int(len(sequence)/chunk_size)
	if nchunk < 4:
		raise ValueError
	chunck = []
	for i in range(nchunk):
		chunck.append(sequence[i*chunk_size:(i+1)*chunk_size])
	return chunck

def cut_kmer(sequence, kmer_size):
	"""prend une séquence et une longueur de séquence k et retourne un générateur de tous les mots de longueur 
	k présents dans cette séquence,  yield kmer
	"""	
	pass

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
	"""prend un dictionnaire ayant pour clé un index de kmer et pour valeur une liste d’identifiant 
	des séquences dont ils proviennent. 
	"""
	pass

def search_mates(kmer_dict, sequence, kmer_size):
	"""prend un dictionnaire ayant pour clé un index de kmer et pour valeur une liste d’identifiant 
	des séquences dont ils proviennent, une séquence et une longueur de kmer: kmer_size. 
	"""
	pass

def get_identity(alignment_list):
	"""prend un alignement (sous forme de liste) et calcule le pourcentage d’identité entre 
	es deux séquences selon la formule: id = nb nucléotides identiqueslongueur de l'alignement
	"""
	pass

def detect_chimera(perc_identity_matrix):
	"""prend une matrice donnant par segment le taux d’identité entre la séquence candidate et 
	deux séquences parentes et retourne un booléen indiquant si la séquence candidate est une 
	chimère (True) ou ne l’est pas (False).
	"""
	pass

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
	"""Fait appel au générateur fourni par dereplication_fulllength et retourne un générateur 
	des séquences non chimérique au format: yield [sequence, count]
	"""
	pass

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
	"""Fait appel à chimera removal et retourne un liste d’OTU, cette liste indiquera pour 
	chaque séquence son occurrence (count).
	"""
	pass

def write_OTU(OTU_list, output_file):
	"""Prend une liste d’OTU et le chemin vers un fichier de sortie et affiche les OTU au format:
	>OTU_{numéro partant de 1} occurrence: {nombre d’occurrence à la déréplication}
	{séquence au format fasta}
	"""
	pass

def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    for liste in dereplication_fulllength(args.amplicon_file,200,3):
    	print(liste)

if __name__ == '__main__':
    main()