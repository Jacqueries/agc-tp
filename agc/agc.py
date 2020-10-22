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
	if amplicon_file.endswith(".gz"):
    file = gzip.open(
	else:
	   file = open(
	file.close()
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
		print(len(sequence))
		raise ValueError
	chunck = []
	for i in range(4):
		chunck.append(sequence[i*chunk_size:(i+1)*chunk_size])
	return chunck

def cut_kmer(sequence, kmer_size):
	"""prend une séquence et une longueur de séquence k et retourne un générateur de tous les mots de longueur 
	k présents dans cette séquence,  yield kmer
	"""
	for i in range(len(sequence)):
		if (i + kmer_size) <= len(sequence):
			yield sequence[i:i+kmer_size]

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
	"""prend un dictionnaire ayant pour clé un index de kmer et pour valeur une liste d’identifiant 
	des séquences dont ils proviennent. 
	"""
	for kmer in cut_kmer(sequence,kmer_size):
		if kmer not in kmer_dict:
			kmer_dict[kmer] = [id_seq]
		else :
			kmer_dict[kmer].append(id_seq)
			kmer_dict[kmer] = list(set(kmer_dict[kmer]).union(set(kmer_dict[kmer])))
	return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
	"""prend un dictionnaire ayant pour clé un index de kmer et pour valeur une liste d’identifiant 
	des séquences dont ils proviennent, une séquence et une longueur de kmer: kmer_size. 
	"""
	return [i[0] for i in Counter([ids for kmer in cut_kmer(sequence, kmer_size) if kmer in kmer_dict for ids in kmer_dict[kmer]]).most_common(8)]


def get_identity(alignment_list):
	"""prend un alignement (sous forme de liste) et calcule le pourcentage d’identité entre 
	es deux séquences selon la formule: id = nb nucléotides identiqueslongueur de l'alignement
	"""
	longueur = len(alignment_list[0])
	identite = 0
	for i in range(longueur):
		if alignment_list[0][i] == alignment_list[1][i]:
			identite +=1	
	return 100*(identite/longueur)

def detect_chimera(perc_identity_matrix):
	"""prend une matrice donnant par segment le taux d’identité entre la séquence candidate et 
	deux séquences parentes et retourne un booléen indiquant si la séquence candidate est une 
	chimère (True) ou ne l’est pas (False).
	"""
	liste_std = 0
	diff = False
	for i in range(len(perc_identity_matrix)):
		liste_std += statistics.stdev([perc_identity_matrix[i][0],perc_identity_matrix[i][1]])
		if i > 0 :
			if perc_identity_matrix[i][0] != perc_identity_matrix[i-1][0] or perc_identity_matrix[i][1] != perc_identity_matrix[i-1][1]:
				diff = True
	return bool((liste_std/i) > 5) and diff

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
	"""Fait appel au générateur fourni par dereplication_fulllength et retourne un générateur 
	des séquences non chimérique au format: yield [sequence, count]
	"""
	kmer_dict = {}
	liste = []
	for i,couple_s_o in enumerate(dereplication_fulllength(amplicon_file,minseqlen,mincount)):
		kmer_dict =  get_unique_kmer(kmer_dict,couple_s_o[0],i,kmer_size)
		liste.append(couple_s_o)
	print(liste)
	print("couple_s_o")
	for i,couple_s_o in enumerate(dereplication_fulllength(amplicon_file,minseqlen,mincount)):
		chunk_cand = get_chunks(couple_s_o[0],chunk_size)
		huit_seq = []
		for chunk in chunk_cand:
			huit_seq.append(search_mates(kmer_dict,chunk,kmer_size))
		commun = huit_seq[0]
		for j in range(len(huit_seq)-1):
			commun = common(commun,huit_seq[j+1])
		if len(commun) < 2:
			print("NO PARENTS DETECTED")
			continue
		else:
			parents = commun[0:2]
		chunk_parents = get_chunks(list(dereplication_fulllength(amplicon_file,minseqlen,mincount))[parents[0]][0],chunk_size)
		chunk_parents.append(list(dereplication_fulllength(amplicon_file,minseqlen,mincount))[parents[1]][0])
		perc_identity_matrix = [[] for c in range(len(chunk_cand))]
		for npar in range(len(chunk_parents)):
			for l,chunk in enumerate(chunk_cand):
				perc_identity_matrix[l].append(get_identity(nw.global_align(chunk, chunk_parents[npar][l],gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),'../agc')) + "/MATCH")))
		if not detect_chimera(perc_identity_matrix):
			yield couple_s_o

def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
	"""Fait appel à chimera removal et retourne un liste d’OTU, cette liste indiquera pour 
	chaque séquence son occurrence (count).
	"""
	otu = []
	# liste_seq = sorted(list(chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)),reverse = False)
	# print(liste_seq)
	for i,seq in enumerate(chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)):
		if i == 0:
			otu.append(seq)
			print(seq)
		else:
			for seq_otu in otu:
				idt = get_identity(nw.global_align(seq_otu[0], seq[0],gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),'../agc')) + "/MATCH"))
				if idt <= 97:
					otu.append(seq)
				else:
					print("Identity (ajouter à l'OTU pre existante)")
	return otu


def write_OTU(OTU_list, output_file):
	"""Prend une liste d’OTU et le chemin vers un fichier de sortie et affiche les OTU au format:
	>OTU_{numéro partant de 1} occurrence: {nombre d’occurrence à la déréplication}
	{séquence au format fasta}
	"""
	with open(output_file, 'w') as output:
		for i,otu in enumerate(OTU_list):
			print()
			output.write(">OTU_{} occurrence: {}\n{}\n".format(i,otu[1],fill(otu[0])))

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    otu_list = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, args.chunk_size, args.kmer_size)
    write_OTU(otu_list,args.output_file)

if __name__ == '__main__':
    main()