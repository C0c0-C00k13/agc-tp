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
import textwrap
from pathlib import Path
from collections import Counter
from typing import Iterator, Dict, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "MAYI Emilie-Jeanne"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["MAYI Emilie-Jeanne"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "MAYI Emilie-Jeanne"
__email__ = "emilie-jeanne.mayi@etu.u-paris.fr"
__status__ = "Developpement"



def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=Path,
                        default=Path("OTU.fasta"), help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed fasta and extract all fasta sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :return: A generator object that provides the Fasta sequences (str).
    """
    with gzip.open(amplicon_file, "rt") as  monfich:
        for line in monfich:
            # print(line)
            if not line.startswith(">") and len(line) >= minseqlen:
                yield line


def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of sequence with a count >= mincount and a length >= minseqlen.
    """
    sequences_fasta = read_fasta(amplicon_file, minseqlen)
    # Conversion du Generator : List
    liste_fasta = list(sequences_fasta)
    # print(len(liste_fasta))

    # Récupération des séquences uniques
    unique_list = set(liste_fasta)
    # print(len(unique_list))
    
    # Comptage d'occurences
    print("Comptage des occurences")
    output = [[sequence, liste_fasta.count(sequence)] for sequence in unique_list]
    for couple_occurence in output:
        if couple_occurence[1]> minseqlen:
            yield couple_occurence
    

def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences in the format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    score = 0.0
    for align1, align2 in zip(alignment_list[0], alignment_list[1]): 
        if align1 == align2:
            score += 1.0
    return score / max(len(alignment_list[0]), len(alignment_list[1]))

def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """

    # Lecture du fichier de sequences
    iterator_counts = dereplication_fulllength(amplicon_file=amplicon_file,\
        minseqlen=minseqlen,mincount=mincount) 
    liste_counts = list(iterator_counts)

    # Creation d'une liste d'OTU : vide
    liste_OTU = []

    # Parcours chaque séquence
    for index_seq, sequence in enumerate(liste_counts):
        is_OTU = True    
        # Comparaison de la séquence au reste des séquences
        for other_sequence in liste_counts[index_seq:]:
            # Ne pas comparer la séquences à elle-même 
            if other_sequence != sequence:
                # Générer l'alignement
                align = nw.global_align(sequence, other_sequence, gap_open=-1,\
                    gap_extend=-1,matrix=\
                        str(Path(amplicon_file).parent / "MATCH")) # Installer le module MATCH
                # Calcul du pourcentage d'identité
                identite = get_identity(align)
                # Vérification des conditions pour annoter OTU (similarité à
                # au moins 97% et nb_sequence > nb_autre_sequence)
                if identite >= 0.97 and sequence[1] > other_sequence[1]:
                    print("OTU")
                else:
                    print("OTU")
                    is_OTU = False        
            # Incrementer la liste d'OTU
            if is_OTU:
                liste_OTU.append(sequence)


def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    pass


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # # Get arguments
    # args = get_arguments()

    # Votre programme ici
    filename = isfile("./data/amplicon.fasta.gz")
    iterator_seq = read_fasta(filename, minseqlen= 5)
    list_seq = list(iterator_seq)
    # print(len(list_seq))
    print("formation de séquence unique")
    # Récupération des séquences uniques
    unique_list = set(list_seq)

    print("Comptage des occurences")
    print(len([[sequence, list_seq.count(sequence)] for sequence in unique_list]))
    
    # list_unique = dereplication_fulllength(FILENAME, 5, 3)
    # print(len(list_unique))






if __name__ == '__main__':
    main()
