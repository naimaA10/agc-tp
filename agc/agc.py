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
import numpy as np
from pathlib import Path
from collections import Counter
from typing import Iterator, Dict, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw
np.int = int

__author__ = "AMMICHE Naïma"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["AMMICHE Naïma"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "AMMICHE Naïma"
__email__ = "naima.ammiche@gmail.com"
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
    sequence = ""
    with gzip.open(amplicon_file, "rt") as  file:
        for line in file :
            if line.startswith(">") : # Si ligne commence par >
                if len(sequence) >= minseqlen :
                    yield sequence
                sequence = ""
            else : # Sinon ajouter la séquence
                sequence += line.strip()          
        if len(sequence) >= minseqlen : # Ajouter la dernière séquence
            yield sequence

def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of sequence with a count >= mincount and a length >= minseqlen.
    """
    seq_occ = {} #[sequence, count]
    for sequence in read_fasta(amplicon_file, minseqlen):
        if sequence not in seq_occ:
            seq_occ[sequence] = 1
        else : 
            seq_occ[sequence] += 1

    seq_occ_sort = sorted(seq_occ, key = seq_occ.get , reverse = True)
    for sequence in seq_occ_sort  :
        if seq_occ[sequence] >= mincount :
            yield [sequence, seq_occ[sequence]]

    

def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences in the format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """

    align = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            align += 1
    return ((align/max(len(alignment_list[0]),len(alignment_list[1]))) *100)


def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int, chunk_size: int, kmer_size: int) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """

    #(list)[sequences, count] dans ordre décroissant
    #(list)[OTU (str), count (int)]
    OTU = []
    f = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    OTU.append(next(f)) #next car générateur

    for seq_count in f:
        ajoute = False  
        for sequence, count in OTU:
            if count > seq_count[1]:
                align = nw.global_align(seq_count[0], sequence, gap_open=-1, gap_extend=-1, matrix=str(Path(__file__).parent / "MATCH"))
                if get_identity(align) <= 97:
                    ajoute = True
        if ajoute :
            OTU.append(seq_count)
    return OTU



def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """ 
    with output_file.open("w") as file :
        for i in range(len(OTU_list)):
            file.write(textwrap.fill(f">OTU_{i+1} occurrence:{OTU_list[i][1]}", width=80))
            file.write("\n")
            file.write(textwrap.fill(f"{OTU_list[i][0]}", width=80))
            file.write("\n")
#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici
    file = args.amplicon_file
    lenmin = args.minseqlen
    min_count = args.mincount
    outputfile = args.output_file

    #print(f"file {args.amplicon_file}\nminseqlen {args.minseqlen}")

    #1
    #print(read_fasta(file,lenmin))
    #for i in read_fasta(file,lenmin):
    #   print(i)

    #3
    np.int = int 
    align = nw.global_align("ATTATATCGCCCG", "ATTCGCGCGATATC", gap_open=-1, gap_extend=-1, matrix=str(Path(__file__).parent / "MATCH"))
    print(align)
   # for i in abundance_greedy_clustering(file, lenmin, min_count, 0, 0):
   #     print(i)


if __name__ == '__main__':
    main()
