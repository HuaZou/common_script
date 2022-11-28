#!/usr/bin/env /ldfssz1/ST_META/share/flow/anaconda3/bin/python3

# ==============================================================================
# gene.filter:
#               1. gene 150bp
#               2. protein 50bp
#
# Authors: ZouHua
#
# Please type "./gene.filter.py -h" for usage help
#
# ==============================================================================


__author__ = ('ZouHua (zouhua@genomics.cn)')
__version__ = '0.2'
__date__ = '20181228-20190108'

import sys
import re
import os
from Bio import SeqIO

try:
    import argparse as ap
except ImportError:
    sys.exit("Unable to find argparse module")

cwd = os.getcwd()   # command running path


def parse_arguments(args):
    """
    parameters input
    """
    parser = ap.ArgumentParser(
        description="DESCRIPTION\n"
        "gene.filter version "+__version__+" ("+__date__+"): \n"
        "Filtering gene whose length is less than 150bp \n"
        "AUTHORS: "+__author__+"\n\n"
        "",
        formatter_class=ap.RawTextHelpFormatter,
        prog="gene.filter.py")
    parser.add_argument(
        '-f', '--infile', metavar='<infile>', type=str,
        help="fasta files\n",
        required=True)
    parser.add_argument(
        '-d', '--dire', metavar='<direction>', type=str,
        help="direction of fasta\n",
        required=True)
    parser.add_argument(
        '-n', '--nucleic', metavar='<nucleic>', type=str,
        help="nucleic fasta output\n",
        required=True)
    parser.add_argument(
        '-p', '--protein', metavar='<protein>', type=str,
        help="protein fasta output\n",
        required=True)
    parser.add_argument(
        '-v', '--version', action='version',
        version="gene.filter version {} ({})".format(__version__, __date__),
        help="Prints the current gene.filter version and exit")
    return parser.parse_args()


def fasta(infile, dire, num):
    """
    consider two types of fasta: nucleic & protein
    """
    fa = {}
    with open(infile, "r") as f:
        for line in f:
            line = line.strip()
            sample = re.split(r'\s+', line)
            name = re.split(r'\/', sample[num])[-2]
            path = "/".join([cwd, dire, name])
            makedir(path)
            fa[sample[num]] = name
    return(fa)


def filter_fasta(fastatype, pathway, suffix, threshold, output):
    """
    filter fasta by threshold
    """
    outf = open(output, "w")

    for key, value in fastatype.items():
        names = value + suffix
        filename = "/".join([cwd, pathway, value, names])
        if os.path.isfile(filename):
            outf.write(value + "\t" + filename + "\n")
        if not os.path.isfile(filename):
            f = open(filename, "w")
            for record in SeqIO.parse(key, "fasta"):
                if len(record.seq) >= threshold:
                    SeqIO.write(record, f, "fasta")
            f.close()
            outf.write(value + "\t" + filename + "\n")

    outf.close()


def makedir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)


def main():
    args = parse_arguments(sys.argv)

    # nucleic fasta
    nfa = fasta(args.infile, args.dire, 8)
    filter_fasta(nfa, args.dire, ".ncleo.150.fa", 150, args.nucleic)

    # protein fasta
    pfa = fasta(args.infile, args.dire, 6)
    filter_fasta(pfa, args.dire, ".aa.50.fa", 50, args.protein)


main()
