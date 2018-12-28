#!/usr/bin/env /ldfssz1/ST_META/share/flow/anaconda3/bin/python3

# ==============================================================================
# gene.filter:
#               1. filter gene whose length is less than 150bp
#
# Authors: ZouHua
#
# Please type "./gene.filter.py -h" for usage help
#
# ==============================================================================


__author__ = ('ZouHua (zouhua@genomics.cn)')
__version__ = '0.1'
__date__ = '28 12 2018'

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
        '-v', '--version', action='version',
        version="gene.filter version {} ({})".format(__version__, __date__),
        help="Prints the current gene.filter version and exit")
    return parser.parse_args()


def filtergene(infile, dire):
    dic = fasta(infile, dire)
    for key, value in dic.items():
        names = value + ".150.nucleo.fa"
        filename = "/".join([cwd, dire, value, names])
        if not os.path.isfile(filename):
            f = open(filename, "w")
            for record in SeqIO.parse(key, "fasta"):
                if len(record.seq) > 150:
                    SeqIO.write(record, f, "fasta")
            f.close()


def fasta(infile, dire):
    fa = {}
    with open(infile, "r") as f:
        for line in f:
            line = line.strip()
            sample = re.split(r'\s+', line)
            name = re.split(r'\/', sample[8])[-2]
            path = "/".join([cwd, dire, name])
            makedir(path)
            fa[sample[8]] = name
    return(fa)


def makedir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)


def main():
    args = parse_arguments(sys.argv)
    filtergene(args.infile, args.dire)


main()
