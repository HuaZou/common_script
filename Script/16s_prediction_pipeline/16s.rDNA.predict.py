#!/usr/bin/env /ldfssz1/ST_META/share/flow/anaconda3/bin/python3

# ==============================================================================
# 16s rDNA predicte:
#               1. predicte 16s rDNA on contigs by rnammer
#
# Authors: ZouHua
#
# Please type "./16s.rDNA.predict.py -h" for usage help
#
# ==============================================================================


__author__ = ('ZouHua (zouhua@genomics.cn)')
__version__ = '0.1'
__date__ = '29 12 2018'

import sys
import re
import os
from Bio import SeqIO

try:
    import argparse as ap
except ImportError:
    sys.exit("Unable to find argparse module")

cwd = os.getcwd()   # command running path
script = "/hwfssz1/ST_OCEAN/USER/zhangpengfan/software/rnamer/rnammer"


def parse_arguments(args):
    """
    parameters input
    """

    parser = ap.ArgumentParser(
        description="DESCRIPTION\n"
        "16s.rDNA.predict version "+__version__+" ("+__date__+"): \n"
        "predicte 16s rDNA on contigs by rnammer \n"
        "AUTHORS: "+__author__+"\n\n"
        "",
        formatter_class=ap.RawTextHelpFormatter,
        prog="16s.rDNA.predict.py")
    parser.add_argument(
        '-f', '--infile', metavar='<infile>', type=str,
        help="fasta files\n",
        required=True)
    parser.add_argument(
        '-d', '--dire', metavar='<direction>', type=str,
        help="direction of fasta\n",
        required=True)
    parser.add_argument(
        '-o', '--out', metavar='<out>', type=str,
        help="output of 16s rDNA prediction\n",
        required=True)     
    parser.add_argument(
        '-v', '--version', action='version',
        version="16s.rDNA.predict version {} ({})".format(__version__, __date__),
        help="Prints the current 16s.rDNA.predict version and exit")
    return parser.parse_args()


def predictrna(infile, dire, out):
    path = "/".join([cwd, dire])
    dic = fasta(infile, dire)
    outf = open(out, "w")
    for key, value in dic.items():
        ssufa = value + "_ssu.fa"
        file1 = "/".join([path, value, ssufa])
        if os.path.isfile(file1):
            print(file1+'\n')
        if not os.path.isfile(file1):
            shell = " ".join([script, "-S bac -m ssu -f", file1, "<", key])
            outf.write(shell + "\n")
    outf.close()


def fasta(infile, dire):
    fa = {}
    with open(infile, "r") as f:
        for line in f:
            line = line.strip()
            sample = re.split(r'\s+', line)
            name = re.split(r'\/', sample[2])[-2]
            path = "/".join([cwd, dire, name])
            makedir(path)
            fa[sample[2]] = name
    return(fa)


def makedir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)


def main():
    args = parse_arguments(sys.argv)
    predictrna(args.infile, args.dire, args.out)


main()
