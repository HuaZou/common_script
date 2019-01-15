#!/usr/bin/env /ldfssz1/ST_META/share/flow/anaconda3/bin/python3

# ==============================================================================
# filter:
#        1. size more than 10M
#        2. trim contigs with 2000bp
#		 3. output result
#
# Authors: ZouHua
#
# Please type "./assembly.filter.py -h" for usage help
#
# ==============================================================================


__author__ = ('ZouHua (zouhua@genomics.cn)')
__version__ = '0.1'
__date__ = '20190108'

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
        "Filtering contigs whose length is less than 150bp \n"
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
        '-o', '--out', metavar='<out>', type=str,
        help="fasta output\n",
        required=True)
    parser.add_argument(
        '-v', '--version', action='version',
        version="assembly filter version {} ({})".format(__version__, __date__),
        help="Prints the current gene.filter version and exit")
    return parser.parse_args()


def contigs_trim(infile, dire, out):
    """
    filter contigs with threshold 2000bp
    """
    outf = open(out, "w")
    contigs = assembly_contigs(infile, dire)
    # key -> trim fa; value -> contigs fa
    for key, value in contigs.items():
        if os.path.isfile(key):
            outf.write(key + "\n")
        if not os.path.isfile(key):
            f = open(key, "w")
            for record in SeqIO.parse(value, "fasta"):
                if len(record.seq) >= 2000:
                    SeqIO.write(record, f, "fasta")
            f.close()
            outf.write(key + "\n")
    outf.close()


def assembly_contigs(infile, dire):
    """
    fasta file exists or size
    """
    trim_fa = {}
    with open(infile, "r") as f:
        for line in f:
            line = line.strip()
            sample = re.split(r'\s+', line)
            if len(sample) > 25:
                # paired
                assembly_dir = sample[21]
                assembly_name = sample[23]
            else:
                # single
                assembly_dir = sample[17]
                assembly_name = sample[19]
            assembly_file = assembly_dir + "/" + assembly_name + ".contigs.fa"
            # file exists or not
            if os.path.isfile(assembly_file):
                assembly_size = file_size(assembly_file)
                size_list = re.split(r'\s+', assembly_size)
                # file size more than 10MB
                if size_list[0] > 10 and size_list[1] in ['MB', 'GB']:
                    dir_new = dire + "_trim_2000"
                    # trim result
                    trim_dir = assembly_dir.replace(dire, dir_new)
                    makedir(trim_dir)
                    trim_file = trim_dir + "/" + assembly_name + ".contigs.2000.fa"
                    trim_fa[trim_file] = assembly_file
    return(trim_fa)


def conver_bytes(num):
    """
    this function will convert bytes into MB.. GB..
    """
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if num < 1024.0:
            return("%3.1f %s" % (num, x))
        num /= 1024.0


def file_size(file_path):
    """
    this function will return the file size
    """
    if os.path.isfile(file_path):
        file_s = os.path.getsize(file_path)
        return(conver_bytes(file_s))


def makedir(path):
    """
    this function will create output folders
    """
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)


def main():
    args = parse_arguments(sys.argv)
    contigs_trim(args.infile, args.dire, args.out)


main()
