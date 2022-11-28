#!/usr/bin/env /ldfssz1/ST_META/share/flow/anaconda3/bin/python3

# ==============================================================================
# Check.remove v1:
#               1. check fasta file size
#               2. fasta into gene
#
# Authors: ZouHua
#
# Please type "./gene.prediction.py -h" for usage help
#
# ==============================================================================


__author__ = ('ZouHua (zouhua@genomics.cn)')
__version__ = '0.2'
__date__ = '20181227-20190109'

import sys
import re
import os

try:
    import argparse as ap
except ImportError:
    sys.exit("Unable to find argparse module")

cwd = os.getcwd()   # command running path
script = "/hwfssz1/ST_OCEAN/USER/zhangpengfan/software/Anaconda/anaconda2/bin/prodigal"


def parse_arguments(args):
    """
    parameters input
    """

    parser = ap.ArgumentParser(
        description="DESCRIPTION\n"
        "gene.prediction version "+__version__+" ("+__date__+"): \n"
        "gene prediction on fasta by prodigal software \n"
        "AUTHORS: "+__author__+"\n\n"
        "",
        formatter_class=ap.RawTextHelpFormatter,
        prog="gene.prediction.py")
    parser.add_argument(
        '-f', '--infile', metavar='<infile>', type=str,
        help="fasta files\n",
        required=True)
    parser.add_argument(
        '-d', '--dire', metavar='<direction>', type=str,
        help="direction of fasta\n",
        required=True)
    parser.add_argument(
        '-o', '--out', metavar='<output>', type=str,
        help="output of gene prediction\n",
        required=True)
    parser.add_argument(
        '-v', '--version', action='version',
        version="gene.prediction version {} ({})".format(__version__, __date__),
        help="Prints the current gene.prediction version and exit")
    return parser.parse_args()


def gene_prediction(infile, dire, out):
    """
    generate shell script for gene prediciton
    """
    outfile = open(out, "w")
    with open(infile, 'r') as f:
        for line in f:
            tmp = line.strip()
            sample = re.split(r'\s+', tmp)
            name = re.split("/", sample[0])[-2]
            size = file_size(sample[0])
            tmp2 = re.split(r'\s+', size)
            # size more than 10 MB
            if tmp2[0] > 10 and tmp2[1] in ['MB', 'GB']:
                genepath = "/".join([cwd, dire, name])
                makedir(genepath)
                filenumber = judge_file(genepath)
                if filenumber != 3:
                    shell = shell_command(genepath, name, sample[0])
                    outfile.write(shell+"\n")
    outfile.close()


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
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)


def judge_file(path):
    """
    determine prodigal step had been done or not by files
    """
    True_list = 0
    files = os.listdir(path)
    for tmp in files:
        filepath = "/".join([path, tmp])
        if os.path.isfile(filepath):
            True_list += 1
    return(True_list)


def shell_command(path, name, fasta):
    path_file = path + "/"
    # specify output file
    file1 = "".join([path_file, name, ".gene"])
    # protein translations
    file2 = "".join([path_file, name, ".pro.aa"])
    # nucleotide sequences
    file3 = "".join([path_file, name, ".nucleo.fa"])
    res = " ".join([script, "-i", fasta, "-o", file1, "-a", file2, "-d", file3,
                    "-p meta -f gff"])
    return(res)


def main():
    args = parse_arguments(sys.argv)
    gene_prediction(args.infile, args.dire, args.out)


main()
