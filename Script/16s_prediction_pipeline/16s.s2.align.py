#!/usr/bin/env /ldfssz1/ST_META/share/flow/anaconda3/bin/python3

# ==============================================================================
# 16s s2 alignment:
#               1. cbind predict result
#               2. rename fasta
#               3. remove redundary fasta
#               4. alignment into E.coli K12
#
# Authors: ZouHua
#
# Please type "./16s.s2.alignment.py -h" for usage help
#
# ==============================================================================


__author__ = ('ZouHua (zouhua@genomics.cn)')
__version__ = '0.1'
__date__ = '20190109'

import sys
import os
from Bio import SeqIO

try:
    import argparse as ap
except ImportError:
    sys.exit("Unable to find argparse module")

vsearch = "/hwfssz1/ST_OCEAN/USER/zhangpengfan/software/vsearch-2.6.0/bin/vsearch"
align = "/hwfssz1/ST_OCEAN/USER/zhangpengfan/software/muscle3.8.31_i86linux64"


def parse_arguments(args):
    """
    1. infile: 16s predict result
    2. E.coli K12 genome fasta
    3. output directory
    4. output script
    """

    parser = ap.ArgumentParser(
        description="DESCRIPTION\n"
        "16s.s2.alignment version "+__version__+" ("+__date__+"): \n"
        "16s alignment on E.coli K12 \n"
        "AUTHORS: "+__author__+"\n\n"
        "",
        formatter_class=ap.RawTextHelpFormatter,
        prog="16s.rDNA.predict.py")
    parser.add_argument(
        '-f', '--infile', metavar='<predict file>', type=str,
        help="fasta files\n",
        required=True)
    parser.add_argument(
        '-g', '--genome', metavar='<genome file>', type=str,
        help="E.coli K12 fasta\n",
        required=True)
    parser.add_argument(
        '-d', '--dir', metavar='<directory>', type=str,
        help="output directory of 16s cbind\n",
        required=True)
    parser.add_argument(
        '-o', '--out', metavar='<out script>', type=str,
        help="output of script\n",
        required=True)
    parser.add_argument(
        '-v', '--version', action='version',
        version="16s.s2.alignment version {} ({})".format(__version__, __date__),
        help="Prints the current 16s.s2.alignment version and exit")
    return parser.parse_args()


def read_file(infile, genome):
    """
    this function is cbind files into a list
    """
    fasta_files = []
    with open(infile, "r") as f:
        for line in f:
            line = line.strip()
            if os.path.isfile(line):
                fasta_files.append(line)
    fasta_files.append(genome)
    return(fasta_files)


def cbind_content(files, out):
    """
    1. cbind contents of files into one file
    2. rename sample fasta by order
    """
    with open(out, "w") as outfile:
        for fname in files:
            sample_name = os.path.basename(fname)
            number = 1
            records = SeqIO.parse(fname, 'fasta')
            for record in records:
                name = sample_name + str(number)
                record.id = name
                record.name = ""
                record.description = ""
                number = number + 1
                SeqIO.write(record, outfile, "fasta")


def makedir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)


def main():
    args = parse_arguments(sys.argv)
    cwd = os.getcwd()   # command running path

    # create output dir
    out_dir = cwd + '/' + args.dir
    makedir(out_dir)
    # file list
    filenames = read_file(args.infile, args.genome)
    # cbind file content
    out_file = out_dir + "/" + "16s.all.rename.fa"
    cbind_content(filenames, out_file)
    # create output script for redundary and alignment
    script_file = out_dir + "/" + args.out
    sf = open(script_file, "w")
    # remove redundary file
    unique_file = cwd + "/" + args.dir + "/16s.all.rename_uniq.fa"
    r_shell = vsearch + " --derep_fulllength " + out_file + " --output " + unique_file
    a_shell = align + " -in " + unique_file + " -out " + args.dir + "/16s.all.rename_uniq.align"
    sf.write(r_shell + "\n" + a_shell + "\n")
    sf.close()

    # run shell
    os.system(r_shell)
    if os.path.isfile(unique_file):
        os.system(a_shell)


main()
