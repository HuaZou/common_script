#!/usr/bin/env /ldfssz1/ST_META/share/flow/anaconda3/bin/python3

# ==============================================================================
# 16s rDNA predicte:
#               1. predicte 16s rDNA on 2Kbp contigs by rnammer
#               2. result output into dire file
#
#       Notice: software will generate lots of temp file in the directory of
#               the running script,so I separate script and output result.
#               On the other hand, you need to run manually script.
#
# Authors: ZouHua
#
# Please type "./16s.s1.predict.py -h" for usage help
#
# ==============================================================================


__author__ = ('ZouHua (zouhua@genomics.cn)')
__version__ = '0.1'
__date__ = '20190109'

import sys
import re
import os
#import subprocess
import shutil

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
        "16s.s1.predict version "+__version__+" ("+__date__+"): \n"
        "predicte 16s rDNA on 2000bp contigs by rnammer \n"
        "AUTHORS: "+__author__+"\n\n"
        "",
        formatter_class=ap.RawTextHelpFormatter,
        prog="16s.s1.predict.py")
    parser.add_argument(
        '-f', '--infile', metavar='<infile>', type=str,
        help="one colum with fa pathway\n",
        required=True)
    parser.add_argument(
        '-d', '--dire', metavar='<direction>', type=str,
        help="output of 16s predict result\n",
        required=True)
    parser.add_argument(
        '-o', '--out', metavar='<out>', type=str,
        help="shell script of 16s\n",
        required=True)
    parser.add_argument(
        '-v', '--version', action='version',
        version="16s.s1.predict version {} ({})".format(__version__, __date__),
        help="Prints the current 16s.s1.predict version and exit")
    return parser.parse_args()


def predict_rna(fa_file, out_dir, out_shell):
    """
    1. 16s output dir pathway and name
    """
    # 16s result output folder
    out_folder = "/".join([cwd, out_dir])
    makedir(out_folder)


    # 16s output files
    # key -> contigs; value -> out prefix
    predict_dict = fasta(fa_file, out_dir)

    # 16s prepa scirpt path
    out_script_folder = out_folder + "/" + "00.run_tmp"
    makedir(out_script_folder)

    out_script_path = out_script_folder + "/" + out_shell
    out_script = open(out_script_path, "w")
    for key, value in predict_dict.items():
        ssufa = value + "_ssu.fa"
        if not os.path.isfile(ssufa) or os.stat(ssufa).st_size == 0:
            shell = " ".join([script, "-S bac -m ssu -f", ssufa, "<", key])
            out_script.write(shell + "\n")
    out_script.close()

    # backup the script into result dir
    backup_file = out_folder + "/" + out_shell
    shutil.copy(out_script_path, backup_file)


def fasta(infile, dire):
    """
    1. one column: more than 2000 bp contigs pathway
    """
    contig_fa = {}
    with open(infile, "r") as f:
        for line in f:
            line = line.strip()
            name = re.split(r'\/', line)[-2]
            prefix = "/".join([cwd, dire, name])
            contig_fa[line] = prefix
    return(contig_fa)


def makedir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)


def main():
    args = parse_arguments(sys.argv)
    predict_rna(args.infile, args.dire, args.out)

    # run script
    #script_path = "/".join([cwd, args.dire, args.out])
    #subprocess.call(['sh', script_path])


main()
