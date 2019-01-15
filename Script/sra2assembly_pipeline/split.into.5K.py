#!/usr/bin/env /ldfssz1/ST_META/share/flow/anaconda3/bin/python3

import sys
import os

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
        description="DESCRIPTION\n")
    parser.add_argument(
        '-f', '--input', metavar='<file>', type=str,
        help="fasta files path\n",
        required=True)
    parser.add_argument(
        '-o', '--output', metavar='<file>', type=str,
        help="splited fasta files path\n",
        required=True)
    return parser.parse_args()


def split_fasta_5K(file_path, out_file):
    """
    split fasta file into 5K files
    outfile contains 5K files path
    """
    with open(file_path, "r") as f_path, open(out_file, "w") as out:
        for path in f_path:
            path = path.strip()
            if os.path.isfile(path):
                split_dir = path + ".split"
                if not os.path.isdir(split_dir):
                    shell = "seqkit split -s 5000 --quiet " + path
                    os.system(shell)
                if os.path.isdir(split_dir):
                    split_files = os.listdir(split_dir)
                    for split_file in split_files:
                        split_file_path = split_dir + "/" + split_file
                        out.write(split_file_path + "\n")


def main():
    args = parse_arguments(sys.argv)
    split_fasta_5K(args.input, args.output)


main()
