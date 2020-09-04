#!/usr/bin/env /ldfssz1/ST_META/share/flow/anaconda3/bin/python3

import sys
import re
import os
import shutil


try:
    import argparse
except ImportError:
    sys.exit("Unable to find argparse module")

cwd = os.path.dirname(os.path.realpath(__file__))


def parse_arguments(args):
    parser = argparse.ArgumentParser(
        description="remove temporary file and director\n",
        formatter_class=argparse.RawTextHelpFormatter,
        prog="remove.redundancy.py")
    parser.add_argument(
        '-f', '--infile', metavar='<infile>', type=str,
        help="the type of script to generate\n",
        required=True)
    parser.add_argument(
        '-o', '--out', metavar='<output>', type=str,
        help="output of fa file\n",
        required=True)
    return parser.parse_args()


def MatchChar(infile, outfa):
    outfile = open(outfa, "w")
    with open(infile, "r") as f:
        for line in f:
            tmp = line.strip()
            sample = re.split(r"\s+", tmp)
            if len(sample) > 25:   # pair
                prefix = "/".join([sample[21], sample[23]])
                dirpath = sample[21]
                prf = sample[23]
            else:                  # single
                prefix = "/".join([sample[16], sample[18]])
                dirpath = sample[16]
                prf = sample[18]
            logfile = prefix + ".log"
            if os.path.exists(logfile):
                with open(logfile, "r") as lf:
                    lflines = lf.readlines()
                for lfline in lflines:
                    lfline = lfline.strip()
                    p = r'ALL DONE'
                    if re.search(p, lfline):
                        if os.path.isdir(dirpath):
                            fa = RemoveFile(dirpath)
                            outfile.write(prf + "\t" + fa + "\n")
            # else:
                # outfile.write(prf + "\t" + "No contigs.fa\n")
    outfile.close()


def RemoveFile(path):
    files = os.listdir(path)
    for tmp in files:
        p1 = r'(done|intermediate_contigs)'
        filepath = "/".join([path, tmp])
        if re.search(p1, tmp):
            if os.path.isdir(filepath):
                shutil.rmtree(filepath)
            else:
                os.remove(filepath)
        if re.search("contigs.fa", tmp):
            fafile = filepath
    return(fafile)


def main():
    args = parse_arguments(sys.argv)
    MatchChar(args.infile, args.out)


main()
