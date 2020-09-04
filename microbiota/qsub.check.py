#!/usr/bin/env /ldfssz1/ST_META/share/flow/anaconda3/bin/python3

# ==============================================================================
# Check.remove v1:
#               1. check qsub status.
#               2. make failed script into one file.
#
#
# Authors: ZouHua
#
# Please type "./qsub.check.py -h" for usage help
#
# ==============================================================================


__author__ = ('ZouHua (zouhua@genomics.cn)')
__version__ = '0.1'
__date__ = '27 12 2018'

import sys
import re
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
        description="DESCRIPTION\n"
        "qsub_check version "+__version__+" ("+__date__+"): \n"
        "qsub_check status\n"
        "AUTHORS: "+__author__+"\n\n"
        "",
        formatter_class=ap.RawTextHelpFormatter,
        prog="qsub.check.py")
    parser.add_argument(
        '-q', '--qsub', metavar='<qsub>', type=str,
        help="qsub status log file\n",
        required=True)
    parser.add_argument(
        '-o', '--out', metavar='<output>', type=str,
        help="output of failed script\n",
        required=True)
    parser.add_argument(
        '-v', '--version', action='version',
        version="QsubCheck version {} ({})".format(__version__, __date__),
        help="Prints the current QsubCheck version and exit")
    return parser.parse_args()


def CheckFile(infile, newfile):
    """
    check qsub file
    """
    logfile = infile + ".log"
    qsblog = "/".join([cwd, logfile])
    qsb = qsblog.replace(".log", "_qsub")
    newf = open(newfile, 'w')
    with open(logfile, 'r') as logF:
        for line in logF:
            tmp = line.strip()
            if re.findall(r'not', tmp):
                sample = re.split(r'\s+', tmp)
                shfile = "/".join([qsb, sample[0]])
                with open(shfile, 'r') as shf:
                    lines = shf.read().splitlines()
                lines = lines[:-1]
				newf.write(lines[0]+"\n")
    newf.close()


def main():
    args = parse_arguments(sys.argv)
    CheckFile(args.qsub, args.out)


main()
