#!/usr/bin/env /ldfssz1/ST_META/share/flow/anaconda3/bin/python3

# ==============================================================================
# gene.trim.check:
#               *. check file exists or not.
#               *. make failed script into one file.
#
#
# Authors: ZouHua
#
# Please type "./gene.trim.check.py -h" for usage help
#
# ==============================================================================


__author__ = ('ZouHua (zouhua@genomics.cn)')
__version__ = '0.1'
__date__ = '20190109'

import sys
import re
import os
import shutil
import time

try:
    import argparse
except ImportError:
    sys.exit("Unable to find argparse module")

cwd = os.getcwd()   # command running path


def parse_arguments(args):
    """
    parameters input
    """

    parser = argparse.ArgumentParser(
        description="DESCRIPTION\n"
        "gene.trim.check version "+__version__+" ("+__date__+"): \n"
        "gene.trim.check status\n"
        "AUTHORS: "+__author__+"\n\n"
        "",
        formatter_class=argparse.RawTextHelpFormatter,
        prog="Check.remove.py")
    parser.add_argument(
        '-f', '--infile', metavar='<infile>', type=str,
        help="the type of script to generate\n",
        required=True)
    parser.add_argument(
        '-t', '--trim', metavar='<trimout>', type=str,
        help="trim out file\n",
        required=True)
    parser.add_argument(
        '-n', '--new', metavar='<update>', type=str,
        help="update file\n",
        required=True)
    parser.add_argument(
        '-v', '--version', action='version',
        version="gene.trim.check version {} ({})".format(__version__, __date__),
        help="Prints the current gene.trim.check version and exit")
    return parser.parse_args()


def 


def main():
    args = parse_arguments(sys.argv)


main()
