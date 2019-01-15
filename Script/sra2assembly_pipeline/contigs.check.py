#!/usr/bin/env /ldfssz1/ST_META/share/flow/anaconda3/bin/python3

# ==============================================================================
# Check.remove v1:
#               1. check qsub status.
#               2. remove redundary files.
#               3. make failed script into one file.
#
#
# Authors: ZouHua
#
# Please type "./Check.remove.py -h" for usage help
#
# ==============================================================================


__author__ = ('ZouHua (zouhua@genomics.cn)')
__version__ = '0.2'
__date__ = '21 12 2018'

import sys
import re
import os
import shutil
import time

try:
    import argparse
except ImportError:
    sys.exit("Unable to find argparse module")

# cwd = os.path.dirname(os.path.realpath(__file__)) # script path
cwd = os.getcwd()   # command running path


def parse_arguments(args):
    """
    parameters input
    """

    parser = argparse.ArgumentParser(
        description="DESCRIPTION\n"
        "Check.remove version "+__version__+" ("+__date__+"): \n"
        "Check task status and remove redundary files\n"
        "AUTHORS: "+__author__+"\n\n"
        "",
        formatter_class=argparse.RawTextHelpFormatter,
        prog="Check.remove.py")
    parser.add_argument(
        '-f', '--infile', metavar='<infile>', type=str,
        help="the type of script to generate\n",
        required=True)
    parser.add_argument(
        '-o', '--out', metavar='<output>', type=str,
        help="output of fa file\n",
        required=True)
    parser.add_argument(
        '-n', '--new', metavar='<update>', type=str,
        help="update file\n",
        required=True)
    parser.add_argument(
        '-v', '--version', action='version',
        version="CheckRemove version {} ({})".format(__version__, __date__),
        help="Prints the current CheckRemove version and exit")
    return parser.parse_args()


def MatchChar(infile, outfa):
    outfile = open(outfa, "w")
    with open(infile, "r") as f:
        for line in f:
            tmp = line.strip()
            sample = re.split(r'\s+', tmp)
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
    """
    remove file
    """
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
                sample2 = re.split(r'\s+', lines[0])
                if len(sample2) > 25:
                    rmdir = sample2[21]
                else:
                    rmdir = sample2[16]
                if os.path.isdir(rmdir):
                    shutil.rmtree(rmdir)
                newf.write(lines[0]+"\n")
    newf.close()


def main():
    args = parse_arguments(sys.argv)
    MatchChar(args.infile, args.out)
    CheckFile(args.infile, args.new)

    localtime = time.asctime(time.localtime(time.time()))
    print("Local current time :", localtime)


main()
