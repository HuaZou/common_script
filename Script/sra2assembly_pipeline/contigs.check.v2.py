#!/usr/bin/env /ldfssz1/ST_META/share/flow/anaconda3/bin/python3

# ==============================================================================
# contigs.check.remove v1:
#               1. check file exists or not.
#               2. remove redundary files.
#               3. make failed script into one file.
#
#
# Authors: ZouHua
#
# Please type "./contigs.check.remove.py -h" for usage help
#
# ==============================================================================


__author__ = ('ZouHua (zouhua@genomics.cn)')
__version__ = '0.3'
__date__ = '20181221-20190109'

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
        "contigs.check.remove version "+__version__+" ("+__date__+"): \n"
        "Check file status and remove redundary files\n"
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
        version="contigs.check.remove version {} ({})".format(__version__, __date__),
        help="Prints the current contigs.check.remove version and exit")
    return parser.parse_args()


def check_contigs(infile, out_contigs, update_script):
    """
    1. check assembly All done in log file again
    2. remove redundary files and folders
    3. assembly fasta result
    """
    outfile = open(out_contigs, "w")
    contigs_fa = assembly_shell(infile, update_script)
    # key -> fa prefix; value -> fa folder
    for key, value in contigs_fa.items():
        # assembly log file
        contig_logfile = key + ".log"
        if os.path.exists(contig_logfile):
            with open(contig_logfile, "r") as logf:
                logf_lines = logf.readlines()
                for line in logf_lines:
                    line = line.strip()
                    p = r'ALL DONE'
                    if re.search(p, line):
                        if os.path.isdir(value):
                            finished_fa = remove_finished_file(value)
                            outfile.write(finished_fa + "\n")
    outfile.close()


def remove_finished_file(contig_dir):
    """
    remove the finished temp files
    """
    files = os.listdir(contig_dir)
    for tmp in files:
        pattern = r'(done|intermediate_contigs)'
        files_dir = "/".join([contig_dir, tmp])
        if re.search(pattern, tmp):
            if os.path.isdir(files_dir):
                shutil.rmtree(files_dir)
            else:
                os.remove(files_dir)
                print("OK")
        if re.search("contigs.fa", tmp):
            fafile = files_dir
    return(fafile)


def assembly_shell(shell_file, update):
    """
    1. assembly file exist or not : suffix is "contigs.fa"
    2. megahit: the dir must creat by this software without former one
    3. output nunfinished assembly script
    4. return finished contigs fasta
    """
    update_file = open(update, "w")
    assembly_fa = {}
    with open(shell_file, "r") as f:
        for line in f:
            tmp = line.strip()
            sample = re.split(r'\s+', tmp)
            if len(sample) > 25:
                # pair
                contigs_prefix = "/".join([sample[21], sample[23]])
                contigs_dir = sample[21]
            else:
                # single
                contigs_prefix = "/".join([sample[17], sample[19]])
                contigs_dir = sample[17]
            # assembly contigs file
            contig_fa_file = contigs_prefix + ".contigs.fa"
            if os.path.exists(contig_fa_file):
                # contigs finished
                # key->fa prefix; value->fa folder
                assembly_fa[contigs_prefix] = contigs_dir
            else:
                # contigs unfinished and output script
                if os.path.isdir(contigs_dir):
                    shutil.rmtree(contigs_dir)
                update_file.write(line)
    return(assembly_fa)


def main():
    args = parse_arguments(sys.argv)
    check_contigs(args.infile, args.out, args.new)

    localtime = time.asctime(time.localtime(time.time()))
    print("Local current time :", localtime)


main()
