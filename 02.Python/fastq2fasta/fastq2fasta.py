#!/usr/bin/env python

import sys
import argparse as ap
from Bio import SeqIO


def parse_arguments(args):
    parser = ap.ArgumentParser(description='convert fq into fa.')
    parser.add_argument(
        '-i', '--input', metavar='<file>', type=str,
        help="fastq file\n")
    parser.add_argument(
        '-o', '--output', metavar='<file>', type=str,
        help="fasta file\n")
    return parser.parse_args()


def fastq2fasta(fastq_file, fasta_file):
    SeqIO.convert(fastq_file, 'fastq', fasta_file, 'fasta')


def main():
    args = parse_arguments(sys.argv)
    fastq2fasta(args.input, args.output)


main()
