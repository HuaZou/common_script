#!/usr/bin/python 

import sys, re, os
import hashlib
import argparse as ap 

def parse_argument(args):
	parser = ap.ArgumentParser(description='put the script matrix into batch mode')
	parser.add_argument(
		'-m', '--matrix', metavar='<script matrix file>', type=str,
		help='one line one script in file', required=True)		 
	parser.add_argument(
		'-o', '--out', metavar='<result>', type=str,
		help='bacthing run mode', required=True)

	return parser.parse_args()


def batch_run(infile, outfile):
    out_f = open(outfile, "w")
    in_f = open(infile, 'r')  
    while True:
        line = in_f.readline()
        line = line.strip()
        out_f.write(line + "&")
        if not line:
            break
    out_f.close()


def main():
    args = parse_argument(sys.argv)
    batch_run(args.matrix, args.out)

main()
