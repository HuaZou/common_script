#!/usr/bin/python3

#----------------------------------------------------------------------------#
# Copyright (c) 2018 Hua Zou (BGI-shenzhen). Allrights reserved              #
# This program Generate origin python/R/perl script with basic funciton      #
# Args :  1 args                                                             #
#   Type:	which type of script to generate								 #
#                                                                            #
# Output: 1 script		                                                     #
#----------------------------------------------------------------------------#


import sys
import os
import re


def parse_arguments(args):
	parser = argparse.ArgumentParser(
		description = "initialize  script",
		formatter_class = argparse.RawTextHelpFormatter,
		prog = "init")
	parser.add_argument('-i','--input', metavar='<input>', 
		help = "input file",
		required=True)
	parser.add_argument('-o','--output', metavar='<output>', 
		help = "output file",
		required=True)
	return parser.parse_args()

def OpenMethod(args):
	f = open(args, "r")
	lines = f.readlines()
	for line in lines:
		samples = line.strip().split()
	f.close()

def PythonFunction(args):
	statements

def main():
	args =  parse_arguments(sys.argv)
main()