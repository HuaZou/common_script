#!/usr/bin/python3

#----------------------------------------------------------------------------#
# Copyright (c) 2018 Hua Zou (BGI-shenzhen). Allrights reserved              #
# This program caculate shannon index  by gene profile table		     #
# Args :  3 args                                                             #
#   gene_table:		profile : rownames->gene ; colnames->sampleID	     #
#   col:	        1,2,3-> reads_pairs,base_abundance,reads_abundance   #	
#   output: 		director and prefix 				     #
#                                                                            #
# Output: one result		                                             #
#	shannon index:  colname is shannon; rownames are types    	     #
#----------------------------------------------------------------------------#

import sys, re, gzip
from math import log
try:
	import argparse
except ImportError:
	sys.exit("Unable to find the argparse module")

def parse_arguments(args):
	parser = argparse.ArgumentParser(
		description = "shannon index caculate")
	parser.add_argument('-p','--profile', metavar='<profile>', 
		help = "profile table rownames->gene ; colnames->sampleID",
		required=True)
	parser.add_argument('-c','--col', metavar='<number>', 
		help = "choose col to calculate",
		required=True)
	parser.add_argument('-o','--output', metavar='<output>', 
		help = "prefix or director",
		required=True)
	return parser.parse_args()

def ShannonIndex(prf, col, output):
	num = int(col)
	shannon = {}
	profile = OpenMethod(prf)
	head = profile.readline().decode("utf-8").strip().split()
	line = profile.readline().decode("utf-8")
	while line:
		line = profile.readline().decode("utf-8")
		sample = line.strip().split()
		if len(sample) != 0 : 
			if num in [1, 2, 3]:
				shannon = CaculateValue(shannon, sample, num, head)
			#else:
			#	for j in [1, 2, 3]:
			#		shannon = CaculateValue(shannon, sample, j, head)	
	profile.close()

	out = open(output, "w")
	out.write("Type\tIndex\n")
	for k, v in shannon.items():
		res = v[1]/v[0] + log(v[0])
		out.write(k + "\t" + str(res) + "\n")
	out.close()

def CaculateValue(Shan, Sample, number, header):
	if Sample[number] != "0":
		count = float(Sample[number])
		index = -float(Sample[number]) * log(float(Sample[number]))
		tmp = [count, index]
		if header[number] in Shan:
			for key in Shan.keys():
				value = [float(i) + float(j) for i, j in zip(Shan[key], tmp)]
				Shan[key] = value
		else:
			Shan[header[number]] = tmp
	return(Shan)

def OpenMethod(infile):
	FILE = gzip.open(infile, "rb") if re.findall(".gz", infile) else open(infile, "r")
	return(FILE)

def main():
	args =  parse_arguments(sys.argv)
	ShannonIndex(args.profile, args.col, args.output)
main()
