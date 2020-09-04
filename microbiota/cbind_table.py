#!/usr/bin/python3

#----------------------------------------------------------------------------#
# Copyright (c) 2018 Hua Zou (BGI-shenzhen). Allrights reserved              #
# This program is used to cbind sample profile into a table 		         #
# Args :  3 args                                                             #
#   one profile:	profile : rownames->gene ; colnames->sampleID			 #
#   column:			taxonomy reference										 #
#   Output: 		director											     #
#                                                                            #
# Output: one result		                                                 #
#	profile table: colnames are sampleID; rownames are taxonomy			     #
#----------------------------------------------------------------------------#

import sys
import re
import gzip

try:
	import argparse
except ImportError:
	sys.exit("Unable to find the argparse module")


def parse_arguments(args):
	parser = argparse.ArgumentParser(
		description="taxonomy annotation")
	parser.add_argument('-p', '--profile', metavar='<table>',
		help="one profile rownames->taxonomy ; colnames->type of abundance",
		required=True)
	parser.add_argument('-c','--cols', metavar='<column>',
		help="columns to be annotated by abundance",
		required=True)
	parser.add_argument('-o','--output', metavar='<output>',
		help="combind profile",
		required=True)
	return parser.parse_args()


def DictFun(input1):
	"""
	dictionary for storage
	"""
	storage = {}
	f = OpenMethod(input1)
	lines = f.readlines()
	for line in lines:
		sample = line.strip().split("\t")
		if len(sample) >= 2:
			# k->sampleID  v-> path
			storage[sample[0]] = sample[1]
	f.close()
	return(storage)


def CaculateValue(anno, dict, sample):
	if anno[sample[0]] in dict:
		# foreach values and add each
		for key in dict.keys():
			value = [str(float(i) + float(j)) for i, j in zip(dict[key], sample[1:])]
			dict[key] = value
	else:
		dict[anno[sample[0]]] = sample[1:]
	return(dict)


def AnnotationProfile(gene, ref, out):
	"""
	extract gene profile and generate annotate files
	"""
	result = {}
	annotation = DictFun(ref)
	profile = OpenMethod(gene)
	head = profile.readline().strip().split()
	lines = profile.readline()

	# Read a Text File Line by Line Using while
	while lines:
		lines = profile.readline()
		sample = lines.strip().split()
		# judge profile geneID in GeneKO
		if len(sample) != 0:
			# geneID in dictionary
			if sample[0] in annotation:
				result = CaculateValue(annotation, result, sample)
			else:
				result = CaculateValue(annotation, result, sample)
	profile.close()
	# output annotation profile
	output = open(out, "w")
	output.write("\t" + "\t".join(head[1:]))
	for k in sorted(result.keys()):
		abundance = "\t".join(result[k])
		output.write("\n" + k + "\t" + abundance)
	output.close()


def OpenMethod(infile):
	result = gzip.open(infile, "rb") if re.findall(".gz", infile) else open(infile, "r")
	return(result)


def main():
	args = parse_arguments(sys.argv)
	AnnotationProfile(args.gene, args.reference, args.output)


main()
