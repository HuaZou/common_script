#!/usr/bin/python3

#----------------------------------------------------------------------------#
# Copyright (c) 2018 Hua Zou (BGI-shenzhen). Allrights reserved              #
# This program generate taxonomy profile by gene profile table		         #
# Args :  3 args                                                             #
#   gene_table:		profile : rownames->gene ; colnames->sampleID			 #
#   reference:		taxonomy reference										 #	
#   Output: 		taxonomy profile table									 #
#                                                                            #
# Output: one result		                                                 #
#	taxonomy table: colnames are sampleID; rownames are taxonomy			 #
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
    parser.add_argument(
        '-p', '--gene', metavar='<table>',
        help="gene profile table rownames->gene ; colnames->sampleID",
        required=True)
    parser.add_argument(
        '-r', '--reference', metavar='<annotation>',
        help="reference for annotation with two columns\n",
        required=True)
    parser.add_argument(
        '-o', '--output', metavar='<output>',
        help="annotation profile",
        required=True)
    return parser.parse_args()


# dictionary for storage
def DictFun(input1):
    storage = {}
    f = OpenMethod(input1)
    lines = f.readlines()
    for line in lines:
        sample = line.strip().split("\t")
		if len(sample) >= 2:
			# k->gene  v->annotation
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


# extract gene profile and generate annotate files
def AnnotationProfile(gene, ref, out):
	result = {}

	annotation = DictFun(ref)
	profile = OpenMethod(gene)
	head = profile.readline().strip().split()
	lines = profile.readline()
	# Read a Text File Line by Line Using While
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
    res = gzip.open(infile, "rb") if re.findall(".gz", infile) else open(infile, "r")
    return(res)


def main():
    args = parse_arguments(sys.argv)
    AnnotationProfile(args.gene, args.reference, args.output)


main()
