#!usr/bin/env python

#----------------------------------------------------------------------------#
# Copyright (c) 2018 Hua Zou (BGI-shenzhen). Allrights reserved              #
# This program Generate origin python/R/perl script with basic funciton      #
# Args :  3 args                                                             #
#   Type:	which type of script to generate								 #
#   packages:	python module or R packages or perl module					 #
#   Output:	which type of script to 										 #
#                                                                            #
# Output: 1 script		                                                     #
#----------------------------------------------------------------------------#


import sys

try:
	import argparse
except ImportError:
	sys.exit("Unable to find the argparse module")

def parse_arguments(args):
	"""
	Parse the arguments from the user 
	"""
	parser = argparse.ArgumentParser(
		description = "initialize  script\n",
		formatter_class = argparse.RawTextHelpFormatter,
		prog = "init")
	# type of script	
	parser.add_argument(
		'-t','--type', metavar='<type>', type=str,
		help = "the type of script to generate\n",
		required=True)
	# module/Library 	
	parser.add_argument(
		'-m','--module', metavar='<module>', type=str,
		help = "module for import or library\n",
		required=True)		
	# output 	
	parser.add_argument(
		'-o','--out', metavar='<output>', type=str,
		help = "output of script\n",
		required=True)
	
	return parser.parse_args()

def Right():
	"""
	To declare Right of script 
	"""

	right1 = "#----------------------------------------------------------------------------#"
	right2 = "# Copyright (c) 2018 Hua Zou (BGI-shenzhen). Allrights reserved              #"
	right3 = "# This program Generate origin python/R/perl script with basic funciton      #"		
	right4 = "# Args :  1 args                                                             #"		
	right5 = "#   Type:	which type of script to generate								 #"		
	right6 = "#                                                                            #"		
	right7 = "# Output: 1 script		                                                     #"		
	right8 = "#----------------------------------------------------------------------------#"

	right = "\n".join([right1, right2, right3, right4, right5, right6, right7, right8])	+ "\n\n\n"	
	return right

def ModuleLibray(kind, packages):
	"""
	import or library in python,perl and R
	"""

	if kind == "python":
		lib = packages.strip().split(",")
		res = ""
		for i in lib:
			res = res + "import "+ i + "\n"
		res = res + "\n\n"
	elif kind == "R":
		lib = packages.strip().split(",")
		res = "pacman::p_load(" + ",".join(lib) + ")\n\n"
	elif kind == "perl":
		lib = packages.strip().split(",")
		res = ""
		for i in lib:
			res = res + "use " + i + "\n"
		res = res + "\n\n"
	return res

def Function(kind):
	"""
	defined function 
	"""

	if kind == "python":
		res = "def PythonFunction(args):" + "\n\t" + "statements" + "\n\n"
	elif kind == "R":
		res = "fun <- function(args){" + "\n\t" + "statements" + "}\n\n"
	elif kind == "perl":
		res = "sub PerlFunction{" + "\n\t" + "statements;\n}" + "\n\n"

	return res	


def PythonScript(Type, Module, Out):
	"""
	generate python script 
	"""

	tmp = Out + ".py"
	try:
		f = open(tmp, "w")
	except IOError:
		print("Error: file can not open python script.")
	
	# envirnment
	f.write("#!/usr/bin/python3" + "\n\n")

	# scrip right declare
	right = Right()
	f.write(right)

	# import module
	module = ModuleLibray(Type, Module)
	f.write(module)

	# def paramater function
	f.write("def parse_arguments(args):" + "\n\t" +
			"parser = argparse.ArgumentParser(" + "\n\t\t" +
			"description = \"initialize  script\"," + "\n\t\t" +
			"formatter_class = argparse.RawTextHelpFormatter," + "\n\t\t" +
			"prog = \"init\")" +"\n\t")
	f.write("parser.add_argument(\'-i\',\'--input\', metavar=\'<input>\', " + "\n\t\t" +
			"help = \"input file\"," + "\n\t\t" +
			"required=True)"+"\n\t")
	f.write("parser.add_argument(\'-o\',\'--output\', metavar=\'<output>\', " + "\n\t\t" +
			"help = \"output file\"," + "\n\t\t" +
			"required=True)"+"\n\t")
	f.write("return parser.parse_args()" +"\n\n")

	# def open function
	f.write("def OpenMethod(args):" + "\n\t" +
			"f = open(args, \"r\")" + "\n\t" +
			"lines = f.readlines()" + "\n\t" +
			"for line in lines:" + "\n\t\t" +
			"samples = line.strip().split()" +"\n\t" +
			"f.close()" + "\n\n")
	
	# def function
	fun = Function(Type)
	f.write(fun)

	# def main function
	f.write("def main():" + "\n\t" +
			"args =  parse_arguments(sys.argv)" + "\n"+
			"main()")
	
	f.close()


def RScript(Type, Module, Out):
	"""
	generate python script 
	"""

	tmp = Out + ".R"
	try:
		f = open(tmp, "w")
	except IOError:
		print("Error: file can not open R script.")
	
	# envirnment
	f.write("#!/usr/bin/R" + "\n\n")

	# scrip right declare
	right = Right()
	f.write(right)
	
	# def library 
	f.write("if (!require(pacman)) {" +"\n\t"+
		"install.packages(\"pacman\", dependencies=T)" +"\n"+
		"}" + "\n")
	module = ModuleLibray(Type, Module)
	f.write(module)
	
	# def args 
	f.write("args <- commandArgs(T)" + "\n\n")

	# def usage 
	f.write("if (length(args) < num) {" +"\n\t"+
		"stop(\"Usage:\")" +"\n" +
		"}"+ "\n\n")
	
	# def args input 
	f.write("args1 <- read.csv(args[1])" + "\n" + 
		"args2 <- args[2]" + "\n\n")

	# def function
	fun = Function(Type)
	f.write(fun)

	f.close()

def PerlScript(Type, Module, Out):
	"""
	generate python script 
	"""	
	tmp = Out + ".pl"
	try:
		f = open(tmp, "w")
	except IOError:
		print("Error: file can not open perl script.")
	
	# envirnment
	f.write("#!/usr/bin/perl" + "\n\n")

	# scrip right declare
	right = Right()
	f.write(right)

	# def use module 
	f.write("use strict;" + "\n" +
		"use warnings;" + "\n" + 
		"use Pod::Usage;" + "\n" +
		"use Getopt::Long;" + "\n")
	module = ModuleLibray(Type, Module)
	f.write(module)	

	# def usages 
	f.write("my $usage = \"Usage: perl $0 -f infile -o outfile\";" + "\n\n")

	# def parameters
	f.write("my ($file, $out)=();" + "\n" +
		"GetOptions(" + "\n\t" + 
		"\'f=s\' 	=>\$file," +"\n\t" + 
		"\'o=s\' 	=>\$out," + "\n\t" +
		"help 		=> sub{ pod2usage($usage); }, " + "\n"+
		")" + "or die($usage);\n\n")
	
	# def open function  
	f.write("if ($file =~ /(\/|)(\S+).gz$/) {" + "\n\t"+
		"open(IN, \"gzip -dc $file | \") or die($! $file);" + "\n"+
		"}else{" + "\n\t" + 
		"open(IN, $file) or die($! $file);" + "\n" + "}" + "\n" +
		"while(<IN>){" + "\n\t"+
		"chomp;"+ "\n\t" + "my @tmp = split;"+"\n}\n\n\n") 

	# def function
	fun = Function(Type)
	f.write(fun)

	f.close()


def main():
	args = parse_arguments(sys.argv)

	if args.type == "python":
		PythonScript(args.type, args.module, args.out)
	elif args.type == "R":
		RScript(args.type, args.module, args.out)
	elif args.type == "perl":
		PerlScript(args.type, args.module, args.out)		

main()