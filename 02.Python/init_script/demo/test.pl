#!/usr/bin/perl

#----------------------------------------------------------------------------#
# Copyright (c) 2018 Hua Zou (BGI-shenzhen). Allrights reserved              #
# This program Generate origin python/R/perl script with basic funciton      #
# Args :  1 args                                                             #
#   Type:	which type of script to generate								 #
#                                                                            #
# Output: 1 script		                                                     #
#----------------------------------------------------------------------------#


use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use test


my $usage = "Usage: perl $0 -f infile -o outfile";

my ($file, $out)=();
GetOptions(
	'f=s' 	=>\$file,
	'o=s' 	=>\$out,
	help 		=> sub{ pod2usage($usage); }, 
)or die($usage);

if ($file =~ /(\/|)(\S+).gz$/) {
	open(IN, "gzip -dc $file | ") or die($! $file);
}else{
	open(IN, $file) or die($! $file);
}
while(<IN>){
	chomp;
	my @tmp = split;
}


sub PerlFunction{
	statements;
}

