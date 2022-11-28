#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';

my $cwd = abs_path;
my ($file, $dir, $out);
GetOptions(
    "f|file:s"  =>  \$file,
    "d|dir:s"   =>  \$dir,
    "o|out:s"   =>  \$out
);
&usage if(!defined $out);

open(IN, $file) or die "can't open $file";
open(OT, "> $out") or die "can't open $out\n";
while(<IN>){
    chomp;
    my @tmp = split("\t", $_);
    my $filename = join(".", $tmp[0], "sra");
    my $path = join("/", $dir, "sra", $filename);
    unless (-e $path){
        print OT "$tmp[0]\n";
    }
}
close(IN);
close(OT);

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <file> -d <dir> -o <out>
options:
    -f|file :[essential] the sra accession list.
    -d|dir  :[essential] the directory of download sra.
    -o|out  :[essential] the update sra accession list.
USAGE
    exit;
}

