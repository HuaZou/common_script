#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';

my $cwd = abs_path;
my ($file, $dir, $prefix);
GetOptions(
    "f|file:s"          => \$file,
    "d|dir:s"           => \$dir,	
    "p|prefix:s"        => \$prefix   
);
print &usage && exit if(!defined $prefix);
system "mkdir -p -m 755 $dir" unless (-d $dir);


open(IN, &OpenM($file)) or die "can't open $file\n";
open(OUT, "> $cwd/$prefix\.sh") or die "can't open outfile\n"; 

<IN>;
while(<IN>){
    chomp;
    my @tmp = split("\t", $_);
    if(-e $tmp[2]){
        if ($tmp[1] eq "PAIRED") {
            print OUT "fastq-dump $tmp[2] -O $dir --gzip --defline-seq '\@\$sn\[_\$rn\]/\$ri' --split-files\n";
        } else {
            print OUT "fastq-dump $tmp[2] -O $dir --gzip --defline-seq '\@\$sn\[_\$rn\]/\$ri'\n";
        }
    }
}
close(IN);
close(OUT);


sub OpenM {
    $_=shift;
    return(($_=~/\.gz$/)?"gzip -dc $_|":"$_")
}

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <input> -d <dir> -p <prefix>
options:
	-f|file		:sample path file (SampleID|mode|fqPath).
	-d|outdir	:output directory path. Conatins the dir. 
	-p|prefix	:prefix of shell script.
USAGE
    exit;
};
