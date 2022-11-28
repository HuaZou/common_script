#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';

my $cwd = abs_path;
my ($file, $dir, $out1, $out2);
GetOptions(
    "f|file:s"    =>  \$file,
    "d|dir:s"     =>  \$dir,
    "o1|out1:s"   =>  \$out1,
    "o2|out2:s"   =>  \$out2
);
&usage if(!defined $out2);

open(IN, $file) or die "can't open $file";
open(OT1, "> $out1") or die "can't open $out1\n";
open(OT2, "> $out2") or die "can't open $out2\n";
while(<IN>){
    chomp;
    my @tmp = split("\t", $_);
    # download failed file
    my $download_dir = join("/", $dir, "sra", $tmp[0]);
    my $failed_file = join(".", $tmp[0], "partial");
    my $succeeded_file = join("/", $download_dir, $tmp[0]);
    my $failed_filepath = join("/", $download_dir, $failed_file);
    my $failed_dir = join("/", $dir, "sra", $tmp[0]);
    unless(-e $succeeded_file){
        print OT2 "$tmp[0]\n";
        if(-e $failed_filepath){
            print OT1 "$tmp[0]\t$failed_filepath\t$failed_dir\n";
        }
        if(-e $download_dir){
            unless(-e $failed_filepath){
                print OT1 "$tmp[0]\t$download_dir\t$download_dir\n";    
            }
        }
    }
}
close(IN);
close(OT1);
close(OT2);

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <file> -d <dir> -o1 <out1> -o2 <out2>
options:
    -f|file     :[essential] the sra accession list.
    -d|dir      :[essential] the directory of download sra.
    -o1|out1    :[essential] the download failed acession list.
    -o2|out2    :[essential] the update sra accession list.
USAGE
    exit;
}

