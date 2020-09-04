#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';

my $cwd = abs_path;
my ($file, $dir);
GetOptions(
    "f|file:s"          => \$file,
    "d|dir:s"           => \$dir
);
print &usage && exit if(!defined $dir);

my $dir_all = "/hwfssz1/ST_AGRIC/P18Z10200N0112/USER/zhangpengfan/projects/metadata";
my $fqdir = join("\/", $dir_all, "00.fastq");
my $cldir = join("\/", $dir_all, "01.cleandata");  

open(IN, &OpenM($file)) or die "can't open $file\n";
<IN>;
while(<IN>){
    chomp;
    my @tmp = split("\t", $_);
    # $tmp[1] -> Run; $tmp[7] -> bioproject; $tmp[9] -> sra path; $tmp[3] -> PE|SE
    my $fqpath = join("/", $fqdir, $dir, $tmp[7], $tmp[1]);
    my $clpath = join("/", $cldir, $dir, $tmp[7], $tmp[1]);
    if ($tmp[3] eq "PAIRED") {
        my $fq1 = join("\_", $fqpath, "1.fastq.gz");
        my $fq2 = join("\_", $fqpath, "2.fastq.gz");
        my $cl1 = join(".", $clpath, "pair.1.fq.gz");
        my $cl2 = join(".", $clpath, "pair.2.fq.gz");
        my $cls = join(".", $clpath, "single.fq.gz");

        my @fqfiles = ($cl1, $cl2, $cls); 
        unlink $fq1 if @fqfiles;
        unlink $fq2 if @fqfiles;
        #print "$fq1\n" if @fqfiles;        
        #print "$fq2\n" if @fqfiles;
    } else {
        my $fqp = join(".", $fqpath, "fastq.gz");
        my $clp = join(".", $clpath, "fq.gz"); 
        unlink $fqp if (-e $clp); 
        #print "$fqp\n" if (-e $clp);        
    }
}
close(IN);



sub OpenM {
    $_=shift;
    return(($_=~/\.gz$/)?"gzip -dc $_|":"$_")
}

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <input> -d <dir>
options:
	-f|file		:sample path file (SampleID|fqID|fqPath).
	-d|outdir	:gz files in directory path to be removed.
USAGE
    exit;
};
