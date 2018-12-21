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

system "mkdir -p -m 777 $dir" unless (-d $dir);
my $script = "/hwfssz1/ST_AGRIC/P18Z10200N0112/USER/zhangpengfan/projects/metadata/script/PQScan.pl"; 

open(IN, &OpenM($file)) or die "can't open $file\n";
open(OUT, "> $cwd/$prefix\.sh") or die "can't open outfile\n"; 

<IN>;
while(<IN>){
    chomp;
    my @tmp = split("\t", $_);
    # $tmp[1] -> Run; $tmp[7] -> bioproject; $tmp[9] -> sra path; $tmp[3] -> PE|SE    
    my $pro = join("/", $cwd, $dir, $tmp[7]);
    system "mkdir -p -m 777 $pro" unless (-d $pro);
    my $prf = join("/", $pro, $tmp[1]);
    if ($tmp[3] eq "PAIRED") {
        print OUT "$script  $prf\_1\.fastq.gz $prf\_1\.sra2fq.stat 33\n";
        print OUT "$script  $prf\_2\.fastq.gz $prf\_2\.sra2fq.stat 33\n";
    } else {
        print OUT "$script  $prf\.fastq.gz $prf\.sra2fq.stat 33\n";
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
	-f|file		:sample path file (SampleID|fqID|fqPath).
	-d|outdir	:output directory path. Conatins the dir. 
	-p|prefix	:prefix of shell script.
USAGE
    exit;
};
