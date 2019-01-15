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
		my $fq1 = join(".", $prf, "pair.1.fq.gz"); 
		my $fq2 = join(".", $prf, "pair.2.fq.gz");
		my $fqs = join(".", $prf, "single.fq.gz");
        print OUT "$script  $fq1 $prf\.fq1.trim.stat 33\n" if (-e $fq1);  
        print OUT "$script  $fq2 $prf\.fq2.trim.stat 33\n" if (-e $fq2);
        print OUT "$script  $fqs $prf\.single.trim.stat 33\n" if (-e $fqs);
    } else {
		my $fq = join(".", $prf, "fq.gz");
        print OUT "$script  $fq $prf\.fq.trim.stat 33\n" if (-e $fq) ;
    }
}


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

