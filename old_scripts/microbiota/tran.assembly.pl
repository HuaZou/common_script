#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';

my $cwd = abs_path;
my ($file, $dir, $prefix, $max_m, $cpu); 
GetOptions(
    "f|file:s"          => \$file,
    "d|dir:s"           => \$dir,
    "p|prefix:s"        => \$prefix,
    "m|max_memoty:i"    => \$max_m,
    "c|cpu_num:i"       => \$cpu,    
);
print &usage && exit if(!defined $max_m);

system "mkdir -p -m 777 $dir" unless (-d $dir);
my $script = "/hwfssz1/ST_OCEAN/USER/zhangpengfan/software/Anaconda/anaconda2/bin/Trinity";
my $fqdir = "/hwfssz1/ST_AGRIC/P18Z10200N0112/USER/zhangpengfan/projects/metadata/01.cleandata/";  

open(IN, &OpenM($file)) or die "can't open $file\n";
open(OUT, "> $cwd/$prefix\.sh") or die "can't open outfile\n"; 

<IN>;
while(<IN>){
    chomp;
    my @tmp = split("\t", $_);
    # $tmp[1] -> Run; $tmp[7] -> bioproject; $tmp[9] -> sra path; $tmp[3] -> PE|SE
    my $fqpath = join("/", $fqdir, $dir, $tmp[7], $tmp[1]);
    my $pro = join("/", $cwd, $dir, $tmp[7]);
    system "mkdir -p -m 777 $pro" unless (-d $pro); 
    if ($tmp[3] eq "PAIRED") {
        #print OUT "$script filterMeta -1 $fqpath\_1.fastq.gz -2 $fqpath\_2.fastq.gz -o $pro -Q 2 -S -L 15 -N 3 -P 0.5 -q 20 -l 60 -R 0.5 -5 0 -c $tmp[1]\n";
    } else {
        print OUT "$script --seqType fq --max_memory $max_m --single $fqpath\.fq.gz --min_contig_length 150 --CPU $cpu --output $fqpath --no_normalize_reads --run_as_paired --no_version_check\n";
    }
}
close(IN);
close(OUT);


sub usage{
	print <<USAGE;
usage:
	perl $0 -f <input> -d <dir> -p <prefix> -m <memory> -c <cpu>
options:
	-f|file		:sample path file (SampleID|fqID|fqPath).
	-d|outdir	:output directory path. Conatins the dir. 
	-p|prefix	:prefix of shell script.
	-m|memory	:max memory for run assemby.
	-c|cpu		:CPU number.
USAGE
    exit;
};

sub OpenM {
    $_=shift;
    return(($_=~/\.gz$/)?"gzip -dc $_|":"$_")
}
