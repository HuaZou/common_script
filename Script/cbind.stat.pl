#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';

my $cwd = abs_path;
my ($file, $dir, $type, $prefix);
GetOptions(
    "f|file:s"          => \$file,
    "d|dir:s"           => \$dir,
    "t|type gz file:s"  => \$type, 	
    "p|prefix:s"        => \$prefix   
);
print &usage && exit if(!defined $prefix);

open(IN, &OpenM($file)) or die "can't open $file\n";
open(OUT, "> $cwd/$prefix\.txt") or die "can't open outfile\n"; 
print OUT "Run\tA\tT\tC\tG\tBaseNum\tReadNum";
<IN>;
while(<IN>){
    chomp;
    my @tmp = split("\t", $_);
    # $tmp[1] -> Run; $tmp[7] -> bioproject; $tmp[9] -> sra path; $tmp[3] -> PE|SE
	my $fqpath = join("\/", $dir, $tmp[7], $tmp[1]); 

	my ($fq1, $fq1res, $fq2, $fq2res, $fqs, $fqsres); 	  
    if ($tmp[3] eq "PAIRED") {
        if ($type eq "sra2fq") {
			# fq1 sra2fq
			$fq1 = join("\_", $fqpath, "1.sra2fq.stat");
			$fq1res = &calculate($fq1); 
			print OUT "\n$tmp[1]\.fq1\t$fq1res";
			# fq2 sra2fq
			$fq2 = join("\_", $fqpath, "2.sra2fq.stat");
			$fq2res = &calculate($fq2); 
			print OUT "\n$tmp[1]\.fq2\t$fq2res";			
		} else {
			# fq1 trim			
			$fq1 = join(".", $fqpath, "fq1.trim.stat");
			$fq1res = &calculate($fq1); 
			print OUT "\n$tmp[1]\.fq1\t$fq1res";
			# fq2 trim
			$fq2 = join(".", $fqpath, "fq2.trim.stat");
			$fq2res = &calculate($fq2); 
			print OUT "\n$tmp[1]\.fq2\t$fq2res";
			# fq single trim
			$fqs = join(".", $fqpath, "single.trim.stat");
			$fqsres = &calculate($fqs); 
			print OUT "\n$tmp[1]\.single\t$fqsres";						
		}
    } else {
		my $fq; 
       if ($type eq "sra2fq") {
			# fq sra2fq		   
			$fq = join(".", $fqpath, "sra2fq.stat");
		} else {
			# fq trim			
			$fq = join(".", $fqpath, "fq.trim.stat");						
		}
		my $fqres = &calculate($fq); 
		print OUT "\n$tmp[1]\.fq1\t$fqres";	
    }
}
close(IN);
close(OUT);

sub calculate { 
	my $fqfile = shift;
	my ($A, $T, $C, $G, $base, $read, $res) = 0;	
	if (-e $fqfile){	
		open(IN2, &OpenM($fqfile)) or die "can't open $fqfile\n";
		my $head = <IN2>;
		my @array = split("\t", $head); 
		while(<IN2>){
			chomp; 
			my @arr = split("\t", $_);
			$A += $arr[1];
			$T += $arr[2];
			$C += $arr[3];
			$G += $arr[4];	
		}
		close(IN2);
		$base = $A+$T+$C+$G;
		$read = $array[0];
		$res = join("\t", $A, $T, $C, $G, $base, $read);	
	} else{
		$res = "NA\tNA\tNA\tNA\tNA\tNA";
	}
	return($res);
}


sub OpenM {
    $_=shift;
    return(($_=~/\.gz$/)?"gzip -dc $_|":"$_")
}

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <input> -d <dir> -t <memory> -p <prefix>
options:
	-f|file		:sample path file (SampleID|fqID|fqPath).
	-d|outdir	:output directory path. Conatins the dir.
	-t|type		:type of gz file.	 
	-p|prefix	:prefix of shell script.
USAGE
    exit;
};
