#!/usr/bin/perl -w
# (c) 2016-2017 Chao IN-HOURSE SHARE ONLY
# fangchao@genomics.cn
# Scan fastQ file and stat position-quality stat table
use warnings;
use strict;
use File::Basename; 

die &usage if @ARGV < 2;
sub usage {
	print <<USAGE;
usage:
se pattern:
	perl $0 <input.fq> <output> <Qual system(33|64)> 
e.g	perl $0 sample.fq sample.pqt 33
USAGE
	exit;
}
sub openM {$_=shift;return(($_=~/\.gz$/)?"gzip -dc $_|":"$_")}

### BODY ###
my ($fq, $out, $Qsys) = @ARGV;
$Qsys ||= 33;

open FQ,&openM($fq) or die "error\n";
open STAT,"> $out",or die "error\n";

my (%STAT, $readNum);
my ($maxP, $maxQ) = (0,0);
while(<FQ>){
	#FQ info
	my ($seq, $num, $qual, $baseNum, $Q, $GC) =();
	chomp;
	chomp($seq = <FQ>);
	chomp($num = <FQ>);
	chomp($qual= <FQ>);
	$readNum ++;
    &scan($seq,$Qsys,$qual);
}
close FQ;

print STAT "$readNum\tA\tT\tC\tG";
for(my $q=0;$q<=$maxQ;$q++){
	print STAT "\t$q";
}
print STAT "\n";
for(my $p=1;$p<=$maxP;$p++){
	print STAT "$p\t$STAT{ATCGN}{$p}{A}\t$STAT{ATCGN}{$p}{T}";
	print STAT "\t$STAT{ATCGN}{$p}{C}\t$STAT{ATCGN}{$p}{G}";
	for(my $q=0;$q<=$maxQ;$q++){
		$STAT{POS}{$p}{$q}||=0;
		print STAT "\t$STAT{POS}{$p}{$q}";
	}
	print STAT "\n";
}

close STAT;
# sub

sub scan {
    my $seq  = shift @_;
    my $sysQ = shift @_;
    my $qual = shift @_;
    my ($s,$p,$l) = (0, 1, length($seq));
    $STAT{info}{baseNum} += $l;
	$maxP = ($maxP < $l)?$l:$maxP;
    while($p<=length($qual)){
		my $p_ = $p -1;
		$s = substr($seq,$p-1,1);
        $_ = substr($qual,$p-1,1);
        $_ = ord($_) - $sysQ;
        $STAT{POS}{$p}{$_} ++;
		$STAT{ATCGN}{$p}{$s} ++;
        $p ++;
		$maxQ = ($maxQ < $_)?$_:$maxQ;
    }
}
