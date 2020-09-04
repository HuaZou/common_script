#!/usr/bin/perl -w
use strict;
################################################################################
#unless(1==@ARGV) {
#    &usage;
#    exit;
#}
################################################################################a
my $user=$ENV{'USER'};
$user=$ARGV[0] if @ARGV==1;
my(@info,$i,@jobs,%hash,%computer);
################################################################################
open IN,"qstat -u $user|" or die "qstat -u $user $!\n";
<IN>;
<IN>;
while(<IN>) {
    chomp;
    last unless(/(\d+)\s+/);
    s/^\s+(\d+)/$1/;
    @info=split;
    $hash{$info[0]}=$info[4];
    push(@jobs,$info[0]);
    if($info[4]=~m/qw/)
    {
        $computer{$info[0]}="non";
    }else{
        $computer{$info[0]}=$info[-2];
    }
}
close IN;
################################################################################
$i=join ",",@jobs;
exit unless($i=~/\d+/);
print "job:\tuser\thard\tinformation\n";
open IN,"qstat -j $i|" or die "qstat -j $i $!\n";
while(<IN>) {
    chomp;
    if(/^job_number:\s+(\d+)/) {
        print $1,"\t";
        $i=$1;
    }elsif(/script_file:\s+(\S+)/){
        #print $1,"\t";
        print "$hash{$i}\t$computer{$i}\t";
        chomp(my $seq=<IN>);
		if($seq=~m/^version/){
			chomp($seq=<IN>);
		}
		if($seq=~m/^project/){
			chomp($seq=<IN>);
		}
		if($seq=~m/^binding/){
                        chomp($seq=<IN>);
                }
		if($seq=~m/^job_type/){
                        chomp($seq=<IN>);
                }
        if($seq=~/^usage/){
            $seq=~s/^usage\s+(\d+):\s+//;
            print "$seq\n";
       }else{print "\n";}
    }
    elsif(/job_name:\s+(\S+)/){
        print $1,"\t";
    }elsif(/owner:\s+(\w+)/) {
        print $1,"\t";
    }elsif(/^hard resource_list:\s+virtual_free=(.*),num_proc=(.*)/) {
        print $1,",p=",$2,"\t";
    }elsif(/^hard resource_list:\s+virtual_free=(.*)/) {
        print $1,"\t";
    }elsif(/^hard resource_list:\s+num_proc=(.*),virtual_free=(.*)/){
		print $2,",p=",$1,"\t";
	}
}
close IN;
################################################################################
=ss
sub usage {
    print STDERR
        "\n
        Description\n
        Start by Fri Feb 18 12:09:28 2011\n
        This script is to view jobs' informations\n
        Usage:  \$perl $0 [user]\n
        Author by zhouyuanjie\@genomics.org.cn\n
        \n"
}
=cut
