#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;
#use Cwd 'abs_path';

#my $cwd = abs_path;
my ($file, $site, $dir, $out, $help, $version);
GetOptions(
    "f|file:s"  =>  \$file,
    "s|site:s"  =>  \$site,
    "d|dir:s"   =>  \$dir,
    "o|out:s"   =>  \$out,
    "h|help:s"  =>  \$help,
    "v|version" =>  \$version
);

system "mkdir -p -m 755 $out" unless(-d $out);

open(IN1, $file) or die "can't open $file";
my %file_name;
while(<IN1>){
    chomp;
    my @tmp = split("\t", $_);
    $file_name{$tmp[0]} = $tmp[1];
}
close(IN1);

foreach my $key (keys %file_name){
    my $scp_file_name = join(".", $key, "scp.sh");
    my $remote_path = join(":", $site, $dir);
    my $shell_command = join(" ", "scp", $file_name{$key}, $remote_path);
    my $scp_status = join("'", $key, "scp was succeed");

    open(OT, "> $out/$scp_file_name") or die "can't open $scp_file_name\n";
    print OT "$shell_command\n";
    print OT "echo \"$scp_status\"\n";
    close(OT);
}

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <file> -s <site> -d <dir> -o <out>
options:
	-f|file	    :[essential] the file path and their names (names|path).
	-s|site     :[essential] the remote server site. 
	-d|dir      :[essential] the file for remote path.
	-o|out      :[essential] the output directory of the shell scripts.
USAGE
    exit;
};

sub version {
    print <<VERSION;
    version:    v1.0
    update:     20200703 - 20200703
    author:     zouhua1\@outlook.com
VERSION
};
