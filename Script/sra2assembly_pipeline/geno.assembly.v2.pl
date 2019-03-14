#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';

my $cwd = abs_path;
my ($file, $dir, $prefix, $help, $version); 
GetOptions(
    "f|file:s"    => \$file,
    "d|dir:s"     => \$dir,
    "p|prefix:s"  => \$prefix,
    "h|help:s"    => \$help,
    "v|version:s" => \$version
);
print &version && exit if(defined $version);
print &usage && exit if((!defined $prefix)||(defined $help));

system "mkdir -p -m 755 $dir" unless (-d $dir);
my $script = "/hwfssz1/ST_OCEAN/PMO/F13ZOOYJSY1313/huangyueying_backup/software/megahit-1.0.3/megahit";
my $fqdir = "/hwfssz1/ST_AGRIC/P18Z10200N0112/USER/zhangpengfan/projects/metadata/01.cleandata";  

my(@bioproject, @mode, %hash);
open(IN, &OpenM($file)) or die "can't open $file\n";
open(OUT, "> $cwd/$prefix\.sh") or die "can't open outfile\n";
<IN>;
while(<IN>){
    chomp;
    # $tmp[0]->Run; $tmp[1]->bioproject; $tmp[2]->PE|SE; $tmp[3]->model; $tmp[4]->number    
    my @tmp = split("\t", $_);
    # clean fq path 
    my $fqpath = join("/", $fqdir, $dir, $tmp[1], $tmp[0]);
    my $pro = join("/", $cwd, $dir, $tmp[1]);
    my $outdir = join("/", $pro, $tmp[0]); 
    if ($tmp[2] =~ /PAIRED/) {
		my $fq1 = join(".", $fqpath, "pair.1.fq.gz"); 
		my $fq2 = join(".", $fqpath, "pair.2.fq.gz");
		my $fqs = join(".", $fqpath, "single.fq.gz");  
        if ($tmp[3] =~ /single/){
            my $outscript = &scriptout($outdir, $tmp[0], $fq1, $fq2, $fqs, $script);
            if(defined $outscript){
                print OUT "$outscript\n";
            }          
        } elsif($tmp[3] =~ /mix/){
            push(@{$hash{$tmp[1]}{"fq1"}{"PE"}}, $fq1); 
            push(@{$hash{$tmp[1]}{"fq2"}{"PE"}}, $fq2);
            push(@{$hash{$tmp[1]}{"fqs"}{"PE"}}, $fqs);
            push(@bioproject, $tmp[1]);
            push(@mode, "PE");            
        }        
    } elsif ($tmp[2] =~ /SINGLE/){
        my $fq = join(".", $fqpath, "fq.gz");
        if ($tmp[3] =~ /single/){
            my $outscript = &scriptout($outdir, $tmp[0], $fq, $script);          
            if(defined $outscript){
                print OUT "$outscript\n";
            }           
        } elsif($tmp[3] =~ /mix/){   
            push(@{$hash{$tmp[1]}{"fq"}{"SE"}}, $fq);
            push(@bioproject, $tmp[1]);
            push(@mode, "SE");              
        }
    }
}
close(IN);

# fearch hash 
my @bio = uniq(@bioproject);
my @mod = uniq(@mode);

foreach my $key1 (0..$#bio){   
    my $pro = join("\/", $cwd, $dir, $bio[$key1]);
    foreach my $key2 (0..$#mod){
        if ($mod[$key2] =~ /SE/){
            if (exists $hash{$bio[$key1]}{"fq"}{"SE"}){
                my $fq = join(",", @{$hash{$bio[$key1]}{"fq"}{"SE"}});
                my $outscript = &scriptout($pro, $bio[$key1], $fq, $script);         
                if(defined $outscript){
                    print OUT "$outscript\n";
                }
            }

        } elsif($mod[$key2] =~ /PE/){
            if (exists $hash{$bio[$key1]}{"fq1"}{"PE"}){
                my $fq1 = join(",", @{$hash{$bio[$key1]}{"fq1"}{"PE"}});
                my $fq2 = join(",", @{$hash{$bio[$key1]}{"fq2"}{"PE"}});
                my $fqs = join(",", @{$hash{$bio[$key1]}{"fqs"}{"PE"}});
                my $outscript = &scriptout($pro, $bio[$key1], $fq1, $fq2, $fqs, $script);          
                if(defined $outscript){
                    print OUT "$outscript\n";
                }
            }
        }
    }
}
close(OUT);

# ####################
# SUB FUNCTION
# ####################
sub OpenM {
    $_=shift;
    return(($_=~/\.gz$/)?"gzip -dc $_|":"$_")
}

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

sub scriptout{
    my($dire, $pref, $fasq, $faq1, $faq2, $faqs, $scpt, $con, $contigs, $res);
    my $len = @_;
    if ($len == 4) {
        $dire = shift @_;
        $pref = shift @_;
        $fasq = shift @_;
        $scpt = shift @_;
        system "mkdir -p -m 755 $dire" unless (-d $dire);
        # contigs file exist or not                
        $con = join("/", $dire, $pref);
        $contigs = join(".", $con, "contigs.fa");
        #unless (-e $contigs){
        $res = join(" ", $scpt, "-r", $fasq, "-m 100000000000 -t 20 --kmin-1pass --min-count 2 --k-min 27 --k-max 87 --k-step 10 -o", $dire, "--out-prefix", $pref, "--min-contig-len 500");
        #}
    } elsif($len == 6){
        $dire = shift @_;
        $pref = shift @_;
        $faq1 = shift @_;
        $faq2 = shift @_;
        $faqs = shift @_;
        $scpt = shift @_;
        system "mkdir -p -m 755 $dire" unless (-d $dire);        
        # contigs file exist or not        
        $con = join("/", $dire, $pref);
        $contigs = join(".", $con, "contigs.fa");
        $res = join(" ", $scpt, "-1", $faq1, "-2", $faq2, "-r", $faqs, "-m 100000000000 -t 20 --kmin-1pass --min-count 2 --k-min 27 --k-max 87 --k-step 10 -o", $dire, "--out-prefix", $pref, "--min-contig-len 500");        
    }
    return($res)
}

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <input> -d <dir> -p <prefix> 
options:
	-f|file		:[essential]sample information (Run|project|library|model).
	-d|outdir	:[essential output directory path. Conatins the dir. 
	-p|prefix	:prefix of shell script.
USAGE
    exit;
};

sub version {
    print <<VERSION;
    version:    v1.1
    update:     20181210-20190102
    author:     zouhua\@genomics.cn
VERSION
};
