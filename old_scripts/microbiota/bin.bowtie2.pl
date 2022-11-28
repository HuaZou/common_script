#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Long;

my ($file, $cfg, $dir, $idx, $prefix, $bam, $rpath, $help, $version); 
GetOptions(
    "f|file:s"    => \$file,
    "c|cfg:s"     => \$cfg,
    "d|dir:s"     => \$dir,
    "i|idx:s"     => \$idx,
    "p|prefix:s"  => \$prefix,
	"b|bam:s"	  => \$bam,
	"r|rpath:s"	  => \$rpath,	
    "h|help:s"    => \$help,
    "v|version:s" => \$version
);
print &version && exit if(defined $version);
print &usage && exit if((!defined $prefix)||(defined $help));

system "mkdir -p -m 755 $dir" unless (-d $dir);
my $index_dir = "00.index";
system "mkdir -p -m 755 $index_dir" unless (-d $index_dir);
my $bam_dir = "02.bam";
system "mkdir -p -m 755 $bam_dir" unless(-d $bam_dir);

my $script = "/hwfssz1/ST_OCEAN/USER/zhangpengfan/software/Anaconda/anaconda2/bin/bowtie2";
my $fqdir = "/zfssz3/ST_AGRIC/P18Z10200N0112/root_micro_db/";  

# contigs 
my %contig;
open(CTG, "< $cfg") or die "can't open $cfg\n";
while(<CTG>){
    chomp;
    my $tmp =(split("\/", $_))[-1];
    if($tmp=~/(\S+)\.contigs/){
        my $name = $1;
        $contig{$name} = $_;
    }
}
close(CTG);

# reads path
my %reads;
open(RP, "< $rpath") or die "can't open $rpath\n";
while(<RP>){
	chomp;
	my @tmp = split("\/", $_);
	my $sample = pop(@tmp);
	my $sid = $1 if $sample =~ /^(\D+\d+)\.[pair|single|fq]/;
	my $path = join("\/", @tmp);
	$reads{$sid} = $path;
}
close(RP);

# bowtie2 
my(@bioproject, @mode, %hash1, %hash2);
open(IN, &OpenM($file)) or die "can't open $file\n";
open(OUT, "> $prefix\.sh") or die "can't open $prefix\n";
open(IDX, "> $idx\.sh") or die "can't open $idx\n";
open(BM, "> $bam\.sh") or die "can't open $bam\n";
<IN>;
while(<IN>){
    chomp;
    #$tmp[0]->Run;$tmp[1]->bioproject;$tmp[2]->PE|SE;$tmp[3]->model;$tmp[4]->number 
    my @tmp = split("\t", $_);
    if ($tmp[2] =~ /PAIRED/) {# PE sequence model
        # clean fq path 
        my $fq1 = join("", $fqdir, join(".", $tmp[0], "pair.1.fq.gz")); 
        my $fq2 = join("", $fqdir, join(".", $tmp[0],"pair.2.fq.gz"));
        my $fqs = join("", $fqdir, join(".", $tmp[0], "single.fq.gz"));
		# judge whether the reads are here or not
		unless(-e $fq1){
			if(exists $reads{$tmp[0]}){
				$fq1 = "$reads{$tmp[0]}"."/"."$tmp[0]".".pair.1.fq.gz";
				$fq2 = "$reads{$tmp[0]}"."/"."$tmp[0]".".pair.2.fq.gz";
				$fqs = "$reads{$tmp[0]}"."/"."$tmp[0]".".single.fq.gz";
			}	
		}
        if ($tmp[3] =~ /single/){
            if(exists $contig{$tmp[0]}){
                my $sampleID = $tmp[0]; 
                my $configpath = $contig{$tmp[0]};
                my $mkidx = &mkindex($configpath, $sampleID);
                if(defined $mkidx){
                    print IDX "$mkidx\n";
                }   
				# bowtie2
                # one project has more than two contigs
                my $index = "00.index/".$sampleID;
                if($tmp[4] > 1){
                    # projectID; sampleID; configID; 
                    $hash1{$tmp[1]}{$tmp[0]}{$configpath} = 1;
                    $hash2{$tmp[1]}{$configpath} = 1;
                } else {
                    my $outscript = &scriptout($dir, $sampleID, $index, $script, $fq1, $fq2, $fqs);
                    if(defined $outscript){
                        print  OUT "$outscript\n";
                    }
                    # sam to bam
                    my $sb = &sam2bam($dir, $tmp[0]);
                    if(defined $sb){
						print BM "$sb\n";
					}
                }              

			}
        } elsif($tmp[3] =~ /mix/){ 
            if(exists $contig{$tmp[1]}){
                my $sampleID = $tmp[0];
                my $configpath = $contig{$tmp[1]};
                my $mkidx = &mkindex($configpath, $tmp[1]);
                if(defined $mkidx){
                    print IDX "$mkidx\n";
                } 
                # bowtie2
                my $index = "00.index/".$tmp[1];
				my $outscript = &scriptout($dir, $sampleID, $index, $script, $fq1, $fq2, $fqs);
                if(defined $outscript){
                    print  OUT "$outscript\n";
                }
                # sam to bam
                my $sb = &sam2bam($dir, $tmp[0]);
                if(defined $sb){
					print BM "$sb\n";
				}
            }          
        }        
    } elsif ($tmp[2] =~ /SINGLE/){ # SE sequence model
        my $fq = join("", $fqdir, join(".", $tmp[0], "fq.gz"));
		unless(-e $fq){
			if(exists $reads{$tmp[0]}){
				$fq = "$reads{$tmp[0]}"."/"."$tmp[0]".".fq.gz";	
			}
		}
        if ($tmp[3] =~ /single/){
            if(exists $contig{$tmp[0]}){
                my $sampleID = $tmp[0];
                my $configpath = $contig{$tmp[0]};
                my $mkidx = &mkindex($configpath, $sampleID);
                if(defined $mkidx){
                    print IDX "$mkidx\n";
                } 
                my $index = "00.index/".$sampleID;           
                my $outscript = &scriptout($dir, $sampleID, $index, $script, $fq);          
                if(defined $outscript){
                    print OUT "$outscript\n";
                } 
                # sam to bam
                my $sb = &sam2bam($dir, $tmp[0]);
                if(defined $sb){
					print BM "$sb\n";
				}
			}
        } elsif($tmp[3] =~ /mix/){   
            if(exists $contig{$tmp[1]}){
                my $sampleID = $tmp[0];
                my $configpath = $contig{$tmp[1]};
                my $mkidx = &mkindex($configpath, $tmp[1]);
                print IDX "$mkidx\n";
                my $index = "00.index/".$tmp[1];
                my $outscript = &scriptout($dir, $sampleID, $index, $script,$fq);
                if(defined $outscript){
                    print  OUT "$outscript\n";
                } 
                # sam to bam
                my $sb = &sam2bam($dir, $tmp[0]);
                if(defined $sb){
					print BM "$sb\n";
				}
			}              
        }
    }

}
close(IN);
close(IDX);

# foreach hash
foreach my $proj(keys %hash1){
    foreach my $samp (keys %{$hash1{$proj}}){
        foreach my $cong (keys %{$hash2{$proj}}){
            my $fq1 = join("", $fqdir, join(".", $samp, "pair.1.fq.gz"));
            my $fq2 = join("", $fqdir, join(".", $samp,"pair.2.fq.gz"));
            my $fqs = join("", $fqdir, join(".", $samp, "single.fq.gz"));
 			unless(-e $fq1){
				if(exists $reads{$samp}){
                	$fq1 = "$reads{$samp}"."/"."$samp".".pair.1.fq.gz";
					$fq2 = "$reads{$samp}"."/"."$samp".".pair.2.fq.gz";
					$fqs = "$reads{$samp}"."/"."$samp".".single.fq.gz";
				}
			}
			my $prf = $1  if $cong =~ m/\_fa\/(\S+)\.contigs/;
            my $sampleID = "$prf"."$samp";
            my $configpath = $cong;
            my $index = "00.index/".$samp;
            my $outscript = &scriptout($dir, $sampleID, $index, $script, $fq1, $fq2, $fqs);
            #print("$outscript\n");
            if(defined $outscript){
                print OUT "$outscript\n";
            }
            my $sb = &sam2bam($dir, $sampleID);
            if(defined $sb){
				print BM "$sb\n";
			}
        }
    }
}
close(OUT);
close(BM);

# ####################
# SUB FUNCTION
# ####################
sub OpenM {
    $_=shift;
    return(($_=~/\.gz$/)?"gzip -dc $_|":"$_")
}

sub scriptout{
    # samout;sample;contig;
    my($dir, $name, $path, $scpt, $fasq, $faq1, $faq2, $faqs, $res);
    my $len = @_;
    if ($len == 5) {
        $dir  = shift @_;
        $name = shift @_;
        $path = shift @_;
        $scpt = shift @_;
        $fasq = shift @_;
        my $sam = join(".", $name, "sam");
        my $sam_out = join("\/", $dir, $sam);              
        #unless(-e $sam_out){
			$res = join(" ", $scpt, "-x", $path, "-U", $fasq, "-S", $sam_out, "--no-unal --no-discordant -p 30 -k 50");
        #}
    } elsif($len == 7){
        $dir  = shift @_;
        $name = shift @_;
        $path = shift @_;
        $scpt = shift @_;
        $faq1 = shift @_;
        $faq2 = shift @_;
        $faqs = shift @_;
        my $sam = join(".", $name, "sam");
		my $sam_out = join("\/", $dir, $sam);
        unless(-e $sam_out){
			$res = join(" ", $scpt, "-x", $path, "-1", $faq1, "-2", $faq2, "-U", $faqs, "-S", $sam_out, "--no-unal --no-discordant -p 30 -k 50");
        }       
    }
    return($res)
}

sub mkindex{
    my($cfp, $sam, $res);
    $cfp = shift @_;
    $sam = shift @_;
    my $index = "00.index/".$sam;
    my $index2 = "$index".".rev.1.bt2";
    unless(-e $index2){
        $res = join(" ", "bowtie2-build", "$cfp", "00.index/", $sam);
    }
    return($res)
}

sub sam2bam{
    my($out, $sap, $sam, $bam_name1, $bam_name2, $res);
    $out = shift @_;
    $sap = shift @_;
    $sam = "$out"."/"."$sap".".sam";
    $bam_name1 = "02.bam"."/"."$sap"."\.tmp.bam";
    $bam_name2 = "02.bam"."/"."$sap"."\.sort.bam";
    my $soft = "sambamba_v0.6.6";
    if(-e $sam){
		unless(-e $bam_name2){
        	$res = join(" ", $soft, "view", $sam, "-S -f bam -t 30 -o", $bam_name1, "\n", $soft, "sort -o", $bam_name2, "-t 30", $bam_name1);
    	}    
    	return($res);
	}
}

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <input> -c <ctg> -i <index> -b <bam> -d <dir> -p <prefix> -r <readpath> 
options:
	-f|file		:[essential]sample information (Run|project|library|model).
	-c|ctg		:[essential] contig more than 2K.
	-i|index	:[essential] make fa into index for bowtie2.
	-b|bam		:[essential] convert sam into bam.
	-d|outdir	:[essential output directory path. Conatins the dir. 
	-p|prefix	:prefix of shell script.
	-r|readpath	:reads path
USAGE
    exit;
};

sub version {
    print <<VERSION;
    version:    v1.1
    update:     20190311
    author:     zouhua\@genomics.cn
VERSION
};
