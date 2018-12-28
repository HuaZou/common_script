#!/usr/bin/perl

use strict;
use warnings;

my ($dir, $out) = @ARGV;
if(@ARGV < 2){
 print "Usage :
 			 perl $0 [directory] [output]\n";
}

my (%hash, @id, @geneID);
### open dir
foreach my $file (glob("$dir/*diff")){
  my $name=$1 if $file =~ /gene_expression\/(\S+)\.diff/;
  push(@id, $name);
  open(IN, &openM($file)) or die "can't open $file $!\n";
  <IN>;
  while(<IN>){
     chomp;
     my @tmp = split("\t", $_);
     my $res =  &calculate(\@tmp);
     my $gene = join("link", $tmp[1], $tmp[3]);
     $hash{$gene}{$name} = $res;
     push(@geneID, $gene);
  }
  close(IN);
}

### output title
open(OUT, "> $out/result.re") or die "can't open out $!\n";
print OUT "GeneID\tLocation";
for(my $i=0; $i < @id; $i++){
    my $v1 = join("_", $id[$i], "Value1");
    my $v2 = join("_", $id[$i], "Value2");
    print OUT "\t$v1\t$v2\tFPKM\tFC";
}
print OUT "\n";

### output content
my %hash2;
@geneID = grep {++$hash2{$_} < 2} @geneID;
@geneID = sort {$a cmp $b} @geneID;
for(my $j=0; $j<@geneID; $j++){
    my @tmp = split("link", $geneID[$j]);
    print OUT "$tmp[0]\t$tmp[1]";
    for(my $i=0; $i<@id; $i++){
    $hash{$geneID[$j]}{$id[$i]} ||= "0\t0\tFALSE\tNA";
       print OUT "\t$hash{$geneID[$j]}{$id[$i]}";
    }
    print OUT "\n";
}
close(OUT);


# open file function
sub openM{
   my $f=shift;
   return(($f=~/gz$/)?"gzip -dc $f |":"$f");
}

# calculate function 
sub calculate{
   my ($arr, $FPKM, $FC, $va);
   $arr = shift;
   if($$arr[7] > 0.5 || $$arr[8] > 0.5){
      $FPKM = "TRUE";
   }else{
      $FPKM = "FALSE";
   }
   if($$arr[7]==0 && $$arr[8]==0){
      $FC = "NA";
   }else{
      $FC = sprintf("%.4f",($$arr[7]/($$arr[7]+$$arr[8])));
   }
   $va = join("\t", $$arr[7], $$arr[8], $FPKM, $FC);
   return ($va);
} 
