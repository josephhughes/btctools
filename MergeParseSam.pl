#!/usr/bin/perl

# use this to parse VarScan outputs and determine overlapping mutations 
# IUPAC
# W	weak	A			T	
# S	strong		C	G	
# M	amino	A	C		
# K	keto			G	T
# R	purine	A		G	
# Y	pyrimidine		C		T
# B	not A (B comes after A)		C	G	T	
# D	not C (D comes after C)	A		G	T
# H	not G (H comes after G)	A	C		T
# V	not T (V comes after T and U)	A	C	G	
# N or -	any base (not a gap)	A	C	G	T	



use strict;
use Getopt::Long; 
use Bio::SeqIO;

my %IUPAC;
$IUPAC{"W"} = "T A";
$IUPAC{"S"} = "C G";
$IUPAC{"M"} = "A C";
$IUPAC{"K"} = "G T";
$IUPAC{"R"} = "A G";
$IUPAC{"Y"} = "C T";
$IUPAC{"B"} = "C G T";
$IUPAC{"D"} = "A G T";
$IUPAC{"H"} = "A C T";
$IUPAC{"V"} = "A C G";
$IUPAC{"N"} = "A C G T";

        # flu
# Ambiguous Amino Acids	3-Letter	1-Letter
# Asparagine or aspartic acid	Asx	B
# Glutamine or glutamic acid	Glx	Z
# Leucine or Isoleucine	Xle	J
# Unspecified or unknown amino acid	Xaa	X
        #--------------------------#
        
    my %c2p;
    $c2p{$_} = "Z" for qw(SAA);
    $c2p{$_} = "J" for qw(MTT MTA);
    $c2p{$_} = "B" for qw(AAT AAC GAC GAT RAT RAC);
    $c2p{$_} = "L" for qw(CTA CTT CTG CTC CTK TTA TTG TTR YTA YTG);
    $c2p{$_} = "R" for qw(CGA CGT CGC CGG AGG AGA AGR MGA MGG);
    $c2p{$_} = "S" for qw(TCA TCG TCT TCC AGT AGC AGY);
    $c2p{$_} = "A" for qw(GCC GCT GCA GCG);
    $c2p{$_} = "G" for qw(GGC GGT GGA GGG);
    $c2p{$_} = "P" for qw(CCA CCT CCG CCC);
    $c2p{$_} = "T" for qw(ACA ACG ACC ACT);
    $c2p{$_} = "V" for qw(GTA GTC GTG GTT);
    $c2p{$_} = "I" for qw(ATT ATC ATY ATA ATW);
    $c2p{$_} = "_" for qw(TAA TAG TAR TGA);
    $c2p{$_} = "*" for qw(TAA TAG TAR TGA);
    $c2p{$_} = "C" for qw(TGT TGC TGY);
    $c2p{$_} = "D" for qw(GAT GAC GAY);
    $c2p{$_} = "E" for qw(GAA GAG GAR);
    $c2p{$_} = "F" for qw(TTT TTC TTY);
    $c2p{$_} = "H" for qw(CAT CAC CAY);
    $c2p{$_} = "K" for qw(AAA AAG AAR);
    $c2p{$_} = "N" for qw(AAT AAC AAY);
    $c2p{$_} = "Q" for qw(CAA CAG CAR);
    $c2p{$_} = "Y" for qw(TAT TAC TAY);
    $c2p{$_} = "M" for qw(ATG);
    $c2p{$_} = "W" for qw(TGG);
    $c2p{$_} = "X" for qw(...);#TYA CAM
    

# input will be a text-tab-delimited file with the names of the parsesam files and corresponding reference
# the alignment of all the genes of interest

# go through the text-tab to find the names of parsesam files

my ($files,$out,$refaligned,%ref,%nuccnt,%refseq,%genes,$UTRfile,%aUTR,%UTR);
&GetOptions(
	    'files:s'      => \$files,#tab delimited format
	    'out:s'         => \$out, # output 
	    'refaligned:s'    => \$refaligned, # reference in fasta format comma separated
	    'UTR:s'  => \$UTRfile, #text-tab file with name of gene and length of UTR
           );

open (UTR,"<$UTRfile")||die "Can't open $UTRfile\n";
while (<UTR>){
  chomp($_);
  my @elements=split(/\t/,$_);
  $aUTR{$elements[3]}{$elements[0]}{"start"}=$elements[1];
  $aUTR{$elements[3]}{$elements[0]}{"stop"}=$elements[2];
}
my (@parsefiles);
my @alignments=split(/,/,$refaligned);# all the different gene alignments
open (FILES,"<$files")|| die "Can't open $files\n";
while(<FILES>){
  chomp($_);
  my @elements=split(/\t/,$_);
  $ref{$elements[0]}=$elements[1]; #adding the name of the ParseSam.txt file as a key
  push(@parsefiles,$elements[0]);
}

my (%length,@seqids,%sequences,%parsefileref);
foreach my $ParseFile (keys %ref){
  open (PARSESAM, "<$ParseFile")||die "Can't open $ParseFile\n";
  while (<PARSESAM>){
    chomp($_);
     my @elements=split(/\t/,$_);#Chr Position Acnt Ccnt Tcnt Gcnt Total
     $nuccnt{$ParseFile}{$elements[0]}{$elements[1]}{"A"}=$elements[2];
     $nuccnt{$ParseFile}{$elements[0]}{$elements[1]}{"C"}=$elements[3];
     $nuccnt{$ParseFile}{$elements[0]}{$elements[1]}{"T"}=$elements[4];
     $nuccnt{$ParseFile}{$elements[0]}{$elements[1]}{"G"}=$elements[5]; 
     $nuccnt{$ParseFile}{$elements[0]}{$elements[1]}{"total"}=$elements[2]+$elements[3]+$elements[4]+$elements[5]; 
     $genes{$elements[0]}{$ParseFile}++; #this will provide the length of the unaligned gene and the list of parsefiles with the same reference
  }
}


foreach my $alignment (@alignments){
  my $seq_in  = Bio::SeqIO->new(-format => 'fasta',-file   => $alignment);
	while( my $seq = $seq_in->next_seq() ) {
	  my $id=$seq->display_id();
	  push (@seqids,$id);
	  my $seq_str=$seq->seq();
	  $sequences{$id}=$seq_str;
	  $UTR{$id}=$aUTR{$alignment};
	  $length{$alignment}=length($seq_str);
	  my @bases = split(//,$seq_str);
	  for (my $i = 0; $i<scalar(@bases); $i++){
	    my $site=$i+1;
	    $refseq{$alignment}{$id}{$site}=$bases[$i];# $id corresponds to Chr above
	  }
	}
}

# go through each alignment and check all mutation possibilities
# store the name of the gene id and the mutation and aamutation for each alignment site
my (%gap,%table,%nuctable);
foreach my $alignment (keys %refseq){
  print "Parsing $alignment...\n";
  # check how many proteins are present 
  my @prot=keys %{$aUTR{$alignment}};
  print "number of proteins @prot: ".scalar(@prot)."\n";
  
  
  for (my $i = 0; $i<$length{$alignment}; $i++){
        printf STDOUT "Alignment position %04d\r",$i;
  	    foreach my $id (keys %{$refseq{$alignment}}){
          my $site = $i + 1;
          $gap{$id}{$site}=$gap{$id}{$site-1};
          if ($refseq{$alignment}{$id}{$site}=~/-/){
            $gap{$id}{$site}++;#corresponds to the number of gaps seen before this site in the alignment
          }
		  foreach my $prot (keys %{$aUTR{$alignment}}){            
		    if ($aUTR{$alignment}{$prot}{"start"}<=$site && $aUTR{$alignment}{$prot}{"stop"}>=$site){
		      #print "Coding region $site\n";
              my ($aamut,$type)=aamut($site,$aUTR{$alignment}{$prot}{"start"},"A",$sequences{$id});
          	  $nuctable{$prot}{$id}{$site}{"A"}{"aamut"}=$aamut;
          	  $nuctable{$prot}{$id}{$site}{"A"}{"type"}=$type;
          	  my ($aamut,$type)=aamut($site,$aUTR{$alignment}{$prot}{"start"},"C",$sequences{$id});
          	  $nuctable{$prot}{$id}{$site}{"C"}{"aamut"}=$aamut;
              $nuctable{$prot}{$id}{$site}{"C"}{"type"}=$type;
              my ($aamut,$type)=aamut($site,$aUTR{$alignment}{$prot}{"start"},"T",$sequences{$id});
              $nuctable{$prot}{$id}{$site}{"T"}{"aamut"}=$aamut;
              $nuctable{$prot}{$id}{$site}{"T"}{"type"}=$type;
              my ($aamut,$type)=aamut($site,$aUTR{$alignment}{$prot}{"start"},"G",$sequences{$id});
              $nuctable{$prot}{$id}{$site}{"G"}{"aamut"}=$aamut;
              $nuctable{$prot}{$id}{$site}{"G"}{"type"}=$type;
            }else{
          	  $nuctable{$prot}{$id}{$site}{"A"}{"aamut"}="NA";
          	  $nuctable{$prot}{$id}{$site}{"A"}{"type"}="NA";
          	  $nuctable{$prot}{$id}{$site}{"C"}{"aamut"}="NA";
              $nuctable{$prot}{$id}{$site}{"C"}{"type"}="NA";
              $nuctable{$prot}{$id}{$site}{"T"}{"aamut"}="NA";
              $nuctable{$prot}{$id}{$site}{"T"}{"type"}="NA";
              $nuctable{$prot}{$id}{$site}{"G"}{"aamut"}="NA";
              $nuctable{$prot}{$id}{$site}{"G"}{"type"}="NA";
            }
          }
	    }
  }
}
my @nucs=qw/A C T G/; 

foreach my $alignment (@alignments){
  print "Printing output file $alignment\n"; 
  foreach my $prot (keys %{$aUTR{$alignment}}){   
    open (OUT,">$out\_$prot\.txt")||die "Can't open $out\_$prot\.txt\n";
    #print OUT "Gene\tSite\t";
    print OUT "Site\t";
    foreach my $id (%{$refseq{$alignment}}){
      my @sameref = keys %{$genes{$id}};
      if (@sameref){
      print OUT "$ref{$sameref[0]} Mutation\t$ref{$sameref[0]} aaMutation\t";
      foreach my $parsefile (@sameref){
        #print OUT "$ParseFile Mutation\t$ParseFile aaMutation\t$ParseFile MutationType\t$ParseFile MutationFreq\t$ParseFile Coverage\t";
        print OUT "$parsefile MutationFreq\t";
      }
      }
    }
    print OUT "\n";
    for (my $i = 0; $i<$length{$alignment}; $i++){
      my $site = $i+1;
      
      foreach my $nuc (@nucs){
        print OUT "$site\t";
        foreach my $id (keys %{$refseq{$alignment}}){
          #print "$site $id Ref base $refseq{$alignment}{$id}{$site}\t>>>$gap{$id}{$site} gaps\t ";
          my @sameref = keys %{$genes{$id}};
          my $nogapsite=$site-$gap{$id}{$site};
          print OUT $refseq{$alignment}{$id}{$site}.$nogapsite.$nuc."\t";
          print OUT $nuctable{$prot}{$id}{$site}{$nuc}{"aamut"}."\t";
          #print "$nogapsite\t$nuc\n";
          foreach my $parsefile (@sameref){
            if ($nuccnt{$parsefile}{$id}{$nogapsite}{"total"}>0){
              print OUT $nuccnt{$parsefile}{$id}{$nogapsite}{$nuc}/$nuccnt{$parsefile}{$id}{$nogapsite}{"total"}."\t";
            }else{
              print OUT "0\t";
            }
          }
        }
        print OUT "\n";
      }
    }
    close(OUT)
  }
}
        


sub aamut{
  my $site=$_[0];
  my $UTR=$_[1];
  my $mut=$_[2];
  my $seq_str=$_[3];
  # what is the AA change
  my $noUTR=$site-$UTR+1;
  my ($aasite,$mod);
  if ($noUTR==1){
    $mod=1;
  }elsif ($noUTR==2){
    $mod=2;
  }elsif ($noUTR>2){
    $mod = $noUTR % 3;
  }
  
  my ($raa,$saa,$rcodon,$scodon);#amino acid, reference codon, subject codon
  if ($mod==0){
      $rcodon=substr($seq_str, $site-3, 3);
      $scodon=$rcodon;
      $scodon=~s/(\w\w)\w/$1$mut/;
      $raa= $c2p{uc($rcodon)};
      $saa= $c2p{uc($scodon)};
      $aasite = ($noUTR)/3;
      #print "No UTR $noUTR modular $mod $aasite\n";

  }if ($mod==1){
      $rcodon=substr($seq_str, $site-1, 3);
      $scodon=$rcodon;
      $scodon=~s/\w(\w\w)/$mut$1/;
      $raa= $c2p{uc($rcodon)};
      $saa= $c2p{uc($scodon)};
      if ($noUTR==1){
        $aasite=1;
      }else{
        $aasite = ($noUTR+3-$mod)/3;
        #print "No UTR $noUTR modular $mod $aasite\n";
      }
  }if ($mod==2){
      $rcodon=substr($seq_str, $site-2, 3);
      $scodon=$rcodon;
      $scodon=~s/(\w)\w(\w)/$1$mut$2/;
      $raa= $c2p{uc($rcodon)};
      $saa= $c2p{uc($scodon)};
      $aasite = ($noUTR+3-$mod)/3;
      #print "No UTR $noUTR modular $mod $aasite\n";
      if ($noUTR==2){
        $aasite=1;
      }else{
        $aasite = ($noUTR+3-$mod)/3;
        #print "No UTR $noUTR modular $mod $aasite\n";
      }

  }
  my $aamut=$raa.$aasite.$saa;
  my $type="";
  if ($raa eq $saa){
    $type="syn";
  }
  if ($raa ne $saa){
    $type="nonsyn";
  }
  #print "$site\t$mod\t$mut\t$rcodon\t$scodon\t$aamut\t$type\n"; 
  return($aamut,$type);
}


