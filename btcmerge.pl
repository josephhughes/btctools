#!/usr/bin/perl

# use this to merge the files from two runs of btcutils 

use strict;
use Getopt::Long; 
use Bio::SeqIO;

# input will be a text-tab-delimited file for the nucleotide frequency, an an AA table
# the alignment of all the genes of interest needs to be provided if they are not aligned against the
# same reference
# if no alignment is provided, then they are considered to be replicates
# go through the text-tab to find the names of parsesam files

my ($files,$out,$refalign,%ref,%nuccnt,%refseq,%genes,%nuctable,%aatable);

&GetOptions(
	    'files:s'      => \$files,#the stubs, comma separated
	    'out:s'         => \$out, # output 
	    'refalign:s'    => \$refalign, # reference in fasta format comma separated, the chr names need to be the same
           );
my (@seqids,%sequences,%length);
my @files=split(/,/,$files);#list of files to merge
my (%chr);
open(OUT,">$out\_entropy.txt")||die "Can't open $out\_entropy.txt\n";
if(!$refalign){
  if (scalar(@files)>2){
    print "Will not be conducting the randomisation or nucletoide frequency test\n";
  }elsif (scalar(@files)==2){
    print "will conduct entropy randomisation base on nucleotide counts\n";
  }else{
    print "expecting two or more input files\n";
  }
}elsif($refalign){
  my (@parsefiles);
  my @alignments=split(/,/,$refalign);# all the different gene alignments
  foreach my $file (@files){
    open (FILE,"<$file\_entropy.txt")|| die "Can't open $file\_entropy.txt\n";
    my $header=<FILE>;
    my @colnames=split(/\t/,$header);
    while(<FILE>){
      chomp($_);
      my @values=split(/\t/,$_);
      for (my $i=0; $i<scalar(@values);$i++){
        $nuctable{$values[1]}{$file}{$values[2]}{$colnames[$i]}=$values[$i];#hash of all the stubs that share a particular gene sequences
        print ">$values[1]< $file $values[2]\n";
        $chr{$values[1]}++;
      }
    }
    close(FILE);
    open (AAFILE,"<$file\_AA.txt")|| die "Can't open $file\_AA.txt\n";
    my $aaheader=<AAFILE>;
    my @aacolnames=split(/\t/,$aaheader);
    while(<AAFILE>){
      chomp($_);
      my @aavalues=split(/\t/,$_);
      for (my $i=0; $i<scalar(@aavalues);$i++){
        $aatable{$file}{$aacolnames[$i]}=$aavalues[$i];
      }
    }
    close(AAFILE);
  }  
  foreach my $alignment (@alignments){
    my $seq_in  = Bio::SeqIO->new(-format => 'fasta',-file   => $alignment);
	while( my $seq = $seq_in->next_seq() ) {
	  my $id=$seq->display_id();
	  push (@seqids,$id);
	  my $seq_str=$seq->seq();
	  $sequences{$id}=$seq_str;
	  $length{$alignment}=length($seq_str);
	  my @bases = split(//,$seq_str);
	  for (my $i = 0; $i<scalar(@bases); $i++){
	    my $site=$i+1;
	    $refseq{$alignment}{$id}{$site}=$bases[$i];# $id corresponds to Chr above
	  }
	}
  }
  for my $id (keys %nuctable){
    for my $file (keys %{$nuctable{$id}}){
      for my $site (keys %{$nuctable{$id}{$file}}){
        for my $info (keys %{$nuctable{$id}{$file}{$site}}){
          print ">$id<\t$file\t$site\t$info\t$nuctable{$id}{$file}{$site}{$info}\n";
        }
      }
    }
  }
        
  
  
  # go through each alignment and check all mutation possibilities
  # store the name of the gene id and the mutation and aamutation for each alignment site
  my (%gap,%table,%nuctable);
  foreach my $alignment (keys %refseq){
    print "Parsing $alignment...\n";
    my $nogapsite;
    for (my $i = 0; $i<$length{$alignment}; $i++){
      
      
      for my $id (keys %{$refseq{$alignment}}){
        my $alpos = $i + 1;
        $gap{$id}{$alpos}=$gap{$id}{$alpos-1};
        $nogapsite=$alpos-($gap{$id}{$alpos-1});
        
        
        #print "$id Alignment position $alpos Number of gaps so far: $gap{$id}{$alpos} Site: $nogapsite\n";
        print OUT "$alpos\t";
	    #print "$alignment $id $refseq{$alignment}{$id}{$alpos}\n";
	    if ($refseq{$alignment}{$id}{$alpos}=~/-/){
	       
	       
	       $gap{$id}{$alpos}++;#corresponds to the number of gaps seen before this site in the alignment
	       print ">$id< $alpos $alignment $refseq{$alignment}{$id}{$alpos} >$gap{$id}{$alpos}< $nogapsite\n";
	       print keys %{$nuctable{$id}},"\n";
	       print keys %chr,"\n";
	       for my $sample (keys %{$nuctable{$id}}){
	         print "$sample\t$id\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"; 
	         
	       }
	       if (%{$nuctable{$id}}){
	         print keys %{$nuctable{$id}},"\n";
	       }else{
	         print "THE VALUE >$id< DOESN T EXISTS\n";
	       
	       }
	    }else{ 
	       print "$id $alpos $alignment $refseq{$alignment}{$id}{$alpos} $gap{$id}{$alpos} $nogapsite\n";
	       for my $file (keys %{$nuctable{$id}}){
	           
	           
	         print OUT "$file\t$id\t$nogapsite\t";
	         for my $sample (keys %{$nuctable{$id}}){
	           for my $info (keys %{$nuctable{$id}{$sample}{$nogapsite}}){
	             print "$nuctable{$id}{$file}{$nogapsite}{$info}\t";
	           }
	         }	        
	       }
	     }
	   }
     }
     print OUT "\n";
  }
}
  
  
  
  