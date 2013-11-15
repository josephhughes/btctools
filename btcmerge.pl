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
# example:
# perl btcmerge.pl -files 1351_LPAIH7N1,test -out refalign -refalign SamTestFiles/HA.fa,SamTestFiles/NA.fa,SamTestFiles/M.fa,SamTestFiles/NS.fa,SamTestFiles/PA.fa,SamTestFiles/PB1.fa,SamTestFiles/PB2.fa,SamTestFiles/NP.fa 


my ($files,$out,$refalign,%ref,%nuccnt,%refseq,%genes);
my $reps;
# add argument to specify whether two files are technical replicates
&GetOptions(
	    'files:s'      => \$files,#the stubs, comma separated
	    'out:s'         => \$out, # output 
	    'refalign:s'    => \$refalign, # reference in fasta format comma separated, the chr names need to be the same
# 	    'rep:i'  => \$reps, #number of replicates for the randomisation test of the entropy 
           );
my (@seqids,%sequences,%length);

my @files=split(/,/,$files);#list of files to merge
my (%chr);
if(!$refalign){
  my (%nuctable,@colnames,%aatable,@aacolnames);
  if (scalar(@files)>2){
    open(OUT,">$out\_multi.txt")||die "Can't open $out\_multi.txt\n";
    open(AAOUT,">$out\_AA_multi.txt")||die "Can't open $out\_AA_multi.txt\n";
    #print "Will NOT be conducting the randomisation or nucletoide frequency test\n";
    print OUT "Chr\tSite\t";
    print AAOUT "Protein\tAAPosition\t";
    foreach my $file (@files){
      open (FILE,"<$file\_entropy.txt")|| die "Can't open $file\_entropy.txt\n";
      my $header=<FILE>;
      chomp($header);
      @colnames=split(/\t/,$header);
      while(<FILE>){
        chomp($_);
        my @values=split(/\t/,$_);
        for (my $i=3; $i<scalar(@values);$i++){
          #print "$values[1] $file $values[2]\n";
          $nuctable{$values[1]}{$values[2]}{$file}{$colnames[$i]}=$values[$i];
          #print "$nuctable{$values[1]}{$values[2]}{$file}{$colnames[$i]} $values[1] $values[2] $file $colnames[$i]\n";
        }
      }
      for (my $i=3; $i<scalar(@colnames);$i++){
        print OUT "$file\_$colnames[$i]\t";
      }
      open (AAFILE,"<$file\_AA.txt")|| die "Can't open $file\_AA.txt\n";
      my $aaheader=<AAFILE>;
      #print "$aaheader\n";
      chomp($aaheader);
      @aacolnames=split(/\t/,$aaheader);
      while(<AAFILE>){
        chomp($_);
        my @values=split(/\t/,$_);
        for (my $j=4; $j<scalar(@values);$j++){
          #print "$values[1] $file $values[2]\n";
          $aatable{$values[2]}{$values[3]}{$file}{$aacolnames[$j]}=$values[$j];
          #print "$aatable{$values[2]}{$values[3]}{$file}{$aacolnames[$j]} $values[2] $values[3] $file $aacolnames[$j]\n";
        }
      }
      for (my $j=4; $j<scalar(@aacolnames);$j++){
        print AAOUT "$file\_$aacolnames[$j]\t";
      }
      close(FILE);
      close(AAFILE);
    }
    print OUT "\n";
    print AAOUT "\n";    
    foreach my $gene (keys %nuctable){
      foreach my $site (sort {$a<=>$b} keys %{$nuctable{$gene}}){
        print OUT "$gene\t$site\t";
        foreach my $sample (@files){
          for (my $j=4; $j<scalar(@colnames);$j++){
            print OUT "$nuctable{$gene}{$site}{$sample}{$colnames[$j]}\t";
          }
        } 
        print OUT "\n";   
      }
    }
    close(OUT);
    foreach my $prot (keys %aatable){
      foreach my $site (sort {$a<=>$b} keys %{$aatable{$prot}}){
        print AAOUT "$prot\t$site\t";
        foreach my $sample (@files){
          for (my $i=4; $i<scalar(@aacolnames);$i++){
            print AAOUT "$aatable{$prot}{$site}{$sample}{$aacolnames[$i]}\t";
          }
        } 
        print AAOUT "\n";   
      }
    }
    close(AAOUT);

  }elsif (scalar(@files)==2){
    open(OUT,">$out\_rep.txt")||die "Can't open $out\_rep.txt\n";
    open(AAOUT,">$out\_AA_rep.txt")||die "Can't open $out\_AA_rep.txt\n";
#    print "will conduct entropy randomisation base on nucleotide counts\n";
    print OUT "Chr\tSite\t";
    print AAOUT "Protein\tAAPosition\t";
    foreach my $file (@files){
      open (FILE,"<$file\_entropy.txt")||die "Can't open $file\_entropy.txt\n";
      my $header=<FILE>;
      chomp($header);
      @colnames=split(/\t/,$header);
      while(<FILE>){
        chomp($_);
        my @values=split(/\t/,$_);
        for (my $i=3; $i<scalar(@values);$i++){#
          $nuctable{$values[1]}{$values[2]}{$file}{$colnames[$i]}=$values[$i];
          #print "$nuctable{$values[1]}{$values[2]}{$file}{$colnames[$i]} $values[1] $values[2] $file $colnames[$i]\n";
        }
      }
      #table column names
      for (my $i=3; $i<scalar(@colnames);$i++){
        print OUT "$file\_$colnames[$i]\t";
      }
      open (AAFILE,"<$file\_AA.txt")|| die "Can't open $file\_AA.txt\n";
      my $aaheader=<AAFILE>;
      print "$aaheader\n";
      chomp($aaheader);
      @aacolnames=split(/\t/,$aaheader);
      while(<AAFILE>){
        chomp($_);
        my @values=split(/\t/,$_);
        for (my $j=4; $j<scalar(@values);$j++){
          #print "$values[1] $file $values[2]\n";
          $aatable{$values[2]}{$values[3]}{$file}{$aacolnames[$j]}=$values[$j];
          #print "$aatable{$values[2]}{$values[3]}{$file}{$aacolnames[$j]} $values[2] $values[3] $file $aacolnames[$j]\n";
        }
      }
      for (my $j=4; $j<scalar(@aacolnames);$j++){
        print AAOUT "$file\_$aacolnames[$j]\t";
      }
      close(AAFILE);
      close(FILE);
    }
    print AAOUT "\n";    
    print OUT "Truths";
#     if ($reps){
#       print OUT "\tHdiff(Hb-Hq)\tNbRandH>Hdiff\tHighestRandH\tP-value";
#     }
    print OUT "\n";
    foreach my $gene (keys %nuctable){
      foreach my $site (sort {$a<=>$b} keys %{$nuctable{$gene}}){
        print OUT "$gene\t$site\t";
        foreach my $sample (@files){
          for (my $i=3; $i<scalar(@colnames);$i++){
            print OUT "$nuctable{$gene}{$site}{$sample}{$colnames[$i]}\t";
          }
        }
        my @OrderOfNucs0=split(//,$nuctable{$gene}{$site}{$files[0]}{"OrderOfNucs"});
        my @OrderOfNucs1=split(//,$nuctable{$gene}{$site}{$files[1]}{"OrderOfNucs"});
        my $ordermatch;
        for (my $i=0; $i<4; $i++){
          if ($OrderOfNucs0[$i] eq $OrderOfNucs1[$i] && $OrderOfNucs0[$i]=~/./ && $OrderOfNucs1[$i]=~/./){
            print OUT "Y";
          }else{
            print OUT "N";
          }
        }
#         if ($reps){
#           my $difshannon=$nuctable{$gene}{$site}{$files[0]}{"entropy(base e)"}-$nuctable{$gene}{$site}{$files[1]}{"entropy(base e)"};
#           my (@bbases,@qbases);
#           @bbases = (@bbases, ('A') x $nuctable{$gene}{$site}{$files[0]}{"Acnt"},('C') x $nuctable{$gene}{$site}{$files[0]}{"Ccnt"},('T') x $nuctable{$gene}{$site}{$files[0]}{"Tcnt"},('G') x $nuctable{$gene}{$site}{$files[0]}{"Gcnt"});
#           @qbases = (@qbases, ('A') x $nuctable{$gene}{$site}{$files[1]}{"Acnt"},('C') x $nuctable{$gene}{$site}{$files[1]}{"Ccnt"},('T') x $nuctable{$gene}{$site}{$files[1]}{"Tcnt"},('G') x $nuctable{$gene}{$site}{$files[1]}{"Gcnt"});
# 	      my @randentropies=randomize(\@bbases,\@qbases);
# 	      #print "@randentropies\n";
# 	      my $max = (sort { $b <=> $a } @randentropies)[0];
# 	      my $randcnt=0;
# 	      foreach my $rand (@randentropies){
# 	        if ($rand>$difshannon){
# 	          $randcnt++;
# 	        }
# 	      }
# 	      my $pvalue=$randcnt/$reps;
# 	      print OUT "\t$difshannon\t$max\t$randcnt\t$pvalue";
# 	    }
        print OUT "\n";    
      }
    }
    close(OUT);
    foreach my $prot (keys %aatable){
      foreach my $site (sort {$a<=>$b} keys %{$aatable{$prot}}){
        print AAOUT "$prot\t$site\t";
        foreach my $sample (@files){
          for (my $i=4; $i<scalar(@aacolnames);$i++){
            print AAOUT "$aatable{$prot}{$site}{$sample}{$aacolnames[$i]}\t";
          }
        } 
        print AAOUT "\n";   
      }
    }
    close(AAOUT);

  }else{
    print "Expecting two or more input files\n";
  }
}elsif($refalign){
# try the following approach
# open the reference alignments. Map alignment position to original reference position for each gene/chr (same as input tables), if the alignment position is a gap, then don't put anything in the hash (or put a -)
# open each input table and put gene and reference position into a hash.
# then go through every position of the alignment mapping hash and if the hash value exists then
# check what the corresponding values are in the table hash.
# if the hash value doesn't exist, then put NAs in that row.


  my (%nuctable,%aatable,%newnuc,%newaa,@aacolnames,@colnames,%gappos,%sharedref,%reftoaln);
  foreach my $file (@files){
    open (FILE,"<$file\_entropy.txt")|| die "Can't open $file\_entropy.txt\n";
    my $header=<FILE>;
    chomp($header);
    # Sample	Chr	Position	RefBase	Coverage	AvQual	Acnt	Apval	Ccnt	Cpval	Tcnt	Tpval	Gcnt	Gpval	entropy(base e)	NonRefCnt	CntTv	CntTs	OrderOfNucs
    @colnames=split(/\t/,$header);
    while(<FILE>){
	  chomp($_);
	  my $str=$_;
	  $str =~ s/\s+$//;
	  my @values=split(/\t/,$_);
	  $sharedref{$values[1]}{$file}++;
	  $nuctable{$file}{$values[1]}{$values[2]}=$str;#In NUCLEOTIDE FILE 
    }
    open (AAFILE,"<$file\_AA.txt")|| die "Can't open $file\_AA.txt\n";
    print "$file\_AA.txt\n";
    my $aaheader=<AAFILE>;
    print "$aaheader\n";
    chomp($aaheader);
    @aacolnames=split(/\t/,$aaheader);
    while(<AAFILE>){
	  chomp($_);
	  my $str=$_;
	  $str =~ s/\s+$//;
	  my @values=split(/\t/,$_);
	  # Sample	Chr	Protein	AAPosition	RefAA	RefSite	RefCodon	FstCodonPos	SndCodonPos	TrdCodonPos	CntNonSyn	CntSyn	NbStop	TopAA	TopAAcnt	SndAA	SndAAcnt	TrdAA	TrdAAcnt	AAcoverage
	  $aatable{$file}{$values[1]}{$values[5]}{$values[2]}=$str;#IN AATABLE 
    }
    close(FILE);
    close(AAFILE);
  }
  my @alignments=split(/,/,$refalign);# all the different gene alignments
  foreach my $alignment (@alignments){
    my $seq_in  = Bio::SeqIO->new(-format => 'fasta',-file   => $alignment);
	while( my $seq = $seq_in->next_seq() ) {
	  my $id=$seq->display_id();
	  push (@seqids,$id);
	  my $seq_str=$seq->seq();
	  $sequences{$id}=$seq_str;
	  $length{$alignment}=length($seq_str);
	  $genes{$alignment}{$id}++;
	  my $gapcnt=0;
	  my @bases = split(//,$seq_str);
	  for (my $i = 0; $i<scalar(@bases); $i++){
	    my $alnpos=$i+1;
	    if ($bases[$i]=~/-/){
	      $gapcnt++;
	    }
	    $refseq{$alignment}{$alnpos}{$id}{"site"}=$alnpos-$gapcnt;# from the alignment position get the real position
	    $refseq{$alignment}{$alnpos}{$id}{"base"}=$bases[$i];
	    $gappos{$id}{$alnpos-$gapcnt}=$bases[$i];
	    foreach my $file (keys %{$sharedref{$id}}){
	      if ($bases[$i]=~/-/){
            $newnuc{$alignment}{$alnpos}{$file}="$file\t$id\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
            #print "$alignment\t$id\t$bases[$i]\t$alnpos\t$gapcnt => NA\n";
	      }elsif($bases[$i]!~/-/ && ($alnpos-$gapcnt)>0){
            if ($nuctable{$file}{$id}{($alnpos-$gapcnt)}=~/\w+/){
              $newnuc{$alignment}{$alnpos}{$file}=$nuctable{$file}{$id}{($alnpos-$gapcnt)};
              if (keys %{$aatable{$file}{$id}{($alnpos-$gapcnt)}}){
                for my $prot (keys %{$aatable{$file}{$id}{($alnpos-$gapcnt)}}){
                  $newaa{$alignment}{$prot}{$alnpos}{$file}=$aatable{$file}{$id}{($alnpos-$gapcnt)}{$prot};# adding the protein name
                  print "Protein $alignment $prot $newaa{$alignment}{$prot}{$alnpos}{$file}\n";
                }
                print "END\n";
              }
            }else{
              $newnuc{$alignment}{$alnpos}{$file}="$file\t$id\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
            }
            #print "Populate\t$bases[$i]\t$id\t$alnpos\t$gapcnt\t\t".($alnpos-$gapcnt).">>>>>>>>>>>>$nuctable{$file}{$id}{($alnpos-$gapcnt)}\n";
	      }elsif(($alnpos-$gapcnt)<1){
	        print "Doesn't exist".$alnpos-$gapcnt."\n";
	      }
	    }
	  }
	}
  }
  open(OUT,">$out\_realign.txt")||die "Can't open $out\_realign.txt\n";
  print OUT "Alignment\tAlignPos\t";
  foreach my $file (@files){
    for (my $j=0; $j<scalar(@colnames);$j++){
      print OUT "$file\_$colnames[$j]\t";
    }
  }
  print OUT "\n";

  for my $alignment (keys %newnuc){
    for my $alnpos (sort {$a<=>$b} keys %{$newnuc{$alignment}}){
      print OUT "$alignment\t$alnpos\t";
      foreach my $file (@files){
        print OUT ">$newnuc{$alignment}{$alnpos}{$file}<\t";
      }
      print OUT "\n";
    }
  }
# from the aa refsite position, check what the corresponding alignment position in the newnuc
# put this in a %newaa
# go through every alignment position of the %newaa and if empty for a particular file, add NAs
  open(AAOUT,">$out\_AA_realign.txt")||die "Can't open $out\_AA_realign.txt\n";      
  print AAOUT "Alignment\tAlignPos\t";
  foreach my $file (@files){
    for (my $j=0; $j<scalar(@aacolnames);$j++){
      print AAOUT "$file\_$aacolnames[$j]\t";
    }
  }
  print AAOUT "\n";
  for my $alignment (keys %newaa){
    for my $prot (keys %{$newaa{$alignment}}){
      for my $alnpos (sort {$a<=>$b} keys %{$newaa{$alignment}{$prot}}){
        print AAOUT "$alignment\t$alnpos\t";
        foreach my $file (@files){
          if ($newaa{$alignment}{$prot}{$alnpos}{$file}=~/.+/){
            print AAOUT ">$newaa{$alignment}{$prot}{$alnpos}{$file}<\t";
          }elsif ($newaa{$alignment}{$prot}{$alnpos}{$file}!~/.+/){
            #print AAOUT "|>$file\tGene\t$prot\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA<|\t";
          }
        }
        print AAOUT "\n";
      }
    }
  }


#   my (%nuctable,@colnames,%aatable,@aacolnames,%sharedref,%proteins);
#   open(OUT,">$out\_realign.txt")||die "Can't open $out\_realign.txt\n";
#   open(AAOUT,">$out\_AA_realign.txt")||die "Can't open $out\_AA_realign.txt\n";
#   print OUT "Chr\tAlignPos\t";
#   #print AAOUT "Protein\tAlignPos\t";
#   foreach my $file (@files){
#     open (FILE,"<$file\_entropy.txt")|| die "Can't open $file\_entropy.txt\n";
#     my $header=<FILE>;
#     chomp($header);
#     @colnames=split(/\t/,$header);
#     while(<FILE>){
# 	  chomp($_);
# 	  my @values=split(/\t/,$_);
# 	  for (my $i=1; $i<scalar(@values);$i++){
# 	    #print "$values[1] $file $values[2]\n";
# 	    $nuctable{$values[1]}{$values[2]}{$file}{$colnames[$i]}=$values[$i];#In NUCLEOTIDE FILE chr position samplename columnname and data
# 	    #print "$nuctable{$values[1]}{$values[2]}{$file}{$colnames[$i]} $values[1] $values[2] $file $colnames[$i]\n";
# 	    $sharedref{$values[1]}{$file}++;
# 	  }
#     }
#     open (AAFILE,"<$file\_AA.txt")|| die "Can't open $file\_AA.txt\n";
#     print "$file\_AA.txt\n";
#     my $aaheader=<AAFILE>;
#     print "$aaheader\n";
#     chomp($aaheader);
#     @aacolnames=split(/\t/,$aaheader);
#     while(<AAFILE>){
# 	  chomp($_);
# 	  my @values=split(/\t/,$_);
# 	  for (my $j=1; $j<scalar(@values);$j++){
# 	    # Chr	Protein	AAPosition	RefAA	RefSite	RefCodon	FstCodonPos	SndCodonPos	TrdCodonPos	CntNonSyn	CntSyn	NbStop	TopAA	TopAAcnt	SndAA	SndAAcnt	TrdAA	TrdAAcnt	AAcoverage
# 	    #print "$values[1] $file $values[2]\n";
# 	    # listed all files that have same chr, protein, refsite 
# 	    $aatable{$values[1]}{$values[2]}{$values[5]}{$file}{$aacolnames[$j]}=$values[$j];#IN AATABLE Chr Protein RefSite samplename columname and data
# 	    #print "$aatable{$values[2]}{$values[3]}{$file}{$aacolnames[$j]} $values[2] $values[3] $file $aacolnames[$j]\n";
# 	    $proteins{$values[1]}{$values[2]}++;
# 	  }
#     }
#     close(FILE);
#     close(AAFILE);
#   }
#   #read in the alignment
#   my (%genes,%gappos,%newtable);
#   my @alignments=split(/,/,$refalign);# all the different gene alignments
#   foreach my $alignment (@alignments){
#     my $seq_in  = Bio::SeqIO->new(-format => 'fasta',-file   => $alignment);
# 	while( my $seq = $seq_in->next_seq() ) {
# 	  my $id=$seq->display_id();
# 	  push (@seqids,$id);
# 	  my $seq_str=$seq->seq();
# 	  $sequences{$id}=$seq_str;
# 	  $length{$alignment}=length($seq_str);
# 	  $genes{$alignment}{$id}++;
# 	  my $gapcnt=0;
# 	  my @bases = split(//,$seq_str);
# 	  for (my $i = 0; $i<scalar(@bases); $i++){
# 	    my $alnpos=$i+1;
# 	    if ($bases[$i]=~/-/){
# 	      $gapcnt++;
# 	    }
# 	    $refseq{$alignment}{$alnpos}{$id}{"site"}=$alnpos-$gapcnt;# from the alignment position get the real position
# 	    $refseq{$alignment}{$alnpos}{$id}{"base"}=$bases[$i];
# 	    $gappos{$id}{$alnpos-$gapcnt}=$bases[$i];
# 	  }
# 	}
#   }
#   foreach my $id (keys %{$genes{$alignments[0]}}){
#     print OUT "Gene\t";
#     foreach my $sample (keys %{$sharedref{$id}}){
#       for (my $j=2; $j<scalar(@colnames);$j++){
#         print OUT "$sample\_$colnames[$j]\t";
#       }
#     }
#   }
#   print OUT "\n";
#   # need to repeat the above for each protein
#   # listed all files that have same chr, protein, refsite 
#   #$aatable{$values[1]}{$values[2]}{$values[5]}{$file}{$aacolnames[$j]}=$values[$j];#Chr Protein RefSite
#   print AAOUT "AlignmentPosition\t";
#   foreach my $id (keys %{$genes{$alignments[0]}}){
#     
#     foreach my $sample (keys %{$sharedref{$id}}){
#       for (my $j=1; $j<scalar(@aacolnames);$j++){
#         print AAOUT "$sample\_$aacolnames[$j]\t";
#       }
#     }
#   }
#   #foreach protein, I need the different chromosome ids that correspond to it
#   #from these different chromosome ids, I can get the different samples
#   #from these samples I can check whether the site corresponds to  
#   print AAOUT "\n";
# 
#   my (%alnpos, %gapcnt,%newaatable);
#   foreach my $gene (keys %refseq){
#     foreach my $alnpos (sort {$a<=>$b} keys %{$refseq{$gene}}){
#       print OUT "$gene\t$alnpos\t";
#       foreach my $id (keys %{$genes{$gene}}){    
#         print OUT "$id\t";
#         my $site=$refseq{$gene}{$alnpos}{$id}{"site"};
#         my $base=$refseq{$gene}{$alnpos}{$id}{"base"};
#         
#         
#         #print "$nuctable{$id}{$site}{$sample}{$colnames[$j]}\t"; 
#         if ($refseq{$gene}{$alnpos}{$id}{"base"}=~/-/){
#           #print OUT "NA\tNA\t";
#           foreach my $prot (keys %{$proteins{$id}}){
#             $gapcnt{$prot}{$id}++;
#             #$alnpos{$prot}{$id}{$alnpos}="NA";
#           
# 	       my $refsite = $alnpos-$gapcnt{$prot}{$id};
#            if (keys %{$aatable{$id}{$prot}{$refsite}}){
#              #$alnpos{$prot}{$id}{$alnpos}=$alnpos-$gapcnt{$prot}{$id};
#              foreach my $sample (keys %{$aatable{$id}{$prot}{$refsite}}){
#                print "$prot $id $alnpos $sample $refsite\n";
#                for (my $j=1; $j<scalar(@aacolnames);$j++){
#                  $newaatable{$prot}{$id}{$alnpos}{$sample}{$aacolnames[$j]}=$aatable{$id}{$prot}{$refsite}{$sample}{$aacolnames[$j]};
#                  
#                }
#              }
#            }
# 		  }
#           foreach my $sample (keys %{$sharedref{$id}}){
#             for (my $j=2; $j<scalar(@colnames);$j++){
#               print OUT "NA\t";
#             }
#           }
#         }else{    
#          # print OUT "$site\t$base\t";
#          foreach my $prot (keys %{$proteins{$id}}){ 
#            #perhaps only add it in the hash if it is the first codon position of an aa (i.e. in the coding region) 
#            # listed all files that have same chr, protein, refsite 
# 	       #$aatable{$values[1]}{$values[2]}{$values[5]}{$file}{$aacolnames[$j]}=$values[$j];#Chr Protein RefSite
# 	       my $refsite = $alnpos-$gapcnt{$prot}{$id};
#            if (keys %{$aatable{$id}{$prot}{$refsite}}){
#              #$alnpos{$prot}{$id}{$alnpos}=$alnpos-$gapcnt{$prot}{$id};
#              foreach my $sample (keys %{$aatable{$id}{$prot}{$refsite}}){
#                for (my $j=1; $j<scalar(@aacolnames);$j++){
#                  $newaatable{$prot}{$alnpos}{$id}{$sample}{$aacolnames[$j]}=$aatable{$id}{$prot}{$refsite}{$sample}{$aacolnames[$j]};
#                }
#              }
#            }
#          }
#           foreach my $sample (keys %{$sharedref{$id}}){
#             for (my $j=2; $j<scalar(@colnames);$j++){
#               if ($nuctable{$id}{$site}{$sample}{$colnames[$j]}=~/./){
#                 print OUT "$nuctable{$id}{$site}{$sample}{$colnames[$j]}\t";
#               }else{
#                 print OUT "NA\t";
#               }
#             }
#           } 
#         }
#       }
#       print OUT "\n";
#       
#     }
#   }
#   foreach my $prot (keys %newaatable){
#     foreach my $alnpos (sort {$a<=>$b} keys %{$newaatable{$prot}}){
#       foreach my $id (keys %{$newaatable{$prot}{$alnpos}}){
#         #print AAOUT "$prot\t$id\t$alnpos\t";
#         print AAOUT "$alnpos\t";
#         foreach my $sample (keys %{$newaatable{$prot}{$alnpos}{$id}}){
#           for (my $j=1; $j<scalar(@aacolnames);$j++){
#           
#             #print "$sample ";
#             print AAOUT "$newaatable{$prot}{$alnpos}{$id}{$sample}{$aacolnames[$j]}\t";
#           }
#         }
#       }
#        print AAOUT "\n";
#     }
#   }

}


sub randomize {
  my ($bground,$query) = @_;
  my (%bbasecnt,%qbasecnt,@diffs,@allbases);
  #print "@$bground\n@$query\n";
  my $nbbground=scalar(@$bground);
  my $nbquery=scalar(@$query);
  push (@allbases,@$bground);
  push (@allbases,@$query);
  for (my $r=0; $r<$reps; $r++){
	  my ($bshannon,$qshannon); #to be in the loop for the number of randomisations
	  # this is with replacement (use shuffle indexes for without replacement)
	  for (my $i=0; $i<$nbbground; $i++){
		my $randombase = $allbases[rand @allbases];
		$bbasecnt{$randombase}++;
	  }
	  for (my $i=0; $i<$nbquery; $i++){
		my $randombase = $allbases[rand @allbases];
		$qbasecnt{$randombase}++;
	  }
	  if ($nbbground>0){
		foreach my $nuc (keys %bbasecnt){
		  my $nucnt=$bbasecnt{$nuc};
		  my $p = $nucnt / $nbbground;
		  if($nucnt > 0){
			$bshannon += -$p*log($p);#natural log, i.e. log base e
		  }
		}   
	  }
	  if ($nbquery>0){
		foreach my $nuc (keys %qbasecnt){
		  my $nucnt=$qbasecnt{$nuc};
		  my $p = $nucnt / $nbquery;
		  if($nucnt > 0){
			$qshannon += -$p*log($p);#natural log, i.e. log base e
		  }
		}   
	  }
	  my $dif=$bshannon-$qshannon;
	  #print "$r\t$dif\n";
	  push(@diffs,$dif);
  }
  return(@diffs);
}

  
  