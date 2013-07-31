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

my ($files,$out,$refalign,%ref,%nuccnt,%refseq,%genes,%aatable);
#my $reps;
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
    open(OUT,">$out\_entropy.txt")||die "Can't open $out\_entropy.txt\n";
    open(AAOUT,">$out\_AA.txt")||die "Can't open $out\_AA.txt\n";
    print "Will NOT be conducting the randomisation or nucletoide frequency test\n";
    open(OUT,">$out\_multi.txt")||die "Can't open $out\_multi.txt\n";
    print "Will conduct entropy randomisation base on nucleotide counts\n";
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
          print "$nuctable{$values[1]}{$values[2]}{$file}{$colnames[$i]} $values[1] $values[2] $file $colnames[$i]\n";
        }
      }
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
    print "will conduct entropy randomisation base on nucleotide counts\n";
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
  my (%nuctable);
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

  
  