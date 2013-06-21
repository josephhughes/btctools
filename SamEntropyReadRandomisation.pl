#!/usr/bin/perl

# use this to obtain the number of each base at each site from a bam file
# you need to give the bam file and the reference fasta
# this randomises the whole read rather than each site as in CompareEntropyFromSam.pl


use strict;
use Getopt::Long; 
use Bio::DB::Sam;

my ($bam, $ref,$help, $out,%basefreq,$i,%allreads,%nbreads,%readcnt,%totalreads);

&GetOptions(
	    'bam:s'      => \$bam,#bam file (binary sam)
	    'ref:s'  =>  \$ref, #reference file in fasta  
        "out=s" => \$out,
           );

my $reps=1000;
print "REPLICATES: $reps\n";
if (($help)&&!($help)||!($bam)||!($ref)){
 print "Usage : perl CompareEntropyFromSam.pl -bam S1_refHPAI_cons_stampy.bam,S2_refHPAI_cons_stampy.bam -ref refHPAI_cons.fa -out S1_basefreq.txt \n";
 print " -bam <txt> - the input two files in bam format\n";
 print " -ref <txt>  - the reference fasta file\n";
 print " -out - the output in text-tab delimited\n";
 print " -help        - Get this help\n";
 exit();
 }
 
my @bam=split(/,/,$bam);
print "$bam[0]\t$bam[1]\n";
foreach my $bam(@bam){
     # high level API
	 my $sam = Bio::DB::Sam->new(-bam  => $bam,
								 -fasta=> $ref);
	my $ins_cnt=0;
	my $nocigs_cnt=0;
	my @targets    = $sam->seq_ids;
	# print "@targets\n";
	foreach my $target (@targets){
	 # print "$target\n";
	 my @alignments = $sam->get_features_by_location(-seq_id => $target);
	 for my $a (@alignments) {
	    #print "$a";
	    my $cigar  = $a->cigar_str;
		# where does the alignment start in the reference sequence
		my $seqid  = $a->seq_id;
		my $start  = $a->start;
		my $query_dna = $a->query->dna; # query sequence bases
		#print "Ref $start Cigar $cigar $ref_dna $query_dna\n";
		if ($cigar=~/I/){
		  $ins_cnt++;
		}
		if ($cigar!~/I|D|N|S|H|P|=|X/){
		  #print "$cigar\n";
		  $nocigs_cnt++;
		  my @bases=split(//,$query_dna);
		  $readcnt{$target}++;
		  $nbreads{$bam}{$target}++;
		  $allreads{$target}{$readcnt{$target}}="$start\t$query_dna";
		  for ($i=0; $i<scalar(@bases); $i++){
			# chromosome site nuc 
			my $site=$i+$start;
			$basefreq{$bam}{$target}{$site}{$bases[$i]}++;
		  }
		}
		my @scores    = $a->qscore;     # per-base quality scores
		my $match_qual= $a->qual;       # quality of the match
	 }
	}
	print "$bam:\nNumber of inserts $ins_cnt\nNumber of sequence with NO cigars $nocigs_cnt\n";
}

my %shannon;#Shannon-Wiener Diversity Index (also known as Shannon's diversity index, the Shannon-Wiener index, the Shannon-Weaver index and the Shannon entropy)
my %simpson;#the Simpson diversity index D = probability that two sequences taken at random from the dataset represent the same sequence
my %entropy;#the shannon entropy as calculated at Los Alamos


my @nuc=qw/A C T G/;
open (OUT, ">$out")||die "can't open $out\n";
print OUT "Sample\tChr\tPosition\tTotal\tAcnt\tCcnt\tTcnt\tGcnt\tentropy(base e)\tentropy(base 2)\tsimpson\n";
foreach my $bam (@bam){
  foreach my $gene (keys %{$basefreq{$bam}}){
    foreach my $site (sort {$a<=>$b} keys %{$basefreq{$bam}{$gene}}){
	  my $cntA = $basefreq{$bam}{$gene}{$site}{"A"};
	  my $cntC = $basefreq{$bam}{$gene}{$site}{"C"};
	  my $cntT = $basefreq{$bam}{$gene}{$site}{"T"};
	  my $cntG = $basefreq{$bam}{$gene}{$site}{"G"};
	  my $total = $cntA + $cntT + $cntC + $cntG;
	  if ($total>0){
	    foreach my $nuc (@nuc){
		  my $nucnt=$basefreq{$bam}{$gene}{$site}{$nuc};
		  my $p = $nucnt / $total;
		  if($nucnt > 0){
		    $shannon{$bam}{$gene}{$site} += -$p*log($p);#natural log, i.e. log base e
			$entropy{$bam}{$gene}{$site} -= ((log_base(2,$p))*$p);
		  }
		  $simpson{$bam}{$gene}{$site} += $p**2;#to the power of two
		}   
	  }   
	  print OUT "$bam\t$gene\t$site\t$total\t$cntA\t$cntC\t$cntT\t$cntG\t$shannon{$bam}{$gene}{$site}\t$entropy{$bam}{$gene}{$site}\t$simpson{$bam}{$gene}{$site}\n";
	}
  }
}

open (ENTROPY, ">Compare$out")||die "can't open entropy_results.txt\n";
print ENTROPY "Gene\tSite\tBackgroundEntropy\tS1Total\tS1Acnt\tS1Ccnt\tS1Tcnt\tS1Gcnt\tS1entropybe\tS1entropyb2\tS1simpson\tQueryEntropy\tS2Total\tS2Acnt\tS2Ccnt\tS2Tcnt\tS2Gcnt\tS2entropybe\tS2entropyb2\tS1simpson\tHdiff(Hb-Hq)\tNbRandH>Hdiff\tHighestRandH\tNbRandH>Hdiff\tP-value\n";

# with replacement select n many reads corresponding to the number of reads in each set
my (%randshannon, %randentropy, %randsimpson,%randif);

for (my $i=0; $i<$reps; $i++){
 my %randfreq=();
 foreach my $bam (keys %nbreads){
  foreach my $gene (keys %{$nbreads{$bam}}){
    for (my $j=0; $j<$nbreads{$bam}{$gene}; $j++){
      my $random = int (rand($readcnt{$gene})+1);
      # print "gene $gene from $bam rep $i random number $j is $random\n";
      my @seqinfo=split(/\t/,$allreads{$gene}{$random});
      # print "INFORMATION $seqinfo[1]\n";
      my @bases=split(//,$seqinfo[1]);
	  for (my $k=0; $k<scalar(@bases); $k++){
		  # chromosome site nuc 
		  my $site=$k+$seqinfo[0];
		  # print "Site $site $bases[$k]\n";
		  $randfreq{$bam}{$gene}{$site}{$bases[$k]}++;
	  }
	}
   }
  }
  foreach my $gene (keys %{$randfreq{$bam[0]}}){
	 #print "$gene ";
	  foreach my $site (sort {$a<=>$b} keys %{$randfreq{$bam[0]}{$gene}}){
	    # print "Site for randfreq $site \n";
		my $bcntA = $randfreq{$bam[0]}{$gene}{$site}{"A"};
		my $bcntC = $randfreq{$bam[0]}{$gene}{$site}{"C"};
		my $bcntT = $randfreq{$bam[0]}{$gene}{$site}{"T"};
		my $bcntG = $randfreq{$bam[0]}{$gene}{$site}{"G"};
		my $btotal = $bcntA + $bcntT + $bcntC + $bcntG;
		# print "Rep $i Number of reads $nbreads{$bam[0]}{$gene} \n total $btotal A:$bcntA C:$bcntC T:$bcntT G:$bcntG\n";
		my $rtotal=$basefreq{$bam[0]}{$gene}{$site}{"A"}+$basefreq{$bam[0]}{$gene}{$site}{"C"}+$basefreq{$bam[0]}{$gene}{$site}{"T"}+$basefreq{$bam[0]}{$gene}{$site}{"G"};
		# print "real $rtotal A:".$basefreq{$bam}{$gene}{$site}{"A"}." C:".$basefreq{$bam}{$gene}{$site}{"C"}." T:".$basefreq{$bam}{$gene}{$site}{"T"}." G:".$basefreq{$bam}{$gene}{$site}{"G"}."\n";
		my $qcntA = $randfreq{$bam[1]}{$gene}{$site}{"A"};
		my $qcntC = $randfreq{$bam[1]}{$gene}{$site}{"C"};
		my $qcntT = $randfreq{$bam[1]}{$gene}{$site}{"T"};
		my $qcntG = $randfreq{$bam[1]}{$gene}{$site}{"G"};
		my $qtotal = $qcntA + $qcntT + $qcntC + $qcntG;
		# print "Rep $i Number of reads $nbreads{$bam[1]}{$gene} \n total $qtotal A:$qcntA C:$qcntC T:$qcntT G:$qcntG\n";
		my $rtotal=$basefreq{$bam[1]}{$gene}{$site}{"A"}+$basefreq{$bam[1]}{$gene}{$site}{"C"}+$basefreq{$bam[1]}{$gene}{$site}{"T"}+$basefreq{$bam[1]}{$gene}{$site}{"G"};
		if ($btotal>0){
		  foreach my $nuc (@nuc){
			my $nucnt=$randfreq{$bam[0]}{$gene}{$site}{$nuc};
			my $p = $nucnt / $btotal;
			# print "frequency $p\n";
			if($nucnt > 0){
			  $randshannon{$bam[0]}{$gene}{$site}{$i} += -$p*log($p);#natural log, i.e. log base e
			  $randentropy{$bam[0]}{$gene}{$site}{$i} -= ((log_base(2,$p))*$p);
			  # print "$randshannon{$bam[0]}{$gene}{$site}{$i}\n";
			}
		    $randsimpson{$bam[0]}{$gene}{$site}{$i} += $p**2;#to the power of two
		  }   
		}  
		if ($qtotal>0){
		  foreach my $nuc (@nuc){
			my $nucnt=$randfreq{$bam[1]}{$gene}{$site}{$nuc};
			my $p = $nucnt / $qtotal;
			# print "frequency $p\n";
			if($nucnt > 0){
			  $randshannon{$bam[1]}{$gene}{$site}{$i} += -$p*log($p);#natural log, i.e. log base e
			  $randentropy{$bam[1]}{$gene}{$site}{$i} -= ((log_base(2,$p))*$p);
			  # print "$randshannon{$bam}{$gene}{$site}{$i}\n";
			}
		    $randsimpson{$bam[1]}{$gene}{$site}{$i} += $p**2;#to the power of two
		  }   
		}   
        my $randif=$randshannon{$bam[0]}{$gene}{$site}{$i}-$randshannon{$bam[1]}{$gene}{$site}{$i};
        my $difshannon=$shannon{$bam[0]}{$gene}{$site}-$shannon{$bam[1]}{$gene}{$site};
        # print "Real dif $difshannon\n";
        # print "Random dif for rep $i $randif\n";
	}
  }
}
foreach my $gene (keys %{$shannon{$bam[0]}}){
  foreach my $site (sort {$a<=>$b} keys %{$shannon{$bam[0]}{$gene}}){
    my $difshannon=$shannon{$bam[0]}{$gene}{$site}-$shannon{$bam[1]}{$gene}{$site};
	my $btotal = $basefreq{$bam[0]}{$gene}{$site}{"A"} + $basefreq{$bam[0]}{$gene}{$site}{"C"} + $basefreq{$bam[0]}{$gene}{$site}{"T"} + $basefreq{$bam[0]}{$gene}{$site}{"G"};
	my $qtotal = $basefreq{$bam[1]}{$gene}{$site}{"A"} + $basefreq{$bam[1]}{$gene}{$site}{"C"} + $basefreq{$bam[1]}{$gene}{$site}{"T"} + $basefreq{$bam[1]}{$gene}{$site}{"G"};

	my @randentropies;
	for (my $i=0; $i<$reps; $i++){
	  my $dif=$randshannon{$bam[0]}{$gene}{$site}{$i}-$randshannon{$bam[1]}{$gene}{$site}{$i};
	  push (@randentropies,$dif);
	}
	# print "@randentropies\n";
	my $max = (sort { $b <=> $a } @randentropies)[0];
	my $randcnt=0;
	foreach my $rand (@randentropies){
	  if ($rand>$difshannon){
	    $randcnt++;
	  }
	  #print "$max\t$difshannon\t$rand\t$randcnt\n";
	}
	
	my $pvalue=$randcnt/$reps;
	print ENTROPY "$gene\t$site\t$shannon{$bam[0]}{$gene}{$site}\t$btotal\t";
	print ENTROPY $basefreq{$bam[0]}{$gene}{$site}{"A"}."\t".$basefreq{$bam[0]}{$gene}{$site}{"C"}."\t".$basefreq{$bam[0]}{$gene}{$site}{"T"}."\t".$basefreq{$bam[0]}{$gene}{$site}{"G"}."\t";
	print ENTROPY "$shannon{$bam[0]}{$gene}{$site}\t$entropy{$bam[0]}{$gene}{$site}\t$simpson{$bam[0]}{$gene}{$site}\t";
	print ENTROPY "$shannon{$bam[1]}{$gene}{$site}\t$qtotal\t";
	print ENTROPY $basefreq{$bam[1]}{$gene}{$site}{"A"}."\t".$basefreq{$bam[1]}{$gene}{$site}{"C"}."\t".$basefreq{$bam[1]}{$gene}{$site}{"T"}."\t".$basefreq{$bam[1]}{$gene}{$site}{"G"}."\t";
	print ENTROPY "$shannon{$bam[1]}{$gene}{$site}\t$entropy{$bam[1]}{$gene}{$site}\t$simpson{$bam[1]}{$gene}{$site}\t";
	print ENTROPY "$difshannon\t$max\t$randcnt\t$pvalue\n";
  }
}



sub log_base {
    my ($base, $value) = @_;
    return log($value)/log($base);
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
	  #print "$dif\n";
	  push(@diffs,$dif);
  }
  return(@diffs);
}

