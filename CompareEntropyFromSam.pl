#!/usr/bin/perl

# use this to obtain the number of each base at each site from a bam file
# you need to give the bam file and the reference fasta
#it then randomizes each site to determine whether there is a significant difference in entropy based
# on the proportion of each bases found in total (in both sets) at that site 
# if you want to randomize the whole read you need to use SamEntropyRandomisationFast.pl

use strict;
use Getopt::Long; 
use Bio::DB::Sam;

my ($bam, $ref,$help, $out,%basefreq,$i);

&GetOptions(
	    'bam:s'      => \$bam,#bam file (binary sam)
	    'ref:s'  =>  \$ref, #reference file in fasta  
        "out=s" => \$out,
           );

my $reps=1000;

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
	 print "@targets\n";
	foreach my $target (@targets){
	 print "$target\n";
	 my @alignments = $sam->get_features_by_location(-seq_id => $target);
	 for my $a (@alignments) {
	
		# where does the alignment start in the reference sequence
		my $seqid  = $a->seq_id;
		my $start  = $a->start;
		my $end    = $a->end;
		my $strand = $a->strand;
		my $cigar  = $a->cigar_str;
		my $paired = $a->get_tag_values('PAIRED');
	
		# where does the alignment start in the query sequence
		my $query_start = $a->query->start;     
		my $query_end   = $a->query->end;
	
		my $ref_dna   = $a->dna;        # reference sequence bases
		my $query_dna = $a->query->dna; # query sequence bases
		#print "Ref $start Cigar $cigar $ref_dna $query_dna\n";
		if ($cigar=~/I/){
		  $ins_cnt++;
		}
		if ($cigar!~/I|D|N|S|H|P|=|X/){
		  #print "$cigar\n";
		  $nocigs_cnt++;
		  my @bases=split(//,$query_dna);
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
		  foreach my $nuc (keys %{$basefreq{$bam}{$gene}{$site}}){
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
print ENTROPY "Gene\tSite\tBackgroundEntropy\tS1Total\tS1Acnt\tS1Ccnt\tS1Tcnt\tS1Gcnt\tS1entropybe\tS1entropyb2\tS1simpson\tQueryEntropy\tS2Total\tS2Acnt\tS2Ccnt\tS2Tcnt\tS2Gcnt\tS2entropybe\tS2entropyb2\tS1simpson\tHdiff(Hb-Hq)\tNbRandH>Hdiff\tHighestRandH\tP-value\n";

foreach my $gene (keys %{$shannon{$bam[0]}}){
  foreach my $site (sort {$a<=>$b} keys %{$shannon{$bam[0]}{$gene}}){
    my (@bbases,@qbases);
    my $difshannon=$shannon{$bam[0]}{$gene}{$site}-$shannon{$bam[1]}{$gene}{$site};
    @bbases = (@bbases, ('A') x $basefreq{$bam[0]}{$gene}{$site}{"A"},('C') x $basefreq{$bam[0]}{$gene}{$site}{"C"},('T') x $basefreq{$bam[0]}{$gene}{$site}{"T"},('G') x $basefreq{$bam[0]}{$gene}{$site}{"G"});
    @qbases = (@qbases, ('A') x $basefreq{$bam[1]}{$gene}{$site}{"A"},('C') x $basefreq{$bam[1]}{$gene}{$site}{"C"},('T') x $basefreq{$bam[1]}{$gene}{$site}{"T"},('G') x $basefreq{$bam[1]}{$gene}{$site}{"G"});
	
	my $btotal = $basefreq{$bam[0]}{$gene}{$site}{"A"} + $basefreq{$bam[0]}{$gene}{$site}{"C"} + $basefreq{$bam[0]}{$gene}{$site}{"T"} + $basefreq{$bam[0]}{$gene}{$site}{"G"};
	my $qtotal = $basefreq{$bam[1]}{$gene}{$site}{"A"} + $basefreq{$bam[1]}{$gene}{$site}{"C"} + $basefreq{$bam[1]}{$gene}{$site}{"T"} + $basefreq{$bam[1]}{$gene}{$site}{"G"};
	#print "$gene $site $difshannon\n $bam[0] @bbases\n$bam[1] @qbases\n";
	my @randentropies=randomize(\@bbases,\@qbases);
	#print "@randentropies\n";
	my $max = (sort { $b <=> $a } @randentropies)[0];
	my $randcnt=0;
	foreach my $rand (@randentropies){
	  if ($rand>$difshannon){
	    $randcnt++;
	  }
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

