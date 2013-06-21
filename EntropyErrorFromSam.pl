#!/usr/bin/perl

# use this to obtain the number of each base at each site from a bam file
# you need to give the bam file and the reference fasta
# it also calculate the avergae per site entropy for the sample
# and if -window is specified calculates a sliding window of the average entropy

use strict;
use Getopt::Long; 
use Bio::DB::Sam;
use Math::CDF;

my ($bam, $ref,$help, $out,%basefreq,$i,$stub,$window,%cumulqual);

&GetOptions(
	    'bam:s'      => \$bam,#bam file (binary sam)
	    'ref:s'  =>  \$ref, #reference file in fasta  
	    'window=i'  => \$window, #window size
        "stub=s" => \$stub,
           );


if (($help)&&!($help)||!($bam)||!($ref)){
 print "Usage : perl EntropyFromSam.pl -bam S1_refHPAI_cons_stampy.bam -ref refHPAI_cons.fa -out S1_basefreq.txt \n";
 print " -bam <txt> - the input a file in bam format\n";
 print " -ref <txt>  - the reference fasta file\n";
 print " -stub - the output in text-tab delimited\n";
 print " -window - the sliding window size\n";
 print " -help        - Get this help\n";
 exit();
 }

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
	my @scores    = $a->qscore;     # per-base quality scores
	my $match_qual= $a->qual;       # quality of the match

	#print "Ref $start Cigar $cigar $ref_dna $query_dna\n";
	if ($cigar=~/I/){
	  $ins_cnt++;
	}
	if ($cigar!~/I|D|N|S|H|P|X/){#|= if there is an = it means exact match so want to parse
	  #print "$cigar\n";
	  $nocigs_cnt++;
	  my @bases=split(//,$query_dna);
	  my $cumP=0;
	  for ($i=0; $i<scalar(@bases); $i++){
		# chromosome site nuc 
		my $site=$i+$start;
		$basefreq{$bam}{$target}{$site}{$bases[$i]}++;
		my $P = 10**(-$scores[$i]/10);
		if ($cumulqual{$bam}{$target}{$site}){
		  $cumulqual{$bam}{$target}{$site}=$cumulqual{$bam}{$target}{$site}+$P;
		}else{
		  $cumulqual{$bam}{$target}{$site}=$P;
		}
		# create a hash for the read information (position of mismatches and motifs upstream and downstream) => Jo's 
		
		# create a hash for the codon and aa information => Joseph
		
	  }
	}
 }
}
print "$bam:\nNumber of inserts $ins_cnt\nNumber of sequence with NO cigars $nocigs_cnt\n";

my %shannon;#Shannon-Wiener Diversity Index (also known as Shannon's diversity index, the Shannon-Wiener index, the Shannon-Weaver index and the Shannon entropy)
my %simpson;#the Simpson diversity index D = probability that two sequences taken at random from the dataset represent the same sequence
my %entropy;#the shannon entropy as calculated at Los Alamos
my %rawshannon;# this calculates the shannon entropy taking into account all variation, i.e. not removing sites below the illumina threshold

my @nuc=qw/A C T G/;
open (OUT, ">$stub\_entropy.txt")||die "can't open $stub\_entropy.txt\n";
# add the reference base in the table
print OUT "Sample\tChr\tPosition\tTotal\tAvQual\tAcnt\tApval\tCcnt\tCpval\tTcnt\tTpval\tGcnt\tGpval\tentropy(base e)\tentropy(base 2)\tsimpson\n";
foreach my $gene (keys %{$basefreq{$bam}}){
  my $nbsites=0;
  my $sumentropy=0;
  my $rawsumentropy=0;
  foreach my $site (sort {$a<=>$b} keys %{$basefreq{$bam}{$gene}}){
    my ($cntA,$cntC,$cntT,$cntG);
	if ($basefreq{$bam}{$gene}{$site}{"A"}){ $cntA = $basefreq{$bam}{$gene}{$site}{"A"}}else{$cntA=0}
	if ($basefreq{$bam}{$gene}{$site}{"C"}){ $cntC = $basefreq{$bam}{$gene}{$site}{"C"}}else{$cntC=0}
	if ($basefreq{$bam}{$gene}{$site}{"T"}){ $cntT = $basefreq{$bam}{$gene}{$site}{"T"}}else{$cntT=0}
	if ($basefreq{$bam}{$gene}{$site}{"G"}){ $cntG = $basefreq{$bam}{$gene}{$site}{"G"}}else{$cntG=0}
	
	my $total = $cntA + $cntT + $cntC + $cntG;
	my $average_p=$cumulqual{$bam}{$gene}{$site}/$total;
	my $p = $average_p/3;
	my %prob;
    $prob{"A"} = 1 - (&Math::CDF::pbinom(($cntA-1), $total, $p));# need to double check with Richard about the -1
    $prob{"C"} = 1 - (&Math::CDF::pbinom(($cntC-1), $total, $p));
    $prob{"T"} = 1 - (&Math::CDF::pbinom(($cntT-1), $total, $p));
    $prob{"G"} = 1 - (&Math::CDF::pbinom(($cntG-1), $total, $p));
 	
	if ($total>0){
	  foreach my $nuc (keys %{$basefreq{$bam}{$gene}{$site}}){
		my $nucnt=$basefreq{$bam}{$gene}{$site}{$nuc};
		my $p = $nucnt / $total;
		if($nucnt > 0){
		  $rawshannon{$bam}{$gene}{$site} += -$p*log($p);#natural log, i.e. log base e
		  $entropy{$bam}{$gene}{$site} -= ((log_base(2,$p))*$p);
		}
		$simpson{$bam}{$gene}{$site} += $p**2;#to the power of two
		if($nucnt > 0 && $prob{$nuc}<0.05){
		  $shannon{$bam}{$gene}{$site} += -$p*log($p);#natural log, i.e. log base e
		}
	  }   
	}   
	$nbsites++;
	$rawsumentropy=$rawsumentropy+$rawshannon{$bam}{$gene}{$site};
	$sumentropy=$sumentropy+$shannon{$bam}{$gene}{$site};
	print OUT "$bam\t$gene\t$site\t$total\t$average_p\t$cntA\t".$prob{"A"}."\t$cntC\t".$prob{"C"}."\t$cntT\t".$prob{"T"}."\t$cntG\t".$prob{"G"}."\t";
	print OUT "$rawshannon{$bam}{$gene}{$site}\t$entropy{$bam}{$gene}{$site}\t$simpson{$bam}{$gene}{$site}\t$shannon{$bam}{$gene}{$site}\n";
  }
  print "Gene $gene Average raw entropy = ".$rawsumentropy/$nbsites."\n";
  print "Gene $gene Average entropy = ".$sumentropy/$nbsites."\n";
}
# Jo's code for output of the read and motif table



# Jo's loop for the aa mutations and position of mismatches in the codon



if ($window){
	open (WINDOW, ">$stub\_$window\.txt")||die "can't open $stub\_$window\.txt\n";
	print WINDOW "Gene\tPosition\tRaw_Average_entropy\tAverage_Entropy\n";
	foreach my $gene (keys %{$rawshannon{$bam}}){
	  print "$gene\n";
	  foreach my $site (sort {$a<=>$b} keys %{$rawshannon{$bam}{$gene}}){
		my $rawsum=0;
		my $sum=0;
		if ($rawshannon{$bam}{$gene}{$site+$window}){
		  for ($i=0; $i<$window; $i++){
			$rawsum=$rawsum+$rawshannon{$bam}{$gene}{$site+$i};
			$sum=$sum+$shannon{$bam}{$gene}{$site+$i};
		  }
		  print WINDOW "$gene\t$site\t".$rawsum/$window."\t".$sum/$window."\n";
		}
	  }
	
	}
	close(WINDOW);
	system("Rscript PlotAvEntropyWindowPerGene.R $stub\_$window\.txt");
}


sub log_base {
    my ($base, $value) = @_;
    return log($value)/log($base);
}
