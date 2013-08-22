#!/usr/bin/perl

# use this to obtain the number of each base at each site from a bam file
# you need to give the bam file and the reference fasta
# it also calculate the average per site entropy for the sample
# 
# produces a table with Sample\tChr\tPosition\tRefBase\tCoverage\t
# AvQual\tAcnt\tApval\tCcnt\tCpval\tTcnt\tTpval\tGcnt\tGpval\tentropy(base e)\tentropy(base 2)\tsimpson\tNonRefCnt
# CntTv\tCntTs\tCodonPos\tCntNonSyn\tCntSyn\tOrderATCG
# example command: ./btcutils -bam SamTestFiles/S1_refHPAI_cons_stampy.bam -ref SamTestFiles/refHPAI_cons.fa -stub out
# perl btcutils.pl -bam SamTestFiles/S1_refHPAI_cons_stampy.bam -ref SamTestFiles/refHPAI_cons.fa -stub out
use strict;
use Getopt::Long; 
use Bio::DB::Sam;
use Math::CDF;
use Bio::SeqIO;

# global variables
my ($bam, $ref,$help, $out,%basefreq,%refbase,$i,%cumulqual,$orfs,%readinfo,%IUPAC,%c2p,%aafreq,%aaorder,%readmis);
# addition of avqual threshold, coverage threshold, pval threshold => need to re-think this
# %readmis used to count the number of mismatches per reads
# my $tpval=0;
# my $tcov=0;
# my $tqual=0;
my $stub="output";

&GetOptions(
	    'bam:s'  => \$bam,#bam file (binary sam)
	    'ref:s'  => \$ref, #reference file in fasta  
	    'orfs:s' => \$orfs, #start stop position for each gene labelled the same way as the ref file, keep in mind that a gene may code for multiple proteins
        "stub:s" => \$stub,
           );

if (($help)&&!($help)||!($bam)||!($ref)){
 print "Usage : perl btcutils.pl -bam input.bam -ref ref.fa -out stub \n";
 print " -bam <txt>  - the input a file in bam format\n";
 print " -ref <txt>  - the reference fasta file\n";
 print " -orfs <txt> - text tab-delimited file with Protein,Beginning,End,Reference(Chr) of the coding sequence [optional: only if you want information about dN/dS and aa frequencies etc...]\n";
 print " -stub <txt> - the output in text-tab delimited\n";
 print " -help        - Get this help\n";
 exit();
}
# open an output file with the motif upstream and downstream of a mismatch
open(MOTIF,">$stub\_motif.fa")||die "Can't open output $stub\_motif.fa\n";
open(LOG,">$stub\_log.txt")||die "Can't open output $stub\_log.txt\n";

# get a hash of the reference sequences with gene, site, refbase
my %refseq;
my $seq_in  = Bio::SeqIO->new(-format => 'fasta',-file   => $ref);
while( my $seq = $seq_in->next_seq() ) {
  my $id=$seq->display_id();
  my $seq_str=$seq->seq();
  my @refbases=split(//,$seq_str);
  for (my $i=0; $i<scalar(@refbases); $i++){
    my $site=1+$i;
    $refseq{$id}{$site}=$refbases[$i];
  }
}

# high level API
my $sam = Bio::DB::Sam->new(-bam  => $bam,
							 -fasta=> $ref);
my @targets    = $sam->seq_ids;
my $ins_cnt=0;
my $nocigs_cnt=0;
my $Ncnt=0;
my %codreg;
if ($orfs){
    open (CODING,"<$orfs")||die "Can't open $orfs\n";
    my $firstLine = <CODING>; 
    if ($firstLine !~/Protein\tBeg\tEnd\tReference/){
      die "Error: The input file $orfs does not have the proper HEADER:\nProtein	Beg	End	Reference\n";
    }else{
      while(<CODING>){
	chomp($_);
	my @elements=split(/\t/,$_);
	$codreg{$elements[3]}{$elements[0]}{"Beg"}=$elements[1];#$codreg{Chr name}{ProteinName}{"Beg"}
	$codreg{$elements[3]}{$elements[0]}{"End"}=$elements[2];
      }
    }
    #load translation table if start and stop specified
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
    # Ambiguous Amino Acids	3-Letter	1-Letter
    # Asparagine or aspartic acid	Asx	B
    # Glutamine or glutamic acid	Glx	Z
    # Leucine or Isoleucine	Xle	J
    # Unspecified or unknown amino acid	Xaa	X
    #--------------------------#
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
    #$c2p{$_} = "_" for qw(TAA TAG TAR TGA);
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
    foreach my $target (@targets){
      if (!keys %{$codreg{$target}}){
        my @orfIDs=keys %codreg;
	    die "Error: The reference identifer $target in the bam file does not match the reference identifier in your $orfs:\n @orfIDs\n";
      }
    }
}
my $increment;
foreach my $target (@targets){
  print "\n$target\n";
  $increment=0;
  my @alignments = $sam->get_features_by_location(-seq_id => $target);
  my $refobj = $sam->segment(-seq_id => $target);
  my $refseq = $refobj->dna;
  my $total = scalar(@alignments);
  my $progress=0;
  for my $a (@alignments) {
	# where does the alignment start in the reference sequence
	$progress++;
	progress_bar( $progress, $total, 25, '=' );
	my $seqid  = $a->seq_id;
	my $name   = $a->display_name;
	my $start  = $a->start; #position in the reference
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
    my $mismatches=0;
	#print "Ref $start Cigar $cigar $ref_dna $query_dna\n";
	if ($cigar=~/I|D|N|S|H|P|X/){
	  $ins_cnt++;
	}
	if ($cigar!~/I|D|N|S|H|P|X/){#|= if there is an = it means exact match so want to parse
	  #print "$cigar\n";
	  $nocigs_cnt++;
	  if ($query_dna=~/N/){
	    $Ncnt++;
	  }else{
	  my @bases=split(//,$query_dna);
	  my @refbases=split(//,$ref_dna);
	  my $cumP=0;
	  my $matches=0;
	  for ($i=0; $i<scalar(@bases); $i++){
		# chromosome site nuc 
		my $site=$i+$start;
		my $readpos=$i+1;
		$basefreq{$bam}{$target}{$site}{$bases[$i]}++;
		$refbase{$bam}{$target}{$site}=$refbases[$i];
		my $P = 10**(-$scores[$i]/10);
		if ($cumulqual{$bam}{$target}{$site}){
		  $cumulqual{$bam}{$target}{$site}=$cumulqual{$bam}{$target}{$site}+$P;
		}else{
		  $cumulqual{$bam}{$target}{$site}=$P;
		}
		# create a hash for the read information (position of mismatches relative to the start position of a read and motifs upstream and downstream of a mismatch)
		$readinfo{$readpos}{"AvQual"}=$readinfo{$readpos}{"AvQual"}+$P;
		if ($refbases[$i]=~/$bases[$i]/i){
		  $readinfo{$readpos}{"CntRef"}++;
		  $readinfo{$readpos}{"AvQualRef"}=$readinfo{$readpos}{"AvQualRef"}+$P;
		}elsif ($refbases[$i]!~/$bases[$i]/i){ 
		  $mismatches++;
		  $readinfo{$readpos}{"CntNonRef"}++;
		  $readinfo{$readpos}{"AvQualNonRef"}=$readinfo{$readpos}{"AvQualNonRef"}+$P;
		  # motif on either side of mismatch
		  my ($motif_start,$motif_end,$endfromref,$startfromref);
		  if (length($query_dna)<=($i+10)){
		    $motif_end=length($query_dna);
		    $endfromref=$i+10+1-length($query_dna);
		  }elsif (length($query_dna)>($i+10)){
		    $motif_end=$i+10;
		  }
		  if ($i-10<0){
		    $motif_start=0;
		    $startfromref=10-$i;
		  }elsif($i-10>=0){
		    $motif_start=$i-10;
		  }
		  if ($start-1+$i-10>=0 && length($refseq)>=$i+10){
			  my $motif = substr ($query_dna,$motif_start,$motif_end+1-$motif_start);
			  my $begmotif = lc(substr ($refseq,($start-1+$i-10),$startfromref));
			  my $endmotif = lc(substr ($refseq,(($start-1+$i-10+21)-$endfromref),$endfromref));
			  $motif = $begmotif.$motif.$endmotif;
			  # print "Motif $begmotif$motif$endmotif\n";
			  my $refmotif = substr ($refseq,($start-1+$i-10),21);
			  # print "Refer $refmotif\n\n";
			  if (length($motif)<21){
			    #the motif is too close to the beginning or end of the reference sequence
			    #print "Motifend $motif_end EndFromRef $endfromref\t StartMotif $motif_start StartFromRef $startfromref\tlength ".length($query_dna)." position $i Length of reference ".length($refseq)." site $site\n";
			  }elsif (length($motif)==21){
			    print MOTIF ">$name\_$mismatches\n$motif\n";
			  }
		  }
		}
		
		# create a hash for the codon and aa information only if information on start and stop are given
		if ($orfs){
		  #$codreg{Chr name}{ProteinName}{"Beg"}
		  foreach my $prot (keys %{$codreg{$target}}){
		    my $noUTR=($site-$codreg{$target}{$prot}{"Beg"})+1;
		    my ($aasite,$mod);
		    if ($noUTR==1){
              $mod=1;
            }elsif ($noUTR>1){
              $mod = $noUTR % 3;
            }
            # check that the site is in a coding region, that it is the first codon position of the coding region and that the read is long enough for final codon
            if ($site>=$codreg{$target}{$prot}{"Beg"} && $site<=$codreg{$target}{$prot}{"End"} && $mod==1 && $i<(scalar(@bases)-2)){
              #print "Coding region for $prot Site $site and Modular $mod\t";
              my $rcodon=$refbases[$i].$refbases[$i+1].$refbases[$i+2];
              my $qcodon=$bases[$i].$bases[$i+1].$bases[$i+2];
              #ignore AA if 
              my $raa= $c2p{uc($rcodon)};
              my $qaa= $c2p{uc($qcodon)};
              #print "$rcodon\t$qcodon\n";
              #nucsite aasite $rcodon $raa $qcodon $qaa $codonposmis
              if ($noUTR==1){
                $aasite=1;
              }else{
                $aasite = ($noUTR+3-$mod)/3;
                #print "No UTR $noUTR modular $mod $aasite\n";
              }
              # the position of the mutation
              # Sample\tChr\tAAPosition\tRefAA\tRefCodon\tCntNonSyn\tCntSyn\tTopAA\tTopAAcnt\tSndAA\tSndAAcnt\tTrdAA\tTrdAAcnt\tAAcoverage\tNbStop
              $aafreq{$bam}{$target}{$prot}{$aasite}{"AAcoverage"}++;#this will provide the coverage
              $aaorder{$bam}{$target}{$prot}{$aasite}{$qaa}++;
              $aafreq{$bam}{$target}{$prot}{$aasite}{"RefAA"}=$raa;
              $aafreq{$bam}{$target}{$prot}{$aasite}{"RefCodon"}=$rcodon;
              $aafreq{$bam}{$target}{$prot}{$aasite}{"RefSite"}=$site;
              my $aamut=$raa.$aasite.$qaa;
              if (uc($rcodon) ne uc($qcodon)){
				if(uc($raa) eq uc($qaa)){
				    my $mm=mismatch_count($rcodon,$qcodon);
				   	$aafreq{$bam}{$target}{$prot}{$aasite}{"syn"}+=$mm;
				}if(uc($raa) ne uc($qaa)){
				    my $mm=mismatch_count($rcodon,$qcodon);
					$aafreq{$bam}{$target}{$prot}{$aasite}{"nonsyn"}+=$mm;
				}if(uc($qaa)=~/\*/){
				    $aafreq{$bam}{$target}{$prot}{$aasite}{"stop"}++;
				}if($refbases[$i] ne $bases[$i]){
				    $aafreq{$bam}{$target}{$prot}{$aasite}{"firstcodonpos"}++;
				}if($refbases[$i+1] ne $bases[$i+1]){
				    $aafreq{$bam}{$target}{$prot}{$aasite}{"secondcodonpos"}++;
			    }if($refbases[$i+2] ne $bases[$i+2]){
			       $aafreq{$bam}{$target}{$prot}{$aasite}{"thirdcodonpos"}++;
				}
              }
            }
		  }
		}
	  }
	}
	$readmis{$mismatches}++;
	}#close the if no bad cigars
	
 }
}
print LOG "$bam:\nNumber of reads with inserts: $ins_cnt\nNumber of reads with Ns: $Ncnt\nNumber of sequence used: $nocigs_cnt\n";
# print out the number of reads with 1 , 2, 3, 4, etc.. mismatches
print LOG "Frequency of mismatches per read (mismatches: number of reads)\n";
for my $mismatchcnt (sort {$a<=>$b} keys %readmis){
  print LOG "$mismatchcnt: $readmis{$mismatchcnt}\n";
} 

my %shannon;#Shannon-Wiener Diversity Index (also known as Shannon's diversity index, the Shannon-Wiener index, the Shannon-Weaver index and the Shannon entropy)

my @nuc=qw/A C T G/;
open (OUT, ">$stub\_entropy.txt")||die "can't open $stub\_entropy.txt\n";
# add the reference base in the table
# Nucleotide Table:
# Sample\tChr\tPosition\tRefBase\tCoverage\tAvQual\tAcnt\tApval\tCcnt\tCpval\tTcnt\tTpval\tGcnt\tGpval\tentropy(base e)\tNonRefCnt
# CntTv\tCntTs\tOrderATCG

print OUT "Sample\tChr\tPosition\tRefBase\tCoverage\tAvQual\tAcnt\tApval\tCcnt\tCpval\tTcnt\tTpval\tGcnt\tGpval\tentropy(base e)\tNonRefCnt\tCntTv\tCntTs\tOrderOfNucs\n";
#foreach my $gene (keys %{$basefreq{$bam}}){
foreach my $gene (keys %refseq){
  my $nbsites=0;
  my $sumentropy=0;
  my $rawsumentropy=0;
  #foreach my $site (sort {$a<=>$b} keys %{$basefreq{$bam}{$gene}}){
  foreach my $site (sort {$a<=>$b} keys %{$refseq{$gene}}){
    if (keys %{$basefreq{$bam}{$gene}{$site}}){# if there is information in the bam file about this site
    my ($cntA,$cntC,$cntT,$cntG);
	if ($basefreq{$bam}{$gene}{$site}{"A"}){ $cntA = $basefreq{$bam}{$gene}{$site}{"A"}}else{$cntA=0}
	if ($basefreq{$bam}{$gene}{$site}{"C"}){ $cntC = $basefreq{$bam}{$gene}{$site}{"C"}}else{$cntC=0}
	if ($basefreq{$bam}{$gene}{$site}{"T"}){ $cntT = $basefreq{$bam}{$gene}{$site}{"T"}}else{$cntT=0}
	if ($basefreq{$bam}{$gene}{$site}{"G"}){ $cntG = $basefreq{$bam}{$gene}{$site}{"G"}}else{$cntG=0}
	
	my $coverage = $cntA + $cntT + $cntC + $cntG;
	my $average_p=$cumulqual{$bam}{$gene}{$site}/$coverage;
	my $p = $average_p/3;
	my %prob;
    $prob{"A"} = 1 - (&Math::CDF::pbinom(($cntA-1), $coverage, $p));# need to double check with Richard about the -1
    $prob{"C"} = 1 - (&Math::CDF::pbinom(($cntC-1), $coverage, $p));
    $prob{"T"} = 1 - (&Math::CDF::pbinom(($cntT-1), $coverage, $p));
    $prob{"G"} = 1 - (&Math::CDF::pbinom(($cntG-1), $coverage, $p));
 	
	if ($coverage>0){
	  foreach my $nuc (keys %{$basefreq{$bam}{$gene}{$site}}){
		my $nucnt=$basefreq{$bam}{$gene}{$site}{$nuc};
		my $p = $nucnt / $coverage;
		if($nucnt > 0){
		  $shannon{$bam}{$gene}{$site} += -$p*log($p);#natural log, i.e. log base e
		}
	  }   
	}   
	$nbsites++;
	$sumentropy=$sumentropy+$shannon{$bam}{$gene}{$site};
	my $refbase=$refbase{$bam}{$gene}{$site};
	my $nonrefcnt=$coverage-($basefreq{$bam}{$gene}{$site}{$refbase});
	my ($Ts,$Tv) = cntTsTv($refbase,$cntA,$cntT,$cntG,$cntC);
    my $NucOrder="";
    foreach my $nuc (sort { $basefreq{$bam}{$gene}{$site}{$b} <=> $basefreq{$bam}{$gene}{$site}{$a} } keys %{$basefreq{$bam}{$gene}{$site}}) {
        $NucOrder=$NucOrder.$nuc;
    }
    #print "$site $refbase $NucOrder\n";	 
	print OUT "$bam\t$gene\t$site\t$refbase\t$coverage\t$average_p\t$cntA\t".$prob{"A"}."\t$cntC\t".$prob{"C"}."\t$cntT\t".$prob{"T"}."\t$cntG\t".$prob{"G"}."\t";
	print OUT "$shannon{$bam}{$gene}{$site}\t$nonrefcnt\t$Ts\t$Tv\t$NucOrder\n";
	}else{#there is no coverage for that site in the bam
	print OUT "$bam\t$gene\t$site\t$refseq{$gene}{$site}\t0\t\t\t\t\t\t\t\t\t\t";
	print OUT "\t\t\t\t\n";
	}
  }
  print LOG "Gene $gene Average entropy = ".$sumentropy/$nbsites."\n";
}
# Read mismatch table:
# ReadPos\tCntNonRef\tCntRef\tTotalCnt\tFreq\tAvQual (Freq=NonRef/TotalCnt)
open (READ, ">$stub\_read.txt")||die "can't open $stub\_read.txt\n";
print READ "ReadPos\tCntRef\tCntNonRef\tTotalCnt\tFreq\tAvQual\tAvQualRef\tAvQualNonRef\n";
foreach my $readpos (sort {$a<=>$b} keys %readinfo){
  print READ "$readpos\t".$readinfo{$readpos}{"CntRef"}."\t".$readinfo{$readpos}{"CntNonRef"}."\t";
  my $total = $readinfo{$readpos}{"CntNonRef"}+$readinfo{$readpos}{"CntRef"};
  my $freq = $readinfo{$readpos}{"CntNonRef"}/$total;
  print READ "$total\t$freq\t";
  print READ $readinfo{$readpos}{"AvQual"}/$total."\t";
  print READ $readinfo{$readpos}{"AvQualRef"}/$readinfo{$readpos}{"CntRef"}."\t";
  print READ $readinfo{$readpos}{"AvQualNonRef"}/$readinfo{$readpos}{"CntNonRef"}."\n";
}
# loop for the aa mutations and position of mismatches in the codon
if ($orfs){
  open (AA, ">$stub\_AA.txt")||die "can't open $stub\_AA.txt\n";
  # the position of the mutation
  print AA "Sample\tChr\tProtein\tAAPosition\tRefAA\tRefSite\tRefCodon\tFstCodonPos\tSndCodonPos\tTrdCodonPos\tCntNonSyn\tCntSyn\tNbStop\tTopAA\tTopAAcnt\tSndAA\tSndAAcnt\tTrdAA\tTrdAAcnt\tAAcoverage\n";
  foreach my $bam (keys %aafreq){
    foreach my $target (keys %{$aafreq{$bam}}){
      foreach my $prot (keys %{$aafreq{$bam}{$target}}){
        foreach my $aasite (sort {$a<=>$b} keys %{$aafreq{$bam}{$target}{$prot}}){
          print AA "$bam\t$target\t$prot\t$aasite\t";
          print AA $aafreq{$bam}{$target}{$prot}{$aasite}{"RefAA"}."\t";
          print AA $aafreq{$bam}{$target}{$prot}{$aasite}{"RefSite"}."\t";
          print AA $aafreq{$bam}{$target}{$prot}{$aasite}{"RefCodon"}."\t";
          print AA $aafreq{$bam}{$target}{$prot}{$aasite}{"firstcodonpos"}."\t";
          print AA $aafreq{$bam}{$target}{$prot}{$aasite}{"secondcodonpos"}."\t";
          print AA $aafreq{$bam}{$target}{$prot}{$aasite}{"thirdcodonpos"}."\t";
          print AA $aafreq{$bam}{$target}{$prot}{$aasite}{"nonsyn"}."\t";
          print AA $aafreq{$bam}{$target}{$prot}{$aasite}{"syn"}."\t";
          print AA $aafreq{$bam}{$target}{$prot}{$aasite}{"stop"}."\t";
           my $topAAs=0;
          foreach my $aa (sort { $aaorder{$bam}{$target}{$prot}{$aasite}{$b} <=> $aaorder{$bam}{$target}{$prot}{$aasite}{$a} } keys %{$aaorder{$bam}{$target}{$prot}{$aasite}}) {
            if ($topAAs<3){# provide the top three AAs and their counts
              print AA "$aa\t$aaorder{$bam}{$target}{$prot}{$aasite}{$aa}\t";
              $topAAs++;
            }
          }
          while ($topAAs<3){# if there are less than three different AAs
            print AA "\t\t";
            $topAAs++;
          }
          print AA $aafreq{$bam}{$target}{$prot}{$aasite}{"AAcoverage"}."\n";
        }
      }
    }
  }
}

# count mismatches between two strings
sub mismatch_count($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }

sub cntTsTv{
  #count transitions A<=>G (purine to purine) and C<=>T *pyrimidine to pyrimidines
  #count transversions A<=>C, A<=>T, G<=>T, G<=>C 
  my ($refbase, $Acnt, $Tcnt, $Gcnt, $Ccnt)= @_;
  my $Ts=0;
  my $Tv=0;
  if ($refbase=~/A/){
    $Ts=$Gcnt;
    $Tv=$Ccnt+$Tcnt;
  }
  if ($refbase=~/T/){
    $Ts=$Ccnt;
    $Tv=$Acnt+$Gcnt;
  }
  if ($refbase=~/G/){
    $Ts=$Acnt;
    $Tv=$Tcnt+$Ccnt;
  }
  if ($refbase=~/C/){
    $Ts=$Tcnt;
    $Tv=$Acnt+$Gcnt;
  }
  #print "$refbase A:$Acnt T:$Tcnt G:$Gcnt C:$Ccnt Transition:$Ts Transversion:$Tv\n";
  return($Ts,$Tv);
}

sub log_base {
    my ($base, $value) = @_;
    return log($value)/log($base);
}

sub progress_bar {
    my ( $got, $total, $width, $char ) = @_;
    if ($got/$total>=$increment){
    $width ||= 25;
    $char  ||= '=';
    my $num_width = length $total;
    local $| = 1;
    $increment=$increment+0.05;
#     printf "|%-${width}s| Processed %${num_width}s reads of %s (%.2f%%)\r", 
#         $char x (($width-1)*$got/$total). '>', $got, $total, 100*$got/+$total;
    printf "|%-${width}s| Processed %.2f%% of reads out of %s\r", 
        $char x (($width+1)*$got/$total). '>', 100*$increment, $total;

    }    
}