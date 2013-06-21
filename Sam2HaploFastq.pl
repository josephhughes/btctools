#!/usr/bin/perl

# use this to obtain the number of each base at each site from a bam file
# you need to give the bam file and the reference fasta

use strict;
use Getopt::Long; 
use Bio::DB::Sam;
use Bio::SeqIO;


my ($bam, $ref,$help, $out,%basefreq,%sitequal,$i,$chr,$start,$end,$threshold);

&GetOptions(
	    'bam:s'   => \$bam,#bam file (binary sam)
	    'ref:s'   =>  \$ref, #reference file in fasta  
        "out:s"   => \$out,
        "chr:s"   => \$chr,
        "start:i" => \$start,
        "end:i"   => \$end,
        "t=i"     => \$threshold,
        
           );

if (($help)&&!($help)||!($bam)||!($ref)){
 print "Usage : perl ParseSam.pl -bam S1_refHPAI_cons_stampy.bam -ref refHPAI_cons.fa -out S1_basefreq.txt \n";
 print " -bam <txt> - the input bam file\n";
 print " -ref <txt>  - the reference fasta file\n";
 print " -out - the output in text-tab delimited\n";
 print " -chr  - the name of the chromosome of interes\n";
 print " -start - start of region of interest\n";
 print " -end - end of region of interest\n";
  print " -t - minum length threshold\n";
 print " -help        - Get this help\n";
 exit();
 }

 # high level API
 my $sam = Bio::DB::Sam->new(-bam  => $bam,
                             -fasta=> $ref,
                             -expand_flags  => 1);
my $ins_cnt=0;
my $nocigs_cnt=0;

my $inref = Bio::SeqIO->new(-file => "$ref" , '-format' => 'fasta');

my $region;
while ( my $refseq = $inref->next_seq() ) {
    $region=lc($refseq->subseq($start,$end));
    
}

open(FASTA,">$out\.fa")||die "Can't open $out\.fa\n";
open(QUAL,">$out\.qual")||die "Can't open $out\.qual\n";

my @pairs = $sam->get_features_by_location( -type   => 'read_pair',
                                            -seq_id => $chr,
                                            -start  => $start,
                                            -end    => $end);

my $nbpairs=0;
my $max_length=0;
my $min_length=length($region);
my $above_thres=0;
my $fullcnt=0;
 for my $pair (@pairs) {
    my $length = $pair->length;   # insert length
    #print "PAIR $length\n";
    my ($first_mate,$second_mate) = $pair->get_SeqFeatures;
    if ($first_mate && $second_mate){
      my $f_start = $first_mate->start;
      my $f_end = $first_mate->end;
      my $f_strand = $first_mate->strand;
      my $f_id = $first_mate->seq_id;      
      my $f_dna = $first_mate->query->dna;
   
        
      #print ">$f_id start=$f_start end=$f_end strand=$f_strand\n$f_dna\n";
      my $s_start = $second_mate->start;
      my $s_end = $second_mate->end;
      my $s_strand = $second_mate->strand;
      my $s_id = $second_mate->seq_id;
      my $s_dna = $second_mate->query->dna;
      #print ">$s_id start=$s_start end=$s_end strand=$s_strand\n$s_dna\n";
      $nbpairs++; 
      my ($f_trimseq,$f_start,$f_end,@f_scores)=&trim($first_mate);
     # print "LEngth\n".scalar(@f_scores);
   #   my @str=split(//,$f_trimseq);
    #  print "\n".scalar(@str);
      #print "$f_start end $f_end\n";
      #print "\nCOUNT $nbpairs\n";
      $f_trimseq='-'x$f_start."$f_trimseq".'-'x$f_end."\n";
      my $f_score='- 'x$f_start.join(" ",@f_scores).' -'x$f_end."\n";
      my ($s_trimseq,$s_start,$s_end,@s_scores)=&trim($second_mate);
      #print "$s_start end $s_end\n";
      $s_trimseq='-'x$s_start."$s_trimseq".'-'x$s_end."\n";
      my $s_score='- 'x$s_start.join(" ",@s_scores).' -'x$s_end."\n";
      print "$region\n$f_trimseq\n$f_score\n$s_trimseq\n$s_score\n";
      
 	  #comparing each site 
 	  my @f_sites=split(//,$f_trimseq);
 	  my @s_sites=split(//,$s_trimseq);
 	  my @r_sites=split(//,$region);
 	  my @f_scores=split(/ /,$f_score);
 	  my @s_scores=split(/ /,$s_score);
 	  my %contig;
 	  my %qual;
 	  my $length=length($region);
 	  for (my $i=0; $i<scalar(@r_sites);$i++){
 	    #print "site $i $r_sites[$i] $f_sites[$i] $f_scores[$i] $s_sites[$i] $s_scores[$i]\n";
 	    if ($f_sites[$i]=~/\-/ && $s_sites[$i]=~/\-/){
 	      $contig{$i}=$r_sites[$i];
 	      $qual{$i}="0";
 	      $length--;
 	    }
 	    if ($f_sites[$i]!~/\-/ && $s_sites[$i]=~/\-/){
  	      $contig{$i}=$f_sites[$i];
 	      $qual{$i}=$f_scores[$i];
 	   }
 	   if ($f_sites[$i]=~/\-/ && $s_sites[$i]!~/\-/){
  	      $contig{$i}=$s_sites[$i];
 	      $qual{$i}=$s_scores[$i];
 	   }
 	   if ($f_sites[$i]!~/\-/ && $s_sites[$i]!~/\-/ && $f_sites[$i]==$s_sites[$i]){
 	      $contig{$i}=$f_sites[$i];
 	      
	      if ($f_scores[$i]>=$s_scores[$i]){
	        $qual{$i}=$f_scores[$i];# change this to do the average
	      }else{
	        $qual{$i}=$s_scores[$i];# change this to do the average
	      }
	   } 
 	   if ($f_sites[$i]!~/\-/ && $s_sites[$i]!~/\-/ && $f_sites[$i]!=$s_sites[$i]){
 	      if ($f_scores[$i]>=$s_scores[$i]){
 	      	$contig{$i}=$f_sites[$i];
	        $qual{$i}=$f_scores[$i];# change this to do the average
	      }else{
	        $qual{$i}=$s_scores[$i];# change this to do the average
	        $contig{$i}=$r_sites[$i];
	      }
	    }
	   
 	  }
      if ($length>$max_length){
        $max_length=$length;
      }
      if($length<$min_length){
        $min_length=$length;
      }
      if($length>=length($region)){
       $fullcnt++;
      }
      if ($length > $threshold){
 	  print FASTA ">seq$nbpairs LENGTH=$length\n"; 
 	  foreach my $pos (sort { $a <=> $b} keys %contig){
 	    print FASTA "$contig{$pos}";
 	  }
 	  print FASTA "\n";
 	  print QUAL ">seq$nbpairs LENGTH=$length\n"; 
 	  foreach my $pos (sort { $a <=> $b} keys %qual){
 	    print QUAL "$qual{$pos} ";
 	  }
 	  print QUAL "\n";
 	  $above_thres++;
      }
    }
 }

print "Total number coverage ".scalar(@pairs)."\n";
print "Number of pairs $nbpairs\n";
print "Maximum length $max_length\nMinimum length $min_length\n";
print "Above threshold $above_thres\n";
print "Full length $fullcnt\n";

sub trim {
  my ($mate)=@_;
  my $m_start = $mate->start;
  my $m_end = $mate->end;
  my $m_strand = $mate->strand;
  my $m_id = $mate->seq_id;      
  my $m_dna = $mate->query->dna;  
  #print "SUBROUTINE\n>$m_id start=$m_start end=$m_end strand=$m_strand\n$m_dna\n";
  my $m_trimdna=$m_dna;
  my @scores    = $mate->qscore; 
  my $start_gaps=0;
  my $end_gaps=0;
    if ($m_start<$start){
      my $start_trim=$start-$m_start;
      $m_trimdna=substr($m_trimdna,$start_trim);
      #print "Trim start $start_trim\n$m_trimdna\n";
      @scores = @scores[ $start_trim .. $#scores ];
    }
    if ($m_end>$end){
      my $end_trim=$m_end-$end;
      $m_trimdna=substr($m_trimdna,0,-$end_trim);
      #$m_trimdna=~ s/\w{$end_trim}$//;
      #print "Trim end $end $end_trim\n$m_trimdna\n";
      @scores = @scores[$#scores - $end_trim .. $#scores];
    }
  if ($start<$m_start){#do padding with gaps
      $start_gaps=$m_start-$start;
      #print "start $start gaps $start_gaps\n";
      #print "Gaps\n$m_trimdna\n";
  }
  if ($m_end<$end){
    $end_gaps=$end-$m_end;
    #print "end $end gaps $end_gaps\n";
  }
  #print "TRIM\n$region\n$m_trimdna\n$start_gaps\n$end_gaps\n";
  return($m_trimdna,$start_gaps,$end_gaps,@scores);

}

