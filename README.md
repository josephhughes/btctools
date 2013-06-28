BTCtools - Beyond the consensus NGS tools
========================================

Single SAM
----------

Input: 

* A single bam file

* A reference in FASTA format

* A text-tab delimited file with coding regions [optional]

Output:

STDOUT (or logfile):

* Entropy per gene

* Number of synonymous and non-syn changes of individual genes

* Number of reads with 1, 2 ,3 mismatches

Four files:
  1. A nucleotide base frequency table with:
    * Sequencing error estimates
    * Per site entropy 
    
  2. An amino-acid table with:
    *  dN/dS
    * stop codons
    
  3. A FASTA file with the sequence around the mismatches:
    * Motif upstream of mutation this would mainly be useful for clonal NGS data, low frequency samples, for looking at sequencing error

  4. Position of mismatches in a read
    * Read mismatch location (position of mismathces in each of the reads): Read position, number of mismatches, nb of reference bases, total number of bases, mismatch frequency

Use table 1 to calculate and plot:

  1. Sequence error estimates => Mutation spectrum = plot distribution of frequency estimates (figure)

  2. GC content versus coverage plot(add reference base in table)

  3. Per site entropy => plot entropy along each gene and plot sliding window of entropy (figure)
  
  4. Per site coverage

Use table 2:

  1. Per aa site dN/dS => plot dN/dS along the length of each gene (figure)
  
Use file 3:

  1. Motif => web-logo plot  (figure)
  
use table 4:

  1. Mismatch location plot from first base (figure)

  2. mutations per codon position (table and figure)


Comparison between two SAM
--------------------------
Sequencing error estimates - consistancy between samples (replicate sample) => run independently the two  samples and output a merged table (merge and/or filter based on threshold)
(merge script)

Same ordering of nucleotides in replicates => sites with consistent and sig mutations between the replicates  

MacDonald-Kreitman 

Randomisation differences - add columns with the entropy randomisation
(randomisation scripts)

Figures from two-samples
------------------------

Scatterplot of mismatches frequencies of rep 1 against rep 2 (coloured by coverage) 

Coverage plots of both

Mutation spectrum from both

Per site entropy for both mirrored on x-axis (after randomisation)



Multiple samples
----------------

Track position for multiple samples over time

Input
-----

One or two SAM file 
Start and stop position for each gene (gene labelled identical to reference in bam) Table: Gene(text)\tStart(digits)\tStop(digits)
Reference fasta file (one or two)

Output
------
Nucleotide Table:
Sample\tChr\tPosition\tRefBase\tCoverage\tAvQual\tAcnt\tApval\tCcnt\tCpval\tTcnt\tTpval\tGcnt\tGpval\tentropy(base e)\tentropy(base 2)\tsimpson\tNonRefCnt
CntTv\tCntTs\tCodonPos\tCntNonSyn\tCntSyn\tOrderATCG

Aa Table:
Sample\tChr\tAAPosition\tRefAA\tRefCodon\tCntNonSyn\tCntSyn\tTopAA\tTopAAcnt\tSndAA\tSndAAcnt\tTrdAA\tTrdAAcnt\tAAcoverage\tNbStop

=>Filtered table would need to parse through the bam again to exclude the sites that have not met the threshold (Version 2)

Fasta:
Fasta sequence of 10 up and down of a mistmatch (padding with the consensus)

Read mismatch table:
ReadPos\tCntNonRef\tCntRef\tTotalCnt\tFreq\tAvQual (Freq=NonRef/TotalCnt)

Merged Table output:
Merged Nucleotide Table, recalculate entropy, randomisation, 
Sample\tChr\tPosition\tRefBase\tCombinedCoverage\tCombinedAvQual\tCombAcnt\tCombApval\tCombCcnt\tCombCpval\tCombTcnt\tCombTpval\tCombGcnt\tCombGpval\tCombentropy(base e)\tCombentropy(base 2)\tCombsimpson\tCombNonRefCnt
NonRefCntSample1/CoverageSample1\tNonRefCntSample2/CoverageSample2\tRejectBasedOnOrder\tRejectBasedOnCombinedBinomial

Merged AA Table



TO DO:
======
Codon based mutations/aa mutations => non-synonymous/synonyous
Mutations at each codon position


Java frontend:
==============
Removing previously created files
Checking progress of parsing
Talking to perl
Saving tables
Saving summaries
Saving Figures
Zooming
Cursor over peaks on graphs






