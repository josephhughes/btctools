BTCtools
========

Beyond the consensus NGS tools
==============================



Single SAM
----------
* Sequencing error estimates

dN/dS (need aa loop)

* entropy

* Per site entropy 

Per site dN/dS (aa loop)

* Minimum frequency threshold table

* GC content (add reference base in table)

Number of synonymous and non-syn changes of individual genes (need  aa loop)

Motif upstream of mutation (read loop, need to provide a fasta file for weblogo) this would mainly be useful for clonal NGS data, low frequency samples, for looking at sequencing error

Introduced stop codons (aa loop)

Read mismatch location (read loop : position of mismathces in each of the reads, get information about how many reads have 1, 2, 3 etcc mutations): Read position, number of mismatches, nb of reference bases, total number of bases, mismatch frequency

vPhaser - haplotype diversity 

Figures and tables from single SAM:
-------------------------
1) Sequence error estimates => Mutation spectrum = plot distribution of frequency estimates (figure)

2) dN/dS => single value per gene (table)

3) entropy => single value per gene (table)

4) Per site entropy => plot entropy along each gene and plot sliding window of entropy (figure)

5) Per aa site dN/dS => plot dN/dS along the length of each gene (figure)

6) Motif => web-logo plot  (figure)

7) Mismatch location plot from first base (figure)

8) mutations per codon position (table and figure)

9) gc versus coverage plot (figure)

10) mutation frequency per site coloured by coverage

11) per site coverage


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






