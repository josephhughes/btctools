btctools
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

Analyses from single SAM:
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
Sequencing error estimates - consistancy between samples (replicate sample) => run independently the two  samples

MacDonald-Kreitman 

Randomisation differences

Multiple samples
----------------

Track position for multiple samples over time


viralbox, vbox, viraltools, …
These are a few suggestions for which I think we already have most of the code.



Comparison between two SAMs
---------------------------
1) Sequence error estimates - consistency between replicates 

2) differences in entropy tested by randomisation => plot of the entropy per site with significant difference in entropy.


Input
-----

One or two SAM file 
Start and stop position for each gene (gene labelled identical to reference in bam) Table: Gene(text)\tStart(digits)\tStop(digits)


TO DO:
Codon based mutations/aa mutations => non-synonymous/synonyous
Mutations at each codon position


Java frontend:
Removing previously created files
Checking progress of parsing
Talking to perl
Saving tables
Saving summaries
Saving Figures
Zooming
Cursor over peaks on graphs






