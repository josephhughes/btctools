btctools
========

Beyond the consensus NGS tools
==============================

Single SAM
----------
Sequencing error estimates

dN/dS

entropy

Per site entropy 

Per site dN/dS

Minimum frequency threshold table

GC content

Number of synonymous and non-syn changes of individual genes

Motif upstream of mutation above a threshold

Introduced stop codons

Read mismatch location


Comparison between two SAM
--------------------------
Sequencing error estimates - consistancy between samples

MacDonald-Kreitman

Randomisation differences

Multiple samples
----------------

Track position for multiple samples over time


viralbox, vbox, viraltools, …
These are a few suggestions for which I think we already have most of the code.

Analyses from single SAM:
-------------------------
1) Sequence error estimates => plot distribution of frequency estimates

2) dN/dS => single value per gene

3) entropy => single value per gene

4) Per site entropy => plot entropy along each gene and plot sliding window of entropy

5) Per site dN/dS => plot dN/dS along the length of each gene

6) Motif => web-logo plot

7) Mismatch location plot from first base


Comparison between two SAMs
---------------------------
1) Sequence error estimates - consistency between replicates 

2) differences in entropy tested by randomisation => plot of the entropy per site with significant difference in entropy.


Input
-----

One or two SAM file 






