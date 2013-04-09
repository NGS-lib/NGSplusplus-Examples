The code samples in this repository require NGS++ to be installed (i.e.: make install) on your computer (see www.ngsplusplus.ca).

How to compile  <a id="Compilation"></a>
-------------
To compile programs, type:
make

Programs description<a id="Description"></a>
-------------
DensityStrand: 

1.  Parse a sorted SAM file to determine the density of sequences on the "+" strand compared to the density of sequences on the "-" strand.
2.  The program will analyze a chromosome at a time and will produce a result file for each chromosome.
3.  The first part of the analysis is to produce a density vector for every position in the current chromosome based on the aligned sequences.
4.  Once the density vector is done, a 200 nucleotides sliding window will move along every positions in the current chromosome and output regions that contains at least 30 sequences.
5.  The output will be a .txt file with 2 columns. Each row contains the score for one region. The first columns is for the sequence count on the "+" strand and the second column is for the sequence count on the "-" strand.

Note: The count is based on the average aligned sequence length for the sequences on the current chromosome. The sequence count for each sliding window region is based on the average lenght, which means that the result will not necessarily be an integer.
