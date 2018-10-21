# CSHL Sequence Homology and Alignment Workshop
 
## CD Search and HMM Profiles Searching

- Search the *E.Coli* protein P45796 using [CD Search](https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi)
- Search the *E.Coli* protein P45796 using [HMMSCAN](https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan) using PFAM
- Search the *E.Coli* protein P45796 using [Psi-Search](https://www.ebi.ac.uk/Tools/sss/psisearch/)

```
>NP_416557.1 GDP-mannose 4,6-dehydratase [Escherichia coli str. K-12 substr. MG1655]
MSKVALITGVTGQDGSYLAEFLLEKGYEVHGIKRRASSFNTERVDHIYQDPHTCNPKFHLHYGDLSDTSN
LTRILREVQPDEVYNLGAMSHVAVSFESPEYTADVDAMGTLRLLEAIRFLGLEKKTRFYQASTSELYGLV
QEIPQKETTPFYPRSPYAVAKLYAYWITVNYRESYGMYACNGILFNHESPRRGETFVTRKITRAIANIAQ
GLESCLYLGNMDSLRDWGHAKDYVKMQWMMLQQEQPEDFVIATGVQYSVRQFVEMAAAQLGIKLRFEGTG
VEEKGIVVSVTGHDAPGVKPGDVIIAVDPRYFRPAEVETLLGDPTKAHEKLGWKPEITLREMVSEMVAND
LEAAKKHSLLKSHGYDVAIALES
```

---
Questions:
1. What domains are identified in using each method?
2. What is the expectation value for these hits?
3. What percent of the protein is "covered" by these domains?
4. How many iteration of psi-search is necessary to find the protein ABNA_ASPNG?
5. What domain(s) match ABNA_ASPNG and P45796? *Hint: See the Visual Output*

## Write a simple global nucleotide alignment program of your own for sequences of similar size (think amplicon sequences) using Python

The input for your program should be the scores for match and mismatches and 2 sequences and allow for an option to calculate the scores for each combination of sequences and their reverse complements ie sequence 1 vs sequence 2, reverse complement sequence 1 vs sequence 2, sequence 1 vs reverse complement of sequence 2 and reverse complement of both sequence 1 and 2.

1. What is the score of the alignment using the scores 1,-1 and the sequences:
- sequence1 = "agtctgtca"
- sequence1 = "gatctctgc"

2. What are the scores of each combinations of reverse complements using the score 2/-2.

## Write a fasta parse to determine percent alignment coverage for each hit

The objective of this exercise is to practice regular expressions using the [original FASTA program output](out.txt).

### Hints:
- Query sequences are denoted with a number and a record separator '>>>' then the name of the query and the length denoted with a number and aa to indicate amino acid
  - Example: 1>>>Query 1 Description Species or Organism information 787 aa
- The alignments section of the file contains information about the alignment and the library hit sequence
  - Each aligned hit begins with >>LibraryName Description with spaces
  - Proceeding lines include information about the alignment including:
    - score =  s-w opt: score
    - bit score = bits: bit score 
    - identify = number% identity
    - alignment lenght= number aa overlap

**Parse this file and create output with the following columns**

- Hit Name [example = SP:XYND_PAEPO]
- Percent Query Coverage [example 32%] ie alignment length/query length (ignore gaps)
- Bit Score [example=130.0]
- E-Value [example=1.6e-28]
- Alignment Length [example=131]
- Query Length example [635]
