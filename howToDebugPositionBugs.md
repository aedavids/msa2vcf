# How to debug position bugs

If you have a position bug in your converted VCF file you can use [https://github.com/vgteam/vg](https://github.com/vgteam/vg) to get more information when you run 'vg construct'


```
$ utils/cleanRun.sh ../data/E_gen.txt ../data/E_gen.reference.fasta
```

The error will look like
```
+ ../vg/bin/vg construct -r ../data/E_gen.reference.fasta -v E_gen.vcf.gz
index file ../data/E_gen.reference.fasta.fai not found, generating...
error:[vg::Constructor] Variant/reference sequence mismatch: TC vs pos: 2299: TT; do your VCF and FASTA coordinates match?
Variant: E*01:01:01:01  2299    .   TC  T   0   PASS    .
zero ind: 2298 1-indexed: 2299
```
## approach 1 cut down the file
### step 1 find the sample in the vcf that has the bug. 
In our example it will be the entry for position 2299

if we look at the the position vector we should be able to figure out what the sample name is. In our example
the sample name is E*01:01:01:05

### step 2 create a msf file to test with
* create a new msf file from the one you used in utils/cleanRun.sh.
* all you need is the header
* next use utils/fetchMSASample.sh to grab the reference line and sample line
```
$ utils/fetchMSASample.sh ../data/E_gen.txt 'E*01:01:01:01' > 'E*01:01:01:01'.sequence
$ utils/fetchMSASample.sh ../data/E_gen.txt 'E*01:01:01:05' > 'E*01:01:01:05'.sequence
```

* cut and past the ref sequence and sample sequence into your test msf file. 
* the MSF file only has a single sample. making it easier to track down the bug

### step 3 make sure the test msf file reproduced the bug
run utils/cleanRun.sh

### Step 4 tracking down the bug
use ShortenMultipleAlignmentSequence.py to trim bases from reference sequence while preserving 
the alignment

```
(HLA_haplotype) $ python ShortenMultipleAlignmentSequence.py --help
usage: n inputFile outputFile [options]

positional arguments:
  n           The number of nucleotides to trim from the front of the
              reference sequenc
  inputFile   This file has 3 lines. The first line is the reference sequence.
              The second line is the reference sequence in multiple sequence
              alignment format. The third line is the sample sequence in
              multiple sequence alignment format.
  outputFile  The trimed version of the input file

optional arguments:
  -h, --help  show this help message and exit
(HLA_haplotype) $ 
```
  
example

```
(HLA_haplotype) $ cat testShorten
ACGTACGT
A|CGTA.CGT
-|-*--G---
(HLA_haplotype) $ python ShortenMultipleAlignmentSequence.py 2 testShorten testShorten.out
(HLA_haplotype) $ cat testShorten.out
GTACGT
GTA.CGT
*--G---

```

### clean up artifacts
```
rm *.vcf *.vg *.log ../data/*.fasta.fia *.dot *.pdf *.tbi
```

## approach 2, try and figure out how to recreate small unit test
issue: lines can be thousands of chars long. Hard to see the difference. Our approach create fasta like file
that we can compare easily

1) use utils/fetchMSASample.sh to get the reference in multiple sequence alignment format and the sample as one
long continuous string. use sed and tr to get the reference sequence as a single continous string

The reference in MSA format
```
$ utils/fetchMSASample.sh ../data/E_gen.txt > 'E*01:01:01:01'.sequence
```

The sample we are trying to debug
```
$ utils/fetchMSASample.sh ../data/E_gen.txt 'E*01:03:02:07.sequence' > 'E*01:03:02:07.sequence' 
```

remove the first line of the reference fasta file, and strip all new lines
```
 sed -e '1d' ../data/E_gen.reference.fasta | tr -d '\n' > E_gen.refernce.sequence
```

2) use fold to create files with  reasonable line lengths. 

```
$ fold -w100 'E*01:01:01:01'.sequence >  'E*01:01:01:01'.sequence.w100
$ fold -w100 'E*01:03:02:07.sequence' > 'E*01:03:02:07.sequence'.w100
$ fold -w100  E_gen.refernce.sequence >  E_gen.refernce.sequence.w100
```

3) use paste to interleave the files
```
(HLA_haplotype) $ cat f1
line1.1
line1.2
line1.3

(HLA_haplotype) $ cat f2
line2.1
line2.2
line2.3


# note f3 has a single \n in it, this will cause white space to be added between blocks
(HLA_haplotype) $ cat f3



(HLA_haplotype) $ !past
paste -d '\n' f1 f2
line1.1
line2.1

line1.2
line2.2

line1.3
line2.3

```

working example
```
AEDWIP we want to to past the ref, msfRef, sample, empty file

(HLA_haplotype) $ paste -d '\n' 'E*01:01:01:01'.sequence.w100 'E*01:03:02:07.sequence'.w100 >'E*01:01:01:01'.sequence.w100-'E*01:03:02:07.sequence'.w100
(HLA_haplotype) $ head 'E*01:01:01:01'.sequence.w100-'E*01:03:02:07.sequence'.w100
CAAAGTGCTGAGATTACAGGCGTGAGCCACCGCGCCCAGCCAGGACTAATTTCTAAGAGTGTGCAGAGATACCGAAACCTAAAAGTTTAAGAACTGCTGA
************************************************************************************----------------
TTGCTGGGAAACTCTGCAGTTTCCCGTTCCTCTCGTAACCTGGTCATGTGTCCTTCTTCCTGGATACTCATGACGCAGACTCAGTTCTCATTCCCAATGG
----------------------------------------------------------------------------------------------------
GTGTCGGGTTTCTAGAGAAGCCAATCAGCGTCGCCACGACTCCCGACTATAAAGTCCCCATCCGGACTCAAGAAGTTCTCAGGACTCA.GAGGCTGGGAT
----------------------------------------------------------------------------------------.-----------
C|ATGGTAGATGGAACCCTCCTTTTACTCCTCTCGGAGGCCCTGGCCCTTACCCAGACCTGGGCGG|GTGAGTGCGGGGTCGGGATGGAAACGGCCTCTA
-|----------------------------------------------------------------|---------------------------------
CCGGGAGTAGAGAGGGGCCGGCCCGGCGGGGGCGAAGGACTCGGGGAGCCGCGCCGGGAGGAGGGTCGGGCCGATCTCAGCCCCTCCTCGCCCCCAG|GC
-------------------------------------------------------------------------------------------------|--
```

4. if you did not use the paste trick to get a new line between blocks use an emacs macro to add white space 
