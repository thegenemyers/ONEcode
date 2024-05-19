# ONEcode applications from Richard Durbin

This directory contains core DNA sequence utilities that interconvert
between ONEcode and fasta and other formats, summarise the
properties of a sequence file, and support extract sequences from a
sequence file.

- **seqconvert** to convert between DNA file types.  This also will
  homopolymer compress (-H "hoco"), and if writing to a 1seq file
  while doing this will store the offsets in the original file of each
  new sequence base position, allowing to unhoco (-U) back to the
  original sequence.  Finally, it can split a scaffold sequence into
  contigs at blocks of n's (-K), and if writing to a 1seq file then
  it stores the gap information in scaffold ('s') objects, and they
  can be rejoined back to the original sequences with -J. Because
  ONEcode files compress all lists, this is all quite space efficient.
- **seqstat** to give information about the data in DNA sequence
  files: fasta[.gz], fastq[.gz], 1seq (if compiled with -DONEIO, as by
  default in this project), SAM/BAM/CRAM (if compiled with -DBAMIO,
  see below), and a custom binary sequence file format (now deprecated
  with 1seq preferred). 
- **seqextract** to extract sequences or parts of sequences from any
  type of sequence file that seqconvert or the seqio library can read.

For each program, running it without any arguments gives usage
information.

Note that this directory contains **copies** of ../../ONElib.[ch]
rather than directly linking against them.  It is therefore
standalone.  This is my standard pattern for using ONElib: I copy a
working version of ONElib.[ch] into the new project to ensure there
are no dependencies and it stays stable.  I can update the copies to
acquire new features whenever I want.

## Building
```
  git clone https://github.com/thegenemyers/ONEcode
  cd APPLICATIONS/Durbin
  make
```

If you want to be able to read SAM/BAM/CRAM files then you need to install htslib in a parallel directory and use Makefile.bam:
```
  pushd ../../..
  git clone https://github.com/samtools/htslib.git
  cd htslib
  autoreconf -i  # Build the configure script and install files it uses
  ./configure    # Optional but recommended, for choosing extra functionality
  make
  make install
  popd
  make clean
  make -f Makefile.bam
```

## Libraries

There are various useful libraries, with header files:

- **ONElib.[hc]** symbolic links up to the versions in ../...
- **seqio.[hc]** supports reading, writing DNA files with a few other basic operations.  Implementation in seqio.c, with dependencies on utils.[hc], libz, ONElib and htslib depending on compile operations.
- **utils.[hc]** some very low level type definitions (e.g. I8 to I64 and U8 to U64), die(), warn(), new(), new0(), and a timing package.  No dependencies beyond normal C run time library.  NB there is a handy fzopen() which will silently open .gz files as standard files, but this depends on funopen() which is not available on all systems.  If this does not compile/link then you will need to undefine WITH_ZLIB to link.
- **array.[ch], dict.[ch], hash.[ch]** respectively implement advanced language style extendable arrays, dictionaries (hashes of strings) and general hashes of basic types (up to 64-bit).

## Synopsis

An example usage pattern is given below.  This uses three files which
all contain the same set of 10 short sequences: small.fa (fasta
), small.seq (text ONEcode), small.1seq (binary ONEcode).

```
> seqstat -b small.fa
onecode file, 10 sequences >= 0, 577 total, 57.70 average, 42 min, 72 max
bases
  a 202 35.0 %
  c 88 15.3 %
  g 98 17.0 %
  t 189 32.8 %
  
>seqconvert -1 small.fa > small.1seq
  
>composition -l small.1seq
fasta file, 10 sequences >= 0, 577 total, 57.70 average, 42 min, 72 max
approximate N50 60
length distribution (quadratic bins)
  40    1
  46    1
  51    2
  57    2
  64    1
  70    3
// quadratic binning scales naturally with the expected standard deviation
under Poisson sampling

> seqconvert -o A.fa.gz small.1seq    // -fa -z are implied by the output filename
reading from file type onecode  with 10 sequences totLen 577
written 10 sequences to file type fasta, total length 577, max length 72
user    0.000156        system  0.001141        max_RSS 1359872 memory 33555094

>gunzip A.fa.gz
>diff A.fa small.fa            // these two files are the same

> seqconvert -1 -H -o B.1seq small.1seq
reading from file type onecode  with 10 sequences totLen 577
written 10 sequences to file type onecode, total length 404, max length 50
user    0.000432        system  0.002380        max_RSS 2998272 memory
16864607

> seqconvert -fa B.1seq
reading from file type onecode  with 10 sequences totLen 404
written 10 sequences to file type fasta, total length 404, max length 50
user    0.000282        system  0.001571        max_RSS 1671168 memory  33555038
>seq1
ctagtagcgatatagtatagtatcatgcgagtgtagat
>seq2
ctactcgagctctatcacagactcgcgtcagactctacgtctctgtat
>seq3
catatctgtcgtatgtagagagtagtagacactcagacgatcagacgt
>seq4
tgagcgagagagatgatagactcgagagctga
>seq5
tatcagcgagtagcgacagcactatatcatatag
>seq6
agagtgatatcatactagacatcacgatagatagtatatctctcgagta
>seq7
gctctgtatatgtctgtactgtgtgatatgctagc
>seq8
cgagatctatacagtatcatactatatatatatgat
>seq9
tagagtgtatagtatcacatcgtatatacacata
>seq10
acatacatatgatgtacactctatagctgatgacgactatatacgagagt

> seqconvert -U -o C.fa B.1seq
> diff C.fa small.fa

> ./seqextract -fa -f seq6 -c 3:10-20 small.1seq
>seq3:10-20
cgtaaatgta
>seq6:0-66
agagtgaatatcattaaactagacattcacgatagaaaattagttaattatttctctttccgagta

```
With larger data sets the 1seq files are not only smaller than compressed fasta, they also read much faster because parsing is trivial (composition reads in each sequence, even though the basic counts are available in the 1seq header).
