# ONEcode sequence utilities

This directory contains core DNA sequence utilities that interconvert
between ONEcode and fasta and other DNA sequence formats, summarise the contents of a sequence file, and flexibly extract a subset of sequences from a sequence file.

- **seqconvert** to convert between DNA file types.  Supported file
  types are: fasta[.gz], fastq[.gz], 1seq (if compiled with -DONEIO, as by
  default in this project), SAM/BAM/CRAM (if compiled with -DBAMIO,
  see below), and a custom binary sequence file format (now deprecated
  with 1seq preferred). Gzipped files are automatically decompressed
  and file types are automatically recognized based on the first
  character of their contents. When specifying output files with -o
  if the file type is not explicitly given then it will be 
  auto-recognized from the filename ending, understanding 
  "fa", "fa.gz", "fasta", "fasta.gz", "fq", "fq.gz", "fastq",
  "fastq.gz", and any extension starting with '1' as a 1seq file. 
  This also will
  homopolymer compress (-H "hoco"), and if writing to a 1seq file
  while doing this will store the offsets in the original file of each
  new sequence base position, allowing to unhoco (-U) back to the
  original sequence.
- **seqstat** to give information about the data in DNA sequence
  files, supporting all file types handled by seqconvert.
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

## Synopsis

An example usage pattern is given below.  This uses three files which
all contain the same set of 10 short sequences: small.fa (fasta), small.seq (text ONEcode), small.1seq (binary ONEcode).

```
> seqstat -b small.fa
fasta file, 10 sequences >= 0, 580 total, 58.00 average, 42 min, 72 max
bases
  a 202 34.8 %
  c 88 15.2 %
  g 98 16.9 %
  n 3  0.5 %
  t 189 32.6 %
  
> seqconvert -1 small.fa > small.1seq
reading from file type fasta
written 10 sequences to file type onecode, total length 580, max length 72
  
> seqstat -l small.1seq
onecode file, 10 sequences >= 0, 580 total, 58.00 average, 42 min, 72 max
approximate N50 60
length distribution (quadratic bins)
  40    1
  51    3
  57    2
  64    1
  70    3
// quadratic binning scales naturally with the expected standard deviation under Poisson sampling

> seqconvert -o A.fa.gz small.1seq    // -fa -z are implied by the output filename
reading from file type onecode  with 10 sequences totLen 580
written 10 sequences to file type fasta, total length 580, max length 72

> gunzip A.fa.gz
> diff A.fa small.fa            // these two files are the same

> seqconvert -t -1 -H -o B.1seq small.1seq
reading from file type onecode  with 10 sequences totLen 580
written 10 sequences to file type onecode, total length 406, max length 50
user    0.000690        system  0.002908        max_RSS 1441792 memory  16786601

> seqconvert -fa B.1seq
reading from file type onecode  with 10 sequences totLen 406
written 10 sequences to file type fasta, total length 406, max length 50
>seq1
ctagtagcgatatagtatagtatcatgcgagtgtagat
>seq2
ctactcgagctctatcacagactcgcgtcagactctacgtctctgtat
>seq3
catatctgtcgtnatgtagagagtagtagacactcagacgatcagacgt
>seq4
tgagcgagagagatgatagactcgagagctga
>seq5
tatcagcgagtagcgacagcactatatcatatag
>seq6
agagtgatatcatactagacatcacgatagatagtatatctctcgagta
>seq7
gctctgtatatgtctgtnactgtgtgatatgctagc
>seq8
cgagatctatacagtatcatactatatatatatgat
>seq9
tagagtgtatagtatcacatcgtatatacacata
>seq10
acatacatatgatgtacactctatagctgatgacgactatatacgagagt

> seqconvert -U -o C.fa B.1seq
reading from file type onecode  with 10 sequences totLen 406
written 10 sequences to file type fasta, total length 580, max length 72

> diff C.fa small.fa

> seqextract -fa -f seq6 -c 3:10-20 small.1seq
>seq3:10-20
cgtnaaatgt
>seq6:0-66
agagtgaatatcattaaactagacattcacgatagaaaattagttaattatttctctttccgagta

```
With larger data sets the 1seq files are not only smaller than compressed fasta, they also read much faster because parsing is trivial (composition reads in each sequence, even though the basic counts are available in the 1seq header).

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
- **array.[ch], dict.[ch], hash.[ch]** respectively implement auto-extending arrays over any type, dictionaries (hashes of strings) and general hashes of basic types (up to 64-bit).

## ONEcode .1seq format

The schema for the ONEcode 1seq format is:

```
  O S 1 3 DNA               sequence: the DNA string
  D I 1 6 STRING            id: (optional) sequence identifier
  D Q 1 6 STRING            quality: Q values (ascii string = q+33)
  D N 3 3 INT 4 CHAR 3 INT  non-acgt base: pos (0-indexed), base, number
```

The DNA type is native encoded in ONEcode, stored in binary as a 2-bit
per base (4 bases per byte) array. Using the ONElib.c library you can
access either this compressed version or the original. Non-acgt bases 
are converted to 'a' but their original values are recorded in N
lines, so that full sequence information is preserved and can be
restored on conversion back to fasta or other formats (and also used
by software accessing sequences in 1seq files directly. 

Richard Durbin
July 2024
Last edit: 14 August 2024
