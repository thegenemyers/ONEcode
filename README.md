# The ONEcode Data Framework

ONEcode is a powerful general data representation framework with a growing collection of associated software. It was originally designed for genomic data in the context of the
[Vertebrate Genomes Project](http://vertebrategenomesproject.org) (VGP).   

Data is represented in a very simple ASCII file format that is easy for both humans and
programs to read and interpret.  Moreover, there is a corresponding compressed and indexed binary
version for each ASCII file so that production tools built with the
ONEcode library `ONElib` are very efficient in time and the
size of the data files they manipulate.  All fields are strongly typed
and a specific collection of data types constituting a *schema* that
is normaly encoded at the top of the file, itself in ONEcode format.

A generic viewer `ONEview` allows one to view ascii and binary files,
and convert between them, also supporting viewing of a subset of
objects in the file.  Another core tool `ONEstat` can provide various
checks including validating files against a given schema in a separate file.

To make the library and command line tools just type `make` in this top level directory:

```
make
make install // copies executables to ~/bin
```

To see a simple annotated ONEcode file and some example usage:

```
cat TEST/small.seq
make test
```

The package has no dependencies on other software.
The `.md` files contain documentation, `ONElib.c` and
`ONElib.h` contain the C code library for developers, and
`ONEview.c` and `ONEstat.c` encode their respective programs.

The subdirectory `SEQUENCE_UTILITIES` contains a set of sequence utilities to interconvert between (compressed) fasta/fastq, ONEcode, and BAM/SAM, report statistics, and flexibly extract sequences. It contains its own README.md file.
 
The documents describing the framework, generic tools, and development library are as follows:

- [Framework Description](https://github.com/thegenemyers/ONE-Code/blob/master/Format-description.md): a description of ONEcode schemas and the data representation framework.

- [Generic Tools](https://github.com/thegenemyers/ONE-Code/blob/master/Generic-tools.md): the man pages for the generic ONEcode command line tools.

- [Library Interface](https://github.com/thegenemyers/ONE-Code/blob/master/Library-interface.md):
a standard C library interface for ONEcode developers.

In addition, specific technical documentation can be found within
individual source files, in particular ```ONElib.h```.

Authors:  Richard Durbin & Gene Myers

Date: October 2022  
Updated: July 2024
