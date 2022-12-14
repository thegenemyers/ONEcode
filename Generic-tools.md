# 1-Code Generic Command Line Tools

### Authors:  Richard Durbin & Gene Myers
### Last Update: October 19, 2022

This document describes generic command line tools for interacting with One-Code files.

#### <code>1. ONEstat [-Hu] [-o \<name>] [-t \<type suffix>] \<input:ONE-file></code>**

ONEstat provides information about a 1-code data file.  Without arguments it validates an ASCII file, including reporting any missing header information, and states how many objects, groups, and lines it contains.  Details of how many lines of each type are present are available in the count '@' header lines output by the -H option.

The -H option calculates a full and correct header from the contents of the file, and writes it out in ASCII.  A header can be added to an ASCII file that is lacking one using ONEview as explained below.

The -u option outputs the number of bytes used by each line type.

The -o option redirects the output to the named file. The default is stdout.

The -t option specifies the file type, and is required if the inspected file is an ascii file without a header, but is not needed for a binary file or an ASCII file with a proper header.

#### <code>2. ONEview [-bhH] [-o \<filename>] [-t \<type suffix>] [-i \<ranges>] [-g \<ranges>] \<input:ONE-file></code>
	
ONEview is the standard utility to extract data from 1-code files and convert between ASCII and binary forms of the format.

The -b option outputs the file in binary.  The default is ASCII.  Note that the binary form is compressed and indexed, and should be the standard form for programmatic access.

The -h option drops the header lines from ASCII output.  It has no effect when writing binary because all binary files automatically generate a full header (much of which is actually written as a footer, since it can only be created once all the data is processed).

The -H option just prints out the header, in ASCII.

The -o option redirects the output to the named file. The default is stdout.

The -t option specifies the file type, and is required if the inspected file is an ASCII file without a header, but is not needed for a binary file or an ASCII file with a proper header.

The -i and -g options make use of the binary file indices to allow random access to arbitrary sets of ojects or groups.  Legal range arguments include "0-10" which outputs the first 10 items, "7" which outputs the eighth item (remember numbering starts at 0), or compound ranges such as "3,5,9,20-33,4" which returns the requested items in the specified order.

It is possible to stream from a binary file to ascii and back from ascii to binary, so a standard pattern is 
```
   ONEview -h <binary-file> | <script operating on ascii> | ONEview -b -t <type> - > <new-binary-file>
```
Note that the script can ignore the header information, which is reconstructed by the second call to ONEview.

Regrettably it is not possible to stream a binary file into ONEview, since reading a binary file requires a seek operation to the end of the file to read footer information before reading the rest of the file.  An intermediate binary file must therefore be created to pass from ASCII to binary and back again, as follows:
```
   ONEview -b -t <type> <ascii-file> > <binary-file>
   ONEview <binary-file> > <new-ascii-file>
```
This pattern has the effect of standardising an ASCII file, and is the recommended way to add a header to an ASCII 1-code file that lacks a header.  Although some format consistency checks will be performed, if you want to fully validate a 1-code file then use ONEstat.
