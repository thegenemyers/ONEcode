# ONEcode Data Files and Schemas

### Authors:  Gene Myers & Richard Durbin
### Last Update: July 26, 2024

## 0. Introduction

This document describes ONEcode data files and the schemas that describe them.  Schemas themselves are defined by a series of ONEcode data lines following a predefined metaschema.

We describe here the ASCII version of a data file.  There is automatically also a binary version encoding the same information that is much more compact and efficient that can be produced either by converting from an ASCII version with ONEview, or produced directly with our C support library.  By convention The binary versions are distinguished by having the same suffix as the ASCII version, but with a "1" prepended, e.g. dataset.1seq.

Each ONEcode file is a container for *objects*, each of which has a *type*. Each object is made of *records*, which are each a line in There is a schema for each file that describes the types that are allowed, and what 

There are *primary* types and *secondary* types which are specializations of a primary type for
particular applications, requiring particular properties or line types to be present by convention (that is not checked or enforced by the ONEcode framework).
Every primary and secondary file type is assigned a short string, that must be used as a filename's suffix, e.g. "seq" as in dataset.seq or "kmer" as in bigdata.kmer.  For example, a primary data type could be sequence, with secondary types k-mer, illumina read-pair, pacbio long read, or contig that are specializations of sequence to particular types of sequences.

###Philsophy

The design of the ONEcode data file format is based on the following principles:

1.	The format should be trivial to parse as input.  All the burden of encoding is placed on the software that produces the formatted file. ONElib deals with most of this burden automatically.

2.	The length of a list should always precede the list (e.g. 3 xxx).  Bracket (.eg [xxx]) or terminator (e.g. xxx0) constructions are not permitted.

3.	The total number of items in the file, and the maximum size and total size of variable-length items,
are specified in the header of a file so that memory allocation can be performed
once at the start by the reader if desired.

4.	There is a generic tool `ONEview` that translates the ASCII format
to a terse, indexable binary form and can invert this, as well as a code-level interface `ONElib` to the binary encoding.

5.	The ASCII form should not be overly verbose. It is not primarily for human consumption, but must be easily interpretable by a human being when and if they need to look at the data.

6.	Complex identifiers or symbolic names for objects are not used.  The n'th occurrence of an
object of a given type is referred to by the number n in any future context within the file. Object numbering starts at 1.  Object number 0 refers to the first data line following the file header.

7.	A reference to an object is not permitted before its definition.

8.	The provenance of a file must be given as the sequence of tools applied to obtain it,
where each tool application is described by its name, version number, command line applied,
and date & time performed. 

Together, these design decisions simplify many things.
For example, an ASCII ONEcode data file can be read by any simple parser that can
(a) read the next character, (b) read the next integer, (c) read the next real number, (d) read the next n symbols, and
(e) skip to the next line, where tokens are separated by a single white space character.
Moreover, any program reading in ONEcode data should never require dynamic memory
allocation to store objects – the size of any object or object collection should be given
before the first object of that type is described.  We realized this by designing the ONEcode
data files so that each begins with a header section that gives the number of objects, maximum
size of an object, and total size of all objects, along with the provenance information of
Principle 8.

So every ASCII ONEcode data file begins with a header segment followed by a data segment in order to abide by Principles 3 and 8.  Headers, relevant for all
file types, are described in Section 1 following this introduction, and the encoding of the data is documented in Section 2.
Schemas, that are themselves in ONEcode format, describe the data content of a data file and are documented in Section 3.
The utility program **ONEstat** validates the format of any given file against its internal or user-supplied schema, and also
can optionally reconstruct a file's header count information given only the data part.

Each line in either the header or data section starts with a single character (the "one-code") which determines the encoding of the remainder of the line. Header lines are specified by *non-alphabetic* characters, and data lines by *alphabetic* letters.  Lines starting with a period '.' are comment lines. An important restriction is that any given data line can contain only one list item.  Consider as an example a hypothetical Illumina read pair file of type .irp which is a specialization of the sequence type .seq.     

```
     1 3 seq 2 1            header start: primary file type, major and minor version
     2 3 irp                optional secondary subtype
     ! 7 VGPpair 3 0.1 ...  provenance line saying how file arose
     . schema
     ~ O P 0                defines object type P for read pair
     ~ G S                  tells us that P lines group S lines
     ~ O S 1 3 DNA          defines object type S for sequence
     . count information
     # P 3                  number of read pairs in file
     # S 6                  number of sequences in file
     @ S 5                  maximum number of bp’s in a read
     + S 26                 total number of bp's in file
     . start of Data Section
     P                      start first read pair 
     S 5 acgta              forward sequence of pair 1
     S 3 gtt                reverse sequence of pair 1
     P                      start of pair 2
     S 4 gcta               sequence 3 paired with ...
     S 5 ggtac              sequence 4
     P                      pair 3
     S 4 atta               sequence 5 paired with ...
     S 5 cctac              sequence 6
```

The header codes in the example are '1', '2', '!', '#', '@', and '+', and the data codes are either 'P' or 'S'.  The sequence of items on the remainder of the line is determined by the code.  Tokens are separated by a single space character.  Variable length lists,
including strings, are preceded by their length in keeping with Principle 2.
Thus, the specific ONEcode and
any subsequent list length items determine when the encoding of information on a line is at an end, so that any additional text on a line is understood to be a comment, as in the example, or auxiliary information providing a mechanism for extensibility.  In addition, a line with a ONEcode of '.' is also considered a comment in its entirety.

The first header line must always be a 1-line confirming the primary file type and specifying
the major and minor version numbers separated by whitespace.
This can optionally be followed immediately by a file subtype line, a 2-line, that gives the
secondary file type.
One or more provenance lines (!-lines) then inform one about how the
particular file came to be.  This is followed by the schema definition (see below). Additional header lines give information about the number of items
in the file (#-lines), the maximum length of lists (@-lines), and the total number of items in
a given list class (+-lines).  The data segment in this simple example then consists of pairs
of reads encoded by P-lines indicating the beginning of a pair, followed by two S-lines
giving the two sequence reads in the pair.

Conceptually ONEcode files are *immutable*, meaning that we do not expect their contents to change.
This means that subsequently in the same file, or more often in future files in a pipeline, we
can refer to objects by their ordinal position in a file, 1...n, not requiring named identifiers
(Principle 6). 

The type and arguments of each data line in any ONEcode data file are specified by a schema.  These schemas are part of each data file so that the file is self documenting.
The one exemption is potentially an ASCII file produced by an external program, in which case its schema can subsequently be embedded with one of the generic tools.  The schema for the example above is:

```
P 3 seq           primary type is .seq
S 3 irp           secondary type is .irp
O P 0
O S 1 3 DNA       objects are DNA sequences given in S-lines
D P 0             sequences are paired by proceeding their S-lines with a P-line
``` 

A schema file consists of predefined 'P', 'S', 'O', and 'D' lines.  A single, leading P-line specifies the primary type file extension and may be followed optionally be an S-line giving the
secondary type file extension.  The O-line specifies that S-lines in the data file are the objects described in the file and they have a single DNA string as an argument.  The D-line specifies that P-lines are auxiliary lines giving additional information, in this example, the line has no arguments.


## 1. ONEcode Headers

Considerable effort is invested on headers in ONEcode in keeping with Principles 3 and 8.
We will specify ONEcode syntax with "casual" context-free grammar rules.
Tokens are between angle-braces ```<>``` and may optionally have a descriptive tag following
a colon after the token name, e.g. ```<string:file_name>```.  The token name can be a predefined
type, such as ```int```, ```real```, or ```char```, in which case the syntax of the token is
the usual lexical syntax for such literals.  The right hand side of a rule can also be a
regular expression, where ```*``` denotes 0-or-more, ```+``` 1-or-more, ```^n``` exactly n
repetitions, ```|``` separates two alternatives, an item between square brackets ```[]``` is
optional, and parentheses may be used to disambiguate precedence.  For example, a string
is defined by the rule:

```
    <string> = <int:n> <char>^n
```

The first line of any header must be the primary file type declaration that has
the synatx:

```
    <primary_type_header> = 1 <string:file_type> <int:major> <int:minor>
```
where the initial ```1``` indicates that this is a "ONEcode" file (as well as this being line 1
&#x1F609;)
and ```<file_type>``` gives the file suffix string.

The initial header line can be followed by an optional subtype line

```
    <subtype_header> = 2 <string:file_subtype>
```
where now file type is one of the suffixes above for a secondary file type.

How a file was produced, or its provenance, we believe is very important and users are encouraged to include a sequence of provenance or !-lines that each record a processing step that
was involved in producing the current file.  Each line contains four strings giving (a) the
program name, (b) the version of that program as a string, (c) the command line that was executed, and (d) the date and time it was run.

```
    <provenance_step> = ! <string:name> <string:version> <string:command> <string:date>
```

The schema for a data file is always embedded within the file save in those occasional cases where the file is created by an external agent.  If present it is listed in the header with each schema line (see Section2 below) preceeded with a ONEcode of '~': 

```
    <schema_line> = ~ ( <primary_type> | <secondary_type|  <data_line> )
```
See Section 2 for the definition of the schema line types.  For example, the schema for our example, when embedded in the header has the form:

```
     ~ P 3 seq           primary type is .seq
     ~ S 3 irp           secondary type is .irp
     ~ O S 1 3 DNA       objects are DNA sequences given in S-lines
     ~ D P 0             sequences are paired by proceeding their S-lines with a P-line
``` 

Next there are three header line types - #, +, and @ - that allow one to specify the number, total size, and maximum size of objects across a file.  These all have the syntax:

```
    <size_header> = [#+@] <symbol:C> <int>
```
\#-lines tell you the number of lines in the file of type C.  For line types that encode a
list of items (recall they can have only one list), such as ```<string>``` (a list of characters) or say a list of restriction map
sites, a +-line tells you the total number of items in all the lists in the file,
e.g. "<code>+ S 26</code>" in the example of the introduction indicates that altogether the sequences in the file total 26 bases.  Similarly, an @-line indicates the length of the largest list that occurs in any given line of the specified type (assuming the line type has a list).
These "limit" lines are only present for data line types that occur in the body and not lines
that occur in the header itself.  There is one unusual case and that is for a list of strings
because strings are themselves lists and so in some sense this data type is a list of lists.
In this case, ```+``` is the sum of all of the list element lengths over all lines of the given type, and ```@``` is the maximum of the sum of the list lengths of each line.

Often the objects in a data file are naturally partitioned into groups, e.g. all the read pairs in a flow-cell lane, or all the read pairs in a "cloud".   ONEcode supports this with the concept
of "group" data lines that by convention have lower case symbols as the line type indicator.
For a group, one would like to know the maximum number of items of a given type in any group
and the maximum size needed to contain all the list objects occurring in any group.  So in the
header, a %-line indicates the maximum number of objects or maximum total size of
list objects within any given group type.  The syntax for these lines is:  

```
    <group_header> = % <symbol:G> [#+] <symbol:C> <int>
```
where G is the group line designator and C is the line type in question.
(The ```@``` reduction is not necessary as it is the same as the global maximum.)

Another important header line type indicates that references in this file are made to objects
in another file. This has the syntax:

```
    <reference_header> = < <string:file_name> <int:nx>
```
All the objects (all of the same type) in the specified file are available and
and ```nx``` indicates the number of these items in the file (hence
the range of reference indices is ```[1,nx)```).  For example a hypothetical alignment file would refer to sequence objects in another file of type sequence.

A related concept is to specify another file whose objects depend on the objects
in the current file.
We denote these with a '>'-line that has the opposite direction to the '<' of the reference line
above.

```
    <forward_header> = > <string:file_name>
```
In this case there is no need to indicate the number of objects in the file, since the current
file will not refer to them.


In summary, every (complete) ONEcode data file begins with a header.  Every header starts with a primary type
line optionally followed by a subtype line and provenance information.  Then ensue a number of size lines for every relevant
data line of the file type.  And finally, at the end, any relevant reference- and forward-lines.  In a rule:

```
    <header> = <primary_type_header> [<subtype_header>] <provenance_step>+ <schema_line>*
                   (<size_header>|<group_header>)+ (<reference_header>|<forward_header>)+
```

## 2. ONEcode Data

Every data line begins with a type designated by a single alphabetic letter.
This is then followed by a number of arguments as defined by the schema for the file.
All tokens are separated by a single space character, and lists are proceeded by an
integer giving the length of the list.  After the last argument, any text to the end of the line is considered a comment.

```
    <data>      = <data_line>*
      <data_line> = <char:ONEcode> <datum>*
```
ONEcode supports arguments of the following types: CHAR, INT, REAL, STRING, DMA, INT\_LIST, REAL\_LIST, and STRING\_LIST.  One does not need to worry about the word-size of integers or reals, they both accommodate the machine maximum of 64-bits but internally are encoded more efficiently when possible.  Note that a STRING is just a CHAR\_LIST, and DNA represents a STRING over the letters acgt where case is ignored.  DNA strings are handled internally in a more efficient way then general strings, hence the distinction.  Finally, note that a STRING\_LIST is subtle in that it is technically a list of lists, e.g. ```3 1 a 5 small 7 example```, requiring clarification for the meaning of any relevant size lines in the header.

```
    <datum>       = <scalar_data> | <list_data>
      <scalar_data> = <char> | <int> | <real>
      <list_data>   = <string> | <list_of(int)> | <list_of(real)> | <list_of(string)> | <dna> 
        <list_of(T)>  = <int:n> <T>^n
          <dna>         = <list_of([acgtACGT]>
```

There is always one line type that specifies the object being encoded in a data file.  All other data lines provide auxiliary information, and are associated with an object by convention, e.g. those lines immediately following an object line, or say the P-lines applying to the next two objects in our running example.

The one special auxiliary line is a group line type, of which there can currently only be one in the current implementation.  By convention a group line begins with a lower-case alphabetic letter, whereas normal data lines begin with a capital letter.  The group lines conceptually partition the file into the segments between them.  If there are data lines before the first group line, then those lines are not in any group.  In the binary ONEcode version, the number of objects in each group is always the first field of the line, and the file can be indexed by group (as well as by object).  In the ASCI version, group size is not necessarily known, in which case the first field has value 0.

## 3. ONEcode Schemas

The ONEcode framework allows one to encode almost any kind of data.  A schema for a primary or secondary file type specifies the type of lines that can be in a data file of that type and their arguments.  The schema applicable to a given file is given in the header with the exception that one can produce an ASCII data file without header or schema and subsequently create the header and associate the schema with the generic ONEcode tools.

Schemas are themselves specified in the ONEcode format (there is a schema for schemas &#x1F609;) with the following ONEcode line types.  A schema always begins with a P-line that has a string argumnt specifying the primary file type suffix extension.

```
    <primary_type> = P <string:file_type>
```

This can then be optionally followed by an S-line that similarly specifies the secondary suffix extension.  In this case, the schema is defining the secondary file type.

```
    <secondary_type> = S <string:file_type>
```

The remaining lines specify data lines by giving their defining ONEcode character followed by the argument types for such a line:

```
    <data_line> = [ODG] <char> <field_list>
```

An O-line specifies that the data line being defined is considered the object of the data file type.  There can only be one such line.  A G-line specifies that the data line being defined is a group line.  While conceptually there is no reason for a limitation, the current implementation of ONEcode allows only one such line ins a schema.  Finally, D-lines define any number of auxiliary data lines that augment the content of the primary objects.

```
      <field_list>   = <int:n> <field>^n
        <field>        = <scalar_field> | <list_field>
          <scalar_field> = '4 CHAR' | '3 INT' | '4 REAL'
          <list_field>   = '6 STRING' | '8 INT_LIST' | '9 REAL_LIST' | '11 STRING_LIST' | '3 DNA' 
```

The fields of a line are specified with the special strings CHAR, INT, etc. corresponding to the data types supported by the ONEcode framework.

In summary a schema is:

```
    <schema> = <primary_type> [<secondary_type] <data_line>*