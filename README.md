# The 1-Code Data Framework

1-Code is a data representation framework with a growing collection of associated software,
initially designed in the context of the
[Vertebrate Genomes Project](http://vertebrategenomesproject.org) (VGP) to operate on all the forms of 
data involved in a DNA sequencing and assembly project.   Data is represented in a very simple ASCII file format that is easy for both humans and
programs to read and interpret.  Moreover, there is a corresponding compressed and indexed binary
version for each ASCII file so that production tools built with the 1-code library are very efficient in time and the
size of the data files they manipulate.  All fields are strongly typed and a specific collection of data types constituting a *schema* is
specified in a schema file, that is itself a 1-code data file.  A generic converter allows one to move between the ASCII and
binary representations of data, and another core tool is provided to validate files against a 
given schema.

To make the library and command line tools just type ```make``` in this top
level directory.  The package has no dependencies on other software.  The .md files contain documentation, ONElib.c & .h contain the C-library for developers, and the rest is for building the two generic tools, ONEview and ONEstat.  The subdirectory `DEVELOPMENT` contains work areas for developers, currently Richard and Gene.

The documents describing the framework, generic tools, and development library are as follows:

- [Framework Description](https://github.com/thegenemyers/ONE-Code/blob/master/Format-description.md): a description of 1-code schema's and the data representation framework.

- [Generic Tools](https://github.com/thegenemyers/ONE-Code/blob/master/Generic-tools.md): the man pages for the generic 1-code command line tools.

- [Library Interface](https://github.com/thegenemyers/ONE-Code/blob/master/Library-interface.md):
a standard C library interface for 1-code developers.

In addition, specific technical documentation can be found within individual source files.

Authors:  Richard Durbin & Gene Myers

Date: October 2022
