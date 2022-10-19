# The One-Code Data Framework

One-Code is a data representation framework with a growing collection of associated software,
initially designed in the context of the
[Vertebrate Genomes Project](http://vertebrategenomesproject.org) (VGP) to operate on all the forms of 
data involved in a DNA sequencing and assembly project.   Data is represented in a very simple ASCII file format that is easy for both humans and
programs to read and interpret.  Moreover, there is a corresponding compressed and indexed binary
version for each ASCII file so that production tools built with the One-Code library are very efficient in time and the
size of the data files they manipulate.  All fields are strongly typed and a specific collection of data types constituting a *schema* is
specified in a schema file, that is itself a One-Code data file.  A generic converter allows one to move between the ASCII and
binary representations of data, and another core tool is provided to validate files against a 
given schema.

To make the library and command line tools just type ```make``` in this top
level directory.  The package has no dependencies on other software.  The .md files contain documentation, ONElib.c & .h contain the C-library for developers, and the rest is for generic tools.  The subdirectory `DEVELOPMENT` contains work areas for developers, currently Richard and Gene.

The documents describing the framework, generic tools, and development library are as follows:

- [Framework Description](https://github.com/ONE-Code/blob/master/Core/Format-description.md): a description of One-Code schema's and the data representation framework.

- [Generic Tools](https://github.com/ONE-Code/blob/master/Core/Generic-tools.md): the man pages for the generic One-Code command line tools.

- [Library Interface](https://github.com/ONE-Code/blob/master/Core/Library-interface.md):
a standard C library interface for One-Code developers.

In addition, specific technical documentation can be found within individual source files.

Authors:  Richard Durbin & Gene Myers

Date: October 2023
