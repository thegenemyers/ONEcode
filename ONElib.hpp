/******************************************************************************************
 *
 *  File: ONElib.hpp
 *    C++ header code for ONE file reading and writing
 *
 *  Authors: Richard Durbin (rd109@cam.ac.uk), Gene Myers (gene.myers@gmail.com)
 *  Copyright (C) Richard Durbin, Gene Myers, 2019-
 *  with considerable help from Regev Schweiger
 *
 *****************************************************************************************/
/*  Last edited: Jun  9 11:07 2023 (rd109) */

#ifndef ONEfile_h
#define ONEfile_h

#include <string>
#include <stdexcept>
#include <vector>
using namespace std ;

namespace C_1F {	// private namespace for C interface
extern "C" { 
#include "ONElib.h"
}
}

class ONEschema
{
 private:
  C_1F::OneSchema *os ;
  
 public:
  ONEschema (const string &text) { os = C_1F::oneSchemaCreateFromText (text.c_str()) ; }
  ~ONEschema () { C_1F::oneSchemaDestroy (os) ; }

  friend class ONEfile ; 
} ;

class ONEfile
{
 private:
  C_1F::OneFile *of ;
  
 public:
  ONEfile (const string &path) // just to open an existing file for reading with one thread
    { of = NULL ;
      of = C_1F::oneFileOpenRead (path.c_str(), 0, 0, 1) ;
      if (of == NULL) { throw runtime_error("failed to open ONEfile") ; } 
    }
  ONEfile (const string &path, int nthreads) // open an existing file for reading with nthreads
    { of = NULL ;
      of = C_1F::oneFileOpenRead (path.c_str(), 0, 0, nthreads) ;
      if (of == NULL) { throw runtime_error("failed to open ONEfile") ; } 
    }
  ONEfile (const string &path, const string &mode, const ONEschema &schema, const string &type, int nthreads) // full version
    { of = NULL ;
      const char* tc = (type.size() > 0) ? type.c_str() : 0 ;
      if (mode[0] == 'r')
	of = C_1F::oneFileOpenRead (path.c_str(), schema.os, tc, nthreads) ;
      else if (mode[0] == 'w' && mode[1] == 'b')
	of = C_1F::oneFileOpenWriteNew (path.c_str(), schema.os, tc, true, nthreads) ;
      else if (mode[0] == 'w' && mode.size() == 1)
	of = C_1F::oneFileOpenWriteNew (path.c_str(), schema.os, tc, false, nthreads) ;
      if (of == NULL) { throw runtime_error("failed to open ONEfile") ; } 
    }
  ONEfile (const string &path, const string &mode, ONEfile &from, int nthreads)
    { of = NULL ;
      if (mode[0] == 'w' && mode[1] == 'b')
	of = C_1F::oneFileOpenWriteFrom (path.c_str(), from.of, true, nthreads) ;
      else if (mode[0] == 'w' && mode.size() == 1)
	of = C_1F::oneFileOpenWriteFrom (path.c_str(), from.of, false, nthreads) ;
      if (of == NULL) { throw runtime_error("failed to open ONEfile") ; } 
    }
  ~ONEfile () { C_1F::oneFileClose (of) ; }

  bool      checkSchema (ONEschema &schema, bool isRequired)
  { return C_1F::oneFileCheckSchema (of, schema.os, isRequired) ; }
  bool      checkSchemaText (const string &text)
    { return C_1F::oneFileCheckSchemaText (of, text.c_str()) ; }

  char      readLine() { return C_1F::oneReadLine (of) ; }

  int64_t   listLength()
    { return ((of)->field[((of)->info[(int)(of)->lineType]->listField)].len & 0xffffffffffffffll) ; }
  int64_t   getInt(int x) { return of->field[x].i ; }
  double    getReal(int x) { return of->field[x].r ; }
  char      getChar(int x) { return of->field[x].c ; }
  int64_t*  getIntList() { return (int64_t*) C_1F::_oneList(of) ; }
  double*   getRealList() { return (double*) C_1F::_oneList(of) ; }
  char*     getDNAchar() { return (char*) _oneList(of) ; }
  uint8_t*  getDNA2bit () { return (uint8_t*) C_1F::_oneCompressedList(of) ; }
  string    getString() { return (char *) C_1F::_oneList(of) ; }
  vector<string> getStringList ()
  { char *s = (char*) C_1F::_oneList(of) ;
      int size = listLength() ;
      vector<string> vs(size) ;
      int i ; for (i = 0 ; i < size ; ++i) { vs[i] = strdup(s) ; s += strlen(s) + 1 ; }
      return vs ;
    }

  char*     getComment() { return C_1F::oneReadComment(of) ; }

  bool      gotoObject(char lineType, int64_t i) { return C_1F::oneGoto (of, lineType, i) ; }

  // routines to write

  void      setInt(int x, int64_t val) { of->field[x].i = val ; }
  void      setReal(int x, double val) { of->field[x].r = val ; }
  void      setChar(int x, char val) { of->field[x].c = val ; }
  void      writeLine(char lineType) { C_1F::oneWriteLine(of, lineType, 0, 0) ; }
  void      writeLine(char lineType, int64_t listLen, void *listBuf)
    { C_1F::oneWriteLine(of, lineType, listLen, listBuf) ; }
  void      writeLine(char lineType, string s)
  { C_1F::oneWriteLine(of, lineType, s.length(), (void*) s.c_str()) ; }
  void      writeLine(char lineType, const vector<string>& slist)
    { int i, size = slist.size() ;
      int totLen = 0 ; for (i = 0 ; i < size ; ++i) totLen += slist[i].length() ;
      char *s, *s0 = new char[totLen+size] ;
      for (i = 0, s = s0 ; i < size ; ++i, s += slist[i].length() + 1)
	strcpy (s, slist[i].c_str()) ;
      C_1F::oneWriteLine(of, lineType, size, s0) ;
    }
  void      writeLineDNA2bit (char lineType, int64_t dnaLen, uint8_t *dnaBuf)
  { C_1F::oneWriteLineDNA2bit (of, lineType, dnaLen, dnaBuf) ; }
  
  void      writeComment (string s) { C_1F::oneWriteComment(of, (char*)"%s", s.c_str()) ; }

  // provenance etc.

  bool addProvenance(string prog, string version, string commandLine)
  { return C_1F::oneAddProvenance(of, prog.c_str(), version.c_str(),
				  (char*)"%s", commandLine.c_str()); }
  bool inheritProvenance (ONEfile source) { return C_1F::oneInheritProvenance (of, source.of) ; }
  bool addReference(string filename, int64_t count)
  { return C_1F::oneAddReference(of, filename.c_str(), count) ; }
  bool inheritReference(ONEfile source) { return C_1F::oneInheritReference (of, source.of) ; }

  // some extra functions to hide readable class attributes

  char      lineType() { return of->lineType ; }
  int64_t   lineNumber() { return of->line ; }
  string    fileName() { return of->fileName ; }
  int64_t   givenCount(char lineType) { return of->info[(int)lineType]->given.count ; }
  int64_t   givenMax(char lineType) { return of->info[(int)lineType]->given.max ; }
  int64_t   givenTotal(char lineType) { return of->info[(int)lineType]->given.total ; }
  int64_t   currentCount(char lineType) { return of->info[(int)lineType]->accum.count ; }
  int64_t   currentMax(char lineType) { return of->info[(int)lineType]->accum.max ; }
  int64_t   currentTotal(char lineType) { return of->info[(int)lineType]->accum.total ; }
} ;

#ifdef TEST_HEADER

// to use this link this file to a filename ending .cpp and compile with -D TEST_HEADER e.g.
// ln -s ONElib.hpp ONEtest.cpp
// gcc -c ONElib.c
// g++ -D TEST_HEADER -o ONEtest ONEtest.cpp ONElib.o

#include <iostream>

static string schemaText = 
  "P 3 seq                 SEQUENCE\n"
  "S 6 segseq              segment sequences - objects are 1:1 with those in seg file\n"
  "S 7 readseq             read sequences\n"
  "O S 1 3 DNA             sequence: the DNA string\n"
  "D I 1 6 STRING          id - sequence identifier; unnecessary for segments\n" ;

int main (int argc, char *argv[])
{
  if (argc < 2) { cerr << "need a filename as argument\n" ; exit (1) ; }
  ONEschema os(schemaText) ;
  ONEfile of(argv[1], "r", os, "", 1) ;

  cout << "opened 1seq file " << string(argv[1]) << " with " << of.givenCount('S') << " sequences\n" ;
  
  while (of.readLine())
    if (of.lineType() == 'S')
      cout << "sequence length " << of.listLength() << "\n" ;
}

#endif

#endif
