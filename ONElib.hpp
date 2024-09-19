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
  C_1F::OneFile *vf ;
  
 public:
  ONEfile (const string &path) // just to open an existing file for reading with one thread
    { vf = NULL ;
      vf = C_1F::oneFileOpenRead (path.c_str(), 0, 0, 1) ;
      if (vf == NULL) { throw runtime_error("failed to open ONEfile") ; } 
    }
  ONEfile (const string &path, int nthreads) // open an existing file for reading with nthreads
    { vf = NULL ;
      vf = C_1F::oneFileOpenRead (path.c_str(), 0, 0, nthreads) ;
      if (vf == NULL) { throw runtime_error("failed to open ONEfile") ; } 
    }
  ONEfile (const string &path, const string &mode, const ONEschema &schema, const string &type, int nthreads) // full version
    { vf = NULL ;
      const char* tc = (type.size() > 0) ? type.c_str() : 0 ;
      if (mode[0] == 'r')
	vf = C_1F::oneFileOpenRead (path.c_str(), schema.os, tc, nthreads) ;
      else if (mode[0] == 'w' && mode[1] == 'b')
	vf = C_1F::oneFileOpenWriteNew (path.c_str(), schema.os, tc, true, nthreads) ;
      else if (mode[0] == 'w' && mode.size() == 1)
	vf = C_1F::oneFileOpenWriteNew (path.c_str(), schema.os, tc, false, nthreads) ;
      if (vf == NULL) { throw runtime_error("failed to open ONEfile") ; } 
    }
  ONEfile (const string &path, const string &mode, ONEfile &from, int nthreads)
    { vf = NULL ;
      if (mode[0] == 'w' && mode[1] == 'b')
	vf = C_1F::oneFileOpenWriteFrom (path.c_str(), from.vf, true, nthreads) ;
      else if (mode[0] == 'w' && mode.size() == 1)
	vf = C_1F::oneFileOpenWriteFrom (path.c_str(), from.vf, false, nthreads) ;
      if (vf == NULL) { throw runtime_error("failed to open ONEfile") ; } 
    }
  ~ONEfile () { C_1F::oneFileClose (vf) ; }

  bool      checkSchemaText (const string &text)
    { return C_1F::oneFileCheckSchemaText (vf, text.c_str()) ; }

  char      readLine() { return C_1F::oneReadLine (vf) ; }

  int64_t   length()
    { return ((vf)->field[((vf)->info[(int)(vf)->lineType]->listField)].len & 0xffffffffffffffll) ; }
  int64_t   getInt(int x) { return vf->field[x].i ; }
  void      setInt(int x, int64_t val) { vf->field[x].i = val ; }
  double    getReal(int x) { return vf->field[x].r ; }
  void      setReal(int x, double val) { vf->field[x].r = val ; }
  char      getChar(int x) { return vf->field[x].c ; }
  void      setChar(int x, char val) { vf->field[x].c = val ; }
  int64_t*  getIntList() { return (int64_t*) C_1F::_oneList(vf) ; }
  double*   getRealList() { return (double*) C_1F::_oneList(vf) ; }
  char*     getDNAchar() { return (char*) _oneList(vf) ; }
  uint8_t*  getDNA2bit () { return (uint8_t*) C_1F::_oneCompressedList(vf) ; }
  string    getString() { return (char *) C_1F::_oneList(vf) ; }
  char*     nextString(char* s) { return s + strlen(s) + 1 ; }

  char*     getComment() { return C_1F::oneReadComment(vf) ; }

  bool      gotoObject(char lineType, int64_t i) { return C_1F::oneGoto (vf, lineType, i) ; }

  // some extra functions to hide readable class attributes

  char      lineType() { return vf->lineType ; }
  int64_t   lineNumber() { return vf->line ; }
  int64_t   givenCount(char lineType) { return vf->info[lineType]->given.count ; }
  int64_t   givenMax(char lineType) { return vf->info[lineType]->given.max ; }
  int64_t   givenTotal(char lineType) { return vf->info[lineType]->given.total ; }
  int64_t   currentCount(char lineType) { return vf->info[lineType]->accum.count ; }
  int64_t   currentMax(char lineType) { return vf->info[lineType]->accum.max ; }
  int64_t   currentTotal(char lineType) { return vf->info[lineType]->accum.total ; }
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
      cout << "sequence length " << of.length() << "\n" ;
}

#endif

#endif
