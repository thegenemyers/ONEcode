/*  File: seqextract.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2023
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Sep 28 00:51 2024 (rd109)
 * Created: Thu Dec  7 22:26:05 2023 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "array.h"
#include "dict.h"
#include "seqio.h"
#include "ONElib.h"

#define VERSION "1.0"

typedef struct {
  U64 k ;
  U64 start, end ;
  bool isRC ;
} Frag ;

int fragOrder (const void *a, const void *b)
{
  Frag *fa = (Frag*)a, *fb = (Frag*)b ;
  if (fa->k < fb->k) return -1 ; else if (fa->k > fb->k) return 1 ;
  if (fa->start < fb->start) return -1 ; else if (fa->start > fb->start) return 1 ;
  if (fa->end < fb->end) return -1 ; else if (fa->end > fb->end) return 1 ;
  return 0 ;
}

DICT *dict ;

void usage (void)
{
  fprintf (stderr, "usage: seqextract <commands> <seqfile>\n") ;
  fprintf (stderr, "  extract sequence fragments from a sequence file:\n") ;
  fprintf (stderr, "    -fa                       output FASTA file\n") ;
  fprintf (stderr, "    -1                        output ONEcode file\n") ;
  fprintf (stderr, "    -I                        output identifiers - only aplies to ONEcode\n") ;
  fprintf (stderr, "    -o outfilename            output file [-] : autorecognizes .1*, .fa, .fa.gz\n") ;
  fprintf (stderr, "    -f id[:start-[end]]    fragment to extract - can do many of these\n") ;
  fprintf (stderr, "    -F fragfile               file of fragments to extract\n") ;
  fprintf (stderr, "    -c count[{:_}start-[end]] fragment to extract by position in file\n") ;
  fprintf (stderr, "    -C countfragfile          file of count fragments to extract\n") ;
  fprintf (stderr, "  add R at end to reverse complement, :R for whole sequence\n") ;
  fprintf (stderr, "    e.g. \"-f id\", \"-f id:R\", \"-f id:10-20\", \"-f id:10-20R\"\n") ;
  fprintf (stderr, "  you can escape colons in identifiers with \\, as in \"-f run5\\:read2\"\n") ;
  fprintf (stderr, "  start and end use 0-based coords with open end (so length = end-start)\n") ;
  fprintf (stderr, "    if no end then go to end of seq, so :0- for whole sequence\n") ;
  fprintf (stderr, "  count is 1-based (since 0 is before the first object)\n") ;
  
  exit (1) ;
}

void parse (char *s, U64 *len, U64 *start, U64 *end, bool *isRC)
{
  char *text = s ; *start = 0 ; *end = 0 ; *isRC = false ;
  while (*s && *s != ':') { if (*s == '\\' && s[1]) ++s ; ++s ; }
  *len = s - text ;
  if (!*s++) return ; // just the identifier given
  if (*s == 'R' && !*++s) { *isRC = true ; return ; }
  if (*s == '-') die ("must have start coord after ':' in %s", text) ;
  *start = strtoll(s, &s, 10) ;
  if (*s++ != '-') die ("must have '-' after start coord in %s", text) ;
  if (*s == 'R' && !*++s) { *isRC = true ; return ; }
  if (!*s) return ; // just the identifier and the start given
  *end = strtoll(s, &s, 10) ;
  if (*s == 'R' && !*++s) { *isRC = true ; return ; }
  if (*s) die ("bad end in %s", text) ;
  if (*end && *start > *end) die ("start %llu > end %llu in %s", *start, *end, text) ;
}

void parseText (char *s, Frag *f)
{
  U64 len ;
  parse (s, &len, &f->start, &f->end, &f->isRC) ;
  char c = s[len] ; s[len] = 0 ; dictAdd (dict, s, &f->k) ; s[len] = c ;
}

void parseCount (char *s, Frag *f)
{
  U64 len ;
  parse (s, &len, &f->start, &f->end, &f->isRC) ;
  char *kEnd ;
  f->k = strtoll (s, &kEnd, 10) ; if (kEnd - s != len) die ("bad count in %s", s) ;
  if (f->k <= 0) die ("count must be > 0 in %s", s) ;
}

void parseFile (char *fileName, Array frags, bool isCount)
{
  FILE *f = fopen (fileName, "r") ; if (!f) die ("failed to open file %s", fileName) ;
  U64 line = 0 ;
  while (!feof (f))
    { char *s = fgetword(f) ;
      if (feof (f)) break ;
      if (isCount) parseCount (s, arrayp(frags, arrayMax(frags), Frag)) ;
      else parseText (s, arrayp (frags, arrayMax(frags), Frag)) ;
      ++line ;
      if (getc(f) != '\n' && !feof(f)) die ("bad end of line %llu in file %s", line, fileName) ;
    }
  fclose (f) ;
}

void writeFrag (SeqIO *siIn, SeqIO *siOut, Frag *f, U64 count)
{
  static char idBuf[1024] ;
  char *seq = sqioSeq(siIn) ;
  if (!f->end) f->end = siIn->seqLen ;
  if (f->isRC)
    { I64 t = f->start ; f->start = siIn->seqLen - f->end ; f->end = siIn->seqLen - t ;
      seq = seqRevComp (seq, siIn->seqLen) ;
    }      
  seq += f->start ;
  if (count)
    { if (siIn->idLen)
	{ if (snprintf (idBuf, 1024, "%s:%llu-%llu", sqioId(siIn), f->start,f->end) >= 1023)
	    die ("buffer overflow in writeFrag - please use smaller identifiers or recompile with bigger buffer") ;
	}
      else sprintf (idBuf, "%llu:%llu-%llu", count, f->start, f->end) ;
      if (f->isRC) strcat (idBuf, "R") ;
      seqIOwrite (siOut, idBuf, 0, f->end - f->start, seq, 0) ;
    }
  else
    seqIOwrite (siOut, 0, 0, f->end - f->start, seq, 0) ;
  if (f->isRC) free (seq - f->start) ;
}

int main (int argc, char *argv[])
{
  storeCommandLine (argc, argv) ;
  --argc ; ++argv ;

  if (sizeof (I64) != sizeof (long long)) die ("I64 size mismatch") ;

  Array textFrags = arrayCreate (256, Frag) ;
  Array countFrags = arrayCreate (256, Frag) ;
  dict = dictCreate (256) ;
  SeqIOtype outType = FASTA ;
  char* outFileName = "-" ;
  bool isWriteIdentifiers = false ;
  
  while (argc > 1)
    if (!strcmp (*argv, "-1")) { outType = ONE ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-I")) { isWriteIdentifiers = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-fa")) { outType = FASTA ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-o")) { outFileName = argv[1] ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-f"))
      { parseText (argv[1], arrayp(textFrags,arrayMax(textFrags),Frag)) ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-c"))
      { parseCount (argv[1], arrayp(countFrags,arrayMax(countFrags),Frag)) ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-F")) { parseFile (argv[1], textFrags, false) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-C")) { parseFile (argv[1], countFrags, true) ; argc -= 2 ; argv += 2 ; }
    else { fprintf (stderr, "unknown option %s\n", *argv) ; usage () ; }
  if (argc != 1) usage () ;

  arraySort (textFrags, fragOrder) ;
  arraySort (countFrags, fragOrder) ;
  Array kStart = arrayCreate (256, I64) ; // build an index of the start of k in textFrags
  array(kStart, 0, I64) = 0 ; // must start like this
  U64 i ;
  for (i = 1 ; i < arrayMax(textFrags) ; ++i)
    if (arrp(textFrags,i,Frag)->k > arrp(textFrags,i-1,Frag)->k) // should be +1
      array(kStart, arrayMax(kStart), I64) = i ;

  SeqIO *siIn = seqIOopenRead (*argv, dna2textConv, 0) ;
  if (!siIn) die ("failed to open input sequence file %s", *argv) ;
  SeqIO *siOut = seqIOopenWrite (outFileName, outType, dna2textConv, 0) ;
  if (!siOut) die ("failed to open output file %s", outFileName) ;
  if (outType == ONE)
    { oneAddProvenance ((OneFile*)siOut->handle, "seqextract", VERSION, getCommandLine()) ;
      oneAddReference ((OneFile*)siOut->handle, *argv, 0) ;
    }
  else
    isWriteIdentifiers = true ;
  U64 count = 0, kt = 0, kc = 0 ;
  U64 k ;
  while (seqIOread (siIn))
    { ++count ;
      while (kc < arrayMax (countFrags) && arrp(countFrags,kc,Frag)->k == count)
	writeFrag (siIn, siOut, arrp(countFrags,kc++,Frag), isWriteIdentifiers ? count : 0) ;
      if (siIn->idLen && dictFind (dict, sqioId(siIn), &k))
	for (i = arr(kStart,k,I64) ; i < arrayMax(textFrags) && arrp(textFrags,i,Frag)->k == k ; ++i)
	  { writeFrag (siIn, siOut, arrp(textFrags,i,Frag), count) ; ++kt ; }
      if (kt == arrayMax (textFrags) && kc == arrayMax (countFrags)) break ;
    }
  seqIOclose (siIn) ;
  seqIOclose (siOut) ;
  arrayDestroy (textFrags) ; arrayDestroy (countFrags) ; arrayDestroy (kStart) ; dictDestroy (dict) ;
		    
  return 0 ;
}
