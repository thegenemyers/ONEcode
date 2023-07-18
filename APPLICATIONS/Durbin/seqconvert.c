/*  File: seqconvert.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description: utility to convert between sequence formats
 * Exported functions:
 * HISTORY:
 * Last edited: May 30 14:02 2023 (rd109)
 * Created: Sun Feb 17 10:23:37 2019 (rd109)
 *-------------------------------------------------------------------
 */

#include "seqio.h"

static void hoco (char *seq, U64 *seqLen, U64 **runLengths, U64 *runLengthSize) ;
static void writeHoco (SeqIO *si, U64 seqLen, U64 runLen, U64 *runLengths) ;
static void convertUnHoco (SeqIO *siIn, SeqIO *siOut) ;

int main (int argc, char *argv[])
{
  storeCommandLine (argc, argv) ;
  --argc ; ++argv ;

  timeUpdate (stderr) ;

  if (!argc || !strcmp(*argv,"-h") || !strcmp(*argv,"--help"))
    { fprintf (stderr, "Usage: seqconvert [-fa|fq|b|1] [-Q T] [-H|U] [-z] [-S] [-o outfile] [infile]\n") ;
      fprintf (stderr, "   autodetects input file type: fasta/q (.gz), binary, ONEcode, BAM/SAM\n") ;
      fprintf (stderr, "   .gz ending outfile name implies gzip compression\n") ;
      fprintf (stderr, "   -fa output as fasta, -fq as fastq, -b as binary, -1 as ONEcode\n") ;
      fprintf (stderr, "      else .fa or .fq in outfile name imply fasta, fastq else binary\n") ;
      fprintf (stderr, "   -Q sets the quality threshold for single bit quals in -b option [0]\n") ;
      fprintf (stderr, "   -S silent - else it reports to stderr on what it is doing\n") ;
      fprintf (stderr, "   -H homopolymer compress (hoco) - saves run lengths if ONEcode\n") ;
      fprintf (stderr, "   -U homopolymer uncompress - only works on ONEcode input\n") ;
      fprintf (stderr, "   NB gzip is not compatible with binary\n") ;
      fprintf (stderr, "   if no infile then use stdin\n") ;
      fprintf (stderr, "   if no -o option then use stdout and -z implies gzip\n");
      exit (0) ;
    }
  
  SeqIOtype type = UNKNOWN ;
  bool isVerbose = true ;
  bool isGzip = false ;
  bool isHoco = false ;
  bool isUnHoco = false ;
  char *inFileName = "-" ;
  char *outFileName = "-z" ;
  int qualThresh = 0 ;
  while (argc)
    { if (!strcmp (*argv, "-fa")) type = FASTA ;
      else if (!strcmp (*argv, "-fq")) type = FASTQ ;
      else if (!strcmp (*argv, "-b")) type = BINARY ;
      else if (!strcmp (*argv, "-1")) type = ONE ;
      else if (!strcmp (*argv, "-Q") && argc >1)
	{ --argc ; ++argv ; qualThresh = atoi (*argv) ; }
      else if (!strcmp (*argv, "-z")) isGzip = true ;
      else if (!strcmp (*argv, "-H")) isHoco = true ;
      else if (!strcmp (*argv, "-U")) isUnHoco = true ;
      else if (!strcmp (*argv, "-o") && argc > 1)
	{ --argc ; ++argv ; outFileName = *argv ; }
      else if (!strcmp (*argv, "-S")) isVerbose = false ;
      else if (argc == 1 && **argv != '-') inFileName = *argv ;
      else die ("unknown option %s - run without arguments for help\n", *argv) ;
      --argc ; ++argv ;
    }

  if (!strcmp(outFileName, "-z") && !isGzip) outFileName = "-" ; /* remove 'z' */
  SeqIO *siOut = seqIOopenWrite (outFileName, type, 0, qualThresh) ;
  if (!siOut) die ("failed to open output file %s", outFileName) ;
  bool isQual = ((siOut->type == BINARY && qualThresh > 0) || siOut->type == FASTQ || siOut->type == ONE) && !isHoco && !isUnHoco ;
  SeqIO *siIn = seqIOopenRead (inFileName, 0, isQual) ;
  if (!siIn) die ("failed to open input file %s", inFileName) ;
  if (isUnHoco && siIn->type != ONE) die ("can only Unhoco ONEcode files") ;
  if (isVerbose)
    { fprintf (stderr, "reading from file type %s", seqIOtypeName[siIn->type]) ;
      if (siIn->type == BINARY || siIn->type == ONE)
	fprintf (stderr, "  with %" PRIu64 " sequences totLen %" PRIu64 "", siIn->nSeq, siIn->totSeqLen) ;
      fprintf (stderr, "\n") ;
    }

  if (isUnHoco)
    convertUnHoco (siIn, siOut) ;
  else
    while (seqIOread (siIn))
      { U64 seqLen = siIn->seqLen ;
	U64 *runLengths, runLengthSize = 0 ;
	if (isHoco) hoco (sqioSeq(siIn), &seqLen, (siOut->type == ONE) ? &runLengths : 0, &runLengthSize) ;
	seqIOwrite (siOut,
		    siIn->idLen ? sqioId(siIn) : 0,
		    siIn->descLen ? sqioDesc(siIn) : 0,
		    seqLen, sqioSeq(siIn),
		    siIn->isQual ? sqioQual(siIn) : 0) ;
	if (isHoco && siOut->type == ONE) writeHoco (siOut, siIn->seqLen, seqLen, runLengths) ;
      }

  if (isVerbose)
    { fprintf (stderr, "written %" PRIu64 " sequences to file type %s, total length %" PRIu64 ", max length %" PRIu64 "\n",
	       siOut->nSeq, seqIOtypeName[siOut->type], siOut->totSeqLen, siOut->maxSeqLen) ;
  
      timeTotal (stderr) ;
    }

  seqIOclose (siIn) ;
  seqIOclose (siOut) ;
}

/****************/

static void hoco (char *seq, U64 *seqLen, U64 **runLengths, U64 *runLengthSize)
{
  char *s = seq, *t = seq ;
  U64  tot = *seqLen ; 
  if (runLengths)
    { if (*runLengthSize < *seqLen)
	{ if (*runLengthSize) free (*runLengths) ;
	  *runLengthSize = *seqLen + 1024 ;
	  *runLengths = new (*runLengthSize, U64) ;
	}
      U64 *r = *runLengths ;
      *t = *s++ ; --tot ; *r = 1 ;
      while (tot--)
	if (*s == *t) { ++s ; ++*r ; }
	else { *++t = *s++ ; r[1] = 1 + *r ; ++r ; }
    }
  else // duplicate code to avoid extra tests inside the loop
    { *t = *s++ ; --tot ;
      while (tot--)
	if (*s == *t) ++s ; 
	else *++t = *s++ ;
    }
  *seqLen = (t - seq) + 1 ;
}

#include "ONElib.h"

static void writeHoco (SeqIO *si, U64 seqLen, U64 runLen, U64 *runLengths)
{
  OneFile *vf = (OneFile*) si->handle ;
  oneInt(vf,0) = seqLen ;
  oneWriteLine (vf, 'H', runLen, runLengths) ;
}


static void storeIdLine (SeqIO *si, OneFile *vf)
{ char *desc = oneReadComment (vf) ;
  si->idLen = oneLen(vf) ;
  if (desc) si->descLen = strlen(desc) ; else si->descLen = 0 ;
  si->idStart = 0 ; strcpy (si->buf, oneString(vf)) ;
  si->descStart = si->idLen+1 ;
  if (desc) strcpy (si->buf+si->descStart, desc) ; else si->buf[si->descStart] = 0 ;
}


static void convertUnHoco (SeqIO *siIn, SeqIO *siOut)
{
  OneFile *vf = (OneFile*) siIn->handle ;
  char *sbuf, *tbuf ;
  I64 sbufSize = 0, tbufSize = 0 ;

  while (vf->lineType == 'S') // we got to an S line in sqioOpenRead()
    { I64 slen = 0, tlen = oneLen(vf) ;
      if (tbufSize < tlen)
	{ if (tbufSize) free (tbuf) ; tbufSize = tlen+4096 ; tbuf = new (tbufSize, char) ; }
      memcpy (tbuf, oneString(vf), tlen) ; // need to make a copy, because the ONElib buffer will be reused
      
      while (oneReadLine (vf) && vf->lineType != 'H' && vf->lineType != 'S')
	if (vf->lineType == 'I') storeIdLine (siIn, vf) ;
      if (vf->lineType == 'H')
	{ assert (oneLen(vf) == tlen) ;
	  slen = oneInt(vf,0) ;
	  if (sbufSize < slen)
	    { if (sbufSize) free (sbuf) ; sbufSize = slen+4096 ; sbuf = new (sbufSize, char) ; }
	  I64 *r = oneIntList(vf) ;
	  I64 n = *r ;
	  char *s = sbuf, *t = tbuf ;
	  while (tlen--)
	    { while (n--) *s++ = *t ;
	      n = r[1] - *r ; ++r ; ++t ;
	    }
	  while (oneReadLine (vf) && vf->lineType != 'S')
	    if (vf->lineType == 'I') storeIdLine (siIn, vf) ;
	}
      if (slen) // then I read a hoco line and made the expanded sequence in sbuf
	seqIOwrite (siOut,
		    siIn->idLen ? sqioId(siIn) : 0,
		    siIn->descLen ? sqioDesc(siIn) : 0,
		    slen, sbuf, 0) ;
      else // just write the input sequence which is still in tbuf
	seqIOwrite (siOut,
		    siIn->idLen ? sqioId(siIn) : 0,
		    siIn->descLen ? sqioDesc(siIn) : 0,
		    tlen, tbuf, 0) ;
    }
  
  if (sbufSize) free (sbuf) ;
  if (tbufSize) free (tbuf) ;
}

/****************/
