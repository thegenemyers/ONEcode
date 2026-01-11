/*  File: seqconvert.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description: utility to convert between sequence formats
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 15 17:27 2025 (rd109)
 * Created: Sun Feb 17 10:23:37 2019 (rd109)
 *-------------------------------------------------------------------
 */

#include "seqio.h"
#include "ONElib.h"

static char *hocoSchemaText ;
static void hoco (char *seq, U64 *seqLen, U64 **runLengths) ;
static void writeHoco (SeqIO *siOut, U64 seqLen, U64 runLen, U64 *runLengths) ;
static void convertUnHoco (SeqIO *siIn, SeqIO *siOut) ;

static char *scaffoldSchemaText ;
static void scaffoldBreak (SeqIO *siOut, char *id, char *desc, U64 seqLen, char *seq, char *qual,
			   int scaffThresh) ; 
static void scaffoldJoin (char *fileName, SeqIO *siOut, bool isVerbose) ; 

int main (int argc, char *argv[])
{
  storeCommandLine (argc, argv) ;
  --argc ; ++argv ;

  timeUpdate (stdout) ;

  if (!argc || !strcmp(*argv,"-h") || !strcmp(*argv,"--help"))
    { fprintf (stderr, "Usage: seqconvert [-fa|fq|b|1] [-t] [-Q T] [-H|U] [-K|J] [-KT T] [-z] [-S] [-R cramRefFile] [-o outfile] [infile]\n") ;
      fprintf (stderr, "   autodetects input file type: fasta/q (.gz), binary, ONEcode, BAM/SAM\n") ;
      fprintf (stderr, "   .gz ending outfile name implies gzip compression\n") ;
      fprintf (stderr, "   -fa : output as fasta, -fq as fastq, -b as binary, -1 as ONEcode\n") ;
      fprintf (stderr, "      else .fa or .fq in outfile name imply fasta, fastq else ONEcode\n") ;
      fprintf (stderr, "   -Q  : sets the quality threshold for single bit quals in -b option [30]\n") ;
      fprintf (stderr, "   -S  : silent - else it reports to stderr on what it is doing\n") ;
      fprintf (stderr, "   -H  : homopolymer compress (hoco) - stores run lengths if ONEcode\n") ;
      fprintf (stderr, "   -U  : homopolymer uncompress - only works on ONEcode input\n") ;
      fprintf (stderr, "   -t  : show time and memory usage\n") ;
      fprintf (stderr, "   -K  : scaffold break sequences at >KT N's - stores breaks if ONEcode\n") ;
      fprintf (stderr, "   -J  : scaffold rejoin - only works on ONEcode input\n") ;
      fprintf (stderr, "   -KT : sets the threshold for scaffold breaking [20]\n") ;
      fprintf (stderr, "   -R refFileName : fasta reference file for cram\n") ;
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
  bool isScaffold = false ;
  bool isJoin = false ;
  bool isTime = false ;
  char *inFileName = "-" ;
  char *outFileName = "-z" ;
  int qualThresh = 30 ;
  int scaffThresh = 20 ;
  
  while (argc)
    { if (!strcmp (*argv, "-fa")) type = FASTA ;
      else if (!strcmp (*argv, "-fq")) type = FASTQ ;
      else if (!strcmp (*argv, "-b")) type = BINARY ;
      else if (!strcmp (*argv, "-1")) type = ONE ;
      else if (!strcmp (*argv, "-Q") && argc > 1)
	{ --argc ; ++argv ; qualThresh = atoi (*argv) ; }
      else if (!strcmp (*argv, "-z")) isGzip = true ;
      else if (!strcmp (*argv, "-H")) isHoco = true ;
      else if (!strcmp (*argv, "-U")) isUnHoco = true ;
      else if (!strcmp (*argv, "-t")) isTime = true ;
      else if (!strcmp (*argv, "-K")) isScaffold = true ;
      else if (!strcmp (*argv, "-J")) isJoin = true ;
      else if (!strcmp (*argv, "-KT") && argc > 1)
	{ --argc ; ++argv ; scaffThresh = atoi (*argv) ; }
      else if (!strcmp (*argv, "-o") && argc > 1)
	{ --argc ; ++argv ; outFileName = *argv ; }
      else if (!strcmp (*argv, "-S")) isVerbose = false ;
      else if (!strcmp (*argv, "-R") && argc > 1)
	{ --argc ; ++argv ; seqIOreferenceFileName (*argv) ; }
      else if (argc == 1 && **argv != '-') inFileName = *argv ;
      else die ("unknown option %s - run without arguments for help\n", *argv) ;
      --argc ; ++argv ;
    }

  if (isHoco && isScaffold) die ("sorry, can't do both scaffold and hoco for now") ;

  SeqIO *siOut ;
  if (!strcmp(outFileName, "-z") && !isGzip) outFileName = "-" ; /* remove 'z' */
  if (type == ONE && isHoco)
    { OneSchema *schema = oneSchemaCreateFromText (hocoSchemaText) ;
      OneFile *vf = oneFileOpenWriteNew (outFileName, schema, "seq", true, 1) ;
      oneSchemaDestroy (schema) ;
      siOut = seqIOadoptOneFile (vf, 0, qualThresh) ;
    }
  else if (type == ONE && isScaffold)
    { OneSchema *schema = oneSchemaCreateFromText (scaffoldSchemaText) ;
      OneFile *vf = oneFileOpenWriteNew (outFileName, schema, "seq", true, 1) ;
      oneSchemaDestroy (schema) ;
      if (!vf) die ("didn't open %s", outFileName) ;
      siOut = seqIOadoptOneFile (vf, 0, qualThresh) ;
      if (!siOut) die ("didn't adopt %s", outFileName) ;
    }
  else
    siOut = seqIOopenWrite (outFileName, type, 0, qualThresh) ;
  if (!siOut) die ("failed to open output file %s", outFileName) ;

  if (isJoin)
    { scaffoldJoin (inFileName, siOut, isVerbose) ;
      goto cleanup ;
    }
  
  bool isQual = ((siOut->type == BINARY && qualThresh > 0) || siOut->type == FASTQ || siOut->type == ONE) && !isHoco && !isUnHoco ;
  SeqIO *siIn = seqIOopenRead (inFileName, 0, isQual) ;
  if (!siIn) die ("failed to open input file %s", inFileName) ;
  if (isUnHoco && siIn->type != ONE) die ("can only Unhoco ONEcode files") ;
  if (isJoin && siIn->type != ONE) die ("can only reJoin ONEcode files") ;
  if (isVerbose)
    { fprintf (stderr, "reading from file type %s", seqIOtypeName[siIn->type]) ;
      if (siIn->type == BINARY || siIn->type == ONE)
	fprintf (stderr, "  with %llu sequences totLen %llu", siIn->nSeq, siIn->totSeqLen) ;
      fprintf (stderr, "\n") ;
    }

  U64 *runLengths ; // needed for hoco - need to declare here outside seqio loop
  if (isUnHoco)
    convertUnHoco (siIn, siOut) ;
  else
    while (seqIOread (siIn))
      { U64 seqLen = siIn->seqLen ;
	if (isScaffold)
	  scaffoldBreak (siOut,
			 siIn->idLen ? sqioId(siIn) : 0,
			 siIn->descLen ? sqioDesc(siIn) : 0,
			 seqLen, sqioSeq(siIn),
			 siIn->isQual ? sqioQual(siIn) : 0,
			 scaffThresh) ;
	else
	  { if (isHoco) hoco (sqioSeq(siIn), &seqLen, (siOut->type == ONE) ? &runLengths : 0) ;
	    seqIOwrite (siOut,
			siIn->idLen ? sqioId(siIn) : 0,
			siIn->descLen ? sqioDesc(siIn) : 0,
			seqLen, sqioSeq(siIn),
			siIn->isQual ? sqioQual(siIn) : 0) ;
	    if (isHoco && siOut->type == ONE) writeHoco (siOut, siIn->seqLen, seqLen, runLengths) ;
	  }
      }

  if (isVerbose)
    { fprintf (stderr, "written %llu sequences to file type %s, total length %llu, max length %llu\n",
	       siOut->nSeq, seqIOtypeName[siOut->type], siOut->totSeqLen, siOut->maxSeqLen) ;
    }

  seqIOclose (siIn) ;
  
 cleanup:
  seqIOclose (siOut) ;
  if (isTime) timeTotal (stdout) ;
}

/****************/

static void hoco (char *seq, U64 *seqLen, U64 **runLengths)
{
  static size_t runLengthSize = 0 ;
  char *s = seq, *t = seq ;
  U64  tot = *seqLen ; 
  if (runLengths)
    { if (runLengthSize < *seqLen)
	{ if (runLengthSize) free (*runLengths) ;
	  runLengthSize = *seqLen + 1024 ;
	  *runLengths = new (runLengthSize, U64) ;
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

static char *hocoSchemaText =
  "1 3 def 1 0  schema for seqconvert to hoco\n"
  ".\n"
  "P 3 seq SEQUENCE\n"
  "O S 1 3 DNA               sequence: the DNA string\n"
  "D I 1 6 STRING            id: (optional) sequence identifier\n"
  "D Q 1 6 STRING            quality: Q values (ascii string = q+33)\n"
  "D N 3 3 INT 4 CHAR 3 INT  non-acgt base\n"
  "D H 2 3 INT 8 INT_LIST    original length, list of run lengths\n" ;

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

/******** a little buffer package *********/

typedef struct { I64 size ; char *buf ; } Buffer ;

static inline Buffer *bufferCreate (I64 size)
{ Buffer *b = new (1, Buffer) ; b->size = size ; b->buf = new (size, char) ; return b ; }

static inline void bufferCheckSize (Buffer *b, int size)
{ if (size > b->size) { free (b->buf) ; b->size += size ; b->buf = new (b->size, char) ; } }

static inline void bufferDestroy (Buffer *b)
{ free (b->buf) ; free (b) ; }

/*******************************************/

static void convertUnHoco (SeqIO *siIn, SeqIO *siOut)
{
  OneFile *vf = (OneFile*) siIn->handle ;
  Buffer *sbuf = bufferCreate (1<<20) ;
  Buffer *tbuf = bufferCreate (1<<20) ;

  while (vf->lineType == 'S') // we got to an S line in sqioOpenRead()
    { I64 slen = 0, tlen = oneLen(vf) ;
      bufferCheckSize (tbuf, tlen) ;
      memcpy (tbuf->buf, oneString(vf), tlen) ; // needed because the ONElib buffer will be reused
      
      while (oneReadLine (vf) && vf->lineType != 'H' && vf->lineType != 'S')
	if (vf->lineType == 'I') storeIdLine (siIn, vf) ;
	else if (vf->lineType == 'N') tbuf->buf[oneInt(vf,0)] = oneChar(vf,1) ; // will only be one
      if (vf->lineType == 'H')
	{ assert (oneLen(vf) == tlen) ;
	  slen = oneInt(vf,0) ;
	  bufferCheckSize (sbuf, slen) ;
	  I64 *r = oneIntList(vf) ;
	  I64 n = *r ;
	  char *s = sbuf->buf, *t = tbuf->buf ;
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
		    slen, sbuf->buf, 0) ;
      else // just write the input sequence which is still in tbuf
	seqIOwrite (siOut,
		    siIn->idLen ? sqioId(siIn) : 0,
		    siIn->descLen ? sqioDesc(siIn) : 0,
		    tlen, tbuf->buf, 0) ;
    }
  
  bufferDestroy (sbuf) ;
  bufferDestroy (tbuf) ;
}

/****************/

// This breaks sequences at >= scaffThresh consecutive non-ACGT characters, and trims such characters
// off the start and end of the sequence.
// If any trimming or breaking has taken place and siOut is a ONEfile then a scaffold object (type 's')
// is created and the non-ACGT characters recorded as n, allowing reconstruction of the sequence with
// scaffoldJoin().
// Within a scaffold, if not a ONEfile, the sequences are named id.<k> where k = 1...

static char *scaffoldSchemaText =
  "1 3 def 1 0  schema for seqconvert to scafffold\n"
  ".\n"
  "P 3 seq SEQUENCE\n"
  "O s 2 3 INT 6 STRING      scaffold: length then names, made of S objects and n lines\n"
  "D g 1 3 INT               gap: length of block of n's in scaffold\n"
  "G S                       scaffolds group sequences\n"
  "O S 1 3 DNA               sequence: the DNA string\n"
  "D I 1 6 STRING            id: (optional) sequence identifier\n"
  "D Q 1 6 STRING            quality: Q values (ascii string = q+33)\n"
  "D N 3 3 INT 4 CHAR 3 INT  non-acgt base\n" ;

static void scaffoldBreak (SeqIO *siOut, char *id, char *desc, U64 seqLen, char *seq, char *qual,
			   int scaffThresh)
{
  OneFile *vf = (OneFile*) siOut->handle ;
  I64      i, firstN1 = 0 ; // firstN1 is 0 if not in a run of Ns, 1 + start of run if in a run of Ns
  bool     isOne = (siOut->type == ONE) ;
  int      k = 0 ;             // number of contig in scaffold

  if (!seqLen) return ;

  static Buffer *idBuf = 0 ;
  if (!idBuf) idBuf = bufferCreate (64) ; // NB never destroyed - small one-time memory leak
  if (id) { bufferCheckSize (idBuf, strlen(id) + 12) ; strcpy (idBuf->buf, id) ; }
  else sprintf (idBuf->buf, "s%lld", siOut->nSeq+1) ;
  char *idTail = idBuf->buf + strlen(idBuf->buf) ;
  
  if (isOne)
    { oneInt(vf,0) = seqLen ;
      oneWriteLine (vf, 's', strlen (idBuf->buf), idBuf->buf) ; // write the scaffold record
      if (desc) oneWriteComment (vf, "%s", desc) ;
    }
  
  for (i = 0 ; i < seqLen ; ++i)
    if (firstN1 && acgtCheck[(int)seq[i]]) // end of a block of non-acgt chars
      { --firstN1 ; 
	if (firstN1 == 0 || i - firstN1 >= scaffThresh)
	  { if (firstN1) // write a sequence up until firstN1
	      { ++k ;
		if (isOne)
		  seqIOwrite (siOut, 0, 0, firstN1, seq, qual) ;
		else
		  { sprintf (idTail, ".%d", k) ;
		    seqIOwrite (siOut, idBuf->buf, desc, firstN1, seq, qual) ;
		  }
	      }
	    if (isOne)
	      { oneInt(vf,0) = i - firstN1 ; oneWriteLine (vf, 'g', 0, 0) ; } // write gap record
	    seq += i ; if (qual) qual += i ; seqLen -= i ; i = 0 ;
	  }
	firstN1 = 0 ;
      }
    else if (!firstN1 && !acgtCheck[(int)seq[i]]) // start of a block of non-acgt chars
      firstN1 = i+1 ; // set the marker; +1 so that it is set if i == 0

  // now have finished checking sequence
  I64 finalNs = 0 ;
  if (firstN1--) { finalNs = seqLen - firstN1 ; seqLen = firstN1 ; }
  if (seqLen) // write a sequence up until end
    { ++k ;
      if (isOne)
	seqIOwrite (siOut, 0, 0, seqLen, seq, qual) ;
      else
	{ sprintf (idTail, ".%d", k) ;
	  seqIOwrite (siOut, idBuf->buf, desc, seqLen, seq, qual) ;
	}
    }
  if (isOne && finalNs)
    { oneInt(vf,0) = finalNs ; oneWriteLine (vf, 'g', 0, 0) ; } // write gap record
}

static void scaffoldJoin (char *inFileName, SeqIO *siOut, bool isVerbose)
{
  OneFile *vfIn = oneFileOpenRead (inFileName, 0, "seq", 1) ;
  if (!vfIn) die ("failed to open OneFile %s to read", inFileName) ;
  if (isVerbose)
    { fprintf (stderr, "reading from file type onecode") ;
      fprintf (stderr, " with %lld scaffolds",  vfIn->info['s']->given.count) ;
      fprintf (stderr, " containing %lld sequences", vfIn->info['S']->given.count) ;
      fprintf (stderr, " with total length %lld\n", vfIn->info['S']->given.total) ;
    }

  Buffer *seqBuf = bufferCreate (1 << 20) ;
  Buffer *idBuf  = bufferCreate (64) ;
  Buffer *dscBuf = bufferCreate (256) ;
  char *dsc ;
  
  if (!oneGoto (vfIn, 's', 1)) die ("can't locate to start of first scaffold") ;
  oneReadLine(vfIn) ;
  while (true)
    { if (!vfIn->lineType) break ; // end of file
      I64 scaffLen = oneInt (vfIn, 0) ;
      bufferCheckSize (seqBuf, scaffLen) ;
      char *s = seqBuf->buf, *S = s ;
      bufferCheckSize (idBuf, oneLen(vfIn)) ; strcpy (idBuf->buf, oneString(vfIn)) ;
      if ((dsc = oneReadComment (vfIn)))
	{ bufferCheckSize (dscBuf, strlen(dsc)) ; strcpy (dscBuf->buf, dsc) ; dsc = dscBuf->buf ; }
      while (oneReadLine(vfIn))
	{ if (vfIn->lineType == 's') break ;
	  if (vfIn->lineType == 'g') // gap
	    { memset (s, 'n', oneInt(vfIn,0)) ; s += oneInt(vfIn,0) ; }
	  else if (vfIn->lineType == 'S') //
	    { memcpy (s, oneDNAchar(vfIn), oneLen(vfIn)) ; S = s ; s += oneLen(vfIn) ; }
	  else if (vfIn->lineType == 'N')
	    S[oneInt(vfIn,0)] = oneChar(vfIn,1) ;
	}
      seqIOwrite (siOut, idBuf->buf, dsc, scaffLen, seqBuf->buf, 0) ;
    }
}

/****************/

