/*  File: seqsample.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2026
 *-------------------------------------------------------------------
 * Description: to sample reads from sequences (potentially a genome)
 * Exported functions:
 * HISTORY:
 * Last edited: May 29 20:40 2026 (rd109)
 * Created: Wed May 27 16:50:55 2026 (rd109)
 *-------------------------------------------------------------------
 */

#include "seqio.h"
#include <math.h>

static int length = 10000 ;
static double coverage = 10.0 ;

void usage (void)
{
  fprintf (stderr, "Usage: seqsample [opts] <filename> [> output.fa]\n") ;
  fprintf (stderr, "  will read [gzipped] fasta, fastq, or bam/sam/cram, 1code.  Use filename '-' for stdin (not 1code binary)\n") ;
  fprintf (stderr, "  options:\n") ;
  fprintf (stderr, "    -L <int> : mean length [%d]\n", length) ;
  fprintf (stderr, "    -C <float> : mean coverage [%f]\n", coverage) ;
  fprintf (stderr, "    -e <float> : SNV error rate [0.0] - use 0.01 for 1%%\n") ;
  fprintf (stderr, "    -d <float> : single base indel error rate [0.0] - use 0.01 for 1%%\n") ;
  fprintf (stderr, "    -o <filename> : output file [-]\n") ;
  fprintf (stderr, "        if filename ends .fa.gz then gzipped fasta\n") ;
  fprintf (stderr, "        if filename ends .1* then writes ONEcode with basic 1seq schema\n") ;
  fprintf (stderr, "    -z : output gzipped fasta to standard output\n") ;
  fprintf (stderr, "    -1 : output 1seq to standard output\n") ;
  exit (0) ;
}

double sampleExponential (double mean) // 
{
  static double uFac = 1.0 / ((I64)1 << 31) ;
  long x = random () ;
  double U = uFac * (1 + x) ; // a random variate between 0 and 1
  double y = - mean * log(U) ;
  // printf ("** sampleExponential uFac %e x %ld U %f y %f\n", uFac, x, U, y) ;
  return y ;
}

char sampleBase (void)
{
  static char acgt[] = "acgt" ;
  return acgt[random() & 0x03] ;
}

bool sampleBool (void)
{ return (random() & 0x01) == 0 ;
}

int main (int argc, char *argv[])
{
  SeqIOtype type = FASTA ;
  double errSNV = 0.0 ;
  double errIndel = 0.0 ;
  char  *outFileName = "-" ;
  U32    seed = 7 ;

  timeUpdate (0) ;
  storeCommandLine (argc, argv) ;
  --argc ; ++argv ;
  if (!argc) usage() ;
  while (argc >= 2 && **argv == '-')
    if (!strcmp (*argv, "-L")) { length = atoi (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-C")) { coverage = atof (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-e")) { errSNV = atof (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-d")) { errIndel = atof (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-o")) { outFileName = argv[1] ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-z")) { outFileName = "-z" ; argc-- ; argv++ ; }
    else if (!strcmp (*argv, "-1")) { type = ONE ; argc-- ; argv++ ; }
    else { warn ("unknown option %s", *argv) ; usage() ; }

  srandom (seed) ;
  
  if (coverage <= 0) die ("coverage %f must be positive", coverage) ;
  double spacing = length / coverage ;
  if (errSNV < 0 || errSNV >= 1.0) die ("errSNV %f must be >= 0 and less than 1", errSNV) ;
  double snvSpacing = errSNV ? 1.0 / errSNV : 0.0 ;
  if (errIndel < 0 || errIndel >= 1.0) die ("errSNV %f must be >= 0 and less than 1", errIndel) ;
  double indelSpacing = errIndel ? 1.0 / errIndel : 0.0 ;
  
  if (!argc) die ("no sequence argument given") ;
  if (argc > 1) die ("more than one sequence argument - first would be %s", *argv) ;

  SeqIO *siIn = seqIOopenRead (*argv, 0, false) ;
  if (!siIn) die ("failed to open sequence file %s", *argv) ;

  { char *s = outFileName + strlen(outFileName) ;
    while (s > outFileName) if (*--s == '.') { if (s[1] == '1') type = ONE ; break ; }
  }
  SeqIO *siOut = seqIOopenWrite (outFileName, type, 0, false) ;
  if (!siOut) die ("failed to open output file %s", outFileName) ;

  int   idMax = 1024 ;
  char *idBuf = new (idMax, char) ;
  while (seqIOread (siIn))
    { int start = 0 ;
      int idLen = 0 ;
      if (sqioId(siIn))
	{ idLen = strlen(sqioId(siIn)) ;
	  while (idLen+11 > idMax)
	    { newFree (idBuf, idMax, char) ; idMax *= 2 ; idBuf = new(idMax, char) ; }
	  strcpy (idBuf, sqioId(siIn)) ;
	}
      while (true)
	{ start += sampleExponential (spacing) ;
	  int len = length ; // could add uncertainty here
	  if (len > (1<<30)) die ("len %'d > 1Gb - too big", len) ;
	  if (start + len > siIn->seqLen) break ;
	  char *seq = sqioSeq(siIn) + start ;
	  sprintf (idBuf+idLen, "+%d", start) ;
	  if (errSNV)
	    { int off = sampleExponential (snvSpacing) ;
	      while (off < len)
		{ seq[start+off] = sampleBase () ;
		  off += sampleExponential (snvSpacing) ;
		}
	    }
	  if (errIndel) // NB keep length
	    { int off = sampleExponential (indelSpacing) ;
	      while (off < len)
		{ if (sampleBool()) // delete
		    if (start + len == siIn->seqLen)
		      { memmove (seq+off, seq+off+1, len-off-1) ; seq[len-1] = sampleBase() ; }
		    else memmove (seq+off, seq+off+1, len-off) ;
		  else // insert
		    { memmove (seq+off+1, seq+off, len-off-1) ; seq[off] = sampleBase() ; }
		  off += sampleExponential (indelSpacing) ;
		}
	    }
	  if (sampleBool()) { seqRevComp (seq, len) ; idBuf[idLen] = '-' ; }
	  seqIOwrite (siOut, idBuf, 0, len, seq, 0) ;
	}
    }
  fprintf (stderr, "read %'llu sequences, total length %'llu\n", siIn->nSeq, siIn->totSeqLen) ;
  seqIOclose (siIn) ;
  fprintf (stderr, "wrote %'llu sequences, total length %'llu", siOut->nSeq, siOut->totSeqLen) ;
  if (siOut->nSeq) fprintf (stderr, ", average %.1f\n", siOut->totSeqLen / (double)siOut->nSeq) ;
  seqIOclose (siOut) ;
  timeTotal (stderr) ;
}
