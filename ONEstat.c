/*****************************************************************************************
 *
 *  File: ONEstat.c
 *    one format validator and header generator
 *
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *
 * HISTORY:
 * Last edited: Sep 19 16:13 2024 (rd109)
 *   * Dec 27 09:20 2019 (gene): style edits
 *   * Created: Thu Feb 21 22:40:28 2019 (rd109)
 *
 ****************************************************************************************/

#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "ONElib.h"
extern void oneFinalizeCounts (OneFile *vf) ; // secret connection into ONElib.c for checking

// forward declarations for utilities at the end of the file
void die (char *format, ...) ;
void timeUpdate (FILE *f) ;
void timeTotal (FILE *f) ;

int main (int argc, char **argv)
{ int        i ;
  char      *fileType = 0 ;
  char      *outFileName = "-" ;
  bool       isHeader = false, isUsage = false, isVerbose = false ;
  char      *schemaFileName = 0 ;
  char      *checkText = 0 ;
  
  timeUpdate (0) ;

  //  Process command line arguments

  --argc ; ++argv ;		/* drop the program name */
  if (argc == 0)
    { fprintf (stderr, "ONEstat [options] onefile\n") ;
      fprintf (stderr, "  -t --type <abc>          file type, e.g. seq - required if no header\n") ;
      fprintf (stderr, "  -S --schema <schema>     schema file - required if not in file\n") ;
      fprintf (stderr, "  -C --check 'schematext'  check for a limited set of features\n") ;
      fprintf (stderr, "  -H --header              output header accumulated from data\n") ;
      fprintf (stderr, "  -o --output <filename>   output to filename\n") ;
      fprintf (stderr, "  -u --usage               byte usage per line type; no other output\n") ;
      fprintf (stderr, "  -v --verbose             else only errors and requested output\n") ;
      fprintf (stderr, "ONEstat aborts on a syntactic parse error with a message.\n") ;
      fprintf (stderr, "Otherwise information is written to stderr about any inconsistencies\n") ;
      fprintf (stderr, "between the header and the data in the body of the file.\n") ;
      fprintf (stderr, "Output is to stdout by default, use -o to overide\n");
      exit (0) ;
    }
  
  while (argc && **argv == '-')
    if (!strcmp (*argv, "-H") || !strcmp (*argv, "--header"))
      { isHeader = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-u") || !strcmp (*argv, "--usage"))
      { isUsage = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-v") || !strcmp (*argv, "--verbose"))
      { isVerbose = true ; --argc ; ++argv ; }
    else if (argc > 1 && (!strcmp (*argv, "-t") || !strcmp (*argv, "--type")))
      { fileType = argv[1] ;
	argc -= 2 ; argv += 2 ;
      }
    else if (argc > 1 && (!strcmp (*argv, "-S") || !strcmp (*argv, "--schema")))
      { schemaFileName = argv[1] ;
	argc -= 2 ; argv += 2 ;
      }
    else if (argc > 1 && (!strcmp (*argv, "-C") || !strcmp (*argv, "--check")))
      { checkText = argv[1] ;
	argc -= 2 ; argv += 2 ;
      }
    else if (argc > 1 && (!strcmp (*argv, "-o") || !strcmp (*argv, "--output")))
      { outFileName = argv[1] ;
	argc -= 2 ; argv += 2 ;
      }
    else die ("unknown option %s - run without arguments to see options", *argv) ;
  
  if (argc != 1)
    die ("need to give a single data file as argument") ;

  //  Open subject file for reading and read header (if present)
  OneSchema *vs = 0 ;
  if (schemaFileName)
    { vs = oneSchemaCreateFromFile (schemaFileName) ;
      if (!vs) die ("failed to read schema file %s", schemaFileName) ;
    }
  OneFile *vf = oneFileOpenRead (argv[0], vs, fileType, 1) ;
  if (!vf) die ("failed to open OneFile %s", argv[0]) ;
  oneSchemaDestroy (vs) ; // no longer needed

  if (isVerbose)
    { if (vf->line == 1)
	fprintf (stderr, "header missing\n") ;
      else
	fprintf (stderr, "read %lld header lines\n", vf->line) ;
    }

  if (checkText)
    oneFileCheckSchemaText (vf, checkText) ;

  vf->isCheckString = true ;

  // if requesting usage, then 

  if (isUsage)
    { I64 usage[128] ; memset (usage, 0, 128*sizeof(I64)) ; 
      off_t u, uLast = ftello (vf->f) ;

      while (oneReadLine (vf))
	{ u = ftello (vf->f) ; usage[(int)vf->lineType] += u-uLast ; uLast = u ; }
      u = ftello (vf->f) ; usage[(int)vf->lineType] += u-uLast ; uLast = u ;

      FILE *f ;
      if (strcmp (outFileName, "-") && !(f = fopen (outFileName, "w")))
	die ("failed to open output file %s", outFileName) ;
      else
	f = stdout ;
      
      for (i = 'A' ; i < 128 ; ++i)
	if (usage[i]) fprintf (f, "usage line type %c bytes %lld\n", (char)i,  usage[i]) ;

      if (f != stdout) fclose (f) ;
     }

  else
    {
      //  Read data portion of file checking syntax and group sizes (if present)

      while (oneReadLine (vf)) ;

      if (isVerbose)
	fprintf (stderr, "read %lld lines from OneFile %s type %s\n",
		 vf->line, argv[0], vf->fileType) ;

      oneFinalizeCounts (vf) ;
    
      //  Check count statistics for each line type versus those in header (if was present)

      { I64 nTotal = 0, nBad = 0, nMissing = 0 ;
	  
#define CHECK(X,Y,Z)							                       \
  if (li->X > 0 && li->X != li->Y)				                               \
    { fprintf (stderr, "header mismatch %s %c: header %lld data %lld\n", Z, i,  li->X,  li->Y) ; \
      nBad += 1 ;							                       \
   } 											       \
 else if (li->Y > 0 && li->X == 0)							       \
   { fprintf (stderr, "header %s line missing for %c, value is %lld\n", Z, i,  li->Y) ;	       \
     nMissing += 1 ;									       \
   } 											       \
 if (li->Y > 0)										       \
   nTotal += 1 ;

	for (i = 0; i < 128; i++)
	  if (((i >= 'A' && i <= 'Z') || (i >= 'a' && i <= 'z')) && vf->info[i] != NULL)
	    { OneInfo *li = vf->info[i] ;
	      CHECK(given.count, accum.count, "count") ;
	      CHECK(given.max, accum.max, "max") ;
	      CHECK(given.total, accum.total, "total") ;
	  }
	if (isVerbose || nBad || nMissing)
	  fprintf (stderr, "expected %lld header content lines, of which %lld bad and %lld missing\n",
		    nTotal,  nBad,  nMissing) ;
      }

  //  Write header if requested - achieved by writing a "from" file with no lines

      if (isHeader)
	{ OneFile *vfOut = oneFileOpenWriteFrom (outFileName, vf, false, 1) ;
	  if (vfOut == NULL)
	    die ("failed to open output file %s", outFileName) ;

	  for (i = 0 ; i < 128 ; i++) // transfer accumulated counts for vfIn to given for vfOut
	    if (vfOut->info[i])
	      vfOut->info[i]->given = vf->info[i]->accum ;
  
	  fflush (vfOut->f) ;
	  oneFileClose (vfOut) ;
	}
    }

  oneFileClose (vf) ;

  if (isVerbose) timeTotal (stderr) ;

  exit (0) ;
}

/********************* utilities *************************/


void die (char *format, ...)
{
  va_list args ;

  va_start (args, format) ;
  fprintf (stderr, "FATAL ERROR: ") ;
  vfprintf (stderr, format, args) ;
  fprintf (stderr, "\n") ;
  va_end (args) ;

  exit (-1) ;
}

/***************** rusage for timing information ******************/

#include <sys/resource.h>
#ifndef RUSAGE_SELF     /* to prevent "RUSAGE_SELF redefined" gcc warning, fixme if this is more intricate */
#define RUSAGE_SELF 0
#endif

#ifdef RUSAGE_STRUCTURE_DEFINITIONS
struct rusage {
  struct timeval ru_utime; /* user time used */
  struct timeval ru_stime; /* system time used */
  long ru_maxrss;          /* integral max resident set size */
  long ru_ixrss;           /* integral shared text memory size */
  long ru_idrss;           /* integral unshared data size */
  long ru_isrss;           /* integral unshared stack size */
  long ru_minflt;          /* page reclaims */
  long ru_majflt;          /* page faults */
  long ru_nswap;           /* swaps */
  long ru_inblock;         /* block input operations */
  long ru_oublock;         /* block output operations */
  long ru_msgsnd;          /* messages sent */
  long ru_msgrcv;          /* messages received */
  long ru_nsignals;        /* signals received */
  long ru_nvcsw;           /* voluntary context switches */
  long ru_nivcsw;          /* involuntary context switches */
};

struct timeval {
  time_t       tv_sec;   /* seconds since Jan. 1, 1970 */
  suseconds_t  tv_usec;  /* and microseconds */
} ;
#endif /* RUSAGE STRUCTURE_DEFINITIONS */

static struct rusage rOld, rFirst ;

void timeUpdate (FILE *f)
{
  static bool isFirst = 1 ;
  struct rusage rNew ;
  int secs, usecs ;

  getrusage (RUSAGE_SELF, &rNew) ;
  if (!isFirst)
    { secs = rNew.ru_utime.tv_sec - rOld.ru_utime.tv_sec ;
      usecs =  rNew.ru_utime.tv_usec - rOld.ru_utime.tv_usec ;
      if (usecs < 0) { usecs += 1000000 ; secs -= 1 ; }
      fprintf (f, "user\t%d.%06d", secs, usecs) ;
      secs = rNew.ru_stime.tv_sec - rOld.ru_stime.tv_sec ;
      usecs =  rNew.ru_stime.tv_usec - rOld.ru_stime.tv_usec ;
      if (usecs < 0) { usecs += 1000000 ; secs -= 1 ; }
      fprintf (f, "\tsystem\t%d.%06d", secs, usecs) ;
      fprintf (f, "\tmax_RSS\t%ld", rNew.ru_maxrss - rOld.ru_maxrss) ;
      fputc ('\n', f) ;
    }
  else
    { rFirst = rNew ;
      isFirst = false ;
    }

  rOld = rNew ;
}

void timeTotal (FILE *f) { rOld = rFirst ; timeUpdate (f) ; }

/********************* end of file ***********************/
