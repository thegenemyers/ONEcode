/*  File: ONEview.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Apr 30 23:43 2024 (rd109)
 * Created: Thu Feb 21 22:40:28 2019 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "ONElib.h"

#include <assert.h>

#include <string.h>		/* strcmp etc. */
#include <stdlib.h>		/* for exit() */

static char *commandLine (int argc, char **argv)
{
  int i, totLen = 0 ;
  for (i = 0 ; i < argc ; ++i) totLen += 1 + strlen(argv[i]) ;
  char *buf = new (totLen, char) ;
  strcpy (buf, argv[0]) ;
  for (i = 1 ; i < argc ; ++i) { strcat (buf, " ") ; strcat (buf, argv[i]) ; }
  return buf ;
}

typedef struct IndexListStruct {
  I64 i0, iN ;
  struct IndexListStruct *next ;
} IndexList ;

static IndexList *parseIndexList (char *s)
{
  IndexList *ol, *ol0 = ol = new0 (1, IndexList) ;
  while (*s)
    { while (*s >= '0' && *s <= '9') ol->i0 = ol->i0*10 + (*s++ - '0') ;
      if (*s == '-')
	{ ++s ; while (*s >= '0' && *s <= '9') ol->iN = ol->iN*10 + (*s++ - '0') ;
	  if (ol->iN <= ol->i0) die ("end index %lld <= start index %lld", ol->iN, ol->i0) ;
	}
      else
	ol->iN = ol->i0 + 1 ;
      if (*s == ',') { ol->next = new0 (1, IndexList) ; ol = ol->next ; ++s ; }
      else if (*s) die ("unrecognised character %c at %s in object list\n", *s, s) ;
    }
  return ol0 ; 
}

static void transferLine (OneFile *vfIn, OneFile *vfOut, size_t *fieldSize)
{ memcpy (vfOut->field, vfIn->field, fieldSize[(int)vfIn->lineType]) ;
  oneWriteLine (vfOut, vfIn->lineType, oneLen(vfIn), oneString(vfIn)) ;
  char *s = oneReadComment (vfIn) ; if (s) oneWriteComment (vfOut, "%s", s) ;
}

int main (int argc, char **argv)
{
  I64 i ;
  char *fileType = 0 ;
  char *outFileName = "-" ;
  char *schemaFileName = 0 ;
  bool  isNoHeader = false, isHeaderOnly = false, isWriteSchema = false, 
    isBinary = false, isVerbose = false ;
  char  indexType ;
  IndexList *objList = 0 ;
  
  timeUpdate (0) ;

  char *command = commandLine (argc, argv) ;
  --argc ; ++argv ;		/* drop the program name */

  if (!argc)
    { fprintf (stderr, "ONEview [options] onefile\n") ;
      fprintf (stderr, "  -t --type <abc>           file type, e.g. seq, aln - required if no header\n") ;
      fprintf (stderr, "  -S --schema <schemafile>      schema file name for reading file\n") ;
      fprintf (stderr, "  -h --noHeader                 skip the header in ascii output\n") ;
      fprintf (stderr, "  -H --headerOnly               only write the header (in ascii)\n") ;
      fprintf (stderr, "  -s --writeSchema              write a schema file based on this file\n") ;
      fprintf (stderr, "  -b --binary                   write in binary (default is ascii)\n") ;
      fprintf (stderr, "  -o --output <filename>        output file name (default stdout)\n") ;
      fprintf (stderr, "  -i --index T x[-y](,x[-y])*   write specified objects/groups of type T\n") ;
      fprintf (stderr, "  -v --verbose                  write commentary including timing\n") ;
      fprintf (stderr, "index only works for binary files; '-i A 0-10' outputs first 10 objects of type A\n") ;
      exit (0) ;
    }
  
  while (argc && **argv == '-')
    if (!strcmp (*argv, "-t") || !strcmp (*argv, "--type"))
      { fileType = argv[1] ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-S") || !strcmp (*argv, "--schema"))
      { schemaFileName = argv[1] ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-h") || !strcmp (*argv, "--noHeader"))
      { isNoHeader = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-H") || !strcmp (*argv, "--headerOnly"))
      { isHeaderOnly = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-s") || !strcmp (*argv, "--writeSchema"))
      { isWriteSchema = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-b") || !strcmp (*argv, "--binary"))
      { isBinary = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-v") || !strcmp (*argv, "--verbose"))
      { isVerbose = true ; --argc ; ++argv ; }
    else if ((!strcmp (*argv, "-o") || !strcmp (*argv, "--output")) && argc >= 2)
      { outFileName = argv[1] ; argc -= 2 ; argv += 2 ; }
    else if ((!strcmp (*argv, "-i") || !strcmp (*argv, "--index")) && argc >= 3)
      { indexType = *argv[1] ; objList = parseIndexList (argv[2]) ; argc -= 3 ; argv += 3 ; }
    else die ("unknown option %s - run without arguments to see options", *argv) ;

  if (isBinary) isNoHeader = false ;
  if (isHeaderOnly) isBinary = false ;
    
  if (argc != 1)
    die ("need a single data one-code file as argument") ;

  OneSchema *vs = 0 ;
  if (schemaFileName && !(vs = oneSchemaCreateFromFile (schemaFileName)))
    die ("failed to read schema file %s", schemaFileName) ;
  OneFile *vfIn = oneFileOpenRead (argv[0], vs, fileType, 1) ; /* reads the header */
  if (!vfIn) die ("failed to open one file %s", argv[0]) ;

  if (objList)
    { if (!vfIn->isBinary)
	die ("%s is ascii - you can only access objects and groups by index in binary files", argv[0]) ;
      if (!vfIn->info[(int)indexType])
	die ("requested index type %c is not present in the schema", indexType) ;
      if (!vfIn->info[(int)indexType]->index)
	die ("no index for line type %c", indexType) ;
    }

  if (isWriteSchema)
    { oneFileWriteSchema (vfIn, outFileName) ; }
  else
    { OneFile *vfOut = oneFileOpenWriteFrom (outFileName, vfIn, isBinary, 1) ;
      if (!vfOut) die ("failed to open output file %s", outFileName) ;

      if (isNoHeader) vfOut->isNoAsciiHeader = true ; // will have no effect if binary

      if (!isHeaderOnly)
	{ oneAddProvenance (vfOut, "ONEview", "0.0", command) ;
      
	  static size_t fieldSize[128] ;
	  for (i = 0 ; i < 128 ; ++i)
	    if (vfIn->info[i]) fieldSize[i] = vfIn->info[i]->nField*sizeof(OneField) ;
      
	  if (objList)
	    while (objList)
	      { if (!oneGoto (vfIn, indexType, objList->i0))
		  die ("can't locate to object %c %lld", indexType, objList->i0 ) ;
		if (!oneReadLine (vfIn))
		  die ("can't read object %c %lld", indexType, objList->i0) ;
		if (objList->i0 == 0 && vfIn->lineType == indexType) ++objList->i0 ;
		while (vfIn->lineType && objList->i0 < objList->iN) // lineType 0 is end of file
		  { transferLine (vfIn, vfOut, fieldSize) ;
		    oneReadLine (vfIn) ;
		    if (vfIn->lineType == '/' && oneChar(vfIn,0) == indexType) // end of object
		      while (oneReadLine (vfIn) && vfIn->lineType != indexType) ;
		    if (vfIn->lineType == indexType) ++objList->i0 ;
		  }
		objList = objList->next ;
	      }
	  else
	    while (oneReadLine (vfIn))
	      transferLine (vfIn, vfOut, fieldSize) ;
	}
      oneFileClose (vfOut) ;
    }
      
  oneFileClose (vfIn) ;
  if (vs) oneSchemaDestroy (vs) ;
  
  free (command) ;
  if (isVerbose) timeTotal (stderr) ;

  exit (0) ;
}

/********************* end of file ***********************/
