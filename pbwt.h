/*  File: pbwt.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2013-
 *-------------------------------------------------------------------
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *-------------------------------------------------------------------
 * Description: header file for pbwt package
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 14 13:53 2015 (rd)
 * added paintSparse function
 * Created: Thu Apr  4 11:02:39 2013 (rd)
 *-------------------------------------------------------------------
 */

#include "utils.h"

const static int pbwtMajorVersion = 3, pbwtMinorVersion = 0 ;

const char *pbwtCommitHash(void);
const char *pbwtHtslibVersionString(void);

/* data types */

typedef unsigned char uchar ;

class PBWT {
	public:
		int N ;			/* number of sites */
		int M ;			/* number of samples */
		char* chrom ;			/* chromosome name */
		Array sites ;			/* array of Site */
		static Array samples ;		/* array of int index into global samples */
		Array yz ;			/* compressed PBWT array of uchar */
		int *aFstart, *aFend ;	/* start and end a[] index arrays for forwards cursor */
		Array zz ;			/* compressed reverse PBWT array of uchar */
		int *aRstart, *aRend ; /* start and end a[] index arrays for reverse cursor */
		/* NB aRend is the lexicographic sort order for the data, and aFend the reverse lex order */
		/* probably it is optimal to have aFstart == aRend and vice versa: to be done */
		Array zMissing ;		/* compressed array of uchar - natural not sort order */
		Array missingOffset ;		/* of long, site index into zMissing, 0 if no missing data at site */
		Array zDosage ;		/* run-length compressed array of uchar in local sort order */
		Array dosageOffset ;		/* of long, site index into zDosage, 0 if no dosage data */
		BOOL  isRefFreq ;		/* some flags for the whole VCF */
		BOOL  isUnphased ;

		/* philosophy is to be lazy about PBWT - only fill items for which we have info */

		struct Site {
			int x ;			/* position on chromosome */
			int varD ;			/* index in variationDict */
			double freq ;			/* frequency */
			double refFreq ;		/* frequency from reference used for last phasing or imputation */
			double imputeInfo ;		/* estimated r^2 from imputation */
		};

		struct Sample {
			int nameD ;			/* index in sampleDict */
			int father ;			/* index into samples */
			int mother ;			/* index into samples */
			int popD ;			/* index in populationDict */
			BOOL isMale ;			/* treat X chromosome as haploid */
			BOOL isFemale ;		/* treat X chromosome as diploid and ignore Y */
		};


		/* pbwtMain.c */

		//extern char *commandLine ;	/* a copy of the command line */
		//extern FILE *logFile ;  /* log file pointer */

		static char * commandLine;
		static FILE * logFile;

		/* pbwtCore.c */

		//extern BOOL isCheck ;		/* when TRUE carry out various checks */
		//extern BOOL isStats ;		/* when TRUE report stats in various places */
		//extern DICT *variationDict ;	/* "xxx|yyy" where variation is from xxx to yyy in VCF */
		static BOOL isCheck;
		static BOOL isStats;
		static DICT *variationDict;
		/* NB using a global DICT for variation means that identical variations use the same string */

		static void pbwtInit (void) ;
		static PBWT *pbwtCreate (int M, int N) ; /* OK to have N == 0 and set p->N later if not known now */
		static void pbwtDestroy (PBWT *p) ;
		PBWT *pbwtSubSites (PBWT *pOld, double fmin, double frac) ;
		PBWT *pbwtSubRange (PBWT *pOld, int start, int end) ;
		static void pbwtBuildReverse (PBWT *p) ;
		static uchar **pbwtHaplotypes (PBWT *p) ;
		PBWT *pbwtSelectSites (PBWT *pOld, Array sites, BOOL isKeepOld) ;
		PBWT *pbwtRemoveSites (PBWT *pOld, Array sites, BOOL isKeepOld) ;

		/* operations to move forwards and backwards in the pbwt using the cursor structure */

		struct PbwtCursor {		/* data structure for moving forwards - doesn't know PBWT */
			int M ;
			Array z ;			/* packed byte array; if zero y needs loading from elsewhere */
			long n ;			/* position in packed byte array */
			BOOL isBlockEnd ;		/* TRUE if n is at end of next block, FALSE if at start */
			uchar *y ;			/* current value in sort order */
			int c ;			/* number of 0s in y */
			int *a ;			/* index back to original order */
			int *d ;			/* location of last match */
			int *u ;			/* number of 0s up to and including this position */
			int *b ;			/* for local operations - no long term meaning */
			int *e ;			/* for local operations - no long term meaning */
			long nBlockStart ;		/* u->n at start of block encoding current u->y */
		};

			static PbwtCursor *pbwtCursorCreate (PBWT *p, BOOL isForwards, BOOL isStart) ;
			static PbwtCursor *pbwtNakedCursorCreate (int M, int *aInit) ;
			static void pbwtCursorDestroy (PbwtCursor *u) ;
			static void pbwtCursorForwardsA (PbwtCursor *u) ; /* algorithm 1 in the manuscript */
			static void pbwtCursorForwardsAPacked (PbwtCursor *u) ; /* faster version, when have read y and set u->nBlockStart */
			static void pbwtCursorBackwardsA (PbwtCursor *u) ; /* undo algorithm 1 */
			static void pbwtCursorForwardsAD (PbwtCursor *u, int k) ; /* algorithm 2 in the manuscript */
			static void pbwtCursorCalculateU (PbwtCursor *x) ;   /* calculate u required for CursorMap */
			static void pbwtCursorForwardsRead (PbwtCursor *u) ; /* move forwards and read (unless at end) */
			static void pbwtCursorForwardsReadAD (PbwtCursor *u, int k) ;
			static void pbwtCursorReadBackwards (PbwtCursor *u) ; /* read and move backwards (unless at start) */
			static void pbwtCursorWriteForwards (PbwtCursor *u) ; /* write then move forwards */
			static void pbwtCursorWriteForwardsAD (PbwtCursor *u, int k) ;
			static void pbwtCursorToAFend (PbwtCursor *u, PBWT *p) ; /* utility to copy final u->a to p->aFend */
			/* basic update operations - inline them to make them tight */
			/* NB run pbwtCursorCalculateU() before pbwtCursorMap() */
			static inline int pbwtCursorMap (PbwtCursor *u, int x, int i)
			{ return x ? u->c + i - u->u[i] : u->u[i] ; }
			static inline int pbwtCursorMapDplus (PbwtCursor *u, int x, int i, int dplus)
			{ for ( ; i < u->M && u->y[i] != x ; ++i) if (u->d[i] > dplus) dplus = u->d[i] ;
				return dplus ;
			}
			static inline int pbwtCursorMapDminus (PbwtCursor *u, int x, int i, int dminus)
			{ for (--i ; i >= 0 && u->y[i] != x ; --i) if (u->d[i] > dminus) dminus = u->d[i] ;
				return dminus ;
			}

		/* low level operations on packed PBWT, argument yzp in these calls */

#define Y_SENTINEL 2			   /* needed to pack efficiently */
		static int pack3 (uchar *yp, int M, uchar *yzp) ; /* pack M values from yp into yzp */
		static int pack3arrayAdd (uchar *yp, int M, Array ayz) ; /* normally use this one */
		static int unpack3 (uchar *yzp, int M, uchar *yp, int *n0) ; /* unpack M values from yzp into yp, return number of bytes used from yzp, if (n0) write number of 0s into *n0 */
		static int packCountReverse (uchar *yzp, int M) ; /* return number of bytes to reverse one position */
		static int extendMatchForwards (uchar *yzp, int M, uchar x, int *f, int *g) ; /* move hit interval f,g) forwards one position, matching x */
		static int extendPackedForwards (uchar *yzp, int M, int *f, uchar *zp) ; /* move f forwards one position */
		static int extendPackedBackwards (uchar *yzp, int M, int *f, int c, uchar *zp) ; /* move f backwards one position - write value into *zp if zp non-zero */

		/* pbwtSample.c */

		static void sampleInit (void) ;
		static void sampleDestroy (void) ;
		static Sample *sample (PBWT *p, int i) ; /* give back Sample structure for sample i from p */
		static int  sampleAdd (char* name, char *father, char *mother, char *pop) ;
		static char* sampleName (Sample *s) ;
		static char* popName (Sample *s) ;	/* give back population name for sample i */
		static PBWT *pbwtSubSample (PBWT *pOld, Array select) ;
		static PBWT *pbwtSubSampleInterval (PBWT *pOld, int start, int Mnew) ;
		static PBWT *pbwtSelectSamples (PBWT *pOld, FILE *fp) ;

		/* pbwtIO.c */

		//extern int nCheckPoint ;	/* if set non-zero write pbwt and sites files every n sites when parsing external files */
		static int nCheckPoint;

		static void pbwtWrite (PBWT *p, FILE *fp) ; /* just writes packed PBWT p->yz */
		static void pbwtWriteSites (PBWT *p, FILE *fp) ;
		static void pbwtWriteSamples (PBWT *p, FILE *fp) ;
		static void pbwtWriteMissing (PBWT *p, FILE *fp) ;
		static void pbwtWriteDosage (PBWT *p, FILE *fp) ;
		static void pbwtWriteReverse (PBWT *p, FILE *fp) ;
		static void pbwtWriteAll (PBWT *p, char *fileNameRoot) ;
		static void pbwtWriteGen (PBWT *p, FILE *fp) ; /* write gen file as for impute etc. */
		static void pbwtWritePhase (PBWT *p, char *filename); /* Write phase file as output by impute and input for chromopainter */
		static PBWT *pbwtRead (FILE *fp) ;
		static Array pbwtReadSitesFile (FILE *fp, char **chrom) ;
		static void pbwtReadSites (PBWT *p, FILE *fp) ;
		static void pbwtReadRefFreq (PBWT *p, FILE *fp) ;
		static Array pbwtReadSamplesFile (FILE *fp) ;
		static void pbwtReadSamples (PBWT *p, FILE *fp) ;
		static void pbwtReadMissing (PBWT *p, FILE *fp) ;
		static void pbwtReadDosage (PBWT *p, FILE *fp) ;
		static void pbwtReadReverse (PBWT *p, FILE *fp) ;
		static PBWT *pbwtReadAll (char *fileNameRoot) ; /* reads .pbwt, .sites, .samples, .missing  */
		static PBWT *pbwtReadMacs (FILE *fp) ;
		static PBWT *pbwtReadVcfq (FILE *fp) ;	/* reduced VCF style file made by vcf query */
		static PBWT *pbwtReadGen (FILE *fp, char *chrom) ;	/* gen file as used by impute2 (unphased) */
		static PBWT *pbwtReadHap (FILE *fp, char *chrom) ;	/* hap file as used by impute2 (unphased) */
		static PBWT *pbwtReadHapLegend (FILE *fp, FILE *lp, char *chrom) ; /* hap and legend file as used by impute2 (phased) */
		PBWT *pbwtReadPhase (FILE *fp) ; /* Li and Stephens PHASE file */
		static void pbwtWriteHaplotypes (FILE *fp, PBWT *p) ;
		static void pbwtWriteTransposedHaplotypes (PBWT *p, FILE *fp) ;
		static void pbwtWriteImputeRef (PBWT *p, char *fileNameRoot) ;
		static void pbwtWriteImputeHapsG (PBWT *p, FILE *fp) ;
		static void pbwtCheckPoint (PbwtCursor *u, PBWT *p) ; /* need cursor to write end index */

		/* pbwtHtslib.c */
		/* all these functions also read and write samples and sites */
		PBWT *pbwtReadVcfGT (char *filename) ;	/* read GTs from vcf/bcf using htslib */
		PBWT *pbwtReadVcfPL (char *filename) ;	/* read PLs from vcf/bcf using htslib */
		// mode: wb=compressed BCF; wbu=uncompressed BCF; wz=compressed VCF; w=uncompressed VCF
		static void pbwtWriteVcf (PBWT *p, char *filename, char *reference_fname, char *mode) ;  /* write vcf/bcf using htslib */

		/* pbwtMatch.c - functions as in Bioinformatics 2014 paper */

		static void matchMaximalWithin (PBWT *p, void (*report)(int, int, int, int)) ;
		static void pbwtLongMatches (PBWT *p, int L) ; /* internal matches longer than L, maximal if L=0 */
		static void matchSequencesNaive (PBWT *p, FILE *fp) ; /* fp is a pbwt file of sequences to match */
		static void matchSequencesIndexed (PBWT *p, FILE *fp) ;
		static void matchSequencesDynamic (PBWT *p, FILE *fp) ;
		static void matchSequencesSweep (PBWT *p, PBWT *q, void (*report)(int, int, int, int)) ;
		static void matchSequencesSweepSparse (PBWT *p, PBWT *q, int nSparse,
				void (*report)(int, int, int, int, BOOL)) ;

		/* pbwtImpute.c */

		void imputeExplore (PBWT *p, int test) ;
		PBWT *phase (PBWT *p, int nSparse) ;
		PBWT *referencePhase (PBWT *p, char *fileNameRoot) ;
		PBWT *referenceImpute (PBWT *p, char *fileNameRoot, int nSparse, double fSparse) ;
		void genotypeCompare (PBWT *p, char *fileNameRoot) ;
		PBWT *imputeMissing (PBWT *p) ;
		PBWT *pbwtCorruptSites (PBWT *pOld, double pSite, double pChange) ;
		PBWT *pbwtCorruptSamples (PBWT *pOld, double pSample, double pChange) ;
		PBWT *pbwtCopySamples (PBWT *pOld, int Mnew, double meanLength) ;
		static void pbwtDosageStore (PBWT *p, double *dosage, int k) ;
		static double *pbwtDosageRetrieve (PBWT *p, PbwtCursor *u, double *dosage, int k) ; 
		/* if arg dosage == 0 then create and return, else fill and return; uses u->y */

		/* pbwtLikelihood.c */

		void pbwtFitAlphaBeta (PBWT *p, int model) ;
		void pbwtLogLikelihoodCopyModel (PBWT *p, double theta, double rho) ;

		/* pbwtPaint.c */

		void paintAncestryMatrix (PBWT *p, char *fileRoot,int chunksperregion) ;
		void paintAncestryMatrixSparse (PBWT *p, char *fileRoot,int chunksperregion,int cutoff) ;

		/* pbwtMerge.c */

		PBWT *pbwtMerge(const char **file_names, int nfiles);

		/* pbwtGeneticMap.c */

		void readGeneticMap (FILE *fp) ;
		double geneticMap (int x) ;
		void pbwt4hapsStats (PBWT *p) ;

};

/******************* end of file *******************/
