#ifndef UTILS_H_
#define UTILS_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "bitLib.h"

// ------------------ Global Variable -------------------------
#define HASHSIZE 2384185791015625

static const int PROTINDX[26] = {0/*A*/, 20/*B*/, 4/*C*/, 3/*D*/, 6/*E*/, 13/*F*/, 7/*G*/, 8/*H*/, 9/*I*/, 21/*J*/, 11/*K*/,10/*L*/, 12/*M*/, 2/*N*/, 23/*O*/, 14/*P*/, 5/*Q*/, 1/*R*/, 15/*S*/, 16/*T*/, 23/*U*/, 19/*V*/, 17/*W*/, 23/*X*/, 18/*Y*/, 22/*Z*/};
static const int PROTSIZE = 24;

static const int NUCLINDX[4] = { 0/*A*/, 1/*C*/, 2/*G*/, 3/*T*/};
static const int NUCLSIZE = 4;

// ------------------ Substitution Matrixes -------------------------
static const int Blosum62PSM[25 * 25] = {
    /*       A,  R,  N,  D,  C,  Q,  E,  G,  H,  I,  L,  K,  M,  F,  P,  S,  T,  W,  Y,  V,  B,  J,  Z,  X,  *        */
    /*A*/    4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1, -1, -1, -4,
    /*R*/   -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1, -2,  0, -1, -4,
    /*N*/   -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  4, -3,  0, -1, -4,
    /*D*/   -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4, -3,  1, -1, -4,
    /*C*/    0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -1, -3, -1, -4,
    /*Q*/   -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0, -2,  4, -1, -4,
    /*E*/   -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1, -3,  4, -1, -4,
    /*G*/    0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -4, -2, -1, -4,
    /*H*/   -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0, -3,  0, -1, -4,
    /*I*/   -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1, 0, -3, -2, -1, -3, -1,  3, -3,  3, -3, -1, -4,
    /*L*/   -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2, 0, -3, -2, -1, -2, -1,  1, -4,  3, -3, -1, -4,
    /*K*/   -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0, -3,  1, -1, -4,
    /*M*/   -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5, 0, -2, -1, -1, -1, -1,  1, -3,  2, -1, -1, -4,
    /*F*/   -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0, 6, -4, -2, -2,  1,  3, -1, -3,  0, -3, -1, -4,
    /*P*/   -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -3, -1, -1, -4,
    /*S*/    1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0, -2,  0, -1, -4,
    /*T*/    0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1, -1, -1, -4,
    /*W*/   -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11,  2, -3, -4, -2, -2, -1, -4,
    /*Y*/   -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1, 3, -3, -2, -2,  2,  7, -1, -3, -1, -2, -1, -4,
    /*V*/    0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3,  2, -2, -1, -4,
    /*B*/   -2, -1,  4,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4, -3,  0, -1, -4,
    /*J*/   -1, -2, -3, -3, -1, -2, -3, -4, -3,  3,  3, -3,  2, 0, -3, -2, -1, -2, -1,  2, -3,  3, -3, -1, -4,
    /*Z*/   -1,  0,  0,  1, -3,  4,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -2, -2, -2,  0, -3,  4, -1, -4,
    /*X*/   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -4,
    /***/   -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1
};

static const char Blosum62ALP[26] = {"ARNDCQEGHILKMFPSTWYVBJZX\0"};
static const int Blosum62INDX[26] = {0/*A*/, 20/*B*/, 4/*C*/, 3/*D*/, 6/*E*/, 13/*F*/, 7/*G*/, 8/*H*/, 9/*I*/, 21/*J*/, 11/*K*/,
			10/*L*/, 12/*M*/, 2/*N*/, 23/*O*/, 14/*P*/, 5/*Q*/, 1/*R*/, 15/*S*/, 16/*T*/, 23/*U*/, 19/*V*/, 17/*W*/, 23/*X*/, 18/*Y*/, 22/*Z*/};
static const int Blosum62ALP_SIZE = 25;

// ------------------ Global Variable -------------------------

// ------------------- Structures -----------------------------
typedef struct DICTIOSEQ {
	char * sequence;
	int l_sequence;
	char * gid;
	int l_gid;
	int strand;
	double lambda;
	double k;
	double h;
} DICTIOSEQ;

typedef struct DICTIONARY {
	struct DICTIOSEQ * sequences;
} DICTIONARY;
// ------------------- Structures -----------------------------

// ------------------ Methods declarations --------------------
int parseFasta (FILE * fParam, int protFlag, unsigned int * tamSec, unsigned char ** formatedSec, int threshold, char * psmName, int fdb_score);

long * fetchingDBwithStrand(char * dbFileName, int nucleotideStrand);

char * readGID_lmc(FILE * fParam);
char * readSeq_lmc(FILE * fParam);

int processLmersProtWithHits(char ** formatedSec , unsigned int * tamSec, char * seq, int lmerl, int qid, int threshold, char * psmName);
int processLmersProtWithScores(char ** formatedSec , unsigned int * tamSec, char * seq, int lmerl, int qid, int threshold, char * psmName);
int processLmersNucl(unsigned char ** formatedSec , unsigned int * tamSec, char * seq, int lmerl, int qid);
int parseSequence(FILE * fParam, char * gid, char * seq, int protFlag, int qid, unsigned int * tamSec, unsigned char ** formatedSec, int threshold, char * psmName, int fdb_score);
unsigned long juan_hash ( char * name, int length);
int getHashIndx(char charParam, const int * indxParam);

long * seedingWithScores(char * lmer, int size, const char * alp, const int alpsize, const int * psm, const int * indx, int threshold, int * neighbour_counter, int ** neighbour_score);
long * seeding(char * lmer, int size, const char * alp, const int alpsize, const int * psm, const int * indx, int threshold, int * neighbour_counter);
int getMatrixScore(const char * alpParam, const int alpsize, const int * psmParam, const int * indxParam, char firstChar, char secondChar);
int getNuclHashIndx(char charParam, const int * indxParam);
// ------------------ Methods declarations --------------------

#endif
