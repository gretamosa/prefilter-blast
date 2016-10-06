#ifndef PARSE_H_
#define PARSE_H_

#define PROT_ALP_SIZE 26

#include "prefilterCuda.h"

typedef struct filter_dictionary {
	char * gid;
	char * seq;
} dictionary;

static const unsigned short int gb_protIndx[26] = {0/*A*/, 20/*B*/, 4/*C*/, 3/*D*/, 6/*E*/, 13/*F*/, 7/*G*/, 8/*H*/, 9/*I*/, 21/*J*/, 11/*K*/, 10/*L*/, 12/*M*/, 2/*N*/, 23/*O*/, 14/*P*/, 5/*Q*/, 1/*R*/, 15/*S*/, 16/*T*/, 23/*U*/, 19/*V*/, 17/*W*/, 23/*X*/, 18/*Y*/, 22/*Z*/};

void parse (char * fn, dictionary ** dict, char ** sec, unsigned long int * ofl, unsigned int * ofn, unsigned ** ofla, char * bi, unsigned short int ws, unsigned short int threshold, char * sm_name, unsigned short int fs);
void parseOut(FILE * fout, dictionary * db, int indx);
void parseSequence(FILE * fParam, const unsigned short int * indx, unsigned short int alpsize, dictionary ** dict, char ** data, unsigned long int * ofl, unsigned int * ofn, unsigned int ** ofla, unsigned short int ws, unsigned short int threshold, char * sm_name, unsigned short int fs);
void parseSequencePerWarp(FILE * fParam, const unsigned short int * indx, unsigned short int alpsize, dictionary ** dict, char ** data, unsigned long int * ofl, unsigned int * ofn, unsigned int ** ofla, unsigned short int ws, unsigned short int threshold, char * sm_name, unsigned short int fs);
char * readGid(FILE * fParam);
char * readSeq(FILE * fParam);

void parseFDB (char * fn, char ** sec, unsigned long int * fl, long int * fn, unsigned int ** flsec, char * bi);

#endif
