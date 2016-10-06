#ifndef PARSE_H_
#define PARSE_H_

#define PROT_ALP_SIZE 26

static const unsigned short int gb_protIndx[26] = {0/*A*/, 20/*B*/, 4/*C*/, 3/*D*/, 6/*E*/, 13/*F*/, 7/*G*/, 8/*H*/, 9/*I*/, 21/*J*/, 11/*K*/, 10/*L*/, 12/*M*/, 2/*N*/, 23/*O*/, 14/*P*/, 5/*Q*/, 1/*R*/, 15/*S*/, 16/*T*/, 23/*U*/, 19/*V*/, 17/*W*/, 23/*X*/, 18/*Y*/, 22/*Z*/};

void parse (char * fn, char ** ofs, unsigned long int * ofl, unsigned int * ofn, unsigned ** ofla, char * bi, unsigned short int ws);
void parseSequence(FILE * fParam, const unsigned short int * indx, unsigned short int alpsize, char ** data, unsigned long int * ofl, unsigned int * ofn, unsigned int ** ofla, unsigned short int ws);
char * readGid(FILE * fParam);
char * readSeq(FILE * fParam);

#endif
