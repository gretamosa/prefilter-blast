#ifndef BITLIB_H_
#define BITLIB_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

//BUENAS
#define NUM_NUCLEOTIDE_POSITIONS 4194304
#define NUCL_BIT_SIZE 524288

#define BYTE_SIZE 8

unsigned char valueToBit(unsigned char * values);

unsigned char maskAnd(unsigned char value, int maskPos);

unsigned char * transformByteToBit(unsigned char * byteValues, unsigned int * lmerIndex, int numLmer);

unsigned char * transformBitToByte(unsigned char * bitValues);

unsigned char fromGuidToResult(unsigned char * bitValues, long int guid);

#endif
