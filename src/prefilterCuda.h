#ifndef PREFILTERCUDA_H_
#define PREFILTERCUDA_H_

#define TH_PER_BLOCK 256
#define NUM_SEC_DB_PROCESS 1000000

#define MISMATCH_PENALTY 1

#define PROTDBSIZE 13824
// *** CHANGE LIMIT DEFINE IF HW CHANGE
#define LIMIT_PROT_SIZE 108555

#define NUCL_BIT_SIZE 524288
// *** CHANGE LIMIT DEFINE IF HW CHANGE
#define LIMIT_NUCL_SIZE 2860

#define BYTE_SIZE 8

#define DEBUG 0
#define LOCAL_TEST 0
#define REDUCTION 0

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
	#define printf(f, ...) ((void)(f, __VA_ARGS__),0)
#endif

extern "C" void prefilterCuda(char * secuenciasQ, char * secuenciasDB, int numSecDB, int numSecQ, unsigned int * tamSecDb, unsigned int * tamSecQ, unsigned long int sizeSecDB, unsigned long int sizeSecQ, 
								short int lmerLength, char flagInicializacion, char * out, 
								unsigned int *** qblock, unsigned short int ** indxQBlock, unsigned int * numBloques, 
								unsigned int ** resultBlock ,unsigned int ** resultBlockId, int gpuCard, int scoreFilter);


#endif
