#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>

#include <cutil_inline.h>
#include "prefilterCuda.h"

typedef struct {
	
	int numBlkSize;
	unsigned long int sizeSecQ;
	
	char * secuenciasQ_h;

	unsigned int * resultBlock_h;
	unsigned int * resultBlockId_h;
	
	char * secuenciasQ_d;

	unsigned int * resultBlock_d;
	unsigned int * resultBlockId_d;
	
	unsigned int * tamSecDb_d;
	char * secuenciasDB_d;
	
    //Stream for asynchronous command execution
    cudaStream_t stream;
    
} GPUdata;

static const int gpuIndexes[5] = {0,1,2,3,4};

static const int PROTSIZE = 24;
static const int NUCLSIZE = 4;
__constant__ int PROTINDX[26] = {0/*A*/, 20/*B*/, 4/*C*/, 3/*D*/, 6/*E*/, 13/*F*/, 7/*G*/, 8/*H*/, 9/*I*/, 21/*J*/, 11/*K*/, 10/*L*/, 12/*M*/, 2/*N*/, 23/*O*/, 14/*P*/, 5/*Q*/, 1/*R*/, 15/*S*/, 16/*T*/, 23/*U*/, 19/*V*/, 17/*W*/, 23/*X*/, 18/*Y*/, 22/*Z*/};
__constant__ int NUCLINDX[26] = {0/*A*/, -1/*B*/, 1/*C*/, -1/*D*/, -1/*E*/, -1/*F*/, 2/*G*/, -1/*H*/, -1/*I*/, -1/*J*/, -1/*K*/, -1/*L*/, -1/*M*/, -1/*N*/, -1/*O*/, -1/*P*/, -1/*Q*/, -1/*R*/, -1/*S*/, 3/*T*/, -1/*U*/, -1/*V*/, -1/*W*/, -1/*X*/, -1/*Y*/, -1/*Z*/};

/*-------PROTOTYPES---------*/

__global__ void  cudaPrefilter(char * secQ, char * secDB, unsigned short int * tamanosDB, unsigned int * result, unsigned int * resultDB, int numDB, int tamanoQ, short int lmerLength);
void checkCUDAError(const char *msg);

/*--------------------------*/

//KERNEL 1
#define USE_SMEM_ATOMICS 0

#ifdef USE_SMEM_ATOMICS
	#ifdef CUDA_NO_SM13_ATOMIC_INTRINSICS
		#error Compilation target does not support shared-memory atomics
	#elif CUDA_NO_SM20_ATOMIC_INTRINSICS
		#error Compilation target does not support shared-memory atomics
	#else
		inline __device__ unsigned long long int suma(unsigned long long int * addr, unsigned long long int data)
		{
			return atomicAdd(addr, data);
		}
	#endif	
#else
	#error Shared-memory atomics not considered
#endif

#ifdef USE_SMEM_ATOMICS
	#ifdef CUDA_NO_SM13_ATOMIC_INTRINSICS
		#error Compilation target does not support shared-memory atomics
	#elif CUDA_NO_SM20_ATOMIC_INTRINSICS
		#error Compilation target does not support shared-memory atomics
	#else
		inline __device__ unsigned int suma2(unsigned int * addr, unsigned int data)
		{
			return atomicAdd(addr, data);
		}
	#endif	
#else
	#error Shared-memory atomics not considered
#endif

#ifdef USE_SMEM_ATOMICS
	#ifdef CUDA_NO_SM13_ATOMIC_INTRINSICS
		#error Compilation target does not support shared-memory atomics
	#elif CUDA_NO_SM20_ATOMIC_INTRINSICS
		#error Compilation target does not support shared-memory atomics
	#else
		inline __device__ unsigned int minimo(unsigned int * addr, unsigned int data)
		{
			return atomicMin(addr, data);
		}
	#endif	
#else
	#error Shared-memory atomics not considered
#endif

/**
* \fn void BlkToQuery (int numBloques, int numSecQ, unsigned int * tamSecQ, unsigned int *** qblock, unsigned short int ** indxQBlock, short int lmerLength)
* \brief Crea un array en el que contiene las relaciones entre las distintas querys y los bloques que las contienen.
* \param[in] numBloques, Numero de bloques totales.
* \param[in] numSecQ, Numero de secuencias totales.
* \param[in] tamSecQ, Array que contiene todos los tamaños de las secuencias.
* \param[in] qblock, Matriz que contiene los indices query que están contenidos en cada bloque.
* \param[in] indxQBlock, Array que contiene el número de secuencias query contenidas en casa bloque.
* \param[in] lmerLength, tamaño del lmer.
* \return .
*/
void BlkToQuery (int numBloques, int numSecQ, unsigned int * tamSecQ, unsigned int *** qblock, unsigned short int ** indxQBlock, short int lmerLength)
{
	int i, j, k;
	int contadorSec = 0;
	int contadorPos = 0;
	int contadorBlk = 0;
	int contadorAux = 0;
	int qStart = 0;

	(* qblock) = (unsigned int **) malloc (numBloques * sizeof(unsigned int *));
	(* indxQBlock) = (unsigned short int *) malloc (numBloques * sizeof(unsigned short int));
	
	for(i = 0; i < numSecQ; i++)
	{
		contadorSec++;
		for(j = 0; j < tamSecQ[i]; j++)
		{
			contadorPos++;
			
			if (contadorPos == TH_PER_BLOCK)
			{
				(* qblock)[contadorBlk] = (unsigned int *) malloc (contadorSec * sizeof(unsigned int));
				//Se inicializan el numero de secuencias de ese bloque.
				(* indxQBlock)[contadorBlk] = contadorSec;
				
				//Se meten las secuencias en la matriz.
				for(k = 0; k < contadorSec; k++)
				{
					contadorAux = qStart + k;
					(* qblock)[contadorBlk][k] = contadorAux;
				}
				
				contadorBlk++;
				contadorPos = 0;
				//Si se acaba ya la secuencia se usa la siguiente.
				if(j+1 ==  tamSecQ[i])
				{	
					contadorSec = 0;
					qStart = i+1;
				}
				else
				{
					contadorSec = 1;
					qStart = i;
				}
			}
		}
	}
	//Al terminar las secuencias se mete el último bloque.
	if(contadorPos > 0)
	{
		(* qblock)[contadorBlk] = (unsigned int *) malloc (contadorSec * sizeof(unsigned int));
		//Se inicializan el numero de secuencias de ese bloque.
		(* indxQBlock)[contadorBlk] = contadorSec;
		
		//Se meten las secuencias en la matriz.
		for(k = 0; k < contadorSec; k++)
		{
			contadorAux = qStart + k;
			(* qblock)[contadorBlk][k] = contadorAux;
		}
	}

	return;	
}

/**
* \fn __device__  unsigned int getHashIndx(char charParam, const int * indxParam)
* \brief Función que calcula el indice de un elemento del lmer.
* \param[in] charParam, Elemento del Lmer.
* \param[in] indxParam, Array de indices.
* \return indice del Lmer.
*/
__device__  unsigned int getHashIndx(char charParam, const int * indxParam) 
{
	unsigned int indx = 0;

	if (charParam >= 65 && charParam <= 90) {
		charParam -= 65;
		indx = indxParam[(int)charParam];
	} else if (charParam == 42) {
		indx = indxParam[23];
	} else {
		indx = indxParam[23];
	}

	return indx;
}

/**
* \fn __device__ int juan_hash (char * name, int length)
* \brief Función que calcula el valor hash de Nucleótidos o Proteínas.
* \param[in] name, Secuencia que contiene el Lmer.
* \param[in] length, Tamaño del lmer
* \return el valor hash del lmer.
*/
__device__ int juan_hash (char * name, int length) {
	short int i=0;
	int h = 0, h_aux = 0;

	for (i=0;i<length;i++) 
	{
		if ((name[i] == '*') || ((name[i] == 'U') && (length == 11)))
		{
			return -1;
		}
		else 
		{
			if (length == 3) 
			{
				h += getHashIndx(name[i], PROTINDX) * pow((double)PROTSIZE,length - (i + 1));
			}
			
			else if (length == 11) 
			{
				h_aux = getHashIndx(name[i], NUCLINDX) * pow((double)NUCLSIZE, length - (i + 1));
				if(h_aux != -1)
				{
					h += h_aux;
				}
				else
				{
					return -1;
				}
			}
		}
	}

	return h;
}

/**
* \fn unsigned char maskAnd(unsigned char value, int maskPos)
* \brief Función que saca de un byte el valor de un bit determinado por la posición.
* \param[in] value, Valor del que se debe sacar el dato.
* \param[in] maskPos, Máscara binaria para obtener el valor correcto.
* \return el valor del bit (1 o 0).
*/
__device__ unsigned int maskAnd(unsigned char value, int maskPos)
{
	unsigned char mask = 0;
	unsigned char result = 0;
	
	mask = pow((double) 2.0, maskPos);
	
	result = value & mask;
	
	if(result > 0)
		return 1;
	else
		return 0;
}

/**
* \fn unsigned char fromGuidToResult(unsigned char * bitValues, long int guid)
* \brief Función que devuelve un valor de la cadena de bits a partir de un GUID.
* \param[in] bitValues, Cadena de bytes comprimidos.
* \param[in] guid, Identificador Hash del valor buscado.
* \return valor buscado.
*/
__device__ unsigned int fromGuidToResult(unsigned char * bitValues, long int guid)
{
	int block = guid / BYTE_SIZE;
	int position = guid % BYTE_SIZE;
	
	return maskAnd(bitValues[block], position);
}

/**
* \fn __device__ unsigned int reduction(unsigned int * sdata, unsigned int tid)
* \brief Realiza una suma de datos usando el paralelismo de CUDA.
* \param[in] sdata, Array con todos los valores a sumar.
* \param[in] tid, Identificador de cada hilo.
* \return resultado de la suma.
*/
__device__ unsigned int reduction(unsigned int * sdata, unsigned int tid)
{	
	sdata[tid] += sdata[(tid + 128) % 256]; __syncthreads();
	sdata[tid] += sdata[(tid + 64) % 256]; __syncthreads();
	sdata[tid] += sdata[(tid + 32) % 256]; __syncthreads();
   
	sdata[tid] += sdata[(tid + 16) % 256];
	sdata[tid] += sdata[(tid +  8) % 256];
	sdata[tid] += sdata[(tid +  4) % 256];
	sdata[tid] += sdata[(tid +  2) % 256];
	sdata[tid] += sdata[(tid +  1) % 256];
	   
	if(tid == 0) 
	{
		return sdata[0];
	}
	
	return 0;
}

/**
* \fn __global__ void  cudaPrefilterProt (char * secQ, char * secDB, unsigned int * tamanosDB, unsigned int * result, unsigned int * resultDB, int numDB, unsigned long int tamanoQ, short int lmerLength, int scoreFilter)
* \brief Kernel Cuda para procesar Proteínas.
* \param[in] secQ, Secuencia con todas las secuencias query.
* \param[in] secDB, Secuencia con todas las secuencias de base de datos.
* \param[in] tamanosDB, Array que contiene todos los tamaños de las secuencias de base de datos.
* \param[in] result, Array que contiene el mejor Score obtenido por cada bloque.
* \param[in] resultDB, Array que contiene el identificador de la secuencia de base de datos con la que cada bloque ha obtenido el mejor Score.
* \param[in] numDB, Número de secuencias de base de datos.
* \param[in] tamanoQ, Tamaño total de todas las secuencias query.
* \param[in] lmerLength, tamaño del lmer.
* \param[in] scoreFilter, flag que discrimina que valor de score se debe devolver.
* \return .
*/
__global__ void  cudaPrefilterProt(char * secQ, char * secDB, unsigned int * tamanosDB, unsigned int * result, unsigned int * resultDB, int numDB, unsigned long int sizeDB, unsigned long int tamanoQ, short int lmerLength, int scoreFilter)
{
	unsigned int i;
	unsigned int scsMax = 0;
	double scsMax2 = 0;
	double pTotal = 0;	
	unsigned int dbIdMax = 0;
	unsigned int tid = threadIdx.x;
	unsigned int X = ((blockIdx.x * blockDim.x) + tid);
	long int guid = 0;
	
#if REDUCTION == 1
	__shared__ unsigned int sdata[256];
#else
	unsigned int scsAct;
#endif	

	__shared__ unsigned int scsHilo;
	__shared__ unsigned int scsHits;
	__shared__ char scsHitsFlag;
	__shared__ unsigned int valid_range;
	
	char letrita[3];
	
#if REDUCTION == 1
	sdata[tid] = 0;
#else
	scsAct = 0;
#endif

	i = 0;
	
	if (X < (tamanoQ - (lmerLength - 1)))
	{
		
		letrita[0] = secQ[X];
		letrita[1] = secQ[X + 1];
		letrita[2] = secQ[X + 2];
		
		valid_range = 257;
		
		if (letrita[2] == '*') {
			minimo(&valid_range, tid);
		}
		
		__syncthreads();
		
		guid = juan_hash(letrita, lmerLength);
		
		//Se recorren todas las db
		for (i = 0; i < numDB; i++)
		{
			scsHilo = 0;
			scsHits = 0;
			
			if (guid > -1) {
				if (guid >= sizeDB) {
#if REDUCTION == 0
					scsAct = 0;
#else
					sdata[tid] = 0;
#endif
				} else {
#if REDUCTION == 0
					scsAct = (unsigned int) secDB[(i*PROTDBSIZE) + guid];
#else
					sdata[tid] = (unsigned int) secDB[(i*PROTDBSIZE) + guid];
#endif
				}
			} else {
#if REDUCTION == 0
				scsAct = 0;
#else
				sdata[tid] = 0;
#endif
			}
			
			if (scsHitsFlag == 0) {
#if REDUCTION == 0
				if (scsAct > 1) {
					scsAct = 0;
				}
#else
				if (sdata[tid] > 1) {
					sdata[tid] = 0;
				}
#endif
			}
			
#if REDUCTION == 0
			// Se realiza la suma de los valores de cada hilo
			if (scsAct > 1)
			{
				// ATOMIC ADD
				scsHitsFlag = 1;
				suma2(&scsHilo, scsAct);
				if (scoreFilter != 2) 
				{
					suma2(&scsHits, 1);
				}
			} else { 
				//ATOMIC ADD
				suma2(&scsHilo, scsAct);
			}
#else
			if (sdata[tid] > 1) {
				scsHitsFlag = 1;
				// REDUCTION
				scsHilo = reduction(sdata, tid);
				if (scoreFilter != 2) {
					suma2(&scsHits, 1);
				}
			} else {
				// REDUCTION
				scsHilo = reduction(sdata, tid);
			}
#endif

			__syncthreads();			
			
			if (tid == 0) 
			{
				pTotal = scsHilo * 1000; 
				
				//if (tamanosDB[i] >= TH_PER_BLOCK)
				//{
				//	pTotal = pTotal / (tamanosDB[i] / (TH_PER_BLOCK * 1.0));
				//}
				
				// *** Elegir del mejor resultado
				if (pTotal >= scsMax2)
				{
					scsMax2 = pTotal;
					dbIdMax = i;
					
					// Caso Hits
					if(scsHitsFlag == 0)
					{
						scsMax = scsHilo;
					}
					else
					{
						// Caso Scores Promedio
						if (scoreFilter == 0)
						{
							scsMax = (scsHilo / scsHits);
						}
						// Caso Raw Scores con missmatch
						else if (scoreFilter == 1) 
						{
							// Mayor número de missmatch que de match
							if (((valid_range - 1) - scsHits) > scsHilo) 
							{
								scsMax = 0;
							} 
							// Mayor número de match que de missmatch
							else 
							{
								scsMax = scsHilo - (((valid_range - 1) - scsHits) * MISMATCH_PENALTY);
							}
						}
						// Caso Raw Scores 
						else if (scoreFilter == 2) 
						{
							scsMax = scsHilo;
						}
					}
				}
			}
		
			__syncthreads();
		}
		
		if(tid == 0)
		{			
			result[blockIdx.x] = scsMax;
			resultDB[blockIdx.x] = dbIdMax;
		}
	}

	return;
}

/**
* \fn __global__ void  cudaPrefilterNucl(char * secQ, char * secDB, unsigned int * tamanosDB, unsigned int * result, unsigned int * resultDB, int numDB, unsigned long int tamanoQ, short int lmerLength)
* \brief Kernel Cuda para procesar Nucleótidos
* \param[in] secQ, Secuencia con todas las secuencias query.
* \param[in] secDB, Secuencia con todas las secuencias de base de datos.
* \param[in] tamanosDB, Array que contiene todos los tamaños de las secuencias de base de datos.
* \param[in] result, Array que contiene el mejor Score obtenido por cada bloque.
* \param[in] resultDB, Array que contiene el identificador de la secuencia de base de datos con la que cada bloque ha obtenido el mejor Score.
* \param[in] numDB, Número de secuencias de base de datos.
* \param[in] tamanoQ, Tamaño total de todas las secuencias query.
* \param[in] lmerLength, tamaño del lmer.
* \return .
*/
__global__ void  cudaPrefilterNucl(char * secQ, char * secDB, unsigned int * tamanosDB, unsigned int * result, unsigned int * resultDB, int numDB, unsigned long int sizeDB, unsigned long int tamanoQ, short int lmerLength)
{
	unsigned int i;
	unsigned int scsAct = 0;
	//unsigned int sumaPosDB = 0;
	unsigned int scsMax = 0;
	double scsMax2 = 0;
	double pTotal = 0;
	unsigned int scsAux = 0;	
	unsigned int dbIdMax = 0;
	unsigned short int tid = threadIdx.x;
	unsigned int X = ((blockIdx.x * blockDim.x) + tid);
	long int guid = 0;
	
	__shared__ unsigned int scsHilo;
	
	char letritaNucl[11];
				
	if (X < (tamanoQ - (lmerLength - 1))) {		
		letritaNucl[0] = secQ[X];
		letritaNucl[1] = secQ[X + 1];
		letritaNucl[2] = secQ[X + 2];
		letritaNucl[3] = secQ[X + 3];
		letritaNucl[4] = secQ[X + 4];
		letritaNucl[5] = secQ[X + 5];
		letritaNucl[6] = secQ[X + 6];
		letritaNucl[7] = secQ[X + 7];
		letritaNucl[8] = secQ[X + 8];
		letritaNucl[9] = secQ[X + 9];
		letritaNucl[10] = secQ[X + 10];
		
		//Se recorren todas las db
		for (i = 0; i < numDB; i++) {
			scsAct = 0;
			scsHilo = 0;
					
			guid = juan_hash (letritaNucl, lmerLength);
			if (guid != -1)
			{
				scsAct = fromGuidToResult((unsigned char *) &secDB[i*NUCL_BIT_SIZE], guid);
			}
		 
			//Se realiza la suma de los valores de cada hilo mediante atomicAdd.
			suma2(&scsHilo, scsAct);
			
			__syncthreads();
			
			if(tid == 0)
			{
				scsAux = scsHilo * 1000;

				if(scsHilo == TH_PER_BLOCK)
				{
					pTotal = scsAux;
				}
				else if(tamanosDB[i] >= TH_PER_BLOCK)
				{
					pTotal = scsAux / (tamanosDB[i] / (TH_PER_BLOCK * 1.0));
				}
				else
				{
					pTotal = scsAux;
				}
				
				if(pTotal >= scsMax2)
				{
					scsMax2 = pTotal;
					scsMax = scsHilo;
					dbIdMax = i;
				}
			}			
			__syncthreads();
		}
		
		if(tid == 0)
		{	
			result[blockIdx.x] = scsMax;
			resultDB[blockIdx.x] = dbIdMax;
		}
	}

	return;
}

/**
* \fn extern "C" void prefilterCuda(char * secuenciasQ, char * secuenciasDB, int numSecDB, int numSecQ, unsigned int * tamSecDb, unsigned int * tamSecQ, unsigned long int sizeSecDB, unsigned long int sizeSecQ, short int lmerLength, char flagInicializacion, char * out, unsigned int *** qblock, unsigned short int ** indxQBlock, unsigned int * numBloques, unsigned int ** resultBlock, unsigned int ** resultBlockId, int gpuCard, int scoreFilter)
* \brief Función que prepara, reserva e inicializa las estructuras y arrays necesarias para la ejecución de CUDA.
* \param[in] secuenciasQ, Secuencia con todas las secuencias query.
* \param[in] secuenciasDB, Secuencia con todas las secuencias de base de datos.
* \param[in] numSecDB, Número de secuencias de base de datos.
* \param[in] numSecQ, Número de secuencias query.
* \param[in] tamSecDb, Array que contiene todos los tamaños de las secuencias de base de datos.
* \param[in] tamSecQ, Array que contiene todos los tamaños de las secuencias de base de datos.
* \param[in] sizeSecDB, Tamaño total de todas las secuencias de la base de datos.
* \param[in] sizeSecQ, Tamaño total de todas las secuencias query.
* \param[in] lmerLength, Tamaño del Lmer.
* \param[in] flagInicializacion, Flag que indica si el sistema ya se ha inicializado (No se usa).
* \param[in] out, Nombre del fichero de salida (No se usa).
* \param[in] qblock, Matriz que contiene los indices query que están contenidos en cada bloque.
* \param[in] indxQBlock, Array que contiene el número de secuencias query contenidas en casa bloque.
* \param[in] numBloques, Número de bloques que se crearán en la ejecución del Kernel CUDA.
* \param[in] resultBlock, Array que contiene el mejor Score obtenido por cada bloque.
* \param[in] resultBlockId, Array que contiene el identificador de la secuencia de base de datos con la que cada bloque ha obtenido el mejor Score.
* \param[in] gpuCard, Número de Gpus a usar (No se usa).
* \param[in] scoreFilter, flag que discrimina que valor de score se debe devolver.
* \return .
*/
extern "C" void prefilterCuda(char * secuenciasQ, char * secuenciasDB, int numSecDB, int numSecQ, unsigned int * tamSecDb, unsigned int * tamSecQ, unsigned long int sizeSecDB, unsigned long int sizeSecQ, 
								short int lmerLength, char flagInicializacion, char * out, unsigned int *** qblock, unsigned short int ** indxQBlock, unsigned int * numBloques,
								unsigned int ** resultBlock, unsigned int ** resultBlockId, int gpuCard, int scoreFilter)
{
	struct timeval iniCuda, finCuda;
	struct timeval ini1, fin1;
	struct timeval ini2, fin2;
	struct timeval ini3, fin3;
	
	int device_count = 0;
	int device_index = 0;
	
	GPUdata paralelData[5];
	int acumulador = 0;
	int acumulador2 = 0;
			
	//Inicialización de la GPU
	cudaGetDeviceCount(&device_count);
	device_count = 4;
	
	if (numSecQ == 1)
		device_count = 1;

	fprintf(stderr,"Device Numbers: %d\n", device_count);
	
	//Cálculo del numero de bloques que se van a ejecutar.
	if ((sizeSecQ % TH_PER_BLOCK) == 0)
		(*numBloques) = sizeSecQ / TH_PER_BLOCK;
	else
		(*numBloques) = (sizeSecQ / TH_PER_BLOCK) + 1;

	fprintf(stderr,"numBloques %d sizeQ %ld sizeDB %ld\n", (*numBloques), sizeSecQ, sizeSecDB);	
	
	(*resultBlock) = (unsigned int *) malloc ((*numBloques) * sizeof(unsigned int));
	(*resultBlockId) = (unsigned int *) malloc ((*numBloques) * sizeof(unsigned int));
	
	acumulador = 0;
	
	//Inicialización del numero de bloques que va a ejecutar cada GPU
	//Inicialización de las posiciones de memoria de los arrays de entrada para cada GPU.
	//Creación de el hilo de envío de cada GPU.
	for (device_index = 0; device_index < device_count; device_index++)
	{
		cudaSetDevice(gpuIndexes[device_index]);
		
		paralelData[device_index].numBlkSize = (*numBloques) / device_count;
		
		if(acumulador + (((*numBloques) / device_count) * TH_PER_BLOCK) > sizeSecQ)
		{
			paralelData[device_index].sizeSecQ = sizeSecQ - acumulador;
		} else {
			paralelData[device_index].sizeSecQ = ((*numBloques) / device_count) * TH_PER_BLOCK;
		}
		
		paralelData[device_index].secuenciasQ_h = secuenciasQ + acumulador;

		paralelData[device_index].resultBlock_h = (*resultBlock) + acumulador2;
		paralelData[device_index].resultBlockId_h = (*resultBlockId) + acumulador2;

#if DEBUG == 1
		fprintf(stderr,"leido desde bloque %d hasta %d, leido desde sizeSecQ %d hasta %d GPU : %d\n", 
				acumulador2, ((*numBloques) / device_count) + acumulador2, acumulador, 
				paralelData[device_index].sizeSecQ + acumulador, device_index);
#endif
		
		acumulador += ((*numBloques) / device_count) * TH_PER_BLOCK;
		acumulador2 += (*numBloques) / device_count;
		
		cudaStreamCreate(&paralelData[device_index].stream);
	}
	
	//El último bloque lo procesa la última GPU (en caso de no ser una división esacta).
	int restoBloques = ((*numBloques) % device_count);
	if (restoBloques != 0)
	{
		paralelData[device_index-1].numBlkSize += restoBloques;
		paralelData[device_index-1].sizeSecQ += sizeSecQ - acumulador;

#if DEBUG == 1		
		fprintf(stderr, "anadidos %d bloques adicionales, GPU: %d desde: %d hasta: %d\n", restoBloques, device_index-1, acumulador, acumulador + (sizeSecQ - acumulador));
#endif
	}
				
	gettimeofday(&iniCuda, NULL);
	
	//Generar relación bloques query
	BlkToQuery ((*numBloques), numSecQ, tamSecQ, qblock, indxQBlock, lmerLength);
	
	if (numSecDB <= NUM_SEC_DB_PROCESS)
	{
		gettimeofday(&ini1, NULL);
		
		//Reserva de memoria de los datos GPU.	
		for (device_index = 0; device_index < device_count; device_index++)
		{
			cudaSetDevice(gpuIndexes[device_index]);
			
			// *** Reserva memoria CUDA -- BASE DE DATOS.
			cudaMalloc((void**) &paralelData[device_index].secuenciasDB_d, sizeSecDB * sizeof(char));
			cudaMalloc((void**) &paralelData[device_index].tamSecDb_d, numSecDB * sizeof(unsigned int));
			checkCUDAError("cuda malloc db");
			
			// *** Reserva de memoria CUDA -- QUERIES Y RESULTADOS
			cudaMalloc((void**) &paralelData[device_index].secuenciasQ_d, paralelData[device_index].sizeSecQ * sizeof(char));
			cudaMalloc((void**) &paralelData[device_index].resultBlock_d, paralelData[device_index].numBlkSize * sizeof(unsigned int));
			cudaMalloc((void**) &paralelData[device_index].resultBlockId_d, paralelData[device_index].numBlkSize * sizeof(unsigned int));
			checkCUDAError("cuda malloc q");
			
			// *** Copia de datos a CUDA	
			// NOTA: Se ha usado el mismo array para copiar la DB en todas las GPU.
			cudaMemcpyAsync(paralelData[device_index].secuenciasDB_d, secuenciasDB,  sizeSecDB * sizeof(char), cudaMemcpyHostToDevice, paralelData[device_index].stream);
			cudaMemcpyAsync(paralelData[device_index].tamSecDb_d, tamSecDb,  numSecDB * sizeof(unsigned int), cudaMemcpyHostToDevice, paralelData[device_index].stream);
			
				// SYNC MODE
				//cudaMemcpy(paralelData[device_index].secuenciasDB_d, secuenciasDB, sizeSecDB * sizeof(char), cudaMemcpyHostToDevice);
				//cudaMemcpy(paralelData[device_index].tamSecDb_d, tamSecDb, numSecDB * sizeof(unsigned int), cudaMemcpyHostToDevice);
			
			checkCUDAError("cuda memcpy Send BD");
			
			cudaMemcpyAsync(paralelData[device_index].secuenciasQ_d, paralelData[device_index].secuenciasQ_h, paralelData[device_index].sizeSecQ * sizeof(char), cudaMemcpyHostToDevice, paralelData[device_index].stream);
			 		 
				// SYNC MODE
				//cudaMemcpy(paralelData[device_index].secuenciasQ_d, paralelData[device_index].secuenciasQ_h, paralelData[device_index].sizeSecQ * sizeof(char), cudaMemcpyHostToDevice);
			
			checkCUDAError("cuda memcpy Send Q");
		}
			
		gettimeofday(&fin1, NULL);
#if DEBUG == 1
		fprintf(stderr,"Time CUDA send: %ld msec ---- ---- \n", (((fin1.tv_sec*1000000)+fin1.tv_usec)-((ini1.tv_sec*1000000)+ini1.tv_usec))/1000);
#endif
		
		gettimeofday(&ini2, NULL);
			
		//Inicialización de los hilos y bloques para CUDA, y lanzamiento del Kernel	
		for (device_index = 0; device_index < device_count; device_index++)
		{
			cudaSetDevice(gpuIndexes[device_index]);
			
			dim3 dimBlock(TH_PER_BLOCK, 1);
			dim3 dimGrid(paralelData[device_index].numBlkSize, 1);

			
			//Kernel
			if(lmerLength == 3)
			{
				cudaPrefilterProt <<< dimGrid, dimBlock, 0, paralelData[device_index].stream >>>
				(paralelData[device_index].secuenciasQ_d, paralelData[device_index].secuenciasDB_d, 
				paralelData[device_index].tamSecDb_d, paralelData[device_index].resultBlock_d, 
				paralelData[device_index].resultBlockId_d, numSecDB, sizeSecDB, paralelData[device_index].sizeSecQ, lmerLength, scoreFilter);
			}
			else if (lmerLength == 11)
			{
				cudaPrefilterNucl <<< dimGrid, dimBlock, 0, paralelData[device_index].stream >>>
				(paralelData[device_index].secuenciasQ_d, paralelData[device_index].secuenciasDB_d, 
				paralelData[device_index].tamSecDb_d, paralelData[device_index].resultBlock_d, 
				paralelData[device_index].resultBlockId_d, numSecDB, sizeSecDB, paralelData[device_index].sizeSecQ, lmerLength);
			}
		}
		
		for (device_index = 0; device_index < device_count; device_index++)
		{
			// *** Puede quitarse en modo release (únicamente afecta en la medida de tiempos no en el rendimiento)
			cudaThreadSynchronize();																		   
			checkCUDAError("kernel invocation");
		}
		
		gettimeofday(&fin2, NULL);
#if DEBUG == 1
		fprintf(stderr,"Time CUDA kernel: %ld msec ---- ---- \n", (((fin2.tv_sec*1000000)+fin2.tv_usec)-((ini2.tv_sec*1000000)+ini2.tv_usec))/1000);
#endif
			
		gettimeofday(&ini3, NULL);
		
		//Recuperación de los resultados.
		for (device_index = 0; device_index < device_count; device_index++)	
		{	
			cudaSetDevice(gpuIndexes[device_index]);
			
			cudaMemcpyAsync(paralelData[device_index].resultBlock_h, paralelData[device_index].resultBlock_d, paralelData[device_index].numBlkSize * sizeof(unsigned int), cudaMemcpyDeviceToHost, paralelData[device_index].stream);
			cudaMemcpyAsync(paralelData[device_index].resultBlockId_h, paralelData[device_index].resultBlockId_d, paralelData[device_index].numBlkSize * sizeof(unsigned int), cudaMemcpyDeviceToHost, paralelData[device_index].stream);
			
			//cudaMemcpy(paralelData[device_index].resultBlock_h, paralelData[device_index].resultBlock_d, paralelData[device_index].numBlkSize * sizeof(unsigned int), cudaMemcpyDeviceToHost);
			//cudaMemcpy(paralelData[device_index].resultBlockId_h, paralelData[device_index].resultBlockId_d, paralelData[device_index].numBlkSize * sizeof(unsigned int), cudaMemcpyDeviceToHost);
			
			checkCUDAError("cuda memcpy Recv");
		}

#if DEBUG == 1
		fprintf(stderr,"Synchronizing Streams\n");
#endif
		acumulador = 0;
		//Se juntan todas los resultados parciales.
		for (device_index = 0; device_index < device_count; device_index++)
		{
			cudaSetDevice(gpuIndexes[device_index]);
			checkCUDAError("cudaSetDevice");
			
			cudaStreamSynchronize(paralelData[device_index].stream);
			checkCUDAError("cudaStreamSynchronize");
			
			// *** Remove to avoid overlapping
			// memcpy((*resultBlock) + acumulador, paralelData[device_index].resultBlock_h, paralelData[device_index].numBlkSize * sizeof(unsigned int));
			// memcpy((*resultBlockId) + acumulador, paralelData[device_index].resultBlockId_h, paralelData[device_index].numBlkSize * sizeof(unsigned int));
			acumulador += paralelData[device_index].numBlkSize;
		}
		
		gettimeofday(&fin3, NULL);
#if DEBUG == 1		
		fprintf(stderr,"Time CUDA recv: %ld msec ---- ---- \n", (((fin3.tv_sec*1000000)+fin3.tv_usec)-((ini3.tv_sec*1000000)+ini3.tv_usec))/1000);
#endif
	}
	// *** El numero de secuencias BD es demasiado grande
	else
	{
		fprintf(stderr,"Demasiadas secuencias de Base de Datos\n");
		exit(0);
	}
		
	//Liberación memoria CUDA
	for (device_index = 0; device_index < device_count; device_index++)
	{
#if DEBUG == 1
		fprintf(stderr,"Liberando recursos del device %d\n",device_index);
#endif
		cudaSetDevice(gpuIndexes[device_index]);
		checkCUDAError("cudaSetDevice");
		
		if (paralelData[device_index].secuenciasQ_d)
			cudaFree(paralelData[device_index].secuenciasQ_d);
		if (paralelData[device_index].resultBlock_d)
			cudaFree(paralelData[device_index].resultBlock_d);
		if (paralelData[device_index].resultBlockId_d)
			cudaFree(paralelData[device_index].resultBlockId_d);
		
		if (paralelData[device_index].tamSecDb_d)
			cudaFree(paralelData[device_index].tamSecDb_d);
		if (paralelData[device_index].secuenciasDB_d)
			cudaFree(paralelData[device_index].secuenciasDB_d);
			
		checkCUDAError("free memory");
		cudaStreamDestroy(paralelData[device_index].stream);
		
		cudaDeviceReset();
	}
	
	gettimeofday(&finCuda, NULL);
	
#if DEBUG == 1
	fprintf(stderr,"Time CUDA: %ld msec ---- ---- \n", (((finCuda.tv_sec*1000000)+finCuda.tv_usec)-((iniCuda.tv_sec*1000000)+iniCuda.tv_usec))/1000);
#endif
	
	return;
}

/**
* \fn void sys_call(char * sysdata)
* \brief Función que realiza una llamada al sistema con los datos pasados.
* \param[in] sysdata, Datos para la llamada al sistema.
* \return .
*/
void report_error(char * sysdata)
{
	char systemCall[200];
	int size = 0;
	
	if (sysdata)
	{
		memset(systemCall, 0, 200 * sizeof(char));
		size = strlen(sysdata);
		memcpy(systemCall, sysdata, size*sizeof(char));
		system((const char *) systemCall);
	}
}

/**
* \fn void checkCUDAError(const char *msg)
* \brief Función que evalua si ha habido algún error de ejecucion en CUDA.
* \param[in] msg, Mensaje devuelto por la librería de CUDA.
* \return .
*/
void checkCUDAError(const char *msg) {
    cudaError_t err = cudaGetLastError();
    int flagError = 0, device_index = 0, device_count = 4;
    char aux[100];
    
    if (cudaErrorPriorLaunchFailure == err)
    {
    	memset(aux, 0, 100 * sizeof(char));
#if LOCAL_TEST == 0
		sprintf(aux , "blast_reporter \"Cuda error: cudaErrorPriorLaunchFailure: %s: %s.\" %c", msg, cudaGetErrorString(err),'\0');
#else
		sprintf(aux , "./blast_reporter \"Cuda error: cudaErrorPriorLaunchFailure: %s: %s.\" %c", msg, cudaGetErrorString(err),'\0');
#endif
		report_error(aux);
		fprintf(stderr, "Cuda error: cudaErrorPriorLaunchFailure: %s: %s.\n", msg, cudaGetErrorString(err) );
        flagError = 1;
	}
	if (cudaErrorLaunchTimeout == err)
    {
    	memset(aux, 0, 100 * sizeof(char));
#if LOCAL_TEST == 0
		sprintf(aux , "blast_reporter \"Cuda error: cudaErrorLaunchTimeout: %s: %s.\" %c", msg, cudaGetErrorString(err),'\0');
#else
		sprintf(aux , "./blast_reporter \"Cuda error: cudaErrorLaunchTimeout: %s: %s.\" %c", msg, cudaGetErrorString(err),'\0');
#endif
		report_error(aux);
		fprintf(stderr, "Cuda error: cudaErrorLaunchTimeout: %s: %s.\n", msg, cudaGetErrorString(err) );
        flagError = 1;
	}
	if (cudaErrorInvalidDeviceFunction == err)
    {
    	memset(aux, 0, 100 * sizeof(char));
#if LOCAL_TEST == 0
		sprintf(aux ,"blast_reporter \"Cuda error: cudaErrorInvalidDeviceFunction: %s: %s.\" %c", msg, cudaGetErrorString(err),'\0');
#else
		sprintf(aux ,"./blast_reporter \"Cuda error: cudaErrorInvalidDeviceFunction: %s: %s.\" %c", msg, cudaGetErrorString(err),'\0');
#endif	
		report_error(aux);
		fprintf(stderr, "Cuda error: cudaErrorInvalidDeviceFunction: %s: %s.\n", msg, cudaGetErrorString(err) );
        flagError = 1;
	}
	if (cudaErrorInvalidValue == err)
    {
    	memset(aux, 0, 100 * sizeof(char));
#if LOCAL_TEST == 0
		sprintf(aux ,"blast_reporter \"Cuda error: cudaErrorInvalidValue: %s: %s.\" %c", msg, cudaGetErrorString(err),'\0');
#else
		sprintf(aux ,"./blast_reporter \"Cuda error: cudaErrorInvalidValue: %s: %s.\" %c", msg, cudaGetErrorString(err),'\0');
#endif
		report_error(aux);
		fprintf(stderr, "Cuda error: cudaErrorInvalidValue: %s: %s.\n", msg, cudaGetErrorString(err) );
        flagError = 1;
	}
	
    if (cudaSuccess != err) {
    	memset(aux, 0, 100 * sizeof(char));
#if LOCAL_TEST == 0
		sprintf(aux ,"blast_reporter \"Cuda error: %s: %s.\" %c", msg, cudaGetErrorString(err),'\0');
#else
		sprintf(aux ,"./blast_reporter \"Cuda error: %s: %s.\" %c", msg, cudaGetErrorString(err),'\0');
#endif
		report_error(aux);
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString(err) );
        flagError = 1;
    }
    
    if (flagError == 1) {
    	//Liberación memoria CUDA
		for (device_index = 0; device_index < device_count; device_index++)
		{
#if DEBUG == 1
			fprintf(stderr,"Liberando recursos del device %d\n",device_index);
#endif
			cudaError_t ret = cudaSetDevice(gpuIndexes[device_index]);
			if (ret != cudaSuccess) {
				fprintf(stderr,"\t%s\n",cudaGetErrorString(ret));
			}
			
			ret = cudaDeviceReset();
			if (ret != cudaSuccess) {
				fprintf(stderr,"\t%s\n",cudaGetErrorString(ret));
			}
		}
		
		exit(-1);
    }
}
