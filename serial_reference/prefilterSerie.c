#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>

#include "util/parse.h"

#define TH_PER_BLOCK 512

char * gb_dbFileName = NULL;
char * gb_blastInstance = NULL;
unsigned short int gb_blastWordSize = 3;
char * gb_queryFileName = NULL;
char * gb_blastSubstitutionMatrix = NULL;
unsigned short int gb_blastThreshold = 11;

unsigned short int parseInputParameters(int argc, char ** argv)
{
	int c = 0;

	while ((c = getopt(argc, argv, "d:p:w:i:m:f:")) != -1)
	{
		switch (c)
		{
			case 'd':
				gb_dbFileName = optarg;
				break;
			case 'p':
				gb_blastInstance = optarg;
				break;
			case 'w':
				gb_blastWordSize = atoi(optarg);
				break;
			case 'f':
				gb_blastThreshold = atoi(optarg);
				break;
			case 'i':
				gb_queryFileName = optarg;
				break;
			case 'm':
				gb_blastSubstitutionMatrix = optarg;
				break;
			default:
				perror("Unknown option\n");
				break;
		}
	}

	if (!gb_dbFileName || !gb_blastInstance || !gb_queryFileName)
	{
		return 1;
	}

	if (!gb_blastSubstitutionMatrix)
	{
		gb_blastSubstitutionMatrix = (char *) calloc(9, sizeof(char));
		memcpy(gb_blastSubstitutionMatrix,"BLOSUM62", 8 * sizeof(char));
	}

	return 0;
}

void prefilterSerie (char * secQ, char * secDB, unsigned int * tamanosQ, unsigned int * tamanosDB, int numQ, int numDB, int sizeQ, int sizeDB, short int lmerLength)
{
	int i, k, m, n;

	unsigned long int scsAct;
	unsigned long int scsMax = 0;
	unsigned int dbIdMax = 0;
	unsigned long int sumaPosDB = 0;
	int indx;
	char valido = 0;
	int numBloques;
	
	unsigned int scsAux = 0, scsAux2 = 0;

	if ((sizeQ % TH_PER_BLOCK) == 0)
		numBloques = sizeQ / TH_PER_BLOCK;
	else
		numBloques = (sizeQ / TH_PER_BLOCK) + 1;
		
	fprintf(stderr,"numBloques %d sizeQ %d sizeDB %d\n", numBloques,sizeQ, sizeDB);	

	//SE RECORREN LOS BLOQUES DE TH_PER_BLOCK LMERS
	for(i = 0; i < numBloques; i++)
	{
		sumaPosDB = 0;
		dbIdMax = 0;
		scsMax = 0;
		//Se recorre cada secuencia DB
		for(k = 0; k < numDB; k++)
		{
			indx = 0;
			scsAct = 0;

			//SE RECORREN LOS LMERS DEL BLOQUE
			while ((indx < TH_PER_BLOCK) && ((i * TH_PER_BLOCK) + indx < (sizeQ - (lmerLength - 1))))
			{
				//Se recorren los Lmer de la sec DB
				for(m = 0; m < (tamanosDB[k] - (lmerLength - 1)); m++)
				{
					valido = 1;
					for (n = 0; n < lmerLength; n++)
					{
						if(secQ[(i * TH_PER_BLOCK) + indx + n] != secDB[sumaPosDB + m + n])
						{
							valido = 0;
						}
					}

					if(valido == 1)
					{
						scsAct++;
					}
				}
				indx++;
			}

			//printf("db %d bloque : %d dbId: %d score : %ld sumapos:%d\n", k, i, dbIdMax, scsAct);
				
			scsAux2 = scsAct * 1000;
			scsAux = (scsAux2 / tamanosDB[k]);
			
			//Se guarda la mejor secuencia DB de contra el bloque
			if(scsAux >= scsMax)
			{				
				scsMax = scsAux;
				dbIdMax = k;
			}

			sumaPosDB += tamanosDB[k];
		}
		
		printf("bloque : %d dbId: %d score : %ld\n", i, dbIdMax, scsMax);
	}

	return;
}

void prefiltering(char * secuenciasQ, char * secuenciasDB, unsigned int * tamanosQ, unsigned int * tamanosDB, int numQ, int numDB, int sizeQ, int sizeDB, short int lmerLength)
{
	struct timeval iniSerie, finSerie;

	gettimeofday(&iniSerie, NULL);
	
	prefilterSerie(secuenciasQ, secuenciasDB, tamanosQ, tamanosDB, numQ, numDB, sizeQ, sizeDB, lmerLength);

	gettimeofday(&finSerie, NULL);

	fprintf(stderr,"Time Serie: %ld msec ---- ---- \n", (((finSerie.tv_sec*1000000)+finSerie.tv_usec)-((iniSerie.tv_sec*1000000)+iniSerie.tv_usec))/1000);

	return;
}

void createSecs(char * secuenciasQ, char * secuenciasDB, int size, int numSecDB, int numSecQ, unsigned short int * tamSecDb)
{
	int i, j;

	//QUERY
	for (i = 0; i < numSecQ; i++)
	{
		for(j = 0; j < size; j++)
		{
			secuenciasQ[(i * size) + j] = 'A';
		}
	}

	//DB
	for(i = 0; i < numSecDB; i++)
	{
		for(j = 0; j < size; j++)
		{
			secuenciasDB[(i * size) + j] = 'A';
		}

		tamSecDb[i] = size;
	}

	return;
}

int main(int argc, char* argv[])
{
	char * qseq = NULL;
	unsigned long int ql = 0;
	unsigned int qn = 0;
	unsigned int * qla = NULL;

	char * dbseq = NULL;
	unsigned long int dl = 0;
	unsigned int dn = 0;
	unsigned int * dla = NULL;


	parseInputParameters(argc, argv);

	parse(gb_dbFileName, &dbseq, &dl, &dn, &dla, gb_blastInstance, gb_blastWordSize);
	parse(gb_queryFileName, &qseq, &ql, &qn, &qla, gb_blastInstance, gb_blastWordSize);

	//fprintf(stderr,"DB (%ld): %s\n", dl, dbseq);
	//fprintf(stderr,"Q (%ld): %s\n", ql, qseq);
		
	prefiltering(qseq, dbseq, qla, dla, qn, dn, ql, dl, gb_blastWordSize);

	if (qseq)
	{
		free((void *)qseq);
	}
	if (qla)
	{
		free((void *)qla);
	}
	if (dbseq)
	{
		free((void *)dbseq);
	}
	if (dla)
	{
		free((void *)dla);
	}

	if(gb_blastSubstitutionMatrix)
		free(gb_blastSubstitutionMatrix);

	return 0;
}
