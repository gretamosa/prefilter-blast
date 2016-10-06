#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>

#include "parse.h"
#include "prefilterCuda.h"

char * gb_dbFileName = NULL;
char * gb_blastInstance = NULL;
unsigned short int gb_blastWordSize = 3;
char * gb_queryFileName = NULL;
char * gb_blastSubstitutionMatrix = NULL;
unsigned short int gb_blastThreshold = 11;
float threshold = 0.5;
char * gb_out = NULL;
unsigned short int filter = 1;
int gpuCard = 0;

unsigned short int parseInputParameters(int argc, char ** argv)
{
	int c = 0;

	while ((c = getopt(argc, argv, "d:p:w:i:m:f:t:o:b:c:")) != -1)
	{
		switch (c)
		{
			case 'd':
				gb_dbFileName = optarg;
				break;
			case 'p':
				gb_blastInstance = optarg;
				if (memcmp(gb_blastInstance, "blastp", 6) == 0)
				{
					gb_blastWordSize = 3;
				}
				else if (memcmp(gb_blastInstance, "blastn", 6) == 0)
				{
					gb_blastWordSize = 11;
				}
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
			case 't':
				threshold = atof(optarg);
				break;
			case 'o':
				gb_out = optarg;
				break;
			case 'b':
				filter = atoi(optarg);
				break;
			case 'c':
				gpuCard = atoi(optarg);
				break;
			default:
				perror("Unknown option\n");
				fprintf(stderr,"-d -> database file name\n-i object file name\n-p blast option {blastp, blastn}\n-w blast word size\n-m blast substitution matrix\n-f blast threshold\n-t filter threshold\n-o output fasta file\n-b filter to be used {1- different aligment, 2- best alignment}\n-c gpu card to be used {0,1}\n");
				break;
		}
	}

	if (!gb_dbFileName || !gb_blastInstance || !gb_queryFileName)
	{
		fprintf(stderr,"-d -> database file name\n-i object file name\n-p blast option {blastp, blastn}\n-w blast word size\n-m blast substitution matrix\n-f blast threshold\n-t filter threshold\n-o output fasta file\n-b filter to be used {1- different aligment, 2- best alignment}\n-c gpu card to be used {0,1}\n");
		return 1;
	}

	if (!gb_blastSubstitutionMatrix)
	{
		gb_blastSubstitutionMatrix = (char *) calloc(9, sizeof(char));
		memcpy(gb_blastSubstitutionMatrix,"BLOSUM62", 8 * sizeof(char));
	}

	return 0;
}

void sys_call(char * sysdata)
{
	char * systemCall = NULL;
	int size = 0;
	
	if (sysdata)
	{
		size = strlen(sysdata);
		
		systemCall = (char *) calloc (size+1, sizeof(char));
		if (!systemCall)
		{
			perror("Unable to allocate memory\n");
			exit(-1);
		}
		memcpy(systemCall, sysdata, size*sizeof(char));
		
		system((const char *) systemCall);
		
		if (systemCall)
		{
			free((void *) systemCall);
		}
	}
}

void renameFile(char * oldFileName, char * newFileName)
{
	if (oldFileName && newFileName)
	{
		char aux[100];
#if DEBUG == 1
		fprintf(stderr, "pasando fichero de %s a %s\n", oldFileName, newFileName);
#endif
		
		sprintf(aux ,"mv %s %s%c", oldFileName,newFileName,'\0');
		sys_call(aux);
	}
	else
	{
		perror("Input/Output filename not specified\n");
		exit(-1);
	}
}

void prefiltering(dictionary * qdict, char * qsec, dictionary * dbdict, char * dbsec, unsigned int numSecDB, unsigned int numSecQ, unsigned int * tamSecDb, unsigned int * tamSecQ, unsigned long int sizeSecDB, unsigned long int sizeSecQ, unsigned short int lmerLength, float threshold, char * out, unsigned short int filter)
{
	int i,j;
	struct timeval iniSerie, finSerie;
	unsigned int ** qblock;
	unsigned short int * indxQBlock;
	unsigned int numBloques = 0;
	char * included;
	unsigned int * resultBlock; 
	unsigned int * resultBlockId;
	FILE * f = NULL;
	char aux [100];

	gettimeofday(&iniSerie, NULL);
	
	prefilterCuda(qsec, dbsec, numSecDB, numSecQ, tamSecDb, tamSecQ, sizeSecDB, sizeSecQ, lmerLength, 1, out, &qblock, &indxQBlock, &numBloques, &resultBlock, &resultBlockId, gpuCard); 
	
	included = (char *) calloc (numSecQ, sizeof(char));
		
	//Filtrar los resultados correctos.	
#if DEBUG == 1	
	fprintf(stderr,"threshold: %.2f nbloques: %d filtro: %d\n", threshold, numBloques, filter);
#endif
		
	if((out[0] == 'd') && (out[1] == 'b') && (out[2] == '/') && (out[3] == 'p') && (out[4] == 'r') && (out[5] == 'e') && (out[6] == '-'))
		sprintf(aux, "db/tmp-%s%c", out+7, '\0');
	else
		sprintf(aux, "db/tmp-%s%c", out, '\0');
	
	f = fopen(aux,"w+");	
	for (i = 0; i < numBloques; i++)
	{
#if DEBUG == 1
		printf("bloque: %d (", i, resultBlockId[i], resultBlock[i]);
#endif		
	
		for (j = 0; j < indxQBlock[i]; j++)
		{
			if (threshold > 1)
			{
#if DEBUG == 1
				printf("%d [ %.2f ] caso scores",qblock[i][j], resultBlock[i] * 1.0);				
#endif
			}
			else 
			{
				if(strlen(qdict[qblock[i][j]].seq) > TH_PER_BLOCK)
				{
#if DEBUG == 1
					printf("%d [ %.2f ] caso 1",qblock[i][j], (resultBlock[i] * 1.0) / TH_PER_BLOCK);				
#endif
				}
				else
				{
#if DEBUG == 1
					printf("%d [ %.2f ] caso 2",qblock[i][j], (resultBlock[i] * 1.0) / (strlen(qdict[qblock[i][j]].seq) - (lmerLength -1)));
#endif
				}
			}
			
#if DEBUG == 1
			printf(") dbId: %d score %d (256 caso 1)(%d caso 2) %s (sizeDB = %d)\n", resultBlockId[i], resultBlock[i], strlen(qdict[qblock[i][j]].seq), qdict[qblock[i][j]].gid, tamSecDb[resultBlockId[i]]);
#endif
		}
		
		//Si el score está por debajo del threshold se añade al nuevo fichero query.
		float promedio = 0.0;
		
		if(filter == 1)
		{
			for(j = 0; j < indxQBlock[i]; j++)
			{
				if(strlen(qdict[qblock[i][j]].seq) > TH_PER_BLOCK)
				{
					promedio = ((resultBlock[i] * 1.0) / TH_PER_BLOCK);
				}
				else
				{
					promedio = (resultBlock[i] * 1.0) / (strlen(qdict[qblock[i][j]].seq) - (lmerLength -1));
				}
				
				if (promedio <= threshold)
				{
					if(included[qblock[i][j]] == 0)
					{
						parseOut(f, qdict, qblock[i][j]);
						included[qblock[i][j]] = 1;
					}
				}
				//Una parte de la secuencia se parece a otra, el promedio es mayor que el threshold (luego no se mete chequean el resto de sus partes)
				else
				{
					included[qblock[i][j]] = 1;
				}
			}
		}		
		//Si el score está por encima del threshold se añade al nuevo fichero query.
		else if(filter == 2)
		{
			for(j = 0; j < indxQBlock[i]; j++)
			{
				if(strlen(qdict[qblock[i][j]].seq) > TH_PER_BLOCK)
				{
					promedio = ((resultBlock[i] * 1.0) / TH_PER_BLOCK);
				}
				else
				{
					promedio = (resultBlock[i] * 1.0) / (strlen(qdict[qblock[i][j]].seq) - (lmerLength -1));
				}
				
				if (promedio >= threshold)
				{
					if(included[qblock[i][j]] == 0)
					{
						parseOut(f, qdict, qblock[i][j]);
						included[qblock[i][j]] = 1;
					}
				}				
			}
		}
		else
		{
			fprintf(stderr, "filtro no reconocido = %d\n", filter);
		}
	}
	
	free (included);
	if(resultBlock)
		free(resultBlock);
	if(resultBlockId)
		free(resultBlockId);	
	if(indxQBlock)
		free(indxQBlock);	
	for(i = 0; i < numBloques; i++)
		free(qblock[i]);
	if(qblock)
		free(qblock);
	
	gettimeofday(&finSerie, NULL);

#if DEBUG == 1
	fprintf(stderr,"Time CUDA prefiltering: %ld msec ---- ---- \n", (((finSerie.tv_sec*1000000)+finSerie.tv_usec)-((iniSerie.tv_sec*1000000)+iniSerie.tv_usec))/1000);
#endif

	fclose(f);

	renameFile(aux, out);
	
	return;
}

int main(int argc, char* argv[])
{
	unsigned long int ql = 0;
	unsigned int qn = 0;
	unsigned int * qla = NULL;
	unsigned long int dl = 0;
	long int dn = 0;
	unsigned int * dla = NULL;
	dictionary * qdict = NULL;
	dictionary * dbdict = NULL;
	char * qsec = NULL;
	char * dbsec = NULL;

	if (parseInputParameters(argc, argv) == 1)
	{
		return 0;
	}

	//meter en el parse las secuencias ducplicadas?
	parseFDB(gb_dbFileName, &dbsec, &dl, &dn, &dla, gb_blastInstance);
	parse(gb_queryFileName, &qdict, &qsec, &ql, &qn, &qla, gb_blastInstance, gb_blastWordSize, gb_blastThreshold, gb_blastSubstitutionMatrix, 1);
	
	prefiltering (qdict, qsec, dbdict, dbsec, dn, qn, dla, qla, dl, ql, gb_blastWordSize, threshold, gb_out, filter);

	if (qdict)
	{
		free((void *)qdict);
	}
	if (qla)
	{
		free((void *)qla);
	}
	
	if (dbdict)
	{
		free((void *)dbdict);
	}
	if (dla)
	{
		free((void *)dla);
	}

	return 0;
}
