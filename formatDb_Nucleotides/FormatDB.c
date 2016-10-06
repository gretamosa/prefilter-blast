#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "util/utils.h"
#include "util/bitLib.h"

// Input file name to parse
char * gb_outputFileName = NULL;

// Input file name to parse
char * gb_inputFileName = NULL;

// Proteins or Nucleotides flag
int gb_proteinsFlag = 0;
int gb_formatDbScoresFlag = 0;
char gb_psmName[20] = "BLOSUM62";
int gb_threshold = 11;

void parseInputParameters(int argc, char ** argv) {
	int c = 0;

	while ((c = getopt(argc, argv, "i:p:t:m:s:")) != -1) {
		switch (c) {
			case 'i':
				gb_inputFileName = optarg;
				break;
			case 'p':
				if (strcmp(optarg,"T") == 0) {
					gb_proteinsFlag = 1;
				} else if (strcmp(optarg,"F") == 0) {
					gb_proteinsFlag = 0;
				} else {
					perror("[FORMATDB_ERROR] Unknown formatdb option\n");
					return;
				}
				break;
			case 's':
				if (strcmp(optarg,"T") == 0) {
					gb_formatDbScoresFlag = 1;
				} else if (strcmp(optarg,"F") == 0) {
					gb_formatDbScoresFlag = 0;
				} else {
					perror("[FORMATDB_ERROR] Unknown formatdb scoring option\n");
					return;
				}
				break;			
			case 't':
				gb_threshold = atoi(optarg);
				break;
			case 'm':
				memcpy(gb_psmName ,optarg, strlen(optarg) * sizeof(char));
				break;
			default:
				perror("[FORMATDB_ERROR] Unknown option\n");
				break;
		}
	}
}

int main(int argc, char ** argv) {
	long int i;
	FILE * dbf = NULL;
	FILE * of = NULL;
	
	unsigned char ** formatedSec = NULL;
	unsigned int * sizeSec = NULL;

	long * fdb = NULL;

	// parse parameters
	parseInputParameters(argc, argv);

	// fetching database
	fdb = fetchingDBwithStrand(gb_inputFileName, (gb_proteinsFlag == 1)? -1 : 0);

	// parse database into hashmap
	dbf = fopen(gb_inputFileName,"r");
	if (!dbf) {
		perror("[FORMATDB_ERROR] Unable to open input file\n");
		return -1;
	}

	// Allocate memory for formated secuences and sizes (incluida la Aux antes del formateo)
	formatedSec = (unsigned char **) malloc (sizeof(unsigned char *) * fdb[1]);
	sizeSec = (unsigned int *) calloc (fdb[1], sizeof(unsigned int));
	
	/*
	for(i = 0; i < fdb[1]; i++)
		formatedSec[i] = (unsigned char *) calloc (NUCL_BIT_SIZE, sizeof(unsigned char));
	*/
		
	//Read the input file and change format.
	if (parseFasta(dbf, gb_proteinsFlag, sizeSec, formatedSec, gb_threshold, gb_psmName, gb_formatDbScoresFlag) == -1) {
		perror("[FORMATDB_ERROR] Unable to parse FASTA file\n");
		if (dbf) {
			fclose(dbf);
		}
		return -1;
	}

	if (dbf) {
		fclose(dbf);
	}

	// write data into binary file
	gb_outputFileName = (char *) calloc (strlen(gb_inputFileName) + 5, sizeof(char));
	if (!gb_outputFileName) 
	{
		perror("[FORMATDB_ERROR] Unable to allocate memory\n");
		return -1;
	}
	
	memcpy(gb_outputFileName, gb_inputFileName, strlen(gb_inputFileName) * sizeof(char));
	strcat(gb_outputFileName, ".fdb");

	of = fopen(gb_outputFileName,"wb");
	if (!of) {
		perror("[FORMATDB_ERROR] Unable to open output file\n");
		return -1;
	}

	//NUMERO DE SECUENCIAS FORMATEADAS.
	fwrite(&fdb[1], sizeof(long int), 1, of);
	// SECUENCIAS.
	for (i = 0; i < fdb[1]; i++)
	{
		fwrite(&sizeSec[i], sizeof(unsigned int), 1, of);
		fwrite(formatedSec[i], sizeof(unsigned char), NUCL_BIT_SIZE, of);
	}
	
	if(sizeSec)
		free(sizeSec);
		
	if(formatedSec)
	{
		for(i = 0; i < fdb[1]; i++)
		{
			if(formatedSec[i])
				free(formatedSec[i]);
		}
		
		free(formatedSec);
	}
	
	if (of) {
		fclose(of);
	}

	if (fdb) {
		free((void *)fdb);
	}

	return 0;
}
