#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>

#include "parse.h"

void parseFDB (char * ofn, char ** sec, unsigned long int * fl, long int * fn, unsigned int ** flsec, char * bi)
{
	int i,j;
	FILE * f = NULL;
	
	f = fopen(ofn,"rb");
	
	fread(fn, sizeof(long int), 1, f);
	//fprintf(stderr,"numero de cadenas %ld\n", (*fn));
	
	// PROTEINS CASE
	if (memcmp(bi,"blastp", 6) == 0)
	{
		(*fl) = (*fn) * PROTDBSIZE;
		
		(*flsec) = (unsigned int *) malloc (sizeof(unsigned int) * (*fn));
		(*sec) = (char *) malloc (sizeof(char) * (*fl));
		
		//Se leen todas las cadenas, longitud - secuencia.
		for(i = 0; i < (*fn); i++)
		{
			fread(&((*flsec)[i]), sizeof(unsigned int), 1, f);
			fread(&((*sec)[i * PROTDBSIZE]), sizeof(char), PROTDBSIZE, f);
		}
	}
	else
	{
		fprintf(stderr,"NOT IMPLEMENTED YET \n");
		exit(0);
	}
	
	fclose(f);
	return;
}

void parse(char * fn, dictionary ** dict, char ** sec, unsigned long int * ofl, unsigned int * ofn, unsigned ** ofla, char * bi, unsigned short int ws, unsigned short int threshold, char * sm_name, unsigned short int fs)
{
	FILE * f = NULL;
	char cAux = 0;

	if (fn && bi)
	{

		if (memcmp(bi,"blastp", 6) == 0) // PROTEINS CASE
		{
			f = fopen(fn,"r");

			do
			{
				parseSequencePerWarp(f, gb_protIndx, PROT_ALP_SIZE, dict, sec, ofl, ofn, ofla, ws, threshold, sm_name, fs);

				cAux = fgetc(f);
				if (cAux != EOF)
				{
					fseek(f,-2, SEEK_CUR);
				}
			} while (cAux != EOF);


			if (f)
			{
				fclose(f);
			}
		}
		else if (memcmp(bi,"blastn", 6) == 0)
		{
			f = fopen(fn,"r");

			do
			{
				parseSequencePerWarp(f, gb_protIndx, PROT_ALP_SIZE, dict, sec, ofl, ofn, ofla, ws, threshold, sm_name, 0);

				cAux = fgetc(f);
				if (cAux != EOF)
				{
					fseek(f,-2, SEEK_CUR);
				}
			} while (cAux != EOF);


			if (f)
			{
				fclose(f);
			}
		}
	}
}

void parseOut(FILE * fout, dictionary * db, int indx)
{
	if (indx != -1)
	{
		fprintf(fout,"%s\n%s\n",db[indx].gid,db[indx].seq);
	}
}

void parseSequence(FILE * fParam, const unsigned short int * indx, unsigned short int alpsize, dictionary ** dict, char ** data, unsigned long int * ofl, unsigned int * ofn, unsigned int ** ofla, unsigned short int ws, unsigned short int threshold, char * sm_name, unsigned short int fs) {
	char * gid = NULL, *seq = NULL;
	unsigned int gl = 0, sl = 0, si = 0;

	char * neighbours = NULL;
	unsigned int neighbours_length = 0;

	if (fParam && indx) {
		// Read GID header
		gid = readGid(fParam);
		gl = strlen(gid);

		// Read sequence
		seq = readSeq(fParam);
		sl = strlen(seq);

		// Include neighbours lmers
		if (fs == 1)
		{
			//seeding(seq, sl, indx, alpsize, ws, threshold, sm_name, &neighbours, &neighbours_length);
		}

		// allocate memory
		si = (*ofl);
		(*ofl) += sl;
		(*ofn)++;

		//DICCIONARIO
		if (!(*dict)) // malloc
		{
			(*dict) = (dictionary *) malloc ((*ofn) * sizeof(dictionary));
			(*ofla) = (unsigned int *) malloc ((*ofn) * sizeof(unsigned int));
		}
		else // realloc
		{
			(*dict) = (dictionary *) realloc ((*dict), (*ofn) * sizeof(dictionary));
			(*ofla) = (unsigned int *) realloc ((*ofla), (*ofn) * sizeof(unsigned int));
		}
		
		//No hace falta hacer callocs, se puede poner \0 al final y usar mallocs.
		(*dict)[(*ofn)-1].gid = (char *) calloc (gl + 1, sizeof(char));
		memcpy((*dict)[(*ofn)-1].gid, gid, gl);
		(*dict)[(*ofn)-1].seq = (char *) calloc (sl + 1, sizeof(char));		
		
		memcpy((*dict)[(*ofn)-1].seq, seq, sl);

		(*ofla)[(*ofn)-1] = sl;
		
		//SECUENCIA TOTAL
		if (*data) {
			(*data) = (char *) realloc ((*data), ((*ofl) + 1) * sizeof(char));
		} else {
			(*data) = (char *) malloc (((*ofl) + 1) * sizeof(char));
		}
		memcpy((*data) + si, seq, sl);
		(*data)[(*ofl)] = '\0';

		if (gid)
		{
			free((void *) gid);
		}
		if (seq)
		{
			free((void *)seq);
		}
	}
}

void parseSequencePerWarp(FILE * fParam, const unsigned short int * indx, unsigned short int alpsize, dictionary ** dict, char ** data, unsigned long int * ofl, unsigned int * ofn, unsigned int ** ofla, unsigned short int ws, unsigned short int threshold, char * sm_name, unsigned short int fs) {
	char * gid = NULL, *seq = NULL;
	unsigned int gl = 0, sl = 0, sl_aux = 0, sl_aux2 = 0, si = 0;

	char * neighbours = NULL;
	unsigned int neighbours_length = 0;

	if (fParam && indx) {		
		// Read GID header
		gid = readGid(fParam);
		gl = strlen(gid);

		// Read sequence
		seq = readSeq(fParam);
		sl_aux2 = sl = strlen(seq);
		
		// ADDED ONE SEQUENCE PER WARP
		sl_aux = 256;
		
		while (sl > sl_aux)
		{
			sl_aux += 256;
		}
		sl = sl_aux;

		// Include neighbours lmers
		if (fs == 1)
		{
			//seeding(seq, sl, indx, alpsize, ws, threshold, sm_name, &neighbours, &neighbours_length);
		}

		// allocate memory
		si = (*ofl);
		(*ofl) += sl;
		(*ofn)++;

		//DICCIONARIO
		if (!(*dict)) // malloc
		{
			(*dict) = (dictionary *) malloc ((*ofn) * sizeof(dictionary));
			(*ofla) = (unsigned int *) malloc ((*ofn) * sizeof(unsigned int));
		}
		else // realloc
		{
			(*dict) = (dictionary *) realloc ((*dict), (*ofn) * sizeof(dictionary));
			(*ofla) = (unsigned int *) realloc ((*ofla), (*ofn) * sizeof(unsigned int));
		}
		
		//No hace falta hacer callocs, se puede poner \0 al final y usar mallocs.
		(*dict)[(*ofn)-1].gid = (char *) calloc (gl + 1, sizeof(char));
		memcpy((*dict)[(*ofn)-1].gid, gid, gl);
		
		// el diccionario no incluye 'X'
		(*dict)[(*ofn)-1].seq = (char *) calloc (sl_aux2 + 1, sizeof(char));
		memcpy((*dict)[(*ofn)-1].seq, seq, sl_aux2 * sizeof(char));

		(*ofla)[(*ofn)-1] = sl;
		
		//SECUENCIA TOTAL
		if (*data) {
			(*data) = (char *) realloc ((*data), (*ofl) * sizeof(char));
		} else {
			(*data) = (char *) malloc ((*ofl) * sizeof(char));
		}
		memset((*data) + si, '*', sl * sizeof(char));
		memcpy((*data) + si, seq, sl_aux2 * sizeof(char));

		if (gid)
		{
			free((void *) gid);
		}
		if (seq)
		{
			free((void *)seq);
		}
	}
}

char * readGid(FILE * fParam) {
	char * str = NULL;
	char cAux = 0;

	unsigned int gc = 0;
	unsigned long int initial_file_pos = 0;
	unsigned short int ret = 0;

	if (fParam) {
		// 1. Obtain the initial position of the file pointer
		initial_file_pos = ftell(fParam);
		// 2. Obtain the length of the GID
		cAux = fgetc(fParam);
		while (cAux != '\n' && cAux != EOF) {
			gc++;
			cAux = fgetc(fParam);
		}
		// 2. Allocate memory
		str = (char *) calloc ((gc + 2), sizeof(char));
		if (!str) {
			perror("Unable to allocate memory\n");
			ret = 1;
		}
		if (ret == 0) {
			// 3. Back to the initial position
			fseek(fParam, initial_file_pos, SEEK_SET);
			// 4. Read the GID
			fgets(str, gc+1, fParam);
		}
		// 5. Jump to sequence position
		fgetc(fParam);
	} else {
		perror("Input File Null Pointer Exception\n");
	}

	return str;
}

char * readSeq(FILE * fParam) {
	char * str = NULL;
	char cAux = 0;

	unsigned int sc = 0, scl = 0, sct = 0;
	unsigned int * sc_ar = NULL;
	unsigned int j=0, i=0;
	unsigned long int initial_file_pos = 0;
	unsigned short int ret = 0;

	if (fParam) {
		// 1. Obtain the initial position of the file pointer
		initial_file_pos = ftell(fParam);
		// 2. Obtain the number of sequence blocks (in FASTA format)
		cAux = fgetc(fParam);
		while (cAux != '>' && cAux != EOF) {
			if (cAux == '\n') {
				scl++;
			}
			cAux = fgetc(fParam);
		}
		// 3. Back to the initial position
		fseek(fParam, initial_file_pos, SEEK_SET);
		// 4. Allocate memory
		sc_ar = (unsigned int *) malloc (scl * sizeof(unsigned int));
		if (!sc_ar) {
			perror("Unable to allocate memory.\n");
			ret = 1;
		}
		if (ret == 0) {
			// 5. Obtain the length of each sequence block
			cAux = fgetc(fParam);
			while (cAux != '>' && cAux != EOF) {
				sc++;
				sct++;
				if (cAux == '\n') {
					sc_ar[j] = sc;
					sc = 0;
					j++;
				}
				cAux = fgetc(fParam);
			}
			// 6. Allocate memory
			str = (char *) calloc ((sct + 1), sizeof(char));
			if (!str) {
				perror("Unable to allocate memory.\n");
				ret = 1;
			}

			if (ret == 0) {
				// 7. Back to the initial position
				fseek(fParam, initial_file_pos, SEEK_SET);
				// 8. Read the GID
				i = 0;
				for (j=0;j<scl;j++) {
					fgets(str+i, sc_ar[j], fParam);
					fgetc(fParam);
					i += sc_ar[j]-1;
				}
			}
			fgetc(fParam);
			// 9. Free memory
			if (sc_ar) {
				free((void *)sc_ar);
			}
		}
	} else {
		perror("Input File Null Pointer Exception\n");
	}

	return str;
}
