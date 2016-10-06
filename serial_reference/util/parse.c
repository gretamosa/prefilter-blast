#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>

#include "parse.h"

void parse (char * fn, char ** ofs, unsigned long int * ofl, unsigned int * ofn, unsigned ** ofla, char * bi, unsigned short int ws) {
	FILE * f = NULL;
	char cAux = 0;

	if (fn && bi) {
		// PROTEINS CASE
		if (memcmp(bi,"blastp", 6) == 0) {
			f = fopen(fn,"r+");

			do {
				parseSequence(f, gb_protIndx, PROT_ALP_SIZE, ofs, ofl, ofn, ofla, ws);

				cAux = fgetc(f);
				if (cAux != EOF) {
					fseek(f,-2, SEEK_CUR);
				}
			} while (cAux != EOF);


			if (f) {
				fclose(f);
			}
		}
	}
}

void parseSequence(FILE * fParam, const unsigned short int * indx, unsigned short int alpsize, char ** data, unsigned long int * ofl, unsigned int * ofn, unsigned int ** ofla, unsigned short int ws) {
	char * gid = NULL, *seq = NULL;
	unsigned int sl = 0, si = 0;

	if (fParam && indx) {
		// Read GID header
		gid = readGid(fParam);
		if (gid) {
			free((void *) gid);
		}

		// Read sequence
		seq = readSeq(fParam);

		// allocate memory
		sl = strlen(seq);
		si = (*ofl);
		(*ofl) += sl;
		(*ofn)++;

		if (*data) {
			(*data) = (char *) realloc ((*data), ((*ofl) + 1) * sizeof(char));
			(*ofla) = (unsigned int *) realloc ((*ofla), (*ofn) * sizeof(unsigned int));
		} else {
			(*data) = (char *) malloc (((*ofl) + 1) * sizeof(char));
			(*ofla) = (unsigned int *) malloc ((*ofn) * sizeof(unsigned int));
		}

		(*ofla)[(*ofn)-1] = sl;
		memcpy((*data) + si, seq, sl);
		(*data)[(*ofl)] = '\0';

		// Process lmer
		//seeding(seq, indx, alpsize, ws, dl, data);

		if (seq) {
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
