#include "utils.h"

// =============================== SUBSTITUTION MATRIX METHODS PATH ===============================
int getMatrixScore(const char * alpParam, const int alpsize, const int * psmParam, const int * indxParam, char firstChar, char secondChar) {
	int indxFirst = 0, indxLast = 0, ret = -9999;

	if (alpParam && psmParam) {
		indxFirst = getHashIndx(firstChar, indxParam);
		indxLast = getHashIndx(secondChar, indxParam);
		indxLast = (alpsize * indxLast) + indxFirst;
		ret = psmParam[indxLast];
	}

	return ret;
}

// =============================== HASHING METHODS PATH ===============================
unsigned long juan_hash (char * name, int length) {
	int i=0;
	unsigned long h =0;

	for (i=0;i<length;i++) {
		if (length == 3) {
			h += getHashIndx(name[i], PROTINDX) * pow(PROTSIZE,length - (i + 1));
		} else if (length == 11) {
			h += getNuclHashIndx(name[i], NUCLINDX) * pow(NUCLSIZE, length - (i + 1));
		}
	}

	return h;
}

int getNuclHashIndx(char charParam, const int * indxParam) {
	int indx = 0;

	if (charParam == 65) //A
	{
		indx = indxParam[0];
	} 
	else if (charParam == 67) //C
	{
		indx = indxParam[1];
	}
	else if (charParam == 71) //G
	{
		indx = indxParam[2];
	}
	else if (charParam == 84) //T
	{
		indx = indxParam[3];
	}
	
	/*
	else if (charParam == 85) // U
	{
		indx = indxParam[4];
	}
	*/
 
	return indx;
}

// TODO: NUCLEOTIDES HASH INDEX FUNCTION
int getHashIndx(char charParam, const int * indxParam) {
	int indx = 0;

	if (charParam >= 65 && charParam <= 90) {
		charParam -= 65;
		indx = indxParam[(int)charParam];
	} else if (charParam == 42) {
		indx = indxParam[26];
	} else {
		indx = indxParam[23];
	}

	return indx;
}

// =============================== FETCHING METHOD =====================================
long * fetchingDBwithStrand(char * dbFileName, int nucleotideStrand) {
	FILE * f = NULL;
	char cAux = 0;
	long cc = 0, sc = 0;
	long * acc = NULL;

	int ret = 0;

	if (dbFileName) {
		acc = (long *) malloc (2 * sizeof(long));
		if (!acc) {
			perror("[GPUBLAST_ERROR] Unable to allocate memory (fetchingDBwithStrand)\n");
			return NULL;
		}
		f = fopen(dbFileName, "r");
		if (!f) {
			perror("[GPUBLAST_ERROR] Unable to open sequence file\n");
			return NULL;
		}
		do  {
			ret = fread(&cAux, sizeof(char), 1, f);
			while (cAux != '\n' && ret != 0) {
				ret = fread(&cAux, sizeof(char), 1, f);
			}
			while (cAux != '>' && ret != 0) {
				if (cAux != '\n') {
					cc++;
				}
				ret = fread(&cAux, sizeof(char), 1, f);
			}
			sc++;
			
			//ESTO PA QUE!!
			if (nucleotideStrand != -1)
			{
				//sc++;
			}
		} while (ret != 0);
		fclose(f);

		acc[0] = cc;
		acc[1] = sc;
	}

	return acc;
}

int parseFasta (FILE * fParam, int protFlag, unsigned int * tamSec, unsigned char ** formatedSec, int threshold, char * psmName, int fdb_score) {
	char * gid = NULL, * seq = NULL;
	char cAux = 0;
	unsigned long qid = 0;

	if (fParam) {
		do {
			if (parseSequence(fParam, gid, seq, protFlag, qid, tamSec, formatedSec, threshold, psmName, fdb_score) == -1) {
				perror("[FORMATDB_ERROR] Error while parsing sequence\n");
				return -1;
			}

			qid++;

			cAux = fgetc(fParam);
			if (cAux != EOF) {
				fseek(fParam,-2, SEEK_CUR);
			}

			if (gid) {
				free((void *)gid);
			}
			if (seq) {
				free((void *)seq);
			}

		} while (cAux != EOF);

		return 0;
	}

	return -1;
}

// =============================== SEEDING SEQUENCES =============================
long * seedingWithScores(char * lmer, int size, const char * alp, const int alpsize, const int * psm, const int * indx, int threshold, int * neighbour_counter, int ** neighbour_score) {
	long * neighbour_list = NULL;
	int * scs = NULL;
	int i=0, j=0, c=0, score=-1;	
	char cAux = 0;
	
	if (lmer) {
		// initialize
		*neighbour_counter = 0;
		
		// allocate memory
		neighbour_list = (long *) malloc ((size * alpsize) * sizeof(long));
		if (!neighbour_list) 
		{
			perror("Unable to allocate memory\n");
			return NULL;
		}
		
		(*neighbour_score) = (int *) malloc ((size * alpsize) * sizeof(int));
		if (!(*neighbour_score)) 
		{
			perror("Unable to allocate memory\n");
			return NULL;
		}
		
		scs = (int *) malloc (size * sizeof(int));
		if (!scs) {
			perror("Unable to allocate memory\n");
			return NULL;
		}
		
		// obtain score equivalent score
		for (i=0;i<size;i++) {
			scs[i] = getMatrixScore(alp, alpsize, psm, indx, lmer[i], lmer[i]);
		}
		
		for (i=0;i<size;i++) 
		{
			c = 0;
			for (j=0;j<size;j++) 
			{
				if (j != i){
					c += scs[j];
				}
			}
						
			for (j=0;j < alpsize - 1;j++) 
			{									
				score = getMatrixScore(alp, alpsize, psm, indx, lmer[i], alp[j]);					
				cAux = lmer[i];
				lmer[i] = alp[j];
				score += c;

				if (score >= threshold) 
				{
					neighbour_list[(*neighbour_counter)] = juan_hash(lmer, size);
					(*neighbour_score)[(*neighbour_counter)] = score;
					(*neighbour_counter)++;
				}
				
				if (cAux != 0) 
				{
					lmer[i] = cAux;
				}
			}
		}
		
		if (scs) {
			free((void *) scs);
		}
	}
	
	return neighbour_list;
}

long * seeding(char * lmer, int size, const char * alp, const int alpsize, const int * psm, const int * indx, int threshold, int * neighbour_counter) {
	long * neighbour_list = NULL;
	int * scs = NULL;
	int i=0, j=0, c=0, score=-1;	
	char cAux = 0;
	
	if (lmer) {
		// initialize
		*neighbour_counter = 0;
		
		// allocate memory
		neighbour_list = (long *) malloc ((size * alpsize) * sizeof(long));
		if (!neighbour_list) 
		{
			perror("Unable to allocate memory\n");
			return NULL;
		}
		
		scs = (int *) malloc (size * sizeof(int));
		if (!scs) {
			perror("Unable to allocate memory\n");
			return NULL;
		}
		
		// obtain score equivalent score
		for (i=0;i<size;i++) {
			scs[i] = getMatrixScore(alp, alpsize, psm, indx, lmer[i], lmer[i]);
		}
		
		for (i=0;i<size;i++) 
		{
			c = 0;
			for (j=0;j<size;j++) 
			{
				if (j != i){
					c += scs[j];
				}
			}
						
			for (j=0;j < alpsize - 1;j++) 
			{									
				score = getMatrixScore(alp, alpsize, psm, indx, lmer[i], alp[j]);					
				cAux = lmer[i];
				lmer[i] = alp[j];
				score += c;

				if (score >= threshold) 
				{
					neighbour_list[(*neighbour_counter)] = juan_hash(lmer, size);
					(*neighbour_counter)++;
				}
				
				if (cAux != 0) 
				{
					lmer[i] = cAux;
				}
			}
		}
		
		if (scs) {
			free((void *) scs);
		}
	}
	
	return neighbour_list;
}

// =============================== PARSING SEQUENCES =============================
int parseSequence(FILE * fParam, char * gid, char * seq, int protFlag, int qid, unsigned int * tamSec, unsigned char ** formatedSec, int threshold, char * psmName, int fdb_score) {
	if (fParam) {
		// Read GID header
		gid = readGID_lmc(fParam);
		if (!gid) {
			perror("[FORMATDB_ERROR] Unable to read GID\n");
			return -1;
		}

		// Read sequence
		seq = readSeq_lmc(fParam);
		if (!seq) {
			perror("[FORMATDB_ERROR] Unable to read Sequence\n");
			if (!gid) {
				free((void *)gid);
			}
			return -1;
		}

		// Process lmer
		if (processLmersNucl(formatedSec, tamSec, seq, 11, qid) == -1) {
			perror("[FORMATDB_ERROR] Unable to process lmer sequence\n");
			if (!gid) {
				free((void *)gid);
			}
			if (!seq) {
				free((void *)seq);
			}
			return -1;
		}
	}

	return 0;
}

int processLmersProtWithScores(char ** formatedSec , unsigned int * tamSec, char * seq, int lmerl, int qid, int threshold, char * psmName) {
	int i=0, j = 0;
	int sl = 0;

	unsigned long guid = 0;
	int neighbour_counter = 0;
	long * neighbour_list = NULL;
	int * neighbour_score = NULL;
	int score = 0;

	if (seq) {
		sl = strlen(seq);
		tamSec[qid] = sl;
		
		for (i = 0; i <= (sl-lmerl); i++) 
		{
			guid = juan_hash(seq+i, lmerl);
			
			//Se calcula el score del LMER
			score = 0;
			for (j=0;j<lmerl;j++) 
			{
				if(memcmp(psmName, "BLOSUM62", 8 * sizeof(char)) == 0)
				{
					score += getMatrixScore(Blosum62ALP, Blosum62ALP_SIZE, Blosum62PSM, Blosum62INDX, seq[i+j], seq[i+j]);
				}
			}
			
			formatedSec[qid][guid] = score;
			
			//Añadir los seeds vecinos de este lmer.
			if(memcmp(psmName, "BLOSUM62", 8 * sizeof(char)) == 0)
			{
				neighbour_list = seedingWithScores(seq+i, lmerl, Blosum62ALP, Blosum62ALP_SIZE, Blosum62PSM, Blosum62INDX, threshold, &neighbour_counter, &neighbour_score);
			}
			
			for(j = 0; j < neighbour_counter; j++)
			{
				//Si el vecino es mejor que el actual, se cambia al valor del vecino (solo se cumplirá cuando el elemento no esté inicializado)
				//Esto consigue que si se tiene el elemento real y el vecino en la secuencia, te quedes siempre con el real.
				if(neighbour_score[j] > formatedSec[qid][neighbour_list[j]])
				{
					formatedSec[qid][neighbour_list[j]] = neighbour_score[j];
				}
			}
			
			if(neighbour_list) 
			{
				free((void *)neighbour_list);
			}				
			if(neighbour_score)
			{
				free((void *)neighbour_score);
			}
		}		

		return 0;
	}

	return -1;
}

int processLmersProtWithHits(char ** formatedSec , unsigned int * tamSec, char * seq, int lmerl, int qid, int threshold, char * psmName) {
	int i=0, j = 0;
	int sl = 0;

	unsigned long guid = 0;
	int neighbour_counter = 0;
	long * neighbour_list = NULL;

	if (seq) {
		sl = strlen(seq);
		tamSec[qid] = sl;
		
		for (i = 0; i <= (sl-lmerl); i++) 
		{
			guid = juan_hash(seq+i, lmerl);			
			formatedSec[qid][guid] = 1;
			
			//Añadir los seeds vecinos de este lmer.
			if(memcmp(psmName, "BLOSUM62", 8 * sizeof(char)) == 0)
			{
				neighbour_list = seeding(seq+i, lmerl, Blosum62ALP, Blosum62ALP_SIZE, Blosum62PSM, Blosum62INDX, threshold, &neighbour_counter);
			}
			
			for(j = 0; j < neighbour_counter; j++)
			{
				formatedSec[qid][neighbour_list[j]] = 1;
			}
			
			if(neighbour_list) 
			{
				free((void *)neighbour_list);
			}
		}		

		return 0;
	}

	return -1;
}

int processLmersNucl(unsigned char ** formatedSec , unsigned int * tamSec, char * seq, int lmerl, int qid) {
	unsigned long i = 0;
	int sl = 0;
	unsigned char * SecBeforeFormat = NULL;
	unsigned long guid = 0;

	if (seq) 
	{
		SecBeforeFormat = (unsigned char *) calloc (NUM_NUCLEOTIDE_POSITIONS, sizeof(unsigned char));
		
		sl = strlen(seq);
		tamSec[qid] = sl;

		for (i = 0; i <= (sl-lmerl); i++) {
			guid = juan_hash(seq+i, lmerl);
			SecBeforeFormat[guid] = 1;
			//printf("se han introducido un elemento en la posicion %ld\n", guid);
		}
		
		// Se formatea la cadena.
		formatedSec[qid] = transformByteToBit(SecBeforeFormat);
		
		/*
		//-------------------------------------------------------------------------------
		//Test de que se ha realizado correctamente
		//Se pone a cero el array.
		unsigned char result = 0;
		memset(SecBeforeFormat, 0, NUM_NUCLEOTIDE_POSITIONS * sizeof(char));
		
		for(i = 0; i < NUM_NUCLEOTIDE_POSITIONS; i++)
		{				
			result = fromGuidToResult(formatedSec[qid], i);
			if(result == 1)
				fprintf(stderr,"se han introducido un elemento en la posicion %ld\n", i);
		}
		//-------------------------------------------------------------------------------
		exit(0);
		*/ 
		
		if(SecBeforeFormat)
			free(SecBeforeFormat);
			
		return 0;
	}

	return -1;
}

// =============================== PARSING SEQUENCES (GID - SEQ)  ==============
char * readGID_lmc(FILE * fParam) {
	char * str = NULL;
	char cAux = 0;

	int gc = 0;
	long int initial_file_pos = 0;

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
			perror("[FORMATDB_ERROR] Unable to allocate memory (readGID_lmc)\n");
			return NULL;
		}
		// 3. Back to the initial position
		fseek(fParam, initial_file_pos, SEEK_SET);
		// 4. Read the GID
		fgets(str, gc+1, fParam);
		// 5. Jump to sequence position
		fgetc(fParam);
	} else {
		perror("[FORMATDB_ERROR] NullPointerException\n");
	}

	return str;
}

char * readSeq_lmc(FILE * fParam) {
	char * str = NULL;
	char cAux = 0;

	int sc = 0, scl = 0, sct = 0;
	int * sc_ar = NULL;
	int j=0, i=0;
	long int initial_file_pos = 0;

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
		sc_ar = (int *) malloc (scl * sizeof(int));
		if (!sc_ar) {
			perror("[FORMATDB_ERROR] Unable to allocate memory (readSeq_lmc)\n");
			return NULL;
		}
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
			perror("[FORMATDB_ERROR] Unable to allocate memory (readSeq_lmc)\n");
			return NULL;
		}
		// 7. Back to the initial position
		fseek(fParam, initial_file_pos, SEEK_SET);
		// 8. Read the GID
		i = 0;
		for (j=0;j<scl;j++) {
			fgets(str+i, sc_ar[j], fParam);
			fgetc(fParam);
			i += sc_ar[j]-1;
		}
		fgetc(fParam);
		// 9. Free memory
		if (sc_ar) {
			free((void *)sc_ar);
		}
	} else {
		perror("[FORMATDB_ERROR] NullPointerException\n");
	}

	return str;
}
