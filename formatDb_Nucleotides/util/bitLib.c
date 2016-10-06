#include "bitLib.h"

// Función que transforma el valor de 8 chars (8 elementos en bytes) a un byte (8 elementos en bits) 
unsigned char valueToBit(unsigned char * values)
{
	int i;
	int base = 1;
	unsigned char bits = 0;
	
	for(i = 0; i < BYTE_SIZE; i++)
	{
		//printf("%d + %d\n", bits, values[i] * base);
		bits += values[i] * base;
		base = base * 2;
	}
	
	//printf("byte valor %d \n", bits);
	
	return bits;
}

// Función que saca de un byte el valor de un bit determinado por la posición.
unsigned char maskAnd(unsigned char value, int maskPos)
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

// Función que transforma una cadena de descomprimidas a comprimidas.
unsigned char * transformByteToBit(unsigned char * byteValues)
{
	int i, j;
	int counter = 0;
	unsigned char byte[BYTE_SIZE];
	unsigned char * bitString = NULL;
	
	bitString = (unsigned char *) calloc (NUCL_BIT_SIZE, sizeof(unsigned char));
	if(!bitString)
	{
		perror("[FormatDb Nucl] Unable to allocate memory (transformByteToBit)\n");
		return NULL;
	}
	
	//Se recorre
	for(i = 0; i < NUCL_BIT_SIZE; i++)
	{
		// Big Endian
		for(j = 0; j < BYTE_SIZE; j++)
		{
			byte[j] = byteValues[counter];
			counter ++;
		}
		
		bitString[i] = valueToBit(byte);
	}
	
	return bitString;
}

// Función que transforma una cadena de comprimidas a descomprimidas.
unsigned char * transformBitToByte(unsigned char * bitValues)
{
	int i, j, counter = 0;
	unsigned char * byteString = NULL;
	
	byteString = (unsigned char *) calloc (NUM_NUCLEOTIDE_POSITIONS, sizeof(unsigned char));
	if(!byteString)
	{
		perror("[FormatDb Nucl] Unable to allocate memory (transformBitToByte)\n");
		return NULL;
	}
	
	for(i = 0; i < NUCL_BIT_SIZE; i++)
	{
		for(j = 0; j < BYTE_SIZE; j++)
		{
			if(((counter % BYTE_SIZE) == 0) && (counter != 0))
			//printf("\t");
			byteString[counter] = maskAnd(bitValues[i], j);
			//printf("%d",byteString[counter]);
			counter++;
		}
	}
		
	//printf("\n");
	
	return byteString;
}

// Función que devuelve un valor de la cadena de bits a partir de un GUID
unsigned char fromGuidToResult(unsigned char * bitValues, long int guid)
{
	int block = guid / BYTE_SIZE;
	int position = guid % BYTE_SIZE;
	//printf("guid %ld, block %d, pos %d -> ", guid, block, position);
	
	return maskAnd(bitValues[block], position);
}
/*
// Pruebita pequeña de uso.
int main(int argc, char ** argv) {
	int i, k, test = 10;
	unsigned char valoresIniciales[NUM_NUCLEOTIDE_POSITIONS];
	unsigned char * valoresFinales = NULL;
	unsigned char * valoresComprimidos = NULL;
	unsigned char result = 0;
	
	//Comprobación de que se traducen correctamente las cadenas en ambos sentidos.
	printf("\nPRIMERA PRUEBA\n");
	
	for(k = 1; k < test; k++)
	{
		printf("\n------------ RONDA %d ------------\n",k);
		
		printf("\nINICIALMENTE\n");
		//Se inicializan los valores Ini y Fin.
		for(i = 0; i < NUM_NUCLEOTIDE_POSITIONS; i++)
		{
			if(((i % BYTE_SIZE) == 0) && (i != 0))
				printf("\t");
				
			if((i % k) == 0)
				valoresIniciales[i] = 1;
			else
				valoresIniciales[i] = 0;
		
			printf("%d",valoresIniciales[i]);
		}
		
		printf("\n");
		
		//Se comprimen los resultados.
		valoresComprimidos = transformByteToBit(valoresIniciales);
		
		printf("\nFINALMENTE\n");
		valoresFinales = transformBitToByte(valoresComprimidos);
		
		if(valoresFinales)
			free((void *)valoresFinales);
			
		//La última no se libera para la siguiente prueba	
		if(k != (test-1))
		{
			if(valoresComprimidos)
				free((void *)valoresComprimidos);
		}
	}
	
	
	//Comprobación de la función que obtiene valores a partir de un GUID.
	printf("\nSEGUNDA PRUEBA\n");
	
	for(i = 0; i < NUM_NUCLEOTIDE_POSITIONS; i++)
	{
		if((i % BYTE_SIZE) == 0)
			printf("\n");
			
		result = fromGuidToResult(valoresComprimidos, i);
		printf("%d\n", result);
	}
	
	if(valoresComprimidos)
		free((void *)valoresComprimidos);

	return 0;
}
*/
