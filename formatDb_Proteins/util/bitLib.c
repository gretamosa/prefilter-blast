#include "bitLib.h"

/**
* \fn unsigned char valueToBit(unsigned char * values)
* \brief Función que transforma el valor de 8 chars (8 elementos en bytes) a un byte (8 elementos en bits).
* \param[in] values, valores a transformar.
* \return byte con los valores comprimidos.
*/
unsigned char valueToBit(unsigned char * values)
{
	int i;
	int base = 1;
	unsigned char bits = 0;
	
	for(i = 0; i < BYTE_SIZE; i++) {
		bits += values[i] * base;
		base = base * 2;
	}
	
	return bits;
}

/**
* \fn unsigned char maskAnd(unsigned char value, int maskPos)
* \brief Función que saca de un byte el valor de un bit determinado por la posición.
* \param[in] value, Valor del que se debe sacar el dato.
* \param[in] maskPos, Máscara binaria para obtener el valor correcto.
* \return el valor del bit (1 o 0).
*/
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

/**
* \fn unsigned char * transformByteToBit(unsigned char * byteValues, unsigned int * lmerIndex, int numLmer)
* \brief Función que transforma una cadena de valores descomprimidos (bytes) a valores comprimidos (bits).
* \param[in] byteValues, Valores en Bytes a comprimir.
* \param[in] lmerIndex, array con los bloques que se deben comprimir.
* \param[in] numLmer, Número de Lmers que se deben comprimir.
* \return cadena ce valores comprimidos.
*/
unsigned char * transformByteToBit(unsigned char * byteValues, unsigned int * lmerIndex, int numLmer)
{
	int i, j;
	//int counter = 0;
	unsigned char byte[BYTE_SIZE];
	unsigned char * bitString = NULL;
	char callFlag[NUCL_BIT_SIZE];
	unsigned int block = 0;
	
	bitString = (unsigned char *) calloc (NUCL_BIT_SIZE, sizeof(unsigned char));
	if(!bitString) {
		perror("[FormatDb Nucl] Unable to allocate memory (transformByteToBit)\n");
		return NULL;
	}
	
	memset(callFlag, 0, sizeof(char) * NUCL_BIT_SIZE);
	
	for (i = 0; i <= numLmer; i++)
	{
		block = lmerIndex[i];
		
		if(callFlag[block] == 0)
		{
			for(j = 0; j < BYTE_SIZE; j++)
			{
				byte[j] = byteValues[(block * BYTE_SIZE) + j];
			}
			
			bitString[block] = valueToBit(byte);
			callFlag[block] = 1;
		}
	}
	
	return bitString;
}

/**
* \fn unsigned char * transformBitToByte(unsigned char * bitValues)
* \brief Función que transforma una cadena de valores comprimidos a valores descomprimidos.
* \param[in] bitValues, Valores en bits a descomprimir.
* \return cadena ce valores descomprimidos.
*/
unsigned char * transformBitToByte(unsigned char * bitValues)
{
	int i, j, counter = 0;
	unsigned char * byteString = NULL;
	
	byteString = (unsigned char *) calloc (NUM_NUCLEOTIDE_POSITIONS, sizeof(unsigned char));
	if (!byteString) {
		perror("[FormatDb Nucl] Unable to allocate memory (transformBitToByte)\n");
		return NULL;
	}
	
	for(i = 0; i < NUCL_BIT_SIZE; i++) {
		for(j = 0; j < BYTE_SIZE; j++)
		{
			if (((counter % BYTE_SIZE) == 0) && (counter != 0)) {
				byteString[counter] = maskAnd(bitValues[i], j);
			}
			
			counter++;
		}
	}
	
	return byteString;
}


/**
* \fn unsigned char fromGuidToResult(unsigned char * bitValues, long int guid)
* \brief Función que devuelve un valor de la cadena de bits a partir de un GUID.
* \param[in] bitValues, Cadena de bytes comprimidos.
* \param[in] guid, Identificador Hash del valor buscado.
* \return valor buscado.
*/
unsigned char fromGuidToResult(unsigned char * bitValues, long int guid)
{
	int block = guid / BYTE_SIZE;
	int position = guid % BYTE_SIZE;
	
	return maskAnd(bitValues[block], position);
}
