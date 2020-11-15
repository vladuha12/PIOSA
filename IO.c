/*
 ============================================================================
 Name        : IO.c
 Author      : Vlad Chebanov
 Version     : FINAL
 Description : Parallel implementation of Sequence Alignment 
 ============================================================================
*/

#include "IO.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defines.h"

int writeToFile(char *fileName, int n, int k)
{
	/*declare a file pointer */
	FILE * infile;

	/*open an existing file for writing */
	infile = fopen(fileName, "a");

	/*quit if the file does not exist */
	CHECK_ALLOC(infile, NULL, EXIT_FAILURE);

	/* writing results to output file */
	fprintf(infile, "n = %d k = %d\n", n, k);
	fclose(infile);
	return EXIT_SUCCESS;
}

char *readFile(char *fileName)
{
	/*declare a file pointer */
	FILE * infile;
	char *buffer;
	long numbytes;

	/*open an existing file for reading */
	infile = fopen(fileName, "r");

	/*quit if the file does not exist */
	CHECK_ALLOC(infile, NULL, NULL);

	/*Get the number of bytes */
	fseek(infile, 0, SEEK_END);
	numbytes = ftell(infile);

	/*reset the file position indicator to
	 the beginning of the file */
	fseek(infile, 0, SEEK_SET);

	/*grab sufficient memory for the
	 buffer to hold the text */
	buffer = (char*) calloc(numbytes, sizeof(char));

	/*memory error */
	CHECK_ALLOC(buffer, NULL, NULL);

	/*copy all the text into the buffer */
	fread(buffer, sizeof(char), numbytes, infile);
	fclose(infile);
	return buffer;
}

/* separate sequences 2 from buffer */
char **seperateSeq2FromBuffer(char *buffer, int maxLength, int ns2)
{
	int i;
	/* allocating 2d array for sequences 2 */
	char **seq2 = (char **) calloc(ns2, sizeof(char*));
	CHECK_ALLOC(seq2, NULL, NULL);
	
	for (i = 0; i < ns2; i++)
	{
		seq2[i] = (char*) calloc(maxLength + 1, sizeof(char));
		CHECK_ALLOC(seq2[i], NULL, NULL);
	}
	/* allocating temp array for buffer because strtok is used */
	char *tempBuffer = (char*) calloc(strlen(buffer) + 1, sizeof(char));
	CHECK_ALLOC(tempBuffer, NULL, NULL);
	
	/* copying buffer to temp array */
	strcpy(tempBuffer, buffer);
	char delim[] = "\n";
	char *ptr = strtok(tempBuffer, delim);
	for (i = 0; i < 3; i++)
		ptr = strtok(NULL, delim);

	/* separating string from temp buffer */
	for (i = 0; i < ns2; i++)
	{
		/* checking that one of the lines is not over the limit */
		if (strlen(ptr) > maxLength){
			printf("Length of sequence 2 is over the size limit. Order of string = %d\n",i + 1);
			return NULL;
		}
		strcpy(seq2[i], ptr);
		ptr = strtok(NULL, delim);
	}
	
	/* freeing temp buffer */
	free(tempBuffer);
	return seq2;
}

/* seperate number of strings from buffer */
int seperateNumOfStringsFromBuffer(char *buffer)
{
	/* allocating temp array for buffer because strtok is used */
	char *tempBuffer = (char*) calloc(strlen(buffer) + 1, sizeof(char));
	CHECK_ALLOC(tempBuffer, NULL, EXIT_FAILURE);
	
	/* copying buffer to temp array */
	strcpy(tempBuffer, buffer);
	char delim[] = "\n";
	/* separating number of string from temp buffer */
	char *ptr = strtok(tempBuffer, delim);
	for (int i = 0; i < 2; i++)
		ptr = strtok(NULL, delim);
	
	/* freeing temp buffer */
	free(tempBuffer);
	return atoi(ptr);
}


/* separate sequence 1 from buffer */
char *seperateSeq1FromBuffer(char *buffer, int maxLength)
{
	/* allocating temp array for buffer because strtok is used */
	char *tempBuffer = (char*) calloc(strlen(buffer) + 1, sizeof(char));
	CHECK_ALLOC(tempBuffer, NULL, NULL);

	/* copying buffer to temp array */
	strcpy(tempBuffer, buffer);
	
	/* allocating array of string 1 to be reterned */
	char *seq1 = (char*) calloc(maxLength + 1, sizeof(char));
	CHECK_ALLOC(seq1, NULL, NULL);

	char delim[] = "\n";
	/* separating sequence 1 from temp buffer */
	char *ptr = strtok(tempBuffer, delim);
	ptr = strtok(NULL, delim);
	if (strlen(ptr) > maxLength){
		printf("Length of sequence 1 is over the size limit\n");
		return NULL;
	}
	/* copying sequence 1 */	
	strcpy(seq1, ptr);

	/* freeing temp buffer */
	free(tempBuffer);
	return seq1;

}

/* separate weights from buffer */
int seperateWeightsFromBuffer(char *buffer, float *w1, float *w2, float *w3,
	float *w4)
{
	/* allocating temp array for buffer because strtok is used */
	char *tempBuffer = (char*) calloc(strlen(buffer), sizeof(char));
	CHECK_ALLOC(tempBuffer, NULL, EXIT_FAILURE);

	/* copying buffer to temp array */
	strcpy(tempBuffer, buffer);
	char delim[] = " ";
	char *ptr = strtok(tempBuffer, delim);

	/* separating weights from temp buffer */
	for (int i = 0; i < 4; i++)
	{
		switch (i)
		{
			case 0:
				*w1 = atof(ptr);
				break;
			case 1:
				*w2 = atof(ptr);
				break;
			case 2:
				*w3 = atof(ptr);
				break;
			case 3:
				*w4 = atof(ptr);
				break;
		}
		ptr = strtok(NULL, delim);
	}
	/* freeing temp buffer */
	free(tempBuffer);
	return EXIT_SUCCESS;
}
