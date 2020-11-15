/*
 ============================================================================
 Name        : IO.h
 Author      : Vlad Chebanov
 Version     : FINAL
 Description : Parallel implementation of Sequence Alignment 
 ============================================================================
*/

#ifndef IO_H_
#define IO_H_

int writeToFile(char *fileName, int n, int k);
char *readFile(char *fileName);
int seperateNumOfStringsFromBuffer(char *buffer);
char **seperateSeq2FromBuffer(char *buffer, int maxLength, int ns2);
int seperateWeightsFromBuffer(char *buffer, float *w1, float *w2, float *w3,
	float *w4);
char *seperateSeq1FromBuffer(char *buffer, int maxLength);

#endif /*IO_H_ */