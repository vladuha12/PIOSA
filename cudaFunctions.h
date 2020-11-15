/*
 ============================================================================
 Name        : cudaFunctions.h
 Author      : Vlad Chebanov
 Version     : FINAL
 Description : Parallel implementation of Sequence Alignment 
 ============================================================================
*/

#ifndef CUDAFUNCTIONS_H_
#define CUDAFUNCTIONS_H_

int countScore(int *arr, int sizeOfArr, float w1, float w2, float w3, float w4, float *maxScore);
int *countChars(char *resultArr, int sizeOfSeq2, int sizeOfMutationSeq);
int computeOnGPU(char *seq1, int sizeOfSeq1 , char **data, int numOfiterations, int sizeOfSeq2, int sizeOfMutatedSeq, float w1, float w2, float w3, float w4, int *n, int *k, float *score);

#endif /*CUDAFUNCTIONS_H_ */
