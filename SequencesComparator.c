/*
 ============================================================================
 Name        : SequencesComparator.c
 Author      : Vlad Chebanov
 Version     : FINAL
 Description : Parallel implementation of Sequence Alignment 
 ============================================================================
*/

#include "SequencesComparator.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include "cudaFunctions.h"
#include <mpi.h>
#include "defines.h"

/* making mutated sequence */
char *mutation(char *seq, int f, int size)
{
	/* allocating space for mutated sequence */
	char *newStr = (char*) calloc(size + 1, sizeof(char));
	CHECK_ALLOC(newStr, NULL, NULL);
	char ch = '-';
	int carry = 0;

	/* copying chars from original sequence to mutated */
	for (int i = 0; i < size + 1; i++)
	{
		if (i == f)
		{
			newStr[i] = ch;
			carry = 1;
		}
		else
			newStr[i] = seq[i - carry];
	}
	return newStr;
}

/* calculate best allignment score */
int scoreCalculator(char *seq1, char *seq2, float w1, float w2, float w3,
	float w4, int *n, int *k)
{
	int lenSeq1 = strlen(seq1);			/* calculating sequence 1 size */
	int lenSeq2 = strlen(seq2);			/* calculating sequence 2 size */
	int carry,i;
	float tempScore, score, tempScore2;		/* allocating temp scores */
	int offset, tempN, tempK, tempN2, tempK2;	/* allocating offset and temp results */
	char *mutatedSeq;				/* allocating pointer to mutated sequence */
	int numOfiterations = lenSeq1 - lenSeq2;	/* calculating maximum number of iterations */
	
	char **mutatedSeqDB = (char **) calloc(lenSeq2, sizeof(char*));		/* allocating muteted sequence 2d array to be sent to CUDA */
	CHECK_ALLOC(mutatedSeqDB, NULL, EXIT_FAILURE);
	for (i = 0; i < lenSeq2; i++)
	{
		mutatedSeqDB[i] = (char*) calloc(lenSeq2 + 1, sizeof(char));
		CHECK_ALLOC(mutatedSeqDB[i], NULL, EXIT_FAILURE);
	}

	/* runnig mutation creation in parallel with OpenMP and adding them to db*/
	omp_set_num_threads(4);
	#pragma omp parallel for
	for (i = 1; i <= lenSeq2; i++)
	{
		mutatedSeq = mutation(seq2, i, lenSeq1);
		memcpy(mutatedSeqDB[i - 1], mutatedSeq, lenSeq2 + 1);
	}
	
	/* initializing clock */
	clock_t t;
	t = clock();

	/* in case of sequence is too large to upload to cuda, separating string to 3 portions */
	if (lenSeq2 > MAX_CUDA_SEQUENCE_LENGTH)
	{
		/* calculating first third of mutations */
		carry = lenSeq2 % 3;
		RUN_AND_CHECK(computeOnGPU(seq1, lenSeq1, mutatedSeqDB, numOfiterations, lenSeq2 / 3, lenSeq2 + 1, w1, w2, w3, w4, n, k, &score));
		printf("Sequence too large for CUDA. Separating string\n");

		/* saving results */
		tempN = *n;
		tempK = *k;
		tempScore = score;

		/* calculating second third of mutations */
		RUN_AND_CHECK(computeOnGPU(seq1, lenSeq1, mutatedSeqDB + lenSeq2 / 3, numOfiterations, lenSeq2 / 3, lenSeq2 + 1, w1, w2, w3, w4, n, k, &score));

		/* saving results */
		tempN2 = *n;
		tempK2 = *k;
		tempScore2 = score;

		/* calculating last third of mutations */
		RUN_AND_CHECK(computeOnGPU(seq1, lenSeq1, mutatedSeqDB + lenSeq2 / 3 * 2, numOfiterations, lenSeq2 / 3 + carry, lenSeq2 + 1, w1, w2, w3, w4, n, k, &score));
		
		/* choosing maximum result from 3 */
		if (tempScore > score)
		{ *n = tempN;
			*k = tempK;
			score = tempScore;
		}
		else
		{
			if (tempScore2 > score)
			{ 	*n = tempN2;
				*k = tempK2 + lenSeq2 / 3 - carry;
				score = tempScore2;
			}
			else
			{ 	*k = *k + lenSeq2 / 3 * 2;
			}
		}
	}

	/* running regular calculation on cuda if sequence isn't too large for cuda */
	else
	{
		RUN_AND_CHECK(computeOnGPU(seq1, lenSeq1, mutatedSeqDB, numOfiterations, lenSeq2, lenSeq2 + 1, w1, w2, w3, w4, n, k, &score));
	}
	
	/* stopping clock */
	t = clock() - t;
	double time_taken = ((double) t) / CLOCKS_PER_SEC;	// in seconds
	printf("CUDA took %f seconds to execute\n", time_taken);

	/* freeing arrays */
	free(mutatedSeq);
	FREE_2D_ARRAY(mutatedSeqDB, lenSeq2);
	return EXIT_SUCCESS;
}
