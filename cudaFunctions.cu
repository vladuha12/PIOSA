/*
 ============================================================================
 Name        : cudaFunctions.cu
 Author      : Vlad Chebanov
 Version     : FINAL
 Description : Parallel implementation of Sequence Alignment 
 ============================================================================
*/

#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "cudaFunctions.h"
#include <cuda.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "defines.h"

/* defining conservative groups for GPU usage */
__device__ char conservativeGroupCUDA[][MAX_CHAR_LENGTH_CONSERVATIVE] = { CONSERVATIVEGROUPCUDA};
__device__ char semiConservativeGroupCUDA[][MAX_CHAR_LENGTH_CONSERVATIVE] = { SEMICONSERVATIVEGROUPCUDA};

/* counting scores for all mutations */
int countScore(int *arr, int sizeOfArr, float w1, float w2, float w3, float w4,
	float *maxScore)
{
	float max = -99999999.9;
	float tempScore;
	int k = 0, i;
	for (i = 0; i < sizeOfArr * WEIGHTS_NUM; i += WEIGHTS_NUM)
	{
		tempScore = arr[i] *w1 - arr[i + 1] *w2 - arr[i + 2] *w3 - arr[i + 3] *w4;
		if (max < tempScore)
		{
			max = tempScore;
			k = i / WEIGHTS_NUM + 1;
		}
	}
	*maxScore = max;
	return k;
}

/* counting signs of all mutations */
int *countChars(char *resultArr, int sizeOfSeq2, int sizeOfMutationSeq)
{
	int *result = (int*) calloc(sizeOfSeq2 * WEIGHTS_NUM, sizeof(int));
	CHECK_ALLOC(result, NULL, NULL);
	int i;
	for (i = 0; i < sizeOfSeq2 * sizeOfMutationSeq; i++)
	{
		switch (resultArr[i])
		{
			case '*':
				result[(i / sizeOfMutationSeq) * WEIGHTS_NUM + 0]++;
				break;
			case ':':
				result[(i / sizeOfMutationSeq) * WEIGHTS_NUM + 1]++;
				break;
			case '.':
				result[(i / sizeOfMutationSeq) * WEIGHTS_NUM + 2]++;
				break;
			case ' ':
				if (i != sizeOfSeq2 *sizeOfMutationSeq - 1)
					result[(i / sizeOfMutationSeq) * WEIGHTS_NUM + 3]++;
				break;
		}
	}
	return result;
}

/* finding similarities in semi conservative group */
__device__ int semiConservativeGroupComareCUDA(char ch1, char ch2)
{
	int i, j, flagCh1, flagCh2, conservativeGroupSize =
		sizeof(semiConservativeGroupCUDA) / (MAX_CHAR_LENGTH_CONSERVATIVE* sizeof(char));
	for (i = 0; i < conservativeGroupSize; i++)
	{
		flagCh1 = 0;
		flagCh2 = 0;
		for (j = 0; j < MAX_CHAR_LENGTH_CONSERVATIVE; j++)
		{
			if (!flagCh1 && semiConservativeGroupCUDA[i][j] == ch1)
				flagCh1 = 1;
			else
			{
				if (!flagCh2 && semiConservativeGroupCUDA[i][j] == ch2)
					flagCh2 = 1;
			}
			if (flagCh1 && flagCh2)
				return 1;
		}
	}
	return 0;
}

/* finding similarities in conservative group */
__device__ int conservativeGroupComareCUDA(char ch1, char ch2)
{
	int i, j, flagCh1, flagCh2, conservativeGroupSize =
		sizeof(conservativeGroupCUDA) / (MAX_CHAR_LENGTH_CONSERVATIVE* sizeof(char));
	for (i = 0; i < conservativeGroupSize; i++)
	{
		flagCh1 = 0;
		flagCh2 = 0;
		for (j = 0; j < MAX_CHAR_LENGTH_CONSERVATIVE; j++)
		{
			if (!flagCh1 && conservativeGroupCUDA[i][j] == ch1)
				flagCh1 = 1;
			else
			{
				if (!flagCh2 && conservativeGroupCUDA[i][j] == ch2)
					flagCh2 = 1;
			}
			if (flagCh1 && flagCh2)
				return 1;
		}
	}
	return 0;
}

/* main cuda kernel function for chars comperation running on all mutation at once */
__global__ void cudaSequenceCompare(char *seq1, char *data, char *result,
	int numElements, int size, int offset)
{
	int i = blockDim.x *blockIdx.x + threadIdx.x;		/* signing ID to every core */

	/* runnig all comparetions in parallel */
	if (i < numElements)
	{
		if (seq1[(i % size) + offset] == data[i])
			result[i] = '*';
		else if (conservativeGroupComareCUDA(seq1[(i % size) + offset],
				data[i]))
		{
			result[i] = ':';
		}
		else if (semiConservativeGroupComareCUDA(seq1[(i % size) + offset],
				data[i]))
		{
			result[i] = '.';
		}
		else
			result[i] = ' ';
	}
}

/* cuda kernel initializer */
int computeOnGPU(char *seq1, int sizeOfSeq1 , char **data, int numOfiterations, int sizeOfSeq2, int sizeOfMutatedSeq, float w1, float w2, float w3, float w4, int *n, int *k, float *score)
{
	int i, kTemp = 0, kMax = 0, nMax = 0;		/* temp result */
	float maxTemp, maxScore;			/* score results */
	int *result;					/* pinter to an array of signs results */

	cudaError_t err = cudaSuccess;			/* cuda error flag */

	/*Host allocations */
	/* initializing 1D array of mutations */
	char *h_Array = (char*) calloc(sizeOfSeq2 *sizeOfMutatedSeq + 1, sizeof(char));
	CHECK_ALLOC(h_Array, NULL, EXIT_FAILURE);

	/* appending 2D array to 1D */
	for (i = 0; i < sizeOfSeq2; i++)
	{
		strncat(h_Array, data[i], sizeOfMutatedSeq);
	}

	/*Allocate memory on GPU to copy the data from the host */
	/* allocating sequence one space on GPU */
	char *d_Seq1;
	err = cudaMalloc((void **) &d_Seq1, sizeof(char) * sizeOfSeq1 + 1);
	CHECK_CUDA_ALLOC(err);

	/* allocating results array space on GPU */
	char *d_Result;
	err = cudaMalloc((void **) &d_Result, sizeof(char) *sizeOfSeq2 *sizeOfMutatedSeq + 1);
	CHECK_CUDA_ALLOC(err);

	/* allocating array of mutations space on GPU */
	char *d_Array;
	err = cudaMalloc((void **) &d_Array, sizeof(char) *sizeOfSeq2 *sizeOfMutatedSeq + 1);
	CHECK_CUDA_ALLOC(err);

	/* copy host strings to device memory */
	/* copying array of mutations from host to GPU memory */
	err = cudaMemcpy(d_Array, h_Array, sizeof(char) *sizeOfSeq2 *sizeOfMutatedSeq + 1, cudaMemcpyHostToDevice);
	CHECK_CUDA_MEMCPY(err);

	/* copying sequence 1 from host to GPU memory */
	err = cudaMemcpy(d_Seq1, seq1, sizeof(char) * sizeOfSeq1 + 1, cudaMemcpyHostToDevice);
	CHECK_CUDA_MEMCPY(err);

	/* initializing cuda */
	int threadsPerBlock = CUDA_THREADS_PER_BLOCK;
	int blocksPerGrid = sizeOfSeq2 * sizeOfSeq2;
		
	/* runnig loop maximum offset iterations possible */
	for (i = 0; i < numOfiterations; i++)
	{
		/* launching GPU kernel functions */
		cudaSequenceCompare <<<blocksPerGrid, threadsPerBlock>>> (d_Seq1, d_Array, d_Result, sizeOfSeq2 *sizeOfMutatedSeq, sizeOfMutatedSeq, i);
		
		/* checking for errors */
		err = cudaGetLastError();
		CHECK_CUDA_KERNEL(err);

		/* copying results from GPU to host */
		err = cudaMemcpy(h_Array, d_Result, sizeof(char) *sizeOfSeq2 *sizeOfMutatedSeq + 1, cudaMemcpyDeviceToHost);
		CHECK_CUDA_MEMCPY(err);

		/* counting signs */
		result = countChars(h_Array, sizeOfSeq2, sizeOfMutatedSeq);

		/* counting score */
		kTemp = countScore(result, sizeOfSeq2, w1, w2, w3, w4, &maxTemp);
	
		/* finding maximum score */
		if (maxScore < maxTemp)
		{
			kMax = kTemp;
			maxScore = maxTemp;
			nMax = i;
		}
		/* freeing result array */
		free(result);
	}

	/* rolling back maximum results */
	*n = nMax;
	*k = kMax;
	*score = maxScore;

	// Free allocated memory on GPU
	CUDA_FREE(d_Array);
	CUDA_FREE(d_Result);
	CUDA_FREE(d_Seq1);

	// Free allocated memory on host
	free(h_Array);

	return 0;
}
