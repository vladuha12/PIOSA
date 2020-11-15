/*
 ============================================================================
 Name        : DNA.c
 Author      : Vlad Chebanov
 Version     : FINAL
 Description : Parallel implementation of Sequence Alignment 
 ============================================================================
*/

#include "defines.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "IO.h"
#include "SequencesComparator.h"

int main(int argc, char *argv[])
{
	int my_rank; 				/* rank of process */
	int p; 					/* number of processes */
	int master = 0; 			/* rank of master */
	int slave = 1; 				/* rank of slave */
	int tag = 0; 				/* tag for messages */
	int i, j;				/* index */
	MPI_Status status; 			/* return status for receive */

	char *seq1, *buffer; 					/* pointer to sequence 1 from file, buffer to send to via MPI */
	char **seq2; 						/* pointer to 2d array of sequence 2 from file */
	int ns2; 						/* number of sequences */
	int *bufferLength = (int*) calloc(1, sizeof(int));	/* buffer length to send via MPI */
	CHECK_ALLOC(bufferLength, NULL, EXIT_FAILURE);	
	float w1, w2, w3, w4; 					/* weights from file */
	int *finalResultPC1;					/* final result from PC 1 */
	int *finalResultPC2;					/* final result from PC 2 */
	int *resultFromPC2;					/* final result from PC 2 recieved from MPI */
	int *n = (int*) calloc(1, sizeof(int));			/* final offset location */
	CHECK_ALLOC(n, NULL, EXIT_FAILURE);
	int *k = (int*) calloc(1, sizeof(int));			/* final mutation location */
	CHECK_ALLOC(k, NULL, EXIT_FAILURE);

	
	MPI_Init(&argc, &argv);					/* start up MPI */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);		/* find out process rank */
	MPI_Comm_size(MPI_COMM_WORLD, &p);			/* find out number of processes */
	
	/* checking that no more than 2 processes are run */
	if (p > 2)
	{
		printf("Please run on 1 or 2 pc's\n");
		return EXIT_FAILURE;
	}

	/* slave portion of calculation */
	if (my_rank != master)
	{	
		/* receiving size of buffer from master */
		MPI_Recv(bufferLength, 1, MPI_INT, master, tag, MPI_COMM_WORLD, &status);
		buffer = (char*) calloc(*bufferLength, sizeof(char));
		CHECK_ALLOC(buffer, NULL, EXIT_FAILURE);
		
		/* receiving buffer from master (the file content) */
		MPI_Recv(buffer, *bufferLength* sizeof(char), MPI_CHAR, master, tag, MPI_COMM_WORLD, &status);
		
		/* separating weights from buffer */
		RUN_AND_CHECK(seperateWeightsFromBuffer(buffer, &w1, &w2, &w3, &w4));
		
		/* separate sequence 1 from buffer */
		seq1 = seperateSeq1FromBuffer(buffer, MAX_LENGTH_SEQ1);
		CHECK_ALLOC(seq1, NULL, EXIT_FAILURE);

		/* seperate number of strings from buffer */
		ns2 = seperateNumOfStringsFromBuffer(buffer);
		CHECK_ALLOC(ns2, 0, EXIT_FAILURE);

		/* separate sequences 2 from buffer */
		seq2 = seperateSeq2FromBuffer(buffer, MAX_LENGTH_SEQ2, ns2);
		CHECK_ALLOC(seq2, NULL, EXIT_FAILURE);
		
		/* allocating array for result from slave */
		finalResultPC2 = (int*) calloc(ns2 * 2, sizeof(int));
		CHECK_ALLOC(finalResultPC2, NULL, EXIT_FAILURE);

		/* calculating second half of strings with sequence compare */
		for (i = ns2 / p; i < ns2; i++)
		{
			RUN_AND_CHECK(scoreCalculator(seq1, seq2[i], w1, w2, w3, w4, n, k));
			finalResultPC2[i * 2] = *n;
			finalResultPC2[i * 2 + 1] = *k;
		}

		/* sending result to master */
		MPI_Send(finalResultPC2, ns2 * 2, MPI_INT, master, tag, MPI_COMM_WORLD);
	}
	/* master portion of calculation */
	else
	{
		printf("Starting calculation with %d pc's\n", p);
		/*reading data from file*/
		buffer = readFile(INPUT_FILE_NAME);
		CHECK_ALLOC(buffer, NULL, EXIT_FAILURE);

		/* separating weights from buffer */
		RUN_AND_CHECK(seperateWeightsFromBuffer(buffer, &w1, &w2, &w3, &w4));

		/* separate sequence 1 from buffer */
		seq1 = seperateSeq1FromBuffer(buffer, MAX_LENGTH_SEQ1);
		CHECK_ALLOC(seq1, NULL, EXIT_FAILURE);

		/* seperate number of strings from buffer */
		ns2 = seperateNumOfStringsFromBuffer(buffer);
		CHECK_ALLOC(ns2, 0, EXIT_FAILURE);
		
		/* separate sequences 2 from buffer */
		seq2 = seperateSeq2FromBuffer(buffer, MAX_LENGTH_SEQ2, ns2);
		CHECK_ALLOC(seq2, NULL, EXIT_FAILURE);

		/* allocating array for result from master */
		finalResultPC1 = (int*) calloc(ns2 * 2, sizeof(int));
		CHECK_ALLOC(finalResultPC1, NULL, EXIT_FAILURE);

		/* if 2 pc's are run */
		if (p > 1)
		{
			/* sending buffer size to slave */
			*bufferLength = strlen(buffer) + 1;
			MPI_Send(bufferLength, 1, MPI_INT, 1, tag,
				MPI_COMM_WORLD);

			/* sending buffer to slave (the file content) */
			MPI_Send(buffer, *bufferLength, MPI_CHAR, 1, tag,
				MPI_COMM_WORLD);
			
			/* calculating first half of strings with sequence compare */
			for (i = 0; i < ns2 / p; i++)
			{
				RUN_AND_CHECK(scoreCalculator(seq1, seq2[i], w1, w2, w3, w4, n, k));
				finalResultPC1[i * 2] = *n;
				finalResultPC1[i * 2 + 1] = *k;
			}

			/* allocating memory for results from master */
			resultFromPC2 = (int*) calloc(ns2 * 2, sizeof(int));
			CHECK_ALLOC(resultFromPC2, NULL, EXIT_FAILURE);

			/* receiving results from slave */
			MPI_Recv(resultFromPC2, ns2 * 2, MPI_INT, slave, tag, MPI_COMM_WORLD, &status);

			/* writing results from master to file */
			for (i = 0; i < ns2 / p; i++)
			{
				RUN_AND_CHECK(writeToFile(OUTPUT_FILE_NAME, finalResultPC1[i * 2],finalResultPC1[i * 2 + 1]));
			}
			
			/* writing results from slave to file */
			for (i = ns2 / p; i < ns2; i++)
			{
				RUN_AND_CHECK(writeToFile(OUTPUT_FILE_NAME, resultFromPC2[i * 2],resultFromPC2[i * 2 + 1]));
			}

			/* freeing results array */
			free(resultFromPC2);
		}

		/* if one pc is run */
		else
		{
			/* calculating all the strings with sequence compare */
			for (i = 0; i < ns2; i++)
			{
				RUN_AND_CHECK(scoreCalculator(seq1, seq2[i], w1, w2, w3, w4, n, k));
				finalResultPC1[i * 2] = *n;
				finalResultPC1[i * 2 + 1] = *k;
			}

			/* writnig results to output file */
			for (i = 0; i < ns2; i++)
			{
				RUN_AND_CHECK(writeToFile(OUTPUT_FILE_NAME, finalResultPC1[i * 2],finalResultPC1[i * 2 + 1]));
			}
		}
		
		/* freeing results */
		free(finalResultPC1);
	}
	
	/* starting freeing */
	free(buffer);
	free(seq1);
	FREE_2D_ARRAY(seq2,ns2);

	/*shutting down MPI */
	MPI_Finalize();

	return EXIT_SUCCESS;
}
