/*
 ============================================================================
 Name        : defines.h
 Author      : Vlad Chebanov
 Version     : FINAL
 Description : Parallel implementation of Sequence Alignment 
 ============================================================================
*/

#ifndef DEFINES_H_
#define DEFINES_H_

#define INPUT_FILE_NAME (char*)"input.txt"
#define OUTPUT_FILE_NAME (char*)"output.txt"
#define MAX_LENGTH_SEQ1 3000
#define MAX_LENGTH_SEQ2 2000
#define MAX_CUDA_SEQUENCE_LENGTH 1000
#define MAX_CHAR_LENGTH_CONSERVATIVE 7
#define CUDA_THREADS_PER_BLOCK 1024
#define WEIGHTS_NUM 4
#define CONSERVATIVEGROUPCUDA "NDEQ", "NEQK", "STA", "MILV", "UHRK", "MHQK", "FYW", "HY", "MILF"
#define SEMICONSERVATIVEGROUPCUDA "SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM"


#define CHECK_ALLOC(A,B,C) if(A==B){printf("Error allocating memory at \"%s\" function\n", __func__);return C;}
#define CHECK_CUDA_ALLOC(A) if (err != cudaSuccess){fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));exit(EXIT_FAILURE);}
#define CHECK_CUDA_MEMCPY(A) if (err != cudaSuccess){fprintf(stderr, "Failed to copy to device memory - %s\n", cudaGetErrorString(err));exit(EXIT_FAILURE);}
#define CHECK_CUDA_KERNEL(A) if (err != cudaSuccess){fprintf(stderr, "Failed to launch cudaSequenceCompare kernel - %s\n", cudaGetErrorString(err));exit(EXIT_FAILURE);}
#define CUDA_FREE(A) if(cudaFree(A) != cudaSuccess){fprintf(stderr, "Failed to free device data - %s\n",cudaGetErrorString(err));exit(EXIT_FAILURE);}
#define FREE_2D_ARRAY(A,B) for(int i = 0; i < B; i++){free(A[i]);}free(A)
#define RUN_AND_CHECK(A) if(A){return EXIT_FAILURE;}

#endif /*DEFINES_H_ */
