/*
 ============================================================================
 Name        : SequencesComparator.h
 Author      : Vlad Chebanov
 Version     : FINAL
 Description : Parallel implementation of Sequence Alignment 
 ============================================================================
*/

#ifndef SEQUENCESCOMPARATOR_H_
#define SEQUENCESCOMPARATOR_H_

int scoreCalculator(char *seq1, char *seq2, float w1, float w2, float w3,
	float w4, int *n, int *k);
char *mutation(char *seq, int i, int size);

#endif /*SEQUENCESCOMPARATOR_H_ */
