# Parallel implementation of Sequence Alignment

Sequence Alignment – a way to estimate a similarity of two strings of letters - is an important field in bioinformatics. Sequence is a string of capital letters including hyphen sign (-). Each letter in the sequence represents DNA, RNA, or protein. Identification of region of similarity of set of Sequences is extremely time consuming. This project deals with a simplified version of Sequence Alignment of two sequences. The purpose of the project is to parallelize the basic algorithm to produce an efficient computation within MPI, OpenMP and CUDA environment.
See full project definition [here](https://github.com/vladuha12/PIOSA/blob/main/Parallel%20implementation%20of%20Sequence%20Alignment%20documentation.pdf).


## Running requirements:
    1) GPU with CUDA API.
    2) "inc" folder with CUDA libraries.
    3) linux operating system.
    4) nvcc CUDA compiler.
    5) mpicc MPI compiler.

## Configurating project properties:
    1) Open "defines.h".
    2) Edit whatever you need.

## Input File format:
    First line - weight coefficients W1, W2, W3, W4.
    Next line – Seq1 (not more than 3000 chars in line).
    Next line – the number NS2 of Sequences Seq2 to check against Seq1.
    Next NS2 lines - Seq2 in each line (not more than 2000 chars in each line).

    Example:
	    2  1.5  1.1  1.3
	    MNMLWVVSGQTYQQLPVDFKTFRQATVGNTQHQTFTFSYPFE
	    3
	    GNTQHQTFTFSYPFE
	    EGQATVGNT
	    MNMLWVVSGQTY

## Output File format:
    This file contains NS2 lines with n and k found for each Sequence Seq2 from the input file, in order
    Sequences Seq2 appear in the input file, where n is an Offset and k is a hyphen location of the Mutant
    Sequence MS(k) with the best Alignment Score.
	
    Example for input from above:
		n = 26 k = 1
		n = 21 k = 9
		n = 0 k = 12

## Building the project:
	1) Open terminal in current folder.
	2) Check that cuda version in makefile is right.
	3) Run "make".

## Cleaning project folder:
    1) Open terminal in current folder.
    2) Run "make clean" ("make clean" will delete all the .o, the main executable and output.txt files).

## Running the project:
	1) Put inside current folder file with data named "input.txt"
	2) Open terminal in current folder.
	3) For one thread run "make run1".
	4) For two threads run "make run2".
	5) For two threads on cluster run "make runOn2" after you create "mf" file with the computers ip's.

## Author
[VC](https://github.com/vladuha12)
Vlad Chebanov 
