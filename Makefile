build:
	mpicxx -fopenmp -c DNA.c -o DNA.o
	mpicxx -fopenmp -c SequencesComparator.c -o SequencesComparator.o
	mpicxx -fopenmp -c IO.c -o IO.o
	nvcc -I./inc -c cudaFunctions.cu -o cudaFunctions.o
	mpicxx -fopenmp -o DNA DNA.o SequencesComparator.o IO.o cudaFunctions.o /usr/local/cuda-9.1/lib64/libcudart_static.a -ldl -lrt

clean:
	rm -f *.o ./DNA output.txt

run1:
	mpiexec -np 1 ./DNA

run2:
	mpiexec -np 2 ./DNA

runOn2:
	mpiexec -np 2 -machinefile  mf  -map-by  node  ./DNA

