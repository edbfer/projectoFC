INCLUDES=/usr/local/cuda-9.1/include
LIBINC=/usr/local/cuda-9.1/lib64
CFLAGS= -I$(INCLUDES) -g -O3
LFLAGS= -L$(LIBINC)
CUDA=

all: omp

omp: CUDA = 0
omp: gpe_omp

gpe_omp: main.o complex_omp.o aux.o matriz.o
	g++ $(CFLAGS) main.o complex.o aux.o matriz.o -o gpe -fopenmp

main.o: main.cpp
	g++ $(CFLAGS) -o main.o -c main.cpp -fopenmp

complex.o: complex.cpp complex.h
	nvcc $(CUDAFLAGS) -x cu -dc complex.cpp -o complex.o

complex_omp.o: complex_omp.cpp complex.h
	g++ $(CFLAGS) -o complex.o -c complex_omp.cpp -fopenmp

aux.o: aux.cpp aux.h
	g++ $(CFLAGS) -o aux.o -c aux.cpp -fopenmp

matriz.o: matriz.cpp matriz.h
	g++ $(CFLAGS) -o matriz.o -c matriz.cpp -fopenmp

clean:
	rm *.o gpe dados/*
