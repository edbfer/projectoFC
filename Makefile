INCLUDES=/usr/local/cuda-9.1/include
LIBINC=/usr/local/cuda-9.1/lib64
CFLAGS= -I$(INCLUDES) -g -O3
CUDAFLAGS= $(CFLAGS) -G -gencode arch=compute_50,code=sm_50 -gencode arch=compute_35,code=sm_35 -ccbin /usr/bin/g++-5 -Xcompiler -fopenmp
LFLAGS= -L$(LIBINC)

all: gpe

gpe: main.o complex.o aux.o matriz.o cuda.o gpu.o
	g++ $(CFLAGS) $(LFLAGS) main.o complex.o aux.o matriz.o cuda.o gpu.o -o gpe -fopenmp -lcudart

main.o: main.cpp
	g++ $(CFLAGS) -o main.o -c main.cpp -fopenmp

complex.o: complex.cpp complex.h
	nvcc $(CUDAFLAGS) -x cu -dc complex.cpp -o complex.o

aux.o: aux.cpp aux.h
	g++ $(CFLAGS) -o aux.o -c aux.cpp -fopenmp

matriz.o: matriz.cpp matriz.h
	g++ $(CFLAGS) -o matriz.o -c matriz.cpp -fopenmp

cuda.o: cuda.cu cuda.cuh
	nvcc $(CUDAFLAGS) -x cu -dc cuda.cu -o cuda.o

gpu.o: cuda.o complex.o
	nvcc $(CUDAFLAGS) -dlink cuda.o complex.o -o gpu.o

clean:
	rm *.o gpe dados/*
