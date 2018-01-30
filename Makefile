CFLAGS=-g -O3

all: gpe

gpe: main.o complex.o aux.o matriz.o
	g++ $(CFLAGS) -o gpe main.o complex.o aux.o matriz.o -fopenmp

main.o: main.cpp
	g++ $(CFLAGS) -o main.o -c main.cpp -fopenmp

complex.o: complex.cpp complex.h
	g++ $(CFLAGS) -o complex.o -c complex.cpp -fopenmp

aux.o: aux.cpp aux.h
	g++ $(CFLAGS) -o aux.o -c aux.cpp -fopenmp

matriz.o: matriz.cpp matriz.h
	g++ $(CFLAGS) -o matriz.o -c matriz.cpp -fopenmp

clean:
	rm *.o gpe dados/*
