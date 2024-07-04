CC = gcc 

TARGET = tsp

CFLAGS =  -g -fopenmp -std=c99 -O3

all:
		$(CC) $(CFLAGS) -o $(TARGET) tsp-omp.c queue.c node.c

clear:
		rm $(TARGET)