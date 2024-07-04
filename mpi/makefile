CC = mpicc

TARGET = tsp

CFLAGS =  -g -std=c99 -O3

all:
		$(CC) $(CFLAGS) -o $(TARGET) tsp-mpi.c queue.c node.c

clear:
		rm $(TARGET)