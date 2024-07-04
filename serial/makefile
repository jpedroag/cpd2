CC = gcc 

TARGET = tsp

CFLAGS = -Wall -g -fopenmp -std=c99 -O3

all:
		$(CC) $(CFLAGS) -o $(TARGET) tsp.c queue.c node.c

clear:
		rm $(TARGET)