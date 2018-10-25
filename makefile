CC=gcc
LIBS= -lm

all: fmm.c point_heap.c update.c quadratic.c
	$(CC) -o fmm fmm.c point_heap.c update.c quadratic.c $(LIBS)
