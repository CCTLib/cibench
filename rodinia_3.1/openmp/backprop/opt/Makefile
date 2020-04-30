# C compiler
CC = gcc
CC_FLAGS_1 = -g -fopenmp -O3 -Ofast -march=native -mtune=native -flto -fprofile-generate

CC_FLAGS_2 = -g -fopenmp -O3 -Ofast -march=native -mtune=native -flto -fprofile-use -fprofile-correction


all: first second 

first: backprop.o1 facetrain.o1 imagenet.o1 backprop_kernel.o1 
	$(CC) $(CC_FLAGS_1) backprop.o facetrain.o imagenet.o backprop_kernel.o -o backprop -lm
	./backprop 6553600
	rm -rf backprop

second: backprop.o2 facetrain.o2 imagenet.o2 backprop_kernel.o2 
	$(CC) $(CC_FLAGS_2) backprop.o facetrain.o imagenet.o backprop_kernel.o -o backprop -lm

%.o: %.[ch]a
	$(CC) $(CC_FLAGS) $< -c

facetrain.o1: facetrain.c backprop.h
	$(CC) $(CC_FLAGS_1) facetrain.c -c
	
backprop.o1: backprop.c backprop.h
	$(CC) $(CC_FLAGS_1) backprop.c -c

backprop_kernel.o1: backprop_kernel.c backprop.h
	$(CC) $(CC_FLAGS_1) backprop_kernel.c -c

imagenet.o1: imagenet.c backprop.h
	$(CC) $(CC_FLAGS_1) imagenet.c -c


facetrain.o2: facetrain.c backprop.h
	$(CC) $(CC_FLAGS_2) facetrain.c -c
	
backprop.o2: backprop.c backprop.h
	$(CC) $(CC_FLAGS_2) backprop.c -c

backprop_kernel.o2: backprop_kernel.c backprop.h
	$(CC) $(CC_FLAGS_2) backprop_kernel.c -c

imagenet.o2: imagenet.c backprop.h
	$(CC) $(CC_FLAGS_2) imagenet.c -c

clean:
	rm -f *.o *~ *.gcda backprop backprop_cuda.linkinfo
