# C compiler
CC = g++
CC_FLAGS_1 = -g -fopenmp -O3 -flto -Ofast -march=native -mtune=native -fprofile-generate
CC_FLAGS_2 = -g -fopenmp -O3 -flto -Ofast -march=native -mtune=native -fprofile-use -fprofile-correction
export LIBRARY_PATH=/usr/lib/gcc/x86_64-redhat-linux/4.8.2/

all: first second

first: 
	$(CC) $(CC_FLAGS_1) srad.cpp -o srad
	./srad 6400 6400 0 3100 0 3100 1 0.5 10
	# For longer execution time:
	#./srad 6400 6400 0 3100 0 3100 1 0.5 50
	rm -f srad
 
second: 
	$(CC) $(CC_FLAGS_2) srad.cpp -o srad



clean:
	rm -f srad
	rm -rf srad.gcda
