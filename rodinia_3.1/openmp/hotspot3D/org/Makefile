CC = gcc -g -O3 -Ofast -march=native -mtune=native -fopenmp
CCFLAGS1 = -flto -fprofile-generate
CCFLAGS2 = -flto -fprofile-use -fprofile-correction
OUTPUT = *.txt

all: first second

first:
	$(CC) -o 3D $(CCFLAGS1) 3D.c -lm
	OMP_NUM_THREADS=1 ./3D 512 8 1000 ../../../data/hotspot3D/power_512x8 ../../../data/hotspot3D/temp_512x8 output.txt 
	rm -rf 3D
second:
	$(CC) -o 3D $(CCFLAGS2) 3D.c -lm 


clean:
	rm -f 3D.gcda 3D $(OUTPUT)
