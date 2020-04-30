include ../Makefile.inc

PGO_OPT1 = -g -O3 -Ofast -march=native -mtune=native -flto -fprofile-generate
PGO_OPT2 = -g -O3 -Ofast -march=native -mtune=native -flto -fprofile-use -fprofile-correction

all: first second

first: larson1 larson-hoard1
	./larson 10 7 8 10000 1000 1 4 4
	LD_PRELOAD=../../src/libhoard.so ../../benchmarks/larson/larson 10 7 8 10000 1000 1 4 4

second: larson2 larson-hoard2


larson1: larson.cpp
	$(CXX) $(CXXFLAGS) $(PGO_OPT1) larson.cpp -o larson -lpthread
larson2: larson.cpp
	$(CXX) $(CXXFLAGS) $(PGO_OPT2) larson.cpp -o larson -lpthread

larson-hoard1: larson.cpp
	$(CXX) $(CXXFLAGS) $(PGO_OPT1) larson.cpp -L../../src -lhoard -o larson-hoard -lpthread
larson-hoard2: larson.cpp
	$(CXX) $(CXXFLAGS) $(PGO_OPT2) larson.cpp -L../../src -lhoard -o larson-hoard -lpthread

clean:
	rm -f *.gcda larson larson-hoard
