#!/bin/sh
cd ../benchmarks
LD_PRELOAD=../src/libhoard.so ../benchmarks/larson/larson 10 7 8 10000 1000 1 1 1
#LD_PRELOAD=../../Hoard-gcc6-op/src/libhoard.so ../benchmarks/larson/larson 10 7 8 10000 1000 1 1 1
