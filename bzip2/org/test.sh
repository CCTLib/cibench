#!/bin/bash

for i in {1..15}

do
	./bzip2 -k input.source
	rm input.source.bz2
done
