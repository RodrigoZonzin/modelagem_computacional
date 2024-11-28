#!/bin/bash
for i in 1 2 3 4 5
do
	python3 gillespie-v1.py $i
	echo "Hello $i"
done
