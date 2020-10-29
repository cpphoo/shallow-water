#!/bin/bash

# bash script to run strong scaling 
make PLATFORM=graphite src/lshallow
make PLATFORM=graphite src/lshallow_parallel

mkdir strong_scaling

problem_size=1000
echo "Running strong scaling with the serial code"
src/lshallow tests.lua dam $problem_size > strong_scaling/serial.txt

for num_thread in 2 3 4 5 6 7 8 9 10
do
    echo "Running strong scaling with "$num_thread" threads"
    OMP_NUM_THREADS=$num_thread src/lshallow_parallel tests.lua dam $problem_size > strong_scaling/threads_$num_thread.txt
done

# python util/compile_strong_scaling.py --dir strong_scaling --num_threads 2 3 4 5