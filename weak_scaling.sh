#!/bin/bash

# bash script to run strong scaling 
make PLATFORM=graphite src/lshallow
make PLATFORM=graphite src/lshallow_parallel

mkdir weak_scaling

problem_size=500
echo "Running weak scaling with the serial code"
src/lshallow tests.lua dam $problem_size > weak_scaling/serial.txt

for num_thread in 2 3 4 5 6 7 8 9 10
do  
    local_problem_size=$(python -c "import math; print(int(math.sqrt($num_thread) * $problem_size))")
    echo "Running weak scaling with problem size "$local_problem_size
    echo "Running weak scaling with the serial code"
    src/lshallow tests.lua dam $local_problem_size > weak_scaling/serial_$num_thread.txt
    echo "Running weak scaling with "$num_thread" threads"
    OMP_NUM_THREADS=$num_thread src/lshallow_parallel tests.lua dam $local_problem_size > weak_scaling/threads_$num_thread.txt
done

python util/compile_weak_scaling.py --dir weak_scaling --num_threads 2 3 4 5 6 7 8 9 10