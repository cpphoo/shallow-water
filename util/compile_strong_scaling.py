import argparse

from parse_output import parse

import os

import matplotlib.pyplot as plt

import numpy as np

def main(args):
    serial_time = parse(os.path.join(args.dir, 'serial.txt'))
    print(serial_time)
    num_threads = sorted(args.num_threads)

    times = [serial_time]
    for num_thread in num_threads:
        times.append(parse(os.path.join(args.dir, f'threads_{num_thread}.txt')))

    
    times = np.array(times)
    speedup = serial_time / times
    print('Speedup: ', speedup)
    plt.plot([1] + num_threads, speedup, '-*', label='Speedup')
    plt.xticks([1] + num_threads)
    plt.ylabel('Speedup')
    plt.xlabel('Number of threads')
    plt.grid()
    plt.legend()

    plt.show()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', required=True, type=str,
                help="directory that have all output files")   
    parser.add_argument('--num_threads', nargs='+', type=int, 
                        help='number of threads')
    args = parser.parse_args()
    main(args) 
