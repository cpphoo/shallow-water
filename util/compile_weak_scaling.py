import argparse

from parse_output import parse

import os

import matplotlib.pyplot as plt

import numpy as np

from sklearn.linear_model import LinearRegression


def main(args):
    serial_time = parse(os.path.join(args.dir, 'serial.txt'))
    print(serial_time)
    num_threads = sorted(args.num_threads)

    parallel_times = [serial_time]
    serial_times = [serial_time]
    for num_thread in num_threads:
        parallel_times.append(
            parse(os.path.join(args.dir, f'threads_{num_thread}.txt')))
        serial_times.append(
            parse(os.path.join(args.dir, f'serial_{num_thread}.txt')))

    parallel_times = np.array(parallel_times)
    serial_times = np.array(serial_times)


    print("Parallel_times: ", parallel_times)
    print("Serial times: ", serial_times)
    speedup = serial_times / parallel_times
    print("Speedup: ", speedup)


    # Fitting a linear regression model 
    X = np.arange(1, 1+ len(speedup)).reshape(-1, 1)
    reg = LinearRegression().fit(X, speedup)

    print(reg.coef_)
    print(reg.intercept_)
    print(num_threads)
    plt.plot([1] + num_threads, speedup, '-*', label='Scaled Speedup')
    
    # plot best fitted line 
    plt.plot(X.flatten(), reg.predict(X).flatten(), '--', label='Best Fitted Line: {:.2f}p+{:.2f}'.format(reg.coef_[0], reg.intercept_))
    plt.xticks([1] + num_threads)
    plt.ylabel('Scaled Speedup')
    plt.xlabel('Number of threads, p')
    plt.grid()
    plt.legend()

    if args.show_plot:
        plt.show()
    else:
        plt.savefig(f'{args.dir}.jpg')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', required=True, type=str,
                        help="directory that have all output files")
    parser.add_argument('--num_threads', nargs='+', type=int,
                        help='number of threads')
    parser.add_argument('--show_plot', action='store_true',
                        help='whether to show the plot')

    args = parser.parse_args()
    main(args)
