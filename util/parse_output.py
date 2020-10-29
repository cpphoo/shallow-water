# Contains a function used to parse the output from running src/lshallow
import argparse

def parse(filename):
    '''
        filename: A .txt file containing the output from running src/lshallow or src/lshallow_parallel

        parse the output from running src/lshallow and return the total runtime
    '''
    with open(filename) as f:
        s = f.readlines()[-1]

    t = float(s.rstrip().split(':')[-1])
    return t

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--filename', required=True, type=str, help='filename')
    args = parser.parse_args()
    parse(args.filename)    