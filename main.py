import argparse
import os
import sys
from mcmc import Data

if __name__ == "__main__":
    parser = argparse.ArgumentParser() # TODO : description
    parser.add_argument("infile", type=str, help="Data file")
    args = parser.parse_args()

    if not os.path.isfile(args.infile):
        print(args.infile, "is not a valid file.")
        sys.exit(1)

    with open(args.infile, "r") as f:
        data = Data.from_json(f.read())

    print("Nothing to do...")