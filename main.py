import argparse
import os
import sys
from mcmc import Data
from mcmc.data import Coefficients
from mcmc.mcmc import MCMC

if __name__ == "__main__":
    parser = argparse.ArgumentParser()  # TODO : description
    parser.add_argument("datafile", type=str, help="Data file")
    parser.add_argument("initfile", type=str, help="Initial coefficients")
    parser.add_argument("-n", type=int, default=10**4, help="Chain length")
    parser.add_argument(
        "-psd", "--prop-sd", type=float, default=1, help="Std for random walk in M.H."
    )
    args = parser.parse_args()

    if not os.path.isfile(args.datafile):
        print(args.datafile, "is not a valid file.")
        sys.exit(1)

    if not os.path.isfile(args.initfile):
        print(args.initfile, "is not a valid file.")
        sys.exit(1)

    with open(args.datafile, "r") as f:
        data = Data.from_json(f.read())

    with open(args.initfile, "r") as f:
        init = Coefficients.from_json(f.read())

    mcmc = MCMC(init, args.n)
    print("Training the model...")
    mcmc.fit(data, args.prop_sd)
    print("Done")
    print("Summary:")
    print(mcmc.summary(burnin=1000))
