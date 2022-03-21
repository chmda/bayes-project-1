import argparse
import os
import sys
from typing import Dict
from mcmc import Data
from mcmc.data import Coefficients
from mcmc.mcmc import MCMC


def export_to_latex(summary: Dict[str, Coefficients]) -> str:
    result = "\\begin{center}\n\\begin{tabular}"
    result += "{|" + "c|" * (1 + 8 + 3 + 1 + 1) + "}\n\\hline\n"
    # add headers
    result += "&".join(
        [
            "metric",
            *[f"beta[{i + 1}]" for i in range(8)],
            *[f"gamma[{i + 1}]" for i in range(3)],
            "phi",
            "theta",
        ]
    )
    result += " \\\\\n"
    result += "\\hline\n"
    # print results
    for metric, coeff in summary.items():
        betas = coeff.beta
        gammas = coeff.gamma
        phi = coeff.phi
        theta = coeff.theta
        f = lambda x: "{:.2e}".format(x)
        result += "&".join([metric, *map(f, betas), *map(f, gammas), f(phi), f(theta)])
        result += " \\\\\n"
        result += "\\hline\n"

    result += "\\end{tabular}\n\\end{center}"
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()  # TODO : description
    parser.add_argument("datafile", type=str, help="Data file")
    parser.add_argument("initfile", type=str, help="Initial coefficients")
    parser.add_argument("-n", type=int, default=10**4, help="Chain length")
    parser.add_argument(
        "-psd", "--prop-sd", type=float, default=1, help="Std for random walk in M.H."
    )
    parser.add_argument("-b", "--burnin", type=int, default=1000)
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

    # print summary
    summary = mcmc.summary(args.burnin)
    table_format = "{:<10} " * (1 + 8 + 3 + 1 + 1)
    print(
        table_format.format(
            "metric",
            *[f"beta[{i + 1}]" for i in range(8)],
            *[f"gamma[{i + 1}]" for i in range(3)],
            "phi",
            "theta",
        )
    )

    for metric, coeff in summary.items():
        betas = coeff.beta
        gammas = coeff.gamma
        phi = coeff.phi
        theta = coeff.theta
        f = lambda x: "{:.2e}".format(x)
        print(
            table_format.format(
                metric, *map(f, betas), *map(f, gammas), f(phi), f(theta)
            )
        )

    with open("report/output.tex", "w+") as f:
        f.write(export_to_latex(summary))
