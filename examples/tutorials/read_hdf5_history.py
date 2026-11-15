"""Print the scalar operating conditions and coefficients in an HDF5 history."""

import argparse

import h5py


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("history_file", nargs="?", default="analysis_history.h5")
    args = parser.parse_args()

    with h5py.File(args.history_file, "r") as history:
        keys = sorted(
            (key for key in history if key.startswith("iteration_")),
            key=lambda key: int(key.rsplit("_", 1)[1]),
        )
        for key in keys:
            iteration = history[key]
            alpha = iteration["operating_conditions/alpha"][()]
            cl = iteration["coefficients/cl"][()]
            cd = iteration["coefficients/cd"][()]
            print(f"{key:12s}  {alpha:6.2f}  {cl:9.5f}  {cd:9.6f}")


if __name__ == "__main__":
    main()
