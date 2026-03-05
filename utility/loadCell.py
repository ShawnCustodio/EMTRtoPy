"""Convert a cell annotation .rda into a pandas DataFrame and print head() to terminal.

Usage:
  python followHTML.py --rda data/cell_annotation_file.rda [--out-csv out.csv] [--n 6]

This script prefers `pyreadr` to read .rda files. If not available it
falls back to using `Rscript` to write a CSV which Python then reads.
No HTML is produced; the DataFrame `head()` is printed to stdout.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd


def read_rda_with_pyreadr(rda_path: Path) -> pd.DataFrame:
    import pyreadr

    res = pyreadr.read_r(rda_path)
    if not res:
        raise RuntimeError(f".rda file {rda_path} contained no objects")
    first_key = next(iter(res))
    df = res[first_key]
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)
    return df


def read_rda_with_rscript(rda_path: Path, out_csv: Path) -> None:
    r_exec = shutil.which("Rscript") or shutil.which("R")
    if not r_exec:
        raise RuntimeError("Rscript or R not found on PATH; cannot convert .rda to CSV")

    expr = (
        "e <- new.env();"
        f"load('{rda_path.as_posix()}', envir=e);"
        "objs <- ls(envir=e);"
        "found <- NULL;"
        "for(n in objs){ if(is.data.frame(get(n, envir=e))){ found <- n; break } }"
        "if(is.null(found)) stop('no data.frame found in rda');"
        "write.csv(get(found, envir=e), file='" + out_csv.as_posix() + "', row.names=FALSE)"
    )

    cmd = [r_exec, "-e", expr]
    subprocess.run(cmd, check=True)


def main(argv: list[str] | None = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    p = argparse.ArgumentParser(description="Load .rda annotation and print head()")
    p.add_argument("--rda", required=True, help="Path to cell_annotation_file.rda")
    p.add_argument("--out-csv", help="Optional output CSV path")
    p.add_argument("--n", type=int, default=6, help="Number of rows to show from head()")
    args = p.parse_args(argv)

    rda_path = Path(args.rda)
    if not rda_path.exists():
        print(f".rda file not found: {rda_path}")
        return 2

    out_csv = Path(args.out_csv) if args.out_csv else rda_path.with_suffix('.csv')

    try:
        try:
            df = read_rda_with_pyreadr(rda_path)
            print("Loaded .rda with pyreadr")
        except Exception as e_py:
            print("pyreadr failed or not available, falling back to Rscript:", e_py)
            read_rda_with_rscript(rda_path, out_csv)
            df = pd.read_csv(out_csv)

        # Optionally write CSV if requested or if R was used
        if args.out_csv or 'pyreadr' in sys.modules:
            df.to_csv(out_csv, index=False)

        # Print head() to terminal
        print(df.head(args.n))
        return 0

    except subprocess.CalledProcessError as e:
        print("Rscript failed:", e)
        return 3
    except Exception as e:
        print("Error:", e)
        return 4


if __name__ == '__main__':
    raise SystemExit(main())

## 2.1 Load cell/sample annotation R to Python