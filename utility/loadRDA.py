import pyreadr
from pathlib import Path
import pandas as pd
import argparse
import sys

def rda_to_csv(rda_path: Path):
    if not rda_path.exists():
        print(f"Error: file not found: {rda_path}")
        sys.exit(1)
        
    result = pyreadr.read_r(rda_path)
    df = next(iter(result.values()))
    df = pd.DataFrame(df)
    csv_path = rda_path.with_suffix(".csv")
    df.to_csv(csv_path, index=False)
    print(f"Converted {rda_path} -> {csv_path}")
    
    
def main():
    parser = argparse.ArgumentParser(description = "Convert a single .rda file to CSV")
    parser.add_argument("rda_file", help="Path to .rda file to convert")
    args = parser.parse_args()
        
    rda_to_csv(Path(args.rda_file))
        
    
if __name__ == "__main__":
    main()
    