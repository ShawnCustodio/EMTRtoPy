import requests
import tempfile
from pathlib import Path
import pandas as pd
import pyreadr

# 1. Download geneExp.rda from Zenodo
url = "https://zenodo.org/records/17353682/files/geneExp.rda"
tmp = tempfile.NamedTemporaryFile(suffix=".rda", delete=False)
destfile = tmp.name
tmp.close()

resp = requests.get(url, stream=True)
resp.raise_for_status()
with open(destfile, "wb") as f:
    for chunk in resp.iter_content(chunk_size=8192):
        if chunk:
            f.write(chunk)

# 2. Load .rda
result = pyreadr.read_r(destfile)
geneExp = result.get("geneExp", next(iter(result.values())))
geneExp = pd.DataFrame(geneExp)

# 3. Preview first 4 rows, first 5 columns
print("dim:", geneExp.shape)
print(geneExp.iloc[:4, :5])

# 4. Optionally save CSV
out_csv = Path("data/geneExp.csv")
geneExp.to_csv(out_csv, index=True)