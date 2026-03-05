import numpy as np
import pandas as pd
from sklearn.preprocessing import scale

# nnPCA implementation
def nsprcomp(X, ncomp=1, center=True, scale_=False, nneg=True, nrestart=5, em_tol=1e-3, em_maxiter=100, verbosity=0):
    """
    Non-negative Sparse PCA (simplified)
    Returns: dict with 'x' as projected components
    """
    n_samples, n_genes = X.shape
    X_scaled = scale(X, axis=0, with_mean=center, with_std=scale_, copy=True)
    
    # handle zero std columns
    sc = X_scaled.std(axis=0)
    if any(sc == 0):
        raise ValueError("Cannot rescale a zero column")

    W = np.zeros((n_genes, ncomp))
    for cc in range(ncomp):
        obj_opt = float('-inf')
        w_opt = None
        for rr in range(nrestart):
            w = np.abs(np.random.normal(size=n_genes)) if nneg else np.random.normal(size=n_genes)
            w = w / np.linalg.norm(w)

            # EM iterations
            for ii in range(em_maxiter):
                z = X_scaled @ w
                obj = float(z.T @ z)
                if obj != 0 and abs(obj - obj_opt) / obj < em_tol:
                    break
                w_star = X_scaled.T @ z / (obj if obj != 0 else 1e-12)
                if nneg:
                    w_star[w_star < 0] = 0
                w = w_star / (np.linalg.norm(w_star) if np.linalg.norm(w_star) != 0 else 1e-12)

            # keep best restart
            obj_cur = float((X_scaled @ w).T @ (X_scaled @ w))
            if obj_cur > obj_opt:
                obj_opt = obj_cur
                w_opt = w.copy()
        W[:, cc] = w_opt

    # project X
    X_scores = X_scaled @ W
    return {'x': X_scores}

# parse GMT file
def parse_gmt(gmt_file):
    genesets = {}
    with open(gmt_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            name = parts[0]
            genes = parts[2:]  # skip first 2 columns
            genesets[name] = genes
    return genesets

# compute PCA1 per gene set
def get_nnPCA_result(geneExp, genes):
    # match genes in expression matrix
    available_genes = [g for g in genes if g in geneExp.columns]
    if len(available_genes) == 0:
        return np.zeros(geneExp.shape[0])
    X = geneExp[available_genes].values
    # drop all-zero columns
    X = X[:, np.any(X != 0, axis=0)]
    if X.shape[1] == 0:
        return np.zeros(geneExp.shape[0])
    result = nsprcomp(X, ncomp=1)
    return result['x'][:, 0]

# main function
def run_nnPCA(geneExp: pd.DataFrame, gmt_file: str, dimension: int = 1):
    print("Parsing GMT file...")
    genesets = parse_gmt(gmt_file)
    print(f"Number of gene sets: {len(genesets)}")

    print("Computing nnPCA scores...")
    result_dict = {}
    for i, (name, genes) in enumerate(genesets.items(), 1):
        print(f"Processing gene set {i}/{len(genesets)}: {name}")
        result_dict[name] = get_nnPCA_result(geneExp, genes)

    scores = pd.DataFrame(result_dict, index=geneExp.index)

    if dimension == 1:
        # return PCA1 only
        return scores
    else:
        return scores  # can be extended later to return more components


def execute_nnPCA_single(geneExp: pd.DataFrame, gmt_file: str, score_name: str = "nnPCA"):
    """
    Compute nnPCA score for the first gene set in the GMT file.
    Returns a DataFrame with index as samples and one column with the score.
    """
    # Parse GMT and take only the first gene set
    genesets = parse_gmt(gmt_file)
    first_name, first_genes = next(iter(genesets.items()))

    # Run nnPCA for that single gene set
    score_values = get_nnPCA_result(geneExp, first_genes)

    # Wrap in DataFrame
    result = pd.DataFrame(score_values, index=geneExp.index, columns=[score_name])
    return result

def execute_nnPCA_multidim(
    geneExp: pd.DataFrame,
    genes: list,
    n_components: int = 2,
    score_prefix: str = "nnPCA"
):
    score_values = get_nnPCA_result(
        geneExp[genes],
        n_components=n_components

    )
    
    columns = [f"{score_prefix}_score" for i in range(n_components)]
    result = pd.DataFrame(score_values, index=geneExp.index, columns = columns)
    
    return result
