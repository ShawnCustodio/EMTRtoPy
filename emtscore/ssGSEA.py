import pandas as pd
import gseapy as gp
import numpy as np


def ssgsea_score(expr_matrix: pd.DataFrame, genes: list, alpha: float =0.25) -> pd.Series:
    genes = [g for g in genes if g in expr_matrix.columns]
    if len(genes) == 0:
        return pd.Series(index=expr_matrix.index, name="ssGSVA", dtype=float)
    
    N = expr_matrix.shape[1]
    scores = []
    
    for sample in expr_matrix.index:
        values = expr_matrix.loc[sample]
        ranked = values.sort_values(ascending=False)
        ranked_genes = ranked.index.values
        ranked_values = ranked.values
        
        
        hits = np.isin(ranked_genes, genes)
        Nh = hits.sum()
        Nm = N - Nh
        
        if Nh == 0:
            scores.append(np.nan)
            continue
        
        weights = np.abs(ranked_values) ** alpha
        sum_hit_weights = weights[hits].sum()
        
        P_hit = np.cumsum(hits * weights / sum_hit_weights)
        P_miss = np.cumsum(~hits / Nm)
        
        running_ES = P_hit - P_miss
        ES = running_ES.sum()
        
        scores.append(ES)
        
    return pd.Series(
        scores, index=expr_matrix.index
    )

def parse_gmt(gmt_file: str) -> dict:
    genesets = {}
    with open(gmt_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            name = parts[0]
            genes = parts[2:]
            genesets[name] = genes
    return genesets




def execute_ssgsva(expr_matrix: pd.DataFrame, gmt_file: str) -> pd.DataFrame:
    genesets = parse_gmt(gmt_file)
    result_dict = {}
    
    for gene_set_name, genes in genesets.items():
        result_dict[gene_set_name] = ssgsea_score(expr_matrix, genes)
        
    return pd.DataFrame(result_dict, index=expr_matrix.index)



def execute_ssgsea_single(expr_matrix: pd.DataFrame, gmt_file: str, score_name: str = "ssGSVA",gene_set_index: int = 0) -> pd.DataFrame:
    genes = []
    with open(gmt_file) as f:
        lines = f.readlines()
        if gene_set_index >= len(lines):
            raise IndexError(f"gene_set_index {gene_set_index} is out of range. GMT has {len(lines)} gene sets.")
        
        parts = lines[gene_set_index].strip().split("\t")
        genes = parts[2:]  
    
    genes = [g for g in genes if g in expr_matrix.columns]
    
    if len(genes) == 0:
        return pd.DataFrame({score_name: [0]*expr_matrix.shape[0]}, index=expr_matrix.index)
    
    scores = ssgsea_score(expr_matrix, genes)
    
    result_df = pd.DataFrame(scores, index=expr_matrix.index, columns=[score_name])
    
    return result_df