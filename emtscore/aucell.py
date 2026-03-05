import numpy as np
import pandas as pd

def aucell_score(expr_matrix: pd.DataFrame, genes: list, top_percent: float=0.05):
    genes = [g for g in genes if g in expr_matrix.columns]
    if len(genes) == 0:
        return pd.DataFrame({"SampleId": expr_matrix.index, "AUCell": np.nan})
    
    n_genes = expr_matrix.shape[1]
    max_rank = int(n_genes * top_percent)
    
    scores = []
    
    for sample in expr_matrix.index:
        values = expr_matrix.loc[sample]
        
        ranked_genes = values.sort_values(ascending=False)
        
        top_genes = ranked_genes.iloc[:max_rank].index
        
        enrichment = len(set(genes).intersection(top_genes)) / len(genes)
        
        scores.append(enrichment)
        
    return pd.DataFrame({
        "SampleID": expr_matrix.index,
        "AUCell": scores
    })
    
  
def parse_gmt(gmt_file: str) -> dict:
    genesets = {}
    with open(gmt_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            name = parts[0]
            genes = parts[2:]
            genesets[name] = genes
    return genesets
    
    
def execute_aucell(expr_matrix: pd.DataFrame, gmt_file: str) -> pd.DataFrame:
    genesets = parse_gmt(gmt_file)
    result_dict = {}
    
    for name, genes in genesets.items():
        auc_df = aucell_score(expr_matrix, genes)
        result_dict[name] = auc_df.set_index("SampleID")["AUCell"]
        
    return pd.DataFrame(result_dict, index=expr_matrix.index)
    

def execute_aucell_single(expr_matrix: pd.DataFrame, gmt_file: str, score_name: str, gene_set_index: int = 0):
    # read only one line from GMT
    with open(gmt_file) as f:
        lines = f.readlines()
        parts = lines[gene_set_index].strip().split("\t")
        genes = parts[2:]  # skip name + NA
    
    auc_df = aucell_score(expr_matrix, genes)
    
    result = auc_df.set_index("SampleID")
    result.columns = [score_name]
    
    return result
