import pandas as pd
import scanpy as sc

def prep_anndata(adata):
    
    def fix_dtypes(adata):
        df = pd.DataFrame(adata.X.A, index=adata.obs_names, columns=adata.var_names)
        df = df.join(adata.obs)
        return sc.AnnData(df[adata.var_names], obs=df.drop(columns=adata.var_names))

    adata = fix_dtypes(adata)
    
    return adata