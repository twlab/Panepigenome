# ==========================================
# data.py — loading and preprocessing
# ==========================================
import numpy as np
import pandas as pd
import dask.dataframe as dd
from sklearn.preprocessing import StandardScaler
import gc
import logging
import matplotlib.pyplot as plt
import seaborn as sns
 
log = logging.getLogger(__name__)
 
 
def load_and_save_parquet(raw_path: str,
                           parquet_path: str,
                           na_rows_path: str) -> None:
    """
    One-time conversion: CSV → parquet
    Only needs to run once
    """
    log.info(f"Reading raw file: {raw_path}")
    df = pd.read_csv(raw_path, header=0, sep='\t')

    mask = df['present'] == 0
    n_replace  = mask.sum()
    df.loc[mask, 'expression'] = -9.97
    log.info(f"Set expression=-9.97 for {n_replace:,} present=0 rows")
 
    # Save NA rows for inspection
    na_rows = df[df['expression'].isna()]
    na_rows.to_csv(na_rows_path, sep='\t', index=False)
    log.info(f"Saved {len(na_rows)} NA rows to {na_rows_path}")
 
    df = df.dropna(axis=0)
    df.to_parquet(parquet_path, compression='snappy')
    log.info(f"Saved parquet: {parquet_path} ({len(df):,} rows)")
 
 
def load_parquet(parquet_path: str) -> pd.DataFrame:
    """Load parquet into pandas"""
    log.info(f"Loading parquet: {parquet_path}")
    df = pd.read_parquet(parquet_path)
    log.info(f"Loaded {len(df):,} rows")
    return df
 
 
def encode_features(df: pd.DataFrame, na_value: float, scaler: StandardScaler = None):
    """
    Encode categorical + numerical features into X_train matrix
 
    Returns:
        X_train: np.ndarray shape (n, 3)
            columns: [gene_encoded, sample_encoded, expression]
        encoders: dict with LabelEncoders for gene and sample
    """

    # sns.kdeplot(data=df, x='expression', hue='present', fill=True, alpha=0.5,common_norm=False)
    # 1. Create the dummies as booleans (True/False) so we know exactly which are present
    group_dummies = pd.get_dummies(df['present'], dtype=bool)

    # 2. Multiply by expression (True becomes the expression value, False becomes 0.0)
    transformed_cols = group_dummies.multiply(df['expression'], axis=0)

    # 3. Safely replace ONLY the cells that were not present with -10
    # The ~ symbol means "NOT", so this says "where the group is NOT present, put -10"
    transformed_cols = transformed_cols.mask(~group_dummies, na_value)

    # if scaler is None:
    #     scaler = StandardScaler()
    #     log.info("Fitting new scaler to transformed columns...")
    # else:
    #     log.info("Using provided scaler to transform columns...")

    transformed_cols = pd.DataFrame(transformed_cols, columns=transformed_cols.columns)
    log.info("No Scaling applied to transformed columns...") 
    col_x = transformed_cols[0]
    col_y = transformed_cols[1]
    col_z = transformed_cols[2]

    # plot col_x, col_y, col_z to check distributions 
    # sns.kdeplot(col_x, label='col_x')
    # sns.kdeplot(col_y, label='col_y')
    # sns.kdeplot(col_z, label='col_z')
    # plt.legend()
    # plt.savefig("figs/encoded_feature_distributions.png", dpi=150, bbox_inches='tight')
    
    
    X_train = np.column_stack((col_x, col_y, col_z)).astype(np.float32)
 
    log.info(f"X_train shape: {X_train.shape}")
    log.info(f"X_train dtype: {X_train.dtype}")
 

    return X_train, scaler
 
 
def save_xtrain(X_train: np.ndarray, path: str) -> None:
    np.savez_compressed(path, train_data=X_train)
    log.info(f"Saved X_train → {path}  shape={X_train.shape}")
 
 
def load_xtrain(path: str) -> np.ndarray:
    data    = np.load(path)
    X_train = data['train_data']
    log.info(f"Loaded X_train: {X_train.shape}")
    return X_train
 
 
def free_memory(*dfs):
    """Delete dataframes and run GC"""
    for obj in dfs:
        del obj
    gc.collect()
    log.info("Memory freed")

def sample_for_training(df:          pd.DataFrame,
                         n:           int   = 1_000_000,
                         random_seed: int   = 42) -> pd.DataFrame:
    """
    Sample n rows for training/validation
    Stratified by 'present' to preserve class balance
    """
    total = len(df)

    if n >= total:
        log.info(f"Requested {n:,} >= total {total:,} — using all rows")
        return df

    # Stratified sample — preserves present=0,1,2 proportions
    frac = n / total

    df_sample = df.sample(frac=frac, random_state=random_seed)

    log.info(f"Sampled {len(df_sample):,} rows from {total:,} "
             f"(stratified by present)")
    log.info(f"Sample value counts:\n"
             f"{df_sample['present'].value_counts().to_string()}")

    return df_sample