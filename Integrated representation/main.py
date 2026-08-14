#!/usr/bin/env python
# ==========================================
# main.py — full pipeline, run end-to-end
# ==========================================
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ["KERAS_BACKEND"]                 = "jax"

import logging
import argparse
import numpy as np
import pandas as pd

import config
from data     import (load_and_save_parquet, load_parquet,
                         encode_features, save_xtrain, load_xtrain, free_memory)
from model    import (build_vae, train_vae,
                       get_bottleneck, save_model)
from evaluate import (check_bottleneck_normality,
                       plot_training_history,
                       plot_compressed_distribution,
                       plot_scatter)

# ── Logging setup ─────────────────────────────────────────
logging.basicConfig(
    level   = logging.INFO,
    format  = '%(asctime)s | %(levelname)s | %(name)s | %(message)s',
    handlers = [
        logging.StreamHandler(),
        logging.FileHandler('pipeline.log')
    ]
)
log = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description='VAE preprocessing pipeline'
    )
    parser.add_argument(
        '--step', type=str,
        choices=['all', 'preprocess', 'encode', 'train', 'infer'],
        default='all',
        help='Which step to run (default: all)'
    )
    parser.add_argument(
        '--skip-csv',  action='store_true',
        help='Skip CSV→parquet conversion (parquet already exists)'
    )
    parser.add_argument(
        '--skip-encode', action='store_true',
        help='Skip feature encoding (X_train.npz already exists)'
    )
    parser.add_argument(
        '--skip-train', action='store_true',
        help='Skip training (load saved model instead)'
    )
    parser.add_argument(
        '--beta', type=float, default=config.BETA,
        help=f'KL weight beta (default: {config.BETA})'
    )
    parser.add_argument(
        '--epochs', type=int, default=config.EPOCHS,
        help=f'Training epochs (default: {config.EPOCHS})'
    )
    return parser.parse_args()


def step_preprocess(skip_csv: bool = False):
    """Step 1: Load raw data → clean → save parquet"""
    log.info("=" * 50)
    log.info("STEP 1: Preprocessing")
    log.info("=" * 50)

    if not skip_csv:
        load_and_save_parquet(
            raw_path    = config.RAW_INPUT,
            parquet_path = config.PARQUET_PATH,
            na_rows_path = config.NA_ROWS_PATH
        )

    df = load_parquet(config.PARQUET_PATH)

    return df


def step_encode(df: pd.DataFrame, skip_encode: bool = False):
    """Step 2: Encode features → X_train matrix"""
    log.info("=" * 50)
    log.info("STEP 2: Feature encoding")
    log.info("=" * 50)

    if skip_encode:
        X_train = load_xtrain(config.XTRAIN_PATH)
        return X_train, df

    X_train, scaler = encode_features(df)
    save_xtrain(X_train, config.XTRAIN_PATH)

    # Free df memory — no longer needed after encoding
    free_memory(df)

    return X_train, scaler


def step_train(X_train: np.ndarray,
               beta:    float = None,
               epochs:  int   = None,
               skip_train: bool = False):
    """Step 3: Build and train VAE"""
    log.info("=" * 50)
    log.info("STEP 3: VAE Training")
    log.info("=" * 50)

    beta   = beta   or config.BETA
    epochs = epochs or config.EPOCHS

    vae, encoder, decoder = build_vae(
        input_dim  = config.INPUT_DIM,
        latent_dim = config.LATENT_DIM,
        hidden_dim = config.HIDDEN_DIM,
        beta       = beta
    )

    if skip_train:
        log.info("Loading saved encoder...")
        from model import load_encoder
        encoder = load_encoder("vae_encoder.keras")
        return None, encoder, None, None

    vae.summary()

    history = train_vae(
        vae              = vae,
        X_train          = X_train,
        epochs           = epochs,
        batch_size       = config.BATCH_SIZE,
        validation_split = config.VALIDATION_SPLIT,
        early_stop_patience = config.EARLY_STOP_PATIENCE,
        reduce_lr_patience  = config.REDUCE_LR_PATIENCE,
        reduce_lr_factor    = config.REDUCE_LR_FACTOR,
        min_lr              = config.MIN_LR
    )

    save_model(vae, encoder, decoder, prefix="vae")
    plot_training_history(history,
                          save_path='figs/training_history.png')

    return vae, encoder, decoder, history


def step_infer(encoder,
               X_test:  np.ndarray,
               df_orig:  pd.DataFrame = None):
    """Step 4: Extract bottleneck → save results"""
    log.info("=" * 50)
    log.info("STEP 4: Inference")
    log.info("=" * 50)

    z_mean_vals, z_log_var_vals, z_sampled = get_bottleneck(
        encoder    = encoder,
        X          = X_test,
        batch_size = config.PREDICT_BATCH_SIZE
    )

    # Normality check
    # check_bottleneck_normality(
    #     z_mean_vals = z_mean_vals,
    #     z_sampled   = z_sampled,
    #     save_path   = 'figs/bottleneck_normality.png'
    # )

    # Attach compressed values to original dataframe
    if df_orig is not None:
        df_orig['compressed_1D'] = z_sampled.flatten()
        
        # Plots
        plot_compressed_distribution(
            df        = df_orig,
            col       = 'compressed_1D',
            group_col = 'present',
            save_path = 'figs/compressed_distribution.png'
        )
        plot_scatter(
            df        = df_orig,
            x_col     = 'expression',
            y_col     = 'compressed_1D',
            group_col = 'present',
            save_path = 'figs/scatter_expression_vs_compressed.png'
        )

        # Save output
        log.info("Saving inferred results...")
        df_orig.to_parquet(config.OUTPUT_PARQUET,
                            compression='snappy')
        log.info(f"Saved results → {config.OUTPUT_PARQUET}")
        # save as CSV for easy inspection
        df_orig.to_csv("inferred_results.csv", index=False)
        log.info(f"Saved results → inferred_results.csv")

        return df_orig

    return z_mean_vals, z_log_var_vals, z_sampled


# ==========================================
# Main
# ==========================================
def main():
    args = parse_args()

    # ── Step 1: Load full dataset ─────────────────────────
    df_full = step_preprocess(skip_csv=args.skip_csv)
    log.info(f"Full dataset: {len(df_full):,} rows")

    # ── Step 2: Sample for training ───────────────────────
    from data import sample_for_training, encode_features, \
                     save_xtrain, load_xtrain, free_memory

    if not args.skip_encode:

        # Sample 1M rows for training
        df_train_sample = sample_for_training(
            df          = df_full,
            n           = config.TRAIN_SAMPLE_SIZE,   # 1_000_000
            random_seed = config.RANDOM_SEED
        )
        log.info(f"Training sample: {len(df_train_sample):,} rows")

        # Encode sample → X_train (1M rows)
        log.info("Encoding training sample using NA value = %f...", config.NA_VALUE)
        X_train, scaler = encode_features(df_train_sample, na_value=config.NA_VALUE)
        save_xtrain(X_train, config.XTRAIN_PATH)

        # Free sample — no longer needed
        free_memory(df_train_sample)

    else:
        X_train  = load_xtrain(config.XTRAIN_PATH)


    log.info(f"X_train shape: {X_train.shape}")   # (1_000_000, 3)

    # ── Step 3: Train on 1M sample ────────────────────────
    vae, encoder, decoder, history = step_train(
        X_train    = X_train,       # ← 1M rows only
        beta       = args.beta,
        epochs     = args.epochs,
        skip_train = args.skip_train
    )

    # Free training data — done with it
    free_memory(X_train)

    # ── Step 4: Infer on ALL rows ─────────────────────────
    if config.INFER_ALL_ROWS:
        log.info(f"Encoding ALL {len(df_full):,} rows for inference...")

        # Encode full dataset using same encoders
        X_full, _ = encode_features(df_full, na_value=config.NA_VALUE, scaler=scaler)
        log.info(f"X_full shape: {X_full.shape}")      # (1B rows, 3)
    else:
        log.info("Skipping encoding of full dataset — using training sample for inference")
        X_full = X_train
        df_full = df_train_sample

    # Run inference on all rows in batches
    df_result = step_infer(
        encoder  = encoder,
        X_test  = X_full,     # ← ALL rows
        df_orig  = df_full
    )

    log.info("Pipeline complete!")

if __name__ == "__main__":
    main()