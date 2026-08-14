# ==========================================
# config.py — all parameters in one place
# ==========================================
import os

# ── Environment ───────────────────────────
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ["KERAS_BACKEND"]                 = "jax"

# ── Data folder ───────────────────────────
DATA_DIR        = "data"          # ← change this

# ── Data paths ────────────────────────────
RAW_INPUT       = os.path.join(DATA_DIR, "input_filtered3_M2.log")
PARQUET_PATH    = os.path.join(DATA_DIR, "input_filtered3_M2.parquet")
NA_ROWS_PATH    = os.path.join(DATA_DIR, "input_filtered3_M2_na_rows.csv")
XTRAIN_PATH     = os.path.join(DATA_DIR, "X_train_compressed.npz")
OUTPUT_PARQUET  = os.path.join(DATA_DIR, "input_filtered3_M2_inferred.parquet")
SCALER_PATH     = os.path.join(DATA_DIR, "scaler.pkl")
MODEL_PREFIX    = os.path.join(DATA_DIR, "vae")

# ── Preprocessing ─────────────────────────
NA_VALUE        = -20
RANDOM_SEED     = 42


# ── Sampling ──────────────────────────────
TRAIN_SAMPLE_SIZE = 1e6
INFER_ALL_ROWS    = True

# ── VAE architecture ──────────────────────
INPUT_DIM       = 3
LATENT_DIM      = 1
HIDDEN_DIM      = 32
BETA            = 10.0

# ── Training ──────────────────────────────
EPOCHS               = 3
BATCH_SIZE           = 32
VALIDATION_SPLIT     = 0.1
LEARNING_RATE        = 1e-3
EARLY_STOP_PATIENCE  = 10
REDUCE_LR_PATIENCE   = 5
REDUCE_LR_FACTOR     = 0.5
MIN_LR               = 1e-6

# ── Inference ─────────────────────────────
PREDICT_BATCH_SIZE = 8192