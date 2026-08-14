# ==========================================
# model.py — VAE architecture
# ==========================================
import keras
import keras.ops as ops
from keras.layers import Input, Dense
from keras.models import Model
from keras.callbacks import EarlyStopping, ReduceLROnPlateau
import numpy as np
import logging

log = logging.getLogger(__name__)


# ==========================================
# Custom Layers
# ==========================================
class SamplingLayer(keras.layers.Layer):
    """
    Reparameterization trick:
    z = μ + σ × ε,  ε ~ N(0,1)
    Allows gradients to flow through sampling
    """
    def __init__(self, name="z", **kwargs):
        super().__init__(name=name, **kwargs)
        self.seed_generator = keras.random.SeedGenerator(seed=1337)

    def call(self, inputs):
        z_mean, z_log_var = inputs
        batch   = ops.shape(z_mean)[0]
        dim     = ops.shape(z_mean)[1]
        epsilon = keras.random.normal(
            shape = (batch, dim),
            seed  = self.seed_generator
        )
        return z_mean + ops.exp(0.5 * z_log_var) * epsilon


class VAELoss(keras.layers.Layer):
    """
    VAE loss as a layer:
    Total = Reconstruction (MSE) + β × KL divergence

    KL divergence pushes bottleneck toward N(0,1)
    Higher beta = stronger normality constraint
    """
    def __init__(self, beta=1.0, **kwargs):
        super().__init__(**kwargs)
        self.beta = beta

    def call(self, inputs):
        x_input, x_decoded, z_mean, z_log_var = inputs

        # Reconstruction loss (MSE)
        recon_loss = ops.mean(
            ops.sum(
                ops.square(x_input - x_decoded),
                axis=-1
            )
        )

        # KL divergence: KL(N(μ,σ²) || N(0,1))
        kl_loss = -0.5 * ops.mean(
            ops.sum(
                1 + z_log_var
                  - ops.square(z_mean)
                  - ops.exp(z_log_var),
                axis=-1
            )
        )

        self.add_loss(recon_loss + self.beta * kl_loss)
        return x_decoded


# ==========================================
# Build VAE
# ==========================================
def build_vae(input_dim:  int   = 3,
              latent_dim: int   = 1,
              hidden_dim: int   = 32,
              beta:       float = 1.0):
    """
    Build VAE with:
    Encoder: input_dim → hidden_dim → z_mean + z_log_var → z
    Decoder: z → hidden_dim → input_dim

    Returns:
        vae:     full model for training
        encoder: model for bottleneck extraction
        decoder: model for generation
    """
    # ── Encoder ───────────────────────────────────────────
    input_layer = Input(shape=(input_dim,), name="encoder_input")
    enc_hidden  = Dense(hidden_dim, activation='relu',
                        name="encoder_hidden")(input_layer)
    z_mean      = Dense(latent_dim, activation='linear',
                        name="z_mean")(enc_hidden)
    z_log_var   = Dense(latent_dim, activation='linear',
                        name="z_log_var")(enc_hidden)
    z           = SamplingLayer(name="z")([z_mean, z_log_var])

    encoder = Model(
        inputs  = input_layer,
        outputs = [z_mean, z_log_var, z],
        name    = "Encoder"
    )

    # ── Decoder ───────────────────────────────────────────
    decoder_input = Input(shape=(latent_dim,), name="decoder_input")
    dec_hidden    = Dense(hidden_dim, activation='relu',
                          name="decoder_hidden")(decoder_input)
    output_layer  = Dense(input_dim,  activation='linear',
                          name="decoder_output")(dec_hidden)

    decoder = Model(
        inputs  = decoder_input,
        outputs = output_layer,
        name    = "Decoder"
    )

    # ── Full VAE ───────────────────────────────────────────
    vae_output = decoder(z)
    vae_output = VAELoss(beta=beta, name="vae_loss")(
                     [input_layer, vae_output, z_mean, z_log_var]
                 )

    vae = Model(
        inputs  = input_layer,
        outputs = vae_output,
        name    = "VAE"
    )
    vae.compile(optimizer=keras.optimizers.Adam())

    log.info("VAE built successfully")
    log.info(f"  input_dim={input_dim}, latent_dim={latent_dim}, "
             f"hidden_dim={hidden_dim}, beta={beta}")

    return vae, encoder, decoder


# ==========================================
# Training
# ==========================================
def train_vae(vae,
              X_train:          np.ndarray,
              epochs:           int   = 10,
              batch_size:       int   = 32,
              validation_split: float = 0.1,
              early_stop_patience: int = 10,
              reduce_lr_patience:  int = 5,
              reduce_lr_factor:  float = 0.5,
              min_lr:            float = 1e-6):
    """Train VAE with early stopping and LR reduction"""

    callbacks = [
        EarlyStopping(
            monitor              = 'val_loss',
            patience             = early_stop_patience,
            restore_best_weights = True,
            verbose              = 1
        ),
        ReduceLROnPlateau(
            monitor  = 'val_loss',
            factor   = reduce_lr_factor,
            patience = reduce_lr_patience,
            min_lr   = min_lr,
            verbose  = 1
        )
    ]

    log.info(f"Training VAE: {len(X_train):,} samples, "
             f"{epochs} epochs, batch_size={batch_size}")

    history = vae.fit(
        X_train, X_train,
        epochs           = epochs,
        batch_size       = batch_size,
        validation_split = validation_split,
        callbacks        = callbacks,
        verbose          = 1
    )

    best_epoch = np.argmin(history.history['val_loss']) + 1
    best_loss  = min(history.history['val_loss'])
    log.info(f"Training complete — best epoch={best_epoch}, "
             f"val_loss={best_loss:.4f}")

    return history


# ==========================================
# Inference
# ==========================================
def get_bottleneck(encoder,
                   X:          np.ndarray,
                   batch_size: int = 8192):
    """
    Extract bottleneck (z_mean, z_log_var, z_sampled)
    in batches to handle large datasets
    """
    log.info(f"Extracting bottleneck for {len(X):,} samples...")

    outputs = encoder.predict(
        X,
        batch_size = batch_size,
        verbose    = 1
    )

    z_mean_vals    = outputs[0]
    z_log_var_vals = outputs[1]
    z_sampled      = outputs[2]

    log.info("Bottleneck statistics:")
    log.info(f"  z_mean:    mean={z_mean_vals.mean():.4f}  "
             f"std={z_mean_vals.std():.4f}  (ideal: 0, 1)")
    log.info(f"  z_sampled: mean={z_sampled.mean():.4f}  "
             f"std={z_sampled.std():.4f}  (ideal: 0, 1)")

    return z_mean_vals, z_log_var_vals, z_sampled


def save_model(vae, encoder, decoder, prefix: str = "vae"):
    """Save all three models"""
    vae.save(f"{prefix}_full.keras")
    encoder.save(f"{prefix}_encoder.keras")
    decoder.save(f"{prefix}_decoder.keras")
    log.info(f"Models saved with prefix '{prefix}'")


def load_encoder(path: str):
    """Load encoder for inference only"""
    encoder = keras.models.load_model(
        path,
        custom_objects = {
            'SamplingLayer': SamplingLayer,
            'VAELoss':       VAELoss
        }
    )
    log.info(f"Encoder loaded from {path}")
    return encoder