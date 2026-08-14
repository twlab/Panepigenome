# ==========================================
# evaluate.py — diagnostics and plots
# ==========================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import logging
from scipy.stats import gaussian_kde

log = logging.getLogger(__name__)


def check_bottleneck_normality(z_mean_vals:    np.ndarray,
                                z_sampled:      np.ndarray,
                                save_path:      str = None):
    """
    Full normality diagnostic for bottleneck:
    - Histogram vs N(0,1)
    - QQ plot
    - Statistical tests
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    z_flat = z_sampled.flatten()

    # 1. Histogram
    axes[0].hist(z_flat, bins=50, density=True,
                 alpha=0.7, color='steelblue',
                 label='bottleneck z')
    x_range = np.linspace(-4, 4, 100)
    axes[0].plot(x_range, stats.norm.pdf(x_range, 0, 1),
                 'r-', linewidth=2, label='N(0,1)')
    axes[0].set_title('Bottleneck Distribution')
    axes[0].legend()

    # 2. QQ plot
    stats.probplot(z_flat, dist='norm', plot=axes[1])
    axes[1].set_title('QQ Plot vs N(0,1)')

    # 3. Test results
    sw_stat, sw_p = stats.shapiro(z_flat[:5000])
    da_stat, da_p = stats.normaltest(z_flat)
    skewness      = stats.skew(z_flat)
    kurtosis      = stats.kurtosis(z_flat)

    text = (
        f"Shapiro-Wilk p: {sw_p:.4f} "
        f"{'✓' if sw_p > 0.05 else '✗'}\n"
        f"D'Agostino p:   {da_p:.4f} "
        f"{'✓' if da_p > 0.05 else '✗'}\n\n"
        f"Mean:     {z_flat.mean():.4f}  (ideal=0)\n"
        f"Std:      {z_flat.std():.4f}   (ideal=1)\n"
        f"Skewness: {skewness:.4f}        (ideal=0)\n"
        f"Kurtosis: {kurtosis:.4f}        (ideal=0)\n\n"
        f"{'✓ NORMAL' if sw_p > 0.05 else '✗ NOT NORMAL — increase beta'}"
    )
    axes[2].text(0.1, 0.5, text,
                 transform    = axes[2].transAxes,
                 fontsize     = 11,
                 verticalalignment = 'center',
                 fontfamily   = 'monospace',
                 bbox         = dict(boxstyle='round',
                                     facecolor='lightblue'))
    axes[2].axis('off')
    axes[2].set_title('Normality Tests')

    plt.suptitle('Bottleneck Normality Check',
                 fontweight='bold', fontsize=13)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        log.info(f"Saved normality plot → {save_path}")
    plt.show()


def plot_training_history(history, save_path: str = None):
    """Plot training and validation loss"""
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(history.history['loss'],     label='train loss')
    ax.plot(history.history['val_loss'], label='val loss')
    ax.set_xlabel('Epoch')
    ax.set_ylabel('Loss')
    ax.set_title('VAE Training Loss')
    ax.legend()

    best_epoch = np.argmin(history.history['val_loss'])
    ax.axvline(best_epoch, color='red', linestyle='--',
               label=f'best epoch={best_epoch+1}')
    ax.legend()

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        log.info(f"Saved training plot → {save_path}")
    plt.show()


def plot_compressed_distribution(df:        pd.DataFrame,
                                  col:       str = 'compressed_1D',
                                  group_col: str = 'present',
                                  save_path: str = None):
    """KDE plot of compressed values by group"""
    fig, ax = plt.subplots(figsize=(10, 6))

    for group in sorted(df[group_col].unique()):
        mask = df[group_col] == group
        vals = df.loc[mask, col].dropna()
        ax.hist(vals, bins=100, density=True,
                alpha=0.4, label=f'{group_col}={group}')
        
        kde = gaussian_kde(vals[:100_000])  # sample for speed
        x   = np.linspace(vals.min(), vals.max(), 500)
        ax.plot(x, kde(x), linewidth=2)
    ax.set_xlabel(col)
    ax.set_ylabel('Density')
    ax.set_title(f'Distribution of {col} by {group_col}')
    ax.legend(title=group_col)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        log.info(f"Saved distribution plot → {save_path}")
    plt.show()


def plot_scatter(df:        pd.DataFrame,
                 x_col:     str,
                 y_col:     str,
                 group_col: str = 'present',
                 n_sample:  int = 100_000,
                 save_path: str = None):
    """Scatter plot of two columns colored by group"""
    df_plot = df.sample(n=min(n_sample, len(df)),
                        random_state=42)

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.scatterplot(data   = df_plot,
                    x      = x_col,
                    y      = y_col,
                    hue    = group_col,
                    alpha  = 0.3,
                    s      = 5,
                    palette = "tab20",
                    ax     = ax)
    ax.set_title(f'{y_col} vs {x_col}')
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        log.info(f"Saved scatter plot → {save_path}")
    plt.show()