
import numpy as np
from scipy.special import erfc, erfcx
from scipy.interpolate import interp1d
from scipy.signal import fftconvolve


# =============================================================================
# Closed-form GaLa implementation
# =============================================================================

def fastGaLa(
    x,      # Energy axis [eV]
    x_0,    # Center of the Gaussian mode [eV]
    sigma,  # Gaussian standard deviation [meV]
    w_h,    # High-energy tail width [meV]
    w_l,    # Low-energy tail width [meV]
    y=0,    # Baseline offset
    A=1     # Peak amplitude
):
    """
    Closed-form implementation of the Gaussian-Asymmetric Laplace (GaLa)
    lineshape.

    Recommended implementation for fitting experimental spectra.

    Returns a peak-normalized spectrum with amplitude A and baseline y.
    """
    # Convert energies to meV and shift to the peak center
    x = np.asarray(x) * 1000
    x0 = x_0 * 1000
    x = x - x0

    # - erfc-based implementation - #
    # L = (
    #     np.exp(x / w_l + sigma**2 / (2 * w_l**2))
    #     * erfc((x + sigma**2 / w_l) / (np.sqrt(2) * sigma))
    # )

    # R = (
    #     np.exp(-x / w_h + sigma**2 / (2 * w_h**2))
    #     * erfc((-x + sigma**2 / w_h) / (np.sqrt(2) * sigma))
    # )

    # - a more numerically stable implementation using erfcx instead of erfc - #
    L = np.exp(-x**2/(2*sigma**2)) * erfcx(
        (x + sigma**2/w_l)/(np.sqrt(2)*sigma)
    )
    R = np.exp(-x**2/(2*sigma**2)) * erfcx(
        (-x + sigma**2/w_h)/(np.sqrt(2)*sigma)
    )

    f = L + R

    if not np.all(np.isfinite(f)) or np.max(f) == 0:
        return np.full_like(x, np.nan, dtype=float)

    return A * f / np.max(f) + y


# =============================================================================
# Basic distributions
# =============================================================================

def Gaussian(x, x_0, sigma, y=0, A=1):
    """Gaussian distribution."""
    return A * np.exp(-(x - x_0)**2 / (2 * sigma**2)) + y


def Laplace(x, x_0, w, y=0, A=1):
    """Laplace distribution."""
    return A * np.exp(-np.abs(x - x_0) / w) + y


def Lorentzian(x, x_0, gamma, y=0, A=1):
    """Lorentzian distribution."""
    return A / (1 + ((x - x_0) / gamma)**2) + y


# =============================================================================
# Asymmetric distributions
# =============================================================================

def asymmetricLaplace(x, w_l, w_h):
    """Asymmetric Laplace distribution with separate low- and high-energy tails."""
    return np.where(
        x < 0,
        Laplace(x, 0, w_l),
        Laplace(x, 0, w_h)
    )


def LorentzianLaplace(x, gamma, w_h):
    """Lorentzian low-energy tail and Laplace high-energy tail."""
    return np.where(
        x < 0,
        Lorentzian(x, 0, gamma),
        Laplace(x, 0, w_h)
    )


# =============================================================================
# Numerical convolution
# =============================================================================

def _convolve(x, sigma, tail):
    """Convolve a Gaussian with a supplied tail distribution."""
    dx = x[1] - x[0]

    G = Gaussian(x, 0, sigma)
    G /= np.sum(G) * dx
    tail = tail / (np.sum(tail) * dx)

    conv = fftconvolve(G, tail, mode="same")
    conv /= np.max(conv)

    return interp1d(
        x,
        conv,
        kind="cubic",
        bounds_error=False,
        fill_value=0.0
    )


# =============================================================================
# Numerical GaLa and GaLx implementations
# =============================================================================

def GaLa(x, x_0, sigma, w_h, w_l, y=0, A=1):
    """
    Numerical Gaussian-Asymmetric Laplace (GaLa) lineshape.

    Uses FFT convolution of a Gaussian and an asymmetric Laplace distribution.
    """
    x_grid = np.linspace(-1000, 1000, 4096)
    tail = asymmetricLaplace(x_grid, w_l, w_h)

    return A * _convolve(x_grid, sigma, tail)((x - x_0) * 1000) + y


def GaLx(x, x_0, sigma, gamma, w_h, y=0, A=1):
    """
    Numerical Gaussian-Lorentzian (GaLx) lineshape.

    Uses FFT convolution of a Gaussian and a Lorentzian-Laplace distribution.
    The model can approximate the non-Laplacian low-energy tail observed
    in 2ML nanoplatelets.
    """
    x_grid = np.linspace(-1000, 1000, 4096)
    tail = LorentzianLaplace(x_grid, gamma, w_h)

    return A * _convolve(x_grid, sigma, tail)((x - x_0) * 1000) + y
