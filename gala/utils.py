import numpy as np


# =============================================================================
# GaLa fitting parameters
# =============================================================================

"""
Default fitting parameters and bounds for scipy.optimize.curve_fit.

The initial guess for x_0 should ideally be set to the spectral maximum
to ensure fast and reliable convergence.

The default initial guess (2.4 eV) is suitable for the example data.
"""

GALA_BOUNDS = (
    [0.5, 1.0, 1.0, 1.0, 0.0, 0.99],
    [3.5, 200, 200, 200, 0.00001, 1.01]
)

GALA_P0 = [2.4, 10, 10, 10, 0.0, 1.0]


# =============================================================================
# Spectral conversions
# =============================================================================

"""
Jacobian conversion functions for transforming wavelength spectra
into photon-energy spectra.

Since E = hc / λ is nonlinear, the intensity must be transformed
using the Jacobian to preserve integrated spectral intensity:

    I_E = I_λ * |dλ/dE| = I_λ * λ² / hc

Wavelengths: nm
Energies:    eV
Intensities: counts/nm -> counts/eV
"""

HC_EV_NM = 1239.841984


def nm_to_eV(wavelength):
    """Convert wavelength in nm to photon energy in eV."""
    return HC_EV_NM / np.asarray(wavelength)


def jacobian_nm_to_eV(wavelength, intensity):
    """
    Convert a wavelength spectrum into an energy spectrum.

    Applies the Jacobian transformation and reverses the arrays
    to return energy in ascending order.

    Parameters
    ----------
    wavelength : array-like
        Wavelength values in nm.
    intensity : array-like
        Spectral intensity per nm.

    Returns
    -------
    energy : ndarray
        Photon energy in eV.
    intensity_eV : ndarray
        Spectral intensity per eV.
    """
    wavelength = np.asarray(wavelength)
    intensity = np.asarray(intensity)

    energy = nm_to_eV(wavelength)
    intensity_eV = intensity * wavelength**2 / HC_EV_NM

    return energy[::-1], intensity_eV[::-1]