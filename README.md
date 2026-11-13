
# GaLa Spectrum Fitting

A lightweight Python toolkit for fitting photoluminescence spectra using the Gaussian-Asymmetric Laplace (GaLa) lineshape.

The GaLa model describes asymmetric emission spectra through the convolution of a Gaussian distribution with asymmetric exponential tails. It provides a compact framework for characterizing spectral broadening and deviations from ideal ensemble emission.

The model is described in:

> L. Luber, D. Benedikt, A. Singldinger, A. S. Urban,  
> *Universal Statistics of Nanocrystal Photoluminescence*.  
> [Preprint, DOI: 10.21203/rs.3.rs-8603272/v1](https://doi.org/10.21203/rs.3.rs-8603272/v1).

---

## Getting started

1. Clone or download this repository.
2. Open a terminal in the repository directory.
3. Install the required dependencies:

   ```bash
   python -m pip install -r requirements.txt
   ```

4. Open `fitGaLa.ipynb` in Jupyter Notebook or JupyterLab.
5. Follow the instructions to import, fit, and visualize an experimental spectrum.

## Repository contents

```text
gala-spectrum-fitting/
├── README.md
├── requirements.txt
├── logo.png
├── fitGaLa.ipynb
└── gala/
    ├── __init__.py
    ├── models.py
    └── utils.py
```

## Input data

Experimental spectra can be provided as CSV or whitespace-delimited TXT files.

The input data should contain two columns:

| Column | Description | Unit |
|---|---|---|
| `wavelength` | Emission wavelength | nm |
| `intensity` | Spectral intensity | counts/nm |

Alternatively, spectra already converted to the energy domain can be supplied using `energy` (eV) and `intensity` (counts/eV).

For wavelength-based input, the notebook applies the Jacobian transformation to convert the spectrum into the energy domain while preserving its integrated intensity.

Spectra are normalized before fitting.

Example data are provided in the `data/` directory.

## Citation

If you use GaLa in your research, please cite the associated GaLa publication:

> L. Luber, D. Benedikt, A. Singldinger, A. S. Urban,  
> *Universal Statistics of Nanocrystal Photoluminescence*.

Please use the final publication DOI once available.

## Questions and issues

For bug reports, questions, or suggestions, please open an issue in this repository.
