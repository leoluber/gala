
# GaLa Lineshape

A lightweight Python toolkit for quantitative analysis of semiconductor nanocrystal ensembles using photoluminescence spectroscopy.

The Gaussian-Asymmetric Laplace (GaLa) framework captures a common statistical structure in ensemble photoluminescence across a wide variety of semiconductor nanocrystals by describing emission as the convolution of a Gaussian contribution with independent low- and high-energy spectral tails.

GaLa transforms spectral lineshapes into quantitative information about the emissive ensemble. Its parameters capture distinct contributions to spectral broadening, can be used to investigate temperature-dependent emission, and, where an appropriate size–energy relationship is available, provide access to relative nanocrystal size dispersion. Systematic deviations between experimental spectra and the fitted lineshape can further reveal hidden emissive subpopulations.


The model is described in:

> L. Luber, D. Benedikt, A. Singldinger, A. S. Urban,  
> *Photoluminescence as a Quantitative Probe of Semiconductor Nanocrystals*.  
> [Nano Lett. 2026; 10.1021/acs.nanolett.6c04076](https://doi.org/10.1021/acs.nanolett.6c04076).

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

Experimental spectra can be provided as CSV or TXT files.

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
