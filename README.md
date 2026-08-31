# The $\kappa$-Source Framework

[![Rust](https://img.shields.io/badge/rust-1.75%2B-orange.svg)](https://www.rust-lang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A high-performance astronomical framework and Rust toolchain for the synthesis, simulation, and hierarchical extraction of **$\kappa$-Sources** in astronomical images.

---

## 🔭 Theoretical Concept

In conventional astronomical source detection, faint sub-threshold components buried under the background noise are lost when using standard thresholding. 

Under the **$\kappa$-Source Formalism**, an astronomical image is conceptualized as the union of $\kappa$-sources ($\kappa = 1, 2, \dots, n$):

$$I(x, y) = \mathcal{N}(0, \sigma_{\text{RMS}}^2) + \bigcup_{\kappa=1}^n \bigcup_{j=1}^{N_\kappa} S_{\kappa, j}(x, y)$$

### Physical Definition
A **$\kappa$-Source** $C_\kappa = \{s_1, s_2, \dots, s_\kappa\}$ is a coherent spatial cluster of $\kappa$ constituent point sources where:

1. **Spatial Coherence**: All $\kappa$ constituent point sources lie within a maximum search radius $R \le R_{\text{search}}$ from their flux-weighted centroid:
   $$\max_{i \in \{1..\kappa\}} \|\mathbf{x}_i - \mathbf{x}_{\text{cen}}\| \le R_{\text{search}}$$

2. **Collective Detection vs Sub-Threshold Condition**:
   - **For $\kappa \ge 2$ (Multi-component $\kappa$-sources)**: Every subcomponent is individually undetected ($F_i < 3.0 \times \sigma_{\text{beam}}$), while their collective integrated flux satisfies the detection threshold:
     $$\sum_{i=1}^\kappa F_i \ge 3.0 \times \sigma_{\text{beam}}$$
   - **For $\kappa = 1$ (1-Sources)**: The single isolated point source meets detection on its own ($F_1 \ge 3.0 \times \sigma_{\text{beam}}$).

---

## 🚀 Features

- **Two-Pass Hierarchical Extraction**:
  - *Pass 1*: Detects and locks genuine 1-sources ($F \ge 3\times\text{RMS}$) to prevent contamination.
  - *Pass 2*: Clusters sub-threshold peaks ($F_i < 3\times\text{RMS}$) within search radius $R_{\text{search}}$ satisfying collective significance.
- **Empirical Auto-PSF Calibration (`--psf-auto`)**: Automatically measures the seeing FWHM from the top 20 brightest stars / $\text{SNR} > 20$ via 2D spatial intensity moments.
- **Session-Timestamped Catalog Outputs**:
  - `<fitsname>_<timestamp>.extracted.fits`: Multi-extension FITS table (`KAPPA_SRCS` and `SOURCES`).
  - `<fitsname>_<timestamp>.extracted.cat` & `.csv`: Extracted $\kappa$-sources table.
  - `<fitsname>_<timestamp>.extracted.subcomponents.cat` & `.csv`: Member subcomponents table (unassociated noise peaks $\kappa=0$ are strictly excluded).
  - `<fitsname>_<timestamp>.extracted.reg`: DS9 color-coded region overlay.
- **Benchmark Cross-Matching Tool (`compare_kappa.py`)**: Computes completeness, positional accuracy, flux conservation ratios, and 4-panel diagnostic reports.
- **Complete PDF Documentation Guide**: Compiled with Typst (`kappa_source_guide.pdf`).

---

## 📦 Installation & Build

```bash
git clone https://github.com/bosscha/kappaSource.git
cd kappaSource
cargo build --release
```

Binaries:
- `target/release/kappa_generate` (`./kappa_generate.bin`)
- `target/release/kappa_extract` (`./kappa_extract.bin`)

---

## 🛠️ Usage

### 1. Generating Mock FITS Images (`kappa_generate.bin`)

```bash
# Generate 50 kappa-sources (1 <= kappa <= 4, radius <= 25 px) convolved with 8.7 px Gaussian PSF:
./kappa_generate.bin -N 50 -M 2048 -k 4 -r 25.0 -s 3.0 --psf gaussian --fwhm 8.7 --output mock_image.fits
```

| Option | Short | Default | Description |
| :--- | :--- | :--- | :--- |
| `--num-kappa` | `-N` | `50` | Number of $\kappa$-sources to produce ($N$) |
| `--size` | `-M` | `2048` | Image dimension $M$ for $M \times M$ grid |
| `--max-kappa` | `-k` | `4` | Maximum multiplicity upper limit ($\kappa \in [1, \kappa_{\max}]$) |
| `--max-radius` | `-r` | `25.0` | Maximum spatial cluster radius in pixels ($R_{\max}$) |
| `--psf` | | `gaussian` | PSF convolution model: `gaussian` or `moffat` |
| `--fwhm` | | `10.0` | PSF Full Width at Half Maximum in pixels |
| `--detection-sigma` | `-s` | `3.0` | Collective detection threshold ($\sum F_i \ge S \times \text{RMS}$) |
| `--noise-sigma` | | `1.0` | Background Gaussian noise RMS |
| `--output` | `-o` | `mock_kappa_image.fits` | Output FITS file path |

---

### 2. Extracting $\kappa$-Sources (`kappa_extract.bin`)

```bash
# Auto-PSF mode (measures seeing directly from bright stars):
./kappa_extract.bin mock_image.fits --psf-auto -s 3.0 -r 25.0

# Known PSF mode:
./kappa_extract.bin mock_image.fits --fwhm 8.7 -s 3.0 -r 25.0
```

| Option | Short | Default | Description |
| :--- | :--- | :--- | :--- |
| `<INPUT>` | | *(required)* | Path to input FITS image (e.g. `test.fits`) |
| `--max-kappa` | `-k` | `0` *(all)* | Upper limit on multiplicity $\kappa$ ($0$ for all) |
| `--detection-sigma` | `-s` | `3.0` | Total collective flux detection threshold ($\ge S \times \text{RMS}$) |
| `--search-radius` | `-r` | `25.0` | Search radius $R_{\text{search}}$ in pixels from centroid |
| `--fwhm` | | `10.0` | Estimated PSF FWHM in pixels (set 0 or use `--psf-auto`) |
| `--psf-auto` | | `false` | Automatically measure PSF FWHM from bright stars |
| `--min-psf-snr` | | `20.0` | Minimum SNR threshold for auto-PSF calibration stars |
| `--psf-samples` | | `20` | Maximum number of brightest sources to sample for auto-PSF |
| `--min-sub-snr` | | `1.2` | Minimum peak SNR for subcomponents in matched filter |
| `--seed-snr` | | `2.2` | Candidate cluster seed threshold on smoothed map (alias: `--peak-snr`) |
| `--output` | `-o` | `<auto>` | Output catalog FITS path (defaults to timestamped name) |

---

### 3. Recovery Benchmark Cross-Matching (`compare_kappa.py`)

```bash
# Compare extracted catalog against ground truth:
python3 compare_kappa.py mock_image.fits
```

---

## 📊 Benchmark Results

Performance on a $2048 \times 2048$ image with 50 injected $\kappa$-sources convolved by an $8.7\text{ px}$ Gaussian PSF:

| Multiplicity | Injected | Recovered (Exact) | Recovered (Auto-PSF) | Mean Offset | Mean Flux Ratio |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **1-sources** | 10 | 8 (80.0%) | 8 (80.0%) | 2.68 px | 1.057× |
| **2-sources** | 19 | 16 (84.2%) | 15 (78.9%) | 8.18 px | 1.098× |
| **3-sources** | 10 | 6 (60.0%) | 6 (60.0%) | 9.13 px | 0.864× |
| **4-sources** | 11 | 7 (63.6%) | 7 (63.6%) | 10.24 px | 2.250× |
| **Total / Median** | **50** | **37 (74.0%)** | **36 (72.0%)** | **6.34 px** | **0.986×** |

---

## 🎨 Visualizing in SAOImage DS9

```bash
ds9 mock_image.fits -regions mock_image_*.extracted.reg -scale zscale -zoom to fit
```

- 🟢 **Green**: 1-sources ($\kappa = 1$)
- 🟡 **Yellow**: 2-sources ($\kappa = 2$)
- 🟠 **Orange**: 3-sources ($\kappa = 3$)
- 🔴 **Red**: $\kappa \ge 4$-sources
- ⚪ **Dashed White**: Constituent subcomponent peaks

---

## 📄 Documentation

For the complete theoretical derivation and mathematical formulation, see [**`kappa_source_guide.pdf`**](kappa_source_guide.pdf).

To recompile:
```bash
typst compile kappa_source_guide.typ kappa_source_guide.pdf
```

---

## 📜 License

MIT License. See [LICENSE](LICENSE) for details.
