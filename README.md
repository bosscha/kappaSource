# The $\kappa$-Source Framework

[![Rust](https://img.shields.io/badge/rust-1.75%2B-orange.svg)](https://www.rust-lang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A high-performance astronomical framework and Rust toolchain for the synthesis, simulation, and hierarchical extraction of **$\kappa$-Sources** in astronomical images.

---

## 🔭 Theoretical Concept

In astronomical source detection, faint sub-threshold components buried under the background noise are typically lost when using standard thresholding. 

Under the **$\kappa$-Source Formalism**, an astronomical image is conceptualized as the union of $\kappa$-sources ($\kappa = 1, 2, \dots, n$):

$$I(x, y) = \mathcal{N}(0, \sigma_{\text{RMS}}^2) + \bigcup_{\kappa=1}^n \bigcup_{j=1}^{N_\kappa} S_{\kappa, j}(x, y)$$

### Physical Definition
A **$\kappa$-Source** $C_\kappa = \{s_1, s_2, \dots, s_\kappa\}$ is a coherent spatial cluster of $\kappa$ constituent point sources where:

1. **Spatial Coherence**: All $\kappa$ constituent point sources lie within a maximum clustering radius $R \le R_{\max}$ from their flux-weighted centroid:
   $$\max_{i \in \{1..\kappa\}} \|\mathbf{x}_i - \mathbf{x}_{\text{cen}}\| \le R_{\max}$$

2. **Collective Detection vs Sub-Threshold Condition**:
   - **For $\kappa \ge 2$ (Multi-component $\kappa$-sources)**: Every subcomponent is individually undetected:
     $$\forall i \in \{1..\kappa\}, \quad F_i < 3.0 \times \sigma_{\text{RMS}}$$
     while their collective integrated emission is detectable:
     $$\sum_{i=1}^\kappa F_i \ge 3.0 \times \sigma_{\text{RMS}}$$
   - **For $\kappa = 1$ (1-Sources)**: The single isolated source meets detection on its own ($F_1 \ge 3.0 \times \sigma_{\text{RMS}}$).

---

## 🚀 Features

- **Hierarchical Extraction**: Sequentially extracts 1-sources, then 2-sources, 3-sources, ..., up to $n$-sources.
- **Pure Rust FITS I/O**: Direct read and write support for multi-extension FITS files without external C library dependencies.
- **Multithreaded Simulation**: Parallel Gaussian noise synthesis and local $5\sigma$ PSF rendering via Rayon.
- **Gaussian Matched Filtering**: 2D separable convolution for robust subcomponent peak detection in noisy backgrounds.
- **Multi-Extension FITS Output**:
  - `PRIMARY`: 2D Float image raster + comprehensive metadata header cards.
  - `SOURCES`: Catalog table of all constituent point sources and subcomponents.
  - `KAPPA_SRCS`: Catalog table of extracted $\kappa$-sources with centroids, total fluxes, radii, and SNR.
- **DS9 Color-Coded Overlay**: Exports companion `.reg` files color-coded by $\kappa$ (🟢 $\kappa=1$, 🟡 $\kappa=2$, 🟠 $\kappa=3$, 🔴 $\kappa \ge 4$).
- **Complete PDF Guide**: Theoretical paper and manual included (`kappa_source_guide.pdf`).

---

## 📦 Installation & Build

Ensure you have a Rust toolchain installed:

```bash
git clone https://github.com/bosscha/kappaSource.git
cd kappaSource
cargo build --release
```

This generates two release binaries:
- `target/release/kappa_generate`
- `target/release/kappa_extract`

---

## 🛠️ Usage

### 1. Generating Mock FITS Images (`kappa_generate`)

Synthesize an $M \times M$ FITS image containing $N$ $\kappa$-sources with $\kappa \le \kappa_{\max}$ and $R \le R_{\max}$:

```bash
# Generate 50 kappa-sources (kappa <= 5, radius <= 25 px) on a 4096 x 4096 grid:
cargo run --release --bin kappa_generate -- \
  -N 50 \
  -M 4096 \
  -k 5 \
  -r 25.0 \
  -s 3.0 \
  --output mock_image.fits
```

#### Generator Options
| Option | Short | Default | Description |
| :--- | :--- | :--- | :--- |
| `--num-kappa` | `-N` | `50` | Number of $\kappa$-sources to produce ($N$) |
| `--size` | `-M` | `4096` | Image dimension $M$ for $M \times M$ grid |
| `--max-kappa` | `-k` | `5` | Maximum multiplicity upper limit ($\kappa \in [1, \kappa_{\max}]$) |
| `--max-radius` | `-r` | `25.0` | Maximum spatial cluster radius in pixels ($R_{\max}$) |
| `--detection-sigma` | `-s` | `3.0` | Collective detection threshold ($\sum F_i \ge S \times \text{RMS}$) |
| `--noise-sigma` | | `1.0` | Background Gaussian noise RMS |
| `--fwhm` | | `10.0` | Gaussian PSF FWHM in pixels |
| `--output` | `-o` | `mock_kappa_image.fits` | Output FITS file path |

---

### 2. Extracting $\kappa$-Sources from a FITS Image (`kappa_extract`)

Ingest any 2D astronomical FITS file and extract the $\kappa$-source hierarchy:

```bash
# Extract kappa-sources with 3xRMS detection threshold and max radius 25 px:
cargo run --release --bin kappa_extract -- test.fits -s 3.0 -r 25.0

# Restrict extraction up to 3-sources (kappa <= 3):
cargo run --release --bin kappa_extract -- test.fits -k 3 -s 3.0 -r 25.0 -o catalog.fits
```

#### Extractor Options
| Option | Short | Default | Description |
| :--- | :--- | :--- | :--- |
| `<INPUT>` | | *(required)* | Path to input FITS image (e.g. `test.fits`) |
| `--max-kappa` | `-k` | `0` *(all)* | Upper limit on multiplicity $\kappa$ ($0$ for all) |
| `--detection-sigma` | `-s` | `3.0` | Total collective flux detection threshold ($\ge S \times \text{RMS}$) |
| `--cluster-radius` | `-r` | `25.0` | Maximum clustering radius in pixels from centroid |
| `--fwhm` | | `10.0` | Estimated PSF FWHM in pixels |
| `--peak-snr` | | `2.5` | Peak SNR threshold for candidate subcomponents |
| `--output` | `-o` | `<in>.extracted.fits` | Output catalog FITS path |

---

## 🎨 Visualizing in SAOImage DS9

Both tools automatically export `.reg` companion files color-coded by $\kappa$:

```bash
ds9 mock_image.fits -regions mock_image.reg -scale zscale -zoom to fit
```

- 🟢 **Green**: 1-sources ($\kappa = 1$)
- 🟡 **Yellow**: 2-sources ($\kappa = 2$)
- 🟠 **Orange**: 3-sources ($\kappa = 3$)
- 🔴 **Red**: $\kappa \ge 4$-sources
- ⚪ **Dashed White**: Individual subcomponent peaks

---

## 🐍 Python / Astropy Inspection

```python
from astropy.io import fits

with fits.open("mock_image.fits") as hdul:
    hdul.info()
    
    # Access extracted kappa-sources catalog table:
    kappa_catalog = hdul["KAPPA_SRCS"].data
    print(f"Total kappa-sources: {len(kappa_catalog)}")
    for row in kappa_catalog[:5]:
        print(f"ID={row['KAPPA_ID']}: kappa={row['KAPPA']}, "
              f"Centroid=({row['CEN_X']:.1f}, {row['CEN_Y']:.1f}), "
              f"Flux={row['TOTAL_FLUX']:.3f}, SNR={row['SNR']:.1f}")
```

---

## 📄 Documentation

For a detailed theoretical derivation and complete parameter formulas, see [**`kappa_source_guide.pdf`**](kappa_source_guide.pdf) (compiled with [Typst](https://typst.app/)).

```bash
# To recompile the PDF guide:
typst compile kappa_source_guide.typ kappa_source_guide.pdf
```

---

## 📜 License

MIT License. See [LICENSE](LICENSE) for details.
