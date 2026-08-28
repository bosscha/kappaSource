#set page(
  paper: "a4",
  margin: (x: 2cm, top: 2.2cm, bottom: 2.2cm),
  header: align(right)[
    #text(size: 8.5pt, fill: luma(120))[
      #smallcaps[The $kappa$-Source Framework] --- Theory, Generation & Extraction
    ]
  ],
  footer: context {
    let page_number = counter(page).get().first()
    let total_pages = counter(page).final().first()
    align(center)[
      #text(size: 8.5pt, fill: luma(120))[
        Page #page_number of #total_pages
      ]
    ]
  }
)

#set text(
  font: "Liberation Serif",
  size: 10.5pt,
  lang: "en"
)

#set heading(numbering: "1.1")
#show heading: it => [
  #v(0.6em)
  #text(fill: rgb("#1a365d"), weight: "bold")[#it]
  #v(0.3em)
]

#show raw.where(block: true): block.with(
  fill: rgb("#f8fafc"),
  inset: 10pt,
  radius: 4pt,
  stroke: 1pt + rgb("#e2e8f0"),
)

#show raw.where(block: false): box.with(
  fill: rgb("#f1f5f9"),
  inset: (x: 3pt, y: 1.5pt),
  radius: 2pt,
)

// --- TITLE BLOCK ---
#align(center)[
  #v(1em)
  #text(size: 20pt, weight: "bold", fill: rgb("#0f172a"))[
    The $kappa$-Source Framework
  ] \
  #v(0.3em)
  #text(size: 13pt, weight: "medium", fill: rgb("#334155"))[
    Theoretical Formulation, Mock Image Generation, and Hierarchical Extraction
  ] \
  #v(0.6em)
  #text(size: 9.5pt, fill: rgb("#64748b"))[
    Standalone High-Performance Rust Implementation | Software Version 0.1.0
  ]
  #v(1em)
  #line(length: 100%, stroke: 1.5pt + rgb("#3b82f6"))
]

#v(0.5em)

// --- ABSTRACT ---
#rect(
  fill: rgb("#f0f7ff"),
  stroke: 1pt + rgb("#bfdbfe"),
  radius: 6pt,
  inset: 12pt,
  width: 100%,
)[
  #text(weight: "bold", fill: rgb("#1e40af"))[Abstract] \
  #v(0.2em)
  This document provides the formal theoretical definition and practical operating instructions for the *$kappa$-Source Framework*. Under this formalism, any astronomical image is conceptualized as the union of $kappa$-sources ($kappa in {1, 2, dots, n}$). We define a $kappa$-source as a coherent cluster of $kappa$ constituent point sources within a spatial radius $R <= R_max$ whose individual flux is below standard detection thresholds ($F_i < 3 "RMS"$), but whose collective integrated flux satisfies $sum F_i >= 3 "RMS"$. We present high-performance Rust tools (`kappa_generate.bin` and `kappa_extract.bin`) to synthesize mock FITS images and perform hierarchical extraction from arbitrary FITS files.
]

#v(1em)

= Theoretical Concept of $kappa$-Sources

== The Fundamental Image Decomposition Principle
In conventional astronomical source detection (e.g., standard thresholding), images are decomposed solely into isolated point-like or extended components exceeding a local detection threshold (typically $3 sigma$ or $5 sigma$). Faint sub-threshold components buried in noise are discarded.

Under the *$kappa$-Source Formalism*, an astronomical image $I(x, y)$ is defined as the union of $kappa$-sources:
$ I(x, y) = cal(N)(0, sigma_text("noise")^2) + union.big_(kappa=1)^n union.big_(j=1)^(N_kappa) S_(kappa, j)(x, y) $
where $kappa$ designates the *multiplicity* (number of constituent point sources) of the source, and $N_kappa$ is the number of $kappa$-sources of order $kappa$.

== Mathematical Definition of a $kappa$-Source
A $kappa$-Source $C_kappa = {s_1, s_2, dots, s_kappa}$ of multiplicity $kappa$ is formed by $kappa$ point sources satisfying two strict criteria:

+ *Spatial Coherence Constraint*: All $kappa$ constituent point sources reside within a maximum clustering radius $R_max$ from their flux-weighted centroid:
  $ max_(i in {1..kappa}) norm(bold(x)_i - bold(x)_text("cen")) <= R_max $

+ *Collective Significance & Sub-threshold Condition*:
  - For $kappa >= 2$ (Multi-component $kappa$-sources): Each constituent subcomponent $i$ is individually below the point-source detection threshold:
    $ forall i in {1..kappa}, quad F_i < 3.0 times sigma_text("RMS") $
    while their collective integrated flux satisfies the detection threshold:
    $ F_text("total") = sum_(i=1)^kappa F_i >= 3.0 times sigma_text("RMS") $
  - For $kappa = 1$ (1-Sources): The single isolated source meets detection on its own:
    $ F_1 >= 3.0 times sigma_text("RMS") $

#table(
  columns: (1.8fr, 2.5fr, 2.5fr, 3.2fr),
  fill: (col, row) => if row == 0 { rgb("#e2e8f0") } else if calc.even(row) { rgb("#f8fafc") } else { white },
  stroke: 0.5pt + rgb("#cbd5e1"),
  align: (center, center, center, left),
  [Multiplicity $kappa$], [Individual Flux $F_i$], [Collective Flux $sum F_i$], [Observational Classification],
  [$kappa = 1$ (1-source)], [$F_1 >= 3 sigma_text("RMS")$], [$F_1 >= 3 sigma_text("RMS")$], [Directly detected point source],
  [$kappa = 2$ (2-source)], [$F_i < 3 sigma_text("RMS")$], [$sum F_i >= 3 sigma_text("RMS")$], [Pair of faint subcomponents],
  [$kappa = 3$ (3-source)], [$F_i < 3 sigma_text("RMS")$], [$sum F_i >= 3 sigma_text("RMS")$], [Triplet of faint subcomponents],
  [$kappa = n$ ($n$-source)], [$F_i < 3 sigma_text("RMS")$], [$sum F_i >= 3 sigma_text("RMS")$], [Cluster of $n$ faint subcomponents],
)

== Hierarchical Extraction Sequence
Extraction proceeds *hierarchically in order of increasing multiplicity*:
1. First, all *1-sources* ($kappa = 1$) are identified and cataloged.
2. Next, all *2-sources* ($kappa = 2$) are identified and cataloged.
3. Next, all *3-sources* ($kappa = 3$), up to $kappa = n$.

For every extracted $kappa$-source, the following physical properties are measured:
- *Flux-weighted Centroid*:
  $ bold(x)_text("cen") = (macron(x), macron(y)) = ( (sum_(i=1)^kappa F_i x_i) / (sum_(i=1)^kappa F_i), (sum_(i=1)^kappa F_i y_i) / (sum_(i=1)^kappa F_i) ) $
- *Characteristic Spatial Extent (Radius)*:
  $ R = max_(i in {1..kappa}) sqrt((x_i - macron(x))^2 + (y_i - macron(y))^2) + "FWHM" / 2 $
- *Signal-to-Noise Ratio (SNR)*:
  $ "SNR" = F_text("total") / sigma_text("RMS") $

#v(1em)

= Multi-Extension FITS Architecture

All data generated and extracted are stored in standard multi-extension astronomical FITS files:

#table(
  columns: (1fr, 2fr, 2fr, 5fr),
  fill: (col, row) => if row == 0 { rgb("#e2e8f0") } else if calc.even(row) { rgb("#f8fafc") } else { white },
  stroke: 0.5pt + rgb("#cbd5e1"),
  align: (center, center, center, left),
  [HDU], [Name], [Type], [Contents],
  [0], [`PRIMARY`], [Image HDU], [2D Float array ($M times M$) with metadata (`NKAPPA`, `KAP_MAX`, `KAP_RAD`, `DET_SIG`, `BG_SIGMA`)],
  [1], [`SOURCES`], [Binary Table], [Point sources catalog (`ID`, `X`, `Y`, `FLUX`, `AMPLITUDE`, `SIGMA`, `FWHM`, `KAPPA_ID`, `KAPPA`)],
  [2], [`KAPPA_SRCS`], [Binary Table], [Extracted $kappa$-sources catalog (`KAPPA_ID`, `KAPPA`, `CEN_X`, `CEN_Y`, `TOTAL_FLUX`, `RADIUS`, `SNR`, `N_MEMBERS`)],
)

#v(1em)

= Practical Guide: Generating Mock FITS Images

The tool `kappa_generate.bin` synthesizes an $M times M$ mock FITS image containing $N$ $kappa$-sources with multiplicity $kappa <= kappa_max$ and cluster radius $R <= R_max$ on top of Gaussian background noise.

== Command-Line Parameters

#table(
  columns: (2.2fr, 1fr, 1.2fr, 5.6fr),
  fill: (col, row) => if row == 0 { rgb("#e2e8f0") } else if calc.even(row) { rgb("#f8fafc") } else { white },
  stroke: 0.5pt + rgb("#cbd5e1"),
  align: (left, center, center, left),
  [Option], [Short], [Default], [Description],
  [`--num-kappa`], [`-N`], [`50`], [Total number of $kappa$-sources to produce ($N$)],
  [`--size`], [`-M`], [`4096`], [Image grid dimension $M$ ($M times M$ pixels)],
  [`--max-kappa`], [`-k`], [`5`], [Maximum multiplicity bound ($1 <= kappa <= kappa_max$)],
  [`--max-radius`], [`-r`], [`25.0`], [Maximum spatial cluster radius in pixels ($R_max$)],
  [`--detection-sigma`], [`-s`], [`3.0`], [Collective flux threshold ($sum F_i >= S times "RMS"$)],
  [`--noise-sigma`], [], [`1.0`], [Background Gaussian noise standard deviation (RMS)],
  [`--fwhm`], [], [`10.0`], [Gaussian PSF Full Width at Half Maximum in pixels],
  [`--output`], [`-o`], [`mock_kappa_image.fits`], [Output FITS file path],
)

== Usage Examples

```bash
# 1. Generate 50 kappa-sources (kappa <= 5, radius <= 25 px) on a 4096 x 4096 grid:
./kappa_generate.bin -N 50 -M 4096 -k 5 -r 25.0 -s 3.0 --output mock_kappa.fits

# 2. Generate 100 kappa-sources with custom noise RMS = 2.0 and detection = 5xRMS:
./kappa_generate.bin -N 100 -M 4096 -k 6 -r 30.0 --noise-sigma 2.0 -s 5.0 -o sim.fits
```

#v(1em)

= Practical Guide: Extracting $kappa$-Sources from FITS

The tool `kappa_extract.bin` ingests any arbitrary 2D FITS file (real or simulated) and extracts the $kappa$-source hierarchy.

== Extraction Algorithm Pipeline
1. *Noise Background Ingestion*: Computes image median and Median Absolute Deviation ($sigma = 1.4826 times "MAD"$).
2. *Gaussian Matched Filtering*: Convolves the 2D image with the PSF kernel to maximize the SNR of point-source features and eliminate single-pixel noise fluctuations.
3. *Subcomponent Peak Finding*: Detects local intensity peaks with sub-pixel quadratic centroid refinement and aperture flux integration.
4. *Hierarchical Graph Clustering*: Groups subcomponents within distance $R_max$, validates $sum F_i >= 3 "RMS"$, verifies individual subcomponent limits, and sorts from $kappa=1$ upwards.
5. *Catalog Export*: Produces `<input>.extracted.fits` and `<input>.extracted.reg` (DS9 region overlay).

== Command-Line Parameters

#table(
  columns: (2.4fr, 1fr, 1.2fr, 5.4fr),
  fill: (col, row) => if row == 0 { rgb("#e2e8f0") } else if calc.even(row) { rgb("#f8fafc") } else { white },
  stroke: 0.5pt + rgb("#cbd5e1"),
  align: (left, center, center, left),
  [Option], [Short], [Default], [Description],
  [`<INPUT>`], [], [*(required)*], [Path to input FITS image (e.g. `test.fits`)],
  [`--max-kappa`], [`-k`], [`0` *(all)*], [Maximum multiplicity upper limit (e.g. `-k 3` restricts to $kappa <= 3$)],
  [`--detection-sigma`], [`-s`], [`3.0`], [Collective flux detection threshold ($>= S times "RMS"$)],
  [`--cluster-radius`], [`-r`], [`25.0`], [Maximum clustering radius in pixels from centroid],
  [`--fwhm`], [], [`10.0`], [Estimated PSF FWHM in pixels],
  [`--peak-snr`], [], [`2.5`], [Matched-filter peak SNR threshold for candidate peaks],
  [`--output`], [`-o`], [`<in>.extracted.fits`], [Output FITS catalog path],
)

== Usage Examples

```bash
# 1. Standard extraction with 3xRMS threshold and 25 px radius:
./kappa_extract.bin test.fits -s 3.0 -r 25.0

# 2. Extract only up to 3-sources (kappa <= 3) with 5xRMS detection threshold:
./kappa_extract.bin test.fits -k 3 -s 5.0 -r 30.0 -o catalog.fits
```

#v(1em)

= Visualization with SAOImage DS9 & Python

== DS9 Color-Coded Overlay
Both `kappa_generate.bin` and `kappa_extract.bin` automatically export companion `.reg` files color-coded by multiplicity hierarchy:
- 🟢 *Green circles*: 1-Sources ($kappa = 1$)
- 🟡 *Yellow circles*: 2-Sources ($kappa = 2$)
- 🟠 *Orange circles*: 3-Sources ($kappa = 3$)
- 🔴 *Red circles*: $kappa >= 4$-Sources
- ⚪ *Dashed white circles*: Constituent subcomponent peaks

Open the image and region overlay with one command:
```bash
ds9 test.fits -regions test.extracted.reg -scale zscale -zoom to fit
```

== Inspecting Catalogs in Python (Astropy)
```python
from astropy.io import fits
import numpy as np

with fits.open("test.extracted.fits") as hdul:
    hdul.info()
    
    # 1. Access extracted kappa-sources catalog
    ktable = hdul["KAPPA_SRCS"].data
    print(f"Total kappa-sources extracted: {len(ktable)}")
    for row in ktable[:5]:
        print(f"ID={row['KAPPA_ID']}: kappa={row['KAPPA']}, Pos=({row['CEN_X']:.1f}, {row['CEN_Y']:.1f}), Flux={row['TOTAL_FLUX']:.3f}, SNR={row['SNR']:.1f}")
```

#v(1.5em)
#line(length: 100%, stroke: 0.5pt + rgb("#cbd5e1"))
#align(center)[
  #text(size: 8pt, fill: rgb("#94a3b8"))[
    Generated by the $kappa$-Source Analysis Toolchain | Rust High-Performance Astronomical Software
  ]
]
