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
  This document provides the formal theoretical definition and practical operating instructions for the *$kappa$-Source Framework*. Under this formalism, any astronomical image is conceptualized as the union of $kappa$-sources ($kappa in {1, 2, dots, n}$). We define a $kappa$-source as a coherent cluster of $kappa$ constituent point sources within a spatial radius $R <= R_max$ whose individual flux is below standard detection thresholds ($F_i < 3 "RMS"$), but whose collective integrated flux satisfies $sum F_i >= 3 "RMS"$. In mock image synthesis, each point source subcomponent is explicitly convolved with the telescope Point Spread Function (PSF). We present high-performance Rust tools (`kappa_generate.bin` and `kappa_extract.bin`) to synthesize mock FITS images and perform hierarchical extraction from arbitrary FITS files.
]

#v(1em)

= Theoretical Concept of $kappa$-Sources

== The Fundamental Image Decomposition Principle
In conventional astronomical source detection (e.g., standard thresholding), images are decomposed solely into isolated point-like or extended components exceeding a local detection threshold (typically $3 sigma$ or $5 sigma$). Faint sub-threshold components buried in noise are discarded.

Under the *$kappa$-Source Formalism*, an astronomical sky model $I_text("sky")(x, y)$ consists of constituent subcomponents that are convolved with the instrument Point Spread Function $"PSF"(x, y)$ and perturbed by background noise:
$ I_text("obs")(x, y) = [ I_text("sky") star "PSF" ](x, y) + cal(N)(0, sigma_text("RMS")^2) $
where the sky is formed by the union of $kappa$-sources:
$ I_text("sky")(x, y) = union.big_(kappa=1)^n union.big_(j=1)^(N_kappa) S_(kappa, j)(x, y) $

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

== Point Spread Function (PSF) Convolution
Each constituent subcomponent $i$ is a point source with flux $F_i$ at position $(x_i, y_i)$:
$ I_(i, text("sky"))(x, y) = F_i delta(x - x_i, y - y_i) $
Convolution with the telescope PSF produces the observed 2D intensity profile:
$ I_(i, text("conv"))(x, y) = [ I_(i, text("sky")) star "PSF" ](x, y) = F_i dot "PSF"(x - x_i, y - y_i) $

Supported PSF models:
- *2D Gaussian PSF*:
  $ "PSF"(r) = 1 / (2 pi sigma_text("PSF")^2) exp(- r^2 / (2 sigma_text("PSF")^2) ), quad sigma_text("PSF") = "FWHM" / (2 sqrt(2 ln 2)) $
- *2D Moffat PSF* (Atmospheric Seeing):
  $ "PSF"(r) = (beta - 1) / (pi alpha^2) [ 1 + r^2 / alpha^2 ]^(-beta), quad alpha = "FWHM" / (2 sqrt(2^(1/beta) - 1)) $

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
  [0], [`PRIMARY`], [Image HDU], [2D Float array ($M times M$) with metadata (`CONVOLV`, `PSF_TYPE`, `PSF_FWHM`, `NKAPPA`, `KAP_MAX`, `KAP_RAD`, `DET_SIG`, `BG_SIGMA`)],
  [1], [`SOURCES`], [Binary Table], [Point sources catalog (`ID`, `X`, `Y`, `FLUX`, `AMPLITUDE`, `SIGMA`, `FWHM`, `KAPPA_ID`, `KAPPA`)],
  [2], [`KAPPA_SRCS`], [Binary Table], [Extracted $kappa$-sources catalog (`KAPPA_ID`, `KAPPA`, `CEN_X`, `CEN_Y`, `TOTAL_FLUX`, `RADIUS`, `SNR`, `N_MEMBERS`)],
)

#v(1em)

= Practical Guide: Generating Mock FITS Images

The tool `kappa_generate.bin` synthesizes an $M times M$ mock FITS image containing $N$ $kappa$-sources with multiplicity $kappa <= kappa_max$ and cluster radius $R <= R_max$ convolved with the telescope PSF on top of Gaussian background noise.

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
  [`--psf`], [], [`gaussian`], [PSF profile model: `gaussian` or `moffat`],
  [`--fwhm`], [], [`10.0`], [PSF Full Width at Half Maximum in pixels],
  [`--moffat-beta`], [], [`4.765`], [Power-law index beta for Moffat PSF],
  [`--detection-sigma`], [`-s`], [`3.0`], [Collective flux threshold ($sum F_i >= S times "RMS"$)],
  [`--noise-sigma`], [], [`1.0`], [Background Gaussian noise standard deviation (RMS)],
  [`--output`], [`-o`], [`mock_kappa_image.fits`], [Output FITS file path],
)

== Usage Examples

```bash
# 1. Generate 50 kappa-sources convolved by Gaussian PSF on a 4096 x 4096 grid:
./kappa_generate.bin -N 50 -M 4096 -k 5 -r 25.0 --psf gaussian --fwhm 10.0 --output mock_kappa.fits

# 2. Generate with Moffat atmospheric seeing PSF (FWHM = 8 px, beta = 4.765):
./kappa_generate.bin -N 60 -M 2048 -k 4 -r 20.0 --psf moffat --fwhm 8.0 -s 3.0 -o moffat_sim.fits
```

#v(1em)

= Practical Guide: Extracting $kappa$-Sources from FITS

The tool `kappa_extract.bin` ingests any arbitrary 2D FITS file (real or simulated) and extracts the $kappa$-source hierarchy.

== Five-Stage Extraction Strategy & Algorithm Pipeline

Extracting $kappa$-sources from a noisy astronomical raster requires solving a dual problem: capturing faint sub-threshold peaks without triggering noise percolation, and testing whether spatial clusters reach statistical detection ($>= 3 "RMS"$). The complete strategy comprises five sequential stages:

+ *Noise Background Characterization ($sigma_text("RMS")$)*:
  Computes background median $mu$ and robust dispersion via the Median Absolute Deviation:
  $ sigma_text("RMS") = 1.4826 times "median"(|I(x, y) - mu|) $

+ *PSF Matched Filtering (Spatial Optimal Filtering)*:
  Convolves the background-subtracted raster with the PSF kernel $I_text("filt") = (I - mu) star "PSF"$. By the Matched Filter Theorem, this enhances real point-source peaks by a factor of $approx sqrt(2 pi sigma_text("PSF")^2) approx 10 times$ while suppressing single-pixel high-frequency noise spikes.

+ *Subcomponent Peak Finding & Aperture Photometry*:
  Detects local maxima in the matched-filter map down to low candidate SNR ($>= 2.5$). Sub-pixel centroids $(x_i, y_i)$ are refined via 2D intensity moments, and total flux $F_i$ is integrated using circular aperture photometry with analytic PSF aperture loss correction.

+ *Proximity Graph Clustering*:
  Constructs an adjacency graph where vertices are candidate peaks and edges connect neighbors with Euclidean separation $d(s_i, s_j) < 1.5 R_max$. Connected components of size $kappa$ are extracted, and their cluster radius $R_text("cluster") <= R_max$ from the centroid is verified.

+ *Physical Validation & Hierarchical Sorting*:
  Tests each cluster against physical significance criteria:
  - *For $kappa = 1$*: $F_1 >= S_text("det") times sigma_text("RMS")$
  - *For $kappa >= 2$*: $(forall i, F_i < S_text("det") times sigma_text("RMS")) "AND" (sum_(i=1)^kappa F_i >= S_text("det") times sigma_text("RMS"))$
  Catalog records are sorted hierarchically starting with 1-sources, then 2-sources, 3-sources, ..., $n$-sources.

== Extraction Strategy Flowchart

```text
┌─────────────────────────────────────────────────────────────────────────┐
│                           Raw 2D FITS Image                             │
└────────────────────────────────────┬────────────────────────────────────┘
                                     │
                                     ▼
┌─────────────────────────────────────────────────────────────────────────┐
│ 1. Background Estimation: Median μ  &  Noise σ_RMS = 1.4826 × MAD       │
└────────────────────────────────────┬────────────────────────────────────┘
                                     │
                                     ▼
┌─────────────────────────────────────────────────────────────────────────┐
│ 2. PSF Matched Filter: I_filt = [ (I - μ) ★ PSF ]  (Optimal SNR Boost)  │
└────────────────────────────────────┬────────────────────────────────────┘
                                     │
                                     ▼
┌─────────────────────────────────────────────────────────────────────────┐
│ 3. Subcomponent Photometry: Peak Finding + Sub-Pixel Centroid + Flux F_i│
└────────────────────────────────────┬────────────────────────────────────┘
                                     │
                                     ▼
┌─────────────────────────────────────────────────────────────────────────┐
│ 4. Proximity Graph Clustering: Connected Components with R ≤ R_max      │
└────────────────────────────────────┬────────────────────────────────────┘
                                     │
                                     ▼
┌─────────────────────────────────────────────────────────────────────────┐
│ 5. Physical Validation:                                                 │
│    • For κ = 1: Flux ≥ 3×RMS (Point Source Detection)                   │
│    • For κ ≥ 2: (Each F_i < 3×RMS) AND (∑ F_i ≥ 3×RMS) (Coherent Sum)   │
└────────────────────────────────────┬────────────────────────────────────┘
                                     │
                                     ▼
┌─────────────────────────────────────────────────────────────────────────┐
│ 6. Hierarchical Catalog Output: 1-src → 2-src → 3-src → ... → n-src     │
│    Multi-Extension FITS Table (.extracted.fits) + DS9 Overlay (.reg)    │
└─────────────────────────────────────────────────────────────────────────┘
```

== Output Timestamped Catalog Files

To prevent accidental overwrites and keep session records organized, `kappa_extract.bin` automatically writes timestamped catalog products using the format `<fitsname>_<YYYYMMDD_HHMMSS>`:

- *FITS Catalog* (`<fitsname>_<timestamp>.extracted.fits`): Binary table HDU with full 32-bit float columns (`KAPPA_ID`, `KAPPA`, `CEN_X`, `CEN_Y`, `TOTAL_FLUX`, `MAX_AMP`, `RADIUS`, `SNR`, `N_MEMBERS`).
- *ASCII Table* (`<fitsname>_<timestamp>.extracted.cat`): Plain text table with metadata header cards and aligned columns.
- *CSV Catalog* (`<fitsname>_<timestamp>.extracted.csv`): Comma-separated format for direct loading into Pandas, Astropy, or Excel.
- *DS9 Region File* (`<fitsname>_<timestamp>.extracted.reg`): Color-coded region overlay.

== Command-Line Parameters

#table(
  columns: (2.4fr, 1fr, 1.2fr, 5.4fr),
  fill: (col, row) => if row == 0 { rgb("#e2e8f0") } else if calc.even(row) { rgb("#f8fafc") } else { white },
  stroke: 0.5pt + rgb("#cbd5e1"),
  align: (left, center, center, left),
  [Option], [Short], [Default], [Description],
  [`<INPUT>`], [], [*(required)*], [Path to input FITS image (e.g. `test.fits`)],
  [`--max-kappa`], [`-k`], [`0` *(all)*], [Maximum multiplicity upper limit (e.g. `-k 3` restricts to $kappa <= 3$)],
  [`--detection-sigma`], [`-s`], [`3.0`], [Collective flux detection threshold ($>= S times "Beam RMS"$)],
  [`--search-radius`], [`-r`], [`25.0`], [Search radius $R_text("search")$ in pixels from centroid],
  [`--fwhm`], [], [`10.0`], [Estimated PSF FWHM in pixels],
  [`--min-sub-snr`], [], [`1.2`], [Minimum candidate peak SNR for subcomponents in matched filter],
  [`--seed-snr`], [], [`2.2`], [Candidate cluster seed threshold on smoothed map (alias: `--peak-snr`)],
  [`--output`], [`-o`], [`<auto>`], [Custom output FITS catalog path (defaults to timestamped name)],
)

== Usage Examples

```bash
# 1. Standard extraction generating timestamped FITS, ASCII .cat, CSV, and .reg files:
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
