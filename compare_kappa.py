#!/usr/bin/env python3
"""
Compare Ground-Truth Injected kappa-Sources with Extracted kappa-Sources.
"""

import os
import glob
import sys
import argparse
import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def find_latest_extracted_file(truth_fits_path):
    stem = os.path.splitext(os.path.basename(truth_fits_path))[0]
    parent = os.path.dirname(truth_fits_path) or "."
    candidates = glob.glob(os.path.join(parent, f"{stem}_*.extracted.fits"))
    if candidates:
        candidates.sort(key=os.path.getmtime, reverse=True)
        return candidates[0]
    generic = os.path.join(parent, f"{stem}.extracted.fits")
    if os.path.exists(generic):
        return generic
    return None

def compare_catalogs(truth_fits_path, extracted_fits_path=None, match_radius_px=25.0, plot_output="comparison_report.png"):
    if extracted_fits_path is None:
        extracted_fits_path = find_latest_extracted_file(truth_fits_path)
        if extracted_fits_path is None:
            print(f"Error: Could not find any extracted catalog matching {truth_fits_path}")
            sys.exit(1)

    print("=" * 80)
    print("KAPPA-SOURCE BENCHMARK: Injected Ground Truth vs. Extracted Catalog")
    print("=" * 80)
    print(f"Ground-Truth File : {truth_fits_path}")
    print(f"Extracted File    : {extracted_fits_path}")
    print(f"Match Radius Cut  : <= {match_radius_px:.1f} pixels")
    print("=" * 80)

    with fits.open(truth_fits_path) as hd_true:
        true_hdr = hd_true[0].header
        true_srcs = hd_true["SOURCES"].data
        true_ksrcs = hd_true["KAPPA_SRCS"].data

    with fits.open(extracted_fits_path) as hd_ext:
        ext_hdr = hd_ext[0].header
        ext_ksrcs = hd_ext["KAPPA_SRCS"].data

    n_true = len(true_ksrcs)
    n_ext = len(ext_ksrcs)
    print(f"Ground Truth kappa-Sources Injected : {n_true}")
    print(f"Total Extracted Candidates          : {n_ext}")
    print("-" * 80)

    # Cross-matching
    matched_true_idx = []
    matched_ext_idx = []
    distances = []
    flux_ratios = []
    kappa_true_list = []
    kappa_ext_list = []

    used_ext_indices = set()

    for t_i, t_row in enumerate(true_ksrcs):
        tx, ty = t_row["CEN_X"], t_row["CEN_Y"]
        tf = t_row["TOTAL_FLUX"]
        tk = t_row["KAPPA"]

        best_e_i = None
        min_dist = float("inf")

        for e_i, e_row in enumerate(ext_ksrcs):
            if e_i in used_ext_indices:
                continue
            ex, ey = e_row["CEN_X"], e_row["CEN_Y"]
            d = np.sqrt((tx - ex)**2 + (ty - ey)**2)
            if d <= match_radius_px and d < min_dist:
                min_dist = d
                best_e_i = e_i

        if best_e_i is not None:
            used_ext_indices.add(best_e_i)
            e_row = ext_ksrcs[best_e_i]
            matched_true_idx.append(t_i)
            matched_ext_idx.append(best_e_i)
            distances.append(min_dist)
            flux_ratios.append(e_row["TOTAL_FLUX"] / tf)
            kappa_true_list.append(tk)
            kappa_ext_list.append(e_row["KAPPA"])

    n_matched = len(matched_true_idx)
    recovery_rate = (n_matched / n_true * 100.0) if n_true > 0 else 0.0
    spurious_count = n_ext - len(used_ext_indices)

    print(f"Matched kappa-Sources Successfully  : {n_matched} / {n_true} ({recovery_rate:.1f}% Completeness)")
    print(f"Spurious / Background Noise Clusts  : {spurious_count}")
    if n_matched > 0:
        print(f"Median Centroid Offset (px)         : {np.median(distances):.2f} px (Mean: {np.mean(distances):.2f} px)")
        print(f"Median Flux Recovery Ratio (Ext/Tr) : {np.median(flux_ratios):.3f} (Mean: {np.mean(flux_ratios):.3f})")
    print("-" * 80)

    # Detailed breakdown by true kappa
    max_k = max(true_ksrcs["KAPPA"]) if n_true > 0 else 1
    print(f"{'Multiplicity':<12} {'Injected':<10} {'Recovered':<10} {'Recovery %':<12} {'Mean Offset (px)':<18} {'Mean Flux Ratio':<16}")
    print("-" * 80)

    for k in range(1, max_k + 1):
        t_mask = (true_ksrcs["KAPPA"] == k)
        k_injected = np.sum(t_mask)
        if k_injected == 0:
            continue

        # Matched for this k
        k_matched_mask = [kt == k for kt in kappa_true_list]
        k_recovered = np.sum(k_matched_mask)
        k_rec_rate = (k_recovered / k_injected * 100.0)

        if k_recovered > 0:
            k_dists = [distances[i] for i, m in enumerate(k_matched_mask) if m]
            k_frs = [flux_ratios[i] for i, m in enumerate(k_matched_mask) if m]
            mean_dist_str = f"{np.mean(k_dists):.2f} px"
            mean_fr_str = f"{np.mean(k_frs):.3f}"
        else:
            mean_dist_str = "N/A"
            mean_fr_str = "N/A"

        print(f"{k:<1}-sources{'':<4} {k_injected:<10} {k_recovered:<10} {k_rec_rate:<12.1f} {mean_dist_str:<18} {mean_fr_str:<16}")

    print("=" * 80)

    # Diagnostic Plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle("kappa-Source Framework: Ground Truth vs Extracted Evaluation", fontsize=14, fontweight="bold")

    # Plot 1: Flux correlation
    ax1 = axes[0, 0]
    if n_matched > 0:
        true_fluxes = [true_ksrcs[i]["TOTAL_FLUX"] for i in matched_true_idx]
        ext_fluxes = [ext_ksrcs[i]["TOTAL_FLUX"] for i in matched_ext_idx]
        sc = ax1.scatter(true_fluxes, ext_fluxes, c=kappa_true_list, cmap="viridis", s=50, alpha=0.8, edgecolors="k")
        f_min = min(min(true_fluxes), min(ext_fluxes)) * 0.8
        f_max = max(max(true_fluxes), max(ext_fluxes)) * 1.2
        ax1.plot([f_min, f_max], [f_min, f_max], "r--", label="1:1 Perfect Recovery")
        ax1.set_xlim(f_min, f_max)
        ax1.set_ylim(f_min, f_max)
        plt.colorbar(sc, ax=ax1, label="True kappa")
    ax1.set_xlabel("Ground Truth Flux")
    ax1.set_ylabel("Extracted Flux")
    ax1.set_title("Total Flux Recovery Correlation")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot 2: Centroid offsets
    ax2 = axes[0, 1]
    if n_matched > 0:
        ax2.hist(distances, bins=15, color="royalblue", edgecolor="black", alpha=0.7)
        ax2.axvline(np.median(distances), color="crimson", linestyle="--", label=f"Median Offset = {np.median(distances):.2f} px")
    ax2.set_xlabel("Centroid Offset (pixels)")
    ax2.set_ylabel("Count")
    ax2.set_title("Centroiding Positional Accuracy")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Plot 3: Multiplicity Confusion Matrix
    ax3 = axes[1, 0]
    if n_matched > 0:
        conf_matrix = np.zeros((max_k, max(max(kappa_ext_list), max_k)), dtype=int)
        for kt, ke in zip(kappa_true_list, kappa_ext_list):
            conf_matrix[kt - 1, ke - 1] += 1
        im = ax3.imshow(conf_matrix, cmap="Blues", origin="lower")
        ax3.set_xticks(range(conf_matrix.shape[1]))
        ax3.set_xticklabels([f"κ={i+1}" for i in range(conf_matrix.shape[1])])
        ax3.set_yticks(range(max_k))
        ax3.set_yticklabels([f"κ={i+1}" for i in range(max_k)])
        for i in range(max_k):
            for j in range(conf_matrix.shape[1]):
                val = conf_matrix[i, j]
                ax3.text(j, i, str(val), ha="center", va="center", color="white" if val > conf_matrix.max()/2 else "black", fontweight="bold")
        plt.colorbar(im, ax=ax3, label="Count")
    ax3.set_xlabel("Extracted kappa")
    ax3.set_ylabel("True kappa")
    ax3.set_title("Multiplicity Confusion Matrix (κ_true vs κ_ext)")

    # Plot 4: Recovery rate by kappa
    ax4 = axes[1, 1]
    k_vals = list(range(1, max_k + 1))
    rec_rates = []
    for k in k_vals:
        t_count = np.sum(true_ksrcs["KAPPA"] == k)
        r_count = sum(1 for kt in kappa_true_list if kt == k)
        rec_rates.append(r_count / t_count * 100.0 if t_count > 0 else 0)
    bars = ax4.bar([f"κ={k}" for k in k_vals], rec_rates, color="teal", edgecolor="black", alpha=0.7)
    for b, r in zip(bars, rec_rates):
        ax4.text(b.get_x() + b.get_width()/2, b.get_height() + 1, f"{r:.1f}%", ha="center", va="bottom", fontweight="bold")
    ax4.set_ylim(0, 115)
    ax4.set_ylabel("Completeness (%)")
    ax4.set_title("Completeness by Multiplicity kappa")
    ax4.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(plot_output, dpi=200)
    print(f"Diagnostic report plots saved to : {plot_output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare Injected vs Extracted kappa-Sources")
    parser.add_argument("truth", help="Path to ground truth mock FITS file (e.g. mock_image.fits)")
    parser.add_argument("extracted", nargs="?", default=None, help="Path to extracted catalog FITS file (default: auto-detect latest timestamped catalog)")
    parser.add_argument("-r", "--radius", type=float, default=25.0, help="Positional match radius in pixels (default: 25.0)")
    parser.add_argument("-p", "--plot", default="comparison_report.png", help="Output PNG path for diagnostic figures")
    args = parser.parse_args()

    compare_catalogs(args.truth, args.extracted, match_radius_px=args.radius, plot_output=args.plot)
