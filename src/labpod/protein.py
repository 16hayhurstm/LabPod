from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


DEFAULT_BSA_CONCS = np.array([2000, 1000, 500, 250, 125, 62.5, 31.25, 0], dtype=float)
ROW_LABELS = list("ABCDEFGH")


@dataclass(frozen=True)
class ProteinQuantResult:
    r2: float
    slope: float
    intercept: float
    equation: str
    standards: pd.DataFrame
    samples: pd.DataFrame
    fig: plt.Figure
    csv_path: Optional[str]


def _find_plate_block(df: pd.DataFrame) -> Tuple[int, int]:
    for label_col in range(min(10, df.shape[1])):
        col = df.iloc[:, label_col].astype(str).str.strip().str.upper()
        for r in range(0, df.shape[0] - 7):
            if col.iloc[r : r + 8].tolist() == ROW_LABELS:
                return r, label_col
    raise ValueError("Could not find consecutive row labels A–H in the sheet.")


def _numeric_block(df: pd.DataFrame, start_row: int, label_col: int) -> pd.DataFrame:
    block = df.iloc[start_row : start_row + 8, label_col + 1 :].copy()
    block = block.apply(pd.to_numeric, errors="coerce")
    while block.shape[1] > 0 and block.iloc[:, -1].isna().all():
        block = block.iloc[:, :-1]
    return block


def _linear_fit(x: np.ndarray, y: np.ndarray) -> Tuple[float, float, float]:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if x.size < 2:
        raise ValueError("Not enough valid points to fit a line (need >=2).")

    slope, intercept = np.polyfit(x, y, 1)
    y_pred = slope * x + intercept
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1.0 - (ss_res / ss_tot if ss_tot != 0 else np.nan)
    return float(slope), float(intercept), float(r2)


def protein_quantifier(
    excel_path: Union[str, Path],
    n_protein_samples: int,
    *,
    sheet: Union[int, str] = 0,
    bsa_concs_ug_ml: Optional[np.ndarray] = None,
    dilution_factor: Union[float, np.ndarray] = 1.0,
    clip_negative: bool = True,
    autosave_csv: bool = True,
    csv_out: Optional[Union[str, Path]] = None,
    show_plot: bool = True,
    print_report: bool = True,
) -> ProteinQuantResult:
    """
    One-line user version:
      protein_quantifier("file.xlsx", n_samples)

    What it does automatically:
      - fits BSA curve (blank-corrected using tube H)
      - prints equation + R²
      - prints a clean concentrations table (final_conc_ug_ml)
      - saves CSV (default on)
      - shows plot (default on)

    Fixed layout:
      - Column 0: A..H labels
      - Next 3 columns: standards triplicates
      - Next column: gap (ignored)
      - Next 3 columns: protein replicates (cols 5–7)
      - Sample i uses row (A,B,C,...) across those 3 protein replicate columns
    """
    excel_path = Path(excel_path)
    if not excel_path.exists():
        raise FileNotFoundError(f"Excel file not found: {excel_path}")

    if not (1 <= n_protein_samples <= 8):
        raise ValueError("n_protein_samples must be between 1 and 8 (rows A..H).")

    concs = DEFAULT_BSA_CONCS if bsa_concs_ug_ml is None else np.asarray(bsa_concs_ug_ml, dtype=float)
    if concs.shape != (8,):
        raise ValueError("bsa_concs_ug_ml must be length 8 for tubes A..H.")

    df = pd.read_excel(excel_path, sheet_name=sheet, header=None, engine="openpyxl")
    start_row, label_col = _find_plate_block(df)
    block = _numeric_block(df, start_row, label_col)

    if block.shape[1] < 7:
        raise ValueError(
            f"Expected at least 7 columns after the A–H label column "
            f"(3 standards + gap + 3 protein replicate cols), but found {block.shape[1]}."
        )

    # ---- Standards ----
    std_reps = block.iloc[:, 0:3].to_numpy(dtype=float)  # (8,3)
    std_mean = np.nanmean(std_reps, axis=1)
    blank = float(std_mean[7])                           # tube H
    std_mean_bc = std_mean - blank

    slope, intercept, r2 = _linear_fit(std_mean_bc, concs)
    equation = f"y = {slope:.4g}x + {intercept:.4g}"

    standards_df = pd.DataFrame({
        "tube": ROW_LABELS,
        "bsa_conc_ug_ml": concs,
        "rep1": std_reps[:, 0],
        "rep2": std_reps[:, 1],
        "rep3": std_reps[:, 2],
        "mean_abs": std_mean,
        "blank_abs": blank,
        "mean_abs_minus_blank": std_mean_bc,
    })

    # ---- Protein samples: fixed replicate columns ----
    protein_rep_block = block.iloc[:, 4:7]              # after label: [0..2]=std, [3]=gap, [4..6]=protein reps
    samp_reps_all = protein_rep_block.to_numpy(dtype=float)  # (8,3)
    samp_reps = samp_reps_all[:n_protein_samples, :]         # (n_samples,3)

    samp_mean = np.nanmean(samp_reps, axis=1)
    samp_mean_bc = samp_mean - blank

    # Raw back-calc
    calc_conc_raw = slope * samp_mean_bc + intercept

    # QC flags
    below_blank_flag = samp_mean_bc < 0
    min_x = float(np.nanmin(std_mean_bc))
    max_x = float(np.nanmax(std_mean_bc))
    extrapolated_low_flag = samp_mean_bc < min_x
    extrapolated_high_flag = samp_mean_bc > max_x

    # Clip negatives (optional)
    if clip_negative:
        used_clipped_flag = calc_conc_raw < 0
        calc_conc = np.maximum(calc_conc_raw, 0.0)
    else:
        used_clipped_flag = np.array([False] * n_protein_samples)
        calc_conc = calc_conc_raw

    # Dilution
    if np.isscalar(dilution_factor):
        dil = np.full(n_protein_samples, float(dilution_factor))
    else:
        dil = np.asarray(dilution_factor, dtype=float)
        if dil.shape != (n_protein_samples,):
            raise ValueError("If dilution_factor is an array, it must be length n_protein_samples.")

    final_conc = calc_conc * dil

    samples_df = pd.DataFrame({
        "sample": [f"Sample_{i+1} (row {ROW_LABELS[i]})" for i in range(n_protein_samples)],
        "rep1": samp_reps[:, 0],
        "rep2": samp_reps[:, 1],
        "rep3": samp_reps[:, 2],
        "mean_abs": samp_mean,
        "mean_abs_minus_blank": samp_mean_bc,
        "calc_conc_ug_ml_raw": calc_conc_raw,
        "calc_conc_ug_ml_used": calc_conc,
        "dilution_factor": dil,
        "final_conc_ug_ml": final_conc,
        "below_blank_flag": below_blank_flag,
        "used_clipped_flag": used_clipped_flag,
        "extrapolated_low_flag": extrapolated_low_flag,
        "extrapolated_high_flag": extrapolated_high_flag,
    })

    # ---- Plot ----
    fig, ax = plt.subplots()
    ax.scatter(std_mean_bc, concs, label="BSA standards")
    x_line = np.linspace(np.nanmin(std_mean_bc), np.nanmax(std_mean_bc), 200)
    ax.plot(x_line, slope * x_line + intercept, label="Linear fit")

    ax.set_title("BSA standard curve (blank-corrected)")
    ax.set_xlabel("Absorbance at 620 nm (mean - blank)")
    ax.set_ylabel("BSA concentration (µg/mL)")

    ax.text(
        0.05, 0.95,
        f"{equation}\nR² = {r2:.4f}",
        transform=ax.transAxes,
        va="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.7),
    )
    ax.legend()

    if show_plot:
        plt.show()

    # ---- Print report (one-line UX) ----
    if print_report:
        print("\nBSA standard curve")
        print(f"  Equation: {equation}")
        print(f"  R²: {r2:.4f}")

        print("\nProtein sample concentrations (µg/mL)")
        report_cols = [
            "sample",
            "final_conc_ug_ml",
            "below_blank_flag",
            "extrapolated_low_flag",
            "extrapolated_high_flag",
        ]
        # Keep terminal output tidy
        print(samples_df[report_cols].to_string(index=False))

    # ---- Autosave CSV ----
    csv_path_str: Optional[str] = None
    if autosave_csv:
        if csv_out is None:
            csv_out = excel_path.with_name(excel_path.stem + "_protein_results.csv")
        csv_out = Path(csv_out)
        samples_df.to_csv(csv_out, index=False)
        csv_path_str = str(csv_out.resolve())
        if print_report:
            print(f"\nSaved results CSV: {csv_path_str}")

    return ProteinQuantResult(
        r2=r2,
        slope=slope,
        intercept=intercept,
        equation=equation,
        standards=standards_df,
        samples=samples_df,
        fig=fig,
        csv_path=csv_path_str,
    )