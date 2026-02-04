import argparse
from .protein import protein_quantifier

def main():
    p = argparse.ArgumentParser(description="LabPod protein quantifier (BSA curve from Excel).")
    p.add_argument("excel_path", help="Path to Excel file (.xlsx)")
    p.add_argument("n_samples", type=int, help="Number of protein samples (rows A..)")
    p.add_argument("--sheet", default=0, help="Sheet name or index (default 0)")
    p.add_argument("--dilution", type=float, default=1.0, help="Dilution factor (e.g. 10 for 1:10)")
    p.add_argument("--no-plot", action="store_true", help="Disable plot display")
    p.add_argument("--no-csv", action="store_true", help="Do not autosave CSV")
    args = p.parse_args()

    protein_quantifier(
        args.excel_path,
        args.n_samples,
        sheet=args.sheet,
        dilution_factor=args.dilution,
        show_plot=not args.no_plot,
        autosave_csv=not args.no_csv,
        print_report=True,
    )