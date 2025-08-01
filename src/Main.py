# ============================
# File: Main.py
# ============================
#!/usr/bin/env python3
"""
Main entry-point for running the parametric SIAMCAT analysis pipeline.

Example
-------
python Main.py \
  --raw-data feat.csv \
  --meta-data meta.csv \
  --feature-type wgs \
  --groups CTR CRC ADA \
  --mapping CTR=0 CRC=1 ADA=2 \
  --class-level t_sgb \
  --filter-mode pass \
  --filter-threshold 0.0 \
  --norm-mode pass \
  --correction none \
  --model rf \
  --output-path ./Result
"""

import argparse
import os
import sys
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from Analysis_Pipeline import analysis_pipeline


def parse_mapping(mapping_list: list[str]) -> dict[str, int]:
    """Parse list of strings like ['CTR=0', 'CRC=1'] into dict{"CTR":0, "CRC":1}."""
    mapping = {}
    for item in mapping_list:
        try:
            label, val = item.split('=')
            mapping[label] = int(val)
        except Exception:
            raise ValueError(f"Invalid mapping entry '{item}'. Use LABEL=INTEGER format.")
    return mapping


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the SIAMCAT pipeline with customisable filtering, "
                    "normalisation, batch-correction, feature-selection and modelling."
    )

    # Required input files
    parser.add_argument("--raw-data", required=True,
                        help="CSV file containing raw feature table (samples × features)")
    parser.add_argument("--meta-data", required=True,
                        help="CSV file containing sample metadata (index = sample IDs)")
    parser.add_argument("--feature-type", required=True,
                        help="Data type tag (e.g. wgs)")
    # Classification level
    parser.add_argument("--class-level", default="t_sgb",
                        help="Classification level tag (default: t_sgb)")

    # Study design: groups and mapping
    parser.add_argument("--groups", nargs='+', required=True,
                        help="Group labels to retain (e.g. CTR CRC ADA)")
    parser.add_argument("--mapping", nargs='+', required=True,
                        help="Mapping of each group label to integer (e.g. CTR=0 CRC=1 ADA=2)")

    # Filtering / normalisation
    parser.add_argument("--filter-mode", default="pass",
                        help="Prevalence metric for filtering (default 'pass')")
    parser.add_argument("--filter-threshold", type=float, default=0.0,
                        help="Threshold for filtering (default 0.0)")
    parser.add_argument("--norm-mode", default="pass",
                        help="Abundance metric for normalisation (default 'pass')")

    parser.add_argument("--no-correction",action="store_false",
                        dest="do_correction",
                        help="Skip correction stage"
                        )
    parser.add_argument("--no-filter", action="store_false",
                        dest="do_filter",
                        help="Skip filtering/normalisation stage")
    parser.set_defaults(do_filter=True)
    parser.add_argument("--no-feature-selection", action="store_false",
                        dest="do_feat_sel",
                        help="Skip feature selection stage")
    parser.set_defaults(do_feat_sel=True)
    parser.add_argument("--no-norm", action="store_false",
                        dest="do_norm",
                        help="Skip normalisation stage")
    parser.set_defaults(do_norm=True)
    parser.add_argument("--no-cv", action="store_false",
                        dest="do_cv",
                        help="Skip cv model stage")
    parser.set_defaults(do_norm=True)
    parser.add_argument("--no-lodo", action="store_false",
                        dest="do_lodo",
                        help="Skip lodo stage")
    parser.set_defaults(do_norm=True)

    # Batch correction and feature selection method
    parser.add_argument("--correction", default="MMUPHin",
                        help="Batch correction method (MMUPHin, Combat, none)")
    parser.add_argument("--feature", default="lefse",
                        help="Column for feature stats (default 'lefse')")


    # Output directory
    parser.add_argument("--output-path", required=True,
                        help="Directory for all outputs")

    args = parser.parse_args()

    # Load data
    raw_df = pd.read_csv(args.raw_data, index_col=0)
    meta_df = pd.read_csv(args.meta_data, index_col=0)
    meta_df=meta_df[meta_df["Study"]!="CHN_WF-CRC"].copy()

    # Validate mapping
    mapping = parse_mapping(args.mapping)
    if set(mapping.keys()) != set(args.groups):
        raise ValueError("Mapping keys must exactly match --groups labels.")

    # Construct thresholds and correction
    thresholds = {args.filter_mode: args.filter_threshold} if args.do_filter else None
    correction = None if args.correction.strip().lower() == "none" else args.correction

    # Invoke pipeline
    analysis_pipeline(
        raw_df=raw_df,
        meta_df=meta_df,
        feature_type=args.feature_type,
        groups=args.groups,
        mapping=mapping,
        class_level=args.class_level,
        filter_mode=args.filter_mode,
        norm_mode=args.norm_mode,
        thresholds=thresholds,
        correction_method=correction,
        do_correction=args.do_correction,
        do_filter=args.do_filter,
        do_norm=args.do_norm,
        do_feature_selection=args.do_feat_sel,
        do_cv_model=args.do_cv,
        do_lodo_model=args.do_lodo,
        feature_method=args.feature,
        output_dir=args.output_path,
    )


if __name__ == "__main__":
    main()
