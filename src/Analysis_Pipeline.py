# ============================
# File: Analysis_Pipeline.py
# ============================
"""
Parametric SIAMCAT analysis pipeline (Python ↔ R).
"""
from __future__ import annotations

import json
import locale
import os
os.environ["LC_CTYPE"] = "en_US.UTF-8"

import platform
import subprocess
from typing import Dict, List, Optional
import pandas as pd
from pathlib import Path
script_dir = Path(__file__).resolve().parent
print(f"Script Directory: {script_dir}")
# ───────────────────────── ensure R_HOME ────────────────────────────────
if "R_HOME" not in os.environ:
    if platform.system() == "Windows":
        os.environ["R_HOME"] = rf"{script_dir}\R-4.4"
        os.environ["PATH"] = os.pathsep.join([
            os.path.join(os.environ["R_HOME"], "bin", "x64"),
            os.path.join(os.environ["R_HOME"], "bin"),
            os.environ.get("PATH", ""),
        ])
    else:
        try:
            r_home = subprocess.check_output(["R", "RHOME"], stderr=subprocess.DEVNULL)
            os.environ["R_HOME"] = r_home.decode().strip()
            os.environ["PATH"] = os.path.join(os.environ["R_HOME"], "bin") + os.pathsep + os.environ.get("PATH", "")
            #os.environ["LD_LIBRARY_PATH"] = os.path.join(R_HOME, "lib") + os.pathsep + os.environ.get("LD_LIBRARY_PATH","")
        except subprocess.SubprocessError as exc:
            raise EnvironmentError("Unable to locate R – is it on PATH?") from exc
import rpy2.rinterface_lib.callbacks as callbacks
import rpy2.rinterface_lib.conversion as conversion_lib
from rpy2 import robjects
from rpy2.robjects import StrVector, default_converter, pandas2ri
from rpy2.robjects.conversion import localconverter
from model_utils import run_cv_models,run_lodo_models
# ─────────────────────────── R console silencer ──────────────────────────
conversion_lib._CCHAR_ENCODING = locale.getpreferredencoding()
callbacks.consolewrite_print = lambda *_, **__: None
callbacks.consolewrite_warnerror = lambda *_, **__: None
if platform.system() != "Windows":
    callbacks.consolewrite_ex = lambda *_, **__: None
# ───────────────────────── helper ───────────────────────────────────────
def _select_feature(groups:List[str],csv_path: str, p_col: str = "lefse") -> List[str]:
    alpha = 0.05 if "CRC" in groups else 0.5
    df = pd.read_csv(csv_path, index_col=0)
    df=df.fillna(1)
    if p_col=="All":
        features = df.index.tolist()
    else:
        features = df.loc[df[p_col] < alpha].index.tolist()
    count = len(features)
    print(f"{p_col} filter to {count} feature")
    return features
# ───────────────────────── main driver ──────────────────────────────────
def analysis_pipeline(
    *,
    raw_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    feature_type: str,
    groups: List[str],
    mapping :Dict[str, int],
    class_level: str,
    filter_mode: str,
    norm_mode: str,
    thresholds: Optional[Dict[str, float]],
    correction_method: Optional[str],
    do_correction:bool=True,
    do_filter: bool = True,
    do_norm:bool=True,
    do_feature_selection: bool = True,
    do_cv_model: bool = True,
    do_lodo_model: bool = True,
    feature_method: str,
    output_dir: str,
) -> None:
    """Execute the complete SIAMCAT-based workflow."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    for r_file in ["siamcat_utils.R", "correct_utils.R", "feature_utils.R"]:
        robjects.r["source"](os.path.join(script_dir, r_file).replace("\\", "/"))
    # ── sample subset ────────────────────────────────────────────────────
    valid_studies = (
        meta_df.groupby("Study")["Group"]
        .apply(lambda s: set(groups).issubset(s))
        .loc[lambda x: x]
        .index
    )
    meta_f = meta_df.loc[meta_df["Study"].isin(valid_studies) & meta_df["Group"].isin(groups)].copy()
    common_idx = meta_f.index.intersection(raw_df.index)

    meta_f = meta_f.loc[common_idx]
    raw_f = raw_df.loc[common_idx]
    meta_f["SampleID"] = meta_f.index

    print("Studies (JSON):", json.dumps(list(valid_studies), ensure_ascii=False))
    with localconverter(default_converter + pandas2ri.converter):
        r_raw = pandas2ri.py2rpy(raw_f)
        r_meta = pandas2ri.py2rpy(meta_f)

    # ── filtering / normalisation ────────────────────────────────────────
    if do_filter:
        if thresholds is None:
            raise ValueError("`thresholds` dict required when filtering enabled")
        print(f"[INFO] Pre-filter shape: {raw_f.shape[0]} samples × {raw_f.shape[1]} features")
        r_thresh = robjects.ListVector({k: robjects.FloatVector([v]) for k, v in thresholds.items()})
        filter_df = robjects.globalenv["run_siamcat_filter_analysis"](
            r_raw, r_meta, StrVector(groups),
           filter_mode, r_thresh
        )
        dim_after=tuple(robjects.r["dim"](filter_df))
        print(f"[INFO] Post-filter shape: {int(dim_after[1])} samples × {int(dim_after[0])} features")
    else:
        filter_df = robjects.r['t'](r_raw)
        print(f"[INFO] No filtering .")

    # ── batch correction ────────────────────────────────────────────────
    if do_correction:
        corr_df = robjects.globalenv["run_correct_analysis"](
            class_level, StrVector(groups), correction_method, filter_df, r_meta
        )
        print(f"Correction successful")
    else:
        corr_df = filter_df
        print(f"No correction")

    # ── feature selection ───────────────────────────────────────────────
    if thresholds is not None:
        thresh_tag = "_".join(f"{k}{v}" for k, v in thresholds.items())
    else:
        thresh_tag = "nofilter"
    feature_path = os.path.join(output_dir, "Feature", feature_type, "_".join(groups), class_level)
    os.makedirs(feature_path, exist_ok=True)
    #feature selection
    if do_feature_selection:
        adj_table_r = robjects.globalenv["run_feature_analysis"](corr_df, r_meta,
                                                                 os.path.join(feature_path,f"{thresh_tag}_{norm_mode}_{correction_method}"))
        with localconverter(default_converter + pandas2ri.converter):
            adj_table = pandas2ri.rpy2py(adj_table_r)
        adj_csv = os.path.join(feature_path,f"{thresh_tag}_{norm_mode}_{correction_method}_feature_adj_table.csv")
        adj_table.to_csv(adj_csv)
        print(f"[INFO] Adjusted feature table saved → {adj_csv}")
    else:
        adj_csv = os.path.join(feature_path, f"{thresh_tag}_{norm_mode}_{correction_method}_feature_adj_table.csv")
    sig = _select_feature(groups,adj_csv, feature_method)
    cols = sorted(set(sig).intersection(raw_f.columns))
    if not cols:
        raise RuntimeError("No significant features found despite selection enabled")
    feat_mat = pandas2ri.rpy2py(corr_df).T[cols]

    with localconverter(default_converter + pandas2ri.converter):
        data_raw = pandas2ri.py2rpy(feat_mat)
        r_meta = pandas2ri.py2rpy(meta_f)
    if do_norm:
        if thresholds is None:
            raise ValueError("`thresholds` dict required when filtering enabled")
        print(f"[INFO] Pre-norm shape: {feat_mat.shape[0]} samples × {feat_mat.shape[1]} features")
        norm_df = robjects.globalenv["run_siamcat_norm_analysis"](
           data_raw, r_meta, StrVector(groups),
           norm_mode
        )
        with localconverter(default_converter + pandas2ri.converter):
            norm_df = pandas2ri.rpy2py(norm_df).T  # samples × features
        print(f"[INFO] Post-norm shape: {norm_df.shape}")
    else:
        norm_df =feat_mat
        print(f"[INFO] No normalization.")

    # ── modelling ───────────────────────────────────────────────────────
    model_df = meta_f[["Group", "Study"]].join(norm_df)
    param_fp = os.path.join(output_dir, "figures/param", feature_type, "_".join(groups), class_level)
    if not os.path.exists(param_fp):
        raise FileNotFoundError(param_fp)

    model_path = os.path.join(output_dir, "Model", feature_type, "_".join(groups), class_level,f"{thresh_tag}_{norm_mode}_{correction_method}",feature_method)
    os.makedirs(model_path, exist_ok=True)
    if do_cv_model:
        run_cv_models(model_df,mapping,param_fp, model_path)
    if do_lodo_model:
        run_lodo_models(model_df, mapping, param_fp, model_path)
    print("[INFO] Full workflow finished successfully.")
