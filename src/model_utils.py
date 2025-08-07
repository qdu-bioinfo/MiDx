# ============================
# File: model_utils.R
# ============================

"""
Utility module for five‑fold cross‑validation of six classifier families
with best‑parameter injection **and automatic ROC‑AUC curve export**.

Updates
-------
* **2025‑06‑19** – Added combined ROC plot per model that overlays all five
  folds plus a bold mean‑ROC curve.  File is saved as ``<model>_roc_all.png``
  in each model's output directory.
* **2025‑06‑19** – Path‑agnostic parameter loading.  You now only need to
  pass the *root* directory that contains the ``<model>_best_params.pkl``
  files (recursively).  The script will locate the correct pickle no matter
  how deep it is nested (e.g. ``Result/figures/param/wgs/CTR_CRC/t_sgb``).

Dependencies
~~~~~~~~~~~~
``pandas scikit-learn imbalanced-learn xgboost matplotlib numpy joblib``
"""
from __future__ import annotations
import ast
import sys
from pathlib import Path
from typing import Dict, List
import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json, codecs
from imblearn.over_sampling import SMOTE
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    auc,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,  # noqa: F401 – kept for potential external use
    roc_curve,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import LabelEncoder, MinMaxScaler
from sklearn.svm import SVC
from xgboost import XGBClassifier

# ───────────────────────── hyper‑parameter helpers ─────────────────────
def _select_best_param(base_model, params_df: pd.DataFrame, model_name: str):
    """Return a *new* model instance with the best params injected.
    ``params_df`` can originate from either a plain ``dict`` that was saved
    with ``pickle.dump(best_dict)`` *or* from a DataFrame that stores that
    dict under the column ``best_params``.  Both formats are supported so
    that historical parameter files remain compatible.
    """
    if "model" not in params_df.columns:
        # legacy format ‑ first row is the param blob itself
        blob = params_df.iloc[0]
        best = blob["best_params"] if "best_params" in blob else blob
    else:
        sel = params_df.loc[params_df["model"] == model_name]
        if sel.empty:
            # No matching row – fall back to the un‑tuned base model
            return base_model
        blob = sel["best_params"].iloc[0]
        if isinstance(blob, str):
            blob = ast.literal_eval(blob)
        best = {
            k: (v[0] if isinstance(v, list) and len(v) == 1 else v)
            for k, v in blob.items()
        }
        if model_name == "knn" and "n_neighbors" in best:
            best["n_neighbors"] = int(round(best["n_neighbors"]))
    return base_model.__class__(**best)
# ───────────────────────── cross‑validation core ───────────────────────
def _cross_validate(
    model_key: str,
    X: pd.DataFrame,
    y: pd.Series,
    mapping :Dict[str, int],
    params_df: pd.DataFrame,
    output_dir: Path,
    n_splits: int = 5,
    random_state: int = 42,
):
    """Run *n_splits* stratified CV with SMOTE + MinMax scaling.

    Saves per‑fold models, scalers, ROC plots and a combined ROC plot with
    mean AUC, plus a CSV with all metrics.  Everything lives inside
    ``output_dir / <model_key>/``.
    """
    model_map = {
        "mlp": (MLPClassifier(), "mlp"),
        "xgb": (
            XGBClassifier(use_label_encoder=False, eval_metric="logloss"),
            "xgb",
        ),
        "rf": (RandomForestClassifier(), "rf"),
        "svm": (SVC(probability=True), "svm"),
        "knn": (KNeighborsClassifier(), "knn"),
        "lasso": (LogisticRegression(max_iter=10_000), "lasso"),
    }

    if model_key not in model_map:
        raise ValueError(f"Unsupported model: {model_key}")

    output_dir.mkdir(parents=True, exist_ok=True)
    base_model, formal_name = model_map[model_key]

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)

    y_enc = y.map(mapping).values

    records = []
    fprs: List[np.ndarray] = []
    tprs: List[np.ndarray] = []
    aucs: List[float] = []
    mean_fpr = np.linspace(0, 1, 101)

    for fold, (tr_idx, te_idx) in enumerate(skf.split(X, y_enc), 1):
        X_tr, X_te = X.iloc[tr_idx], X.iloc[te_idx]
        y_tr, y_te = y_enc[tr_idx], y_enc[te_idx]

        # ── balancing + scaling ────────────────────────────────────────
        X_res, y_res = SMOTE(random_state=random_state).fit_resample(X_tr, y_tr)
        scaler = MinMaxScaler()
        X_res = scaler.fit_transform(X_res)
        X_te_scaled = scaler.transform(X_te)

        # ── bestparams injection ─────────────────────────────────────
        model = _select_best_param(base_model, params_df, formal_name)
        model.fit(X_res, y_res)

        # ── predictions & metrics ─────────────────────────────────────
        if hasattr(model, "predict_proba"):
            scores = model.predict_proba(X_te_scaled)[:, 1]
        else:  # SVM without probability=True (but we set it True above)
            scores = model.decision_function(X_te_scaled)
        preds = model.predict(X_te_scaled)

        fpr, tpr, _ = roc_curve(y_te, scores)
        fold_auc = auc(fpr, tpr)
        fprs.append(fpr)
        tprs.append(tpr)
        aucs.append(fold_auc)

        # collect numeric metrics
        records.append({
            "model": model_key,
            "fold": fold,
            "accuracy": accuracy_score(y_te, preds),
            "roc_auc": fold_auc,
            "precision": precision_score(y_te, preds),
            "recall": recall_score(y_te, preds),
            "f1": f1_score(y_te, preds),
        })

    # ── combined ROC plot ─────────────────────────────────────────────
    plt.figure()
    for i, (fpr, tpr) in enumerate(zip(fprs, tprs), 1):
        plt.plot(fpr, tpr, lw=1, alpha=0.6, label=f"fold {i} (AUC={aucs[i-1]:.3f})")

    # mean ROC on common grid
    interp_tprs = [np.interp(mean_fpr, fprs[i], tprs[i]) for i in range(n_splits)]
    mean_tpr = np.mean(interp_tprs, axis=0)
    std_tpr = np.std(interp_tprs, axis=0)
    tpr_upper = np.minimum(mean_tpr + std_tpr, 1)
    tpr_lower = np.maximum(mean_tpr - std_tpr, 0)

    # shaded band ±1 SD
    plt.fill_between(mean_fpr, tpr_lower, tpr_upper, color="tab:blue", alpha=0.2, label="±1 SD")
    plt.plot(mean_fpr, mean_tpr, lw=2.5, color="blue", label=f"mean ROC (AUC={auc(mean_fpr, mean_tpr):.3f})")

    plt.plot([0, 1], [0, 1], lw=1, ls="--")
    plt.xlim([0, 1])
    plt.ylim([0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"{formal_name} – 5‑fold ROC Curves")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(output_dir / f"{model_key}_roc_all.pdf", dpi=300)
    plt.close()

    pd.DataFrame(records).to_csv(output_dir / f"{model_key}_cv_results.csv", index=False)
    print(f"[{model_key}] 5 fold CV finished results saved to {output_dir}, Mean AUROC={auc(mean_fpr, mean_tpr):.3f}")

def _lodo_validate(
    model_key: str,
    X: pd.DataFrame,
    y: pd.DataFrame,  # y should have both ["Group", "Study"]
    mapping: Dict[str, int],
    params_df: pd.DataFrame,
    output_dir: Path,
    random_state: int = 42,
):
    """
    Run LODO (Leave-One-Dataset-Out) validation across different Study groups.
    Saves per-study ROC curves and a combined ROC plot, plus metrics CSV.
    """
    model_map = {
        "mlp": (MLPClassifier(), "mlp"),
        "xgb": (
            XGBClassifier(use_label_encoder=False, eval_metric="logloss"),
            "xgb",
        ),
        "rf": (RandomForestClassifier(), "rf"),
        "svm": (SVC(probability=True), "svm"),
        "knn": (KNeighborsClassifier(), "knn"),
        "lasso": (LogisticRegression(max_iter=10_000), "lasso"),
    }

    if model_key not in model_map:
        raise ValueError(f"Unsupported model: {model_key}")

    output_dir.mkdir(parents=True, exist_ok=True)
    base_model, formal_name = model_map[model_key]

    y_label = y["Group"].map(mapping).values
    studies = y["Study"].values
    unique_studies = np.unique(studies)

    records = []
    fprs: List[np.ndarray] = []
    tprs: List[np.ndarray] = []
    aucs: List[float] = []
    mean_fpr = np.linspace(0, 1, 101)

    for study in unique_studies:
        test_mask = (studies == study)
        X_tr, X_te = X[~test_mask], X[test_mask]
        y_tr, y_te = y_label[~test_mask], y_label[test_mask]

        # ── balancing + scaling ────────────────────────────────────────
        X_res, y_res = SMOTE(random_state=random_state).fit_resample(X_tr, y_tr)
        scaler = MinMaxScaler()
        X_res = scaler.fit_transform(X_res)
        X_te_scaled = scaler.transform(X_te)

        # ── bestparams injection ─────────────────────────────────────
        model = _select_best_param(base_model, params_df, formal_name)
        model.fit(X_res, y_res)

        # ── predictions & metrics ─────────────────────────────────────
        if hasattr(model, "predict_proba"):
            scores = model.predict_proba(X_te_scaled)[:, 1]
        else:
            scores = model.decision_function(X_te_scaled)
        preds = model.predict(X_te_scaled)

        fpr, tpr, _ = roc_curve(y_te, scores)
        fold_auc = auc(fpr, tpr)
        fprs.append(fpr)
        tprs.append(tpr)
        aucs.append(fold_auc)

        records.append({
            "model": model_key,
            "study": study,
            "n_test": len(y_te),
            "accuracy": accuracy_score(y_te, preds),
            "roc_auc": fold_auc,
            "precision": precision_score(y_te, preds),
            "recall": recall_score(y_te, preds),
            "f1": f1_score(y_te, preds),
        })

    # ── combined ROC plot ─────────────────────────────────────────────
    plt.figure()
    for i, (fpr, tpr) in enumerate(zip(fprs, tprs)):
        plt.plot(fpr, tpr, lw=1, alpha=0.6, label=f"{unique_studies[i]} (AUC={aucs[i]:.3f})")

    interp_tprs = [np.interp(mean_fpr, fprs[i], tprs[i]) for i in range(len(unique_studies))]
    mean_tpr = np.mean(interp_tprs, axis=0)
    std_tpr = np.std(interp_tprs, axis=0)
    tpr_upper = np.minimum(mean_tpr + std_tpr, 1)
    tpr_lower = np.maximum(mean_tpr - std_tpr, 0)

    plt.fill_between(mean_fpr, tpr_lower, tpr_upper, color="tab:blue", alpha=0.2, label="±1 SD")
    plt.plot(mean_fpr, mean_tpr, lw=2.5, color="blue", label=f"mean ROC (AUC={auc(mean_fpr, mean_tpr):.3f})")

    plt.plot([0, 1], [0, 1], lw=1, ls="--")
    plt.xlim([0, 1])
    plt.ylim([0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"{formal_name} – LODO ROC Curves")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(output_dir / f"{model_key}_roc_lodo.pdf", dpi=300)
    plt.close()

    pd.DataFrame(records).to_csv(output_dir / f"{model_key}_lodo_results.csv", index=False)
    print(f"[{model_key}] LODO validation finished results saved to {output_dir}, Mean AUROC={auc(mean_fpr, mean_tpr):.3f}")

def _save_models(
    model_key: str,
    X: pd.DataFrame,
    y: pd.Series,
    mapping :Dict[str, int],
    params_df: pd.DataFrame,
    output_dir: Path,
    random_state: int = 42,
):
    """Run *n_splits* stratified CV with SMOTE + MinMax scaling.

    Saves per‑fold models, scalers, ROC plots and a combined ROC plot with
    mean AUC, plus a CSV with all metrics.  Everything lives inside
    ``output_dir / <model_key>/``.
    """
    model_map = {
        "mlp": (MLPClassifier(), "mlp"),
        "xgb": (
            XGBClassifier(use_label_encoder=False, eval_metric="logloss"),
            "xgb",
        ),
        "rf": (RandomForestClassifier(), "rf"),
        "svm": (SVC(probability=True), "svm"),
        "knn": (KNeighborsClassifier(), "knn"),
        "lasso": (LogisticRegression(max_iter=10_000), "lasso"),
    }
    if model_key not in model_map:
        raise ValueError(f"Unsupported model: {model_key}")

    output_dir.mkdir(parents=True, exist_ok=True)
    base_model, formal_name = model_map[model_key]
    y_enc = y.map(mapping).values
    X_res, y_res = SMOTE(random_state=random_state).fit_resample(X, y_enc)
    scaler = MinMaxScaler()
    X_res = scaler.fit_transform(X_res)
    # ── bestparams injection ─────────────────────────────────────
    model = _select_best_param(base_model, params_df, formal_name)
    model.fit(X_res, y_res)
    #Feature
    feature_cols=list(X.columns)
    # ── artefacts (model + scaler) ────────────────────────────────
    prefix = output_dir / f"{model_key}"
    joblib.dump(model, prefix.with_suffix(".ckpt"))
    joblib.dump(scaler, prefix.with_name(prefix.name + "_scaler.pkl"))

    json.dump(
        feature_cols,
        codecs.open(prefix.with_name(prefix.name + "_features.json"), "w", "utf-8"),
        indent=2,
        ensure_ascii=False,
    )
# ───────────────────────── parameter loader ───────────────────────────
def _load_best_param(param_root: Path, model_key: str) -> pd.DataFrame:
    """Recursively locate ``<model_key>_best_params.pkl`` under *param_root*.

    Returns the content converted to a pandas ``DataFrame`` that contains
    two columns: ``model`` and ``best_params`` (dict).  This keeps the
    downstream logic simple and uniform.
    """
    matches = list(param_root.rglob(f"{model_key}_best_params.pkl"))
    if not matches:
        raise FileNotFoundError(
            f"✖ Dont find '{model_key}_best_params.pkl' (search root dir: {param_root})"
        )
    pkl = matches[0]  # usually there is only one – take the first match
    try:
        params = pd.read_pickle(pkl)
    except Exception:
        params = joblib.load(pkl)  # fallback for old joblib‑saved files
    if isinstance(params, dict):
        params = pd.DataFrame([{"model": model_key, "best_params": params}])
    return params
# ───────────────────────── orchestrator ───────────────────────────────
def run_cv_models(df_or_path,mapping,param_root: str, output_dir: str):
    """High‑level API used by the CLI entry‑point."""
    df = pd.read_csv(df_or_path) if isinstance(df_or_path, (str, Path)) else df_or_path
    if "Group" not in df.columns:
        raise ValueError("Input dataframe must contain a 'Group' column")
    y = df["Group"]
    X = df.drop(columns=["Group", "Study"], errors="ignore")

    root = Path(param_root)
    for key in ["xgb", "rf", "knn", "lasso", "mlp","svm"]:
        params_df = _load_best_param(root, key)
        _cross_validate(key, X, y, mapping,params_df, Path(output_dir) / key)
        _save_models(key, X, y, mapping,params_df, Path(output_dir) / key)

def run_lodo_models(df_or_path,mapping,param_root: str, output_dir: str):
    """High‑level API used by the CLI entry‑point."""
    df = pd.read_csv(df_or_path) if isinstance(df_or_path, (str, Path)) else df_or_path
    if "Group" not in df.columns:
        raise ValueError("Input dataframe must contain a 'Group' column")

    y = df[["Group","Study"]]
    X = df.drop(columns=["Group", "Study"], errors="ignore")

    root = Path(param_root)
    for key in ["xgb", "rf", "knn", "lasso", "mlp","svm"]:
        params_df = _load_best_param(root, key)
        _lodo_validate(key, X, y, mapping,params_df, Path(output_dir) / key)


# ───────────────────────── CLI entry‑point ─────────────────────────────
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "Usage:\n"
            "  python model_utils_corrected.py <meta_feature.csv> <param_root_dir> <output_dir>\n\n"
            "Example:\n"
            "  python model_utils_corrected.py meta_feature.csv Result/figures/param Result/cv_output"
        )
        sys.exit(1)

    meta_csv, param_root, out_dir = sys.argv[1:]
    run_cv_models(meta_csv, param_root, out_dir)
