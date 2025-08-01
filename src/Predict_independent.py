"""
Parametric SIAMCAT analysis pipeline (Python ↔ R).
"""
from __future__ import annotations

import json
import locale
import os
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from onedal.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve, auc
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import MinMaxScaler
from tqdm import tqdm
from xgboost import XGBClassifier

os.environ["LC_CTYPE"] = "en_US.UTF-8"
import platform
import subprocess
from typing import Dict, List, Optional
import pandas as pd
from model_utils import _select_best_param,_load_best_param

# ───────────────────────── ensure R_HOME ────────────────────────────────
if "R_HOME" not in os.environ:
    if platform.system() == "Windows":
        os.environ["R_HOME"] = r"D:\software\R-4.4"
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
        except subprocess.SubprocessError as exc:
            raise EnvironmentError("Unable to locate R – is it on PATH?") from exc

import rpy2.rinterface_lib.callbacks as callbacks
import rpy2.rinterface_lib.conversion as conversion_lib
from rpy2 import robjects
from rpy2.robjects import StrVector, default_converter, pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter


# ─────────────────────────── R console silencer ──────────────────────────
conversion_lib._CCHAR_ENCODING = locale.getpreferredencoding()
callbacks.consolewrite_print = lambda *_, **__: None
callbacks.consolewrite_warnerror = lambda *_, **__: None
if platform.system() != "Windows":
    callbacks.consolewrite_ex = lambda *_, **__: None

def select_feature(groups:List[str],csv_path: str, p_col: str = "lefse") -> List[str]:
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

def analysis_pipeline():
    importr('base')
    importr('utils')

    feature_type = "16s"
    class_level, groups = 'otu', ['CTR', 'ADA']
    filter_thre = 0.0001
    norm_method = "log.std"
    feature_method = "wilcoxon"
    model_name = "xgb"
    thresholds = {'abundance': filter_thre}

    script_dir = os.path.dirname(os.path.abspath(__file__))
    siamcat_script = os.path.join(script_dir, 'siamcat_utils.R').replace('\\', '/')
    correct_script = os.path.join(script_dir, "correct_utils.R").replace("\\", '/')
    feature_script = os.path.join(script_dir, "feature_utils.R").replace("\\", '/')

    robjects.r['source'](siamcat_script)
    robjects.r['source'](correct_script)
    robjects.r['source'](feature_script)

    data_dir = os.path.abspath(os.path.join(script_dir, '..', 'Result', 'data', feature_type, 'Raw'))
    meta_dir = os.path.abspath(os.path.join(script_dir, '..', 'Result', 'data', feature_type))

    raw_df = pd.read_csv(os.path.join(data_dir, class_level, 'Raw_feature_finally.csv'), index_col=0)
    meta_df = pd.read_csv(os.path.join(meta_dir, 'Raw_meta_finally.csv'), index_col=0)
    group_str = "_".join(groups)
    feature_dir = os.path.abspath(os.path.join(script_dir, '..', 'Result', 'Feature', feature_type, group_str, class_level))
    param_dir=os.path.abspath(os.path.join(script_dir, '..', 'Result', 'figures','param'))

    flags = meta_df.groupby('Study')['Group'].apply(lambda g: set(groups).issubset(g))
    valid = flags[flags].index
    meta_f = meta_df[meta_df['Study'].isin(valid) & meta_df['Group'].isin(groups)]
    common = meta_f.index.intersection(raw_df.index)
    meta_f, raw_f = meta_f.loc[common], raw_df.loc[common]
    meta_f["SampleID"] = meta_f.index

    r_thresh = robjects.ListVector({k: robjects.FloatVector([v]) for k, v in thresholds.items()})
    filter_mode = list(thresholds.keys())[0]

    r_siamcat = robjects.globalenv["run_siamcat_filter_analysis"]
    with localconverter(default_converter + pandas2ri.converter):
        r_raw = pandas2ri.py2rpy(raw_f)
        r_meta = pandas2ri.py2rpy(meta_f)
        filter_df = r_siamcat(r_raw, r_meta, StrVector(groups), filter_mode, r_thresh)

    sig_feature = select_feature(groups,
        os.path.join(feature_dir, f'abundance{filter_thre}_{norm_method}_no_correct_feature_adj_table.csv'),
                                 feature_method)
    cols = sorted(set(sig_feature).intersection(raw_f.columns))
    if not cols:
        raise RuntimeError("No significant features found despite selection enabled")
    feat_mat = filter_df.T[cols]

    with localconverter(default_converter + pandas2ri.converter):
        data_raw_r = pandas2ri.py2rpy(feat_mat)
        r_meta2 = pandas2ri.py2rpy(meta_f)
    r_siamcat_norm = robjects.globalenv["run_siamcat_norm_analysis"]
    norm_r_df = r_siamcat_norm(data_raw_r, r_meta2, StrVector(groups), norm_method)
    with localconverter(default_converter + pandas2ri.converter):
        norm_df = pandas2ri.rpy2py(norm_r_df).T

    meta_feature = meta_f[["Group", "Study"]].join(norm_df[cols])
    meta_train = meta_feature[meta_feature["Study"] != "CHN_WF-CRC"].copy()
    meta_test = meta_feature[meta_feature["Study"] == "CHN_WF-CRC"].copy()

    raw_train = meta_train.drop(columns=["Group","Study"])
    raw_test = meta_test.drop(columns=["Group","Study"])
    y_train = meta_train["Group"].map({"CTR": 0, "ADA": 1}).values
    y_test = meta_test["Group"].map({"CTR": 0, "ADA": 1}).values

    print(f"[INFO] Train samples: {meta_train.shape[0]}, Test samples: {meta_test.shape[0]}")
    print(f"[INFO] Features used: {len(cols)}")

    param_path = os.path.join(param_dir,feature_type, "_".join(groups),class_level)
    root = Path(param_path)

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

    base_model, formal_name = model_map[model_name]
    params_df = _load_best_param(root, model_name)
    model = _select_best_param(base_model, params_df, formal_name)
    fpr_list = []
    tpr_list = []
    auc_list = []
    mean_fpr = np.linspace(0, 1, 100)

    for i in tqdm(range(100)):
        scaler = MinMaxScaler()
        raw_train = scaler.fit_transform(raw_train)
        raw_test = scaler.transform(raw_test)
        model.fit(raw_train, y_train)
        y_proba = model.predict_proba(raw_test)[:, 1]

        fpr, tpr, _ = roc_curve(y_test, y_proba)
        interp_tpr = np.interp(mean_fpr, fpr, tpr)
        interp_tpr[0] = 0.0
        interp_tpr[-1] = 1.0

        fpr_list.append(mean_fpr)
        tpr_list.append(interp_tpr)
        auc_list.append(auc(mean_fpr, interp_tpr))

    mean_tpr = np.mean(tpr_list, axis=0)
    std_tpr = np.std(tpr_list, axis=0)
    mean_auc = np.mean(auc_list)
    std_auc = np.std(auc_list)

    plt.figure()
    plt.plot(mean_fpr, mean_tpr, color='blue', lw=2,
             label=f'Mean ROC (AUC = {mean_auc:.4f} ± {std_auc:.4f})')
    plt.fill_between(mean_fpr, mean_tpr - std_tpr, mean_tpr + std_tpr,
                     color='blue', alpha=0.2, label='± 1 std. dev.')

    plt.plot([0, 1], [0, 1], 'k--', lw=1, label='Random')
    plt.xlim([0, 1])
    plt.ylim([0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Mean ROC Curve over 100 Runs')
    plt.legend(loc="lower right")
    plt.grid(True)
    plt.tight_layout()
    #plt.savefig(f"../Result/Model/{feature_type}/{group_str}/roc_curve_mean.pdf", dpi=300)
    plt.show()

if __name__ == '__main__':
    analysis_pipeline()
