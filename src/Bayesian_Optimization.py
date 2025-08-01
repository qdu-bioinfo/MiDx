# ============================
# File:Bayesian_Optimization.R
# ============================

import os
import pickle
import sys

import xgboost as xgb
import numpy as np
from bayes_opt import BayesianOptimization
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn.preprocessing import MinMaxScaler
from xgboost import XGBClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.linear_model import Lasso
from dataprocess import split_data, select_feature
from sklearn.linear_model import LogisticRegression
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='xgboost')
xgb.set_config(verbosity=0)
warnings.filterwarnings('ignore')

SEED = 42

new_path = os.path.abspath(os.path.join(__file__, "../../"))
sys.path[1] = new_path
# Define paths for figures, features, and data
param_path = sys.path[1] + '\\Result\\figures\\param\\'

# -------------------------------------------------------------------- #
#                           Utility Functions                          #
# -------------------------------------------------------------------- #
def change_group(meta_feature, group_name):
    """
    Modify 'Group' in-place for binary classification.
    E.g. 'CTR_ADA' => if x == 'ADA' => 1 else 0
    """
    if group_name == 'CTR_ADA':
        meta_feature['Group'] = meta_feature['Group'].apply(lambda x: 1 if x == "ADA" else 0)
    elif group_name == "CTR_CRC":
        meta_feature['Group'] = meta_feature['Group'].apply(lambda x: 1 if x == "CRC" else 0)
    else:
        meta_feature['Group'] = meta_feature['Group'].apply(lambda x: 1 if x == "CRC" else 0)
    return meta_feature

# -------------------------------------------------------------------- #
#                           K-Fold Functions                           #
# -------------------------------------------------------------------- #
def get_kfold(data, meta_group, model_type, **params):
    aucs = []
    tprs = []
    mean_fpr = np.linspace(0, 1, 100)
    splitor = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)

    for train_index, test_index in splitor.split(data, meta_group):
        X_train, X_test = data.iloc[train_index].values, data.iloc[test_index].values
        y_train, y_test = meta_group.iloc[train_index].values, meta_group.iloc[test_index].values

        scaler = MinMaxScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        if model_type == 'rf':
            # RandomForest
            clf = RandomForestClassifier().set_params(**params)
            probas = clf.fit(X_train, y_train).predict_proba(X_test)

        elif model_type == 'xgb':
            # XGBoost
            clf = XGBClassifier().set_params(**params)
            probas = clf.fit(X_train, y_train).predict_proba(X_test)

        elif model_type == 'svm':
            # Support Vector Machine
            clf = SVC(probability=True).set_params(**params)
            probas = clf.fit(X_train, y_train).predict_proba(X_test)

        elif model_type == 'knn':
            # KNN
            clf = KNeighborsClassifier().set_params(**params)
            probas = clf.fit(X_train, y_train).predict_proba(X_test)

        elif model_type == 'mlp':
            # MLP
            clf = MLPClassifier().set_params(**params)
            probas = clf.fit(X_train, y_train).predict_proba(X_test)
        elif model_type == 'lasso':
            # MLP
            clf = LogisticRegression(**params)
            probas = clf.fit(X_train, y_train).predict_proba(X_test)

        else:
            raise ValueError("Invalid model type")

        fpr, tpr, thresholds = roc_curve(y_test, probas[:, 1])
        roc_auc = auc(fpr, tpr)
        print(f"{model_type.upper()} - AUC: {roc_auc:.3f}")
        aucs.append(roc_auc)
        tprs.append(np.interp(mean_fpr, fpr, tpr))

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    return mean_auc, mean_tpr, mean_fpr, std_auc, tprs

# -------------------------------------------------------------------- #
#                       Get K-Fold AUC Helpers                         #
# -------------------------------------------------------------------- #
def get_kfold_auc_op(data, meta_group, model_name, **params):
    """
    A thin wrapper that returns the mean AUC from get_kfold
    (used for Bayesian Optimization).
    """
    mean_auc, _, _, _, _ = get_kfold(data, meta_group, model_name, **params)
    return mean_auc

# -------------------------------------------------------------------- #
#                   Bayesian Optimization (RF, XGB, SVM, KNN, MLP)     #
# -------------------------------------------------------------------- #
def bayesian_optimise_rf(X, y, clf_kfold, n_iter, init_points=5):
    """
    Bayesian optimization for RandomForest hyperparameters
    """
    def rf_crossval(n_estimators, max_features, max_depth, max_samples):
        return clf_kfold(
            data=X,
            meta_group=y,
            model_name="rf",
            n_estimators=int(n_estimators),
            max_samples=max(min(max_samples, 0.999), 0.1),
            max_features=max(min(max_features, 0.999), 0.1),
            max_depth=int(max_depth),
            bootstrap=True
        )

    optimizer = BayesianOptimization(
        f=rf_crossval,
        pbounds={
            "n_estimators": (100, 500),
            "max_features": (0.1, 0.999),
            "max_samples": (0.1, 0.999),
            "max_depth": (1, 10)
        },
        random_state=SEED
    )
    optimizer.maximize(n_iter=n_iter, init_points=init_points)
    print("RandomForest Final result:", optimizer.max)
    return optimizer.max

def bayesian_optimise_xgb(X, y, clf_kfold, n_iter, init_points=5):
    """
    Bayesian optimization for XGBoost hyperparameters
    """
    def xgb_crossval(n_estimators, max_depth, learning_rate, subsample, colsample_bytree):
        return clf_kfold(
            data=X,
            meta_group=y,
            model_name="xgb",
            n_estimators=int(n_estimators),
            max_depth=int(max_depth),
            learning_rate=learning_rate,
            subsample=subsample,
            colsample_bytree=colsample_bytree,
            use_label_encoder=False,
            eval_metric='logloss'
        )
    optimizer = BayesianOptimization(
        f=xgb_crossval,
        pbounds={
            "n_estimators": (100, 1000),
            "max_depth": (3, 10),
            "learning_rate": (0.01, 0.3),
            "subsample": (0.2, 1.0),
            "colsample_bytree": (0.2, 1.0)
        },
        random_state=SEED
    )
    optimizer.maximize(n_iter=n_iter, init_points=init_points)
    print("XGBoost Final result:", optimizer.max)
    return optimizer.max

def bayesian_optimise_svm(X, y, clf_kfold, n_iter, init_points=5):
    """
    Bayesian optimization for SVM (RBF kernel) with log-scale for C, gamma.
    """
    def svm_crossval(log_C, log_gamma):
        C = 10 ** log_C
        gamma = 10 ** log_gamma
        return clf_kfold(
            data=X,
            meta_group=y,
            model_name="svm",
            C=C,
            gamma=gamma,
            kernel='rbf'
        )

    optimizer = BayesianOptimization(
        f=svm_crossval,
        pbounds={
            "log_C":    (-3, 3),
            "log_gamma":(-3, 3),
        },
        random_state=SEED
    )
    optimizer.maximize(n_iter=n_iter, init_points=init_points)

    # —— 新增：把 log_* 转回 real_* ——
    best = optimizer.max  # dict with keys 'target' 和 'params'
    raw = best['params']
    converted = {
        'C':     10 ** raw['log_C'],
        'gamma': 10 ** raw['log_gamma'],
        'kernel':'rbf'
    }
    best['params'] = converted

    print("SVM Final result (converted):", best)
    return best


def bayesian_optimise_knn(X, y, clf_kfold, n_iter, init_points=5):
    """
    Bayesian optimization for KNN hyperparameters.
    Example: n_neighbors in [1..50], p in [1..2].
    """
    def knn_crossval(n_neighbors, p):
        return clf_kfold(
            data=X,
            meta_group=y,
            model_name="knn",
            n_neighbors=int(n_neighbors),
            p=int(p)
        )
    optimizer = BayesianOptimization(
        f=knn_crossval,
        pbounds={
            "n_neighbors": (1, 50),
            "p": (1, 2),
        },
        random_state=SEED
    )
    optimizer.maximize(n_iter=n_iter, init_points=init_points)
    print("KNN Final result:", optimizer.max)
    return optimizer.max


def bayesian_optimise_mlp(X, y, clf_kfold, n_iter=20, init_points=5):
    """
    Bayesian optimization for MLP hyperparameters.
    log_alpha/log_lr → alpha/learning_rate_init 转换。
    """
    hidden_candidates = [(32,), (64,), (32, 32), (64, 32)]

    def mlp_crossval(h_idx, log_alpha, log_lr, max_iter):
        hidden = hidden_candidates[int(round(h_idx))]
        alpha = 10 ** log_alpha
        lr = 10 ** log_lr
        return clf_kfold(
            data=X,
            meta_group=y,
            model_name="mlp",
            hidden_layer_sizes=hidden,
            alpha=alpha,
            learning_rate_init=lr,
            max_iter=int(round(max_iter))
        )

    pbounds = {
        'h_idx': (0, len(hidden_candidates) - 1),
        'log_alpha': (-5, -1),  # alpha ∈ [1e-5, 1e-1]
        'log_lr': (-4, -1),  # lr ∈ [1e-4, 1e-1]
        'max_iter': (100, 1000),
    }

    optimizer = BayesianOptimization(
        f=mlp_crossval,
        pbounds=pbounds,
        random_state=SEED
    )
    optimizer.maximize(init_points=init_points, n_iter=n_iter)

    best = optimizer.max  # dict: {'target': ..., 'params': {...}}
    raw = best['params']
    converted = {
        'hidden_layer_sizes': hidden_candidates[int(round(raw['h_idx']))],
        'alpha': 10 ** raw['log_alpha'],
        'learning_rate_init': 10 ** raw['log_lr'],
        'max_iter': int(round(raw['max_iter']))
    }
    best['params'] = converted

    print("MLP Final result (converted):", best)
    return best
def bayesian_optimise_lasso(X, y, clf_kfold, n_iter=20, init_points=5):
    """
    Bayesian optimization for L1-regularized LogisticRegression (lasso).
    log_C → C 转换。
    """
    def lasso_crossval(log_C):
        C = 10 ** log_C
        return clf_kfold(
            data=X,
            meta_group=y,
            model_name="lasso",
            C=C,
            penalty='l1',
            solver='saga',
            max_iter=1000
        )

    optimizer = BayesianOptimization(
        f=lasso_crossval,
        pbounds={'log_C': (-4, 2)},  # C ∈ [1e-4, 1e2]
        random_state=SEED
    )
    optimizer.maximize(init_points=init_points, n_iter=n_iter)

    best = optimizer.max
    raw = best['params']
    converted = {
        'C': 10 ** raw['log_C'],
        'penalty': 'l1',
        'solver': 'saga',
        'max_iter': 1000
    }
    best['params'] = converted

    print("Lasso (Logistic L1) Final result (converted):", best)
    return best

# -------------------------------------------------------------------- #
#             Main Pipeline: finally_all_model_result                  #
# -------------------------------------------------------------------- #
def finally_all_model_result():
    # ------------------ Data Loading / Preprocessing ------------------ #
    for analysis_level in ["t_sgb"]:
    #analysis_level = "species"
        for data_type in ["wgs"]:
            group_name = "CTR_CRC"
            #data_type = "Raw"
            meta_feature= split_data(analysis_level,group_name,data_type)
            # Filter meta_feature on specific studies
            study_groups = meta_feature.groupby("Study")["Group"].agg(lambda vals: set(vals.unique()))
            studies_of_interest = [study for study, groups in study_groups.items()
                                   if {"CTR", "CRC"}.issubset(groups)]

            # 4. 过滤出有效 Study，并只保留 CTR/CRC 样本
            meta_feature = meta_feature.loc[
                (meta_feature["Study"].isin(studies_of_interest)) &
                (meta_feature["Group"].isin(["CTR", "CRC"]))
                ]
            data = meta_feature.iloc[:, :-2]
            y_data = change_group(meta_feature[["Group"]], group_name)

            # ------------------ Helper to Save Best Params ------------------ #
            def save_best_param(tune_result, model_name):
                """
                Convert known hyperparams to int if they exist,
                then save them in a pickle file for reference.
                """

                if 'n_estimators' in tune_result['params']:
                    tune_result['params']['n_estimators'] = int(tune_result['params']['n_estimators'])
                if 'max_depth' in tune_result['params']:
                    tune_result['params']['max_depth'] = int(round(tune_result['params']['max_depth']))
                if model_name=='knn':
                    if 'n_neighbors' in tune_result['params']:
                        tune_result['params']['n_neighbors']=int(round(tune_result['params']['n_neighbors']))
                    if 'p' in tune_result['params']:
                        tune_result['params']['p'] = int(round(tune_result['params']['p']))

                best_param = tune_result['params']
                filename = f"{param_path}/{data_type}/{group_name}/{analysis_level}/{model_name}_best_params.pkl"

                directory = os.path.dirname(filename)
                if not os.path.exists(directory):
                    os.makedirs(directory)
                with open(filename, 'wb') as f:
                    pickle.dump(best_param, f)
                print(f"Saved best params for {model_name}: {best_param}")

            # ------------------ 1) RandomForest ------------------ #
            tune_result_rf = bayesian_optimise_rf(data, y_data, get_kfold_auc_op, n_iter=100, init_points=5)
            save_best_param(tune_result_rf, "rf")

            #------------------ 2) XGBoost ------------------ #
            tune_result_xg = bayesian_optimise_xgb(data, y_data, get_kfold_auc_op, n_iter=100, init_points=5)
            save_best_param(tune_result_xg, "xgb")

            # ------------------ 3) SVM ------------------ #
            tune_result_svm = bayesian_optimise_svm(data, y_data, get_kfold_auc_op, n_iter=100, init_points=5)
            save_best_param(tune_result_svm, "svm")

            # ------------------ 4) KNN ------------------ #
            tune_result_knn = bayesian_optimise_knn(data, y_data, get_kfold_auc_op, n_iter=100, init_points=5)
            save_best_param(tune_result_knn, "knn")

            # ------------------ 5) MLP ------------------ #
            tune_mlp = bayesian_optimise_mlp(data, y_data, get_kfold_auc_op, n_iter=100, init_points=5)
            save_best_param(tune_mlp, "mlp")
            # ------------------ 5) lasso ------------------ #
            tune_lasso = bayesian_optimise_lasso(data, y_data, get_kfold_auc_op, n_iter=100, init_points=5)
            save_best_param(tune_lasso, "lasso")

if __name__ == "__main__":
    finally_all_model_result()
