#!/usr/bin/env python

import joblib
import optuna
import functools
import pandas as pd
import numpy as np
from optuna.samplers import TPESampler
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import GridSearchCV, train_test_split, StratifiedKFold, cross_val_score
from sklearn.metrics import r2_score, mean_absolute_error


def choose_options(trial, X_train_std, y_train, n_fold, stratify_labels):
    """Obtain the best params"""
    # l1ration & alpha
    l1_ratio = trial.suggest_float("l1_ratio", low=0, high=1, step=0.05) # l1_ratio：弹性网混合比例（0=纯Ridge，1=纯Lasso），连续空间 0~1, 步长0.05，兼顾精度与效率
    alpha = trial.suggest_float("alpha", low=1e-4, high=10, log=True) # 1e-4 ~ 10，对数刻度

    # current params
    print(f"Optuna Trial {trial.number}: l1_ratio={l1_ratio:.2f}, alpha={alpha:.4f}")
    
    # cross-validation: cross_val_score
    cv = StratifiedKFold(n_splits = n_fold, shuffle = True, random_state = 42)
    
    mae_scores = []
    for train_idx, val_idx in cv.split(X_train_std, stratify_labels):
        X_train_cv, X_val_cv = X_train_std[train_idx], X_train_std[val_idx]
        y_train_cv, y_val_cv = y_train[train_idx], y_train[val_idx]

        # model init
        model = ElasticNet(l1_ratio = l1_ratio, alpha = alpha, max_iter = 10000, random_state = 42)
        model.fit(X_train_cv, y_train_cv)
        y_val_pred = model.predict(X_val_cv)
        mae = mean_absolute_error(y_val_cv, y_val_pred)
        mae_scores.append(mae)

        # 因 sklearn 评分函数需“越大越好”，故用负 MAE，后续取反得到真实 MAE
        #mae_scores = -cross_val_score(model, X_train_std, y = stratify_labels, cv = cv, scoring = "neg_mean_absolute_error")
    
    # 为了让超参数优化器（Optuna）能基于模型在交叉验证中的平均性能来判断超参数的优劣，从而找到泛化能力最强的参数组合
    print(f"Optuna Trial {trial.number}: Cross-Validation mean MAE={np.mean(mae_scores):.4f}|Fold MAE：{mae_scores}")

    return np.mean(mae_scores)


### prepare data
# DMPs
dmp_df = pd.read_csv("age_related_DMPs.txt", sep="\t")
dmp_cpgs = dmp_df["CpG"].tolist()

# beta matrix
beta_df = pd.read_csv("CpG.value.txt", sep="\t", index_col=0)
beta_df = beta_df.T # row: sample; col: cpg
beta_df = beta_df[dmp_cpgs] 

# pheno
pheno_df = pd.read_csv("pheno.info.txt", sep="\t", index_col=0)
common_samples = list(set(beta_df.index) & set(pheno_df.index)) 
beta_df = beta_df.loc[common_samples]
pheno_df = pheno_df.loc[common_samples]

y = pheno_df["Age"].values
covar = pheno_df["Gender2"].values.reshape(-1,1)
X = np.hstack([beta_df.values, covar]) 

covariate_names = ["Gender"]
feature_names = dmp_cpgs + covariate_names
with open("feature_names.txt", "w") as f:
    f.write("\n".join(feature_names))

with open("covariate_names.txt", "w") as f:
    f.write("\n".join(covariate_names))

# 分层抽样
stratify_labels = pheno_df["Age_group"].values.ravel()
#print(stratify_labels) 

#print(pheno_df["Age_group"].value_counts())
X_train, X_test, y_train, y_test, stratify_labels_train, stratify_labels_test = train_test_split(X, y, stratify_labels, test_size=0.2, random_state=42, stratify=stratify_labels)
# print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)

# missing value
X_train_imputed = pd.DataFrame(X_train).fillna(pd.DataFrame(X_train).mean())
X_test_imputed = pd.DataFrame(X_test).fillna(pd.DataFrame(X_train).mean())  # 测试集用训练集均值插补

# scale
scaler = StandardScaler()
X_train_std = scaler.fit_transform(X_train_imputed)
X_test_std = scaler.transform(X_test_imputed) # 用训练集的平均值和标准差来标准化测试集

# scaler.mean_ & scaler.var_
with open("scaler_params.pkl", "wb") as f:
    joblib.dump(scaler, f)

# train feature mean
train_features_mean = pd.DataFrame(X_train).mean()

### train model
# start Optuna logging
optuna.logging.set_verbosity(optuna.logging.INFO)

# run optuna TPE 采样器：贝叶斯优化的核心，比默认采样器更高效;优化目标，最小化MAE
study = optuna.create_study(study_name="biological_age_elastic_net",
    sampler = TPESampler(seed=42),
    direction="minimize" 
    )

# 固定 n_fold，生成新的无额外参数的函数
object_params = functools.partial(choose_options, 
    X_train_std=X_train_std,
    y_train=y_train, 
    n_fold=5, 
    stratify_labels=stratify_labels_train)

study.optimize(object_params, n_trials = 50, n_jobs=1, show_progress_bar = True)

# best result
best_params = study.best_params
best_mae = study.best_value
print(f"Best params: l1_ratio={best_params['l1_ratio']}, alpha={best_params['alpha']}")

# restore all results of Optuna, including all options and mae
study_results = study.trials_dataframe()
study_results.to_csv("optuna_search_results.csv", index=False)

# use best params to init model 
best_model = ElasticNet(l1_ratio = best_params["l1_ratio"], 
    alpha = best_params["alpha"], 
    max_iter = 10000, 
    random_state = 42)

best_model.fit(X_train_std, y_train)

y_pred_test = best_model.predict(X_test_std)
test_mae = mean_absolute_error(y_test, y_pred_test)
test_r2 = r2_score(y_test, y_pred_test)
print(f"test dataset: MAE={test_mae}, R2={test_r2}")

# restore model
with open("elastic_net_model.pkl", "wb") as f:
    joblib.dump(best_model, f)

# model package
model_package = {
	"model": best_model,
	"scaler": scaler,
	"train_cpg_names": dmp_cpgs,
	"train_features_mean": train_features_mean
}

joblib.dump(model_package, "bio_age_prediction_model.pkl")
### FastAPI SDK: support other dataset to prediction