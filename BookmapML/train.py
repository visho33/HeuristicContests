import lightgbm as lgb
import pandas as pd
import joblib

# =====================
# GLOBAL CONFIGURATION
# =====================
params = {
    "objective": "binary",
    "metric": "binary_logloss",  # internal metric (not prioritized)
    "boosting_type": "gbdt",
    "verbosity": -1,
    "learning_rate": 0.014122056151496445,        # stable learning rate
    "num_leaves": 23,             # small trees
    "max_depth": 6,               # shallow depth
    "min_data_in_leaf": 210,      # large leaves to regularize
    "feature_fraction": 0.5874960909666207,      # use feature subsets
    "bagging_fraction": 0.6915622199541117,      # use data subsets
    "bagging_freq": 1,
    "lambda_l1": 5.100814388893182,             # moderate L1 regularization
    "lambda_l2": 2.615051404932521,             # moderate L2 regularization
    "min_split_gain": 0.30460133528111405,       # avoid trivial splits
    "num_threads": -1
}

# Fixed number of boosting rounds (no early stopping)
num_boost_round = 200

# =====================
# LOAD AND PREPARE DATA
# =====================
df = pd.read_csv("train_data.csv")
features = [col for col in df.columns if col not in {"id", "p1", "p2", "p3"}]
X = df[features]
y1, y2, y3 = df["p1"], df["p2"], df["p3"]

# =====================
# FUNCTION TO TRAIN WITHOUT VALIDATION
# =====================
def train_full(X, y, params, model_name):
    dtrain = lgb.Dataset(X, label=y)
    print(f"\nTraining {model_name} (no validation)...")
    model = lgb.train(
        params,
        dtrain,
        num_boost_round=num_boost_round
    )
    joblib.dump(model, f"{model_name}.pkl")
    print(f"Model {model_name} saved.")
    return model

# =====================
# TRAIN THE THREE MODELS (T1, T2, T3)
# =====================
model1 = train_full(X, y1, params, "model_T1")
model2 = train_full(X, y2, params, "model_T2")
model3 = train_full(X, y3, params, "model_T3")
