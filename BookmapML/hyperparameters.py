import optuna
import subprocess
import os
import pandas as pd
import joblib
import lightgbm as lgb
from sklearn.model_selection import train_test_split
import re
import shutil

# =====================
# GLOBAL CONFIGURATION
# =====================
TRAIN_DATA = "train_data.csv"
VALID_IN = "validation/in"
VALID_GT = "validation/gt"
VALID_OUT = "validation/out"
VALID_BIN = "validation/out2"

NUM_TRIALS = 200  # number of hyperparameter configurations to try

# =====================
# FUNCTION TO TRAIN 3 MODELS (T1, T2, T3)
# =====================
def train_models(X, y1, y2, y3, params, num_rounds):
    dtrain1 = lgb.Dataset(X, label=y1)
    dtrain2 = lgb.Dataset(X, label=y2)
    dtrain3 = lgb.Dataset(X, label=y3)

    m1 = lgb.train(params, dtrain1, num_boost_round=num_rounds)
    m2 = lgb.train(params, dtrain2, num_boost_round=num_rounds)
    m3 = lgb.train(params, dtrain3, num_boost_round=num_rounds)

    joblib.dump(m1, "model_T1.pkl")
    joblib.dump(m2, "model_T2.pkl")
    joblib.dump(m3, "model_T3.pkl")

# =====================
# FUNCTION TO RUN SCORING PIPELINE
# =====================
def evaluate_score():
    # Clean previous outputs
    if os.path.exists(VALID_OUT):
        shutil.rmtree(VALID_OUT)
    if os.path.exists(VALID_BIN):
        shutil.rmtree(VALID_BIN)
    os.makedirs(VALID_OUT, exist_ok=True)
    os.makedirs(VALID_BIN, exist_ok=True)

    # Run the prediction pipeline
    subprocess.run([
        "python3", "predict2.py",
        "--indir", VALID_IN,
        "--gtdir", VALID_GT,
        "--outdir", VALID_OUT
    ], check=True)

    # Binarize predictions
    subprocess.run([
        "python3", "binarize_output.py",
        "--indir", VALID_OUT,
        "--outdir", VALID_BIN,
        "--threshold", "0.5"
    ], check=True)

    # Compute score
    result = subprocess.run([
        "python3", "scorer.py",
        "-g", VALID_GT,
        "-o", VALID_BIN
    ], capture_output=True, text=True, check=True)

    # Extract "Final score" from output
    match = re.search(r"Final score:\s*(\d+)", result.stdout)
    if match:
        return int(match.group(1))
    return 0

# =====================
# OPTUNA OBJECTIVE FUNCTION
# =====================
def objective(trial):
    # Define hyperparameter search space
    params = {
        "objective": "binary",
        "metric": "binary_logloss",
        "boosting_type": "gbdt",
        "verbosity": -1,
        "learning_rate": trial.suggest_float("learning_rate", 0.005, 0.05, log=True),
        "num_leaves": trial.suggest_int("num_leaves", 4, 64),
        "max_depth": trial.suggest_int("max_depth", 4, 16),
        "min_data_in_leaf": trial.suggest_int("min_data_in_leaf", 50, 500),
        "feature_fraction": trial.suggest_float("feature_fraction", 0.3, 0.8),
        "bagging_fraction": trial.suggest_float("bagging_fraction", 0.3, 0.8),
        "bagging_freq": trial.suggest_int("bagging_freq", 1, 5),
        "lambda_l1": trial.suggest_float("lambda_l1", 3.0, 10.0),
        "lambda_l2": trial.suggest_float("lambda_l2", 1.5, 10.0),
        "min_split_gain": trial.suggest_float("min_split_gain", 0.0, 0.5),
        "num_threads": -1
    }

    # Also optimize number of boosting rounds
    num_rounds = trial.suggest_int("num_boost_round", 10, 500)

    # Load training data
    df = pd.read_csv(TRAIN_DATA)
    features = [col for col in df.columns if col not in {"id", "p1", "p2", "p3"}]
    X = df[features]
    y1, y2, y3 = df["p1"], df["p2"], df["p3"]

    # Train models with current hyperparameters
    train_models(X, y1, y2, y3, params, num_rounds)

    # Evaluate models using external scorer
    score = evaluate_score()

    # Optuna minimizes by default; negate score to maximize it
    return -score

# =====================
# RUN OPTIMIZATION
# =====================
if __name__ == "__main__":
    study = optuna.create_study(direction="minimize")  # minimize negative score
    study.optimize(objective, n_trials=NUM_TRIALS)

    best_params = study.best_params
    best_score = -study.best_value
    print("\n===== BEST PARAMETERS FOUND =====")
    print(best_params)
    print(f"Best Final score on validation: {best_score}")
