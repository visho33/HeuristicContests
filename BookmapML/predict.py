import pandas as pd
import joblib
import argparse
import os
from collections import defaultdict
import numpy as np
from ground_truth_generator import FileReader, OrderBook, Event, parse_feed_line
from typing import List
from datetime import datetime, timedelta

def parse_start_time(header_line):
    """Converts a header line to seconds since midnight (hh:mm:ss -> seconds)."""
    parts = header_line.strip().split()
    time_part = parts[2]  # format hh:mm:ss
    h, m, s = map(int, time_part.split(":"))
    return h * 3600 + m * 60 + s

def parse_file(event_lines):
    """Parses an .in file and reconstructs the order history."""
    orders = {}
    for line in event_lines:
        parts = line.strip().split()
        if not parts:
            continue
        μs_offset = int(parts[0])
        event_type = int(parts[1])
        order_id = parts[2]

        if event_type == 1 and μs_offset >= 0:  # new order
            price = int(parts[3])
            quantity = int(parts[4])
            orders[order_id] = {
                "price": price,
                "quantity": quantity,
                "created_at": μs_offset,
                "modifications": 0,
                "executions": 0,
                "cancelled": False,
                "fully_executed": False
            }
        elif event_type == 2 and order_id in orders:  # modification
            orders[order_id]["price"] = int(parts[3])
            orders[order_id]["quantity"] = int(parts[4])
            orders[order_id]["modifications"] += 1
        elif event_type == 3 and order_id in orders:  # partial execution
            qty = int(parts[3])
            orders[order_id]["quantity"] -= qty
            orders[order_id]["executions"] += 1
        elif event_type == 4 and order_id in orders:  # cancellation
            orders[order_id]["cancelled"] = True
        elif event_type == 5 and order_id in orders:  # full execution
            orders[order_id]["fully_executed"] = True

    return orders

def extract_features(order):
    """Extracts features for a single order as a list."""
    return [
        order["price"],
        order["quantity"],
        order["modifications"],
        order["executions"],
        int(order["cancelled"]),
        int(order["fully_executed"]),
        order["created_at"]
    ]

def process_file(file_path):
    """Processes a single .in file and returns a list of extracted features."""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    event_lines = lines[1:]  # skip header
    orders = parse_file(event_lines)
    features = [extract_features(order) for order in orders.values()]
    return features

def process_directory_parallel(directory_path, n_jobs=-1):
    """Processes all .in files in a directory using parallel execution."""
    files = [os.path.join(directory_path, f) for f in os.listdir(directory_path) if f.endswith('.in')]
    results = joblib.Parallel(n_jobs=n_jobs)(
        joblib.delayed(process_file)(file_path) for file_path in files
    )
    # Flatten the list of lists into a single list of feature vectors
    all_features = [feature for file_features in results for feature in file_features]
    return all_features

def save_features(features, output_path):
    """Saves the features to a CSV file."""
    columns = ["price", "quantity", "modifications", "executions", "cancelled", "fully_executed", "created_at"]
    df = pd.DataFrame(features, columns=columns)
    df.to_csv(output_path, index=False)

def main():
    parser = argparse.ArgumentParser(description="Extract features from .in order files")
    parser.add_argument("directory", type=str, help="Directory containing .in files")
    parser.add_argument("output", type=str, help="Output CSV path")
    parser.add_argument("--n_jobs", type=int, default=-1, help="Number of parallel jobs (default: -1 = all cores)")
    args = parser.parse_args()

    # Extract features in parallel
    features = process_directory_parallel(args.directory, n_jobs=args.n_jobs)
    # Save to CSV
    save_features(features, args.output)

if __name__ == "__main__":
    main()
