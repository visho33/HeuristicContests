import os
import pandas as pd
from datetime import datetime
from collections import defaultdict
from pathlib import Path
from typing import List
import numpy as np

from ground_truth_generator import FileReader, Event  # Make sure this is defined


def extract_features(in_file_path: str, gt_file_path: str):
    # Load ground truth labels
    labels = {}
    with open(gt_file_path, 'r') as f:
        n = int(f.readline())
        for _ in range(n):
            tokens = f.readline().strip().split()
            order_id = int(tokens[0])
            p1, p2, p3 = map(int, tokens[1:])
            labels[order_id] = (p1, p2, p3)

    # Read events from the .in file
    with FileReader(in_file_path) as reader:
        events: List[Event] = list(reader.get_events(datetime.max))

    if not events:
        return pd.DataFrame()  # No data available

    # Session-wide metadata
    start_time = events[0].timestamp
    end_time = events[-1].timestamp
    session_hour = start_time.hour
    session_weekday = start_time.weekday()

    # Global bid/ask prices and quantities
    all_bids = [e.price for e in events if e.type == 1 and e.side == "BID"]
    all_asks = [e.price for e in events if e.type == 1 and e.side == "ASK"]
    bid_qtys = [e.number for e in events if e.type == 1 and e.side == "BID"]
    ask_qtys = [e.number for e in events if e.type == 1 and e.side == "ASK"]

    total_bid_qty = sum(bid_qtys)
    total_ask_qty = sum(ask_qtys)
    bid_proportion = total_bid_qty / (total_bid_qty + total_ask_qty) if (total_bid_qty + total_ask_qty) > 0 else 0

    max_bid = max(all_bids) if all_bids else 0
    min_ask = min(all_asks) if all_asks else 0
    median_bid = np.median(all_bids) if all_bids else 0
    median_ask = np.median(all_asks) if all_asks else 0

    # Group events by order ID
    order_events = defaultdict(list)
    for e in events:
        order_events[e.id].append(e)

    rows = []

    for order_id, evs in order_events.items():
        if order_id not in labels:
            continue  # Skip non-labeled orders

        evs.sort(key=lambda e: e.timestamp)
        create_ev = next((e for e in evs if e.type == 1), None)
        if not create_ev:
            continue

        final_ev = evs[-1]
        price = final_ev.price if final_ev.price != -1 else create_ev.price
        qty = final_ev.number if final_ev.number != -1 else create_ev.number
        side = create_ev.side
        is_bid = 1 if side == "BID" else 0

        timestamps = [e.timestamp for e in evs]
        deltas = [(t2 - t1).total_seconds() for t1, t2 in zip(timestamps, timestamps[1:])]

        modification_count = sum(e.type in {2, 3} for e in evs)
        transaction_count = sum(e.type == 5 for e in evs)  # full executions
        execution_sizes = [e.number for e in evs if e.type == 5]
        partial_execs = [e for e in evs if e.type == 3]

        # Time since last modification or execution
        time_since_last_mod = (end_time - max((e.timestamp for e in evs if e.type in {2, 3}), default=start_time)).total_seconds()
        time_since_last_exec = (end_time - max((e.timestamp for e in evs if e.type in {3, 5}), default=start_time)).total_seconds()

        # Distance to the best price on the same side
        price_distance = 0
        if is_bid and max_bid > 0:
            price_distance = max_bid - price
        elif not is_bid and min_ask > 0:
            price_distance = price - min_ask

        # Construct feature row
        row = {
            "id": order_id,

            # Order price and quantity
            "final_price": price,
            "final_quantity": qty,

            # Time offset from session start
            "offset": (create_ev.timestamp - start_time).total_seconds(),

            # Session info
            "start_time_sec": start_time.hour * 3600 + start_time.minute * 60 + start_time.second,
            "session_weekday": session_weekday,

            # Order activity
            "num_events": len(evs),
            "modification_count": modification_count,
            "transaction_count": transaction_count,
            "avg_time_between_events": np.mean(deltas) if deltas else 0,

            # Execution metrics
            "was_partially_executed": 1 if partial_execs else 0,
            "mean_exec_size": np.mean(execution_sizes) if execution_sizes else 0,
            "modification_ratio": modification_count / max(len(evs), 1),
            "execution_ratio": transaction_count / max(len(evs), 1),

            # Order type info
            "is_bid": is_bid,

            # Price-relative metrics
            "price_above_median_bid": price - median_bid if is_bid else 0,
            "price_below_median_ask": median_ask - price if not is_bid else 0,
            "price_distance_best": price_distance,

            # Time since recent activity
            "time_since_last_mod": time_since_last_mod,
            "time_since_last_exec": time_since_last_exec,

            # Bid-side dominance
            "bid_proportion": bid_proportion,

            # Target labels
            "p1": labels[order_id][0],
            "p2": labels[order_id][1],
            "p3": labels[order_id][2],
        }

        rows.append(row)

    return pd.DataFrame(rows)


if __name__ == "__main__":
    dfs = []
    for fname in os.listdir("train/in/"):
        if fname.endswith(".in"):
            base = fname.replace(".in", "")
            df = extract_features(
                os.path.join("train/in", fname),
                os.path.join("train/gt", base + ".gt")
            )
            dfs.append(df)
            print(f"{fname}: {len(df)} rows")
    df_all = pd.concat(dfs, ignore_index=True)
    df_all.to_csv("train_data3.csv", index=False)

