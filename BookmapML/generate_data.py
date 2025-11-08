import argparse
from datetime import datetime, timedelta

MICROSECONDS_IN_SECOND = 1_000_000

class Order:
    def __init__(self, id_, price, qty, side):
        self.id = id_
        self.price = price
        self.qty = qty
        self.side = side
        self.active = True

    def modify(self, price, qty):
        self.price = price
        self.qty = qty

    def execute_partial(self, qty_left):
        self.qty = qty_left

    def cancel(self):
        self.active = False

    def execute_full(self):
        self.active = False

def add_seconds_to_time(hhmmss_str, seconds):
    """Add seconds to a time string in HH:MM:SS format."""
    t = datetime.strptime(hhmmss_str, "%H:%M:%S")
    t += timedelta(seconds=seconds)
    return t.strftime("%H:%M:%S")

def generate_test_in(input_path, output_path, skip_s):
    start_us = skip_s * MICROSECONDS_IN_SECOND
    end_us = (skip_s + 900) * MICROSECONDS_IN_SECOND  # 15 minutes (900s)

    order_book = {}   # id -> Order
    events_in_block = []

    with open(input_path, 'r') as f:
        header = f.readline().strip()
        n_events, day, start_time_str = header.split()

        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            try:
                offset = int(parts[0])
            except ValueError:
                continue

            event_type = int(parts[1])
            order_id = parts[2]

            # --- Process all events before skip_s to simulate the book state ---
            if offset == -1 or offset < start_us:
                if event_type == 1:
                    price = int(parts[3])
                    qty = int(parts[4])
                    side = parts[5]
                    order_book[order_id] = Order(order_id, price, qty, side)
                elif event_type == 2 and order_id in order_book:
                    price = int(parts[3])
                    qty = int(parts[4])
                    order_book[order_id].modify(price, qty)
                elif event_type == 3 and order_id in order_book:
                    qty_left = int(parts[3])
                    order_book[order_id].execute_partial(qty_left)
                elif event_type == 4 and order_id in order_book:
                    order_book[order_id].cancel()
                    del order_book[order_id]
                elif event_type == 5 and order_id in order_book:
                    order_book[order_id].execute_full()
                    del order_book[order_id]
                continue

            # --- Events inside the 15-minute block ---
            if start_us <= offset <= end_us:
                # Normalize offset to start from 0
                parts[0] = str(offset - start_us)
                events_in_block.append(" ".join(parts))
            elif offset > end_us:
                # Events are sorted, so we can stop early
                break

    # --- Prepare live orders at the start of the block as -1 entries ---
    initial_state_lines = []
    for order in order_book.values():
        if order.active:
            initial_state_lines.append(f"-1 1 {order.id} {order.price} {order.qty} {order.side}")

    # --- Construct final output file ---
    all_lines = initial_state_lines + events_in_block
    total_count = len(all_lines)
    new_start_time = add_seconds_to_time(start_time_str, skip_s)

    with open(output_path, 'w') as fout:
        fout.write(f"{total_count} {day} {new_start_time}\n")
        for line in all_lines:
            fout.write(line + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a test-style .in file from a long input file")
    parser.add_argument('--input', required=True, help='Path to the original long .txt file')
    parser.add_argument('--output', required=True, help='Output .in file path')
    parser.add_argument('--skip_s', type=int, required=True, help='Start time (in seconds) of the 15-minute block')

    args = parser.parse_args()
    generate_test_in(args.input, args.output, args.skip_s)
