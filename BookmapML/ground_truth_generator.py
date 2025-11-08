from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import Iterator
from argparse import ArgumentParser

LEN_S     = 15 * 60
T1_S      = 180
T2_S      = 300
T3_S      = 600

@dataclass
class Event:
    timestamp: datetime
    microseconds_from_start: int
    type: int
    id: int
    price: int = -1
    number: int = -1
    side: str = ""

class OrderInfo:
    def __init__(self, event: Event) -> None:
        assert event.type == 1
        assert event.side in {"BID", "ASK"}
        self.id = event.id

        self.created_timestamp = event.timestamp
        self.last_modified_by_trader_timestamp = event.timestamp
        self.last_modified_by_creator_timestamp = event.timestamp

        self.is_executed = False
        self.is_cancelled = False

        self.price = event.price
        self.number = event.number
        self.side = event.side
        
        # self.events = [event]

    def add_event(self, event: Event) -> None:
        assert self.id == event.id
        assert self.get_last_timestamp() <= event.timestamp

        # self.events += [event]

        if event.type == 2:
            self.price = event.price
            self.number = event.number
            self.price = event.price
            self.last_modified_by_creator_timestamp = event.timestamp
        elif event.type == 3:
            assert self.number > event.number

            self.number = event.number
            self.last_modified_by_trader_timestamp = event.timestamp
        elif event.type == 4:
            self.is_cancelled = True
            self.last_modified_by_creator_timestamp = event.timestamp
        elif event.type == 5:
            self.is_executed = True
            self.last_modified_by_trader_timestamp = event.timestamp
        else:
            raise Exception(f"invalid event type value: {event.type}")

    def get_last_timestamp(self) -> datetime:
        return max(
            self.last_modified_by_creator_timestamp, 
            self.last_modified_by_trader_timestamp,
        )

    def is_active(self) -> bool:
        return not (self.is_executed or self.is_cancelled)
    
    def get_ground_truth_value(self, observation_start_time: datetime, observation_end_time: datetime) -> int:
        if observation_start_time is not None and self.created_timestamp < observation_start_time: return -1

        if self.last_modified_by_creator_timestamp > observation_end_time: return 0
        if self.is_executed: return -1
        if self.is_cancelled: return 0

        return 1
    

class OrderBook:
    def __init__(self) -> None:
        self.orders: dict[int, OrderInfo] = dict()
        self._are_orders_locked = False
        self.ids = set()

    def add_event(self, event: Event) -> None:
        if event.type == 1:
            assert event.id not in self.orders
            if not self._are_orders_locked:
                assert event.id not in self.ids
                self.ids.add(event.id)

                self.orders[event.id] = OrderInfo(event)
        elif event.type in {2, 3, 4, 5}:
            if not self._are_orders_locked or event.id in self.orders:
                assert event.id in self.orders
                assert self.orders[event.id].is_active()

                self.orders[event.id].add_event(event)

                if event.type in {4, 5} and not self._are_orders_locked: 
                    self.orders.pop(event.id)
        else:
            raise Exception(f"Invalid event type value: {event.type}")
        
    def add_events(self, events: list[Event]) -> None:
        number = 0
        for event in events:
            number += 1
            self.add_event(event)

        print(f"{events[-1].timestamp if len(events) > 0 else ''} Added {number} events!")
        
    def lock_orders(self) -> None:
        self._are_orders_locked = True


def parse_feed_line(line: str) -> Event:
    tokens = line.split()
    microseconds = int(tokens[0])
    timestamp = datetime(1970, 1, 1) + timedelta(microseconds=microseconds)
    type = int(tokens[1])
    id = int(tokens[2])
    if type in {1, 2}:
        price = int(tokens[3])
        number = int(tokens[4])
    elif type in {3}:
        price = -1
        number = int(tokens[3])
    else:
        price = -1
        number = -1

    if type == 1:
        side = tokens[5]
    else:
        side = ""
    
    return Event(timestamp, microseconds, type, id, price, number, side)


class FileReader:
    def __init__(self, file_path: str):
        self._file_path = file_path

    def __enter__(self) -> 'FileReader':
        self._file = open(self._file_path, "r")

        self._cur_line = 0
        self._events_left = int(self._file.readline().split()[0])
        self._next_event = self._read_event()
        self._last_event = self._next_event
        return self
    
    def __exit__(self, exc_type, exc_value, traceback) -> bool:
        if self._file:
            self._file.close()
        if exc_value is not None:
            print(f"Error after reading line {self._cur_line}")
        return False

    def _read_event(self) -> Event:
        if self._events_left == 0:
            return None
        self._cur_line += 1
        self._events_left -= 1
        return parse_feed_line(self._file.readline())
    
    def peek_next_event(self) -> Event:
        if self._next_event is None:
            self._next_event = self._read_event()
        return self._next_event
    
    def get_next_event(self) -> Event:
        self._last_event = self.peek_next_event()
        self._next_event = self._read_event()
        if self._next_event is not None:
            assert self._last_event.timestamp <= self._next_event.timestamp
        return self._last_event
    
    def get_next_event_timestamp(self) -> datetime:
        return self._next_event.timestamp
    
    def get_last_event_timestamp(self) -> datetime:
        return self._last_event.timestamp

    def get_events(self, time_to: datetime) -> Iterator[Event]:
        while self.peek_next_event() is not None and self.peek_next_event().timestamp < time_to:
            next_event = self.get_next_event()
            if next_event is None:
                break
            yield next_event

    def has_more_event(self) -> bool:
        return self.peek_next_event() is not None


def print_ground_truth_file(values: dict[int, list[int]], file_path):
    with open(file_path, "w") as f:
        f.write(f'{len(values)}\n')
        for k, v in values.items():
            f.write(f'{k} {" ".join([str(x) for x in v])}\n')


class GroundTruth:
    def __init__(self):
        self._orders: dict[int, list[int]] = dict()

    def add(self, id: int, value: int) -> None:
        if id not in self._orders:
            self._orders[id] = []
        self._orders[id] += [value]

    def get_values(self, id: int) -> list[int]:
        return self._orders[id]
    
    def print(self, file_path: str) -> None:
        lines = [
            f'{key} {" ".join(str(item) for item in values)}\n'
            for key, values in self._orders.items()
            if values[0] != -1
        ]
        with open(file_path, "w") as f:
            f.write(f'{len(lines)}\n')
            f.writelines(lines)


class GroundTruthGenerator:
    def __init__(self, file_path: str):
        self._file_reader = FileReader(file_path)
        self._buffer: list[Event] = []
        self._order_book = OrderBook()

    def __enter__(self) -> 'GroundTruthGenerator':
        self._file_reader.__enter__()
        self._cur_timestamp = self._file_reader.get_last_event_timestamp()
        return self

    def __exit__(self, *args, **kwargs) -> bool:
        return self._file_reader.__exit__(*args, **kwargs)
    
    def _take_from_buffer(self, time_to: timedelta) -> list[Event]:
        count = 0
        for event in self._buffer:
            if event.timestamp < time_to: 
                count += 1

        result = self._buffer[:count]
        self._buffer = self._buffer[count:]
        return result
    
    def _add_to_buffer(self, time_to: datetime):
        self._buffer += list(self._file_reader.get_events(time_to))

    def get_cur_timestamp(self) -> datetime:
        return self._cur_timestamp
    
    def take(self, delta: timedelta):
        self._cur_timestamp = self._cur_timestamp + delta
        self._add_to_buffer(self._cur_timestamp)

        self._order_book.add_events(self._take_from_buffer(self._cur_timestamp))
        return self._file_reader.has_more_event()

    def generate_ground_truth(
            self, 
            ts: list[timedelta],
            *, 
            observation_start_time: datetime = None
            ) -> GroundTruth:
        self._add_to_buffer(self._cur_timestamp + ts[-1]) # t[1] < t[2] < t[3]
        order_book_copy = deepcopy(self._order_book)
        order_book_copy.lock_orders()
        result = GroundTruth()

        buffer_id = 0
        for t in ts:
            while buffer_id < len(self._buffer) and self._buffer[buffer_id].timestamp < self._cur_timestamp + t:
                order_book_copy.add_event(self._buffer[buffer_id])
                buffer_id += 1

            for order in order_book_copy.orders.values():
                result.add(order.id, order.get_ground_truth_value(observation_start_time, self._cur_timestamp))

        return result


def create_ground_truth_file(
        input_file_path: str, 
        ground_truth_file_path: str,
        skip: timedelta, 
        take: timedelta, 
        t1: timedelta, 
        t2: timedelta, 
        t3: timedelta,
        ):
    
    with GroundTruthGenerator(input_file_path) as generator:
        generator.take(skip)
        ts = generator.get_cur_timestamp()
        generator.take(take)
        generator.generate_ground_truth([t1, t2, t3], observation_start_time=ts).print(ground_truth_file_path)


def main():
    parser = ArgumentParser()
    parser.add_argument("input_path", type=str)
    parser.add_argument("output_path", type=str)
    parser.add_argument("--skip_s", type=int)
    parser.add_argument("--take_s", type=int, default=LEN_S)
    parser.add_argument("--t1_s", type=int, default=T1_S)
    parser.add_argument("--t2_s", type=int, default=T2_S)
    parser.add_argument("--t3_s", type=int, default=T3_S)

    args = parser.parse_args()

    create_ground_truth_file(
        args.input_path,
        args.output_path,
        timedelta(seconds=args.skip_s),
        timedelta(seconds=args.take_s),
        timedelta(seconds=args.t1_s),
        timedelta(seconds=args.t2_s),
        timedelta(seconds=args.t3_s),
    )

if __name__ == '__main__':
    main()