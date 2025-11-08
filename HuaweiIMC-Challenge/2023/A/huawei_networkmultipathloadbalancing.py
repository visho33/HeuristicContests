from collections import deque
from dataclasses import dataclass, field, fields
from typing import Any, Deque, Dict, List, Set, Tuple, Optional
from solution_common.message import Message, Request, SwitchStatsInfo
from solution_common.solution import Solution
from random import random, randint, seed as randomseed

@dataclass
class NodeInfo:
    bw_in: int
    bw_out: int
    size: int
    level: int
    """
    Level 0: Server
    Level 1: Access
    Level 2: Aggregate
    Level 3: Core
    """
    node_id: int

@dataclass
class SwitchInfo:
    id: int
    buf_usage: int = 0

    def encode(self) -> List[int]:
        out = []
        for f in fields(self.__class__):
            v = getattr(self, f.name)
            if f.type is int:
                out.append(v)
            elif f.type is List[int]:
                out.append(len(v))
                out.extend(v)
        return out
    
    @classmethod
    def decode(cls, raw: List[int]) -> "SwitchInfo":
        raw = list(raw)
        raw.reverse()
        out = []
        for f in fields(cls):
            if f.type is int:
                out.append(raw.pop())
            elif f.type is List[int]:
                n = raw.pop()
                out.append([raw.pop() for _ in range(n)])
        return cls(*out)
        
@dataclass
class AccessInfo(SwitchInfo):
    request_refill: List[int] = field(default_factory=list)
    sent_up: List[int] = field(default_factory=list)

@dataclass
class AggregateInfo(SwitchInfo):
    request_refill: List[int] = field(default_factory=list)

@dataclass
class CoreInfo(SwitchInfo):
    pass

@dataclass
class ControllerInfo(SwitchInfo):
    pass


class UserSolution(Solution):
    N: int
    K: int
    T: int
    adj: List[List[int]]
    up_dst: List[int]
    dn_dst: List[int]

    nodes_info: List[NodeInfo]
    messages: List[Message]
    nxt: List[Set[int]]
    nxt1: List[int]
    parent: List[Set[int]]

    request_refill: List[int]
    nxt_ask_up: int
    nxt_dst_up: int
    nxt_dst_dn: int

    info_to_send: SwitchInfo


    def __init__(self, node_id: int, bw_in: int, bw_out: int, size: int, level: int, graph: List[List[int]],
                 nodes_info: List[Tuple[int, int, int, int, int]]):
        super().__init__(node_id, bw_in, bw_out, size, level, graph, nodes_info)
        self.nodes_info = [NodeInfo(*nodeinfo) for nodeinfo in nodes_info]

        randomseed(19)

        self.N = len(graph)  # Number of nodes
        self.adj = [[v for v in range(self.N) if graph[u][v]] for u in range(self.N)]
        self.messages = []   # List of active messages in the node
        self.T = 0

        self._last_result = []
        self._request_msgs_received = 0
        self._request_msgs_failed = 0
        self._last_requests = []

        # BFS
        self.nxt = [set() for _ in range(self.N)]  # self.nxt[v] is the set of adjacent nodes that lie in a shortest path to v
        self.parent = [set() for _ in range(self.N)]
        D = [10**9] * self.N
        D[self.node_id] = 0
        queue: Deque[int] = deque()
        queue.append(self.node_id)
        while queue:
            u = queue.popleft()
            for v in self.adj[u]:
                if D[v] >= D[u] + 1:
                    if D[v] > D[u] + 1:
                        D[v] = D[u] + 1
                        queue.append(v)
                    self.parent[v].add(u)
                    if u == self.node_id:
                        self.nxt[v].add(v)
                    else:
                        self.nxt[v] |= self.nxt[u]
        for v, nxt in enumerate(self.nxt):
            for n in nxt:
                assert self.graph[self.node_id][n] == 1
        self.nxt1 = [next(iter(nxtset)) if nxtset else -1 for nxtset in self.nxt]

        self.request_refill = []
        self.nxt_ask_up = 0
        self.nxt_dst_up = 0
        self.nxt_dst_dn = 0

        # Compute half of the K number of the K-fat tree
        self.K = 0
        for u in range(self.N):
            if self.nodes_info[u].level == 3:
                self.K = len(self.adj[u]) // 2
                break
        self.up_dst = [v for v in self.adj[self.node_id] if self.nodes_info[v].level > self.level]
        self.dn_dst = [v for v in self.adj[self.node_id] if self.nodes_info[v].level < self.level]


    # Add request to node and separate it into messages
    def add_request_list(self, request_list: List[Request]) -> None:
        if self.level == 1:
            self.info_to_send = AccessInfo(self.node_id)

        self._request_msgs_received = 0
        self._request_msgs_failed = 0

        self._last_requests = request_list
        for req in request_list:

            # Only accept request if we have enough space in buffer
            if self.size - len(self.messages) < req.data_size:
                self._request_msgs_failed += req.data_size
                continue
            self._request_msgs_received += req.data_size
            
            for i in range(req.data_size):
                msg = Message(self.node_id, -1, req.target_node_id, req.request_id, i, req.begin_time)
                self.messages.append(msg)


    # Receives info packages from neighbors and returns messages to send this round
    def ask_round_solution(self, neighbor_info_list: List[SwitchStatsInfo]) -> List[Message]:
        if self.node_id == -1:
            self.info_to_send = ControllerInfo(self.node_id)
        elif self.level == 2:
            self.info_to_send = AggregateInfo(self.node_id)
        elif self.level == 3:
            self.info_to_send = CoreInfo(self.node_id)

        # Sort by priority 
        def prio(msg: Message):
            age = self.T - msg.request_begin_time
            cnt = sum(1 for msg2 in self.messages if msg2.request_id == msg.request_id)
            limit = 27 if self.N > 200 else 10
            if age >= limit:
                return (1, cnt, msg.request_begin_time)
            else:
                return (0, cnt, msg.request_begin_time)
        self.messages.sort(key=prio)

        # List of messages to send
        to_send: List[Message] = []

        if self.node_id == -1:
            # CONTROLLER

            # Figure out buffer usage of my neighbors
            buf_usage = [0 for u in range(self.N)]
            for neighbor in neighbor_info_list:
                info = SwitchInfo.decode(neighbor.info)
                if info.id != -1:
                    buf_usage[info.id] = info.buf_usage
            
            # if any(buf_usage) and not hasattr(self, "checked_first_fail"):
            #     self.checked_first_fail = True
            #     if any(failed and self.nodes_info[u].level == 3 for u, failed in enumerate(buf_usage)):
            #         self.kill = True

            # if any(buf_usage[u] >= self.nodes_info[u].size and self.nodes_info[u].level == 1 for u in range(self.N)):
            #     self.kill = True

            # buf_dn = sum(buf_usage[u] for u in range(self.N))
            # if buf_dn >= bs_limit:
            #     self.kill = True

            return []
        elif self.level == 1:
            # ACCESS SWITCH

            # Prefer sending up as many messages as possible
            sent_to = [0 for u in range(self.N)]
            for msg in self.messages:
                if len(to_send) < self.bw_out and self.goes_up(msg):
                    choices = []
                    for v1 in self.up_dst:
                        core = next(core for core in self.adj[v1] if self.nodes_info[core].level == 3)
                        for aggr2 in self.parent[next(iter(self.parent[msg.target_node_id]))]:
                            if self.graph[core][aggr2]:
                                v2 = aggr2
                                break
                        else:
                            raise RuntimeError("aggr2 not found")
                        if sent_to[v1] >= self.nodes_info[v1].bw_out // self.K or sent_to[v2] >= self.nodes_info[v2].bw_out // self.K:
                            continue
                        w = min(self.nodes_info[v1].bw_out, self.nodes_info[v2].bw_out)
                        choices.append(((v1, v2), w ** 2.5))
                    if choices:
                        v1, v2 = weighted_random(choices)
                        msg.to_node_id = v1
                        to_send.append(msg)
                        sent_to[v1] += 1
                        if v1 != v2:
                            sent_to[v2] += 1
            self.info_to_send.sent_up = [sent_to[v] for v in self.up_dst]

            # try:
            #     if len(to_send) / min(sum(self.goes_up(msg) for msg in self.messages), self.bw_out) <= bs_limit:
            #         self.kill = True
            # except ZeroDivisionError:
            #     pass

            # If there is bandwidth remaining, send down
            buffer_dn = 0
            for msg in self.messages:
                if self.goes_down(msg):
                    if len(to_send) < self.bw_out:
                        assert self.graph[self.node_id][msg.target_node_id]
                        msg.to_node_id = msg.target_node_id
                        to_send.append(msg)
                    else:
                        buffer_dn += 1
            buffer_dn += sum(self.request_refill)
            max_buffer_dn = 3.4*self.bw_out if self.T <= 180 else self.size

            # If next tick we may not have enough messages to fill bandwidth, request a buffer refill
            up_bandwidths = [self.nodes_info[v].bw_out for v in self.up_dst]
            max_can_recv = sum(up_bandwidths) // self.K
            request_refill = round(max(min(max_buffer_dn - buffer_dn, max_can_recv), 0) * 1.8)
            request_refill = split_proportionally(request_refill, [round(u ** 1.2) for u in up_bandwidths])
            self.request_refill = request_refill
            self.info_to_send.request_refill = request_refill
        elif self.level == 2:
            # AGGREGATION SWITCH
            our_ids = [u for u in self.adj[self.dn_dst[0]] if self.nodes_info[u].level == self.level]
            requested_to_send_dn = {id: 0 for id in self.dn_dst}
            dn_sent_up = {id: {u: 0 for u in our_ids} for id in self.dn_dst}
            for neighbor in neighbor_info_list:
                id = neighbor.info[0]
                if id == -1:
                    pass
                elif self.nodes_info[id].level == 1:
                    # From access switch
                    info: AccessInfo = AccessInfo.decode(neighbor.info)
                    # requested_to_send_dn[id] = sum(neighbor.info[2:2+self.K])
                    requested_to_send_dn[id] = info.request_refill[our_ids.index(self.node_id)]
                    dn_sent_up[id] = {our_ids[i]: info.sent_up[i] for i in range(self.K)}

            # Split the amount to send down
            for id in self.dn_dst:
                # refill_sum = requested_to_send_dn[id]
                # proportions = [self.nodes_info[u].bw_out - dn_sent_up[id][u] for u in our_ids]
                # refill_by_aggr = split_proportionally(refill_sum, proportions, self.dn_dst.index(id))
                # requested_to_send_dn[id] = refill_by_aggr[our_ids.index(self.node_id)]
                pass
            
            # Send messages up to core
            for msg in self.messages:
                if self.goes_up(msg) and len(to_send) < self.bw_out:
                    msg.to_node_id = self.up_dst[self.nxt_dst_up]
                    self.nxt_dst_up = (self.nxt_dst_up + 1) % self.K
                    to_send.append(msg)

            # Send messages down if requested
            sent_down = 0
            pending_dn: Dict[int, List[Message]] = {id: [] for id in self.dn_dst}
            for msg in reversed(self.messages):
                nxt_id = self.nxt1[msg.target_node_id]
                if nxt_id in pending_dn:
                    pending_dn[nxt_id].append(msg)
            idle = 0
            while idle < self.K + 2:
                id = self.dn_dst[self.nxt_dst_dn]
                self.nxt_dst_dn = (self.nxt_dst_dn + 1) % self.K
                if len(to_send) < self.bw_out and pending_dn[id] and requested_to_send_dn[id]:
                    msg = pending_dn[id].pop()
                    msg.to_node_id = id
                    requested_to_send_dn[id] -= 1
                    to_send.append(msg)
                    sent_down += 1
                    idle = 0
                idle += 1

            # Count messages in the down buffer
            buffer_dn = 0
            for msg in self.messages:
                if self.goes_down(msg):
                    buffer_dn += 1
            buffer_dn -= sent_down
            buffer_dn += sum(self.request_refill)
            max_buffer_dn = self.size - self.bw_out

            # If the buffer is overfull, send any excess down messages to the core
            # This may happen if access switches send us too many messages that go directly down to another access switch
            sent_id_set: Set[Tuple[int, int]] = {(msg.request_id, msg.message_id) for msg in to_send}
            unsent_dn = [msg for msg in reversed(self.messages) if (msg.request_id, msg.message_id) not in sent_id_set and self.goes_down(msg)]
            while buffer_dn > max_buffer_dn and len(to_send) < self.bw_out and unsent_dn:
                msg = unsent_dn.pop()
                msg.to_node_id = self.up_dst[self.nxt_dst_up]
                self.nxt_dst_up = (self.nxt_dst_up + 1) % self.K
                to_send.append(msg)
                buffer_dn -= 1

            # Request messages from core switches to keep down buffer full
            up_bandwidths = [self.nodes_info[v].bw_out for v in self.up_dst]
            max_can_recv = sum(up_bandwidths) // self.K
            request_refill = round(max(min(max_buffer_dn - buffer_dn, max_can_recv), 0) * .89)
            request_refill = split_proportionally(request_refill, [u + self.T for u in up_bandwidths])
            self.request_refill = request_refill
            self.info_to_send.request_refill = request_refill
        elif self.level == 3:
            # CORE SWITCH
            requested_to_send_dn = {id: 0 for id in self.adj[self.node_id] if self.nodes_info[id].level == 2}
            for neighbor in neighbor_info_list:
                id = neighbor.info[0]
                if id == -1:
                    pass
                elif self.nodes_info[id].level == 2:
                    # Refill request from an aggregate switch
                    info: AggregateInfo = AggregateInfo.decode(neighbor.info)
                    dn_neighbors = [v for v in self.adj[id] if self.nodes_info[v].level == self.level]
                    my_index = dn_neighbors.index(self.node_id)
                    requested_to_send_dn[id] = info.request_refill[my_index]

            # Send refill messages down
            pending_dn: Dict[int, List[Message]] = {id: [] for id in self.dn_dst}
            for msg in reversed(self.messages):
                nxt_id = self.nxt1[msg.target_node_id]
                if nxt_id in pending_dn:
                    pending_dn[nxt_id].append(msg)
            idle = 0
            while idle < 2*self.K + 2:
                id = self.dn_dst[self.nxt_dst_dn]
                self.nxt_dst_dn = (self.nxt_dst_dn + 1) % (2*self.K)
                if len(to_send) < self.bw_out and pending_dn[id] and requested_to_send_dn[id]:
                    msg = pending_dn[id].pop()
                    msg.to_node_id = id
                    requested_to_send_dn[id] -= 1
                    to_send.append(msg)
                    idle = 0
                idle += 1

        # print(f"buffer usage for node {self.node_id} at t = {self.T}: {len(self.messages)}/{self.size}, {len(to_send)} to be sent now")

        # Remove sent messages from self.messages
        ans_id_set = set((msg.request_id, msg.message_id) for msg in to_send)
        self.messages = [msg for msg in self.messages if (msg.request_id, msg.message_id) not in ans_id_set]

        for msg in to_send:
            msg.from_node_id = self.node_id
            assert self.graph[msg.from_node_id][msg.to_node_id] == 1

        return to_send
        


    # Results of round and calculate info package to neighbors
    def next_round(self, result: List[Tuple[Message, bool]]) -> SwitchStatsInfo:
        # Add received messages and failed sends to self.messages
        result.sort(key=lambda mo: (mo[0].request_id, mo[0].message_id))
        for mes, ok in result:
            if ok:
                if mes.to_node_id == self.node_id:
                    self.messages.append(mes)
            else:
                if mes.from_node_id == self.node_id:
                    self.messages.append(mes)
                    print("physically failed to send message (%s, %s) from node %s to node %s" % (mes.request_id, mes.message_id, self.node_id, mes.to_node_id))
                    # if self.nodes_info[mes.to_node_id].level == 3 and self.nodes_info[mes.from_node_id].level == 2:
                    #     info = SwitchStatsInfo()
                    #     info.info.append(int(1e10))
                    #     return info
        self._last_result = result
        
        self.info_to_send.buf_usage = len(self.messages)

        # Get info to send to neighbors
        ans = SwitchStatsInfo()
        ans.info = self.info_to_send.encode()

        # Finally, next round
        self.T += 1

        if hasattr(self, "kill"):
            info = SwitchStatsInfo()
            info.info = [int(1e10)]
            return info

        return ans

    # Extra data to place in the log
    def get_extra_log(self) -> Dict:
        return {
            "extra": {
                "usedBuf": len(self.messages),
                "usedOutUp": f"{sum(ok and self.nodes_info[msg.to_node_id].level > self.level for msg, ok in self._last_result)}/{sum(self.nodes_info[msg.to_node_id].level > self.level for msg, ok in self._last_result)}",
                "usedOutDn": f"{sum(ok and self.nodes_info[msg.to_node_id].level < self.level for msg, ok in self._last_result)}/{sum(self.nodes_info[msg.to_node_id].level < self.level for msg, ok in self._last_result)}",
                "usedInUp": sum(self.nodes_info[msg.from_node_id].level > self.level for msg, ok in self._last_result),
                "usedInDn": sum(self.nodes_info[msg.from_node_id].level < self.level for msg, ok in self._last_result),
                "reqRecv": self._request_msgs_received,
                "reqFail": self._request_msgs_failed,
            },
        }
    
    def goes_down(self, msg: Message) -> bool:
        return self.nodes_info[self.nxt1[msg.target_node_id]].level < self.level
    
    def goes_up(self, msg: Message) -> bool:
        return self.nodes_info[self.nxt1[msg.target_node_id]].level > self.level


def weighted_random(candidates: List[Tuple[Any, float]]) -> Optional[Any]:
    """
    Sample a random integer from `candidates`, using the `float` value as a random weight.
    """
    total = sum(weight for cand, weight in candidates)
    r = random() * total
    acc = 0
    for cand, weight in candidates:
        acc += weight
        if r < acc:
            return cand
    # No nonzero-weight candidates
    return None

def ffs(x: int) -> int:
    return (x & -x).bit_length() - 1

def split_proportionally(x: int, prop: List[int], i: Optional[int] = None) -> List[int]:
    total = sum(prop)
    if total == 0:
        y = [x // len(prop) for p in prop]
    else:
        y = [x * p // total for p in prop]
    x -= sum(y)
    if i is None:
        i = randint(0, len(prop) - 1)
    while x:
        y[i] += 1
        x -= 1
        i = (i+1) % len(prop)
    return y
