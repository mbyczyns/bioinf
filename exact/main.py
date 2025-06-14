import copy
import sys
import time
from typing import Optional
import xmltodict

sys.setrecursionlimit(1_000_000_000)

strength_scheme = {"A": "W", "T": "W", "C": "S", "G": "S"}
type_scheme = {"A": "R", "G": "R", "C": "Y", "T": "Y"}

class FragmentGroup:
    def __init__(self, data):
        self.kind = data["@pattern"]
        self.frags = data["cell"]

class DNAInfo:
    def __init__(self, doc):
        meta = doc["dna"]
        self.size = int(meta["@length"])
        self.start = meta["@start"]
        self.ws_fragments = FragmentGroup(meta["probe"][0])
        self.ry_fragments = FragmentGroup(meta["probe"][1])
        base_count = self.size - len(self.ws_fragments.frags[0]) + 1
        self.expected_mismatches = (base_count - len(self.ws_fragments.frags)) + (base_count - len(self.ry_fragments.frags))

class AssemblyTrack:
    def __init__(self, scheme, anchor, frags):
        self.map = scheme
        self.base = self._encode(anchor)
        self.history = [self.base]
        self.depth = [0]
        self.used = {frag: False for frag in frags}
        self.used[self.base] = True

    def _encode(self, seq):
        return ''.join(self.map[nuc] for nuc in seq[:-1]) + seq[-1]

    def total_length(self):
        return sum(len(self.history[0]) - d for d in self.depth)

    @staticmethod
    def blend(ws_seq, ry_seq):
        merged = ""
        for w, r in zip(ws_seq, ry_seq):
            if w + r == "SR": merged += "G"
            elif w + r == "SY": merged += "C"
            elif w + r == "WR": merged += "A"
            elif w + r == "WY": merged += "T"
            elif w == r: merged += w
            else: raise Exception(f"Bad pair: {w}, {r}")
        return merged

def read_file(path):
    with open(path, "r", encoding="utf-8") as f:
        return DNAInfo(xmltodict.parse(f.read().encode("utf-8")))

def overlap_len(a, b, limit):
    for i in range(limit - 1, 0, -1):
        if a.endswith(b[:i]):
            return i
    return 0

def suggest(ws_track, ry_track, edge_log):
    ws_copy = copy.deepcopy(ws_track)
    ry_copy = copy.deepcopy(ry_track)
    found = []

    ws_prev = ws_copy.history[-1]
    ry_prev = ry_copy.history[-1]
    temp_w = ws_prev[:-1] + strength_scheme[ry_prev[-1]]
    temp_r = ry_prev[:-1] + type_scheme[ws_prev[-1]]

    for ws_key in ws_copy.used:
        if not ws_copy.used[ws_key]:
            for ry_key in ry_copy.used:
                if not ry_copy.used[ry_key] and ws_key[-1] == ry_key[-1] and (ws_prev, ws_key) not in edge_log:
                    ovl1 = overlap_len(temp_w, ws_key, len(ws_key))
                    ovl2 = overlap_len(temp_r, ry_key, len(ry_key))
                    val = ovl1 if ovl1 == ovl2 and ovl1 > 0 else 0
                    found.append((ws_key, ry_key, val))

    return tuple(sorted(found, key=lambda t: t[2], reverse=True))

def apply_candidate(options, ws_track, ry_track, info, edge_log):
    ws = copy.deepcopy(ws_track)
    ry = copy.deepcopy(ry_track)

    for w, r, d in options:
        ws.history.append(w)
        ry.history.append(r)
        ws.depth.append(d)
        ry.depth.append(d)

        if ws.total_length() > info.size:
            ws.history.pop(); ws.depth.pop()
            ry.history.pop(); ry.depth.pop()
            continue

        edge_log.append((ws.history[-2], w))
        ws.used[w] = True
        ry.used[r] = True
        return ws, ry

    raise Exception("No candidates left")

def error_metric(ws_track):
    return sum((len(ws_track.base) - 1) - x for x in ws_track.depth[1:])

def explore(ws_track, ry_track, info, results):
    candidates = suggest(ws_track, ry_track, edge_memory)
    if not candidates:
        curr_len = ws_track.total_length()
        results.append((ws_track, ry_track, curr_len))
        if curr_len == info.size and error_metric(ws_track) == info.expected_mismatches:
            return

    try:
        new_ws, new_ry = apply_candidate(candidates, ws_track, ry_track, info, edge_memory)
        explore(new_ws, new_ry, info, results)
        return
    except:
        curr_len = ws_track.total_length()
        if curr_len == info.size and error_metric(ws_track) == info.expected_mismatches:
            results.append((ws_track, ry_track, curr_len))
        elif not results or curr_len > results[0][2]:
            results.clear()
            results.append((ws_track, ry_track, curr_len))

    if len(ws_track.history) == 1:
        return

    ws_copy = copy.deepcopy(ws_track)
    ry_copy = copy.deepcopy(ry_track)
    ws_copy.used[ws_copy.history[-1]] = False
    ry_copy.used[ry_copy.history[-1]] = False
    ws_copy.history.pop(); ws_copy.depth.pop()
    ry_copy.history.pop(); ry_copy.depth.pop()
    explore(ws_copy, ry_copy, info, results)

# -------------- Entry Point --------------

if len(sys.argv) != 2:
    print("Usage: python main.py <file.xml>")
    sys.exit()

input_file = sys.argv[1]
data = read_file(input_file)

t0 = time.time()

ws_main = AssemblyTrack(strength_scheme, data.start, data.ws_fragments.frags)
ry_main = AssemblyTrack(type_scheme, data.start, data.ry_fragments.frags)

edge_memory: list[tuple[str, str]] = []
results: list[tuple[AssemblyTrack, AssemblyTrack, int]] = []

explore(ws_main, ry_main, data, results)

match = None
for ws, ry, length in results:
    if error_metric(ws) == data.expected_mismatches and length == data.size:
        match = (ws, ry, length)
        break
if not match:
    match = sorted(results, key=lambda s: (s[2], -error_metric(s[0])), reverse=True)[0]

output = ""
for w, r, d in zip(match[0].history, match[1].history, match[0].depth):
    piece = AssemblyTrack.blend(w, r)
    output += piece[d - len(w):]

print(output)
print(f"Duration: {time.time() - t0:.3f}s")
