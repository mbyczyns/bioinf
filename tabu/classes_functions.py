import requests
import xmltodict
import copy
from typing import Optional


class FragmentGroup:
    def __init__(self, json_xml: dict):
        self.pattern = json_xml["@pattern"]
        self.cells = json_xml["cell"]  

    def __repr__(self):
        return f"Pattern: {self.pattern}, Cells: {self.cells[:3]}"


class DNAInfo:

    def __init__(self, json: dict):
        self.length = int(json["dna"]["@length"])
        self.start = json["dna"]["@start"]
        self.ws_probe = FragmentGroup(json["dna"]["probe"][0])
        self.ry_probe = FragmentGroup(json["dna"]["probe"][1])
        sis = self.length - len(self.ws_probe.cells[0]) + 1
        self.sqne = (sis - len(self.ws_probe.cells)) + (sis - len(self.ry_probe.cells))

    def __repr__(self):
        return f"Length: {self.length}, Start: {self.start}, Probes:[ {self.ws_probe}, {self.ry_probe} ]"


class AssemblyTrack:
    @staticmethod
    def connect_ws_ry(oligo_ws: str, oligo_ry: str) -> str:
        connected = ""
        for nucleotide_ws, nucleotide_ry in zip(oligo_ws, oligo_ry):
            temp = nucleotide_ws + nucleotide_ry
            if temp == "SR":
                connected += "G"
            elif temp == "SY":
                connected += "C"
            elif temp == "WR":
                connected += "A"
            elif temp == "WY":
                connected += "T"
            elif nucleotide_ws == nucleotide_ry:
                connected += nucleotide_ws
            else:
                raise ValueError(
                    f"Invalid nucleotides: {nucleotide_ws}, {nucleotide_ry}"
                )

        return connected

    def __init__(self, dict_convertion: dict, oligo: str, cells: dict):
        self.dict_convertion = dict_convertion
        self.start_converted = self.convert_oligo(oligo)  

        self.cells_dict = {}  
        for cell in cells:
            self.cells_dict[cell] = False
        self.cells_dict[self.start_converted] = True

        self.path = [self.start_converted]
        self.depth = [0]
        self.perfect_overlap = len(self.start_converted) - 1

    def update_depth(self) -> None:
        self.depth = [0]
        for i in range(1, len(self.path)):
            self.depth.append(
                check_overlap(
                    self.path[i - 1][:-1] + self.dict_convertion[self.path[i - 1][-1]],
                    self.path[i],
                    len(self.start_converted),
                )
            )

    def convert_oligo(self, oligo: str) -> str:
        half = ""
        for i in range(len(oligo) - 1):
            half += self.dict_convertion[oligo[i]]
        half += oligo[-1]

        return half

    def get_tmp_length_solution(self) -> int:
        return sum([len(self.path[0]) - depth for depth in self.depth])

    def not_used_oligos(self) -> list[str]:
        return [oligo for oligo in self.cells_dict if not self.cells_dict[oligo]]

    def __repr__(self) -> str:
        return f"Start: {self.start_converted} Path: {self.path} Depth: {self.depth}"

nucleotide_to_weak_strong = {
    "A": "W",
    "T": "W",
    "C": "S",
    "G": "S",
}

nucleotide_to_purine_pyrimidine = {
    "A": "R",
    "G": "R",
    "C": "Y",
    "T": "Y",
}


def fetch_test_data(filename: str,) -> DNAInfo:
    file = open(filename, "r", encoding="utf-8")
    content = file.read().encode("utf-8")
    if not content:
        raise requests.exceptions.RequestException("Failed to fetch test data")
    data = xmltodict.parse(content)
    return DNAInfo(data)


def check_overlap(oligo1: str, oligo2: str, probe: int) -> int:
    for offset in range(probe - 1, 0, -1):
        if oligo1[probe - offset :] == oligo2[:offset]:
            return offset
    return 0

def first_nonzero_overlap_pair(ws: AssemblyTrack, ry: AssemblyTrack) -> Optional[tuple[str, str, int]]:
    last_added_path_ws = ws.path[-1]
    last_added_path_ry = ry.path[-1]
    tmp_last_ws = (
        last_added_path_ws[:-1] + nucleotide_to_weak_strong[last_added_path_ws[-1]]
    )
    tmp_last_ry = (
        last_added_path_ry[:-1]
        + nucleotide_to_purine_pyrimidine[last_added_path_ry[-1]]
    )

    candidates: list[tuple[str, str, int]] = list()

    for vertex_ws in ws.cells_dict: 
        if not ws.cells_dict[vertex_ws]:
            for vertex_ry in ry.cells_dict:
                if not ry.cells_dict[vertex_ry]:
                    if (
                        vertex_ws[-1] == vertex_ry[-1]
                    ):
                        overlap_ws, overlap_ry = check_overlap(
                            tmp_last_ws, vertex_ws, len(vertex_ws)
                        ), check_overlap(tmp_last_ry, vertex_ry, len(vertex_ry))
                        if overlap_ws == overlap_ry and overlap_ws > 0:
                            candidates.append(
                                (vertex_ws, vertex_ry, overlap_ws)
                            )  
    try:
        return tuple(sorted(candidates, key=lambda x: x[2], reverse=True))[0]
    except IndexError:
        return None


def greedy(ws: AssemblyTrack, ry: AssemblyTrack, r: DNAInfo) -> tuple[AssemblyTrack, AssemblyTrack]:
    ws = copy.deepcopy(ws)
    ry = copy.deepcopy(ry)

    reconstructed_dna_length = len(ws.start_converted)

    while True:
        pair = first_nonzero_overlap_pair(
            ws, ry
        )  

        if pair is None:
            break
        reconstructed_dna_length += len(pair[0]) - pair[2]
        if reconstructed_dna_length > r.length:
            break

        ws.cells_dict[pair[0]] = True
        ry.cells_dict[pair[1]] = True
        ws.path.append(pair[0])
        ry.path.append(pair[1])
        ws.depth.append(pair[2])
        ry.depth.append(pair[2])

    return ws, ry