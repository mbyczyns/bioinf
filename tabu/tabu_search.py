import copy
import random
from enum import Enum
import time
from typing import Optional
from classes_functions import AssemblyTrack, DNAInfo


class Operation(Enum):
    ADD_OLIGO = 1
    REMOVE_OLIGO = 2
    REMOVE_CLUSTER = 3
    MOVE_OLIGO = 4
    MOVE_CLUSTER = 5


class TabuSearch:

    @staticmethod
    def assemble_sequence(ws_track: AssemblyTrack, ry_track: AssemblyTrack) -> str:
        sequence = ""
        for ws_oligo, ry_oligo, depth_val in zip(ws_track.path, ry_track.path, ws_track.depth):
            merged = AssemblyTrack.connect_ws_ry(ws_oligo, ry_oligo)
            sequence += merged[depth_val - len(ws_oligo):]
        return sequence

    @staticmethod
    def calculate_negative_mismatches(ws_track: AssemblyTrack) -> int:
        total = 0
        for depth_val in ws_track.depth[1:]:
            total += ws_track.perfect_overlap - depth_val
        return total

    def __init__(self, max_tabu_length, iterations_limit, neighbours_count):
        self.max_tabu_length = max_tabu_length
        self.tabu_ws_list = []
        self.tabu_ry_list = []
        self.iterations_limit = iterations_limit
        self.neighbours_count = neighbours_count
        self.ignore_tabu_flag = False

    def add_move(self, ws_move, ry_move):
        self.tabu_ws_list.append(ws_move)
        self.tabu_ry_list.append(ry_move)

    def check_tabu(self, candidate_move):
        if self.ignore_tabu_flag:
            return False
        return (candidate_move in self.tabu_ws_list) or (candidate_move in self.tabu_ry_list)

    def neighbour_insert_oligo(
        self,
        ws_track: AssemblyTrack,
        ry_track: AssemblyTrack,
        dna_info: DNAInfo,
        available_ws_oligos,
        available_ry_oligos,
    ) -> Optional[tuple[AssemblyTrack, AssemblyTrack, str, str]]:

        insert_positions = []
        for pos in range(1, len(ws_track.path) + 1):
            if pos != len(ws_track.path) and (
                ws_track.depth[pos - 1] == ws_track.perfect_overlap
                and ws_track.depth[pos] == ws_track.perfect_overlap
            ):
                continue
            insert_positions.append(pos)

        if insert_positions:
            random.shuffle(available_ws_oligos)
            for candidate_ws in available_ws_oligos:
                for candidate_ry in available_ry_oligos:
                    if candidate_ws[-1] == candidate_ry[-1]:
                        for position in insert_positions:
                            new_ws = copy.deepcopy(ws_track)
                            new_ry = copy.deepcopy(ry_track)
                            new_ws.path.insert(position, candidate_ws)
                            new_ry.path.insert(position, candidate_ry)
                            new_ws.update_depth()
                            new_ry.update_depth()

                            if (
                                new_ws.depth == new_ry.depth
                                and len(TabuSearch.assemble_sequence(new_ws, new_ry))
                                <= dna_info.length
                            ):
                                new_ws.cells_dict[candidate_ws] = True
                                new_ry.cells_dict[candidate_ry] = True
                                return (new_ws, new_ry, candidate_ws, candidate_ry)

        return None

    def neighbour_remove_oligo(
        self,
        ws_track: AssemblyTrack,
        ry_track: AssemblyTrack,
        dna_info: DNAInfo,
        skip_tabu=False,
    ) -> Optional[tuple[AssemblyTrack, AssemblyTrack, str, str]]:

        removable_positions = []
        for pos in range(1, len(ws_track.path)):
            if (
                pos != len(ws_track.path) - 1
                and (
                    ws_track.depth[pos] == ws_track.perfect_overlap
                    and ws_track.depth[pos + 1] == ws_track.perfect_overlap
                )
            ) or (self.check_tabu(ws_track.path[pos]) and not skip_tabu):
                continue
            removable_positions.append(pos)

        random.shuffle(removable_positions)
        original_ws = copy.deepcopy(ws_track)
        original_ry = copy.deepcopy(ry_track)

        for pos_to_remove in removable_positions:
            ws_copy = copy.deepcopy(original_ws)
            ry_copy = copy.deepcopy(original_ry)
            removed_ws = ws_copy.path.pop(pos_to_remove)
            removed_ry = ry_copy.path.pop(pos_to_remove)
            ws_copy.update_depth()
            ry_copy.update_depth()

            if ws_copy.depth == ry_copy.depth and len(TabuSearch.assemble_sequence(ws_copy, ry_copy)) <= dna_info.length:
                ws_copy.cells_dict[removed_ws] = False
                ry_copy.cells_dict[removed_ry] = False
                return (ws_copy, ry_copy, removed_ws, removed_ry)
        return None

    def neighbour_shift_oligo(
        self,
        ws_track: AssemblyTrack,
        ry_track: AssemblyTrack,
        dna_info: DNAInfo,
    ) -> Optional[tuple[AssemblyTrack, AssemblyTrack]]:

        candidate_positions = []
        insertion_positions = []

        for pos in range(1, len(ws_track.path)):
            if (
                pos != len(ws_track.path) - 1
                and (
                    ws_track.depth[pos] == ws_track.perfect_overlap
                    and ws_track.depth[pos + 1] == ws_track.perfect_overlap
                )
            ) or self.check_tabu(ws_track.path[pos]):
                continue
            candidate_positions.append(pos)

        random.shuffle(candidate_positions)
        original_ws = copy.deepcopy(ws_track)
        original_ry = copy.deepcopy(ry_track)

        for remove_pos in candidate_positions:
            ws_copy = copy.deepcopy(original_ws)
            ry_copy = copy.deepcopy(original_ry)
            oligo_ws = ws_copy.path.pop(remove_pos)
            oligo_ry = ry_copy.path.pop(remove_pos)
            ws_copy.update_depth()
            ry_copy.update_depth()

            if not insertion_positions:
                for insert_pos in range(1, len(ws_copy.path) + 1):
                    if insert_pos == remove_pos or (
                        insert_pos != len(ws_copy.path)
                        and ws_copy.depth[insert_pos - 1] == ws_copy.perfect_overlap
                        and ws_copy.depth[insert_pos] == ws_copy.perfect_overlap
                    ):
                        continue
                    insertion_positions.append(insert_pos)
                random.shuffle(insertion_positions)

            for insert_pos in insertion_positions:
                ws_copy.path.insert(insert_pos, oligo_ws)
                ry_copy.path.insert(insert_pos, oligo_ry)
                ws_copy.update_depth()
                ry_copy.update_depth()

                if (
                    remove_pos != insert_pos
                    and ws_copy.depth == ry_copy.depth
                    and len(TabuSearch.assemble_sequence(ws_copy, ry_copy)) <= dna_info.length
                ):
                    return (ws_copy, ry_copy)
                ws_copy.path.pop(insert_pos)
                ry_copy.path.pop(insert_pos)

        return None

    def neighbour_remove_cluster(
        self,
        ws_track: AssemblyTrack,
        ry_track: AssemblyTrack,
        dna_info: DNAInfo,
    ) -> Optional[tuple[AssemblyTrack, AssemblyTrack]]:

        ws_copy = copy.deepcopy(ws_track)
        ry_copy = copy.deepcopy(ry_track)

        cluster_start_found = False
        cluster_indices = []

        for i in range(1, len(ws_copy.path) - 1):
            if ws_copy.depth[i] == ws_copy.perfect_overlap:
                if not cluster_start_found and ws_copy.depth[i + 1] == ws_copy.perfect_overlap:
                    cluster_indices.append(i)
                    cluster_start_found = True
            else:
                cluster_start_found = False

        if cluster_indices:
            chosen_cluster_idx = random.choice(cluster_indices)
            while chosen_cluster_idx < len(ws_copy.path) and ws_copy.depth[chosen_cluster_idx] == ws_copy.perfect_overlap:
                removed_ws = ws_copy.path.pop(chosen_cluster_idx)
                removed_ry = ry_copy.path.pop(chosen_cluster_idx)
                ws_copy.depth.pop(chosen_cluster_idx)
                ry_copy.depth.pop(chosen_cluster_idx)
                ws_copy.cells_dict[removed_ws] = False
                ry_copy.cells_dict[removed_ry] = False

            ws_copy.update_depth()
            ry_copy.update_depth()

            if ws_copy.depth == ry_copy.depth and len(TabuSearch.assemble_sequence(ws_copy, ry_copy)) <= dna_info.length:
                return (ws_copy, ry_copy)

        return None

    def neighbour_shift_cluster(
        self,
        ws_track: AssemblyTrack,
        ry_track: AssemblyTrack,
        dna_info: DNAInfo,
    ) -> Optional[tuple[AssemblyTrack, AssemblyTrack]]:

        ws_copy = copy.deepcopy(ws_track)
        ry_copy = copy.deepcopy(ry_track)

        cluster_start = False
        cluster_indices = []

        for i in range(1, len(ws_copy.path) - 1):
            if ws_copy.depth[i] == ws_copy.perfect_overlap:
                if not cluster_start and ws_copy.depth[i + 1] == ws_copy.perfect_overlap:
                    cluster_indices.append(i)
                    cluster_start = True
            else:
                cluster_start = False

        cluster_ws = []
        cluster_ry = []

        if cluster_indices:
            cluster_idx = random.choice(cluster_indices)
            while cluster_idx < len(ws_copy.path) and ws_copy.depth[cluster_idx] == ws_copy.perfect_overlap:
                cluster_ws.append(ws_copy.path.pop(cluster_idx))
                cluster_ry.append(ry_copy.path.pop(cluster_idx))
                ws_copy.depth.pop(cluster_idx)
                ry_copy.depth.pop(cluster_idx)

            ws_copy.update_depth()
            ry_copy.update_depth()

            valid_positions = []
            for pos in range(1, len(ws_copy.path) + 1):
                if (
                    pos != len(ws_copy.path)
                    and ws_copy.depth[pos - 1] == ws_copy.perfect_overlap
                    and ws_copy.depth[pos] == ws_copy.perfect_overlap
                ) or pos == cluster_idx:
                    continue
                valid_positions.append(pos)

            random.shuffle(valid_positions)
            original_ws = copy.deepcopy(ws_copy)
            original_ry = copy.deepcopy(ry_copy)

            for insert_pos in valid_positions:
                ws_copy = copy.deepcopy(original_ws)
                ry_copy = copy.deepcopy(original_ry)
                cluster_ws_copy = cluster_ws.copy()
                cluster_ry_copy = cluster_ry.copy()

                while cluster_ws_copy:
                    ws_copy.path.insert(insert_pos, cluster_ws_copy.pop())
                    ry_copy.path.insert(insert_pos, cluster_ry_copy.pop())

                ws_copy.update_depth()
                ry_copy.update_depth()

                if ws_copy.depth == ry_copy.depth and len(TabuSearch.assemble_sequence(ws_copy, ry_copy)) <= dna_info.length:
                    return (ws_copy, ry_copy)

        return None

    def generate_all_neighbours(
        self,
        ws_track: AssemblyTrack,
        ry_track: AssemblyTrack,
        dna_info: DNAInfo,
    ) -> tuple[tuple[AssemblyTrack, AssemblyTrack] | tuple[AssemblyTrack, AssemblyTrack, str, str], ...]:

        generated_neighbours = []
        seen_paths = {str(ws_track.path)}

        available_ws_oligos = [oligo for oligo in ws_track.not_used_oligos() if not self.check_tabu(oligo)]
        available_ry_oligos = [oligo for oligo in ry_track.not_used_oligos() if not self.check_tabu(oligo)]

        for _ in range(self.neighbours_count):
            chosen_operation = random.choice(list(Operation))
            candidate_neighbour = None

            match chosen_operation:
                case Operation.ADD_OLIGO:
                    candidate_neighbour = self.neighbour_insert_oligo(
                        ws_track,
                        ry_track,
                        dna_info,
                        available_ws_oligos,
                        available_ry_oligos,
                    )
                case Operation.REMOVE_OLIGO:
                    candidate_neighbour = self.neighbour_remove_oligo(ws_track, ry_track, dna_info)
                case Operation.REMOVE_CLUSTER:
                    candidate_neighbour = self.neighbour_remove_cluster(ws_track, ry_track, dna_info)
                case Operation.MOVE_OLIGO:
                    candidate_neighbour = self.neighbour_shift_oligo(ws_track, ry_track, dna_info)
                case Operation.MOVE_CLUSTER:
                    candidate_neighbour = self.neighbour_shift_cluster(ws_track, ry_track, dna_info)

            if candidate_neighbour and str(candidate_neighbour[0].path) not in seen_paths:
                seen_paths.add(str(candidate_neighbour[0].path))
                generated_neighbours.append(candidate_neighbour)

        return tuple(generated_neighbours)

    def score_length(self, ws_track: AssemblyTrack, ry_track: AssemblyTrack):
        return len(ws_track.path)

    def score_compactness(self, ws_track: AssemblyTrack, ry_track: AssemblyTrack):
        return len(ws_track.path) / len(TabuSearch.assemble_sequence(ws_track, ry_track))

    def search_best_solution(
        self,
        ws_track: AssemblyTrack,
        ry_track: AssemblyTrack,
        dna_info: DNAInfo,
        initial_solution: tuple[AssemblyTrack, AssemblyTrack],
    ):

        start_time = time.time()
        current_best = copy.deepcopy(initial_solution)
        fallback_solution: Optional[tuple[AssemblyTrack, AssemblyTrack, int, int]] = None

        for _ in range(self.iterations_limit):
            current_time = time.time()
            if current_time - start_time > 60:
                return (current_best[0], current_best[1])

            best_neighbour_info: Optional[
                tuple[int, tuple[AssemblyTrack, AssemblyTrack] | tuple[AssemblyTrack, AssemblyTrack, str, str]]
            ] = None

            neighbours = self.generate_all_neighbours(current_best[0], current_best[1], dna_info)

            for neighbour in neighbours:
                if len(neighbour) == 2:
                    neighbour_score = self.score_length(neighbour[0], neighbour[1])
                else:
                    neighbour_score = self.score_compactness(neighbour[0], neighbour[1])

                if best_neighbour_info is None or neighbour_score > best_neighbour_info[0]:
                    best_neighbour_info = (neighbour_score, neighbour)

            if best_neighbour_info is None:
                break

            best_score, best_neighbour = best_neighbour_info

            if len(best_neighbour) == 2:
                current_best = best_neighbour
            else:
                fallback_solution = (best_neighbour[0], best_neighbour[1], best_neighbour[2], best_neighbour[3])
                continue

        if fallback_solution is not None:
            ws, ry, ws_oligo, ry_oligo = fallback_solution
            self.add_move(ws_oligo, ry_oligo)
            return (ws, ry)

        return (current_best[0], current_best[1])
