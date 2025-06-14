import sys
import time
from tabu import TabuSearch
from classes_functions import (
    fetch_test_data,
    nucleotide_to_weak_strong,
    nucleotide_to_purine_pyrimidine,
    greedy,
    AssemblyTrack,
    DNAInfo,
)


def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: python main.py <xml_file>")
        return
    filename = sys.argv[1]
    dna_info: DNAInfo = fetch_test_data(filename)
    start_time = time.time()

    ws_track = AssemblyTrack(nucleotide_to_weak_strong, dna_info.start, dna_info.ws_probe.cells)
    ry_track = AssemblyTrack(nucleotide_to_purine_pyrimidine, dna_info.start, dna_info.ry_probe.cells)

    tabu_search = TabuSearch(max_tabu_length=3, iterations_limit=100, neighbours_count=100)

    initial_solution = greedy(ws_track, ry_track, dna_info)
    solution = tabu_search.search_best_solution(ws_track, ry_track, dna_info, initial_solution)

    end_time = time.time()

    print(TabuSearch.assemble_sequence(solution[0], solution[1]))
    print("Czas: ", end_time - start_time)


if __name__ == "__main__":
    main()
