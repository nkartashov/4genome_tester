__author__ = 'nikita_kartashov'

from os import path
from sys import argv

from bg.bg_io import GRIMMReader
from scandir import walk

from .statistics import get_dcj_distance_two_genomes


PAIRS_TO_CHECK = ('A', 'Left'), ('B', 'Left'), ('C', 'Right'), ('D', 'Right')


def validate_datasets(root_directory):
    root_directory = path.abspath(root_directory)

    for root, dirs, files in walk(root_directory):
        for f in files:
            block_file_path = path.join(root, f)
            _, e1, e2 = map(int, path.basename(path.dirname(path.dirname(block_file_path))).split('_'))
            with open(block_file_path) as block_file:
                breakpoint_graph = GRIMMReader.get_breakpoint_graph(block_file)
                for genomes in PAIRS_TO_CHECK:
                    leaf_distance = get_dcj_distance_two_genomes(breakpoint_graph, genomes)
                    if leaf_distance != e2:
                        print('Leaf - inner node distance differs in file {0}, expected={1}, real={2}'.
                              format(block_file_path, e2, leaf_distance))
                inner_node_distance = get_dcj_distance_two_genomes(breakpoint_graph, ('Left', 'Right'))
                if inner_node_distance != e1:
                    print('Inner node distance differs in file {0}, expected={1}, real={2}'.
                          format(block_file_path, e1, inner_node_distance))


if __name__ == '__main__':
    if len(argv) != 2:
        print('Need only path to the root dir, containing datasets')
        exit(1)

    root_directory = argv[1]
    validate_datasets(root_directory)

