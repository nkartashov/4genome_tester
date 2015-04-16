__author__ = 'nikita_kartashov'

from sys import argv

from bg.bg_io import GRIMMReader

from statistics import get_distribution_metric, get_simple_paths_metric


if __name__ == '__main__':
    if len(argv) < 2:
        print('No blocks supplied')
        exit(1)
    block_path = argv[1]
    with open(block_path) as block_file:
        breakpoint_graph = GRIMMReader.get_breakpoint_graph(block_file)
        print(get_distribution_metric(breakpoint_graph))
        print(get_simple_paths_metric(breakpoint_graph))