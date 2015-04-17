__author__ = 'nikita_kartashov'

from sys import argv

from bg.bg_io import GRIMMReader

from iterate_trees import compare_metric_results

if __name__ == '__main__':
    if len(argv) < 2:
        print('No blocks supplied')
        exit(1)
    block_path = argv[1]
    with open(block_path) as block_file:
        breakpoint_graph = GRIMMReader.get_breakpoint_graph(block_file)
        print(compare_metric_results(breakpoint_graph))