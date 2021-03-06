__author__ = 'nikita_kartashov'

import sys
from sys import argv
from os import path, walk, listdir
from itertools import chain
import logging as log
import multiprocessing as mp
from ast import literal_eval

from bg.bg_io import GRIMMReader


PACKAGES_USED = ('graph', 'metrics', 'output')

for package in PACKAGES_USED:
    sys.path.append(path.abspath(package))

from .metric_runner import compare_metric_results, METRICS
from .output.stdout_printer import StdOutPrinter

BLOCK_FILE_NAME = 'blocks.txt'
CORRECT_TREE_FILE_NAME = 'correct_tree.newick'


def run_metrics_on_block_file(block_path, full_correct_tree_file_name):
    correct_tree = read_correct_tree(full_correct_tree_file_name)
    with open(block_path) as block_file:
        breakpoint_graph = GRIMMReader.get_breakpoint_graph(block_file)
        return compare_metric_results(breakpoint_graph, correct_tree)


TREE_NODES = ['A', 'B', 'C', 'D']


def read_correct_tree(full_correct_tree_file_name):
    try:
        with open(full_correct_tree_file_name) as correct_tree_file:
            unparsed_tree = correct_tree_file.readline().strip().strip(';')
            for tree_node in TREE_NODES:
                unparsed_tree.replace(tree_node, "'{0}'".format(tree_node))
            return literal_eval(unparsed_tree)
    except FileNotFoundError:
        log.error('Has not found correct tree file {0}'.format(full_correct_tree_file_name))
        exit(2)


def run_metrics_on_block_folder(block_folder_path):
    for root_path, directory_names, file_names in walk(block_folder_path):
        for block_file_name in (name for name in file_names if name == BLOCK_FILE_NAME):
            full_name = path.join(root_path, block_file_name)
            full_correct_tree_file_name = path.join(root_path, CORRECT_TREE_FILE_NAME)
            yield run_metrics_on_block_file(full_name, full_correct_tree_file_name)


def reduce_run_results(run_results):
    result_sum = [0] * METRICS.metric_number()
    result_number = 0
    for result in run_results:
        for i, e in enumerate(result):
            result_sum[i] += e
        result_number += 1

    return result_sum, result_number


def run_computation_on_folder(block_folder_path):
    run_results = run_metrics_on_block_folder(block_folder_path)
    metric_results, file_count = reduce_run_results(run_results)
    return list(result * 1.0 / file_count for result in metric_results)


def setup_logging():
    root = log.getLogger()
    root.setLevel(log.DEBUG)
    handler = log.FileHandler(path.abspath('out.log'))
    root.addHandler(handler)


def main():
    if len(argv) < 2:
        print('No block folder supplied')
        exit(1)
    input_folder = path.abspath(argv[1])
    if not path.exists(input_folder):
        print("Path {0} doesn't exist".format(input_folder))
        exit(1)

    if not path.isdir(input_folder):
        print("Path {0} is not a directory path".format(input_folder))
        exit(1)

    folder_filterer = lambda _: True

    if len(argv) == 3:
        folder_filterer = lambda folder: folder.startswith(argv[2])

    folder_header = ('run', 'e1', 'e2')

    setup_logging()
    max_width = max(map(len, METRICS.metric_annotations()))

    with StdOutPrinter() as printer:
        printer.write_header(chain(folder_header, METRICS.metric_annotations()), max_width)
        folders_to_work_on = [f for f in listdir(input_folder) if
                              folder_filterer(f) and path.isdir(path.join(input_folder, f))]
        parallel_pool = mp.Pool()
        folder_results = parallel_pool.map(run_computation_on_folder,
                                           [path.join(input_folder, f) for f in folders_to_work_on])
        for folder, folder_result in zip(folders_to_work_on, folder_results):
            printer.write_row(folder.split('_'), folder_result, max_width)
            # log.info('Finished directory {0}'.format(folder))


if __name__ == '__main__':
    main()