__author__ = 'nikita_kartashov'

from sys import argv
from os import path, walk, listdir
from itertools import chain
import logging as log

from bg.bg_io import GRIMMReader

from .metric_runner import compare_metric_results, METRIC_NUMBER, METRIC_ANNOTATIONS
from .output.stdout_printer import StdOutPrinter


def run_metrics_on_block_file(block_path):
    with open(block_path) as block_file:
        breakpoint_graph = GRIMMReader.get_breakpoint_graph(block_file)
        return compare_metric_results(breakpoint_graph)


def run_metrics_on_block_folder(block_folder_path):
    for root_path, directory_names, file_names in walk(block_folder_path):
        for f in file_names:
            full_name = path.join(root_path, f)
            yield run_metrics_on_block_file(full_name)


def reduce_run_results(run_results):
    result_sum = [0] * METRIC_NUMBER
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

    folder_header = ('run', 'e1', 'e2')

    setup_logging()
    max_width = max(map(len, METRIC_ANNOTATIONS))

    with StdOutPrinter() as printer:
        printer.write_header(chain(folder_header, METRIC_ANNOTATIONS), max_width)
        for folder in (f for f in listdir(input_folder) if path.isdir(path.join(input_folder, f))):
            folder_result = run_computation_on_folder(path.join(input_folder, folder))
            printer.write_row(folder.split('_'), folder_result, max_width)
            log.info('Finished directory {0}'.format(folder))


if __name__ == '__main__':
    main()