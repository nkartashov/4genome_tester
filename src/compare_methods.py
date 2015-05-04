__author__ = 'nikita_kartashov'

from sys import argv
from os import path, walk
from functools import reduce
import csv

from bg.bg_io import GRIMMReader

from metric_runner import compare_metric_results, METRIC_NUMBER, METRIC_ANNOTATIONS


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
    def sum_every_coordinate(left_acc, right):
        left, i = left_acc
        return (l + r for l, r in zip(left, right)), i + 1

    return reduce(sum_every_coordinate, run_results, ([0 for _ in range(METRIC_NUMBER)], 0))


def output_results(results, csv_output_file_path):
    if csv_output_file_path is not None:
        with open(csv_output_file_path, 'w') as csv_output:
            csv_writer = csv.writer(csv_output)
            csv_writer.writerow(METRIC_ANNOTATIONS)
            csv_writer.writerow(list(map(str, results)))
    else:
        print('\t'.join(METRIC_ANNOTATIONS))
        print('\t'.join(map(str, results)))


def run_computation_on_folder(block_folder_path, csv_output_file_path=None, result_handler=output_results):
    run_results = run_metrics_on_block_folder(block_folder_path)
    metric_results, file_count = reduce_run_results(run_results)
    overall_result = list(result * 1.0 / file_count for result in metric_results)
    result_handler(overall_result, csv_output_file_path)


def main():
    if len(argv) < 2:
        print('No block folder supplied')
        exit(1)
    block_folder_path = argv[1]
    csv_output_file_path = None
    if len(argv) == 3:
        csv_output_file_path = argv[2]
    if not path.exists(block_folder_path):
        print("Path {0} doesn't exist".format(block_folder_path))
        exit(1)

    if not path.isdir(block_folder_path):
        print("Path {0} is not a directory path".format(block_folder_path))
        exit(1)
    run_computation_on_folder(block_folder_path, csv_output_file_path)


if __name__ == '__main__':
    main()