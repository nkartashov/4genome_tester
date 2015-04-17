__author__ = 'nikita_kartashov'

from statistics import get_distribution_metric, \
    get_simple_paths_metric, \
    get_bp_distance_metric, \
    get_dcj_distance_metric

A, B, C, D = ['A', 'B', 'C', 'D']

METRICS = [get_distribution_metric,
           get_simple_paths_metric,
           get_bp_distance_metric,
           get_dcj_distance_metric]

METRIC_ANNOTATIONS = ['Distribution', 'Simple Paths', 'S_BP', 'S_DCJ']

TREES = [((A, B), (C, D)),
         ((A, C), (B, D)),
         ((A, D), (C, B))]

RIGHT_TREE = ((A, B), (C, D))


# If we have m methods and n trees then function returns score matrix of m lines and n columns
def run_metrics(breakpoint_graph):
    return (((metric(breakpoint_graph, tree), tree) for tree in TREES) for metric in METRICS)


def compare_metric_results(breakpoint_graph, right_tree=RIGHT_TREE):
    metric_results = run_metrics(breakpoint_graph)
    return list(right_tree == min(score_tuple)[1] for score_tuple in metric_results)