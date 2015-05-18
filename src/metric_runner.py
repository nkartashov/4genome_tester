__author__ = 'nikita_kartashov'

from src.graph.statistics import get_distribution_metric, \
    get_simple_paths_metric, \
    get_bp_distance_metric, \
    get_dcj_distance_metric, \
    get_ca_metric, \
    get_mca_metric, \
    get_cumulative_metric_batch

from .metrics.metrics import Metrics

A, B, C, D = ['A', 'B', 'C', 'D']

ANNOTATED_SINGLE_METRICS = ((get_distribution_metric, 'D'),  # Distribution
                            (get_simple_paths_metric, 'SP'),  # Simple Paths
                            (get_bp_distance_metric, 'S_BP'),
                            (get_dcj_distance_metric, 'S_DCJ'),
                            (get_ca_metric, 'S_CA'),
                            (get_mca_metric, 'S_MCA'))

ANNOTATED_BATCH_METRICS = ((get_cumulative_metric_batch, 'MCA+Cylinder+Bag'),)

METRICS = Metrics(ANNOTATED_SINGLE_METRICS, ANNOTATED_BATCH_METRICS)

TOPOLOGIES = [((A, B), (C, D)),
              ((A, C), (B, D)),
              ((A, D), (C, B))]

CORRECT_TOPOLOGY = ((A, B), (C, D))


# If we have m methods and n trees then function returns score matrix of m lines and n columns
# def run_metrics(breakpoint_graph):
# return (((metric(breakpoint_graph, topology), topology) for topology in TOPOLOGIES) for metric in METRICS)


def compare_metric_results(breakpoint_graph, right_tree=CORRECT_TOPOLOGY):
    metric_results = METRICS.run_metrics(breakpoint_graph, TOPOLOGIES)

    def decide_if_right(score_tuple):
        scored_trees = list(score_tuple)
        min_score = min(scored_trees)[0]
        trees_with_min_score = list(tree for score, tree in scored_trees if score == min_score)
        return int(len(trees_with_min_score) == 1 and trees_with_min_score[0] == right_tree)

    return (decide_if_right(score_tuple) for score_tuple in metric_results)