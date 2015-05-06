__author__ = 'nikita_kartashov'

from operator import itemgetter

from .statistics import get_distribution_metric, \
    get_simple_paths_metric, \
    get_bp_distance_metric, \
    get_dcj_distance_metric, \
    get_cylinder_pattern_metric, \
    get_ca_metric, \
    get_mca_metric

A, B, C, D = ['A', 'B', 'C', 'D']

ANNOTATED_METRICS = ((get_distribution_metric, 'Distribution'),
                     (get_simple_paths_metric, 'Simple_Paths'),
                     (get_cylinder_pattern_metric, 'Cylinder_pattern'),
                     (get_bp_distance_metric, 'S_BP'),
                     (get_dcj_distance_metric, 'S_DCJ'),
                     (get_ca_metric, 'S_CA'),
                     (get_mca_metric, 'S_MCA'))

METRICS = tuple(map(itemgetter(0), ANNOTATED_METRICS))
METRIC_NUMBER = len(METRICS)
METRIC_ANNOTATIONS = tuple(map(itemgetter(1), ANNOTATED_METRICS))



TREES = [((A, B), (C, D)),
         ((A, C), (B, D)),
         ((A, D), (C, B))]

RIGHT_TREE = ((A, B), (C, D))


# If we have m methods and n trees then function returns score matrix of m lines and n columns
def run_metrics(breakpoint_graph):
    return (((metric(breakpoint_graph, tree), tree) for tree in TREES) for metric in METRICS)


def compare_metric_results(breakpoint_graph, right_tree=RIGHT_TREE):
    metric_results = run_metrics(breakpoint_graph)

    def decide_if_right(score_tuple):
        scored_trees = list(score_tuple)
        min_score = min(scored_trees)[0]
        trees_with_min_score = list(tree for score, tree in scored_trees if score == min_score)
        return int(len(trees_with_min_score) == 1 and trees_with_min_score[0] == right_tree)

    return (decide_if_right(score_tuple) for score_tuple in metric_results)