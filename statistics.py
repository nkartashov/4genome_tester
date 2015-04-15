__author__ = 'nikita_kartashov'

from collections import Counter


def get_distribution(breakpoint_graph):
    return Counter(frozenset(edge.multicolor.colors) for edge in breakpoint_graph.edges())