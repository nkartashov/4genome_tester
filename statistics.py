__author__ = 'nikita_kartashov'

from collections import Counter, defaultdict

ALL_GENOMES = frozenset(['A', 'B', 'C', 'D'])


def normalize_multicolor(multicolor):
    first_color = frozenset(multicolor)
    second_color = ALL_GENOMES.difference(first_color)
    if first_color < second_color:
        return first_color, second_color
    else:
        return second_color, first_color


def get_distribution_metric(breakpoint_graph):
    return Counter(normalize_multicolor(edge.multicolor.colors) for edge in breakpoint_graph.edges())


def get_simple_paths_metric(breakpoint_graph):
    result = defaultdict(lambda: 0)

    def edge_vertices(edge):
        return [edge.vertex1, edge.vertex2]

    def is_vertex_simple(vertex):
        return len(list(breakpoint_graph.get_edges_by_vertex(vertex))) == 2

    simple_vertices = set(filter(is_vertex_simple, breakpoint_graph.nodes()))

    def is_edge_simple(edge):
        return all(vertex in simple_vertices for vertex in edge_vertices(edge))

    for multicolor in (normalize_multicolor(edge.multicolor.colors) for edge in breakpoint_graph.edges() if
                       is_edge_simple(edge)):
        result[multicolor] += 1

    return result