__author__ = 'nikita_kartashov'

from collections import Counter, defaultdict

from hashable_edge import HashableEdge
from branch import compute_tree_score_with_branches


ALL_GENOMES = frozenset(['A', 'B', 'C', 'D'])
DEFAULT_TOPOLOGY = (('A', 'B'), ('C', 'D'))


def normalize_multicolor(multicolor):
    first_color = frozenset(multicolor)
    second_color = ALL_GENOMES - first_color
    if len(first_color) != len(second_color):
        if len(first_color) < len(second_color):
            return first_color, second_color
        else:
            return second_color, first_color
    else:
        if first_color < second_color:
            return first_color, second_color
        else:
            return second_color, first_color


def get_distribution_metric(breakpoint_graph, tree_topology=DEFAULT_TOPOLOGY):
    distribution = \
        Counter(normalize_multicolor(edge.multicolor.colors) for edge in breakpoint_graph.edges())
    # Score is negative, so we can compare metrics
    NEGATIVE = -1
    return NEGATIVE * compute_tree_score_with_branches(distribution.items(), tree_topology)


def edge_vertices(edge):
    return [edge.vertex1, edge.vertex2]


def get_simple_paths_metric(breakpoint_graph, tree_topology=DEFAULT_TOPOLOGY):
    result = defaultdict(lambda: 0)

    def is_vertex_simple(vertex):
        return len(list(breakpoint_graph.get_edges_by_vertex(vertex))) == 2

    simple_vertices = set(filter(is_vertex_simple, breakpoint_graph.nodes()))

    def is_edge_simple(edge):
        return all(vertex in simple_vertices for vertex in edge_vertices(edge))

    for multicolor in \
            (normalize_multicolor(edge.multicolor.colors) for edge in breakpoint_graph.edges() if
             is_edge_simple(edge)):
        result[multicolor] += 1

    # Score is negative, so we can compare metrics
    NEGATIVE = -1
    return NEGATIVE * compute_tree_score_with_branches(result.items(), tree_topology)


def get_distance_by_additive_metric(breakpoint_graph, metric, tree_topology=DEFAULT_TOPOLOGY):
    return sum(metric(breakpoint_graph, pair_genomes)
               for pair_genomes in tree_topology)


# Calculates S_BP as in Wei Xu paper, %tree_topology% is a tuple, of tuples of genomes,
# like (('A', 'B'), ('C', 'D')) denotes AB|CD
def get_bp_distance_metric(breakpoint_graph, tree_topology=DEFAULT_TOPOLOGY):
    return int(get_distance_by_additive_metric(
        breakpoint_graph, get_bp_distance_two_genomes, tree_topology))


# Calculates S_DCJ as in Wei Xu paper, %tree_topology% is a tuple, of tuples of genomes,
# like (('A', 'B'), ('C', 'D')) denotes AB|CD
def get_dcj_distance_metric(breakpoint_graph, tree_topology=DEFAULT_TOPOLOGY):
    return int(get_distance_by_additive_metric(
        breakpoint_graph, get_dcj_distance_two_genomes, tree_topology))


# Calculates d_BP between two genomes, %genomes% is a tuple of them
# like ('A', 'B')
def get_bp_distance_two_genomes(breakpoint_graph, genomes):
    block_number = len(list(breakpoint_graph.nodes())) / 2

    def get_adjacencies(genome):
        return frozenset(HashableEdge(edge) for edge in breakpoint_graph.edges() if genome in edge.multicolor.colors)

    left_adjacencies, right_adjacencies = map(get_adjacencies, genomes)
    return block_number - len(left_adjacencies & right_adjacencies)


# Calculates d_DCJ between two genomes, %genomes% is a tuple of them
# like ('A', 'B')
def get_dcj_distance_two_genomes(breakpoint_graph, genomes):
    def is_consistent_with_genomes(edge):
        return any(genome in edge.multicolor.colors for genome in genomes)

    def has_end_unvisited(edge):
        return any(vertex not in visited for vertex in edge_vertices(edge))

    def should_be_visited(edge):
        return is_consistent_with_genomes(edge) and has_end_unvisited(edge)

    block_number = len(list(breakpoint_graph.nodes())) / 2
    connected_components = 0
    visited = set()
    for node in breakpoint_graph.nodes():
        if node in visited:
            continue
        connected_components += 1
        next_node = node

        while next_node is not None:
            vertex = next_node
            next_node = None
            visited.add(vertex)
            for edge in filter(should_be_visited, breakpoint_graph.get_edges_by_vertex(vertex)):
                def unvisited_end():
                    return edge.vertex1 if edge.vertex1 not in visited else edge.vertex2

                next_node = unvisited_end()

    return block_number - connected_components
