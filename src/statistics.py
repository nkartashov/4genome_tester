__author__ = 'nikita_kartashov'

from collections import Counter, defaultdict
from math import ceil

from bg import Multicolor
from hashable_edge import HashableEdge
from branch import compute_tree_score_with_branches
from breakpoint_graph_extensions import normalize_multicolor, \
    edge_vertices, is_edge_simple, get_vertex_colored_neighbours, get_vertex_sized_neighbours, \
    traverse_node_starting_in_color


ALL_GENOMES = frozenset(['A', 'B', 'C', 'D'])
DEFAULT_TOPOLOGY = (('A', 'B'), ('C', 'D'))
NEGATIVE = -1


def get_distribution_metric(breakpoint_graph, tree_topology):
    distribution = \
        Counter(normalize_multicolor(edge.multicolor.colors, ALL_GENOMES) for edge in breakpoint_graph.edges())
    # Score is negative, so we can compare metrics
    return NEGATIVE * compute_tree_score_with_branches(distribution.items(), tree_topology)


def get_simple_paths_metric(breakpoint_graph, tree_topology):
    result = defaultdict(lambda: 0)

    for multicolor in \
            (normalize_multicolor(edge.multicolor.colors, ALL_GENOMES) for edge in breakpoint_graph.edges() if
             is_edge_simple(breakpoint_graph, edge)):
        result[multicolor] += 1

    # Score is negative, so we can compare metrics

    return NEGATIVE * compute_tree_score_with_branches(result.items(), tree_topology)


def get_distance_by_additive_metric(breakpoint_graph, metric, tree_topology):
    return sum(metric(breakpoint_graph, pair_genomes)
               for pair_genomes in tree_topology)


# Calculates S_BP as in Wei Xu paper, %tree_topology% is a tuple, of tuples of genomes,
# like (('A', 'B'), ('C', 'D')) denotes AB|CD
def get_bp_distance_metric(breakpoint_graph, tree_topology):
    return int(get_distance_by_additive_metric(
        breakpoint_graph, get_bp_distance_two_genomes, tree_topology))


# Calculates S_DCJ as in Wei Xu paper, %tree_topology% is a tuple, of tuples of genomes,
# like (('A', 'B'), ('C', 'D')) denotes AB|CD
def get_dcj_distance_metric(breakpoint_graph, tree_topology):
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
    # Number of connected components when concerning two genomes is equal to the number of alternating cycles
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


def find_cylinder_patterns(breakpoint_graph):
    cylinder_patterns = dict()
    # Iterate on nodes to check all start nodes
    for start_node in breakpoint_graph.nodes():
        # Check every combination of 2 color edge and 1 color edge
        for double_vertex, double_edge in get_vertex_sized_neighbours(breakpoint_graph, start_node, size=2):
            for single_vertex, single_edge in get_vertex_sized_neighbours(breakpoint_graph, start_node, size=1):
                # Take vertices on other side of the edges
                double_edge_color = double_edge.multicolor
                single_edge_color = single_edge.multicolor
                start_to_single_to_double = frozenset(get_vertex_colored_neighbours(breakpoint_graph,
                                                                                    single_vertex,
                                                                                    double_edge_color,
                                                                                    with_edges=False))
                if len(start_to_single_to_double) == 0:
                    continue

                start_to_double_to_new_single = frozenset(get_vertex_sized_neighbours(breakpoint_graph,
                                                                                      double_vertex,
                                                                                      size=1,
                                                                                      with_edges=False))
                if len(start_to_double_to_new_single) == 0:
                    continue

                final_vertices = (final_vertex for final_vertex in
                                  start_to_single_to_double.intersection(start_to_double_to_new_single) if
                                  breakpoint_graph.get_edge_by_two_vertices(final_vertex,
                                                                            double_vertex).multicolor != single_edge_color)

                for pattern in \
                        (frozenset([start_node, single_vertex, double_vertex, final_vertex]) for final_vertex in
                         final_vertices):
                    cylinder_patterns[pattern] = double_edge_color.colors
    return cylinder_patterns


def get_cylinder_pattern_score(breakpoint_graph, topology):
    cylinder_patterns = find_cylinder_patterns(breakpoint_graph)

    def get_score_on_topology(topology, double_color):
        # If both colors are on one side, then give 1 else 0 because
        # cylinder pattern favors double colors on one side
        return int(any(all(c in topology[i] for c in double_color) for i in range(len(topology))))

    return NEGATIVE * sum(get_score_on_topology(topology, double_color) for double_color in cylinder_patterns.values())


def get_size_of_alternating_structures(breakpoint_graph, colors, modifier=lambda x: x,
                                       get_size_of_paths_instead_of_cycles=True):
    colors = tuple(Multicolor(*color) for color in colors)
    color1, color2 = colors

    def alternate_color(color):
        if color == color1:
            return color2
        else:
            return color1

    visited = set()
    resulting_length = 0
    for node in breakpoint_graph.nodes():
        if node in visited:
            continue
        for color in colors:
            resulting_length += modifier(traverse_node_starting_in_color(
                breakpoint_graph, node, color,
                alternator=alternate_color,
                visited=visited, traverse_paths_instead_of_cycles=get_size_of_paths_instead_of_cycles))
    return resulting_length


def get_ca_metric(breakpoint_graph, tree_topology):
    def halver(value):
        return ceil(value * 1.0 / 2)

    return NEGATIVE * get_size_of_alternating_structures(breakpoint_graph, tree_topology, halver)


def get_mca_metric(breakpoint_graph, tree_topology):
    def cycle_specific_halver(value):
        return value / 2 - 1

    ca_score = get_ca_metric(breakpoint_graph, tree_topology)
    cycles_length = get_size_of_alternating_structures(breakpoint_graph, tree_topology, cycle_specific_halver,
                                                       get_size_of_paths_instead_of_cycles=False)
    return ca_score + NEGATIVE * cycles_length


if __name__ == '__main__':
    from bg import BreakpointGraph, Multicolor
    from networkx import MultiGraph

    def test_cylinder():
        multigraph = MultiGraph()
        multigraph.add_nodes_from(range(4))
        breakpoint_graph = BreakpointGraph(multigraph)
        double_color = ['A', 'B']
        breakpoint_graph.add_edge(0, 1, Multicolor(*double_color))
        breakpoint_graph.add_edge(2, 3, Multicolor(*double_color))
        breakpoint_graph.add_edge(1, 2, Multicolor('C'))
        breakpoint_graph.add_edge(0, 3, Multicolor('D'))
        assert (get_cylinder_pattern_score(breakpoint_graph, DEFAULT_TOPOLOGY) == -1)

    def test_paths():
        multigraph = MultiGraph()
        multigraph.add_nodes_from(range(4))
        breakpoint_graph = BreakpointGraph(multigraph)
        first_color = ['A', 'C']
        second_color = ['B', 'D']
        topology = (first_color, second_color)
        breakpoint_graph.add_edge(0, 1, Multicolor(*first_color))
        breakpoint_graph.add_edge(2, 3, Multicolor(*first_color))
        breakpoint_graph.add_edge(1, 2, Multicolor(*second_color))
        breakpoint_graph.add_edge(0, 3, Multicolor(*second_color))
        print(get_size_of_alternating_structures(breakpoint_graph, topology, get_size_of_paths_instead_of_cycles=False))

    test_cylinder()
    test_paths()
