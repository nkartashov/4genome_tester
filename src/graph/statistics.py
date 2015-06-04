__author__ = 'nikita_kartashov'

from collections import Counter, defaultdict
from math import ceil
from itertools import chain
from functools import reduce

from src.graph.cached_statistic import CachedStatistic
from src.graph.hashable_edge import HashableEdge
from src.graph.branch import compute_tree_score_with_branches
from src.graph.breakpoint_graph_extensions import multicolor_to_normalized_split, \
    edge_vertices, is_edge_simple, get_vertex_colored_neighbours, get_vertex_sized_neighbours, \
    traverse_node_starting_in_color, get_vertex_coloured_sized_neighbours


A = 'A'
B = 'B'
C = 'C'
D = 'D'
ALL_GENOMES = frozenset([A, B, C, D])
NEGATIVE = -1


def get_distribution_metric(breakpoint_graph, tree_topology):
    """
    Finds the distribution metric value of a given topology assuming given BP graph
    :param breakpoint_graph: given BP graph
    :param tree_topology: given topology in the form (('A', 'B'), ('C', 'D'))
    :return: minimized score (i.e. smaller the score the better)
    """
    distribution = Counter(
        multicolor_to_normalized_split(edge.multicolor.colors, ALL_GENOMES)
        for edge in breakpoint_graph.edges())
    # Score is negative, so we can compare metrics
    return NEGATIVE * compute_tree_score_with_branches(distribution.items(), tree_topology)


def get_simple_paths_metric(breakpoint_graph, tree_topology):
    """
    Finds the simple paths metric value of a given topology assuming given BP graph
    :param breakpoint_graph: given BP graph
    :param tree_topology: given topology in the form (('A', 'B'), ('C', 'D')) which denotes AB|CD
    :return: minimized score (i.e. smaller the score the better)
    """
    result = defaultdict(lambda: 0)

    for multicolor in \
            (multicolor_to_normalized_split(edge.multicolor.colors, ALL_GENOMES)
             for edge in breakpoint_graph.edges() if is_edge_simple(breakpoint_graph, edge)):
        result[multicolor] += 1

    # Score is negative, so we can compare metrics
    return NEGATIVE * compute_tree_score_with_branches(result.items(), tree_topology)


def get_distance_by_additive_metric(breakpoint_graph, metric, tree_topology):
    """
    Uses an additive metric to score both pairs comprising a topology, adds them, yielding
     a score for the whole topology
    :param breakpoint_graph: given BP graph to score against
    :param metric: additive metric
    :param tree_topology: topology to be scored in the form (('A', 'B'), ('C', 'D')) which denotes AB|CD
    :return: sum of the scores on the pairs of genomes
    """
    return sum(metric(breakpoint_graph, pair_genomes)
               for pair_genomes in tree_topology)


def get_bp_distance_metric(breakpoint_graph, tree_topology):
    """
    Computes S_BP as per Wei Xu paper as an additive metric
    :param breakpoint_graph: given BP graph to score against
    :param tree_topology: topology to be scored in the form (('A', 'B'), ('C', 'D')) which denotes AB|CD
    :return: minimized score (i.e. smaller the score the better)
    """
    return int(get_distance_by_additive_metric(
        breakpoint_graph, get_bp_distance_two_genomes, tree_topology))


def get_dcj_distance_metric(breakpoint_graph, tree_topology):
    """
    Computes S_DCJ as per Wei Xu paper as an additive metric
    :param breakpoint_graph: given BP graph to score against
    :param tree_topology: topology to be scored in the form (('A', 'B'), ('C', 'D')) which denotes AB|CD
    :return: minimized score (i.e. smaller the score the better)
    """
    return int(get_distance_by_additive_metric(
        breakpoint_graph, get_dcj_distance_two_genomes, tree_topology))


def get_bp_distance_two_genomes(breakpoint_graph, genomes: tuple):
    """
    Computes d_BP between two genomes
    :param breakpoint_graph: given BP graph to score against
    :param genomes: tuple of genome names like ('A', 'B')
    :return: BP score
    """
    block_number = len(list(breakpoint_graph.nodes())) / 2

    def get_adjacencies(genome):
        return frozenset(HashableEdge(edge) for edge in breakpoint_graph.edges() if genome in edge.multicolor.colors)

    left_adjacencies, right_adjacencies = map(get_adjacencies, genomes)
    return block_number - len(left_adjacencies & right_adjacencies)


def get_dcj_distance_two_genomes(breakpoint_graph, genomes):
    """
    Computes d_DCJ between two genomes
    :param breakpoint_graph: given BP graph to score against
    :param genomes: tuple of genome names like ('A', 'B')
    :return: DCJ score
    """

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


def get_size_of_alternating_structures(breakpoint_graph, colors, modifier=lambda x: x,
                                       get_size_of_paths_instead_of_cycles=True):
    """
    Looks for structures (paths or cycles) with alternating colours, applies modifier on their length then sums them
    :param breakpoint_graph: given BP graph to look in
    :param colors: pair of alternating colors
    :param modifier: function to transform the length of the found structures
    :param get_size_of_paths_instead_of_cycles: True if you look for paths, otherwise you look for cycles
    :return: sum of the modified lengths of found structures
    """
    colors = tuple(Multicolor(*color) for color in colors)
    color1, color2 = colors

    def alternate_color(color):
        return color2 if color == color1 else color1

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


@CachedStatistic
def get_ca_metric(breakpoint_graph, tree_topology):
    def halver(value):
        return ceil(value * 1.0 / 2)

    return NEGATIVE * get_size_of_alternating_structures(breakpoint_graph, tree_topology, halver)


@CachedStatistic
def get_mca_metric(breakpoint_graph, tree_topology):
    def cycle_specific_halver(value):
        return value / 2 - 1

    ca_score = get_ca_metric(breakpoint_graph, tree_topology)
    cycles_length = get_size_of_alternating_structures(breakpoint_graph, tree_topology, cycle_specific_halver,
                                                       get_size_of_paths_instead_of_cycles=False)
    return ca_score + NEGATIVE * cycles_length


def get_mca_metric_batch(breakpoint_graph, topologies):
    return ((get_mca_metric(breakpoint_graph, topology), topology) for topology in topologies)


def find_cylinder_patterns(breakpoint_graph):
    """
    Looks for cylinder patterns in the given BP graph. Cylinder is two coloured edges on top,
    two with the same colours on the bottom, two with different colors on the sides
    :param breakpoint_graph: given BP graph to look in
    :return: dictionary, keys are sets of vertices, values are double colours on top & bottom
    """
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


def find_bag_patterns(breakpoint_graph):
    """
    Looks for bag patterns in the given BP graph. Bag is two coloured edges on top,
    edge with one of the colours from the top on the bottom, two with different colors on the sides
    :param breakpoint_graph: given BP graph to look in
    :return: dictionary, keys are sets of vertices, values are double colours on top & bottom
    """
    bag_patterns = dict()
    # Iterate on nodes to check all start nodes
    for start_node in breakpoint_graph.nodes():
        # Check every combination of 2 color edge and 1 color edge
        for double_vertex, double_edge in get_vertex_sized_neighbours(breakpoint_graph, start_node, size=2):
            for single_vertex, single_edge in get_vertex_sized_neighbours(breakpoint_graph, start_node, size=1):
                # Take vertices on other side of the edges
                double_edge_color = double_edge.multicolor
                single_edge_color = single_edge.multicolor
                for final_vertex, final_edge in \
                        get_vertex_coloured_sized_neighbours(breakpoint_graph,
                                                             single_vertex,
                                                             size=1,
                                                             colors=double_edge_color.colors):

                    final_vertex_neighbours = get_vertex_sized_neighbours(breakpoint_graph,
                                                                          final_vertex,
                                                                          size=1,
                                                                          with_edges=False)
                    if double_vertex in final_vertex_neighbours and \
                                    breakpoint_graph.get_edge_by_two_vertices(double_vertex,
                                                                              final_vertex).multicolor != single_edge_color:
                        bag_patterns[frozenset(
                            [start_node, single_vertex, double_vertex, final_vertex])] = double_edge_color.colors
    return bag_patterns


def find_diamond_patterns(breakpoint_graph):
    diamond_patterns = dict()
    # Iterate on nodes to check all start nodes
    for start_node in breakpoint_graph.nodes():
        # Check every combination of 2 color edge and 1 color edge
        for first_vertex, first_edge in get_vertex_sized_neighbours(breakpoint_graph, start_node, size=1):
            for second_vertex, second_edge in ((v, e) for v, e in
                                               get_vertex_sized_neighbours(breakpoint_graph, start_node, size=1)
                                               if e.multicolor != first_edge.multicolor):
                for third_vertex, third_edge in ((v, e) for v, e in
                                                 get_vertex_sized_neighbours(breakpoint_graph, second_vertex, size=1)
                                                 if e.multicolor != first_edge.multicolor and
                                                                 e.multicolor != second_edge.multicolor):
                    last_edge = breakpoint_graph.get_edge_by_two_vertices(third_vertex, first_vertex)
                    if last_edge.multicolor != first_edge.multicolor and \
                                    last_edge.multicolor != second_edge.multicolor and \
                                    last_edge.multicolor != third_edge.multicolor:
                        diamond_patterns[frozenset([start_node, first_vertex, second_vertex, third_vertex])] = (
                            first_edge.multicolor + second_edge.multicolor).colors
    return diamond_patterns


def get_score_on_topology_favouring(topology, colour_to_favour):
    """
    Get score of 1 if two colors that are favoured are on one side
    :param topology: given topology in the form (('A', 'B'), ('C', 'D'))
    :param colour_to_favour: pair of colours (to be together)
    :return: 0 if colors are on different sides, 1 otherwise
    """
    return int(any(all(c in topology[i] for c in colour_to_favour) for i in range(len(topology))))


def get_cylinder_pattern_metric_batch(breakpoint_graph, topologies):
    """
    Look for cylinder patterns, then score every topology from given by using favouring colours
    :param breakpoint_graph: given BP graph to score against
    :param topologies: tuple of topologies, each of which is in the form (('A', 'B'), ('C', 'D'))
    :return: tuple of topologies with scores
    """
    cylinder_patterns = find_cylinder_patterns(breakpoint_graph)
    return ((NEGATIVE *
             sum(get_score_on_topology_favouring(topology, double_color)
                 for double_color in cylinder_patterns.values()), topology)
            for topology in topologies)


def get_bag_pattern_metric_batch(breakpoint_graph, topologies):
    """
    Look for bag patterns, then score every topology from given by using favouring colours
    :param breakpoint_graph: given BP graph to score against
    :param topologies: tuple of topologies, each of which is in the form (('A', 'B'), ('C', 'D'))
    :return: tuple of topologies with scores
    """
    bag_patterns = find_bag_patterns(breakpoint_graph)
    return ((NEGATIVE *
             sum(get_score_on_topology_favouring(topology, double_color)
                 for double_color in bag_patterns.values()), topology)
            for topology in topologies)


def get_diamond_pattern_metric_batch(breakpoint_graph, topologies):
    """
    Look for diamond patterns, then score every topology from given by using favouring colours
    :param breakpoint_graph: given BP graph to score against
    :param topologies: tuple of topologies, each of which is in the form (('A', 'B'), ('C', 'D'))
    :return: tuple of topologies with scores
    """
    diamond_patterns = find_diamond_patterns(breakpoint_graph)
    return ((NEGATIVE *
             sum(get_score_on_topology_favouring(topology, double_color)
                 for double_color in diamond_patterns.values()), topology)
            for topology in topologies)


PATTERN_METRICS = (get_cylinder_pattern_metric_batch, get_bag_pattern_metric_batch, get_diamond_pattern_metric_batch)


def transpose(l):
    return tuple(zip(*l))


def get_cumulative_metric_batch(breakpoint_graph, topologies):
    mca = get_mca_metric_batch(breakpoint_graph, topologies)
    pattern_metrics = (metric(breakpoint_graph, topologies) for metric in PATTERN_METRICS)

    def reducer(acc, new_value):
        return acc[0] + new_value[0], acc[1]

    return (reduce(reducer, topology_row) for topology_row in
            transpose(tuple(map(tuple, chain((mca,), pattern_metrics)))))


if __name__ == '__main__':
    from bg import BreakpointGraph, Multicolor
    from networkx import MultiGraph

    def test_cylinder():
        multigraph1 = MultiGraph()
        multigraph1.add_nodes_from(range(4))
        breakpoint_graph1 = BreakpointGraph(multigraph1)
        double_color = ['A', 'B']
        breakpoint_graph1.add_edge(0, 1, Multicolor(*double_color))
        breakpoint_graph1.add_edge(2, 3, Multicolor(*double_color))
        breakpoint_graph1.add_edge(1, 2, Multicolor('C'))
        breakpoint_graph1.add_edge(0, 3, Multicolor('D'))
        assert (len(find_cylinder_patterns(breakpoint_graph1)) == 1)
        multigraph2 = MultiGraph()
        multigraph2.add_nodes_from(range(4))
        breakpoint_graph2 = BreakpointGraph(multigraph2)
        double_color = ['A', 'B']
        breakpoint_graph2.add_edge(0, 1, Multicolor(*double_color))
        breakpoint_graph2.add_edge(2, 3, Multicolor(*double_color))
        breakpoint_graph2.add_edge(1, 2, Multicolor('C'))
        breakpoint_graph2.add_edge(0, 3, Multicolor('C'))
        assert (len(find_cylinder_patterns(breakpoint_graph2)) == 0)

    def test_bag():
        multigraph1 = MultiGraph()
        multigraph1.add_nodes_from(range(4))
        breakpoint_graph1 = BreakpointGraph(multigraph1)
        double_color = ['A', 'B']
        breakpoint_graph1.add_edge(0, 1, Multicolor(*double_color))
        breakpoint_graph1.add_edge(2, 3, Multicolor('A'))
        breakpoint_graph1.add_edge(1, 2, Multicolor('C'))
        breakpoint_graph1.add_edge(0, 3, Multicolor('D'))
        assert (len(find_bag_patterns(breakpoint_graph1)) == 1)
        multigraph2 = MultiGraph()
        multigraph2.add_nodes_from(range(4))
        breakpoint_graph2 = BreakpointGraph(multigraph2)
        double_color = ['A', 'B']
        breakpoint_graph2.add_edge(0, 1, Multicolor(*double_color))
        breakpoint_graph2.add_edge(2, 3, Multicolor('B'))
        breakpoint_graph2.add_edge(1, 2, Multicolor('C'))
        breakpoint_graph2.add_edge(0, 3, Multicolor('C'))
        assert (len(find_bag_patterns(breakpoint_graph2)) == 0)

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
        assert (get_size_of_alternating_structures(breakpoint_graph, topology) == 3)

    def test_cycles():
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
        assert (get_size_of_alternating_structures(breakpoint_graph,
                                                   topology,
                                                   get_size_of_paths_instead_of_cycles=False) == 4)

    def test_diamond():
        multigraph = MultiGraph()
        multigraph.add_nodes_from(range(4))
        breakpoint_graph = BreakpointGraph(multigraph)
        first_color = ['A', 'C']
        second_color = ['B', 'D']
        topology = (first_color, second_color)
        breakpoint_graph.add_edge(0, 1, Multicolor(*[A]))
        breakpoint_graph.add_edge(2, 3, Multicolor(*[B]))
        breakpoint_graph.add_edge(1, 2, Multicolor(*[C]))
        breakpoint_graph.add_edge(0, 3, Multicolor(*[D]))
        assert (len(find_diamond_patterns(breakpoint_graph)) == 1)

    test_cylinder()
    test_bag()
    test_paths()
    test_cycles()
    test_diamond()