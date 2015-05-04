__author__ = 'nikita_kartashov'


def normalize_multicolor(multicolor, all_genomes):
    first_color = frozenset(multicolor)
    second_color = all_genomes - first_color
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


def edge_vertices(edge):
    return [edge.vertex1, edge.vertex2]


def vertex_multidegree(breakpoint_graph, vertex):
    return len(list(breakpoint_graph.get_edges_by_vertex(vertex)))


def is_vertex_simple(breakpoint_graph, vertex):
    return vertex_multidegree(breakpoint_graph, vertex) == 2


def is_edge_simple(breakpoint_graph, edge):
    return all(is_vertex_simple(breakpoint_graph, vertex) for vertex in edge_vertices(edge))


def filter_n_sized_colors(n, edges):
    return (edge for edge in edges if len(edge.multicolor.colors) == n)


def get_vertex_neighbours(breakpoint_graph, vertex, edges=None, with_edges=True):
    def neighbour_getter(edge):
        neighbour = next(edge_vertex for edge_vertex in edge_vertices(edge) if edge_vertex != vertex)
        return (neighbour, edge) if with_edges else neighbour

    if edges is None:
        edges = breakpoint_graph.get_edges_by_vertex(vertex)

    return map(neighbour_getter, edges)


def get_vertex_predicate_neighbours(breakpoint_graph, vertex, predicate=lambda _: True, edges=None, with_edges=True):
    if edges is None:
        edges = breakpoint_graph.get_edges_by_vertex(vertex)

    return get_vertex_neighbours(breakpoint_graph, vertex, filter(predicate, edges), with_edges)


def get_vertex_sized_neighbours(breakpoint_graph, vertex, size, edges=None, with_edges=True):
    return get_vertex_predicate_neighbours(breakpoint_graph,
                                           vertex,
                                           lambda edge: len(edge.multicolor.colors) == size,
                                           edges,
                                           with_edges)


def get_vertex_colored_neighbours(breakpoint_graph, vertex, color, edges=None, with_edges=True):
    return get_vertex_predicate_neighbours(breakpoint_graph,
                                           vertex,
                                           lambda edge: edge.multicolor == color,
                                           edges,
                                           with_edges)
