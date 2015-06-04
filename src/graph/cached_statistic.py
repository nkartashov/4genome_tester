__author__ = 'nikita_kartashov'


class CachedStatistic(object):
    def __init__(self, f):
        self._f = f
        self._last_bp_graph = None
        self._last_topology = None
        self._cached = None

    def __call__(self, *args, **kwargs):
        breakpoint_graph, topology = args
        if breakpoint_graph != self._last_bp_graph or \
                        topology != self._last_topology:
            self._last_bp_graph = breakpoint_graph
            self._last_topology = topology
            self._cached = self._f(*args, **kwargs)
        return self._cached