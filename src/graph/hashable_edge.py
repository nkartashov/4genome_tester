__author__ = 'nikita_kartashov'


class HashableEdge(object):
    def __init__(self, edge):
        self.__edge = edge

    def __eq__(self, other):
        return self.__edge == other.__edge

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        v1, v2 = self.__edge.vertex1, self.__edge.vertex2
        return v1.__hash__() ^ v2.__hash__()