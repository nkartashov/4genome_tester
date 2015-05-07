__author__ = 'nikita_kartashov'

from operator import itemgetter
from itertools import chain


class Metrics(object):
    def __init__(self, single_metrics, batch_metrics):
        """
        Constructs Metrics object, which handles all the metrics
        :param single_metrics: annotated tuple of metrics which
        cannot reuse info on different topologies
        :param batch_metrics: annotated tuple of metrics which
        CAN reuse info on different topologies
        :return: the resulting object
        """
        self._single_metrics = tuple(map(itemgetter(0), single_metrics))
        self._batch_metrics = tuple(map(itemgetter(0), batch_metrics))
        self._metric_number = len(single_metrics) + len(batch_metrics)
        self._metric_annotations = tuple(
            chain(*(map(itemgetter(1), metrics) for metrics in (single_metrics, batch_metrics))))

    def metric_number(self):
        return self._metric_number

    def metric_annotations(self):
        return self._metric_annotations

    def run_metrics(self, breakpoint_graph, topologies):
        return chain(*((runner(breakpoint_graph, topologies)
                        for runner in (self._run_single_metrics, self._run_batch_metrics))))

    def _run_single_metrics(self, breakpoint_graph, topologies):
        return (((metric(breakpoint_graph, topology), topology)
                 for topology in topologies)
                for metric in self._single_metrics)

    def _run_batch_metrics(self, breakpoint_graph, topologies):
        return (metric(breakpoint_graph, topologies) for metric in self._batch_metrics)