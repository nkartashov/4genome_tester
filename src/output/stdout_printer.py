__author__ = 'nikita_kartashov'

from itertools import chain

from .printer import Printer



class StdOutPrinter(Printer):
    def __init__(self):
        super().__init__()

    def write_header(self, header, width):
        print('\t'.join(word.ljust(width) for word in header))

    def write_row(self, prefix, result, width):
        print('\t'.join(word.ljust(width) for word in chain(prefix, map(str, result))))

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass