from typing import NamedTuple

class Span(NamedTuple):
    start: int
    end: int
    missed_cleavages: int