import pytest
from peptacular import sequence, digest, fragment, mass, protein, score, spans
from peptacular.term import modification, residue

modules = [
    sequence,
    digest,
    fragment,
    mass,
    protein,
    score,
    spans,
    modification,
    residue
]


@pytest.mark.parametrize("module", modules)
def test_doctests(module):
    import doctest
    doctest.testmod(module)
