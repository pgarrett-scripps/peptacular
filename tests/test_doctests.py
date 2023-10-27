import pytest
from peptacular import sequence, digest, fragment, mass, protein, score, spans, term

modules = [
    sequence,
    digest,
    fragment,
    mass,
    protein,
    score,
    spans,
    term,
]


@pytest.mark.parametrize("module", modules)
def test_doctests(module):
    import doctest
    doctest.testmod(module)
