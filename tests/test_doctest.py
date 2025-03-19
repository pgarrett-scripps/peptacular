import pytest
import doctest
import peptacular as pt

modules = [
    pt.chem_calc,
    pt.chem_util,
    pt.combinatoric,
    pt.sequence_funcs,
    pt.mod_builder,
    pt.proforma_parser,
    pt.proforma_dataclasses,
    pt.input_convert,
    pt.randomizer,
    pt.digestion,
    pt.fragmentation,
    pt.glycan,
    pt.isotope,
    pt.mass_calc,
    pt.score,
    pt.util,
]


@pytest.mark.parametrize("module", modules)
def test_doctests(module):
    result = doctest.testmod(module)
    assert result.failed == 0, f"Doctests failed in {module.__name__}"
