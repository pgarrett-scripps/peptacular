"""
Randomizer for ProForma annotations.
"""

from functools import lru_cache
from random import choice, randint, random
from typing import TYPE_CHECKING

from ..elements.lookup import ELEMENT_LOOKUP
from ..mods import PSIMOD_LOOKUP, UNIMOD_LOOKUP

if TYPE_CHECKING:
    from ..mods import PsimodInfo, UnimodInfo
    from .annotation import Interval, ProFormaAnnotation

# Constants
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
DEFAULT_MOD_PROBABILITY = 0.05
DEFAULT_INTERVAL_PROBABILITY = 0.3
MULTIPLE_MOD_THRESHOLD = 0.97  # 97% chance of single mod, 3% chance of multiple
MULTIPLE_AA_THRESHOLD = 0.9  # 90% chance of single AA in static mod


def get_random_psimod() -> "PsimodInfo":
    """Get a random PSI-MOD modification with mass and composition."""
    return PSIMOD_LOOKUP.choice(
        require_monoisotopic_mass=True, require_composition=True
    )


def get_random_unimod() -> "UnimodInfo":
    """Get a random Unimod modification with mass and composition."""
    return UNIMOD_LOOKUP.choice(
        require_monoisotopic_mass=True, require_composition=True
    )


def get_random_mod_component(require_composition: bool = True) -> str:
    """Generate a random modification component string.

    Args:
        require_composition: If True, only name and accession types are allowed (no mass-only).
                            If False, mass-based modifications can also be generated.
    """
    mod_list = choice([get_random_psimod, get_random_unimod])
    m = mod_list()

    # When composition is required, exclude mass-only modifications
    mod_type = choice(
        ["name", "accession"] if require_composition else ["name", "accession", "mass"]
    )

    match mod_type:
        case "name":
            return m.as_tag_name(include_cv=choice([True, False])).serialize()
        case "accession":
            return m.as_tag_accession().serialize()
        case "mass":
            return m.as_tag_mass(
                monoisotopic=True, include_cv=choice([True, False])
            ).serialize()
        case _:
            raise ValueError(f"Invalid mod type: {mod_type}")


def get_random_mod_dict(
    mod_probability: float, require_composition: bool = True
) -> dict[str, int]:
    """Generate a dictionary of random modifications.

    Args:
        mod_probability: Probability of adding each modification.
        require_composition: If True, only modifications with composition are allowed (no mass-only).
    """
    mod_dict: dict[str, int] = {}
    while random() < mod_probability:
        mod_component = get_random_mod_component(
            require_composition=require_composition
        )
        # should be 1 most of the time, but occasionally have multiple mods
        mod_count = 1 if random() < MULTIPLE_MOD_THRESHOLD else randint(2, 3)
        mod_dict[mod_component] = mod_dict.get(mod_component, 0) + mod_count

    return mod_dict


@lru_cache(maxsize=1)
def _get_specific_isotopes() -> tuple:
    """Cache the list of non-monoisotopic isotopes."""
    return tuple(
        elem for elem in ELEMENT_LOOKUP.values() if elem.is_monoisotopic is not None
    )


def generate_random_isotope_mod() -> str:
    """Generate a random isotope modification string."""
    isotope = choice(_get_specific_isotopes())
    return str(isotope)


def generate_isotope_mod_dict(mod_probability: float) -> dict[str, int]:
    """Generate a dictionary with a single random isotope modification."""
    if random() < mod_probability:
        mod_component = generate_random_isotope_mod()
        return {mod_component: 1}
    return {}


def generate_static_mod(require_composition: bool = True) -> str:
    """Generate a static modification string with target amino acids.

    Args:
        require_composition: If True, only modifications with composition are allowed (no mass-only).
    """
    mod_component = get_random_mod_component(require_composition=require_composition)

    # should be 1 most of the time, but occasionally have multiple mods
    number_of_aa = 1 if random() < MULTIPLE_AA_THRESHOLD else randint(2, 3)

    # aa string should be a sequence of random amino acids
    aa_chars = [choice(AMINO_ACIDS) for _ in range(number_of_aa)]
    aa_string = ",".join(aa_chars)

    return f"[{mod_component}]@{aa_string}"


def generate_static_mods_dict(
    mod_probability: float, require_composition: bool = True
) -> dict[str, int]:
    """Generate a dictionary with a single static modification.

    Args:
        mod_probability: Probability of generating a static modification.
        require_composition: If True, only modifications with composition are allowed (no mass-only).
    """
    if random() < mod_probability:
        mod_component = generate_static_mod(require_composition=require_composition)
        return {mod_component: 1}
    return {}


def generate_random_intervals(
    seq_length: int,
    interval_probability: float = DEFAULT_INTERVAL_PROBABILITY,
    mod_probability: float = DEFAULT_MOD_PROBABILITY,
    require_composition: bool = True,
) -> list["Interval"]:
    """Generate random non-overlapping intervals for a sequence.

    Args:
        seq_length: Length of the sequence
        interval_probability: Probability of creating an interval at each position
        mod_probability: Probability of adding modifications to each interval
        require_composition: If True, only modifications with composition are allowed (no mass-only)

    Returns:
        List of non-overlapping Interval objects
    """
    from .annotation import Interval

    intervals: list[Interval] = []
    i = 0

    while i < seq_length - 1:  # Need at least 2 positions for an interval
        # Decide if we should start an interval at this position
        if random() < interval_probability:
            start = i
            max_end = seq_length - 1

            # Random end position between start+1 and max_end (end must be > start)
            end = randint(start + 1, max_end)

            ambiguous = choice([True, False])
            mods = get_random_mod_dict(
                mod_probability, require_composition=require_composition
            )

            interval = Interval(start=start, end=end, ambiguous=ambiguous, mods=mods)
            intervals.append(interval)

            # Move past this interval to avoid overlap
            i = end + 1
        else:
            i += 1

    return intervals


def random_charge_state(min_charge: int = 0, max_charge: int = 5) -> int:
    """Generate a random charge state."""
    return randint(min_charge, max_charge)


CHARGE_ADDUCTS = ["H:z+1", "Na:z+1", "K:z+1", "NH4:z+1"]


def random_charge_adduct() -> dict[str, int]:
    """Generate a random charge adduct dictionary."""
    num_adducts = randint(1, 3)
    d: dict[str, int] = {}
    for _ in range(num_adducts):
        adduct_mult = randint(1, 2)
        adduct = choice(CHARGE_ADDUCTS)
        d[adduct] = d.get(adduct, 0) + adduct_mult
    return d


def generate_random_proforma_annotation(
    min_length: int = 6,
    max_length: int = 20,
    mod_probability: float = DEFAULT_MOD_PROBABILITY,
    include_internal_mods: bool = True,
    include_nterm_mods: bool = True,
    include_cterm_mods: bool = True,
    include_labile_mods: bool = True,
    include_unknown_mods: bool = True,
    include_isotopic_mods: bool = True,
    include_static_mods: bool = True,
    generate_intervals: bool = True,
    include_charge: bool = True,
    require_composition: bool = True,
) -> "ProFormaAnnotation":
    from .annotation import ProFormaAnnotation

    length = randint(min_length, max_length)
    sequence_chars = [choice(AMINO_ACIDS) for _ in range(length)]

    internal_mods: dict[int, dict[str, int]] = {}
    if include_internal_mods:
        for i in range(length):
            mod_dict = get_random_mod_dict(
                mod_probability, require_composition=require_composition
            )
            if mod_dict:
                internal_mods[i] = mod_dict

    # Determine charge
    charge: int | dict[str, int] | None = None
    if include_charge:
        if random() < 0.5:
            charge = random_charge_state(0, 5)
        else:
            charge = random_charge_adduct()

    return ProFormaAnnotation(
        sequence="".join(sequence_chars),
        nterm_mods=get_random_mod_dict(
            mod_probability, require_composition=require_composition
        )
        if include_nterm_mods
        else None,
        cterm_mods=get_random_mod_dict(
            mod_probability, require_composition=require_composition
        )
        if include_cterm_mods
        else None,
        labile_mods=get_random_mod_dict(
            mod_probability, require_composition=require_composition
        )
        if include_labile_mods
        else None,
        unknown_mods=get_random_mod_dict(
            mod_probability, require_composition=require_composition
        )
        if include_unknown_mods
        else None,
        internal_mods=internal_mods if internal_mods else None,
        isotope_mods=generate_isotope_mod_dict(mod_probability)
        if include_isotopic_mods
        else None,
        static_mods=generate_static_mods_dict(
            mod_probability, require_composition=require_composition
        )
        if include_static_mods
        else None,
        intervals=generate_random_intervals(
            length,
            interval_probability=mod_probability,
            require_composition=require_composition,
        )
        if generate_intervals
        else None,
        charge=charge,
    )
