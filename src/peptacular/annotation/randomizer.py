"""
Randomizer for ProForma annotations.
"""

from random import choice, randint

from ..elements.lookup import ELEMENT_LOOKUP

from ..mods import PSIMOD_LOOKUP, UNIMOD_LOOKUP
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .annotation import ProFormaAnnotation, Interval


def get_random_psimod():
    while True:
        m = PSIMOD_LOOKUP.choice()
        if m.monoisotopic_mass is not None and m.composition is not None:
            return m

def get_random_unimod():
    while True:
        m = UNIMOD_LOOKUP.choice()
        if m.monoisotopic_mass is not None and m.composition is not None:
            return m

def get_random_mod_component() -> str:
    mod_list = choice([get_random_psimod, get_random_unimod])
    m = mod_list()

    mod_type = choice(['name', 'accession'])

    match mod_type:
        case 'name':
            return m.as_tag_name(include_cv=choice([True, False])).serialize()
        case 'accession':
            return m.as_tag_accession().serialize()
        case 'mass':
            return m.as_tag_mass(monoisotopic=True, 
                                 include_cv=True).serialize()
        case _:
            raise ValueError("Invalid mod type")
        
    raise ValueError("Invalid mod type")
    
def get_random_mod_dict(mod_probability: float) -> dict[str, int]:

    mod_dict: dict[str, int] = {}
    while True:
        if randint(0, 100) / 100 < mod_probability:
            mod_component = get_random_mod_component()
            # should be 1 most of the time, but occasionally have multiple mods
            mod_count = 1 if randint(0, 100) < 97 else randint(2, 3)
            mod_dict[mod_component] = mod_count
        else:
            break
        
    return mod_dict

        
elem_infos = ELEMENT_LOOKUP.values()
specific_isotopes = [elem for elem in elem_infos if elem.is_monoisotopic == False]

def generate_random_isotop_mod() -> str:
    isotope = choice(specific_isotopes)
    isotope_str = f"{isotope}"
    return isotope_str


def generate_isotope_mod_dict(mod_probability: float) -> dict[str, int]:

    while True:
        mod_dict: dict[str, int] = {}
        if randint(0, 100) / 100 < mod_probability:
            mod_component = generate_random_isotop_mod()
            mod_dict[mod_component] = 1
            return mod_dict
        else:
            return {}

def generate_static_mod() -> str:
    mod_component = get_random_mod_component()

    # should be 1 most of the time, but occasionally have multiple mods
    number_of_aa = 1 if randint(0, 10) < 9 else randint(2, 3)

    # aa string should be a sequence of random amino acids
    aa_chars = [choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(number_of_aa)]
    aa_string = ','.join(aa_chars)

    return f"[{mod_component}]@{aa_string}"

def generate_static_mods_dict(mod_probability: float) -> dict[str, int]:

    while True:
        mod_dict: dict[str, int] = {}
        if randint(0, 100) / 100 < mod_probability:
            mod_component = generate_static_mod()
            mod_dict[mod_component] = 1
            return mod_dict
        else:
            return {}


def generate_random_intervals(seq_length: int, interval_probability: float = 0.3, mod_probability: float = 0.5) -> list["Interval"]:
    """Generate random non-overlapping intervals for a sequence.
    
    Args:
        seq_length: Length of the sequence
        interval_probability: Probability of creating an interval at each position
        mod_probability: Probability of adding modifications to each interval
        
    Returns:
        List of non-overlapping Interval objects
    """
    from .annotation import Interval
    
    intervals: list[Interval] = []
    i = 0
    
    while i < seq_length:
        # Decide if we should start an interval at this position
        if randint(0, 100) / 100 < interval_probability:
            start = i
            # Determine the end of this interval (must be > start and < seq_length)
            max_end = seq_length - 1
            # Need at least 2 positions for an interval (start < end)
            if start >= max_end:
                # Not enough space for a valid interval, skip
                i += 1
                continue
            
            # Random end position between start+1 and max_end (end must be > start)
            end = randint(start + 1, max_end)
            
            ambiguous = choice([True, False])
            mods = get_random_mod_dict(mod_probability)
            
            interval = Interval(start=start, end=end, ambiguous=ambiguous, mods=mods)
            intervals.append(interval)
            
            # Move past this interval to avoid overlap
            i = end + 1
        else:
            i += 1
    
    return intervals


def random_charge_state(min_charge: int = 0, max_charge: int = 5) -> int:
    return randint(min_charge, max_charge)


charge_adducts = ['H:z+1', 'Na:z+1', 'K:z+1', 'NH4:z+1']
def random_charge_adduct() -> dict[str, int]:
    num_adducts = randint(1, 3)
    d: dict[str, int] = {}
    for _ in range(num_adducts):
        adduct_mult = randint(1, 2)
        adduct = choice(charge_adducts)
        d[adduct] = adduct_mult
    return d
        

def generate_random_proforma_annotation(min_length: int = 6, 
                                        max_length: int = 20, 
                                        mod_probability: float = 0.05,
                                        include_internal_mods: bool = True,
                                        include_nterm_mods: bool = True,
                                        include_cterm_mods: bool = True,
                                        include_labile_mods: bool = True,
                                        include_unknown_mods: bool = True,
                                        include_isotopic_mods: bool = True,
                                        include_static_mods: bool = True,
                                        generate_intervals: bool = True,
                                        ) -> "ProFormaAnnotation":
    
    from .annotation import ProFormaAnnotation

    length = randint(min_length, max_length)
    sequence_chars = [choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(length)]

    internal_mods: dict[int, dict[str, int]] = {}
    for i in range(length):
        mod_dict = get_random_mod_dict(mod_probability)
        if mod_dict:
            internal_mods[i] = mod_dict


    return ProFormaAnnotation(
        sequence=''.join(sequence_chars),
        nterm_mods=get_random_mod_dict(mod_probability) if include_nterm_mods else None,
        cterm_mods=get_random_mod_dict(mod_probability) if include_cterm_mods else None,
        labile_mods=get_random_mod_dict(mod_probability) if include_labile_mods else None,
        unknown_mods=get_random_mod_dict(mod_probability) if include_unknown_mods else None,
        internal_mods=internal_mods if include_internal_mods else None,
        isotope_mods=generate_isotope_mod_dict(mod_probability) if include_isotopic_mods else None,
        static_mods=generate_static_mods_dict(mod_probability) if include_static_mods else None,
        intervals=generate_random_intervals(length, interval_probability=mod_probability) if generate_intervals else None,
        charge=random_charge_state(0, 5) if randint(0, 1) == 1 else random_charge_adduct(),
    )

