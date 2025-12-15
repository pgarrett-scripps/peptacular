
from typing import IO, Any
import warnings
import peptacular as pt


def calculate_mass(composition: dict[str, int], monoisotopic: bool = True) -> float:
    """Calculate mass from elemental composition using peptacular's element lookup"""
    from peptacular.elements import ELEMENT_LOOKUP
    
    mass = 0.0
    for element_symbol, count in composition.items():
        # Use (symbol, None) to get the most abundant/monoisotopic isotope
        element_info = ELEMENT_LOOKUP[(element_symbol, None)]
        if monoisotopic:
            mass += element_info.mass * count
        else:
            mass += element_info.average_mass * count
    
    return mass


def format_composition_string(composition: dict[str, int]) -> str:
    """Format composition as a string like C2H3NO"""
    if not composition:
        return ""
    
    # Sort by element symbol for consistency
    parts: list[str] = []
    for element in sorted(composition.keys()):
        count = composition[element]
        if count == 1:
            parts.append(element)
        else:
            parts.append(f"{element}{count}")
    
    return "".join(parts)

def read_obo(file: IO[str]) -> list[dict[str, Any]]:
    file.seek(0)

    info: dict[str, str] = {}
    elems: list[dict[str, Any]] = []
    skip: bool = False
    d: dict[str, Any] | None = None

    for line in file:
        line = line.rstrip()
        if line == "":
            continue

        if line.startswith("[Typedef]"):
            skip = True
            continue

        if line.startswith("[Term]"):
            skip = False
            if d is not None:
                elems.append(d)
            d = {}
            continue

        if d is None:
            key, value = line.split(": ", 1)
            info[key] = value
            continue

        if skip:
            continue

        if line:
            key, value = line.split(": ", 1)

            if key not in d:
                d[key] = [value]
            else:
                d[key].append(value)

    if d is not None:
        elems.append(d)

    return elems

def get_id_and_name(term: dict[str, Any]) -> tuple[str, str]:
    term_id = term.get("id", [])
    term_name = term.get("name", [])

    if len(term_id) > 1:
        warnings.warn(f"Multiple ids for {term_name} {term_id}")
        term_id = term_id[0]
    elif len(term_id) == 1:
        term_id = term_id[0]
    else:
        raise ValueError("Entry name is None")

    if len(term_name) > 1:
        warnings.warn(f"Multiple names for {term_id} {term_name}")
        term_name = term_name[0]
    elif len(term_name) == 1:
        term_name = term_name[0]
    else:
        raise ValueError("Entry id is None")

    return term_id, term_name


def is_obsolete(term: dict[str, Any]) -> bool:
    is_obsolete = term.get("is_obsolete", ["false"])[0]
    if is_obsolete in ["true", "false"]:
        is_obsolete = bool(is_obsolete == "true")
    else:
        raise ValueError(f"Invalid is_obsolete value: {is_obsolete}")

    return is_obsolete



def calculate_composition_mass(composition: dict[pt.ElementInfo, int], monoisotopic: bool = True) -> float:
    total_mass = 0.0
    for elem, count in composition.items():
        total_mass += elem.mass * count
    return total_mass