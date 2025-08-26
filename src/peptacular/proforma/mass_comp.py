from __future__ import annotations
from collections import Counter
from typing import TYPE_CHECKING
import warnings

from .dclasses.modlist import ModList
from ..chem.chem_constants import AVERAGE_AA_MASSES, MONOISOTOPIC_AA_MASSES
from ..chem.chem_util import chem_mass
from ..mass_calc import adjust_mass, adjust_mz, mod_mass
from ..util import parse_static_mods
from ..chem.chem_calc import (
    parse_charge_adducts_comp,
    parse_mod_delta_mass_only,
    apply_isotope_mods_to_composition,
    estimate_comp,
    mod_comp,
)
from ..errors import AmbiguousAminoAcidError, UnknownAminoAcidError
from ..constants import (
    AA_COMPOSITIONS,
    FRAGMENT_ION_BASE_CHARGE_ADDUCTS,
    NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS,
    IonType,
    IonTypeLiteral,
)

if TYPE_CHECKING:
    from .annotation import ProFormaAnnotation


def _pop_delta_mass_mods(annotation: ProFormaAnnotation, 
                         inplace: bool = True) -> float:
    """
    Pop the delta mass modifications. This leaves only modifications which
    can have their composition calculated.

    :param annotation: The ProForma annotation.
    :type annotation: ProFormaAnnotationBase

    :return: The delta mass.
    :rtype: float

    .. code-block:: python

        >>> mods = ProFormaAnnotationBase(sequence='', nterm_mods = [42.0, -20.0])
        >>> mods.pop_delta_mass_mods()
        22.0
        >>> mods
        ProFormaAnnotationBase(sequence=)

    """
    if inplace is False:
        # Create a copy of the annotation to modify
        return _pop_delta_mass_mods(annotation=annotation.copy(), inplace=True)

    if annotation.has_static_mods:
        raise ValueError(
            "Static mods are present. Cannot pop delta mass mods when static mods exist."
        )

    if annotation.has_isotope_mods:
        raise ValueError(
            "Isotope mods are present. Cannot pop delta mass mods when isotope mods exist."
        )

    def _pop_delta_mass_from_mod_list(mod_list: ModList) -> float:
        delta_mass = 0.0

        # Iterate backwards to safely remove items during iteration
        for i in range(len(mod_list) - 1, -1, -1):
            mod = mod_list[i]
            val = parse_mod_delta_mass_only(mod)
            if val is not None:
                delta_mass += val
                mod_list.remove(mod)

        return delta_mass

    delta_mass = 0.0
    delta_mass += _pop_delta_mass_from_mod_list(
        mod_list=annotation.get_labile_mod_list()
    )
    delta_mass += _pop_delta_mass_from_mod_list(
        mod_list=annotation.get_unknown_mod_list()
    )
    delta_mass += _pop_delta_mass_from_mod_list(
        mod_list=annotation.get_nterm_mod_list()
    )
    delta_mass += _pop_delta_mass_from_mod_list(
        mod_list=annotation.get_cterm_mod_list()
    )

    # Intervals
    for interval in annotation.get_interval_list():
        for i in range(len(interval.mods) - 1, -1, -1):
            mod = interval.mods[i]
            val = parse_mod_delta_mass_only(mod=mod)
            if val is not None:
                delta_mass += val
                interval.mods.pop(i)

    # Internal mods
    for _, mods in annotation.get_internal_mod_dict().items():
        delta_mass += _pop_delta_mass_from_mod_list(mod_list=mods)

    return delta_mass


def _get_sequence_composition(
    sequence: str,
) -> Counter[str]:
    # Get the composition of the base sequence
    sequence_composition: Counter[str] = Counter()
    for aa in sequence:
        try:
            aa_comp = AA_COMPOSITIONS[aa]
        except KeyError as err:
            raise UnknownAminoAcidError(aa) from err
        sequence_composition.update(aa_comp)

    return sequence_composition


def _sequence_comp(
    annotation: ProFormaAnnotation,
    ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
    isotope: int = 0,
    use_isotope_on_mods: bool = False,
) -> dict[str, int | float]:

    # If charge is not provided, set it to 0
    charge = 0
    if annotation.has_charge:
        charge = annotation.charge

    if charge is None:
        raise ValueError("Charge state is not set on the annotation.")

    # if charge_adducts is not provided, set it to None
    charge_adducts = None
    if annotation.has_charge_adducts:
        if len(annotation.get_charge_adduct_list()) != 1:
            raise ValueError("Multiple charge adducts found; expected exactly one.")
        charge_adducts = annotation.get_charge_adduct_list()[0]

    if charge_adducts is None:
        if ion_type in (IonType.PRECURSOR, IonType.NEUTRAL):
            charge_adducts = f"{charge}H+"
        else:
            charge_adducts = (
                f"{charge-1}H+,{FRAGMENT_ION_BASE_CHARGE_ADDUCTS[ion_type]}"
            )

    if ion_type not in (IonType.PRECURSOR, IonType.NEUTRAL):
        if charge == 0:
            warnings.warn(
                "Calculating the comp of a fragment ion with charge state 0. Fragment ions should have a "
                "charge state greater than 0 since the neutral form doesnt exist."
            )

    if annotation.contains_mass_ambiguity():
        raise AmbiguousAminoAcidError(
            "B/Z",
            "Cannot determine the composition of a sequence with ambiguous amino acids.",
        )

    # Get the composition of the base sequence
    sequence_composition: Counter[str] = _get_sequence_composition(
        sequence=annotation.sequence
    )
    sequence_composition.update(NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[ion_type])
    sequence_composition.update(parse_charge_adducts_comp(adducts=charge_adducts))

    mod_composition: Counter[str] = Counter()
    for unknown_mod in annotation.get_unknown_mod_list():
        mod_composition.update(mod_comp(mod=unknown_mod))

    for interval in annotation.get_interval_list():
        if interval.has_mods:
            for interval_mod in interval.mods:
                mod_composition.update(mod_comp(mod=interval_mod))

    if annotation.has_labile_mods and ion_type == IonType.PRECURSOR:
        for labile_mod in annotation.get_labile_mod_list():
            mod_composition.update(mod_comp(mod=labile_mod))

    for nterm_mod in annotation.get_nterm_mod_list():
        mod_composition.update(mod_comp(mod=nterm_mod))

    for cterm_mod in annotation.get_cterm_mod_list():
        mod_composition.update(mod_comp(mod=cterm_mod))

    for _, internal_mods in annotation.get_internal_mod_dict().items():
        for internal_mod in internal_mods:
            mod_composition.update(mod_comp(mod=internal_mod))

    if annotation.has_static_mods:
        static_map = parse_static_mods(mods=annotation.static_mods)

        n_term_mod = static_map.get("N-Term")
        if n_term_mod is not None:
            for m in n_term_mod:
                mod_composition.update(mod_comp(mod=m))

        c_term_mod = static_map.get("C-Term")
        if c_term_mod is not None:
            for m in c_term_mod:
                mod_composition.update(mod_comp(mod=m))

        for aa, mod in static_map.items():
            if aa in ["N-Term", "C-Term"]:
                continue

            aa_count = annotation.sequence.count(aa)
            for m in mod:
                comp = mod_comp(mod=m)
                mod_composition.update({k: v * aa_count for k, v in comp.items()})

    mod_composition["n"] = mod_composition.get("n", 0) + isotope

    # Apply isotopic mods
    if annotation.has_isotope_mods:
        if use_isotope_on_mods:
            sequence_composition = apply_isotope_mods_to_composition(
                sequence_composition, annotation.isotope_mods  # type: ignore
            )  # type: ignore
            mod_composition = apply_isotope_mods_to_composition(
                mod_composition, annotation.isotope_mods  # type: ignore
            )  # type: ignore
        else:
            sequence_composition = apply_isotope_mods_to_composition(
                sequence_composition, annotation.isotope_mods  # type: ignore
            )  # type: ignore

    composition: dict[str, int | float] = {}
    for k, v in sequence_composition.items():
        composition[k] = composition.get(k, 0) + v

    for k, v in mod_composition.items():
        composition[k] = composition.get(k, 0) + v

    composition = {k: v for k, v in composition.items() if v != 0}

    return composition


def mass(
    annotation: ProFormaAnnotation,
    ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotope: int = 0,
    loss: float = 0.0,
    use_isotope_on_mods: bool = False,
    precision: int | None = None,
) -> float:

    if annotation.contains_mass_ambiguity():
        raise AmbiguousAminoAcidError(
            aa=",".join(annotation.get_residue_ambiguity_residues()),
            msg="Cannot determine the mass of a sequence with ambiguous amino acids: {annotation.sequence}",
        )

    annotation = annotation.copy()
    annotation.condense_static_mods(inplace=True)

    # more complex mass calculation (based on chem composition)
    if annotation.has_isotope_mods:
        peptide_composition, delta_mass = _comp_with_delta_mass(
            annotation=annotation,
            ion_type=ion_type,
            isotope=isotope,
            use_isotope_on_mods=use_isotope_on_mods,
        )

        return (
            chem_mass(
                formula=peptide_composition, monoisotopic=monoisotopic, precision=precision
            )
            + delta_mass
        )

    if ion_type not in (IonType.PRECURSOR, IonType.NEUTRAL):
        if annotation.charge is None or annotation.charge == 0:
            warnings.warn(
                "Calculating the mass of a fragment ion with charge state 0. Fragment ions should have a "
                "charge state greater than 0 since the neutral mass doesnt exist."
            )

    m: float = 0.0
    """
    if annotation.has_static_mods():
        static_map = parse_static_mods(annotation.static_mods)

        n_term_mod = static_map.get("N-Term")
        if n_term_mod is not None:
            m += sum(mod_mass(m, monoisotopic, precision=None) for m in n_term_mod)

        c_term_mod = static_map.get("C-Term")
        if c_term_mod is not None:
            m += sum(mod_mass(m, monoisotopic, precision=None) for m in c_term_mod)

        for aa, mod in static_map.items():

            if aa in ["N-Term", "C-Term"]:
                continue

            aa_count = annotation.sequence.count(aa)
            m += sum(mod_mass(m, monoisotopic, precision=None) for m in mod) * aa_count
    """
    try:
        m += sum(
            MONOISOTOPIC_AA_MASSES[aa] if monoisotopic else AVERAGE_AA_MASSES[aa]
            for aa in annotation.sequence
        )
    except KeyError as err:
        raise UnknownAminoAcidError(str(err)) from err

    # Apply labile mods
    if ion_type == IonType.PRECURSOR:
        for mod in annotation.get_labile_mod_list():
            m += mod_mass(mod)

    # Apply Unknown mods
    for mod in annotation.get_unknown_mod_list():
        m += mod_mass(mod)

    # Apply N-term mods
    for mod in annotation.get_nterm_mod_list():
        m += mod_mass(mod)

    # Apply intervals
    for interval in annotation.get_interval_list():
        for mod in interval.mods:
            m += mod_mass(mod)

    # apply internal mods
    for _, mods in annotation.get_internal_mod_dict().items():
        for mod in mods:
            m += mod_mass(mod)

    # apply C-term mods
    for mod in annotation.get_cterm_mod_list():
        m += mod_mass(mod)

    return adjust_mass(
        base_mass=m,
        charge=annotation.charge,
        ion_type=ion_type,
        monoisotopic=monoisotopic,
        isotope=isotope,
        loss=loss,
        charge_adducts=annotation.get_charge_adduct_list().data,
        precision=precision,
    )


def mz(
    annotation: ProFormaAnnotation,
    ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotope: int = 0,
    loss: float = 0.0,
    precision: int | None = None,
    use_isotope_on_mods: bool = False,
) -> float:

    m = mass(
        annotation=annotation,
        ion_type=ion_type,
        monoisotopic=monoisotopic,
        isotope=isotope,
        loss=loss,
        use_isotope_on_mods=use_isotope_on_mods,
        precision=precision,
    )

    return adjust_mz(base_mass=m, charge=annotation.charge, precision=precision)


def comp(
    annotation: ProFormaAnnotation,
    ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
    estimate_delta: bool = False,
    isotope: int = 0,
    use_isotope_on_mods: bool = False,
) -> dict[str, int | float]:

    composition, delta_mass = _comp_with_delta_mass(
        annotation,
        ion_type=ion_type,
        isotope=isotope,
        use_isotope_on_mods=use_isotope_on_mods,
    )

    if delta_mass != 0:

        if estimate_delta is False:
            raise ValueError(
                f"Non-zero delta mass ({delta_mass}) encountered without estimation enabled for "
                f"sequence '{annotation.serialize(include_plus=True, precision=5)}'."
            )

        delta_mass_comp = estimate_comp(delta_mass, annotation.isotope_mods)

        # Combine the compositions of the sequence and the delta mass
        for element in delta_mass_comp:
            composition.setdefault(element, 0)
            composition[element] += delta_mass_comp[element]

    return composition


def _comp_with_delta_mass(
    annotation: ProFormaAnnotation,
    ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
    isotope: int = 0,
    use_isotope_on_mods: bool = False,
) -> tuple[dict[str, int | float], float]:
    # returns a tuple containing the composition of the peptide, and the sum of all delta mass modifications, 
    # which cannot be directly converted to a composition dictionary

    annotation = annotation.copy()

    # condenses static mods
    annotation.condense_static_mods(inplace=True)
    delta_mass = _pop_delta_mass_mods(
        annotation
    )  # sum delta mass mods and remove them from the annotation
    peptide_composition = _sequence_comp(
        annotation, ion_type, isotope, use_isotope_on_mods
    )

    if use_isotope_on_mods is True and delta_mass != 0:

        warnings.warn(
            "use_isotope_on_mods=True and delta_mass != 0. Cannot apply isotopic modifications to the delta mass."
        )

    return peptide_composition, delta_mass


def condense_to_mass_mods(
    annotation: ProFormaAnnotation,
    include_plus: bool = False,
    inplace: bool = True,
    use_isotope_on_mods: bool = False,
    estimate_delta: bool = False,
) -> ProFormaAnnotation:

    if inplace is False:
        # Create a copy of the annotation to modify
        return condense_to_mass_mods(
            annotation.copy(),
            include_plus=include_plus,
            inplace=True,
            use_isotope_on_mods=use_isotope_on_mods,
            estimate_delta=estimate_delta,
        )

    labile_mods = annotation.pop_labile_mods()
    isotope_mods = annotation.isotope_mods

    # fix this so that isotopes get applied
    n_term_mods = annotation.pop_nterm_mods()
    c_term_mods = annotation.pop_cterm_mods()
    n_term_mods_mass, c_term_mods_mass = None, None

    if n_term_mods:
        n_term_annot = ProFormaAnnotation(
            sequence="",
            nterm_mods=n_term_mods,
            isotope_mods=isotope_mods,
        )
        n_term_mods_mass = mass(n_term_annot, ion_type="n")

    if c_term_mods:
        c_term_annot = ProFormaAnnotation(
            sequence="",
            cterm_mods=c_term_mods,
            isotope_mods=isotope_mods,
        )
        c_term_mods_mass = mass(c_term_annot, ion_type="n")

    # Split into segments and process each one
    segments = list(annotation.split())
    stripped_segments = [seg.strip(inplace=False) for seg in segments]

    # Calculate mass differences
    for i, (segment, stripped) in enumerate(zip(segments, stripped_segments)):
        # Calculate mass difference
        mass_diff = mass(segment) - mass(stripped)
        if abs(mass_diff) > 1e-3:  # Only add if difference is significant
            annotation.set_internal_mods_at_index(index=i, mods=mass_diff, inplace=True)

    if labile_mods:
        labile_mods_mass = sum(mod_mass(mod) for mod in labile_mods)
        annotation.set_labile_mods(mods=labile_mods_mass, inplace=True)

    if n_term_mods_mass is not None:
        annotation.set_nterm_mods(mods=n_term_mods_mass, inplace=True)

    if c_term_mods_mass is not None:
        annotation.set_cterm_mods(mods=c_term_mods_mass, inplace=True)

    annotation.pop_isotope_mods()

    return annotation
