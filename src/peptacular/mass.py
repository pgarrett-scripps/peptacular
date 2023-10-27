"""
mass.py is a simple module for computing the m/z and mass of an amino acid sequence, plus it works with peptacular's
modification notation!
"""

from peptacular.constants import MONO_ISOTOPIC_ATOMIC_MASSES, AVERAGE_ATOMIC_MASSES, AVERAGE_AA_MASSES, \
    MONO_ISOTOPIC_AA_MASSES, ION_ADJUSTMENTS, UWPR_MONO_ISOTOPIC_ATOMIC_MASSES, UWPR_AVERAGE_AA_MASSES, \
    UWPR_AVERAGE_ATOMIC_MASSES, UWPR_MONO_ISOTOPIC_AA_MASSES
from peptacular.sequence import get_modifications, strip_modifications
from peptacular.util import validate_ion_type


def calculate_mass(sequence: str, charge: int = 0, ion_type: str = 'y', monoisotopic: bool = True,
                   uwpr_mass: bool = False) -> float:
    """
    Calculate the mass of an amino acid 'sequence'.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param charge: The charge state, default is [0].
    :type charge: int
    :param ion_type: The ion type, default is [y].
    :type ion_type: str
    :param monoisotopic: If True, use monoisotopic mass else use average mass, defaults to [True].
    :type monoisotopic: bool
    :param uwpr_mass: If true, uses uwpr masses. If false, uses Pyteomics masses, defaults to [False].
    :type uwpr_mass: bool

    :raise ValueError: If the ion type is not one of 'a', 'b', 'c', 'x', 'y', or 'z'.

    :return: Mass of the peptide sequence.
    :rtype: float

    .. code-block:: python

        # Calculate the mass of a peptide sequence.
        >>> calculate_mass('PEPTIDE')
        799.3599640267099

        # Calculate the b-ion mass of a peptide sequence.
        >>> calculate_mass('PEPTIDE', ion_type='b')
        781.34939934301

        # Calulate the average mass of a peptide sequence.
        >>> calculate_mass('PEPTIDE', monoisotopic=False)
        799.8280646837

        # Calculate the mass of a peptide sequence with a charge of 2.
        >>> calculate_mass('PEPTIDE', charge=2)
        801.3745169602499

        # Calcualte the mass of a modified peptide sequence.
        >>> calculate_mass('PE(3.14)PTIDE[80]')
        882.4999640267099

    """

    validate_ion_type(ion_type=ion_type)

    # Select mass set based on monoisotopic flag
    atomic_masses, aa_masses = (MONO_ISOTOPIC_ATOMIC_MASSES, MONO_ISOTOPIC_AA_MASSES) \
        if monoisotopic is True else (AVERAGE_ATOMIC_MASSES, AVERAGE_AA_MASSES)

    if uwpr_mass is True:
        atomic_masses, aa_masses = (UWPR_MONO_ISOTOPIC_ATOMIC_MASSES, UWPR_MONO_ISOTOPIC_AA_MASSES) \
            if monoisotopic is True else (UWPR_AVERAGE_ATOMIC_MASSES, UWPR_AVERAGE_AA_MASSES)

    # Parse modifications and strip them from sequence
    mods = get_modifications(sequence=sequence)
    stripped_sequence = strip_modifications(sequence=sequence)

    # Calculate mass
    mass = sum(aa_masses[aa] for aa in stripped_sequence)
    mass += sum(float(value) for value in mods.values())
    mass += (charge * atomic_masses['PROTON'])

    return mass + ION_ADJUSTMENTS[ion_type]


def calculate_mz(sequence: str, charge: int = 0, ion_type: str = 'y',
                 monoisotopic: bool = True, uwpr_mass: bool = False) -> float:
    """
    Calculate the m/z (mass-to-charge ratio) of an amino acid 'sequence'.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param charge: The charge state, default is [0].
    :type charge: int
    :param ion_type: The ion type, default is [y].
    :type ion_type: str
    :param monoisotopic: If True, use monoisotopic mass else use average mass, defaults to [True].
    :type monoisotopic: bool
    :param uwpr_mass: If true, uses uwpr masses. If false, uses Pyteomics masses, defaults to [False].
    :type uwpr_mass: bool

    :raise ValueError: If the ion type is not one of 'a', 'b', 'c', 'x', 'y', or 'z'.

    :return: m/z of the peptide sequence.
    :rtype: float

    .. code-block:: python

        # Calculate the m/z of a peptide sequence.
        >>> calculate_mz('PEPTIDE', charge = 1)
        800.3672404934799

        # Calculate the b-ion m/z of a peptide sequence.
        >>> calculate_mz('PEPTIDE', charge = 1, ion_type='b')
        782.35667580978

        # Calulate the average m/z of a peptide sequence.
        >>> calculate_mz('PEPTIDE', charge = 1, monoisotopic=False)
        800.83534115047

        # Calculate the m/z of a peptide sequence with a charge of 2.
        >>> calculate_mz('PEPTIDE', charge=2)
        400.68725848012497

        # Calcualte the m/z of a modified peptide sequence.
        >>> calculate_mz('PE(3.14)PTIDE[80]', charge=2)
        442.25725848012496

    """

    mass = calculate_mass(sequence=sequence, charge=charge, ion_type=ion_type,
                          monoisotopic=monoisotopic, uwpr_mass=uwpr_mass)
    return mass if charge == 0 else mass / charge


def valid_mass_sequence(sequence: str):
    """
    Check if a sequence is a valid mass-based sequence, where all modifications are either of type int or float.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: True if the sequence is a valid mass sequence, False otherwise.
    :rtype: bool

    .. code-block:: python

        # Check if a sequence is a valid mass sequence.
        >>> valid_mass_sequence('PEPTIDE')
        True
        >>> valid_mass_sequence('PEPTIDE[80]')
        True
        >>> valid_mass_sequence('PEPTIDE[80.0]')
        True
        >>> valid_mass_sequence('PEPTIDE(Ox)PEPTIDE')
        False

    """

    stripped_sequence = strip_modifications(sequence)
    mods = get_modifications(sequence)

    for aa in stripped_sequence:
        if aa not in MONO_ISOTOPIC_AA_MASSES:
            return False

    for _, mod in mods.items():
        if isinstance(mod, float) is False and isinstance(mod, int) is False:
            return False

    return True
