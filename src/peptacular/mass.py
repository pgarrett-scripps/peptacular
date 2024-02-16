"""
mass.py is a simple module for computing the m/z and mass of an amino acid sequence.
"""
from typing import Union, List, Dict

from peptacular.constants import AVERAGE_AA_MASSES, MONOISOTOPIC_AA_MASSES, MONOISOTOPIC_ION_ADJUSTMENTS, \
    NEUTRON_MASS, PROTON_MASS, AVERAGE_ION_ADJUSTMENTS
from peptacular.sequence import get_modifications, strip_modifications, is_sequence_valid
from peptacular.util import validate_ion_type


def calculate_mass(sequence: str, charge: int = 0, ion_type: str = 'y', monoisotopic: bool = True,
                   isotope: int = 0, loss: float = 0.0, aa_masses: Dict = None) -> float:
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
    :param isotope: The isotope number, defaults to [0].
    :type isotope: int
    :param loss: The loss, defaults to [0.0].
    :type loss: float
    :param aa_masses: A dictionary of amino acid masses, defaults to [None].
    :type aa_masses: Dict

    :raise ValueError: If the ion type is not one of 'a', 'b', 'c', 'x', 'y', or 'z'.

    :return: Mass of the peptide sequence.
    :rtype: float

    .. code-block:: python

        # Calculate the mass of a peptide sequence.
        >>> calculate_mass('PEPTIDE')
        799.359964305

        # Calculate the b-ion mass of a peptide sequence.
        >>> calculate_mass('PEPTIDE', ion_type='b')
        781.3493996049999

        # Calulate the average mass of a peptide sequence.
        >>> calculate_mass('PEPTIDE', monoisotopic=False)
        799.8225199999999

        # Calculate the mass of a peptide sequence with a charge of 2.
        >>> calculate_mass('PEPTIDE', charge=2)
        801.37451723876

        # Calcualte the mass of a modified peptide sequence.
        >>> calculate_mass('PE(3.14)PTIDE[80]', charge=2)
        884.5145172387599

    """

    validate_ion_type(ion_type=ion_type)

    if aa_masses is None:
        aa_masses = MONOISOTOPIC_AA_MASSES if monoisotopic else AVERAGE_AA_MASSES
    else:
        aa_masses = {**MONOISOTOPIC_AA_MASSES, **aa_masses} if monoisotopic else {**AVERAGE_AA_MASSES, **aa_masses}

    # Parse modifications and strip them from sequence
    mods = get_modifications(sequence=sequence)
    stripped_sequence = strip_modifications(sequence=sequence)

    # Calculate mass
    mass = 0
    mass += sum(aa_masses[aa] for aa in stripped_sequence)  # Add amino acids
    mass += sum(float(value) for value in mods.values())  # Add modifications
    mass += (charge * PROTON_MASS)  # Add charge
    mass += isotope * NEUTRON_MASS + loss  # Add isotope and loss
    mass += MONOISOTOPIC_ION_ADJUSTMENTS[ion_type] if monoisotopic else AVERAGE_ION_ADJUSTMENTS[ion_type]
    return mass


def calculate_mz(sequence: str, charge: int = 0, ion_type: str = 'y',
                 monoisotopic: bool = True, isotope: int = 0,
                 loss: Union[List[float], float] = 0.0, aa_masses: Dict = None) -> float:
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
    :param isotope: The isotope number, defaults to [0].
    :type isotope: int
    :param loss: The loss, defaults to [0.0].
    :type loss: float
    :param aa_masses: A dictionary of amino acid masses, defaults to [None].
    :type aa_masses: Dict

    :raise ValueError: If the ion type is not one of 'a', 'b', 'c', 'x', 'y', or 'z'.

    :return: m/z of the peptide sequence.
    :rtype: float

    .. code-block:: python

        # Calculate the m/z of a peptide sequence.
        >>> calculate_mz('PEPTIDE', charge = 1)
        800.3672407718799

        # Calculate the b-ion m/z of a peptide sequence.
        >>> calculate_mz('PEPTIDE', charge = 1, ion_type='b')
        782.3566760718799

        # Calulate the average m/z of a peptide sequence.
        >>> calculate_mz('PEPTIDE', charge = 1, monoisotopic=False)
        800.8297964668799

        # Calculate the m/z of a peptide sequence with a charge of 2.
        >>> calculate_mz('PEPTIDE', charge=2)
        400.68725861938

        # Calcualte the m/z of a modified peptide sequence.
        >>> calculate_mz('PE(3.14)PTIDE[80]', charge=2)
        442.25725861937997

    """

    mass = calculate_mass(sequence=sequence, charge=charge, ion_type=ion_type,
                          monoisotopic=monoisotopic, isotope=isotope, loss=loss,
                          aa_masses=aa_masses)
    return mass if charge == 0 else mass / charge


def valid_mass_sequence(sequence: str) -> bool:
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

    if is_sequence_valid(sequence) is False:
        return False

    mods = get_modifications(sequence)

    for _, mod in mods.items():
        if isinstance(mod, float) is False and isinstance(mod, int) is False:
            return False

    return True
