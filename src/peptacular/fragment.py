from dataclasses import dataclass
from typing import List, Generator, Union

from peptacular.constants import PROTON_MASS
from peptacular.mass import calculate_mz
from peptacular.sequence import strip_modifications
from peptacular.term import strip_c_term, strip_n_term
from peptacular.util import is_forward


@dataclass(frozen=True)
class Fragment:
    """
    Represents a fragment resulting from breaking a sequence at specific points.

    :ivar str sequence: The sequence of the fragment.
    :ivar float mz: The mass-to-charge (m/z) ratio of the fragment.
    :ivar int charge: The charge state of the fragment.
    :ivar str ion_type: The type of ion formed by the fragmentation.
    :ivar int number: The fragment ion number, for abc ions this refers to the number of residues starting from the
                      N-terminus, for xyz ions this refers to the number of residues starting from the C-terminus.
    :ivar bool internal: True if the fragment is an internal fragment, False if it is a terminal fragment.
    :ivar int parent_number: The number of the parent fragment from which the internal fragment is derived. For
    terminal fragments (internal=False) the number and parent number should always be equal.

    .. code-block:: python

        >>> fragment = Fragment("PEPTIDE", 123.456, 1, "b", 2, False, 2)
        >>> fragment.sequence
        'PEPTIDE'
        >>> fragment.mz
        123.456
        >>> fragment.charge
        1
        >>> fragment.ion_type
        'b'
        >>> fragment.number
        2
        >>> fragment.internal
        False
        >>> fragment.parent_number
        2

    """

    sequence: str
    mz: float
    charge: int
    ion_type: str
    number: int
    internal: bool
    parent_number: int

    def __post_init__(self):
        """
        Validates the fragment after initialization.
        """
        if self.internal is False and self.number != self.parent_number:
            raise ValueError("Parent number must be equal to number for terminal fragments.")

    @property
    def start(self):
        """
        Returns the start index of the fragment sequence in the parent sequence.

        :return: Start index of the fragment.
        :rtype: int

        :raises ValueError: If the ion type is invalid.
        """
        if self.ion_type in 'abc':
            if self.internal is False:
                return 0
            return self.parent_number - self.number

        if self.ion_type in 'xyz':
            return -self.parent_number

        raise ValueError(f"Invalid ion type: {self.ion_type}")

    @property
    def end(self):
        """
        Returns the end index of the fragment sequence in the parent sequence.

        :return: End index of the fragment.
        :rtype: int

        :raises ValueError: If the ion type is invalid.
        """
        if self.ion_type in 'abc':
            return self.parent_number

        if self.ion_type in 'xyz':
            if self.internal is False:
                return None
            return self.start + self.number

        raise ValueError(f"Invalid ion type: {self.ion_type}")

    @property
    def mass(self):
        """
        Returns the mass of the fragment.

        :return: Mass of the fragment.
        :rtype: float
        """
        return self.mz * self.charge

    @property
    def neutral_mass(self):
        """
        Returns the neutral mass of the fragment.

        :return: Neutral mass of the fragment.
        :rtype: float
        """
        return self.mass - self.charge * PROTON_MASS

    @property
    def label(self):
        """
        Returns the label of the fragment, e.g., b2, y3i, etc.

        :return: Label of the fragment.
        :rtype: str
        """
        return f"{'+' * self.charge}" \
               f"{self.ion_type}" \
               f"{self.parent_number}" \
               f"{'i' if self.internal else ''}" \
               f"{self.number if self.internal else ''}"


def create_fragments(sequence: str, ion_types: Union[List[str], str], charges: Union[List[int], int],
                     monoisotopic: bool = True, internal: bool = False) -> List[Fragment]:
    """
    Generate a list of fragments for a given peptide sequence.

    :param sequence: The peptide sequence (can be modified).
    :type sequence: str
    :param ion_types: A list of ion types to consider, e.g., ['b', 'y'], or a single ion type, e.g., 'b'.
    :type ion_types: Union[List[str], str]
    :param charges: A list of charge states for the fragments, or a single charge state.
    :type charges: Union[List[int], int]
    :param monoisotopic: If True, use monoisotopic masses. If False, use average masses. Default is [True].
    :type monoisotopic: bool
    :param internal: If True, include internal fragments. If False, only include terminal fragments. Default is [False].
    :type internal: bool

    :return: A list of Fragment objects representing the generated fragments.
    :rtype: List[Fragment]

    .. code-block:: python

        >>> frags = create_fragments("[1.0]P(2.0)E(3.0)[4.0]", 'y', 1, monoisotopic=True, internal=False)
        >>> len(frags)
        2

        >>> fragments[0].sequence
        '[1.0]P(2.0)E(3.0)[4.0]'

        >>> fragments[0].mz
        255.11319808729002

        >>> fragments[1].sequence
        'E(3.0)[4.0]'

        >>> fragments[1].mz
        155.06043423844

    """

    if isinstance(ion_types, str):
        ion_types = [ion_types]

    if isinstance(charges, int):
        charges = [charges]

    fragments = []
    for ion_type in ion_types:
        seq_length = len(strip_modifications(sequence))
        for number, fragment_sequence in enumerate(create_fragment_sequences(sequence=sequence, ion_type=ion_type)):

            frag_number = seq_length - number
            for charge in charges:
                mz = calculate_mz(fragment_sequence, charge=charge, ion_type=ion_type, monoisotopic=monoisotopic)
                fragments.append(Fragment(sequence=fragment_sequence, mz=mz, charge=charge, ion_type=ion_type,
                                          number=frag_number, internal=False, parent_number=frag_number))
            if internal is True:
                fragment_seq_len = len(strip_modifications(fragment_sequence))
                for internal_number, internal_fragment_sequence in enumerate(
                        create_fragment_internal_sequences(sequence=fragment_sequence, ion_type=ion_type)):

                    if internal_number == 0:  # skip first since its already created
                        continue

                    internal_frag_number = fragment_seq_len - internal_number

                    for charge in charges:
                        mz = calculate_mz(internal_fragment_sequence, charge=charge, ion_type=ion_type,
                                          monoisotopic=monoisotopic)
                        fragments.append(
                            Fragment(sequence=internal_fragment_sequence, mz=mz, charge=charge, ion_type=ion_type,
                                     number=internal_frag_number, internal=True, parent_number=frag_number))

    return fragments


def create_fragment_sequences(sequence: str, ion_type: str) -> List[str]:
    """
    Generate fragment sequences for a given peptide sequence based on the specified ion type.

    :param sequence: The peptide sequence.
    :type sequence: str
    :param ion_type: The type of ion for which fragments are generated.
    :type ion_type: str

    :raise ValueError: If the ion type is not one of 'a', 'b', 'c', 'x', 'y', or 'z'.

    :return: The fragment sequences.
    :rtype: List[str]

    .. code-block:: python

        >>> create_fragment_sequences("PEPTIDE", "a")
        ['PEPTIDE', 'PEPTID', 'PEPTI', 'PEPT', 'PEP', 'PE', 'P']

        >>> create_fragment_sequences("PEPTIDE", "y")
        ['PEPTIDE', 'EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E']

        >>> create_fragment_sequences("PEPTIDE", "o")
        Traceback (most recent call last):
        ValueError: Ion type o is invalid.

    """

    return list(sequence_trimmer(sequence, forward=is_forward(ion_type)))


def create_fragment_internal_sequences(sequence: str, ion_type: str) -> List[str]:
    """
    Generate internal fragment sequences for a given peptide sequence based on the specified ion type.

    :param sequence: The peptide sequence.
    :type sequence: str
    :param ion_type: The type of ion for which internal fragments are generated.
    :type ion_type: str

    :raise ValueError: If the ion type is not one of 'a', 'b', 'c', 'x', 'y', or 'z'.

    :return: The internal fragment sequences.
    :rtype: List[str]

    .. code-block:: python

        >>> create_fragment_internal_sequences("PEPTIDE", "a")
        ['PEPTIDE', 'EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E']

        >>> create_fragment_internal_sequences("PEPTIDE", "y")
        ['PEPTIDE', 'PEPTID', 'PEPTI', 'PEPT', 'PEP', 'PE', 'P']

        >>> create_fragment_internal_sequences("PEPTIDE", "o")
        Traceback (most recent call last):
        ValueError: Ion type o is invalid.

    """

    return list(sequence_trimmer(sequence, forward=not is_forward(ion_type)))


def sequence_trimmer(sequence: str, forward: bool) -> Generator[str, None, None]:
    """
    Generates amino acid sequences either from the front or the back.

    :param sequence: The peptide sequence.
    :type sequence: str
    :param forward: If True, generate sequences from the front. If False, generate sequences from the back.
    :type forward: bool

    :return: A generator of sequences.
    :rtype: Generator[str, None, None]

    .. code-block:: python

        >>> list(sequence_trimmer("PEPTIDE", True))
        ['PEPTIDE', 'EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E']

        >>> list(sequence_trimmer("PEPTIDE", False))
        ['PEPTIDE', 'PEPTID', 'PEPTI', 'PEPT', 'PEP', 'PE', 'P']

        # Works with modifications
        >>> list(sequence_trimmer("[Acetyl]PE(3.0)PTIDE", True))
        ['[Acetyl]PE(3.0)PTIDE', 'E(3.0)PTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E']

        >>> list(sequence_trimmer("PEP(1.0)TIDE[Amide]", False))
        ['PEP(1.0)TIDE[Amide]', 'PEP(1.0)TID', 'PEP(1.0)TI', 'PEP(1.0)T', 'PEP(1.0)', 'PE', 'P']

    """

    if forward is True:
        return _trim_from_start(sequence)

    return _trim_from_end(sequence)


def _trim_from_end(sequence: str) -> Generator[str, None, None]:
    """
    This function yields sub-sequences progressively by removing the last amino acid or modification
    from the given sequence until the sequence is empty. If the last character in the sequence is a
    closing parenthesis (indicating a modification), the entire modification is removed. This function
    generates the sequences associated with abc ions.

    :param sequence: The initial amino acid sequence.
    :type sequence: str

    :yield: The sequence with the back amino acid (and modification) removed.
    :rtype: Generator[str, None, None]

    .. code-block:: python

        >>> list(_trim_from_end('PET'))
        ['PET', 'PE', 'P']
        
        >>> list(_trim_from_end('[Acetyl]P(3.14)ET[Amide]'))
        ['[Acetyl]P(3.14)ET[Amide]', '[Acetyl]P(3.14)E', '[Acetyl]P(3.14)']

    """

    while sequence:
        yield sequence
        sequence = strip_c_term(sequence)


def _trim_from_start(sequence: str) -> Generator[str, None, None]:
    """
    This function yields sub-sequences progressively by removing the first amino acid or modification
    from the given sequence until the sequence is empty. If the first character(s) in the sequence
    indicate a modification (enclosed in parentheses), the entire modification is removed. This function
    generates the sequences associated with xyz ions.

    :param sequence: The initial sequence.
    :type sequence: str

    :yield: The sequence with the front amino acid (and modification) removed.
    :rtype: Generator[str, None, None]

    .. code-block:: python

        >>> list(_trim_from_start('PET'))
        ['PET', 'ET', 'T']

        >>> list(_trim_from_start('[Acetyl]P(3.14)ET[Amide]'))
        ['[Acetyl]P(3.14)ET[Amide]', 'ET[Amide]', 'T[Amide]']

    """
    while sequence:
        yield sequence
        sequence = strip_n_term(sequence)


def calculate_fragment_mz_values(sequence: str, types=('b', 'y'), max_charge=1, monoisotopic=True) \
        -> Generator[float, None, None]:
    """
    Generates fragments of a given amino acid sequence based on the specified ion types and charge states.

    :param sequence: The amino acid sequence to be fragmented.
    :type sequence: str
    :param types: A tuple containing the types of ions to be considered in the fragmentation. Defaults to ('b', 'y').
    :type types: tuple
    :param max_charge: The maximum charge to consider for the fragments. Defaults to 1.
    :type max_charge: int
    :param monoisotopic: Indicates whether to calculate monoisotopic mass. Defaults to True.
    :type monoisotopic: bool

    :raise ValueError: If the ion type is not one of 'a', 'b', 'c', 'x', 'y', or 'z'.

    :yield: The calculated m/z for the fragment.
    :rtype: float
    """

    for ion_type in types:
        for charge in range(1, max_charge + 1):
            yield from calculate_fragment_mz_series(sequence=sequence, ion_type=ion_type, charge=charge,
                                                    monoisotopic=monoisotopic)


def calculate_fragment_mz_series(sequence: str, ion_type='y', charge=1, monoisotopic=True) \
        -> Generator[float, None, None]:
    """
    Generates fragment series based on ion type and charge.

    :param sequence: The amino acid sequence to be fragmented.
    :type sequence: str
    :param ion_type: The type of ions to be considered in the fragmentation. Defaults to 'y'.
    :type ion_type: str
    :param charge: The charge to consider for the fragments. Defaults to 1.
    :type charge: int
    :param monoisotopic: Indicates whether to calculate monoisotopic mass. Defaults to True.
    :type monoisotopic: bool

    :raise ValueError: If the ion type is not one of 'a', 'b', 'c', 'x', 'y', or 'z'.

    :yield: The calculated m/z for the fragment.
    :rtype: float

    .. code-block:: python

        >>> list(calculate_fragment_mz_series("TIDE", ion_type="y", charge=1))
        [477.21911970781, 376.17144123940005, 263.08737726227, 148.06043423844]

        >>> list(calculate_fragment_mz_series("TIDE", ion_type="b", charge=2))
        [230.10791574544, 165.58661920145502, 108.07314768954, 51.531115700975]

        >>> list(calculate_fragment_mz_series("[1.0]T(-1.0)IDE(2.0)[-2.0]", ion_type="y", charge=1))
        [477.21911970781, 376.17144123940005, 263.08737726227, 148.06043423844]

        >>> list(calculate_fragment_mz_series("[Acetyl]TIDE", ion_type="b", charge=2))
        Traceback (most recent call last):
        ...
        ValueError: could not convert string to float: 'Acetyl'

    """

    # Loop through the sequence to generate fragments and yield m/z
    for sub_sequence in sequence_trimmer(sequence, forward=is_forward(ion_type)):
        yield calculate_mz(sequence=sub_sequence, charge=charge, ion_type=ion_type, monoisotopic=monoisotopic)
