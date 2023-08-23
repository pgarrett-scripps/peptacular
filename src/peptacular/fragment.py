from dataclasses import dataclass
from typing import List, Generator, Union

from peptacular.constants import PROTON_MASS
from peptacular.mass import calculate_mz
from peptacular.sequence import calculate_sequence_length
from peptacular.term import strip_c_term, strip_n_term
from peptacular.util import is_forward


@dataclass(frozen=True)
class Fragment:
    """
    A dataclass for representing a peptide fragment ion.

    :ivar str sequence: The fragment sequence, including any modifications.
    :ivar float mz: The fragment mass-to-charge (m/z) ratio
    :ivar int charge: The fragment charge state.
    :ivar str ion_type: The fragment ion type.
    :ivar int number: The fragment ion number, for abc ions this refers to the number of residues starting from the
                      N-terminus, for xyz ions this refers to the number of residues starting from the C-terminus.
    :ivar bool internal: True if the fragment is an internal fragment, False if it is a terminal fragment.
    :ivar int parent_number: The number of the parent fragment from which the internal fragment is derived. For
                             terminal fragments (internal=False) the number and parent number should always be equal.

    .. code-block:: python

        >>> frag = Fragment("PE(3.1415)PTIDE", 123.456, 1, "b", 2, False, 2)
        >>> frag.sequence
        'PE(3.1415)PTIDE'
        >>> frag.mz
        123.456
        >>> frag.charge
        1
        >>> frag.ion_type
        'b'
        >>> frag.number
        2
        >>> frag.internal
        False
        >>> frag.parent_number
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


def build_fragments(sequence: str, ion_types: Union[List[str], str], charges: Union[List[int], int],
                    monoisotopic: bool = True, internal: bool = False) -> List[Fragment]:
    """
    Builds all Fragment objects or a given input 'sequence'.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param ion_types: A list of ion types to consider, e.g., ['b', 'y'], or a single ion type, e.g., 'b'.
    :type ion_types: Union[List[str], str]
    :param charges: A list of charge states for the fragments, or a single charge state.
    :type charges: Union[List[int], int]
    :param monoisotopic: If True, use monoisotopic masses. If False, use average masses, default is [True].
    :type monoisotopic: bool
    :param internal: If True, include internal fragments. If False, only include terminal fragments, default is [False].
    :type internal: bool

    :return: List of Fragment objects.
    :rtype: List[Fragment]

    .. code-block:: python

        # Create a list of y ions for a modified peptide.
        >>> frags = build_fragments("[1.0]P(2.0)E(3.0)[4.0]", 'y', 1, monoisotopic=True, internal=False)
        >>> len(frags)
        2

        >>> frags[0].sequence
        '[1.0]P(2.0)E(3.0)[4.0]'

        >>> frags[0].mz
        255.11319808729002

        >>> frags[1].sequence
        'E(3.0)[4.0]'

        >>> frags[1].mz
        155.06043423844

    """

    if isinstance(ion_types, str):
        ion_types = [ion_types]

    if isinstance(charges, int):
        charges = [charges]

    fragments = []
    for ion_type in ion_types:
        seq_length = calculate_sequence_length(sequence)
        for number, fragment_sequence in enumerate(build_fragment_sequences(sequence=sequence, ion_type=ion_type)):

            frag_number = seq_length - number
            for charge in charges:
                mz = calculate_mz(fragment_sequence, charge=charge, ion_type=ion_type, monoisotopic=monoisotopic)
                fragments.append(Fragment(sequence=fragment_sequence, mz=mz, charge=charge, ion_type=ion_type,
                                          number=frag_number, internal=False, parent_number=frag_number))
            if internal is True:
                fragment_seq_len = calculate_sequence_length(fragment_sequence)
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


def build_fragment_sequences(sequence: str, ion_type: str) -> List[str]:
    """
    Builds all fragment subsequences for the input `sequence`.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param ion_type: The type of ion for which fragments are generated.
    :type ion_type: str

    :raise ValueError: If the ion type is not one of 'a', 'b', 'c', 'x', 'y', or 'z'.

    :return: The fragment sequences.
    :rtype: List[str]

    .. code-block:: python

        >>> build_fragment_sequences("PEPTIDE", "a")
        ['PEPTIDE', 'PEPTID', 'PEPTI', 'PEPT', 'PEP', 'PE', 'P']

        >>> build_fragment_sequences("PEPTIDE", "y")
        ['PEPTIDE', 'EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E']

        >>> build_fragment_sequences("PEPTIDE", "o")
        Traceback (most recent call last):
        ValueError: Ion type o is invalid.

    """

    return list(sequence_trimmer(sequence, forward=is_forward(ion_type)))


def create_fragment_internal_sequences(sequence: str, ion_type: str) -> List[str]:
    """
   Builds all internal fragment subsequences for the input `sequence`.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param ion_type: The type of ion for which internal fragments are generated.
    :type ion_type: str

    :raise ValueError: If the ion type is not one of 'a', 'b', 'c', 'x', 'y', or 'z'.

    :return: The internal fragment subsequences.
    :rtype: List[str]

    .. code-block:: python

        >>> create_fragment_internal_sequences("PEPTIDE", "a")
        ['PEPTIDE', 'EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E']

        >>> create_fragment_internal_sequences("PEPTIDE", "y")
        ['PEPTIDE', 'PEPTID', 'PEPTI', 'PEPT', 'PEP', 'PE', 'P']

        >>> create_fragment_internal_sequences("PEPTIDE", "o")
        Traceback (most recent call last):
            ...
        ValueError: Ion type o is invalid.

    """

    return list(sequence_trimmer(sequence, forward=not is_forward(ion_type)))


def sequence_trimmer(sequence: str, forward: bool) -> Generator[str, None, None]:
    """
    A generator that yields all subsequences of the input `sequence`, by trimming from the front or back.

    :param sequence: The amino acid sequence, which can include modifications.
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

    :param sequence: The amino acid sequence, which can include modifications.
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

    :param sequence: The amino acid sequence, which can include modifications.
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


def fragment(sequence: str, ion_types: Union[str, List[str]] = 'y', charges: Union[int, List[int]] = 1,
             monoisotopic: bool = True, internal: bool = False) -> Generator[float, None, None]:
    """
    Generates fragment ion mz values.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param ion_types: A list of ion types to consider, e.g., ['b', 'y'], or a single ion type, e.g., 'b'.
    :type ion_types: Union[List[str], str]
    :param charges: A list of charge states for the fragments, or a single charge state.
    :type charges: Union[List[int], int]
    :param monoisotopic: Indicates whether to calculate monoisotopic mass, defaults to [True].
    :type monoisotopic: bool
    :param internal: Indicates whether to calculate internal fragments, defaults to [False].
    :type internal: bool

    :raise ValueError: If the ion type is not one of 'a', 'b', 'c', 'x', 'y', or 'z'.

    :yield: The calculated m/z for the fragment.
    :rtype: float

    .. code-block:: python

        >>> list(fragment("TIDE", ion_types="y", charges=1))
        [477.21911970781, 376.17144123940005, 263.08737726227, 148.06043423844]

        >>> list(fragment("TIDE", ion_types="b", charges=2))
        [230.10791574544, 165.58661920145502, 108.07314768954, 51.531115700975]

        >>> list(fragment("[1.0]T(-1.0)IDE(2.0)[-2.0]", ion_types="y", charges=1))
        [477.21911970781, 376.17144123940005, 263.08737726227, 148.06043423844]

        >>> list(fragment("[Acetyl]TIDE", ion_types="b", charges=2))
        Traceback (most recent call last):
        ...
        ValueError: could not convert string to float: 'Acetyl'

        # Can generate internal fragments too, but they are not very useful.
        # Here is the order of fragments generated:
        # [PEP, PE, P,  PEP,  E,  P]
        # [y3, y3i, y3i, y2, y2i, y1]
        >>> list(fragment("PEP", ion_types="y", charges=1, internal=True))
        [342.16596193614004, 245.11319808729002, 116.07060499932, 245.11319808729002, 148.06043423844, 116.07060499932]

    """

    if isinstance(ion_types, str):
        ion_types = [ion_types]

    if isinstance(charges, int):
        charges = [charges]

    def _fragment(t: str, c: int):
        # Loop through the sequence to generate fragments and yield m/z
        for sub_sequence in sequence_trimmer(sequence, forward=is_forward(t)):
            yield calculate_mz(sequence=sub_sequence, charge=c, ion_type=t, monoisotopic=monoisotopic)
            if internal is True:
                for i, internal_sub_sequence in enumerate(
                        create_fragment_internal_sequences(sequence=sub_sequence, ion_type=t)):
                    if i == 0:  # skip first since its already created
                        continue
                    yield calculate_mz(sequence=internal_sub_sequence, charge=c, ion_type=t,
                                       monoisotopic=monoisotopic)

    for ion_type in ion_types:
        for charge in charges:
            yield from _fragment(ion_type, charge)
