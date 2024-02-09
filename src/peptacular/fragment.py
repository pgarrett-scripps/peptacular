"""
fragment.py contains functions for.... guess what... fragmenting amino acid sequences!

There are two main methods for generating fragment ions within this module:
- build_fragments() - Generates all fragment objects.
- fragment() - Generates all fragment ion mass values.

fragment() is the recommended way to generate fragment ions, since it is faster and requires less memory than
build_fragments(). Still build_fragments is useful if you want to generate Fragment objects, which contain much more
information about the fragment ions.
"""

from dataclasses import dataclass
from functools import lru_cache, cached_property
from itertools import chain
from typing import List, Generator, Union, Callable, Dict

from peptacular.constants import PROTON_MASS, ION_ADJUSTMENTS
from peptacular.mass import calculate_mz, calculate_mass
from peptacular.sequence import calculate_sequence_length, split_sequence
from peptacular.term.residue import strip_c_term_residue, strip_n_term_residue
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

        >>> frag = Fragment("PE(3.1415)PTIDE", 1, "b", 2, False, 2, True)
        >>> frag.sequence
        'PE(3.1415)PTIDE'
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
        >>> frag.monoisotopic
        True

    """

    sequence: str
    charge: int
    ion_type: str
    number: int
    internal: bool
    parent_number: int
    monoisotopic: bool
    isotope: int = 0
    loss: float = 0.0
    parent_sequence: str = ''
    aa_masses: Dict = None

    def __post_init__(self):
        """
        Validates the fragment after initialization.
        """

        if self.internal is False and self.number != self.parent_number:
            raise ValueError("Parent number must be equal to number for terminal fragments.")

    @cached_property
    def _start(self):
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

    @cached_property
    def start(self):

        if self.ion_type in 'I':
            return 0

        if self.parent_sequence != '' and self._start < 0:
            return self._start + calculate_sequence_length(self.parent_sequence)

        return self._start

    @cached_property
    def _end(self):
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

    @cached_property
    def end(self):

        if self.ion_type in 'I':
            return 0

        if self.parent_sequence != '' and self._end is None:
            return calculate_sequence_length(self.parent_sequence)

        if self.parent_sequence != '' and self._end < 0:
            return self._end + calculate_sequence_length(self.parent_sequence)

        return self._end

    @cached_property
    def mass(self):
        """
        Returns the mass of the fragment.

        :return: Mass of the fragment.
        :rtype: float
        """

        return calculate_mass(sequence=self.sequence, charge=self.charge, ion_type=self.ion_type,
                              monoisotopic=self.monoisotopic, uwpr_mass=False, isotope=self.isotope,
                              loss=self.loss, aa_masses=self.aa_masses)

    @cached_property
    def neutral_mass(self):
        """
        Returns the neutral mass of the fragment.

        :return: Neutral mass of the fragment.
        :rtype: float
        """

        return calculate_mass(sequence=self.sequence, charge=0, ion_type=self.ion_type,
                              monoisotopic=self.monoisotopic, uwpr_mass=False, isotope=self.isotope,
                              loss=self.loss, aa_masses=self.aa_masses)

    @cached_property
    def mz(self):
        """
        Returns the mass-to-charge (m/z) ratio of the fragment.

        :return: Mass-to-charge (m/z) ratio of the fragment.
        :rtype: float
        """

        return calculate_mz(sequence=self.sequence, charge=self.charge, ion_type=self.ion_type,
                            monoisotopic=self.monoisotopic, uwpr_mass=False, isotope=self.isotope,
                            loss=self.loss, aa_masses=self.aa_masses)

    @cached_property
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
               f"{self.number if self.internal else ''}" \
               f"{'(' + str(self.loss) + ')' if self.loss != 0.0 else ''}" \
               f"{'*' * self.isotope if self.isotope > 0 else ''}"

    def __iter__(self):
        # Include regular attributes
        for key, value in self.__dict__.items():
            yield key, value

        # Explicitly include cached properties
        yield 'start', self.start
        yield 'end', self.end
        yield 'mass', self.mass
        yield 'neutral_mass', self.neutral_mass
        yield 'mz', self.mz
        yield 'label', self.label

    # Method to convert the object into a dictionary including cached properties
    def to_dict(self):
        return dict(self)


def build_fragments(sequence: str, ion_types: Union[List[str], str], charges: Union[List[int], int],
                    monoisotopic: bool = True, internal: bool = False, isotopes: Union[List[int], int] = 0,
                    losses: Union[List[float], float] = 0.0, aa_masses: Dict = None,
                    immonium: bool = False) -> List[Fragment]:
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
    :param isotopes: A list of isotope offsets to consider, or a single isotope offset, default is [0].
    :type isotopes: Union[List[int], int]
    :param losses: A list of neutral losses to consider, or a single neutral loss, default is [0.0].
    :type losses: Union[List[float], float]
    :param aa_masses: A dictionary of amino acid masses, default is None.
    :type aa_masses: Dict
    :param immonium: If True, include immonium ions, default is [False].
    :type immonium: bool

    :return: List of Fragment objects.
    :rtype: List[Fragment]

    .. code-block:: python

        # Create a list of y ions for a modified peptide.
        >>> frags = build_fragments("[1.0]P(2.0)E(3.0)[4.0]", 'y', 1, monoisotopic=True, internal=False)
        >>> len(frags)
        2

        >>> frags[0].sequence
        '[1.0]P(2.0)E(3.0)[4.0]'

        >>> frags[1].sequence
        'E(3.0)[4.0]'

    """

    if isinstance(ion_types, str):
        ion_types = [ion_types]

    if isinstance(charges, int):
        charges = [charges]

    if isinstance(isotopes, int):
        isotopes = [isotopes]

    if isinstance(losses, float) or isinstance(losses, int):
        losses = [losses]

    if len(isotopes) == 0 or len(losses) == 0 or len(ion_types) == 0 or len(charges) == 0:
        raise ValueError("No isotopes, losses, ion types, or charges provided.")

    def _fragment(t: str, c: int):
        # Loop through the sequence to generate fragments and yield m/z
        seq_length = calculate_sequence_length(sequence)
        for number, frag in enumerate(build_fragment_sequences(sequence, t)):
            frag_number = seq_length - number
            for iso in isotopes:
                for loss in losses:
                    yield Fragment(sequence=frag, charge=c, ion_type=t, number=frag_number,
                                   internal=False, parent_number=frag_number, monoisotopic=monoisotopic,
                                   isotope=iso, loss=loss, parent_sequence=sequence, aa_masses=aa_masses)

            if internal is True:
                fragment_seq_len = calculate_sequence_length(frag)
                for internal_number, internal_frag in \
                        enumerate(build_fragment_sequences(sequence=frag, ion_type=t, internal=True)):

                    if internal_number == 0:  # skip first since its already created
                        continue

                    internal_frag_number = fragment_seq_len - internal_number
                    for iso in isotopes:
                        for loss in losses:
                            yield Fragment(sequence=internal_frag, charge=c, ion_type=t,
                                           number=internal_frag_number, internal=True, parent_number=frag_number,
                                           monoisotopic=monoisotopic, isotope=iso, loss=loss, parent_sequence=sequence,
                                           aa_masses=aa_masses)

    immonium_frags = []
    if immonium:
        for aa in set(split_sequence(sequence)):
            immonium_frag = Fragment(sequence=aa, charge=1, ion_type='I', number=0, internal=False, parent_number=0,
                                     monoisotopic=monoisotopic, isotope=0, loss=0.0, parent_sequence=sequence,
                                     aa_masses=aa_masses)
            immonium_frags.append(immonium_frag)

    fragments = list(chain.from_iterable(_fragment(t, c) for t in ion_types for c in charges)) + immonium_frags
    return fragments


def build_fragment_sequences(sequence: str, ion_type: str, internal: bool = False) -> List[str]:
    """
    Creates a list of fragment sequences for a given input `sequence` and `ion_type`.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param ion_type: The ion type, e.g., 'b', 'y', 'a', 'c', 'x', or 'z'.
    :type ion_type: str
    :param internal: If True, build internal fragments. If False, build normal fragments, default is [False].
    :type internal: bool

    :return: List of fragment sequences.
    :rtype: List[str]

    .. code-block:: python

        >>> build_fragment_sequences("PEPTIDE", 'y')
        ['PEPTIDE', 'EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E']

        >>> build_fragment_sequences("PEPTIDE", 'b', True)
        ['PEPTIDE', 'EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E']

        >>> build_fragment_sequences("PEPTIDE", 'b')
        ['PEPTIDE', 'PEPTID', 'PEPTI', 'PEPT', 'PEP', 'PE', 'P']

        >>> build_fragment_sequences("PEPTIDE", 'y', True)
        ['PEPTIDE', 'PEPTID', 'PEPTI', 'PEPT', 'PEP', 'PE', 'P']

    """
    if internal is True:
        return sequence_trimmer(sequence, forward=not is_forward(ion_type))
    return sequence_trimmer(sequence, forward=is_forward(ion_type))


@lru_cache(maxsize=2)
def sequence_trimmer(sequence: str, forward: bool) -> List[str]:
    """
    A generator that yields all subsequences of the input `sequence`, by trimming from the front or back.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param forward: If True, generate sequences from the front. If False, generate sequences from the back.
    :type forward: bool

    :return: The trimmed sequeunces.
    :rtype: List[str]

    .. code-block:: python

        >>> sequence_trimmer("PEPTIDE", True)
        ['PEPTIDE', 'EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E']

        >>> sequence_trimmer("PEPTIDE", False)
        ['PEPTIDE', 'PEPTID', 'PEPTI', 'PEPT', 'PEP', 'PE', 'P']

        # Works with modifications
        >>> sequence_trimmer("[Acetyl]PE(3.0)PTIDE", True)
        ['[Acetyl]PE(3.0)PTIDE', 'E(3.0)PTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E']

        >>> sequence_trimmer("PEP(1.0)TIDE[Amide]", False)
        ['PEP(1.0)TIDE[Amide]', 'PEP(1.0)TID', 'PEP(1.0)TI', 'PEP(1.0)T', 'PEP(1.0)', 'PE', 'P']

    """

    def trim_from_end(seq: str) -> Generator[str, None, None]:
        while seq:
            yield seq
            seq = strip_c_term_residue(seq)

    def trim_from_start(seq: str) -> Generator[str, None, None]:
        while seq:
            yield seq
            seq = strip_n_term_residue(seq)

    if forward is True:
        return list(trim_from_start(sequence))

    return list(trim_from_end(sequence))


# TODO: Add new fragment support to fragment function
def fragment(sequence: str, ion_types: Union[str, List[str]] = 'y', charges: Union[int, List[int]] = 1,
             monoisotopic: bool = True, internal: bool = False, mz: bool = True) -> List[float]:
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
    :param mz: Indicates whether to return m/z values or mass values, defaults to [True].
    :type mz: bool

    :raise ValueError: If the ion type is not one of 'a', 'b', 'c', 'x', 'y', or 'z'.

    :return: The fragemnt ion m/z values.
    :rtype: List[float]

    .. code-block:: python

        # Fragments are returned from largest to smallest.
        >>> fragment("TIDE", ion_types="y", charges=1)
        [477.21911970781, 376.17144123940005, 263.08737726227, 148.06043423844]
        >>> fragment("TIDE", ion_types="b", charges=2)
        [230.10791574544, 165.58661920145502, 108.07314768954, 51.531115700975]

        # When using multiple ion types and charge states the fragments will be generated first by ion_types
        # and second by charge state. For example, the following will generate the following fragments:
        # [y2+, y1+, y2++, y1++, b2+, b1+, b2++, b1++]
        >>> fragment("PE", ion_types=["y", "b"], charges=[1,2])[:4] # First 4 fragments
        [245.11319808729002, 148.06043423844, 123.06023727703, 74.53385535260499]
        >>> fragment("PE", ion_types=["y", "b"], charges=[1,2])[4:] # Last 4 fragments
        [227.10263340359, 98.06004031562, 114.05495493517999, 49.533658391195004]

        # Can also include modifications in the sequence.
        >>> fragment("[1.0]T(-1.0)IDE(2.0)[-2.0]", ion_types="y", charges=1)
        [477.21911970781, 376.17144123940005, 263.08737726227, 148.06043423844]

        # If the sequence contains a modification that is not a number, an error will be raised.
        >>> fragment("[Acetyl]TIDE", ion_types="b", charges=2)
        Traceback (most recent call last):
        ...
        ValueError: could not convert string to float: 'Acetyl'

        # Can generate internal fragments too, but they are typically not very useful.
        # order: sequences -> [PEP, PE, P,  PEP,  E,  P], ions -> [y3, y3i, y3i, y2, y2i, y1]
        >>> fragment("PEP", ion_types="y", charges=1, internal=True)
        [342.16596193614004, 245.11319808729002, 116.07060499932, 245.11319808729002, 148.06043423844, 116.07060499932]

    """

    if isinstance(ion_types, str):
        ion_types = [ion_types]

    if isinstance(charges, int):
        charges = [charges]

    if internal is True:
        return _slow_fragment(sequence=sequence, ion_types=ion_types, charges=charges, monoisotopic=monoisotopic,
                              mz=mz, internal=internal)

    return _fast_fragment(sequence=sequence, ion_types=ion_types, charges=charges, monoisotopic=monoisotopic, mz=mz)


def _slow_fragment(sequence: str, ion_types: List[str], charges: List[int],
                   monoisotopic: bool, mz: bool, internal: bool) -> List[float]:
    """
    Original fragment function that is slow, but works with internal fragments.
    """

    def _fragment(t: str, c: int, f: Callable):
        # Loop through the sequence to generate fragments and yield m/z
        for sub_sequence in build_fragment_sequences(sequence, ion_type=t):
            if mz is True:
                yield f(sequence=sub_sequence, charge=c, ion_type=t, monoisotopic=monoisotopic)
            if internal is True:
                for i, internal_sub_sequence in enumerate(
                        build_fragment_sequences(sequence=sub_sequence, ion_type=t, internal=True)):
                    if i == 0:  # skip first since its already created
                        continue
                    yield f(sequence=internal_sub_sequence, charge=c, ion_type=t, monoisotopic=monoisotopic)

    func = calculate_mz if mz else calculate_mass
    fragments = list(chain.from_iterable(_fragment(t, c, func) for t in ion_types for c in charges))

    return fragments


def _fast_fragment(sequence: str, ion_types: List[str], charges: List[int],
                   monoisotopic: bool, mz: bool) -> List[float]:
    """
    Faster version of _slow_fragment that does not support internal fragments.
    """
    aas = split_sequence(sequence)

    masses = [calculate_mass(aa, 0, 'b', monoisotopic) for aa in aas]

    fragments = []

    for ion_type in ion_types:
        for charge in charges:

            if ion_type in ['a', 'b', 'c']:
                cum_sum = [sum(masses[:i + 1]) for i in range(len(masses))]
            elif ion_type in ['x', 'y', 'z']:
                cum_sum = [sum(masses[-(i + 1):]) for i in range(len(masses))]
            else:
                raise ValueError(f"Invalid ion type: {ion_type}")

            cum_sum = [x + charge * PROTON_MASS + ION_ADJUSTMENTS[ion_type] for x in cum_sum]

            if mz and charge > 0:
                cum_sum = [x / charge for x in cum_sum]

            fragments.extend(cum_sum[::-1])

    return fragments
