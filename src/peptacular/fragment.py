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
from functools import cached_property
from typing import List, Generator, Union

from peptacular.constants import FORWARD_ION_TYPES, BACKWARD_ION_TYPES, INTERNAL_ION_TYPES, TERMINAL_ION_TYPES, \
    IMMONIUM_ION_TYPES
from peptacular.mass import calculate_mz, calculate_mass
from peptacular.sequence.sequence import calculate_sequence_length, pop_modifications, span_to_sequence
from peptacular.spans import build_non_enzymatic_spans, build_right_semi_spans, build_left_semi_spans, Span
from peptacular.term.modification import pop_term_modifications

ChargeType = Union[List[int], int]
IsotopeType = Union[List[int], int]
LossType = Union[List[float], float]
IonTypeType = Union[List[str], str]


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

        >>> frag = Fragment(1, "b", 0, 7, True, 0, 0.0, 'PE[3.1415]PTIDES')
        >>> frag.sequence
        'PE[3.1415]PTIDE'
        >>> frag.charge
        1
        >>> frag.ion_type
        'b'
        >>> frag.internal
        False
        >>> frag.monoisotopic
        True

    """

    charge: int
    ion_type: str
    start: int
    end: int
    monoisotopic: bool
    isotope: int
    loss: float
    parent_sequence: str

    @cached_property
    def internal(self):
        return self.start != 0 and self.end != calculate_sequence_length(self.parent_sequence)

    @cached_property
    def sequence(self):
        span = (self.start, self.end, 0)
        if self.ion_type in IMMONIUM_ION_TYPES:  # Immonium fragments shouldn't contain term sequence
            sequence, _, _ = pop_term_modifications(self.parent_sequence)
            return span_to_sequence(sequence, span)
        return span_to_sequence(self.parent_sequence, span)

    @cached_property
    def mass(self):
        """
        Returns the mass of the fragment.

        :return: Mass of the fragment.
        :rtype: float
        """

        return calculate_mass(sequence=self.sequence, charge=self.charge, ion_type=self.ion_type,
                              monoisotopic=self.monoisotopic, isotope=self.isotope,
                              loss=self.loss)

    @cached_property
    def neutral_mass(self):
        """
        Returns the neutral mass of the fragment.

        :return: Neutral mass of the fragment.
        :rtype: float
        """

        return calculate_mass(sequence=self.sequence, charge=0, ion_type=self.ion_type,
                              monoisotopic=self.monoisotopic, isotope=self.isotope,
                              loss=self.loss)

    @cached_property
    def mz(self):
        """
        Returns the mass-to-charge (m/z) ratio of the fragment.

        :return: Mass-to-charge (m/z) ratio of the fragment.
        :rtype: float
        """

        return calculate_mz(sequence=self.sequence, charge=self.charge, ion_type=self.ion_type,
                            monoisotopic=self.monoisotopic, isotope=self.isotope,
                            loss=self.loss)

    @cached_property
    def label(self):
        """
        Returns the label of the fragment, e.g., b2, y3i, etc.

        :return: Label of the fragment.
        :rtype: str
        """

        if self.ion_type in FORWARD_ION_TYPES:
            number = self.end
        elif self.ion_type in BACKWARD_ION_TYPES:
            number = calculate_sequence_length(self.parent_sequence) - self.start
        elif self.ion_type in INTERNAL_ION_TYPES:
            number = f'{self.start}-{self.end}'
        elif self.ion_type == 'I':
            number = self.start
        else:
            raise ValueError('Wrong Ion Type')

        return f"{'+' * self.charge}" \
               f"{self.ion_type}" \
               f"{number}" \
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
        yield 'internal', self.internal
        yield 'sequence', self.sequence

    # Method to convert the object into a dictionary including cached properties
    def to_dict(self):
        return {
            key: value for key, value in self
        }


def _build_fragments(spans: List[Span], ion_types: List[str], charges: List[int], losses: List[float],
                     isotopes: List[int], monoisotopic: bool, sequence: str) -> Generator[Fragment, None, None]:
    for span in spans:
        for ion_type in ion_types:
            for iso in isotopes:
                for loss in losses:
                    for c in charges:
                        yield Fragment(charge=c, ion_type=ion_type,
                                       start=span[0], end=span[1], monoisotopic=monoisotopic, isotope=iso,
                                       loss=loss, parent_sequence=sequence)


def build_internal_fragments(sequence: str, ion_types: IonTypeType, charges: ChargeType,
                             monoisotopic: bool = True, isotopes: IsotopeType = 0,
                             losses: LossType = 0.0) -> List[Fragment]:
    """
    .. code-block:: python

        >>> [f.sequence for f in build_internal_fragments('TIDE', 'by', 1)]
        ['I', 'ID', 'D']

        >>> [round(f.mz, 3) for f in build_internal_fragments('TIDE', 'by', 1)]
        [131.094, 246.121, 133.037]

        >>> [f.label for f in build_internal_fragments('TIDE', 'by', 1)]
        ['+by1-2', '+by1-3', '+by2-3']

    """

    if isinstance(ion_types, str):
        ion_types = [ion_types]

    if isinstance(charges, int):
        charges = [charges]

    if isinstance(isotopes, int):
        isotopes = [isotopes]

    if isinstance(losses, float) or isinstance(losses, int):
        losses = [losses]

    if len(isotopes) == 0 or len(losses) == 0 or len(charges) == 0:
        raise ValueError("No isotopes, losses, ion types, or charges provided.")

    stripped_sequence, mods = pop_modifications(sequence)
    spans = build_non_enzymatic_spans((0, len(stripped_sequence), 0))
    internal_spans = [span for span in spans if span[0] != 0 and span[1] != len(stripped_sequence)]
    return list(_build_fragments(internal_spans, ion_types, charges, losses, isotopes, monoisotopic, sequence))


def build_immonium_fragments(sequence: str, charges: ChargeType, monoisotopic: bool = True, isotopes: IsotopeType = 0,
                             losses: LossType = 0.0) -> List[Fragment]:
    """
    Build immonium ions for a given sequence.

    :param sequence: Peptide sequence.
    :type sequence: str
    :param charges: Charges to consider.
    :type charges: Union[List[int], int]
    :param monoisotopic: Whether to use monoisotopic masses.
    :type monoisotopic: bool
    :param isotopes: Isotopes to consider.
    :type isotopes: Union[List[int], int]
    :param losses: Losses to consider.
    :type losses: Union[List[float], float]

    :return: List of immonium ions.
    :rtype: List[Fragment]

    .. code-block:: python

        >>> [f.sequence for f in build_immonium_fragments('TIDE', 1)]
        ['T', 'I', 'D', 'E']

        >>> [round(f.mz,3) for f in build_immonium_fragments('TIDE', 1)]
        [74.06, 86.096, 88.039, 102.055]

        >>> [f.label for f in build_immonium_fragments('TIDE', 1)]
        ['+I0', '+I1', '+I2', '+I3']

    """
    if isinstance(charges, int):
        charges = [charges]

    if isinstance(isotopes, int):
        isotopes = [isotopes]

    if isinstance(losses, float) or isinstance(losses, int):
        losses = [losses]

    sequence, _, _ = pop_term_modifications(sequence)
    spans = [(i, i + 1, 0) for i in range(len(sequence))]
    return list(_build_fragments(spans, ['I'], charges, losses, isotopes, monoisotopic, sequence))


def _build_forward_fragments(sequence: str, ion_types: List[str], charges: List[int], monoisotopic: bool,
                             isotopes: List[int], losses: List[float]) -> List[Fragment]:
    start_span = (0, calculate_sequence_length(sequence), 0)
    spans = [start_span] + build_left_semi_spans(start_span)
    return list(_build_fragments(spans, ion_types, charges, losses, isotopes, monoisotopic, sequence))


def _build_backward_fragments(sequence: str, ion_types: List[str], charges: List[int], monoisotopic: bool,
                              isotopes: List[int], losses: List[float]) -> Union[float, List[Fragment]]:
    start_span = (0, calculate_sequence_length(sequence), 0)
    spans = [start_span] + build_right_semi_spans(start_span)
    return list(_build_fragments(spans, ion_types, charges, losses, isotopes, monoisotopic, sequence))


def build_terminal_fragments(sequence: str, ion_types: IonTypeType, charges: ChargeType, monoisotopic: bool = True,
                             isotopes: IsotopeType = 0, losses: LossType = 0.0) -> List[Fragment]:
    """

    .. code-block:: python

        >>> frags = build_terminal_fragments('TIDE', 'b', 2)
        >>> [round(f.mass,3) for f in frags]
        [460.216, 331.173, 216.146, 103.062]

        >>> [round(f.mz,3) for f in build_terminal_fragments("TIDE", ion_types="y", charges=1)]
        [477.219, 376.171, 263.087, 148.06]

        >>> [round(f.mz,3) for f in build_terminal_fragments("TIDE", ion_types="b", charges=2)]
        [230.108, 165.587, 108.073, 51.531]

        >>> [round(f.mz,3) for f in build_terminal_fragments("T[10]IDE", ion_types="y", charges=1)]
        [487.219, 376.171, 263.087, 148.06]

    """

    if isinstance(ion_types, str):
        ion_types = [ion_types]

    if isinstance(charges, int):
        charges = [charges]

    if isinstance(isotopes, int):
        isotopes = [isotopes]

    if isinstance(losses, float) or isinstance(losses, int):
        losses = [losses]

    forward_ions = [ion for ion in ion_types if ion in FORWARD_ION_TYPES]
    backward_ions = [ion for ion in ion_types if ion in BACKWARD_ION_TYPES]

    frags = []
    frags.extend(_build_forward_fragments(sequence, forward_ions, charges, monoisotopic, isotopes, losses))
    frags.extend(_build_backward_fragments(sequence, backward_ions, charges, monoisotopic, isotopes, losses))

    return frags


def build_fragments(sequence: str, ion_types: IonTypeType, charges: ChargeType, monoisotopic: bool = True,
                    isotopes: IsotopeType = 0, losses: LossType = 0.0) -> List[Fragment]:
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
    :param isotopes: A list of isotope offsets to consider, or a single isotope offset, default is [0].
    :type isotopes: Union[List[int], int]
    :param losses: A list of neutral losses to consider, or a single neutral loss, default is [0.0].
    :type losses: Union[List[float], float]

    :return: List of Fragment objects.
    :rtype: List[Fragment]

    .. code-block:: python

        # Create a list of y ions for a modified peptide.
        >>> frags = build_fragments("[1.0]-P[2.0]E[3.0]-[4.0]", 'y', 1, monoisotopic=True)
        >>> len(frags)
        2

        >>> frags[0].sequence
        '[1.0]-P[2.0]E[3.0]-[4.0]'

        >>> frags[1].sequence
        'E[3.0]-[4.0]'

    """

    if isinstance(ion_types, str):
        ion_types = [ion_types]

    if isinstance(charges, int):
        charges = [charges]

    if isinstance(isotopes, int):
        isotopes = [isotopes]

    if isinstance(losses, float) or isinstance(losses, int):
        losses = [losses]

    terminal_fragment_types = [i for i in ion_types if i in TERMINAL_ION_TYPES]
    internal_fragment_types = [i for i in ion_types if i in INTERNAL_ION_TYPES]
    immonium = 'I' in ion_types

    frags = []
    if terminal_fragment_types:
        frags.extend(build_terminal_fragments(sequence, terminal_fragment_types, charges,
                                              monoisotopic, isotopes, losses))

    if internal_fragment_types:
        frags.extend(build_internal_fragments(sequence, internal_fragment_types, charges,
                                              monoisotopic, isotopes, losses))

    if immonium:
        frags.extend(build_immonium_fragments(sequence, charges, monoisotopic, isotopes, losses))

    return frags
