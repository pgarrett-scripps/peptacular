"""
fragment.py contains functions for fragmenting peptides
"""
from dataclasses import dataclass
from functools import cached_property
from typing import List, Generator, Union, Dict

from peptacular.proforma.proforma import ProFormaAnnotation
from peptacular.constants import FORWARD_ION_TYPES, BACKWARD_ION_TYPES, INTERNAL_ION_TYPES, TERMINAL_ION_TYPES, \
    NEUTRON_MASS
from peptacular.mass import mz, mass, comp, adjust_mass, adjust_mz
from peptacular.sequence.sequence import sequence_length, span_to_sequence, sequence_to_annotation, split
from peptacular.spans import build_non_enzymatic_spans, build_right_semi_spans, build_left_semi_spans, Span
from peptacular.types import IonTypeType, IsotopeType, LossType, ChargeType


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
    parent_sequence: Union[str, ProFormaAnnotation]
    mass: float
    neutral_mass: float
    mz: float
    sequence: str
    unmod_sequence: str
    internal: bool

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
            number = len(self.unmod_sequence) - self.start
        elif self.ion_type in INTERNAL_ION_TYPES:
            number = f'{self.start}-{self.end}'
        elif self.ion_type == 'i':
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
            if key == 'parent_sequence' and isinstance(value, ProFormaAnnotation):
                yield key, value.serialize()
            else:
                yield key, value

        # Explicitly include cached properties
        yield 'label', self.label

    # Method to convert the object into a dictionary including cached properties
    def to_dict(self):

        return {
            key: value for key, value in self
        }


def _build_fragments(spans: List[Span],
                     ion_types: List[str],
                     charges: List[int],
                     losses: List[float],
                     isotopes: List[int],
                     monoisotopic: bool,
                     sequence: Union[str, ProFormaAnnotation],
                     return_type: str,
                     mass_components: List[float] = None) -> Union[Generator[Fragment, None, None], Generator[float, None, None]]:

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    if mass_components is None:
        components = annotation.split()
        mass_components = [mass(sequence=component, charge=0, ion_type='n', monoisotopic=monoisotopic) for component in components]

    for span in spans:
        base_mass = adjust_mass(sum(mass_components[span[0]:span[1]]), charge=0, ion_type='n', monoisotopic=monoisotopic)
        base_annotation = annotation.slice(span[0], span[1])
        base_sequence = base_annotation.serialize()
        base_unmod_sequence = base_annotation.sequence
        for ion_type in ion_types:
            base_fragment_mass = adjust_mass(base_mass, charge=0, ion_type=ion_type, monoisotopic=monoisotopic)
            for iso in isotopes:
                for loss in losses:
                    for c in charges:
                        fragment_neutral_mass = adjust_mass(base_fragment_mass, charge=0, ion_type='n', monoisotopic=monoisotopic, isotope=iso, loss=loss)
                        fragment_mass = adjust_mass(fragment_neutral_mass, charge=c, ion_type='n', monoisotopic=monoisotopic)
                        fragment_mz = adjust_mz(fragment_mass, charge=c)

                        if return_type == 'fragment':
                            yield Fragment(charge=c, ion_type=ion_type, start=span[0], end=span[1], monoisotopic=monoisotopic,
                                           isotope=iso, loss=loss, parent_sequence=annotation, mass=fragment_mass,
                                           neutral_mass=fragment_neutral_mass, mz=fragment_mz, sequence=base_sequence,
                                           unmod_sequence=base_unmod_sequence, internal=span[0] != 0 and span[1] != len(annotation))

                        elif return_type == 'mass':
                            yield fragment_mass

                        elif return_type == 'mz':
                            yield fragment_mz


def _get_internal_fragments(sequence: Union[str, ProFormaAnnotation],
                            ion_types: IonTypeType,
                            charges: ChargeType,
                            monoisotopic: bool = True,
                            isotopes: IsotopeType = 0,
                            losses: LossType = 0.0,
                            return_type: str = 'fragment',
                            mass_components: List[float] = None) -> Union[List[float], List[Fragment]]:
    """
    .. code-block:: python

        >>> [frag.sequence for frag in _get_internal_fragments('TIDE', 'by', 1)]
        ['I', 'ID', 'D']

        >>> [round(frag.mz, 3) for frag in _get_internal_fragments('TIDE', 'by', 1)]
        [131.094, 246.121, 133.037]

        >>> [frag.label for frag in _get_internal_fragments('TIDE', 'by', 1)]
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

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    spans = build_non_enzymatic_spans((0, len(annotation), 0))
    internal_spans = [span for span in spans if span[0] != 0 and span[1] != len(annotation)]
    return list(_build_fragments(internal_spans, ion_types, charges, losses, isotopes, monoisotopic, annotation, return_type, mass_components))


def _get_immonium_fragments(sequence: Union[str, ProFormaAnnotation], charges: ChargeType, monoisotopic: bool = True,
                            isotopes: IsotopeType = 0,
                            losses: LossType = 0.0,
                            return_type: str = 'fragment',
                            mass_components: List[float] = None) -> Union[List[float], List[Fragment]]:
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

        >>> [frag.sequence for frag in _get_immonium_fragments('TIDE', 1)]
        ['T', 'I', 'D', 'E']

        >>> [round(frag.mz,3) for frag in _get_immonium_fragments('TIDE', 1)]
        [74.06, 86.096, 88.039, 102.055]

        >>> [frag.label for frag in _get_immonium_fragments('TIDE', 1)]
        ['+i0', '+i1', '+i2', '+i3']

    """
    if isinstance(charges, int):
        charges = [charges]

    if isinstance(isotopes, int):
        isotopes = [isotopes]

    if isinstance(losses, float) or isinstance(losses, int):
        losses = [losses]

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    spans = [(i, i + 1, 0) for i in range(len(annotation))]
    return list(_build_fragments(spans, ['i'], charges, losses, isotopes, monoisotopic, annotation, return_type, mass_components))


def _get_forward_fragments(sequence: Union[str, ProFormaAnnotation],
                           ion_types: List[str],
                           charges: List[int],
                           monoisotopic: bool,
                           isotopes: List[int],
                           losses: List[float],
                           return_type: str = 'fragment',
                           mass_components: List[float] = None) -> Union[List[float], List[Fragment]]:
    start_span = (0, sequence_length(sequence), 0)
    spans = [start_span] + build_left_semi_spans(start_span)
    return list(_build_fragments(spans, ion_types, charges, losses, isotopes, monoisotopic, sequence, return_type, mass_components))


def _get_backward_fragments(sequence: Union[str, ProFormaAnnotation],
                            ion_types: List[str],
                            charges: List[int],
                            monoisotopic: bool,
                            isotopes: List[int],
                            losses: List[float],
                            return_type: str = 'fragment',
                           mass_components: List[float] = None)-> Union[List[float], List[Fragment]]:
    start_span = (0, sequence_length(sequence), 0)
    spans = [start_span] + build_right_semi_spans(start_span)
    return list(_build_fragments(spans, ion_types, charges, losses, isotopes, monoisotopic, sequence, return_type, mass_components))


def _get_terminal_fragments(sequence: Union[str, ProFormaAnnotation],
                            ion_types: IonTypeType,
                            charges: ChargeType,
                            monoisotopic: bool = True,
                            isotopes: IsotopeType = 0,
                            losses: LossType = 0.0,
                            return_type: str = 'fragment',
                            mass_components: List[float] = None) -> Union[List[float], List[Fragment]]:
    """

    .. code-block:: python

        >>> frags = _get_terminal_fragments('TIDE', 'b', 2)
        >>> [round(frag.mass,3) for frag in frags]
        [460.216, 331.173, 216.146, 103.062]

        >>> [round(frag.mz,3) for frag in _get_terminal_fragments("TIDE", ion_types="y", charges=1)]
        [477.219, 376.171, 263.087, 148.06]

        >>> [round(frag.mz,3) for frag in _get_terminal_fragments("TIDE", ion_types="b", charges=2)]
        [230.108, 165.587, 108.073, 51.531]

        >>> [round(frag.mz,3) for frag in _get_terminal_fragments("T[10]IDE", ion_types="y", charges=1)]
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

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    forward_ions = [ion for ion in ion_types if ion in FORWARD_ION_TYPES]
    backward_ions = [ion for ion in ion_types if ion in BACKWARD_ION_TYPES]

    if mass_components is None:
        components = annotation.split()
        mass_components = [mass(sequence=component, charge=0, ion_type='n', monoisotopic=monoisotopic) for component in components]

    frags = []
    frags.extend(_get_forward_fragments(annotation, forward_ions, charges, monoisotopic, isotopes, losses, return_type, mass_components))
    frags.extend(_get_backward_fragments(annotation, backward_ions, charges, monoisotopic, isotopes, losses, return_type, mass_components))

    return frags


def fragment(sequence: Union[str, ProFormaAnnotation],
             ion_types: IonTypeType,
             charges: ChargeType,
             monoisotopic: bool = True,
             isotopes: IsotopeType = 0,
             losses: LossType = 0.0,
             return_type: str = 'fragment') -> Union[List[float], List[Fragment]]:
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
        >>> frags = fragment("[1.0]-P[2.0]E[3.0]-[4.0]", 'y', 1, monoisotopic=True)
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

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    if annotation.contains_sequence_ambiguity():
        raise ValueError("Ambiguous sequence")

    terminal_fragment_types = [i for i in ion_types if i in TERMINAL_ION_TYPES]
    internal_fragment_types = [i for i in ion_types if i in INTERNAL_ION_TYPES]
    immonium = 'i' in ion_types

    components = annotation.split()
    mass_components = [mass(sequence=component, charge=0, ion_type='n', monoisotopic=monoisotopic) for component in
                       components]

    frags = []
    if terminal_fragment_types:
        frags.extend(_get_terminal_fragments(sequence, terminal_fragment_types, charges,
                                             monoisotopic, isotopes, losses, return_type, mass_components))

    if internal_fragment_types:
        frags.extend(_get_internal_fragments(sequence, internal_fragment_types, charges,
                                             monoisotopic, isotopes, losses, return_type, mass_components))

    if immonium:
        frags.extend(_get_immonium_fragments(sequence, charges, monoisotopic, isotopes, losses, return_type, mass_components))

    return frags
