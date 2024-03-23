"""
fragment.py contains functions for fragmenting peptides
"""
from dataclasses import dataclass
from functools import cached_property
from typing import List, Generator, Union

from peptacular.proforma.proforma import ProFormaAnnotation
from peptacular.constants import FORWARD_ION_TYPES, BACKWARD_ION_TYPES, INTERNAL_ION_TYPES, TERMINAL_ION_TYPES
from peptacular.mass import mass, adjust_mass, adjust_mz
from peptacular.sequence.sequence import sequence_length, sequence_to_annotation
from peptacular.spans import build_non_enzymatic_spans, build_right_semi_spans, build_left_semi_spans, Span


@dataclass(frozen=True)
class Fragment:
    """
    A dataclass for representing a peptide fragment ion.
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
    def label(self) -> str:
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
                     annotation: ProFormaAnnotation,
                     return_type: str,
                     mass_components: List[float]) -> \
        Union[Generator[Fragment, None, None], Generator[float, None, None]]:
    """
    Builds fragments for a given sequence.
    """

    for span in spans:
        base_mass = adjust_mass(sum(mass_components[span[0]:span[1]]), charge=0, ion_type='n',
                                monoisotopic=monoisotopic)
        base_annotation = annotation.slice(span[0], span[1])
        base_sequence = base_annotation.serialize()
        base_unmod_sequence = base_annotation.sequence
        for ion_type in ion_types:
            base_fragment_mass = adjust_mass(base_mass, charge=0, ion_type=ion_type, monoisotopic=monoisotopic)
            for iso in isotopes:
                for loss in losses:
                    for c in charges:
                        fragment_neutral_mass = adjust_mass(base_mass=base_fragment_mass, charge=0, ion_type='n',
                                                            monoisotopic=monoisotopic, isotope=iso, loss=loss)
                        fragment_mass = adjust_mass(base_mass=fragment_neutral_mass, charge=c, ion_type='n',
                                                    monoisotopic=monoisotopic)
                        fragment_mz = adjust_mz(base_mass=fragment_mass, charge=c)

                        if return_type == 'fragment':
                            yield Fragment(charge=c, ion_type=ion_type, start=span[0], end=span[1],
                                           monoisotopic=monoisotopic,
                                           isotope=iso, loss=loss, parent_sequence=annotation, mass=fragment_mass,
                                           neutral_mass=fragment_neutral_mass, mz=fragment_mz, sequence=base_sequence,
                                           unmod_sequence=base_unmod_sequence,
                                           internal=span[0] != 0 and span[1] != len(annotation))

                        elif return_type == 'mass':
                            yield fragment_mass

                        elif return_type == 'mz':
                            yield fragment_mz


def _get_internal_fragments(annotation: ProFormaAnnotation,
                            ion_types: List[str],
                            charges: List[int],
                            monoisotopic: bool,
                            isotopes: List[int],
                            losses: List[float],
                            return_type: str,
                            mass_components: List[float]) -> Union[List[float], List[Fragment]]:
    """
    Build internal fragments for a given sequence.
    """

    spans = build_non_enzymatic_spans((0, len(annotation), 0))
    internal_spans = [span for span in spans if span[0] != 0 and span[1] != len(annotation)]
    return list(_build_fragments(internal_spans, ion_types, charges, losses, isotopes, monoisotopic, annotation,
                                 return_type, mass_components))


def _get_immonium_fragments(annotation: ProFormaAnnotation,
                            charges: List[int],
                            monoisotopic: bool,
                            isotopes: List[int],
                            losses: List[float],
                            return_type: str,
                            mass_components: List[float]) -> Union[List[float], List[Fragment]]:
    """
    Build immonium ions for a given sequence.
    """
    spans = [(i, i + 1, 0) for i in range(len(annotation))]
    return list(_build_fragments(spans, ['i'], charges, losses, isotopes, monoisotopic, annotation, return_type,
                                 mass_components))


def _get_forward_fragments(annotation: ProFormaAnnotation,
                           ion_types: List[str],
                           charges: List[int],
                           monoisotopic: bool,
                           isotopes: List[int],
                           losses: List[float],
                           return_type: str,
                           mass_components: List[float]) -> Union[List[float], List[Fragment]]:
    """
    Build forward fragments for a given sequence.
    """

    start_span = (0, sequence_length(annotation), 0)
    spans = [start_span] + build_left_semi_spans(start_span)
    return list(_build_fragments(spans, ion_types, charges, losses, isotopes, monoisotopic, annotation, return_type,
                                 mass_components))


def _get_backward_fragments(annotation: ProFormaAnnotation,
                            ion_types: List[str],
                            charges: List[int],
                            monoisotopic: bool,
                            isotopes: List[int],
                            losses: List[float],
                            return_type: str,
                            mass_components: List[float]) -> Union[List[float], List[Fragment]]:
    """
    Build backward fragments for a given sequence.
    """

    start_span = (0, sequence_length(annotation), 0)
    spans = [start_span] + build_right_semi_spans(start_span)
    return list(_build_fragments(spans, ion_types, charges, losses, isotopes, monoisotopic, annotation, return_type,
                                 mass_components))


def _get_terminal_fragments(annotation: ProFormaAnnotation,
                            ion_types: List[str],
                            charges: List[int],
                            monoisotopic: bool,
                            isotopes: List[int],
                            losses: List[float],
                            return_type: str,
                            mass_components: List[float]) -> Union[List[float], List[Fragment]]:
    """
    Build terminal fragments for a given sequence.
    """

    forward_ions = [ion for ion in ion_types if ion in FORWARD_ION_TYPES]
    backward_ions = [ion for ion in ion_types if ion in BACKWARD_ION_TYPES]

    frags = []
    frags.extend(_get_forward_fragments(annotation, forward_ions, charges, monoisotopic, isotopes, losses,
                                        return_type, mass_components))
    frags.extend(_get_backward_fragments(annotation, backward_ions, charges, monoisotopic, isotopes, losses,
                                         return_type, mass_components))

    return frags


def fragment(sequence: Union[str, ProFormaAnnotation],
             ion_types: Union[List[str], str],
             charges: Union[List[int], int],
             monoisotopic: bool = True,
             isotopes: Union[List[int], int] = 0,
             losses: Union[List[float], float] = 0.0,
             return_type: str = 'fragment',
             _mass_components: List[float] = None) -> Union[List[float], List[Fragment]]:
    """
    Builds all Fragment objects or a given input 'sequence'.

    :param sequence: The amino acid sequence or ProForma annotation.
    :type sequence: str | ProFormaAnnotation
    :param ion_types: A list of ion types to consider, e.g., ['b', 'y'], or a single ion type, e.g., 'b'.
    :type ion_types: List[str] | str
    :param charges: A list of charge states for the fragments, or a single charge state.
    :type charges: List[int] | int
    :param monoisotopic: If True, use monoisotopic masses. If False, use average masses, default is [True].
    :type monoisotopic: bool
    :param isotopes: A list of isotope offsets to consider, or a single isotope offset, default is [0].
    :type isotopes: List[int] | int
    :param losses: A list of neutral losses to consider, or a single neutral loss, default is [0.0].
    :type losses: List[float] | float
    :param return_type: The type of data to return, either 'fragment', 'mass', or 'mz', default is 'fragment'.
    :type return_type: str
    :param _mass_components: The split mass components. Used for more efficient fragmenting.
    :type _mass_components: List[float]

    :return: List of Fragment objects or a list of masses or m/z values.
    :rtype: List[Fragment] | List[float]

    .. code-block:: python

        # Create a list of y ions for a modified peptide.
        >>> frags = fragment("[1.0]-P[2.0]E[3.0]-[4.0]", 'y', 1, monoisotopic=True)
        >>> len(frags)
        2

        >>> frags[0].sequence
        '[1.0]-P[2.0]E[3.0]-[4.0]'

        >>> frags[1].sequence
        'E[3.0]-[4.0]'

        >>> frags = fragment('TIDE', ['b'], [2])
        >>> [round(frag.mass,3) for frag in frags]
        [460.216, 331.173, 216.146, 103.062]

        >>> [round(frag.mz,3) for frag in fragment('TIDE', ion_types=["y"], charges=[1])]
        [477.219, 376.171, 263.087, 148.06]

        >>> [round(frag.mz,3) for frag in fragment('TIDE', ion_types=["b"], charges=[2])]
        [230.108, 165.587, 108.073, 51.531]

        >>> [round(frag.mz,3) for frag in fragment('T[10]IDE', ion_types=["y"], charges=[1])]
        [487.219, 376.171, 263.087, 148.06]

        >>> [frag.sequence for frag in fragment('TIDE', ion_types=["i"], charges=[1])]
        ['T', 'I', 'D', 'E']

        >>> [round(frag.mz,3) for frag in fragment('TIDE', ion_types=["i"], charges=[1])]
        [74.06, 86.096, 88.039, 102.055]

        >>> [frag.label for frag in fragment('TIDE', ion_types=["i"], charges=[1])]
        ['+i0', '+i1', '+i2', '+i3']

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

    if _mass_components is None:
        components = annotation.split()
        _mass_components = [mass(sequence=component, charge=0, ion_type='n', monoisotopic=monoisotopic) for component in
                            components]

    frags = []
    if terminal_fragment_types:
        frags.extend(_get_terminal_fragments(annotation, terminal_fragment_types, charges,
                                             monoisotopic, isotopes, losses, return_type, _mass_components))

    if internal_fragment_types:
        frags.extend(_get_internal_fragments(annotation, internal_fragment_types, charges,
                                             monoisotopic, isotopes, losses, return_type, _mass_components))

    if immonium:
        frags.extend(_get_immonium_fragments(annotation, charges, monoisotopic, isotopes, losses, return_type,
                                             _mass_components))

    return frags


class Fragmenter:
    """
    A class for building peptide fragments. Stores annotation and mass components for more efficient fragmenting.
    """

    def __init__(self, sequence: Union[str, ProFormaAnnotation], monoisotopic: bool = True):
        if isinstance(sequence, str):
            self.annotation = sequence_to_annotation(sequence)
        else:
            self.annotation = sequence

        self.monoisotopic = monoisotopic

        self.components = self.annotation.split()
        self.mass_components = [mass(sequence=component, charge=0, ion_type='n', monoisotopic=self.monoisotopic) for
                                component in self.components]

    def fragment(self,
                 ion_types: Union[List[str], str],
                 charges: Union[List[int], int],
                 isotopes: Union[List[int], int] = 0,
                 losses: Union[List[float], float] = 0.0,
                 return_type: str = 'fragment') -> Union[List[float], List[Fragment]]:
        """
        Builds all Fragment objects or a given input 'sequence'.
        """
        return fragment(self.annotation, ion_types, charges, self.monoisotopic, isotopes, losses, return_type,
                        self.mass_components)
