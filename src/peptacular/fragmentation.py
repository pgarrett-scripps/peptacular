"""
fragmentation.py contains functions for fragmenting peptides
"""
import itertools
import re
from dataclasses import dataclass
from functools import cached_property
from typing import List, Union, Literal, Optional, Tuple, Set

from peptacular.proforma.proforma_parser import ProFormaAnnotation
from peptacular.constants import FORWARD_ION_TYPES, BACKWARD_ION_TYPES, INTERNAL_ION_TYPES, TERMINAL_ION_TYPES
from peptacular.mass_calc import mass, adjust_mass, adjust_mz
from peptacular.sequence.sequence_funcs import sequence_length, sequence_to_annotation
from peptacular.spans import build_non_enzymatic_spans, build_right_semi_spans, build_left_semi_spans, Span


def get_number(ion_type: str, len_sequence: int, start: int, end: int) -> str:
    """
    Returns the number of the fragment, e.g., 2 for b2, 3 for y3, etc.

    :return: Number of the fragment.
    :rtype: str
    """

    if ion_type in FORWARD_ION_TYPES:
        number = end
    elif ion_type in BACKWARD_ION_TYPES:
        number = len_sequence - start
    elif ion_type in INTERNAL_ION_TYPES:
        number = f'{start}-{end}'
    elif ion_type == 'i':
        number = start
    else:
        raise ValueError('Wrong Ion Type')

    return number


def get_label(ion_type: str, charge: int, number: str, loss: float, isotope: int) -> str:
    """
    Returns the label of the fragment, e.g., b2, y3i, etc.

    :return: Label of the fragment.
    :rtype: str
    """

    return f"{'+' * charge}" \
           f"{ion_type}" \
           f"{number}" \
           f"{'(' + str(loss) + ')' if loss != 0.0 else ''}" \
           f"{'*' * isotope if isotope > 0 else ''}"


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
    def number(self) -> str:
        """
        Returns the number of the fragment, e.g., 2 for b2, 3 for y3, etc.

        :return: Number of the fragment.
        :rtype: str
        """

        return get_number(self.ion_type, len(self.parent_sequence), self.start, self.end)

    @cached_property
    def label(self) -> str:
        """
        Returns the label of the fragment, e.g., b2, y3i, etc.

        :return: Label of the fragment.
        :rtype: str
        """

        return get_label(self.ion_type, self.charge, self.number, self.loss, self.isotope)

    def __iter__(self):
        # Include regular attributes
        for key, value in self.__dict__.items():
            if key == 'parent_sequence' and isinstance(value, ProFormaAnnotation):
                yield key, value.serialize()
            else:
                yield key, value

        # Explicitly include cached properties
        yield 'label', self.label
        yield 'number', self.number

    def to_dict(self):
        """
        Convert the Fragment object to a dictionary, including cached properties.
        """

        return dict(self)


FragmentReturnType = Literal["fragment", "mass", "mz", "label", "mass-label", "mz-label"]
FRAGMENT_RETURN_TYPING = Union[List[Fragment], List[float], List[str], List[Tuple[float, str]]]


def get_losses(sequence: str, losses: List[Tuple[str, float]], max_losses: int) -> Set[float]:
    """
    Returns a set of applicable losses for a given sequence.

    :param sequence: The sequence to check for losses.
    :type sequence: str
    :param losses: A list of losses to consider.
    :type losses: List[Tuple[str, float]]
    :param max_losses: The maximum number of losses to consider.
    :type max_losses: int

    :return: A set of applicable losses.
    :rtype: Set[float]

    . code-block:: python

        >>> get_losses('AA', [('A', -10), ('A', -5)], 1)
        {0.0, -5, -10}

        >>> get_losses('AA', [('A', -10), ('A', -5)], 2)
        {0.0, -15, -10, -5, -20}

        >>> get_losses('AA', [('A', -10), ('A', -5)], 3)
        {0.0, -15, -10, -25, -5, -20}

    """
    applicable_losses = []
    for restr, loss in losses:
        for _ in re.findall(restr, sequence):
            applicable_losses.append(loss)

    # if max_losses is > 1, find all combinations of losses
    loss_combinations = set()
    if max_losses > 1:
        for loss_count in range(2, max_losses + 1):
            for loss_combination in list(itertools.combinations(applicable_losses, loss_count)):
                loss_combinations.add(sum(loss_combination))

    applicable_losses = set(applicable_losses) | loss_combinations

    if 0.0 not in applicable_losses:
        applicable_losses.add(0.0)

    return applicable_losses


def _build_fragments(spans: List[Span],
                     ion_types: List[str],
                     charges: List[int],
                     losses: List[Tuple[str, float]],
                     isotopes: List[int],
                     monoisotopic: bool,
                     annotation: ProFormaAnnotation,
                     return_type: str,
                     mass_components: List[float],
                     max_losses: int,
                     precision: Optional[int] = None) -> FRAGMENT_RETURN_TYPING:
    """
    Builds fragments for a given sequence.
    """
    frags = []
    for span in spans:
        base_mass = adjust_mass(sum(mass_components[span[0]:span[1]]), charge=0, ion_type='n',
                                monoisotopic=monoisotopic)
        base_annotation = annotation.slice(span[0], span[1])
        base_sequence = base_annotation.serialize()
        base_unmod_sequence = base_annotation.sequence

        applicable_losses = get_losses(sequence=base_unmod_sequence, losses=losses, max_losses=max_losses)

        for ion_type in ion_types:
            # base_fragment_mass = adjust_mass(base_mass, charge=0, ion_type=ion_type, monoisotopic=monoisotopic)
            for iso in isotopes:
                for loss in applicable_losses:
                    for c in charges:
                        fragment_neutral_mass = adjust_mass(base_mass=base_mass, charge=0, ion_type=ion_type,
                                                            monoisotopic=monoisotopic, isotope=iso, loss=loss)
                        fragment_mass = adjust_mass(base_mass=base_mass, charge=c, ion_type=ion_type,
                                                    monoisotopic=monoisotopic, precision=precision,
                                                    isotope=iso, loss=loss)
                        fragment_mz = adjust_mz(base_mass=fragment_mass, charge=c, precision=precision)

                        if return_type == 'fragment':

                            frags.append(Fragment(charge=c, ion_type=ion_type, start=span[0], end=span[1],
                                                  monoisotopic=monoisotopic,
                                                  isotope=iso, loss=loss, parent_sequence=annotation,
                                                  mass=fragment_mass,
                                                  neutral_mass=fragment_neutral_mass, mz=fragment_mz,
                                                  sequence=base_sequence,
                                                  unmod_sequence=base_unmod_sequence,
                                                  internal=span[0] != 0 and span[1] != len(annotation)))

                        elif return_type == 'label':
                            number = get_number(ion_type, len(base_unmod_sequence), span[0], span[1])
                            frags.append(get_label(ion_type, c, number, loss, iso))

                        elif return_type == 'mass':
                            frags.append(fragment_mass)

                        elif return_type == 'mz':
                            frags.append(fragment_mz)

                        elif return_type == 'mass-label':
                            number = get_number(ion_type, len(base_unmod_sequence), span[0], span[1])
                            frags.append((fragment_mass, get_label(ion_type, c, number, loss, iso)))

                        elif return_type == 'mz-label':
                            number = get_number(ion_type, len(base_unmod_sequence), span[0], span[1])
                            frags.append((fragment_mz, get_label(ion_type, c, number, loss, iso)))
    return frags


def _get_internal_fragments(annotation: ProFormaAnnotation,
                            ion_types: List[str],
                            charges: List[int],
                            monoisotopic: bool,
                            isotopes: List[int],
                            losses: List[Tuple[str, float]],
                            return_type: FragmentReturnType,
                            mass_components: List[float],
                            precision: Optional[int],
                            max_losses: int) -> FRAGMENT_RETURN_TYPING:
    """
    Build internal fragments for a given sequence.
    """

    spans = list(build_non_enzymatic_spans((0, len(annotation), 0)))
    internal_spans = [span for span in spans if span[0] != 0 and span[1] != len(annotation)]
    return list(_build_fragments(spans=internal_spans,
                                 ion_types=ion_types,
                                 charges=charges,
                                 losses=losses,
                                 isotopes=isotopes,
                                 monoisotopic=monoisotopic,
                                 annotation=annotation,
                                 return_type=return_type,
                                 mass_components=mass_components,
                                 max_losses=max_losses,
                                 precision=precision))


def _get_immonium_fragments(annotation: ProFormaAnnotation,
                            charges: List[int],
                            monoisotopic: bool,
                            isotopes: List[int],
                            losses: List[Tuple[str, float]],
                            return_type: str,
                            mass_components: List[float],
                            precision: Optional[int],
                            max_losses: int) -> FRAGMENT_RETURN_TYPING:
    """
    Build immonium ions for a given sequence.
    """
    spans = [(i, i + 1, 0) for i in range(len(annotation))]
    return list(_build_fragments(spans=spans,
                                 ion_types=['i'],
                                 charges=charges,
                                 losses=losses,
                                 isotopes=isotopes,
                                 monoisotopic=monoisotopic,
                                 annotation=annotation,
                                 return_type=return_type,
                                 mass_components=mass_components,
                                 max_losses=max_losses,
                                 precision=precision))


def _get_forward_fragments(annotation: ProFormaAnnotation,
                           ion_types: List[str],
                           charges: List[int],
                           monoisotopic: bool,
                           isotopes: List[int],
                           losses: List[Tuple[str, float]],
                           return_type: str,
                           mass_components: List[float],
                           precision: Optional[int],
                           max_losses: int) -> FRAGMENT_RETURN_TYPING:
    """
    Build forward fragments for a given sequence.
    """

    start_span = (0, sequence_length(annotation), 0)
    spans = [start_span] + list(build_left_semi_spans(start_span))
    return list(_build_fragments(spans=spans,
                                 ion_types=ion_types,
                                 charges=charges,
                                 losses=losses,
                                 isotopes=isotopes,
                                 monoisotopic=monoisotopic,
                                 annotation=annotation,
                                 return_type=return_type,
                                 mass_components=mass_components,
                                 max_losses=max_losses,
                                 precision=precision))


def _get_backward_fragments(annotation: ProFormaAnnotation,
                            ion_types: List[str],
                            charges: List[int],
                            monoisotopic: bool,
                            isotopes: List[int],
                            losses: List[Tuple[str, float]],
                            return_type: FragmentReturnType,
                            mass_components: List[float],
                            precision: Optional[int],
                            max_losses: int) -> FRAGMENT_RETURN_TYPING:
    """
    Build backward fragments for a given sequence.
    """

    start_span = (0, sequence_length(annotation), 0)
    spans = [start_span] + list(build_right_semi_spans(start_span))
    return list(_build_fragments(spans=spans,
                                 ion_types=ion_types,
                                 charges=charges,
                                 losses=losses,
                                 isotopes=isotopes,
                                 monoisotopic=monoisotopic,
                                 annotation=annotation,
                                 return_type=return_type,
                                 mass_components=mass_components,
                                 max_losses=max_losses,
                                 precision=precision))


def _get_terminal_fragments(annotation: ProFormaAnnotation,
                            ion_types: List[str],
                            charges: List[int],
                            monoisotopic: bool,
                            isotopes: List[int],
                            losses: List[Tuple[str, float]],
                            return_type: FragmentReturnType,
                            mass_components: List[float],
                            precision: Optional[int],
                            max_losses: int) -> FRAGMENT_RETURN_TYPING:
    """
    Build terminal fragments for a given sequence.
    """

    forward_ions = [ion for ion in ion_types if ion in FORWARD_ION_TYPES]
    backward_ions = [ion for ion in ion_types if ion in BACKWARD_ION_TYPES]

    frags = []
    frags.extend(_get_forward_fragments(annotation=annotation,
                                        ion_types=forward_ions,
                                        charges=charges,
                                        monoisotopic=monoisotopic,
                                        isotopes=isotopes,
                                        losses=losses,
                                        return_type=return_type,
                                        mass_components=mass_components,
                                        precision=precision,
                                        max_losses=max_losses))
    frags.extend(_get_backward_fragments(annotation=annotation,
                                         ion_types=backward_ions,
                                         charges=charges,
                                         monoisotopic=monoisotopic,
                                         isotopes=isotopes,
                                         losses=losses,
                                         return_type=return_type,
                                         mass_components=mass_components,
                                         precision=precision,
                                         max_losses=max_losses))

    return frags


def fragment(sequence: Union[str, ProFormaAnnotation],
             ion_types: Union[List[str], str],
             charges: Union[List[int], int],
             monoisotopic: bool = True,
             isotopes: Union[List[int], int] = 0,
             water_loss: bool = False,
             ammonia_loss: bool = False,
             losses: Optional[Union[List[Tuple[str, float]], Tuple[str, float]]] = None,
             max_losses: int = 1,
             return_type: FragmentReturnType = 'fragment',
             precision: Optional[int] = None,
             _mass_components: Optional[List[float]] = None) -> FRAGMENT_RETURN_TYPING:
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
    :param water_loss: If True, consider water loss, default is False.
    :type water_loss: bool
    :param ammonia_loss: If True, consider ammonia loss, default is False.
    :type ammonia_loss: bool
    :param losses: A list of neutral losses to consider, or a single neutral loss, default is [0.0].
    :type losses: List[float] | float
    :param max_losses: The maximum number of losses to consider, default is 1.
    :type max_losses: int
    :param return_type: The type of data to return, either 'fragment', 'mass', or 'mz', default is 'fragment'.
    :type return_type: str
    :param precision: The number of decimal places to round the masses or m/z values to, default is None.
    :type precision: int
    :param _mass_components: The split mass components. Used for more efficient fragmenting.
    :type _mass_components: List[float]

    :return: List of Fragment objects or a list of masses or m/z values.
    :rtype: List[Fragment] | List[float]

    .. code-block:: python

        # By default a fragment object is returned
        >>> len(fragment("[1.0]-P[2.0]E[3.0]-[4.0]", 'y', 1))
        2
        >>> fragment("[1.0]-P[2.0]E[3.0]-[4.0]", 'y', 1)[0].sequence
        '[1.0]-P[2.0]E[3.0]-[4.0]'
        >>> fragment("[1.0]-P[2.0]E[3.0]-[4.0]", 'y', 1)[1].sequence
        'E[3.0]-[4.0]'

        # Charge 1 B-Ions mass
        >>> fragment(sequence='TIDE', ion_types='b', charges=1, return_type='mass', precision=3)
        [459.209, 330.166, 215.139, 102.055]

        # Charge 1 B-Ions mz
        >>> fragment(sequence='TIDE', ion_types='b', charges=1, return_type='mz', precision=3)
        [459.209, 330.166, 215.139, 102.055]

        # Charge 2 B-Ions mz
        >>> fragment(sequence='TIDE', ion_types='b', charges=2, return_type='mz', precision=3)
        [230.108, 165.587, 108.073, 51.531]

        # Charge 1 Y-Ions mass
        >>> fragment(sequence='T[10]IDE', ion_types='y', charges=1, return_type='mass', precision=3)
        [487.219, 376.171, 263.087, 148.06]

        # Get fragment objects by default
        >>> fragments = fragment(sequence='TIDE', ion_types="i", charges=1, return_type='fragment')
        >>> list(map(lambda frag: frag.sequence, fragments))
        ['T', 'I', 'D', 'E']

        # Immonium ions
        >>> fragment(sequence='TIDE', ion_types="i", charges=1, return_type='mz', precision=3)
        [74.06, 86.096, 88.039, 102.055]

        # Internal fragment ions and 'mz-label' return type
        >>> fragment(sequence='TIDE', ion_types="by", charges=1, return_type='mz-label', precision=3)
        [(114.091, '+by1-2'), (229.118, '+by1-3'), (116.034, '+by2-3')]

        # regex loss losses
        >>> fragment(sequence='TIDE', ion_types="b", charges=2, return_type='mz', precision=3, losses=('E', -10))
        [230.108, 225.108, 165.587, 108.073, 51.531]

        # Multiple losses (when max_losses == 1, only one loss is considered at a time)
        >>> l = [('A', -10), ('A', -5)]
        >>> fragment('AA', "b", 1, return_type='mz', precision=3, losses=l)
        [143.082, 138.082, 133.082, 72.044, 67.044, 62.044]
        >>> fragment('AA', "b", 1, return_type='label', precision=3, losses=l)
        ['+b2', '+b2(-5)', '+b2(-10)', '+b1', '+b1(-5)', '+b1(-10)']

        # Multiple losses (when max_losses > 1, all unique combinations of losses are considered)
        >>> l = [('A', -10), ('A', -5)]
        >>> fragment('AA', "b", 1, return_type='mz', precision=3, losses=l, max_losses=2)
        [143.082, 128.082, 133.082, 138.082, 123.082, 72.044, 57.044, 67.044, 62.044]
        >>> fragment('AA', "b", 1, return_type='label', precision=3, losses=l, max_losses=2)
        ['+b2', '+b2(-15)', '+b2(-10)', '+b2(-5)', '+b2(-20)', '+b1', '+b1(-15)', '+b1(-5)', '+b1(-10)']

        # Water and ammonia losses
        >>> fragment('AQE', "b", 1, return_type='mz', precision=3, water_loss=True, ammonia_loss=True)
        [329.146, 311.135, 312.119, 200.103, 183.076, 72.044]
        >>> fragment('AQE', "b", 1, return_type='label', precision=3, water_loss=True, ammonia_loss=True)
        ['+b3', '+b3(-18.01056)', '+b3(-17.02655)', '+b2', '+b2(-17.02655)', '+b1']

    """

    if not isinstance(ion_types, List):
        ion_types = [ion_types]

    if not isinstance(charges, List):
        charges = [charges]

    if not isinstance(isotopes, List):
        isotopes = [isotopes]

    if losses is None:
        losses = []
    elif not isinstance(losses, List):
        losses = [losses]

    if water_loss:
        losses.append(('[STED]', -18.01056))

    if ammonia_loss:
        losses.append(('[RKNQ]', -17.02655))

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    annotation.pop_labile_mods()

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
                                             monoisotopic, isotopes, losses, return_type,
                                             _mass_components, precision, max_losses))

    if internal_fragment_types:
        frags.extend(_get_internal_fragments(annotation, internal_fragment_types, charges,
                                             monoisotopic, isotopes, losses, return_type,
                                             _mass_components, precision, max_losses))

    if immonium:
        frags.extend(_get_immonium_fragments(annotation, charges, monoisotopic, isotopes, losses, return_type,
                                             _mass_components, precision, max_losses))

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
                 water_loss: bool = False,
                 ammonia_loss: bool = False,
                 losses: Optional[Union[List[Tuple[str, float]], Tuple[str, float]]] = None,
                 max_losses: int = 1,
                 return_type: FragmentReturnType = 'fragment',
                 precision: int = None,
                 ) -> FRAGMENT_RETURN_TYPING:
        """
        Builds all Fragment objects or a given input 'sequence'.
        """
        return fragment(sequence=self.annotation,
                        ion_types=ion_types,
                        charges=charges,
                        monoisotopic=self.monoisotopic,
                        isotopes=isotopes,
                        water_loss=water_loss,
                        ammonia_loss=ammonia_loss,
                        losses=losses,
                        max_losses=max_losses,
                        return_type=return_type,
                        precision=precision,
                        _mass_components=self.mass_components)
