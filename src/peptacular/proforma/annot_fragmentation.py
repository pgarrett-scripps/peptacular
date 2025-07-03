from dataclasses import dataclass
from functools import cached_property
import itertools
from typing import List, Literal, Optional, Set, Tuple, Union

import regex as re

from ..constants import (
    BACKWARD_ION_TYPES,
    FORWARD_ION_TYPES,
    INTERNAL_ION_TYPES,
    TERMINAL_ION_TYPES,
)
from ..spans import (
    build_left_semi_spans,
    build_non_enzymatic_spans,
    build_right_semi_spans,
)

from ..utils2 import get_label, get_number
from ..mass_calc import adjust_mass, adjust_mz
from ..proforma_dataclasses import Span
from .annot_manipulation import ProFormaAnnotationManipulation


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
    parent_sequence: 'ProFormaAnnotationFragmentation'
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

        return get_number(
            self.ion_type, len(self.parent_sequence), self.start, self.end
        )

    @cached_property
    def label(self) -> str:
        """
        Returns the label of the fragment, e.g., b2, y3i, etc.

        :return: Label of the fragment.
        :rtype: str
        """

        return get_label(
            self.ion_type, self.charge, self.number, self.loss, self.isotope
        )

    def __iter__(self):
        # Include regular attributes
        for key, value in self.__dict__.items():
            yield key, value

        # Explicitly include cached properties
        yield "label", self.label
        yield "number", self.number

    def to_dict(self):
        """
        Convert the Fragment object to a dictionary, including cached properties.
        """

        return dict(self)
    


FragmentReturnType = Literal[
    "fragment", "mass", "mz", "label", "mass-label", "mz-label"
]
FRAGMENT_RETURN_TYPING = Union[
    List[Fragment], List[float], List[str], List[Tuple[float, str]]
]



class ProFormaAnnotationFragmentation(ProFormaAnnotationManipulation):
    """
    Fragmentation Methods
    """

    def get_losses(
        self, losses: List[Tuple[str, float]], max_losses: int
    ) -> Set[float]:
        applicable_losses = []
        for restr, loss in losses:
            for _ in re.findall(restr, self.stripped_sequence):
                applicable_losses.append(loss)

        # if max_losses is > 1, find all combinations of losses
        loss_combinations = set()
        if max_losses > 1:
            for loss_count in range(2, max_losses + 1):
                for loss_combination in list(
                    itertools.combinations(applicable_losses, loss_count)
                ):
                    loss_combinations.add(sum(loss_combination))

        applicable_losses = set(applicable_losses) | loss_combinations

        if 0.0 not in applicable_losses:
            applicable_losses.add(0.0)

        return applicable_losses

    def _build_fragments(
        self,
        spans: List[Span],
        ion_types: List[str],
        charges: List[int],
        losses: List[Tuple[str, float]],
        isotopes: List[int],
        monoisotopic: bool,
        return_type: str,
        mass_components: List[float],
        max_losses: int,
        precision: Optional[int] = None,
    ) -> FRAGMENT_RETURN_TYPING:
        """
        Builds fragments for a given sequence.
        """
        frags = []
        for span in spans:
            base_mass = adjust_mass(
                sum(mass_components[span[0] : span[1]]),
                charge=0,
                ion_type="n",
                monoisotopic=monoisotopic,
            )
            base_annotation = self.slice(span[0], span[1], inplace=False)
            base_sequence = base_annotation.serialize()
            base_unmod_sequence = base_annotation.sequence

            applicable_losses = base_annotation.get_losses(
                losses=losses, max_losses=max_losses
            )

            for ion_type in ion_types:
                # base_fragment_mass = adjust_mass(base_mass, charge=0, ion_type=ion_type, monoisotopic=monoisotopic)
                for iso in isotopes:
                    for loss in applicable_losses:
                        for c in charges:
                            fragment_neutral_mass = adjust_mass(
                                base_mass=base_mass,
                                charge=0,
                                ion_type=ion_type,
                                monoisotopic=monoisotopic,
                                isotope=iso,
                                loss=loss,
                            )
                            fragment_mass = adjust_mass(
                                base_mass=base_mass,
                                charge=c,
                                ion_type=ion_type,
                                monoisotopic=monoisotopic,
                                precision=precision,
                                isotope=iso,
                                loss=loss,
                            )
                            fragment_mz = adjust_mz(
                                base_mass=fragment_mass, charge=c, precision=precision
                            )

                            if return_type == "fragment":

                                frags.append(
                                    Fragment(
                                        charge=c,
                                        ion_type=ion_type,
                                        start=span[0],
                                        end=span[1],
                                        monoisotopic=monoisotopic,
                                        isotope=iso,
                                        loss=loss,
                                        parent_sequence=self,
                                        mass=fragment_mass,
                                        neutral_mass=fragment_neutral_mass,
                                        mz=fragment_mz,
                                        sequence=base_sequence,
                                        unmod_sequence=base_unmod_sequence,
                                        internal=span[0] != 0 and span[1] != len(self),
                                    )
                                )

                            elif return_type == "label":
                                number = get_number(
                                    ion_type, len(base_unmod_sequence), span[0], span[1]
                                )
                                frags.append(get_label(ion_type, c, number, loss, iso))

                            elif return_type == "mass":
                                frags.append(fragment_mass)

                            elif return_type == "mz":
                                frags.append(fragment_mz)

                            elif return_type == "mass-label":
                                number = get_number(
                                    ion_type, len(base_unmod_sequence), span[0], span[1]
                                )
                                frags.append(
                                    (
                                        fragment_mass,
                                        get_label(ion_type, c, number, loss, iso),
                                    )
                                )

                            elif return_type == "mz-label":
                                number = get_number(
                                    ion_type, len(base_unmod_sequence), span[0], span[1]
                                )
                                frags.append(
                                    (
                                        fragment_mz,
                                        get_label(ion_type, c, number, loss, iso),
                                    )
                                )
        return frags

    def _get_internal_fragments(
        self,
        ion_types: List[str],
        charges: List[int],
        monoisotopic: bool,
        isotopes: List[int],
        losses: List[Tuple[str, float]],
        return_type: FragmentReturnType,
        mass_components: List[float],
        precision: Optional[int],
        max_losses: int,
    ) -> FRAGMENT_RETURN_TYPING:
        """
        Build internal fragments for a given sequence.
        """

        spans = list(build_non_enzymatic_spans((0, len(self), 0)))
        internal_spans = [
            span for span in spans if span[0] != 0 and span[1] != len(self)
        ]
        return list(
            self._build_fragments(
                spans=internal_spans,
                ion_types=ion_types,
                charges=charges,
                losses=losses,
                isotopes=isotopes,
                monoisotopic=monoisotopic,
                return_type=return_type,
                mass_components=mass_components,
                max_losses=max_losses,
                precision=precision,
            )
        )

    def _get_immonium_fragments(
        self,
        charges: List[int],
        monoisotopic: bool,
        isotopes: List[int],
        losses: List[Tuple[str, float]],
        return_type: str,
        mass_components: List[float],
        precision: Optional[int],
        max_losses: int,
    ) -> FRAGMENT_RETURN_TYPING:
        """
        Build immonium ions for a given sequence.
        """
        spans = [(i, i + 1, 0) for i in range(len(self))]
        return list(
            self._build_fragments(
                spans=spans,
                ion_types=["i"],
                charges=charges,
                losses=losses,
                isotopes=isotopes,
                monoisotopic=monoisotopic,
                return_type=return_type,
                mass_components=mass_components,
                max_losses=max_losses,
                precision=precision,
            )
        )

    def _get_forward_fragments(
        self,
        ion_types: List[str],
        charges: List[int],
        monoisotopic: bool,
        isotopes: List[int],
        losses: List[Tuple[str, float]],
        return_type: str,
        mass_components: List[float],
        precision: Optional[int],
        max_losses: int,
    ) -> FRAGMENT_RETURN_TYPING:
        """
        Build forward fragments for a given sequence.
        """

        start_span = (0, len(self), 0)
        spans = [start_span] + list(build_left_semi_spans(start_span))
        return list(
            self._build_fragments(
                spans=spans,
                ion_types=ion_types,
                charges=charges,
                losses=losses,
                isotopes=isotopes,
                monoisotopic=monoisotopic,
                return_type=return_type,
                mass_components=mass_components,
                max_losses=max_losses,
                precision=precision,
            )
        )

    def _get_backward_fragments(
        self,
        ion_types: List[str],
        charges: List[int],
        monoisotopic: bool,
        isotopes: List[int],
        losses: List[Tuple[str, float]],
        return_type: FragmentReturnType,
        mass_components: List[float],
        precision: Optional[int],
        max_losses: int,
    ) -> FRAGMENT_RETURN_TYPING:
        """
        Build backward fragments for a given sequence.
        """

        start_span = (0, len(self), 0)
        spans = [start_span] + list(build_right_semi_spans(start_span))
        return list(
            self._build_fragments(
                spans=spans,
                ion_types=ion_types,
                charges=charges,
                losses=losses,
                isotopes=isotopes,
                monoisotopic=monoisotopic,
                return_type=return_type,
                mass_components=mass_components,
                max_losses=max_losses,
                precision=precision,
            )
        )

    def _get_terminal_fragments(
        self,
        ion_types: List[str],
        charges: List[int],
        monoisotopic: bool,
        isotopes: List[int],
        losses: List[Tuple[str, float]],
        return_type: FragmentReturnType,
        mass_components: List[float],
        precision: Optional[int],
        max_losses: int,
    ) -> FRAGMENT_RETURN_TYPING:
        """
        Build terminal fragments for a given sequence.
        """

        forward_ions = [ion for ion in ion_types if ion in FORWARD_ION_TYPES]
        backward_ions = [ion for ion in ion_types if ion in BACKWARD_ION_TYPES]

        frags = []
        frags.extend(
            self._get_forward_fragments(
                ion_types=forward_ions,
                charges=charges,
                monoisotopic=monoisotopic,
                isotopes=isotopes,
                losses=losses,
                return_type=return_type,
                mass_components=mass_components,
                precision=precision,
                max_losses=max_losses,
            )
        )
        frags.extend(
            self._get_backward_fragments(
                ion_types=backward_ions,
                charges=charges,
                monoisotopic=monoisotopic,
                isotopes=isotopes,
                losses=losses,
                return_type=return_type,
                mass_components=mass_components,
                precision=precision,
                max_losses=max_losses,
            )
        )

        return frags

    def fragment(
        self,
        ion_types: Union[List[str], str],
        charges: Union[List[int], int],
        monoisotopic: bool = True,
        isotopes: Union[List[int], int] = 0,
        water_loss: bool = False,
        ammonia_loss: bool = False,
        losses: Optional[Union[List[Tuple[str, float]], Tuple[str, float]]] = None,
        max_losses: int = 1,
        return_type: FragmentReturnType = "fragment",
        precision: Optional[int] = None,
        _mass_components: Optional[List[float]] = None,
    ) -> FRAGMENT_RETURN_TYPING:

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
            losses.append(("[STED]", -18.01056))

        if ammonia_loss:
            losses.append(("[RKNQ]", -17.02655))

        if self.contains_sequence_ambiguity():
            raise ValueError("Ambiguous sequence")

        terminal_fragment_types = [i for i in ion_types if i in TERMINAL_ION_TYPES]
        internal_fragment_types = [i for i in ion_types if i in INTERNAL_ION_TYPES]
        immonium = "i" in ion_types

        if _mass_components is None:
            components = self.split()
            _mass_components = [
                c.add_charge(0).mass(ion_type="n", monoisotopic=monoisotopic)
                for c in components
            ]

        frags = []
        if terminal_fragment_types:
            frags.extend(
                self._get_terminal_fragments(
                    terminal_fragment_types,
                    charges,
                    monoisotopic,
                    isotopes,
                    losses,
                    return_type,
                    _mass_components,
                    precision,
                    max_losses,
                )
            )

        if internal_fragment_types:
            frags.extend(
                self._get_internal_fragments(
                    internal_fragment_types,
                    charges,
                    monoisotopic,
                    isotopes,
                    losses,
                    return_type,
                    _mass_components,
                    precision,
                    max_losses,
                )
            )

        if immonium:
            frags.extend(
                self._get_immonium_fragments(
                    charges,
                    monoisotopic,
                    isotopes,
                    losses,
                    return_type,
                    _mass_components,
                    precision,
                    max_losses,
                )
            )

        return frags
