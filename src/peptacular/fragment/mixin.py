from typing import Generator, Literal, Mapping, Sequence, overload

from ..constants import IonType, IonTypeLiteral
from .core import fragment, get_losses
from .types import (
    FRAGMENT_RETURN_TYPING,
    Fragment,
    FragmentableAnnotation,
    FragmentLiteral,
    FragmentReturnLiteral,
    FragmentReturnType,
    LabelLiteral,
    MassLabelLiteral,
    MassLiteral,
    MzLabelLiteral,
    MzLiteral,
)


class FragmenterMixin(FragmentableAnnotation):
    """Mixin class providing fragmentation functionality."""
    __slots__ = ()

    @overload
    def fragment(
        self: FragmentableAnnotation,
        ion_types: Sequence[IonTypeLiteral | IonType] | IonTypeLiteral | IonType,
        charges: Sequence[int] | int,
        monoisotopic: bool = True,
        *,
        isotopes: Sequence[int] | int = 0,
        water_loss: bool = False,
        ammonia_loss: bool = False,
        losses: None | Mapping[str, Sequence[float]] = None,
        max_losses: int = 1,
        precision: None | int = None,
        _mass_components: None | Sequence[float] = None,
        return_type: FragmentLiteral | Literal[FragmentReturnType.FRAGMENT],
    ) -> Generator[Fragment, None, None]: ...

    @overload
    def fragment(
        self: FragmentableAnnotation,
        ion_types: Sequence[IonTypeLiteral | IonType] | IonTypeLiteral | IonType,
        charges: Sequence[int] | int,
        monoisotopic: bool = True,
        *,
        isotopes: Sequence[int] | int = 0,
        water_loss: bool = False,
        ammonia_loss: bool = False,
        losses: None | Mapping[str, Sequence[float]] = None,
        max_losses: int = 1,
        precision: None | int = None,
        _mass_components: None | Sequence[float] = None,
        return_type: MassLiteral | Literal[FragmentReturnType.MASS],
    ) -> Generator[float, None, None]: ...

    @overload
    def fragment(
        self: FragmentableAnnotation,
        ion_types: Sequence[IonTypeLiteral | IonType] | IonTypeLiteral | IonType,
        charges: Sequence[int] | int,
        monoisotopic: bool = True,
        *,
        isotopes: Sequence[int] | int = 0,
        water_loss: bool = False,
        ammonia_loss: bool = False,
        losses: None | Mapping[str, Sequence[float]] = None,
        max_losses: int = 1,
        precision: None | int = None,
        _mass_components: None | Sequence[float] = None,
        return_type: MzLiteral | Literal[FragmentReturnType.MZ],
    ) -> Generator[float, None, None]: ...

    @overload
    def fragment(
        self: FragmentableAnnotation,
        ion_types: Sequence[IonTypeLiteral | IonType] | IonTypeLiteral | IonType,
        charges: Sequence[int] | int,
        monoisotopic: bool = True,
        *,
        isotopes: Sequence[int] | int = 0,
        water_loss: bool = False,
        ammonia_loss: bool = False,
        losses: None | Mapping[str, Sequence[float]] = None,
        max_losses: int = 1,
        precision: None | int = None,
        _mass_components: None | Sequence[float] = None,
        return_type: LabelLiteral | Literal[FragmentReturnType.LABEL],
    ) -> Generator[str, None, None]: ...

    @overload
    def fragment(
        self: FragmentableAnnotation,
        ion_types: Sequence[IonTypeLiteral | IonType] | IonTypeLiteral | IonType,
        charges: Sequence[int] | int,
        monoisotopic: bool = True,
        *,
        isotopes: Sequence[int] | int = 0,
        water_loss: bool = False,
        ammonia_loss: bool = False,
        losses: None | Mapping[str, Sequence[float]] = None,
        max_losses: int = 1,
        precision: None | int = None,
        _mass_components: None | Sequence[float] = None,
        return_type: MassLabelLiteral | Literal[FragmentReturnType.MASS_LABEL],
    ) -> Generator[tuple[float, str], None, None]: ...

    @overload
    def fragment(
        self: FragmentableAnnotation,
        ion_types: Sequence[IonTypeLiteral | IonType] | IonTypeLiteral | IonType,
        charges: Sequence[int] | int,
        monoisotopic: bool = True,
        *,
        isotopes: Sequence[int] | int = 0,
        water_loss: bool = False,
        ammonia_loss: bool = False,
        losses: None | Mapping[str, Sequence[float]] = None,
        max_losses: int = 1,
        precision: None | int = None,
        _mass_components: None | Sequence[float] = None,
        return_type: MzLabelLiteral | Literal[FragmentReturnType.MZ_LABEL],
    ) -> Generator[tuple[float, str], None, None]: ...

    @overload
    def fragment(
        self: FragmentableAnnotation,
        ion_types: Sequence[IonTypeLiteral | IonType] | IonTypeLiteral | IonType,
        charges: Sequence[int] | int,
        monoisotopic: bool = True,
        *,
        isotopes: Sequence[int] | int = 0,
        water_loss: bool = False,
        ammonia_loss: bool = False,
        losses: None | Mapping[str, Sequence[float]] = None,
        max_losses: int = 1,
        precision: None | int = None,
        _mass_components: None | Sequence[float] = None,
        return_type: (
            FragmentReturnLiteral | FragmentReturnType
        ) = FragmentReturnType.FRAGMENT,
    ) -> Generator[FRAGMENT_RETURN_TYPING, None, None]: ...

    def fragment(
        self: FragmentableAnnotation,
        ion_types: Sequence[IonTypeLiteral | IonType] | IonTypeLiteral | IonType,
        charges: Sequence[int] | int,
        monoisotopic: bool = True,
        *,
        isotopes: Sequence[int] | int = 0,
        water_loss: bool = False,
        ammonia_loss: bool = False,
        losses: None | Mapping[str, Sequence[float]] = None,
        max_losses: int = 1,
        precision: None | int = None,
        _mass_components: None | Sequence[float] = None,
        return_type: (
            FragmentReturnLiteral | FragmentReturnType
        ) = FragmentReturnType.FRAGMENT,
    ) -> Generator[FRAGMENT_RETURN_TYPING, None, None]:
        return fragment(
            annotation=self,
            ion_types=ion_types,
            charges=charges,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            water_loss=water_loss,
            ammonia_loss=ammonia_loss,
            losses=losses,
            max_losses=max_losses,
            precision=precision,
            _mass_components=_mass_components,
            return_type=return_type,
        )

    def get_losses(
        self: FragmentableAnnotation,
        losses: Mapping[str, Sequence[float]],
        max_losses: int,
    ) -> set[float]:
        return get_losses(annotation=self, losses=losses, max_losses=max_losses)
