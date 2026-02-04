from collections import Counter
from collections.abc import Mapping
from typing import Any, Literal, cast

from tacular import (
    ELEMENT_LOOKUP,
    FRAGMENT_ION_LOOKUP,
    ElementInfo,
    FragmentIonInfo,
    IonType,
    IonTypeProperty,
)

from ..proforma_components import (
    ChargedFormula,
    GlobalChargeCarrier,
)


class Fragment:
    def __init__(
        self,
        ion_type: IonType | None,
        position: int | tuple[int, int] | str | None,
        mass: float,
        monoisotopic: bool,
        charge_state: int,
        charge_adducts: Mapping[str, int] | None = None,
        isotopes: Mapping[str, int] | int | None = None,
        deltas: Mapping[str | float, int] | None = None,
        composition: Mapping[ElementInfo, int] | None = None,
        sequence: str | None = None,
    ) -> None:
        self.ion_type: IonType | None = ion_type
        self.position: int | tuple[int, int] | str | None = position
        self.mass: int | float = mass
        self.monoisotopic: bool = monoisotopic
        self.charge_state: int = charge_state
        # None means protonated
        self._charge_adducts: Mapping[str, int] | None = charge_adducts
        # int means 13C count
        self._isotopes: Mapping[str, int] | int | None = isotopes
        self._losses: Mapping[str | float, int] | None = deltas
        # Optional composition cache
        self.composition: Counter[ElementInfo] | None = composition
        self.sequence: str | None = sequence

    @property
    def mz(self) -> float:
        return self.mass / self.charge_state if self.charge_state != 0 else self.mass

    @property
    def neutral_mass(self) -> float:
        # subract adduct masses
        total_adduct_mass = 0.0
        for adduct in self.charge_adducts:
            total_adduct_mass += adduct.get_mass(self.monoisotopic) * adduct.occurance
        return self.mass - total_adduct_mass

    @property
    def charge_adducts(self) -> tuple[GlobalChargeCarrier, ...]:
        if self._charge_adducts is not None:
            adducts = []
            for adduct_name, count in self._charge_adducts.items():
                adducts.append(
                    GlobalChargeCarrier(
                        ChargedFormula.from_string(adduct_name, require_formula_prefix=False),
                        count,
                    )
                )
            return tuple(adducts)

        if self.charge_state != 0:
            return (GlobalChargeCarrier.charged_proton(self.charge_state),)

        return ()

    @property
    def is_protonated(self) -> bool:
        if self._charge_adducts is None and self.charge_state != 0:
            return True
        return False

    @property
    def isotopes(self) -> Mapping[ElementInfo, int]:
        if self._isotopes is not None:
            if isinstance(self._isotopes, int):
                return {ELEMENT_LOOKUP["13C"]: self._isotopes}

            if isinstance(self._isotopes, dict):
                return {ELEMENT_LOOKUP[elem]: count for elem, count in self._isotopes.items()}

        return {}

    @property
    def is_c13(self) -> bool:
        if self._isotopes is not None and isinstance(self._isotopes, int):
            return True
        return False

    @property
    def losses(self) -> Mapping[ChargedFormula | float, int]:
        if self._losses is not None:
            losses = {}
            for loss_name, count in self._losses.items():
                if isinstance(loss_name, float | int):
                    losses[loss_name] = count
                    continue
                loss_formula = ChargedFormula.from_string(loss_name, require_formula_prefix=False)
                losses[loss_formula] = count
            return losses
        return {}

    def asdict(self) -> dict[str, Any]:
        return {
            "ion_type": self.ion_type,
            "position": self.position,
            "mass": self.mass,
            "charge_state": self.charge_state,
            "monoisotopic": self.monoisotopic,
            "charge_adducts": self._charge_adducts,
            "isotopes": self._isotopes,
            "losses": self._losses,
        }

    def __str__(self) -> str:
        return (
            f"Fragment(ion_type={self.ion_type}, position={self.position}, "
            f"mass={self.mass}, charge_state={self.charge_state}, "
            f"monoisotopic={self.monoisotopic}, "
            f"charge_adducts={self._charge_adducts}, isotopes={self._isotopes}, "
            f"losses={self._losses})"
        )

    def __repr__(self) -> str:
        return self.__str__()
