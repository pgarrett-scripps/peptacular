from dataclasses import dataclass
from typing import Any, Counter

from tacular import (
    ELEMENT_LOOKUP,
    ElementInfo,
    IonType,
)

from ..constants import ELECTRON_MASS
from ..proforma_components import (
    ChargedFormula,
    GlobalChargeCarrier,
)


@dataclass(frozen=True)
class Fragment:
    ion_type: IonType | None
    position: (
        int | tuple[int, int] | str | None
    )  # int -> terminal, tuple -> internal, str -> immonium, None -> precursor
    mass: float  # not basemass (mass with adjustments + charge/adducts)
    charge_state: int

    # just for record-keeping (already applied to mass)
    monoisotopic: bool
    sequence: str | None = None
    _charge_adducts: tuple[str, ...] | None = None
    _isotopes: dict[str, int] | int | None = None
    _losses: dict[str | float, int] | None = None
    _comp: Counter[ElementInfo] | None = None  # optional composition cache

    @property
    def mz(self) -> float:
        return self.mass / self.charge_state if self.charge_state != 0 else self.mass

    @property
    def composition(self) -> Counter[ElementInfo] | None:
        return self._comp

    @property
    def neutral_mass(self) -> float:
        """Calculate neutral mass (remove charge carriers)"""
        charge_adducts = self.charge_adducts
        adduct_mass = 0.0
        if charge_adducts is not None:
            for adduct in charge_adducts:
                m = adduct.get_mass(monoisotopic=self.monoisotopic)
                if m is None:
                    raise ValueError(f"Mass not available for charge carrier: {adduct}")
                adduct_mass += m
        neutral_mass = self.mass - adduct_mass

        # fix electron mass
        electron_mass_correcttion = self.charge_state * ELECTRON_MASS
        # positive charges had electron mass subtracted, negative charges had electron mass added
        neutral_mass += electron_mass_correcttion

        return neutral_mass

    @property
    def charge_adducts(self) -> tuple[GlobalChargeCarrier, ...] | None:
        if self._charge_adducts is None and self.charge_state != 0:
            if self.charge_state > 0:
                return (GlobalChargeCarrier.charged_proton(self.charge_state),)
            else:
                return None

        if self._charge_adducts is None:
            return None

        adducts: list[GlobalChargeCarrier] = []
        for adduct_str in self._charge_adducts:
            adducts.append(GlobalChargeCarrier.from_string(adduct_str))
        return tuple(adducts)

    @property
    def isotopes(self) -> dict[ElementInfo, int] | None:
        match self._isotopes:
            case dict() as iso_dict:
                result: dict[ElementInfo, int] = {}
                for key, value in iso_dict.items():
                    element = ELEMENT_LOOKUP[key] if isinstance(key, str) else key
                    result[element] = value
                return result
            case int() as iso_count:
                return {ELEMENT_LOOKUP["13C"]: iso_count}
            case None:
                return None
            case _:
                raise TypeError(f"Invalid isotopes type: {type(self._isotopes)}")

    @property
    def losses(self) -> dict[ChargedFormula | float, int] | None:
        if self._losses is None:
            return None

        result: dict[ChargedFormula | float, int] = {}
        for key, value in self._losses.items():
            if isinstance(key, str):
                formula = ChargedFormula.from_string(key, require_formula_prefix=False)
                if formula.charge:
                    raise ValueError(f"Loss formula cannot have charge: {key}")
                result[formula] = value
            elif isinstance(key, float):
                result[key] = value
            else:
                raise TypeError(f"Invalid loss key type: {type(key)}")
        return result

    def __repr__(self) -> str:
        return (
            f"Fragment(ion_type={self.ion_type}, position={self.position}, "
            f"mass={self.mass:.4f}, charge_state={self.charge_state}, "
            f"monoisotopic={self.monoisotopic}, sequence={self.sequence})"
        )

    def asdict(self) -> dict[str, Any]:
        return {
            "ion_type": self.ion_type,
            "position": self.position,
            "mass": self.mass,
            "charge_state": self.charge_state,
            "monoisotopic": self.monoisotopic,
            "sequence": self.sequence,
            "charge_adducts": self._charge_adducts,
            "isotopes": self._isotopes,
            "losses": self._losses,
            "composition": {str(elem): count for elem, count in self._comp.items()}
            if self._comp
            else None,
        }
