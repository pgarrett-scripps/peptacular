from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache, cached_property
from typing import Self, cast

from tacular import ELEMENT_LOOKUP, NEUTRAL_DELTA_LOOKUP, ElementInfo

from ..proforma_components import ChargedFormula, GlobalChargeCarrier


@dataclass(frozen=True)
class IsotopeInfo:
    data: tuple[tuple[ElementInfo, int], ...]

    @staticmethod
    def from_input(
        isotopes: int | dict[ElementInfo | str, int] | None,
    ) -> Self:
        return cast(Self, get_isotopes(isotopes))

    @cached_property
    def _get_adjustments(self) -> tuple[tuple[ElementInfo, int], ...]:
        data = []
        for iso_elem_info, count in self.data:
            mono_elem_info = ELEMENT_LOOKUP[iso_elem_info.symbol]
            data.append((mono_elem_info, -count))
            data.append((iso_elem_info, count))
        return tuple(data)

    def adjust_composition(self, base_comp: Counter[ElementInfo]) -> None:
        for elem_info, count in self._get_adjustments:
            base_comp[elem_info] += count

        if any(v < 0 for v in base_comp.values()):
            negative_counts = {str(k): v for k, v in base_comp.items() if v < 0}

            raise ValueError(f"Isotopic adjustment resulted in negative element counts: {negative_counts}")

    @property
    def composition(self) -> Counter[ElementInfo]:
        composition: Counter[ElementInfo] = Counter()
        for elem_info, count in self._get_adjustments:
            composition[elem_info] += count
        return composition

    @cached_property
    def monoisotopic_mass_delta(self) -> float:
        delta = 0.0
        for elem_info, count in self._get_adjustments:
            delta += elem_info.get_mass(monoisotopic=True) * count
        return delta

    @cached_property
    def average_mass_delta(self) -> float:
        delta = 0.0
        for elem_info, count in self._get_adjustments:
            delta += elem_info.get_mass(monoisotopic=False) * count
        return delta

    def get_mass_delta(self, monoisotopic: bool = False) -> float:
        if monoisotopic:
            return self.monoisotopic_mass_delta
        return self.average_mass_delta

    def to_fragment_mapping(self) -> Mapping[str, int] | int | None:
        if len(self.data) == 1:
            elem_info, count = self.data[0]
            if str(elem_info) == "13C":
                return count
        return {str(elem_info): count for elem_info, count in self.data}


@dataclass(frozen=True)
class ChargeCarrierInfo:
    adducts: tuple[GlobalChargeCarrier, ...]

    @cached_property
    def charge(self) -> int:
        return sum(adduct.get_charge() for adduct in self.adducts)

    def get_mass(self, monoisotopic: bool = False) -> float:
        if monoisotopic:
            return self.monoisotopic_mass
        return self.average_mass

    @cached_property
    def monoisotopic_mass(self) -> float:
        return sum(adduct.get_mass(monoisotopic=True) for adduct in self.adducts)

    @cached_property
    def average_mass(self) -> float:
        return sum(adduct.get_mass(monoisotopic=False) for adduct in self.adducts)

    @cached_property
    def composition(self) -> Counter[ElementInfo]:
        composition: Counter[ElementInfo] = Counter()
        for adduct in self.adducts:
            composition += adduct.get_composition()
        return composition

    def adjust_composition(self, base_comp: Counter[ElementInfo]) -> None:
        for adduct in self.adducts:
            adduct_comp = adduct.get_composition()
            for elem_info, count in adduct_comp.items():
                base_comp[elem_info] += count

        if any(v < 0 for v in base_comp.values()):
            raise ValueError(f"Charge carrier adjustment resulted in negative element counts: {base_comp}")

    @staticmethod
    def from_input(
        charge: int | str | GlobalChargeCarrier | tuple[GlobalChargeCarrier | str, ...] | None,
    ) -> Self:
        carriers = handle_charge_input(charge)
        return cast(Self, ChargeCarrierInfo(carriers))

    def to_fragment_mapping(self, has_internal_charge: bool = False) -> Mapping[str, int] | None:
        composition: Counter[ElementInfo] = self.composition
        if not composition:
            return None
        if len(composition) == 1:
            elem_info, count = next(iter(composition.items()))
            if str(elem_info) == "H" and not has_internal_charge:
                return None  # None type means infer protonation. Only do this is no internal charge.

        return {str(elem_info): count for elem_info, count in composition.items()}


@dataclass(frozen=True)
class DeltaInfo:
    deltas: Mapping[ChargedFormula | float, int]

    @property
    def has_floats(self) -> bool:
        return any(isinstance(k, float) for k in self.deltas.keys())

    @cached_property
    def monoisotopic_mass_delta(self) -> float:
        delta = 0.0
        for key, count in self.deltas.items():
            if isinstance(key, ChargedFormula):
                delta += key.get_mass(monoisotopic=True) * count
            else:
                delta += key * count
        return delta

    @cached_property
    def average_mass_delta(self) -> float:
        delta = 0.0
        for key, count in self.deltas.items():
            if isinstance(key, ChargedFormula):
                delta += key.get_mass(monoisotopic=False) * count
            else:
                delta += key * count
        return delta

    def get_mass_delta(self, monoisotopic: bool = False) -> float:
        if monoisotopic:
            return self.monoisotopic_mass_delta
        return self.average_mass_delta

    @property
    def composition(self) -> Counter[ElementInfo]:
        if self.has_floats:
            raise ValueError("Cannot get composition when deltas include floats")
        composition: Counter[ElementInfo] = Counter()
        for key, count in self.deltas.items():
            if isinstance(key, ChargedFormula):
                key_comp = key.get_composition()
                for elem_info, elem_count in key_comp.items():
                    composition[elem_info] += elem_count * count
        return composition

    def adjust_composition(self, base_comp: Counter[ElementInfo]) -> None:
        if self.has_floats:
            raise ValueError("Cannot adjust composition when deltas include floats")
        for key, count in self.deltas.items():
            if isinstance(key, ChargedFormula):
                key_comp = key.get_composition()
                for elem_info, elem_count in key_comp.items():
                    base_comp[elem_info] += elem_count * count

        if any(v < 0 for v in base_comp.values()):
            raise ValueError(f"Delta adjustment resulted in negative element counts: {base_comp}")

    @staticmethod
    def from_input(
        deltas: str | ChargedFormula | float | dict[str | ChargedFormula | float, int] | None,
    ) -> Self:
        normalized: dict[ChargedFormula | int | float, int] = _handle_delta_input(deltas)
        for key in normalized.keys():
            if isinstance(key, ChargedFormula) and key.is_charged:
                raise ValueError("Delta formulas must be neutral (charge=0)")

        return cast(Self, DeltaInfo(normalized))

    def __add__(self, other: "DeltaInfo") -> "DeltaInfo":
        """Combine two DeltaInfo objects (add deltas together)."""
        combined: dict[ChargedFormula | float, int] = dict(self.deltas)
        for key, count in other.deltas.items():
            combined[key] = combined.get(key, 0) + count

        # Remove zero counts
        combined = {k: v for k, v in combined.items() if v != 0}

        return DeltaInfo(combined)

    def __sub__(self, other: "DeltaInfo") -> "DeltaInfo":
        """Subtract one DeltaInfo from another."""
        result: dict[ChargedFormula | float, int] = dict(self.deltas)
        for key, count in other.deltas.items():
            result[key] = result.get(key, 0) - count

        # Remove zero counts
        result = {k: v for k, v in result.items() if v != 0}

        return DeltaInfo(result)

    def __neg__(self) -> "DeltaInfo":
        """Negate all deltas (flip signs)."""
        return DeltaInfo({k: -v for k, v in self.deltas.items()})

    def __mul__(self, scalar: int) -> "DeltaInfo":
        """Multiply all delta counts by a scalar."""
        return DeltaInfo({k: v * scalar for k, v in self.deltas.items()})

    def __rmul__(self, scalar: int) -> "DeltaInfo":
        """Allow scalar * DeltaInfo."""
        return self.__mul__(scalar)

    def __str__(self) -> str:
        parts = []
        for key, count in self.deltas.items():
            if isinstance(key, ChargedFormula):
                part = f"Formula:{str(key)}"
            else:
                part = str(key)
            if count != 1:
                part += f" x{count}"
            parts.append(part)
        return " + ".join(parts) if parts else "No Delta"

    def to_fragment_mapping(self) -> Mapping[str | float, int] | None:
        if not self.deltas:
            return None
        result: dict[str | float, int] = {}
        for key, count in self.deltas.items():
            if isinstance(key, ChargedFormula):
                result[key.serialize(include_formula_prefix=False)] = count
            else:
                result[key] = count
        return result


def _handle_delta_input(
    deltas: str | ChargedFormula | float | dict[str | ChargedFormula | float, int] | None,
) -> dict[ChargedFormula | float, int]:
    # Dict - recursively normalize each key
    if deltas is None:
        return {}
    if isinstance(deltas, dict):
        result: dict[ChargedFormula | float, int] = {}
        for key, count in deltas.items():
            # Recurse on each key to normalize it
            normalized = _handle_delta_input(key)  # Returns single-item dict
            for k, v in normalized.items():
                result[k] = result.get(k, 0) + (v * count)
        return result

    # Single items - base cases
    if isinstance(deltas, str):
        try:
            return {ChargedFormula.from_composition(NEUTRAL_DELTA_LOOKUP[deltas].composition): 1}
        except KeyError:
            pass

        return {ChargedFormula.from_string(deltas, require_formula_prefix=False): 1}

    if isinstance(deltas, (ChargedFormula, float, int)):
        return {deltas: 1}

    raise TypeError(f"Invalid delta type: {type(deltas)}")


def handle_charge_input(
    charge: int | str | GlobalChargeCarrier | tuple[GlobalChargeCarrier | str, ...] | None,
) -> tuple[GlobalChargeCarrier, ...]:
    if isinstance(charge, int):
        adduct: GlobalChargeCarrier = get_charge_adducts(charge)
        return (adduct,)
    elif isinstance(charge, str):
        adduct = GlobalChargeCarrier.from_string(charge)
        return (adduct,)
    elif isinstance(charge, GlobalChargeCarrier):
        return (charge,)
    elif isinstance(charge, tuple):
        adducts = []
        for c in charge:
            adducts.extend(handle_charge_input(c))
        return tuple(adducts)
    elif charge is None:
        return ()
    else:
        raise TypeError(f"Invalid charge type: {type(charge)}")


@cache
def get_charge_adducts(charge_state: int) -> GlobalChargeCarrier:
    return GlobalChargeCarrier.charged_proton(charge_state)


C13: ElementInfo = ELEMENT_LOOKUP["13C"]


def get_isotopes(
    isotopes: int | dict[ElementInfo | str, int] | None,
) -> IsotopeInfo:
    if isotopes is None:
        return _get_isotopes(0)

    if isinstance(isotopes, dict):
        elem_infos: list[tuple[ElementInfo, int]] = []
        for elem, count in isotopes.items():
            if isinstance(elem, str):
                elem_info: ElementInfo = ELEMENT_LOOKUP[elem]
            else:
                elem_info: ElementInfo = elem
            elem_infos.append((elem_info, count))

        return _get_isotopes(tuple(elem_infos))

    if isinstance(isotopes, int):
        if isotopes < 0:
            raise ValueError("Isotope count cannot be negative")
        return _get_isotopes(isotopes)

    raise TypeError(f"Invalid isotope type: {type(isotopes)}")


@cache
def _get_isotopes(
    isotopes: int | tuple[tuple[ElementInfo, int], ...],
) -> IsotopeInfo:
    if isinstance(isotopes, tuple):
        return IsotopeInfo(isotopes)

    if isotopes == 0:
        return IsotopeInfo(())

    return IsotopeInfo(((C13, isotopes),))


@cache
def _get_losses(
    losses: tuple[tuple[ElementInfo, int], ...],
) -> Mapping[ChargedFormula, int]:
    composition: Counter[ElementInfo] = Counter()
    for element, count in losses:
        composition[element] += count
    charged_formula = ChargedFormula.from_composition(composition)
    return {charged_formula: 1}


def get_losses(
    losses: Counter[ElementInfo],
) -> Mapping[ChargedFormula, int] | None:
    return _get_losses(tuple(losses.items())) if losses else None
