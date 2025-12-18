from collections import Counter
from dataclasses import dataclass
from typing import Any, Literal, Sequence, TypeVar, overload

from ..annotation.mod import Mods

from ..fragment.neutral_deltas.data import NeutralDelta, NeutralDeltaLiteral

from ..fragment.ion_types.dclass import FragmentIonInfo, IonTypeProperty
from ..fragment.mzpaf import (
    INTERNAL_MASS_DIFFS,
    FragmentAnnotation,
    ImmoniumIon,
    InternalFragment,
    IsotopeSpecification,
    MassError,
    PeptideIon,
    PrecursorIon,
)

from ..fragment.neutral_deltas.lookup import NEUTRAL_DELTA_LOOKUP

from ..fragment import FRAGMENT_ION_LOOKUP
from ..proforma_components.comps import ChargedFormula, GlobalChargeCarrier
from ..fragment import IonType, IonTypeLiteral
from ..constants import C13_NEUTRON_MASS, ELECTRON_MASS, PROTON_MASS
from ..elements.dclass import ElementInfo
from ..elements.lookup import ELEMENT_LOOKUP

H_ELEMENT_INFO = ELEMENT_LOOKUP["H"]

"""
charge can be int | dict[str, int] | tuple[GlobalChargeCarrier, ...]
"""


@dataclass(frozen=True)
class Fragment:
    ion_type: IonType
    position: (
        int | tuple[int, int] | str | None
    )  # int -> terminal, tuple -> internal, str -> immonium, None -> precursor
    mass: float  # not basemass (mass with adjustments + charge/adducts)
    charge_state: int
    monoisotopic: bool  # just for record-keeping (already applied to mass)
    sequence: str | None = None  # just for record-keeping
    _charge_adducts: tuple[str, ...] | None = (
        None  # just for record-keeping (already applied to mass)
    )
    _isotopes: dict[str, int] | int | None = (
        None  # just for record-keeping (already applied to mass)
    )
    _losses: dict[str | float, int] | None = (
        None  # just for record-keeping (already applied to mass)
    )
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
                if m is None:  # type: ignore
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

    def to_mzpaf(
        self,
        include_annotation: bool = True,
        confidence: float | None = None,
        mass_error: float | None = None,
        mass_error_type: Literal["ppm", "da"] = "ppm",
    ) -> FragmentAnnotation:
        """Convert fragment to mzPAF format string."""

        ion_info: FragmentIonInfo = FRAGMENT_ION_LOOKUP[self.ion_type]

        ion = None
        internal_loss = None
        match ion_info.properties:
            case IonTypeProperty.FORWARD | IonTypeProperty.BACKWARD:
                if not isinstance(self.position, int):
                    position = -1
                else:
                    position = self.position

                # some peptacular ion_types mapt to different mzPAF ion types
                match ion_info.ion_type:
                    case IonType.W_VALINE:
                        ion_type = IonType.W
                    case IonType.D_VALINE:
                        ion_type = IonType.D
                    case IonType.WB_ISOLEUCINE | IonType.WB_THREONINE:
                        ion_type = IonType.WB
                    case IonType.WA_ISOLEUCINE | IonType.WA_THREONINE:
                        ion_type = IonType.WA
                    case IonType.DB_ISOLEUCINE | IonType.DB_THREONINE:
                        ion_type = IonType.DB
                    case IonType.DA_ISOLEUCINE | IonType.DA_THREONINE:
                        ion_type = IonType.DA
                    case _:
                        ion_type = ion_info.ion_type

                ion = PeptideIon(
                    series=ion_type,
                    position=position,
                    sequence=self.sequence if include_annotation else None,
                )
            case IonTypeProperty.INTERNAL:
                if ion_info.id == IonType.IMMONIUM:
                    if not isinstance(self.position, str):
                        position = "?"
                    else:
                        position = self.position

                    ion = ImmoniumIon(amino_acid=position)

                else:
                    if not isinstance(self.position, tuple) or len(self.position) != 2:
                        start = -1
                        end = -1
                    else:
                        start, end = self.position

                    internal_ion_key = tuple(list(ion_info.ion_type.value))
                    if internal_ion_key not in INTERNAL_MASS_DIFFS:
                        raise ValueError(
                            f"Internal ion type {ion_info.ion_type} not supported in mzPAF."
                        )
                    internal_loss = INTERNAL_MASS_DIFFS[internal_ion_key]  # type: ignore

                    ion = InternalFragment(
                        start_position=start,
                        end_position=end,
                        sequence=self.sequence if include_annotation else None,
                    )
            case IonTypeProperty.INTACT:
                if ion_info.ion_type == IonType.PRECURSOR:
                    ion = PrecursorIon()
                else:
                    raise ValueError(
                        f"Cannot convert intact ion type {ion_info.id} to mzPAF."
                    )
            case _:
                pass

        if ion is None:
            raise ValueError(
                f"Cannot convert fragment with ion type {ion_info.id} to mzPAF."
            )

        # handle isotope
        iso = None
        match self._isotopes:
            case dict():
                iso: list[IsotopeSpecification] = []
                for element, count in self.isotopes.items():
                    iso.append(IsotopeSpecification(element=str(element), count=count))
            case int() as iso_count:
                iso: list[IsotopeSpecification] = [
                    IsotopeSpecification(count=iso_count)
                ]
            case None:
                pass
            case _:
                raise TypeError(f"Invalid isotopes type: {type(self.isotopes)}")

        # handle losses
        loss = None
        match self.losses:
            case dict() as loss_tup:
                loss: list[str] = []
                for charged_formula in loss_tup:
                    if isinstance(charged_formula, float):
                        # add str with +/- sign
                        loss.append(f"{charged_formula:+.6f}")
                    elif isinstance(charged_formula, ChargedFormula):
                        loss.append(str(charged_formula))
                    else:
                        raise TypeError(
                            f"Invalid loss type in tuple: {type(charged_formula)}"
                        )
            case None:
                pass
            case _:
                raise TypeError(f"Invalid losses type: {type(self.losses)}")

        if internal_loss is not None:
            if loss is None:
                loss = []
            loss.append(internal_loss)

        if self._charge_adducts is None:
            charge_state = self.charge_state
            adducts = None
        else:
            charge_state = self.charge_state
            adducts = [a.to_mz_paf() for a in self.charge_adducts]

        mz_paf_mass_error = None
        if mass_error is not None:
            mz_paf_mass_error = MassError(value=mass_error, unit=mass_error_type)

        return FragmentAnnotation(
            ion_type=ion,
            neutral_losses=loss,
            isotope=iso,
            adducts=adducts,
            charge=charge_state,
            mass_error=mz_paf_mass_error,
            confidence=confidence,
        )

    def __str__(self) -> str:
        return (
            self.to_mzpaf().__str__() + f" | Mass: {self.mass:.4f} | m/z: {self.mz:.4f}"
        )


def _handle_charge_input_mass(
    charge: int
    | str
    | GlobalChargeCarrier
    | tuple[GlobalChargeCarrier | str, ...]
    | Mods[GlobalChargeCarrier],
    monoisotopic: bool,
) -> tuple[float, int]:
    if isinstance(charge, Mods):
        mass = charge.get_mass(monoisotopic=monoisotopic)
        charge_state = charge.get_charge()

        if charge_state is None:
            raise ValueError("Charge is not defined for this modification.")
        return mass, charge_state

    base_mass = 0.0
    if isinstance(charge, int):
        charge_value = charge
        base_mass += charge_value * PROTON_MASS
    elif isinstance(charge, str) or isinstance(charge, GlobalChargeCarrier):
        if isinstance(charge, str):
            charge = GlobalChargeCarrier.from_string(charge)
        m = charge.get_mass(monoisotopic=monoisotopic)
        base_mass += m
        charge_value = charge.get_charge()
    elif isinstance(charge, tuple):  # type: ignore
        charge_value = 0
        for global_charge_carrier in charge:
            if isinstance(global_charge_carrier, str):
                global_charge_carrier = GlobalChargeCarrier.from_string(
                    global_charge_carrier
                )
            m = global_charge_carrier.get_mass(monoisotopic=monoisotopic)
            base_mass += m
            charge_value += global_charge_carrier.get_charge()
    else:
        raise TypeError(f"Invalid charge type: {type(charge)}")

    return base_mass, charge_value


def handle_charge_input_comp(
    charge: int | str | GlobalChargeCarrier | tuple[GlobalChargeCarrier | str, ...],
) -> tuple[Counter[ElementInfo], int]:
    base_comp = Counter[ElementInfo]()
    total_charge = 0
    if isinstance(charge, int):
        charge_value = charge
        base_comp[H_ELEMENT_INFO] += charge_value
        total_charge = charge_value
    elif isinstance(charge, str) or isinstance(charge, GlobalChargeCarrier):
        if isinstance(charge, str):
            charge = GlobalChargeCarrier.from_string(charge)
        comp = charge.get_composition()
        for element, elem_count in comp.items():
            base_comp[element] += elem_count
        total_charge = charge.get_charge()
    elif isinstance(charge, tuple):  # type: ignore
        for global_charge_carrier in charge:
            if isinstance(global_charge_carrier, str):
                global_charge_carrier = GlobalChargeCarrier.from_string(
                    global_charge_carrier
                )
            comp = global_charge_carrier.get_composition()
            for element, elem_count in comp.items():
                base_comp[element] += elem_count

            total_charge += global_charge_carrier.get_charge()

    return base_comp, total_charge


def adjust_mass_mz(
    base: float | Counter[ElementInfo],
    charge: int
    | str
    | GlobalChargeCarrier
    | tuple[GlobalChargeCarrier | str, ...]
    | None = None,
    ion_type: IonType | IonTypeLiteral | FragmentIonInfo = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotopes: int | dict[str | ElementInfo, int] | None = None,
    neutral_deltas: dict[ChargedFormula, int] | None = None,
    position: int | tuple[int, int] | str | None = None,
    parent_sequence: str | None = None,
) -> Fragment:
    """Adjust base mass by charge carriers and ion type."""

    base_mass = 0.0
    if isinstance(base, Counter):
        for elem, count in base.items():
            base_mass += elem.get_mass(monoisotopic=monoisotopic) * count
    else:
        base_mass = base

    fragment_isotope_notation: dict[str, int] | int | None = None
    if isotopes is not None:
        if isinstance(isotopes, int):
            base_mass += C13_NEUTRON_MASS * isotopes
            fragment_isotope_notation = isotopes
        else:
            for element, count in process_isotopes(isotopes).items():
                mono_element = ELEMENT_LOOKUP[element.symbol]
                mass_diff = element.get_mass(monoisotopic) - mono_element.get_mass(
                    monoisotopic
                )
                base_mass += mass_diff * count

                elem = str(element)
                if fragment_isotope_notation is None:
                    fragment_isotope_notation = {}
                if elem not in fragment_isotope_notation:
                    fragment_isotope_notation[elem] = 0
                fragment_isotope_notation[elem] += count

    # Apply losses
    fragment_loss_notation: dict[str | float, int] | None = None
    if neutral_deltas is not None:
        loss_comp, float_losses = process_losses(neutral_deltas)
        for element, count in loss_comp.items():
            base_mass -= element.get_mass(monoisotopic) * count
        fragment_loss_notation = {}
        for element, count in loss_comp.items():
            elem_str = str(element)
            if elem_str not in fragment_loss_notation:
                fragment_loss_notation[elem_str] = 0
            fragment_loss_notation[elem_str] += count
        for float_loss, count in float_losses.items():
            if fragment_loss_notation is None:
                fragment_loss_notation = {}
            if float_loss not in fragment_loss_notation:
                fragment_loss_notation[float_loss] = 0
            fragment_loss_notation[float_loss] += count

    charge = charge or 0
    total_charge: int = 0
    if charge:
        mass_adjustment, charge_state = _handle_charge_input_mass(
            charge=charge, monoisotopic=monoisotopic
        )
        base_mass += mass_adjustment
        total_charge += charge_state

    fragment_annotation_charge_adducts: tuple[str, ...] | None = None
    if charge:
        match charge:
            case int():
                fragment_annotation_charge_adducts = None
            case GlobalChargeCarrier() | str():
                fragment_annotation_charge_adducts = (str(charge),)
            case tuple():
                tmp_fragment_annotation_charge_adducts: list[str] = []
                for c in charge:
                    if isinstance(c, GlobalChargeCarrier):
                        c = str(c)
                    tmp_fragment_annotation_charge_adducts.append(c)
                fragment_annotation_charge_adducts = tuple(
                    tmp_fragment_annotation_charge_adducts
                )
            case Mods():
                fragment_annotation_charge_adducts = tuple(
                    str(mod) for mod in charge._mods
                )  # type: ignore
            case _:
                raise TypeError(f"Invalid charge type: {type(charge)}")

    ion_info = (
        FRAGMENT_ION_LOOKUP[ion_type]
        if not isinstance(ion_type, FragmentIonInfo)
        else ion_type
    )
    base_mass += ion_info.get_mass(monoisotopic=monoisotopic)
    return Fragment(
        ion_type=ion_info.ion_type,
        position=position,
        mass=base_mass,
        charge_state=total_charge,
        monoisotopic=monoisotopic,
        sequence=parent_sequence,
        _charge_adducts=fragment_annotation_charge_adducts,
        _isotopes=fragment_isotope_notation,
        _losses=fragment_loss_notation,
    )


def adjust_comp(
    base_comp: Counter[ElementInfo],
    charge: int
    | GlobalChargeCarrier
    | str
    | tuple[GlobalChargeCarrier | str, ...]
    | None = None,
    ion_type: IonType | IonTypeLiteral | FragmentIonInfo = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotopes: int | dict[str | ElementInfo, int] | None = None,
    neutral_deltas: dict[ChargedFormula, int] | None = None,
    inplace: bool = True,
    isotope_map: dict[ElementInfo, ElementInfo] | None = None,
    position: int | tuple[int, int] | str | None = None,
    parent_sequence: str | None = None,
) -> Fragment:
    """Adjust base composition by charge carriers and ion type, returning a Fragment object."""

    if not inplace:
        base_comp = base_comp.copy()

    # Build fragment notation for isotopes
    fragment_isotope_notation: dict[str, int] | int | None = None
    if isotopes is not None:
        fragment_isotope_notation = {}
        for element, count in process_isotopes(isotopes).items():
            mono_element = ELEMENT_LOOKUP[element.symbol]

            if mono_element not in base_comp or base_comp[mono_element] < count:
                raise ValueError(
                    f"Insufficient {mono_element} for isotope adjustment. Composition: { {str(c): i for c, i in base_comp.items()} }, Isotopes: {isotopes}"
                )
            base_comp[mono_element] -= count
            base_comp[element] += count

            elem = str(element)
            if elem not in fragment_isotope_notation:
                fragment_isotope_notation[elem] = 0
            fragment_isotope_notation[elem] += count

        if isinstance(isotopes, int):
            base_comp[ELEMENT_LOOKUP["13C"]] += isotopes
            base_comp[ELEMENT_LOOKUP["C"]] -= isotopes

            fragment_isotope_notation = isotopes

    # Build fragment notation for losses
    fragment_loss_notation: dict[str | float, int] | None = None
    if neutral_deltas is not None:
        loss_comp, float_losses = process_losses(neutral_deltas)

        if float_losses:
            raise ValueError("Cannot adjust composition with float losses.")

        base_comp += loss_comp

        # Build notation from the composition-based losses
        # Only used for record-keeping in Fragment object
        fragment_loss_notation = {}
        for element, count in loss_comp.items():
            elem_str = str(element)
            if elem_str not in fragment_loss_notation:
                fragment_loss_notation[elem_str] = 0
            fragment_loss_notation[elem_str] += count

    ion_info = (
        FRAGMENT_ION_LOOKUP[ion_type]
        if not isinstance(ion_type, FragmentIonInfo)
        else ion_type
    )

    base_comp += ion_info.composition

    # correct for global isotopes
    if isotope_map:
        for original_element, replaced_element in isotope_map.items():
            if original_element in base_comp:
                count = base_comp.pop(original_element)
                base_comp[replaced_element] += count

    # Charge adjustments should not be effected by isotope mapping
    # Build fragment notation for charge adducts
    fragment_annotation_charge_adducts: tuple[str, ...] | None = None
    total_charge: int = 0
    if charge:
        comp_adjustment, charge_value = handle_charge_input_comp(charge=charge)
        base_comp += comp_adjustment
        total_charge = charge_value

        match charge:
            case int():
                fragment_annotation_charge_adducts = None
            case GlobalChargeCarrier() | str():
                fragment_annotation_charge_adducts = (str(charge),)
            case tuple():
                tmp_fragment_annotation_charge_adducts: list[str] = []
                for c in charge:
                    if isinstance(c, GlobalChargeCarrier):
                        c = str(c)
                    tmp_fragment_annotation_charge_adducts.append(c)
                fragment_annotation_charge_adducts = tuple(
                    tmp_fragment_annotation_charge_adducts
                )
            case Mods():
                fragment_annotation_charge_adducts = tuple(
                    str(mod) for mod in charge._mods
                )
            case _:
                raise TypeError(f"Invalid charge type: {type(charge)}")

    # Validate no negative counts
    if any(count < 0 for count in base_comp.values()):
        raise ValueError(f"Negative element counts after adjustments: {base_comp}")

    # Calculate mass from final composition
    base_mass = 0.0
    for elem, count in base_comp.items():
        base_mass += elem.get_mass(monoisotopic=monoisotopic) * count

    # Correct for electron mass based on charge
    if total_charge != 0:
        base_mass -= total_charge * ELECTRON_MASS

    return Fragment(
        ion_type=ion_info.ion_type,
        position=position,
        mass=base_mass,
        charge_state=total_charge,
        monoisotopic=monoisotopic,
        sequence=parent_sequence,
        _charge_adducts=fragment_annotation_charge_adducts,
        _isotopes=fragment_isotope_notation,
        _losses=fragment_loss_notation,
        _comp=base_comp,
    )


def comp_frag(
    comp: Counter[ElementInfo],
    charge: int | tuple[GlobalChargeCarrier | str, ...] | None = None,
    ion_type: IonType | IonTypeLiteral | FragmentIonInfo = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotopes: int | dict[str | ElementInfo, int] | None = None,
    custom_losses: str
    | ChargedFormula
    | float
    | dict[str | ChargedFormula | float, int]
    | None = None,
    position: int | tuple[int, int] | str | None = None,
    parent_sequence: str | None = None,
) -> Fragment:
    mass: float = 0
    for elem, count in comp.items():
        mass += elem.get_mass(monoisotopic=monoisotopic) * count

    return adjust_mass_mz(
        base=mass,
        charge=charge,
        ion_type=ion_type,
        monoisotopic=monoisotopic,
        isotopes=isotopes,
        neutral_deltas=custom_losses,
        position=position,
        parent_sequence=parent_sequence,
    )


def process_isotopes(
    isotopes: int | dict[str | ElementInfo, int],
) -> dict[ElementInfo, int]:
    """Convert isotope input to normalized dict."""
    if isinstance(isotopes, int):
        return {ELEMENT_LOOKUP["13C"]: isotopes}

    result: dict[ElementInfo, int] = {}
    for key, value in isotopes.items():
        element = ELEMENT_LOOKUP[key] if isinstance(key, str) else key
        result[element] = value
    return result


def process_losses(
    losses: str | ChargedFormula | float | dict[str | ChargedFormula | float, int],
) -> tuple[Counter[ElementInfo], dict[float, int]]:
    """Convert loss input to composition counter."""
    if isinstance(losses, str):
        return NEUTRAL_DELTA_LOOKUP[losses].composition, {}
    if isinstance(losses, ChargedFormula):
        return losses.get_composition(), {}
    if isinstance(losses, float):
        return Counter(), {losses: 1}

    # dict case
    total: Counter[ElementInfo] = Counter()
    float_losses: dict[float, int] = {}
    for key, count in losses.items():
        if isinstance(key, str):
            loss_comp = NEUTRAL_DELTA_LOOKUP[key].composition
        elif isinstance(key, ChargedFormula):
            loss_comp = key.get_composition()
        elif isinstance(key, float):
            float_losses[key] = count
            continue
        else:
            raise TypeError(f"Invalid key type for loss: {type(key)}")

        for element, elem_count in loss_comp.items():
            total[element] += elem_count * count
    return total, float_losses


T = TypeVar("T", float, Counter[Any])


@overload
def cumsum(numbers: Sequence[float], reverse: bool = False) -> list[float]: ...


@overload
def cumsum(
    numbers: Sequence[Counter[Any]], reverse: bool = False
) -> list[Counter[Any]]: ...


def cumsum(
    numbers: Sequence[float] | Sequence[Counter[Any]], reverse: bool = False
) -> list[float] | list[Counter[Any]]:
    """Compute cumulative sum of a list of numbers or Counters."""
    match numbers:
        case [Counter(), *_]:
            # Counter case
            total: Counter[Any] = Counter()
            result: list[Counter[Any]] = []
            for counter in reversed(numbers) if reverse else numbers:
                total = total + counter  # type: ignore
                result.append(total.copy())
            return result

        case [float() | int(), *_] | []:
            # Numeric case
            total_num = 0.0
            result_num: list[float] = []
            for number in reversed(numbers) if reverse else numbers:
                total_num += number  # type: ignore
                result_num.append(total_num)
            return result_num

        case _:
            raise TypeError(
                f"cumsum expects sequence of float or Counter, got {type(numbers[0])}"
            )


# Define rules: ion_type -> (position, required_aas, excluded_aas, specific_ion_map)
FRAGMENT_RULES: Any = {
    IonType.D: ("end", None, {"G", "A", "P", "I", "T"}, {"V": IonType.D_VALINE}),
    IonType.DA: (
        "end",
        {"I", "T"},
        None,
        {"I": IonType.DA_ISOLEUCINE, "T": IonType.DA_THREONINE},
    ),
    IonType.DB: (
        "start",
        {"I", "T"},
        None,
        {"I": IonType.DB_ISOLEUCINE, "T": IonType.DB_THREONINE},
    ),
    IonType.D_VALINE: ("end", {"V"}, None, None),
    IonType.DA_THREONINE: ("end", {"T"}, None, None),
    IonType.DA_ISOLEUCINE: ("end", {"I"}, None, None),
    IonType.DB_THREONINE: ("start", {"T"}, None, None),
    IonType.DB_ISOLEUCINE: ("start", {"I"}, None, None),
    IonType.W: ("start", None, {"G", "A", "P", "I", "T"}, {"V": IonType.W_VALINE}),
    IonType.WA: (
        "start",
        {"I", "T"},
        None,
        {"I": IonType.WA_ISOLEUCINE, "T": IonType.WA_THREONINE},
    ),
    IonType.WB: (
        "end",
        {"I", "T"},
        None,
        {"I": IonType.WB_ISOLEUCINE, "T": IonType.WB_THREONINE},
    ),
    IonType.W_VALINE: ("start", {"V"}, None, None),
    IonType.WA_THREONINE: ("start", {"T"}, None, None),
    IonType.WA_ISOLEUCINE: ("start", {"I"}, None, None),
    IonType.WB_THREONINE: ("end", {"T"}, None, None),
    IonType.WB_ISOLEUCINE: ("end", {"I"}, None, None),
}


def can_fragment_sequence(sequence: str, ion_type: IonType | IonTypeLiteral) -> IonType:
    """Check if a sequence can produce a fragment of the given ion type."""

    if not isinstance(ion_type, IonType):
        ion_type = IonType(ion_type)

    first_aa = sequence[0]
    last_aa = sequence[-1]

    if ion_type not in FRAGMENT_RULES:
        return ion_type

    position, required, excluded, specific_map = FRAGMENT_RULES[ion_type]
    aa = last_aa if position == "end" else first_aa

    if excluded and aa in excluded:
        raise ValueError(
            f"{ion_type.name} fragments cannot be produced from sequences {position}ing in {aa}."
        )

    if required and aa not in required:
        raise ValueError(
            f"{ion_type.name} fragments can only be produced from sequences {position}ing in {', or '.join(required)}."
        )

    if specific_map and aa in specific_map:
        return specific_map[aa]

    return ion_type


def combine_loss_types(
    sequence: str,
    losses: NeutralDelta
    | NeutralDeltaLiteral
    | Sequence[NeutralDelta | NeutralDeltaLiteral]
    | None = None,
    custom_losses: str
    | ChargedFormula
    | float
    | dict[str | ChargedFormula | float, int]
    | None = None,
) -> dict[ChargedFormula | str | float, int]:
    all_losses: dict[ChargedFormula | str | float, int] = {}
    if custom_losses is not None:
        if isinstance(custom_losses, (str, ChargedFormula, float)):
            all_losses[custom_losses] = 1
        elif isinstance(custom_losses, dict):
            all_losses.update(custom_losses)
        else:
            raise TypeError(
                "custom_losses must be str, ChargedFormula, float, or dict."
            )

    if losses is not None:
        if isinstance(losses, (NeutralDelta, str)):
            losses = [losses]
        for loss in losses:
            loss_info = NEUTRAL_DELTA_LOOKUP[loss]
            loss_sites = loss_info.calculate_loss_sites(sequence)
            if loss_sites != 0:
                all_losses[loss_info.charged_formula] = loss_sites

    return all_losses
