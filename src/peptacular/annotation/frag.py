from dataclasses import dataclass
from typing import Any, Counter, Literal

from ..constants import ELECTRON_MASS
from ..elements import ELEMENT_LOOKUP, ElementInfo
from ..fragment import FRAGMENT_ION_LOOKUP, IonType
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
    UnknownIon,
)
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

    def to_mzpaf(
        self,
        include_annotation: bool = True,
        confidence: float | None = None,
        mass_error: float | None = None,
        mass_error_type: Literal["ppm", "da"] = "ppm",
    ) -> FragmentAnnotation:
        """Convert fragment to mzPAF format string."""

        ion = None
        internal_loss = None
        if self.ion_type is None:
            ion = UnknownIon()
        else:
            ion_info: FragmentIonInfo = FRAGMENT_ION_LOOKUP[self.ion_type]

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
                        if (
                            not isinstance(self.position, tuple)
                            or len(self.position) != 2
                        ):
                            start = -1
                            end = -1
                        else:
                            start, end = self.position

                        internal_ion_key = tuple(list(ion_info.ion_type.value))
                        if internal_ion_key not in INTERNAL_MASS_DIFFS:
                            raise ValueError(
                                f"Internal ion type {ion_info.ion_type} not supported in mzPAF."
                            )
                        internal_loss = INTERNAL_MASS_DIFFS[internal_ion_key]

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
                isotopes = self.isotopes
                if isotopes is None:
                    raise ValueError("Isotopes property is None despite dict type.")
                for element, count in isotopes.items():
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
                        loss.append(
                            f"-{charged_formula.serialize(include_formula_prefix=False)}"
                        )
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
            pzpaf_adducts = None
        else:
            charge_state = self.charge_state
            _adducts = self.charge_adducts
            if _adducts is None:
                raise RuntimeError(
                    "Charge adducts is None despite _charge_adducts being set."
                )
            pzpaf_adducts = [a.to_mz_paf() for a in _adducts]

        mz_paf_mass_error = None
        if mass_error is not None:
            mz_paf_mass_error = MassError(value=mass_error, unit=mass_error_type)

        return FragmentAnnotation(
            ion_type=ion,
            neutral_losses=loss,
            isotope=iso,
            adducts=pzpaf_adducts,
            charge=charge_state,
            mass_error=mz_paf_mass_error,
            confidence=confidence,
        )

    def __str__(self) -> str:
        return (
            self.to_mzpaf().__str__() + f" | Mass: {self.mass:.4f} | m/z: {self.mz:.4f}"
        )

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
