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
                        ChargedFormula.from_string(
                            adduct_name, require_formula_prefix=False
                        ),
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
                return {
                    ELEMENT_LOOKUP[elem]: count
                    for elem, count in self._isotopes.items()
                }

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
                loss_formula = ChargedFormula.from_string(
                    loss_name, require_formula_prefix=False
                )
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

    def to_mzpaf(
        self,
        include_annotation: bool = True,
        confidence: float | None = None,
        mass_error: float | None = None,
        mass_error_type: Literal["ppm", "da"] = "ppm",
    ) -> "PafAnnotation":
        """Convert fragment to mzPAF format string."""
        import paftacular as pft

        ion = None
        internal_loss = None
        if self.ion_type is None:
            ion = pft.UnknownIon()
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

                    ion = pft.PeptideIon(
                        series=pft.IonSeries(ion_type),
                        position=position,
                        sequence=self.sequence if include_annotation else None,
                    )
                case IonTypeProperty.INTERNAL:
                    if ion_info.id == IonType.IMMONIUM:
                        if not isinstance(self.position, str):
                            position = "?"
                        else:
                            position = self.position

                        ion = pft.ImmoniumIon(amino_acid=pft.AminoAcids(position))

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
                        if internal_ion_key not in pft.INTERNAL_MASS_DIFFS:
                            raise ValueError(
                                f"Internal ion type {ion_info.ion_type} not supported in mzPAF."
                            )
                        internal_loss: None | str = pft.INTERNAL_MASS_DIFFS[
                            internal_ion_key
                        ]
                        ion = pft.InternalFragment(
                            start_position=start,
                            end_position=end,
                            sequence=self.sequence if include_annotation else None,
                        )
                case IonTypeProperty.INTACT:
                    if ion_info.ion_type == IonType.PRECURSOR:
                        ion = pft.PrecursorIon()
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
        isotopes: list[pft.IsotopeSpecification] = []
        match self._isotopes:
            case dict():
                _isotopes: Mapping[ElementInfo, int] = self.isotopes
                if _isotopes is None:
                    raise ValueError("Isotopes property is None despite dict type.")
                for element, count in _isotopes.items():
                    isotopes.append(
                        pft.IsotopeSpecification(element=str(element), count=count)
                    )
            case int() as iso_count:
                isotopes.append(pft.IsotopeSpecification(count=iso_count))
            case None:
                pass
            case _:
                raise TypeError(f"Invalid isotopes type: {type(self.isotopes)}")

        # handle losses
        losses: list[pft.NeutralLoss] = []
        match frag_losses := self.losses:
            case dict():
                for loss, count in frag_losses.items():
                    if isinstance(loss, float):
                        # add str with +/- sign
                        nloss = pft.NeutralLoss(count=count, base_mass=loss)
                        losses.append(nloss)
                    elif isinstance(loss, ChargedFormula):
                        print(loss.to_mz_paf()[1:])
                        paf_formula = loss.to_mz_paf()
                        sign = paf_formula[0]
                        if sign not in ("+", "-"):
                            raise ValueError(
                                f"Invalid formula sign in loss: {paf_formula}"
                            )
                        mult = 1 if sign == "+" else -1

                        nloss = pft.NeutralLoss(
                            count=count * mult, base_formula=paf_formula[1:]
                        )  #   skip +/- sign
                        losses.append(nloss)
                    else:
                        raise TypeError(f"Invalid loss type in tuple: {type(loss)}")
            case None:
                pass
            case _:
                raise TypeError(f"Invalid losses type: {type(self.losses)}")

        # Add internal loss if applicable
        if internal_loss is not None:
            losses.append(pft.NeutralLoss(count=1, base_formula=internal_loss[1:]))

        pzpaf_adducts: list[pft.Adduct] = []
        if self._charge_adducts is not None:
            _adducts: tuple[GlobalChargeCarrier, ...] = self.charge_adducts
            pzpaf_adducts = [
                pft.Adduct(count=a.occurance, base_formula=a.to_mz_paf()[2:])
                for a in _adducts
            ]

        mz_paf_mass_error = None
        if mass_error is not None:
            mz_paf_mass_error = pft.MassError(value=mass_error, unit=mass_error_type)

        return pft.PafAnnotation(
            ion_type=ion,
            neutral_losses=tuple(losses),
            isotopes=tuple(isotopes),
            adducts=tuple(pzpaf_adducts),
            charge=self.charge_state,
            mass_error=mz_paf_mass_error,
            confidence=confidence,
        )

    @staticmethod
    def from_mzpaf(
        paf: "PafAnnotation" | str, mass: float, monoisotopic: bool = True
    ) -> "Fragment":
        """Create Fragment from mzPAF PafAnnotation."""
        import paftacular as pft

        if isinstance(paf, str):
            paf = pft.parse_single(paf)

        paf: pft.PafAnnotation = cast(pft.PafAnnotation, paf)

        ion_type = None
        position = None
        structural_loss_indices = set()

        match paf.ion_type:
            case pft.PeptideIon() as ion:
                try:
                    ion_type = IonType(ion.series.value)
                except ValueError:
                    pass
                position = ion.position

            case pft.ImmoniumIon() as ion:
                ion_type = IonType.IMMONIUM
                position = ion.amino_acid

            case pft.InternalFragment() as ion:
                position = (ion.start_position, ion.end_position)

                # Default to BY if no structural loss matches
                ion_type = IonType.BY

                # Check for structural losses
                candidates = [
                    it
                    for it in IonType
                    if (info := FRAGMENT_ION_LOOKUP.get(it))
                    and info.properties == IonTypeProperty.INTERNAL
                    and it != IonType.IMMONIUM
                ]

                # We prioritize finding a matching structural loss
                for it in candidates:
                    key = tuple(list(it.value))
                    diff = pft.INTERNAL_MASS_DIFFS.get(key)

                    if diff:
                        sign = diff[0]
                        formula = diff[1:]
                        target_count = 1 if sign == "+" else -1

                        for i, loss in enumerate(paf.neutral_losses):
                            if (
                                loss.base_formula == formula
                                and loss.count == target_count
                            ):
                                ion_type = it
                                structural_loss_indices.add(i)
                                break
                    if structural_loss_indices:
                        break

            case pft.PrecursorIon():
                ion_type = IonType.PRECURSOR

            case pft.UnknownIon():
                pass

        # Isotopes
        isotopes = {}
        c13_count = 0
        only_c13 = True
        for iso in paf.isotopes:
            elem = iso.element if iso.element else "13C"
            if elem == "13C":
                c13_count += iso.count
            else:
                only_c13 = False
            isotopes[elem] = isotopes.get(elem, 0) + iso.count

        final_isotopes = (
            c13_count
            if only_c13 and c13_count > 0
            else (isotopes if isotopes else None)
        )

        # Losses
        deltas = {}
        for i, loss in enumerate(paf.neutral_losses):
            if i in structural_loss_indices:
                continue

            if loss.base_mass is not None:
                # Store as float loss
                # Sign handling: loss.count usually positive for loss?
                # In to_mzpaf: count=count, base_mass=loss.
                # If loss was -18.0 (float), it's a mass loss.
                # Here loss.base_mass is positive usually?
                # pft.NeutralLoss doc: "base_mass: float"
                # If parsed from string "-18.0", base_mass is 18.0?
                # Not sure. Assuming base_mass is the magnitude.
                # count indicates direction?
                # to_mzpaf: count=count (int).
                # If count is 1, and base_mass is 18.0.
                # Fragment expects `deltas` where values are counts.
                # And keys are masses/formulas.
                deltas[loss.base_mass] = abs(loss.count)
            elif loss.base_formula is not None:
                deltas[loss.proforma_formula] = abs(loss.count)

        # Adducts
        adducts = {}
        for adduct in paf.adducts:
            adducts[adduct.proforma_formula] = abs(adduct.count)

        return Fragment(
            ion_type=ion_type,
            position=position,
            mass=mass,
            monoisotopic=monoisotopic,
            charge_state=paf.charge,
            charge_adducts=adducts if adducts else None,
            isotopes=final_isotopes,
            deltas=deltas if deltas else None,
            sequence=getattr(paf.ion_type, "sequence", None),
        )

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
