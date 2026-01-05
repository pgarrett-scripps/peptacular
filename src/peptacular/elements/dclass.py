from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class ElementInfo:
    """
    Class to store information about an element isotope
    """

    number: int
    mass_number: int | None
    symbol: str
    mass: float
    abundance: float | None
    average_mass: float
    is_monoisotopic: bool | None

    def __hash__(self) -> int:
        return hash(str(self))

    def __eq__(self, other: object) -> bool:
        if isinstance(other, str):
            return str(self) == other
        if isinstance(other, ElementInfo):
            return (self.number, self.mass_number) == (other.number, other.mass_number)
        return NotImplemented

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, ElementInfo):
            return NotImplemented
        mn = self.mass_number if self.mass_number is not None else 0
        on = other.mass_number if other.mass_number is not None else 0
        return (self.number, mn) < (other.number, on)

    def __le__(self, other: object) -> bool:
        if not isinstance(other, ElementInfo):
            return NotImplemented
        mn = self.mass_number if self.mass_number is not None else 0
        on = other.mass_number if other.mass_number is not None else 0
        return (self.number, mn) <= (other.number, on)

    def __gt__(self, other: object) -> bool:
        if not isinstance(other, ElementInfo):
            return NotImplemented
        mn = self.mass_number if self.mass_number is not None else 0
        on = other.mass_number if other.mass_number is not None else 0
        return (self.number, mn) > (other.number, on)

    def __ge__(self, other: object) -> bool:
        if not isinstance(other, ElementInfo):
            return NotImplemented
        mn = self.mass_number if self.mass_number is not None else 0
        on = other.mass_number if other.mass_number is not None else 0
        return (self.number, mn) >= (other.number, on)

    @property
    def neutron_count(self) -> int:
        if self.mass_number is None:
            raise ValueError("Mass number is None, cannot calculate neutron count")
        return self.mass_number - self.number

    @property
    def proton_count(self) -> int:
        return self.number

    @property
    def is_radioactive(self) -> bool:
        return self.abundance == 0.0

    def __str__(self) -> str:
        if self.mass_number is None:
            return f"{self.symbol}"
        return f"{self.mass_number}{self.symbol}"

    def get_mass(self, monoisotopic: bool = True) -> float:
        """Get the mass of this element isotope"""
        return self.mass if monoisotopic else self.average_mass

    def to_dict(self, float_precision: int = 6) -> dict[str, object]:
        """Convert the ElementInfo to a dictionary"""
        return {
            "number": self.number,
            "symbol": self.symbol,
            "mass_number": self.mass_number,
            "mass": round(self.mass, float_precision),
            "abundance": self.abundance,
            "average_mass": round(self.average_mass, float_precision),
        }

    def __repr__(self) -> str:
        return (
            f"ElementInfo(number={self.number}, symbol={self.symbol}, mass_number={self.mass_number}, "
            f"mass={self.mass}, abundance={self.abundance}, average_mass={self.average_mass}, "
            f"is_monoisotopic={self.is_monoisotopic})"
        )

    def update(self, **kwargs: object) -> "ElementInfo":
        """Return a new ElementInfo with updated fields"""
        # Since we use slots=True, we need to get fields manually
        current_values: dict[str, object] = {
            "number": self.number,
            "symbol": self.symbol,
            "mass_number": self.mass_number,
            "mass": self.mass,
            "abundance": self.abundance,
            "average_mass": self.average_mass,
            "is_monoisotopic": self.is_monoisotopic,
        }
        return self.__class__(**{**current_values, **kwargs})  # type: ignore
