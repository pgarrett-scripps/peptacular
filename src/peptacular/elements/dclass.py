from dataclasses import dataclass, field


@dataclass(frozen=True, slots=True, order=True)
class ElementInfo:
    """
    Class to store information about an element isotope
    """

    number: int = field(compare=True)
    mass_number: int | None = field(compare=True)
    symbol: str = field(compare=False)
    mass: float = field(compare=False)
    abundance: float | None = field(compare=False)
    average_mass: float = field(
        compare=False
    )  # Average mass for this element (all isotopes)
    is_monoisotopic: bool | None = field(compare=False)

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
        return self.__class__(**{**current_values, **kwargs})
