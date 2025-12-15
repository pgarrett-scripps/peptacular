from dataclasses import dataclass, field


@dataclass(frozen=True, slots=True, order=True)
class ElementInfo:
    """
    Class to store information about an element isotope
    """

    number: int = field(compare=True)
    mass_number: int = field(compare=True)
    symbol: str = field(compare=False)
    mass: float = field(compare=False)
    abundance: float = field(compare=False)
    average_mass: float = field(
        compare=False
    )  # Average mass for this element (all isotopes)
    is_monoisotopic: bool = field(compare=False)

    @property
    def neutron_count(self) -> int:
        return self.mass_number - self.number

    @property
    def proton_count(self) -> int:
        return self.number

    @property
    def is_radioactive(self) -> bool:
        return self.abundance == 0.0

    def __str__(self) -> str:
        if self.is_monoisotopic:
            return f"{self.symbol}"
        return f"{self.mass_number}{self.symbol}"

    def get_mass(self, monoisotopic: bool = True) -> float:
        """Get the mass of this element isotope"""
        return self.mass if monoisotopic else self.average_mass

    def to_dict(self, float_format: str = "{:.6f}") -> dict[str, object]:
        """Convert the ElementInfo to a dictionary"""
        return {
            "number": self.number,
            "symbol": self.symbol,
            "mass_number": self.mass_number,
            "mass": float_format.format(self.mass),
            "abundance": self.abundance,
            "average_mass": float_format.format(self.average_mass),
        }

    def __repr__(self) -> str:
        return (
            f"ElementInfo(number={self.number}, symbol={self.symbol}, mass_number={self.mass_number}, "
            f"mass={self.mass}, abundance={self.abundance}, average_mass={self.average_mass}, "
            f"is_monoisotopic={self.is_monoisotopic})"
        )
