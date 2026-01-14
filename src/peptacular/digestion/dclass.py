import re
from dataclasses import dataclass
from functools import cached_property


@dataclass(frozen=True)  # Cannot use slots and cached_property together
class ProteaseInfo:
    """Information about a protease enzyme"""

    id: str
    name: str
    full_name: str
    regex: str

    @cached_property
    def pattern(self) -> re.Pattern[str]:
        """Compiled regex pattern for the protease"""
        return re.compile(self.regex)

    def to_dict(self) -> dict[str, object]:
        """Convert the ProteaseInfo to a dictionary"""
        return {
            "id": self.id,
            "name": self.name,
            "full_name": self.full_name,
            "regex": self.regex,
        }
