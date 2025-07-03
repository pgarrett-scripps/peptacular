from dataclasses import dataclass
from typing import *

from .annot import ProFormaAnnotation


@dataclass
class MultiProFormaAnnotation:
    """
    A multi proforma annotation
    """

    annotations: List[ProFormaAnnotation]
    connections: List[bool]

    def __post_init__(self):
        if len(self.annotations) - 1 != len(self.connections):
            raise ValueError(
                "Number of connections must be one less than the number of annotations."
            )

    def serialize(
        self, include_plus: bool = False, precision: Optional[int] = None
    ) -> str:
        """
        Convert the multi annotation to a proforma string.

        :return: The serialized multi annotation
        :rtype: str
        """
        seq = ""
        for i, annotation in enumerate(self.annotations):
            seq += annotation.serialize(include_plus=include_plus, precision=precision)
            if i != len(self.annotations) - 1:
                connection = self.connections[i]
                if connection is True:
                    seq += r"\\"
                else:
                    seq += r"+"

        return seq
