from dataclasses import dataclass

from peptacular.proforma.parser import ProFormaParser
from .annotation import ProFormaAnnotation


@dataclass
class MultiProFormaAnnotation:
    """
    A multi proforma annotation
    """

    annotations: list[ProFormaAnnotation]
    connections: list[bool]

    def __post_init__(self):
        if len(self.annotations) - 1 != len(self.connections):
            raise ValueError(
                "Number of connections must be one less than the number of annotations."
            )

    def serialize(
        self, include_plus: bool = False, precision: int | None = None
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

    @staticmethod
    def parse(sequence: str) -> "MultiProFormaAnnotation":
        """
        Parse a ProForma sequence string into a ProFormaAnnotation object.

        :param sequence: The ProForma sequence string to parse.
        :type sequence: str

        :return: A ProFormaAnnotation object representing the parsed sequence.
        :rtype: ProFormaAnnotation
        """
        # Implementation goes here
        annots: list[ProFormaAnnotation] = []
        connections: list[bool] = []
        for prof_parser, connection in ProFormaParser(sequence).parse():
            annot = ProFormaAnnotation(
                sequence="".join(prof_parser.amino_acids),
                isotope_mods=prof_parser.isotope_mods,
                static_mods=prof_parser.static_mods,
                labile_mods=prof_parser.labile_mods,
                unknown_mods=prof_parser.unknown_mods,
                nterm_mods=prof_parser.nterm_mods,
                cterm_mods=prof_parser.cterm_mods,
                internal_mods=prof_parser.internal_mods,
                intervals=prof_parser.intervals,
                charge=prof_parser.charge,
                charge_adducts=prof_parser.charge_adducts,
            )
            annots.append(annot)
            if connection is not None:
                connections.append(connection)

        if len(annots) <= 1:
            raise ValueError("No connections found in the sequence.")

        if len(connections) != len(annots) - 1:
            raise ValueError(
                f"Number of connections must be one less than the number of annotations, got {len(connections)} connections and {len(annots)} annotations."
            )

        return MultiProFormaAnnotation(annotations=annots, connections=connections)
