
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
        connections: list[bool | None] = []
        for prof_parser, connection in ProFormaParser(sequence).parse():
            
            annot = ProFormaAnnotation(
                sequence="".join(prof_parser.amino_acids),
                isotope_mods=prof_parser.isotope_mods.copy(),
                static_mods=prof_parser.static_mods.copy(),
                labile_mods=prof_parser.labile_mods.copy(),
                unknown_mods=prof_parser.unknown_mods.copy(),
                nterm_mods=prof_parser.nterm_mods.copy(),
                cterm_mods=prof_parser.cterm_mods.copy(),
                internal_mods=prof_parser.internal_mods.copy(),
                intervals=prof_parser.intervals.copy(),
                charge=prof_parser.charge,
                charge_adducts=prof_parser.charge_adducts.copy(),
            )
            annots.append(annot)
            connections.append(connection)

        if len(annots) <= 1:
            raise ValueError("No connections found in the sequence.")

        # Filter out None values from connections
        non_none_connections = [conn for conn in connections if conn is not None]
        if len(non_none_connections) != len(connections):
            raise ValueError("Some connections are undefined in the sequence.")
        
        assert len(non_none_connections) == len(annots) - 1, "Mismatch between annotations and connections."

        return MultiProFormaAnnotation(annotations=annots, connections=non_none_connections)