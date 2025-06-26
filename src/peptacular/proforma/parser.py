from typing import *

from ..constants import AMINO_ACIDS
from ..errors import ProFormaFormatError
from ..proforma_dataclasses import (
    ACCEPTED_INTERVAL_INPUT,
    ACCEPTED_MOD_INPUT,
    Interval,
    Mod,
    fix_dict_of_mods,
    fix_intervals_input,
    fix_list_of_mods,
)
from ..util import _validate_single_mod_multiplier
from .annotation import ProFormaAnnotation
from .multi_annotation import MultiProFormaAnnotation


def _is_unmodified(proforma_sequence: str) -> bool:
    """
    Check if a proforma sequence is unmodified

    :param proforma_sequence: The proforma sequence
    :type proforma_sequence: str

    :return: True if the sequence is unmodified
    :rtype: bool
    """
    return all(c in AMINO_ACIDS for c in proforma_sequence)


class _ProFormaParser:
    """
    A proforma sequence parser
    """

    def __init__(self, proforma_sequence: str):
        self.sequence = proforma_sequence
        self.position = 0
        self.length = len(proforma_sequence)
        self._amino_acids: List[str] = []
        self._isotope_mods = None
        self._static_mods = None
        self._labile_mods = None
        self._unknown_mods = None
        self._nterm_mods = None
        self._cterm_mods = None
        self._internal_mods = None
        self._charge = None
        self._charge_adducts = None
        self._intervals = None
        self._current_connection = None

    def parse(self) -> Generator[Tuple[ProFormaAnnotation, bool], None, None]:
        """
        Parse the proforma sequence, yielding annotations and connections

        :return: A generator of annotations and connections
        :rtype: Generator[(ProFormaAnnotation, bool), None, None]
        """
        while not self._end_of_sequence():
            self._parse_sequence_start()
            self._parse_sequence_middle()
            self._parse_sequence_end()

            yield self._get_result(), self._current_connection

            if not self._end_of_sequence():
                self._reset_sequence()

    @property
    def _unmod_sequence(self) -> str:
        return "".join(self._amino_acids)

    def _get_result(self) -> ProFormaAnnotation:
        return ProFormaAnnotation(
            sequence=self._unmod_sequence,
            isotope_mods=self._isotope_mods,
            static_mods=self._static_mods,
            labile_mods=self._labile_mods,
            unknown_mods=self._unknown_mods,
            nterm_mods=self._nterm_mods,
            cterm_mods=self._cterm_mods,
            internal_mods=self._internal_mods,
            intervals=self._intervals,
            charge=self._charge,
            charge_adducts=self._charge_adducts,
        )

    def _reset_sequence(self) -> None:
        self._amino_acids = []
        self._isotope_mods = None
        self._static_mods = None
        self._labile_mods = None
        self._unknown_mods = None
        self._nterm_mods = None
        self._cterm_mods = None
        self._internal_mods = None
        self._charge = None
        self._charge_adducts = None
        self._intervals = None

    @_validate_single_mod_multiplier
    def _add_static_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._static_mods is None:
            self._static_mods = []
        if isinstance(mod, list):
            self._static_mods.extend(mod)
        else:
            self._static_mods.append(mod)

    @_validate_single_mod_multiplier
    def _add_isotope_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._isotope_mods is None:
            self._isotope_mods = []
        if isinstance(mod, list):
            self._isotope_mods.extend(mod)
        else:
            self._isotope_mods.append(mod)

    def _add_labile_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._labile_mods is None:
            self._labile_mods = []
        if isinstance(mod, list):
            self._labile_mods.extend(mod)
        else:
            self._labile_mods.append(mod)

    def _add_unknown_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._unknown_mods is None:
            self._unknown_mods = []
        if isinstance(mod, list):
            self._unknown_mods.extend(mod)
        else:
            self._unknown_mods.append(mod)

    def _add_nterm_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._nterm_mods is None:
            self._nterm_mods = []
        if isinstance(mod, list):
            self._nterm_mods.extend(mod)
        else:
            self._nterm_mods.append(mod)

    def _add_cterm_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._cterm_mods is None:
            self._cterm_mods = []
        if isinstance(mod, list):
            self._cterm_mods.extend(mod)
        else:
            self._cterm_mods.append(mod)

    def _add_internal_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._internal_mods is None:
            self._internal_mods = {}
        position = len(self._amino_acids) - 1
        if position not in self._internal_mods:
            self._internal_mods[position] = []
        if isinstance(mod, list):
            self._internal_mods[position].extend(mod)
        else:
            self._internal_mods[position].append(mod)

    def _add_interval(self, interval: Union[Interval, List[Interval]]) -> None:
        if self._intervals is None:
            self._intervals = []
        if isinstance(interval, list):
            self._intervals.extend(interval)
        else:
            self._intervals.append(interval)

    @_validate_single_mod_multiplier
    def _add_charge_adducts(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._charge_adducts is None:
            self._charge_adducts = []
        if isinstance(mod, list):
            self._charge_adducts.extend(mod)
        else:
            self._charge_adducts.append(mod)

    def _parse_sequence_start(self) -> None:
        """
        Parse the start of the sequence, up until the first amino acid or interval
        """

        while not self._end_of_sequence():

            cur = self._current()
            if cur in AMINO_ACIDS or cur == "(":  # End of start sequence
                return
            if cur == "[":  # N-term or unknown mods
                mods = self._parse_modifications("[", "]")
                next_char = self._parse_char()
                if next_char == "-":
                    self._add_nterm_mod(mods)
                elif next_char == "?":
                    self._add_unknown_mod(mods)
                else:
                    raise ProFormaFormatError(
                        f"Expected '-' or '?, but got {cur}",
                        self.position,
                        self.sequence,
                    )
            elif cur == "<":  # Global mods
                for mod in self._parse_modifications("<", ">"):
                    if "@" in mod.val:  # Static mod
                        try:
                            self._add_static_mod(mod)
                        except (
                            ValueError
                        ) as err:  # re-raise error with position, and sequence
                            raise ProFormaFormatError(
                                err, self.position, self.sequence
                            ) from err
                    else:  # Isotope mod
                        try:
                            self._add_isotope_mod(mod)
                        except (
                            ValueError
                        ) as err:  # re-raise error with position, and sequence
                            raise ProFormaFormatError(
                                err, self.position, self.sequence
                            ) from err
            elif cur == "{":  # Labile mods
                self._add_labile_mod(self._parse_modification("{", "}"))
            else:
                raise ProFormaFormatError(
                    r"Expected amino acid, '[', '{', or '<' but got " + cur,
                    self.position,
                    self.sequence,
                )

    def _parse_sequence_middle(self) -> None:
        """
        Parse the middle of the sequence, up until the end of the sequence defined by '/' or '+'
        """
        dummy_interval = None
        while not self._end_of_sequence():
            cur = self._current()
            if cur in AMINO_ACIDS:  # Amino acid
                self._amino_acids.append(self._parse_char())
            elif cur == "[":  # mods for the previous amino acid
                self._add_internal_mod(self._parse_modifications("[", "]"))
            elif cur == "-":  # cterm mods (end of sequence)
                self._skip(1)
                self._add_cterm_mod(self._parse_modifications("[", "]"))
                return
            elif cur in ("/", "+"):  # charge ( end of sequence)
                return
            elif cur == "(":  # Interval start
                if dummy_interval is not None:
                    raise ProFormaFormatError(
                        "Overlapping intervals!", self.position, self.sequence
                    )
                dummy_interval = [len(self._amino_acids), None, False, None]
                self._skip(1)
            elif cur == ")":  # Interval end
                if dummy_interval is None:
                    raise ProFormaFormatError(
                        "Interval ended without starting!", self.position, self.sequence
                    )
                dummy_interval[1] = len(self._amino_acids)

                self._skip(1)
                if not self._end_of_sequence() and self._current() == "[":
                    dummy_interval[3] = self._parse_modifications("[", "]")

                self._add_interval(
                    Interval(
                        start=dummy_interval[0],
                        end=dummy_interval[1],
                        ambiguous=dummy_interval[2],
                        mods=dummy_interval[3],
                    )
                )
                dummy_interval = None

            elif cur == "?":  # unknown mods
                if dummy_interval is None:
                    raise ProFormaFormatError(
                        "Unknown mod outside of interval", self.position, self.sequence
                    )

                # interval is ambiguous
                dummy_interval[2] = True
                self._skip(1)
            else:
                raise ProFormaFormatError(
                    f"Expected either '[', '(', '?', '-', or '/' but got: {cur}",
                    self.position,
                    self.sequence,
                )

    def _parse_sequence_end(self) -> None:
        """
        Parse the end of the sequence, up until the end of the sequence or start of the next sequence
        """
        while not self._end_of_sequence():
            cur = self._current()
            if cur == "/":  # charge
                self._skip(1)

                # check for // (crosslink)
                if not self._end_of_sequence() and self._current() == "/":
                    self._skip(1)
                    self._current_connection = True
                    return

                self._charge = self._parse_integer()

                # check for charge adducts
                if not self._end_of_sequence() and self._current() == "[":
                    self._add_charge_adducts(self._parse_modifications("[", "]"))

            elif cur == "+":  # next sequence
                self._skip(1)
                self._current_connection = False
                return
            else:
                raise ProFormaFormatError(
                    f"Invalid sequence: expected '/' or '+' but got {cur}",
                    self.position,
                    self.sequence,
                )

    def _parse_char(self) -> str:
        # Assuming any character not a '[' or ']' is an amino acid for simplicity
        aa = self._current()
        self.position += 1
        return aa

    def _parse_modifications(
        self, opening_bracket="[", closing_bracket="]"
    ) -> List[Mod]:
        """
        Parses modifications from the sequence starting with the current position. The function will continue parsing
        until it reaches the end of the sequence or the there are no more sequential modifications.
        """
        mods = []
        while not self._end_of_sequence():
            if self._current() == opening_bracket:
                mod = self._parse_modification(opening_bracket, closing_bracket)
                mods.append(mod)
            else:
                break

        return mods

    def _parse_modification(self, opening_bracket="[", closing_bracket="]") -> Mod:
        """
        Parses a single modification from the sequence starting with the current position.
        """
        self.position += 1
        start = self.position
        bracket_depth = 1
        while not self._end_of_sequence() and bracket_depth > 0:
            if self.sequence[self.position] == opening_bracket:
                bracket_depth += 1
            elif self.sequence[self.position] == closing_bracket:
                bracket_depth -= 1
            self.position += 1

        if bracket_depth != 0:
            msg = f"Unmatched {opening_bracket} at position {self.position}"
            raise ProFormaFormatError(msg, self.position, self.sequence)

        mod = self.sequence[start : self.position - 1]

        multiplier = 1
        if not self._end_of_sequence() and self._peek() == "^":
            self.position += 1
            multiplier_start = self.position
            while not self._end_of_sequence() and self._peek().isdigit():
                self.position += 1
            multiplier = int(self.sequence[multiplier_start : self.position])

        return Mod(mod, multiplier)

    def _parse_integer(self) -> int:
        start = self.position
        digit_count = 0
        while not self._end_of_sequence():
            if self._peek().isdigit():
                digit_count += 1
                self.position += 1
            elif digit_count == 0 and self._peek() in ["+", "-"]:
                self.position += 1
            else:
                break
        return int(self.sequence[start : self.position])

    def _current(self) -> str:
        return self.sequence[self.position]

    def _peek(self) -> str:
        return self.sequence[self.position] if not self._end_of_sequence() else None

    def _skip(self, n) -> None:
        self.position += n

    def _end_of_sequence(self) -> bool:
        return self.position >= self.length


def parse(sequence: str) -> Union[ProFormaAnnotation, MultiProFormaAnnotation]:
    """
    Parses a ProForma sequence string and returns its corresponding annotation object.

    Note that the function's behavior and the type of object returned depend on the structure of the input sequence.
    Single sequences result in ProFormaAnnotation objects, while multi-sequences result in MultiProFormaAnnotation
    objects.

    :param sequence: The sequence to parse.
    :type sequence: str

    :raises ProFormaFormatError: If the sequence is not valid.

    :return: Either a ProFormaAnnotation or a MultiProFormaAnnotation, based on the input
    :rtype: Union[ProFormaAnnotation, MultiProFormaAnnotation]

    .. python::

        Parsing a simple peptide sequence:
        >>> isinstance(parse('PEPTIDE'), ProFormaAnnotation)
        True

        Parsing a sequence with multiple peptides or complex modifications:
        >>> isinstance(parse('PEPTIDE+PEPTIDE'), MultiProFormaAnnotation)
        True


    """
    if _is_unmodified(sequence) is True:
        return ProFormaAnnotation(sequence=sequence)

    annotations_connections = list(_ProFormaParser(sequence).parse())
    annotations = [annotation for annotation, _ in annotations_connections]
    annotations_connections = [connection for _, connection in annotations_connections]

    if len(annotations) == 1:
        return annotations[0]

    return MultiProFormaAnnotation(annotations, annotations_connections[:-1])


def serialize(
    annotation: Union[ProFormaAnnotation, MultiProFormaAnnotation],
    include_plus: bool = False,
) -> str:
    """
    Serializes a ProForma annotation or multiple ProForma annotations into a single string representation.

    :param annotation: Either a ProFormaAnnotation or a MultiProFormaAnnotation.
    :type annotation: Union[ProFormaAnnotation, MultiProFormaAnnotation]

    :return: A string representation of the ProForma annotation.
    :rtype: str

    . python::

        Serializing a simple ProForma annotation:
        >>> serialize(ProFormaAnnotation(sequence='PEPTIDE'))
        'PEPTIDE'

        >>> pfa1 = ProFormaAnnotation(sequence='PEPTIDE')
        >>> pfa2 = ProFormaAnnotation(sequence='PEPTIDE')

        Serializing a MultiProFormaAnnotation with chimeric connections:
        >>> multi_annotation = MultiProFormaAnnotation([pfa1, pfa2], [False])
        >>> serialize(multi_annotation)
        'PEPTIDE+PEPTIDE'

        Serializing a MultiProFormaAnnotation with crosslink connections:
        >>> multi_annotation = MultiProFormaAnnotation([pfa1, pfa2], [True])
        >>> p = serialize(multi_annotation)
        >>> p == r'PEPTIDE\\\PEPTIDE'
        True

    """

    return annotation.serialize(include_plus)


def create_annotation(
    sequence: str,
    isotope_mods: Optional[ACCEPTED_MOD_INPUT] = None,
    static_mods: Optional[ACCEPTED_MOD_INPUT] = None,
    labile_mods: Optional[ACCEPTED_MOD_INPUT] = None,
    unknown_mods: Optional[ACCEPTED_MOD_INPUT] = None,
    nterm_mods: Optional[ACCEPTED_MOD_INPUT] = None,
    cterm_mods: Optional[ACCEPTED_MOD_INPUT] = None,
    internal_mods: Optional[Dict[int, ACCEPTED_MOD_INPUT]] = None,
    intervals: Optional[ACCEPTED_INTERVAL_INPUT] = None,
    charge: Optional[int] = None,
    charge_adducts: Optional[ACCEPTED_MOD_INPUT] = None,
) -> ProFormaAnnotation:
    """
    Create a ProFormaAnnotation from a sequence and modifications

    .. code-block:: python

        >>> create_annotation('PEPTIDE', static_mods=['Carbamidomethyl'])
        ProFormaAnnotation(sequence=PEPTIDE, static_mods=[Mod('Carbamidomethyl', 1)])

    """

    isotope_mods = fix_list_of_mods(isotope_mods) if isotope_mods is not None else None
    static_mods = fix_list_of_mods(static_mods) if static_mods is not None else None
    labile_mods = fix_list_of_mods(labile_mods) if labile_mods is not None else None
    unknown_mods = fix_list_of_mods(unknown_mods) if unknown_mods is not None else None
    nterm_mods = fix_list_of_mods(nterm_mods) if nterm_mods is not None else None
    cterm_mods = fix_list_of_mods(cterm_mods) if cterm_mods is not None else None
    internal_mods = (
        fix_dict_of_mods(internal_mods) if internal_mods is not None else None
    )
    intervals = fix_intervals_input(intervals) if intervals is not None else None
    charge_adducts = (
        fix_list_of_mods(charge_adducts) if charge_adducts is not None else None
    )

    return ProFormaAnnotation(
        sequence=sequence,
        isotope_mods=isotope_mods,
        static_mods=static_mods,
        labile_mods=labile_mods,
        unknown_mods=unknown_mods,
        nterm_mods=nterm_mods,
        cterm_mods=cterm_mods,
        internal_mods=internal_mods,
        intervals=intervals,
        charge=charge,
        charge_adducts=charge_adducts,
    )


def create_multi_annotation(
    annotations: List[ProFormaAnnotation], connections: List[bool]
) -> MultiProFormaAnnotation:
    """
    Create a MultiProFormaAnnotation from a list of annotations and connections

    :param annotations: The list of annotations
    :type annotations: List[ProFormaAnnotation]

    :raises ValueError: The number of connections should be one less than the number of annotations

    :param connections: The list of connections
    :type connections: List[bool]

    .. code-block:: python

        >>> annotation = create_multi_annotation([create_annotation('PEP'), create_annotation('TIDE')], [False])
        >>> annotation.annotations[0]
        ProFormaAnnotation(sequence=PEP)
        >>> annotation.annotations[1]
        ProFormaAnnotation(sequence=TIDE)
        >>> annotation.connections
        [False]
    """

    if len(annotations) != len(connections) + 1:
        raise ValueError(
            "The number of connections should be one less than the number of annotations"
        )

    return MultiProFormaAnnotation(annotations=annotations, connections=connections)
