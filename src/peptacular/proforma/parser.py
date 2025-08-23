from collections.abc import Generator
from typing import Self

from .dclasses import Mod, Interval, IntervalList, ModList, ModDict
from ..constants import AMINO_ACIDS
from ..errors import ProFormaFormatError
from ..util import validate_single_mod_multiplier

from .characters import ProformaChar as PC


def _is_unmodified(proforma_sequence: str) -> bool:
    """
    Check if a proforma sequence is unmodified

    :param proforma_sequence: The proforma sequence
    :type proforma_sequence: str

    :return: True if the sequence is unmodified
    :rtype: bool
    """
    return all(c in AMINO_ACIDS for c in proforma_sequence)


class ProFormaParser:
    """
    A proforma sequence parser
    """

    def __init__(self, proforma_sequence: str):
        self.sequence: str = proforma_sequence
        self.position: int = 0
        self.length: int = len(proforma_sequence)
        self.amino_acids: list[str] = []
        self.isotope_mods: ModList = ModList()
        self.static_mods: ModList = ModList()
        self.labile_mods: ModList = ModList()
        self.unknown_mods: ModList = ModList()
        self.nterm_mods: ModList = ModList()
        self.cterm_mods: ModList = ModList()
        self.internal_mods: ModDict = ModDict()
        self.charge: int | None = None
        self.charge_adducts: ModList = ModList()
        self.intervals: IntervalList = IntervalList()
        self.current_connection: bool | None = None

    def parse(self) -> Generator[tuple[Self, bool | None], None, None]:
        """
        Parse the proforma sequence, yielding annotations and connections

        :return: A generator of annotations and connections
        :rtype: Generator[(ProFormaAnnotation, bool), None, None]
        """

        if _is_unmodified(self.sequence) is True:
            self.amino_acids = list(self.sequence)
            yield (self, None)
            return


        while not self._end_of_sequence():
            self._parse_sequence_start()
            self._parse_sequence_middle()
            self._parse_sequence_end()

            yield (self, self.current_connection)

            if not self._end_of_sequence():
                self._reset_sequence()

    @property
    def _unmod_sequence(self) -> str:
        return "".join(self.amino_acids)

    def _reset_sequence(self) -> None:
        self.amino_acids = []
        self.isotope_mods.clear()
        self.static_mods.clear()
        self.labile_mods.clear()
        self.unknown_mods.clear()
        self.nterm_mods.clear()
        self.cterm_mods.clear()
        self.internal_mods.clear()
        self.charge = None
        self.charge_adducts.clear()
        self.intervals.clear()

    @validate_single_mod_multiplier
    def _add_static_mod(self, mod: Mod | list[Mod]) -> None:
        if self.static_mods.has_mods:
            self.static_mods.clear()
        if isinstance(mod, list):
            self.static_mods.extend(mod)
        else:
            self.static_mods.append(mod)

    @validate_single_mod_multiplier
    def _add_isotope_mod(self, mod: Mod | list[Mod]) -> None:
        if self.isotope_mods.has_mods:
            self.isotope_mods.clear()
        if isinstance(mod, list):
            self.isotope_mods.extend(mod)
        else:
            self.isotope_mods.append(mod)

    def _add_labile_mod(self, mod: Mod | list[Mod]) -> None:
        if self.labile_mods.has_mods:
            self.labile_mods.clear()
        if isinstance(mod, list):
            self.labile_mods.extend(mod)
        else:
            self.labile_mods.append(mod)

    def _add_unknown_mod(self, mod: Mod | list[Mod]) -> None:
        if not self.unknown_mods.has_mods:
            self.unknown_mods.clear()
        if isinstance(mod, list):
            self.unknown_mods.extend(mod)
        else:
            self.unknown_mods.append(mod)

    def _add_nterm_mod(self, mod: Mod | list[Mod]) -> None:
        if not self.nterm_mods.has_mods:
            self.nterm_mods.clear()
        if isinstance(mod, list):
            self.nterm_mods.extend(mod)
        else:
            self.nterm_mods.append(mod)

    def _add_cterm_mod(self, mod: Mod | list[Mod]) -> None:
        if self.cterm_mods.has_mods:
            self.cterm_mods.clear()
        if isinstance(mod, list):
            self.cterm_mods.extend(mod)
        else:
            self.cterm_mods.append(mod)

    def _add_internal_mod(self, mod: Mod | list[Mod]) -> None:
        if self.internal_mods.has_mods:
            self.internal_mods.clear()
        position = len(self.amino_acids) - 1
        if position not in self.internal_mods:
            self.internal_mods[position] = []
        if isinstance(mod, list):
            self.internal_mods[position].extend(mod)
        else:
            self.internal_mods[position].append(mod)

    def _add_interval(self, interval: Interval | list[Interval]) -> None:
        if self.intervals.has_intervals:
            self.intervals.clear()
        if isinstance(interval, list):
            self.intervals.extend(interval)
        else:
            self.intervals.append(interval)

    @validate_single_mod_multiplier
    def _add_charge_adducts(self, mod: Mod | list[Mod]) -> None:
        if self.charge_adducts.has_mods:
            self.charge_adducts.clear()
        if isinstance(mod, list):
            self.charge_adducts.extend(mod)
        else:
            self.charge_adducts.append(mod)

    def _parse_sequence_start(self) -> None:
        """
        Parse the start of the sequence, up until the first amino acid or interval
        """

        while not self._end_of_sequence():

            cur = self._current()
            if cur in AMINO_ACIDS or cur == PC.INTERVAL_START:  # End of start sequence
                return
            if cur == PC.MOD_START:  # N-term or unknown mods
                mods = self._parse_modifications(PC.MOD_START, PC.MOD_END)
                next_char = self._parse_char()
                if next_char == PC.TERM_MOD:
                    self._add_nterm_mod(mods)
                elif next_char == PC.UNKNOWN:
                    self._add_unknown_mod(mods)
                else:
                    raise ProFormaFormatError(
                        f"Expected '{PC.TERM_MOD}' or '{PC.UNKNOWN}', but got '{next_char}'",
                        self.position,
                        self.sequence,
                    )
            elif cur == PC.STATIC_MOD_START:  # Global mods
                for mod in self._parse_modifications(PC.STATIC_MOD_START, PC.STATIC_MOD_END):
                    if PC.STATIC_SEP in mod.val:  # type: ignore (val will always be str)
                        try:
                            self._add_static_mod(mod)
                        except (
                            ValueError
                        ) as err:  # re-raise error with position, and sequence
                            raise ProFormaFormatError(
                                str(err), self.position, self.sequence
                            ) from err
                    else:  # Isotope mod
                        try:
                            self._add_isotope_mod(mod)
                        except (
                            ValueError
                        ) as err:  # re-raise error with position, and sequence
                            raise ProFormaFormatError(
                                str(err), self.position, self.sequence
                            ) from err
            elif cur == PC.LABILE_MOD_START:  # Labile mods
                self._add_labile_mod(self._parse_modification(PC.LABILE_MOD_START, PC.LABILE_MOD_END))
            else:
                raise ProFormaFormatError(
                    f"Expected amino acid, '{PC.MOD_START}', '{PC.LABILE_MOD_START}', or '{PC.STATIC_MOD_START}' but got '{cur}'",
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
                self.amino_acids.append(self._parse_char())
            elif cur == PC.MOD_START:  # mods for the previous amino acid
                self._add_internal_mod(self._parse_modifications(PC.MOD_START, PC.MOD_END))
            elif cur == PC.TERM_MOD:  # cterm mods (end of sequence)
                self._skip(1)
                self._add_cterm_mod(self._parse_modifications(PC.MOD_START, PC.MOD_END))
                return
            elif cur in (PC.CHIMERIC, PC.CONNECTED):  # charge ( end of sequence)
                return
            elif cur == PC.INTERVAL_START:  # Interval start
                if dummy_interval is not None:
                    raise ProFormaFormatError(
                        "Overlapping intervals!", self.position, self.sequence
                    )
                dummy_interval = Interval(
                    start=len(self.amino_acids),
                    end=len(self.amino_acids),
                    ambiguous=False,
                    mods=None
                )
                self._skip(1)
            elif cur == ")":  # Interval end
                if dummy_interval is None:
                    raise ProFormaFormatError(
                        "Interval ended without starting!", self.position, self.sequence
                    )
                dummy_interval.end = len(self.amino_acids)  # type: ignore (dummy_interval will never be None here)
                

                self._skip(1)
                if not self._end_of_sequence() and self._current() == PC.MOD_START:
                    dummy_interval.mods += self._parse_modifications(PC.MOD_START, PC.MOD_END)

                self._add_interval(
                    dummy_interval
                )
                dummy_interval = None

            elif cur == PC.UNKNOWN:  # unknown mods
                if dummy_interval is None:
                    raise ProFormaFormatError(
                        "Unknown mod outside of interval", self.position, self.sequence
                    )

                # interval is ambiguous
                dummy_interval.ambiguous = True
                self._skip(1)
            else:
                raise ProFormaFormatError(
                    f"Expected either '{PC.MOD_START}', '{PC.INTERVAL_START}', '{PC.UNKNOWN}', '{PC.TERM_MOD}', or '{PC.CONNECTED}' but got: {cur}",
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
                if not self._end_of_sequence() and self._current() == PC.CONNECTED:
                    self._skip(1)
                    self.current_connection = True
                    return

                self.charge = self._parse_integer()

                # check for charge adducts
                if not self._end_of_sequence() and self._current() == PC.MOD_START:
                    self._add_charge_adducts(self._parse_modifications(PC.MOD_START, PC.MOD_END))

            elif cur == PC.CHIMERIC:  # next sequence
                self._skip(1)
                self.current_connection = False
                return
            else:
                raise ProFormaFormatError(
                    f"Invalid sequence: expected '{PC.CONNECTED}' or '{PC.CHIMERIC}' but got {cur}",
                    self.position,
                    self.sequence,
                )

    def _parse_char(self) -> str:
        # Assuming any character not a '[' or ']' is an amino acid for simplicity
        aa = self._current()
        self.position += 1
        return aa

    def _parse_modifications(
        self, opening_bracket: str = PC.MOD_START, closing_bracket: str = PC.MOD_END
    ) -> list[Mod]:
        """
        Parses modifications from the sequence starting with the current position. The function will continue parsing
        until it reaches the end of the sequence or the there are no more sequential modifications.
        """
        mods: list[Mod] = []
        while not self._end_of_sequence():
            if self._current() == opening_bracket:
                mod = self._parse_modification(opening_bracket, closing_bracket)
                mods.append(mod)
            else:
                break

        return mods

    def _parse_modification(
        self, opening_bracket: str = PC.MOD_START, closing_bracket: str = PC.MOD_END
    ) -> Mod:
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
        if not self._end_of_sequence() and self._peek() == PC.MULTI:
            self.position += 1
            multiplier_start = self.position
            while not self._end_of_sequence() and self._peek().isdigit():  # type: ignore (Will always be valid)
                self.position += 1
            multiplier = int(self.sequence[multiplier_start : self.position])

        return Mod(mod, multiplier)

    def _parse_integer(self) -> int:
        start = self.position
        digit_count = 0
        while not self._end_of_sequence():
            if self._peek().isdigit(): # type: ignore (Will always be valid)
                digit_count += 1
                self.position += 1
            elif digit_count == 0 and self._peek() in [PC.PLUS, PC.MINUS]:
                self.position += 1
            else:
                break
        return int(self.sequence[start : self.position])

    def _current(self) -> str:
        return self.sequence[self.position]

    def _peek(self) -> str | None:
        return self.sequence[self.position] if not self._end_of_sequence() else None

    def _skip(self, n: int) -> None:
        self.position += n

    def _end_of_sequence(self) -> bool:
        return self.position >= self.length

