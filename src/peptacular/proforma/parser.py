from collections.abc import Generator
from typing import Self

from ..funcs import convert_type

from .dclasses import Mod, Interval, IntervalList, ModList, ModDict
from ..constants import AMINO_ACIDS
from ..errors import ProFormaFormatError
from ..util import validate_single_mod_multiplier

from .characters import ProformaChar as PC

end_chars = {PC.CHIMERIC, PC.CONNECTED, PC.CHARGE_SEP}


def _is_unmodified(proforma_sequence: str) -> bool:
    """
    Check if a proforma sequence is unmodified

    :param proforma_sequence: The proforma sequence
    :type proforma_sequence: str

    :return: True if the sequence is unmodified
    :rtype: bool
    """
    return not any(c not in AMINO_ACIDS for c in proforma_sequence)


class ProFormaParser:
    """
    A proforma sequence parser
    """

    def __init__(self, proforma_sequence: str):
        self.sequence: str = proforma_sequence
        self.position: int = 0
        self.length: int = len(proforma_sequence)
        self.amino_acids: list[str] = []
        self.isotope_mods: ModList | None = None
        self.static_mods: ModList | None = None
        self.labile_mods: ModList | None = None
        self.unknown_mods: ModList | None = None
        self.nterm_mods: ModList | None = None
        self.cterm_mods: ModList | None = None
        self.internal_mods: ModDict | None = None
        self.charge: int | None = None
        self.charge_adducts: ModList | None = None
        self.intervals: IntervalList | None = None
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
                self.current_connection = None

    @property
    def _unmod_sequence(self) -> str:
        return "".join(self.amino_acids)

    def _reset_sequence(self) -> None:
        self.amino_acids = []
        self.isotope_mods = None
        self.static_mods = None
        self.labile_mods = None
        self.unknown_mods = None
        self.nterm_mods = None
        self.cterm_mods = None
        self.internal_mods = None
        self.charge = None
        self.charge_adducts = None
        self.intervals = None

    @validate_single_mod_multiplier
    def _add_static_mod(self, mod: Mod) -> None:
        if self.static_mods is None:
            self.static_mods = ModList()
        self.static_mods.append(mod)

    @validate_single_mod_multiplier
    def _add_isotope_mod(self, mod: Mod) -> None:
        if self.isotope_mods is None:
            self.isotope_mods = ModList()
        self.isotope_mods.append(mod)

    def _add_labile_mod(self, mods: list[Mod]) -> None:
        if self.labile_mods is None:
            self.labile_mods = ModList()
        self.labile_mods.extend(mods)

    def _add_unknown_mod(self, mods: list[Mod]) -> None:
        if self.unknown_mods is None:
            self.unknown_mods = ModList()
        self.unknown_mods.extend(mods)

    def _add_nterm_mod(self, mods: list[Mod]) -> None:
        if self.nterm_mods is None:
            self.nterm_mods = ModList()
        self.nterm_mods.extend(mods)

    def _add_cterm_mod(self, mods: list[Mod]) -> None:
        if self.cterm_mods is None:
            self.cterm_mods = ModList()
        self.cterm_mods.extend(mods)

    def _add_internal_mod(self, mods: list[Mod]) -> None:
        if self.internal_mods is None:
            self.internal_mods = ModDict()
        position = len(self.amino_acids) - 1
        self.internal_mods.setdefault(position, mods)

    @validate_single_mod_multiplier
    def _add_charge_adducts(self, mods: list[Mod]) -> None:
        if self.charge_adducts is None:
            self.charge_adducts = ModList()
        self.charge_adducts.extend(mods)

    def _add_interval(self, interval: Interval) -> None:
        if self.intervals is None:
            self.intervals = IntervalList()
        self.intervals.append(interval)

    def _parse_sequence_start(self) -> None:
        """Parse the start of the sequence, up until the first amino acid or interval."""
        seq = self.sequence
        seq_len = len(seq)

        while self.position < seq_len:
            cur = seq[self.position]

            # End of start sequence
            if cur in AMINO_ACIDS or cur == PC.INTERVAL_START:
                return

            if cur == PC.MOD_START:
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

            elif cur == PC.STATIC_MOD_START:
                for mod in self._parse_modifications(
                    PC.STATIC_MOD_START, PC.STATIC_MOD_END
                ):
                    if PC.STATIC_SEP in mod.val:  # type: ignore
                        try:
                            self._add_static_mod(mod)
                        except ValueError as err:
                            raise ProFormaFormatError(
                                str(err), self.position, self.sequence
                            ) from err
                    else:
                        try:
                            self._add_isotope_mod(mod)
                        except ValueError as err:
                            raise ProFormaFormatError(
                                str(err), self.position, self.sequence
                            ) from err

            elif cur == PC.LABILE_MOD_START:
                self._add_labile_mod(
                    self._parse_modifications(PC.LABILE_MOD_START, PC.LABILE_MOD_END)
                )

            else:
                raise ProFormaFormatError(
                    f"Expected amino acid, '{PC.MOD_START}', '{PC.LABILE_MOD_START}', or '{PC.STATIC_MOD_START}' but got '{cur}'",
                    self.position,
                    self.sequence,
                )

    def _parse_sequence_middle(self) -> None:
        """Parse the middle of the sequence."""
        dummy_interval = None
        seq = self.sequence
        seq_len = len(seq)

        while self.position < seq_len:
            cur = seq[self.position]

            # Early termination check
            if cur in end_chars:
                return

            if cur in AMINO_ACIDS:
                self.amino_acids.append(self._parse_char())

            elif cur == PC.MOD_START:
                self._add_internal_mod(
                    self._parse_modifications(PC.MOD_START, PC.MOD_END)
                )

            elif cur == PC.TERM_MOD:
                self._skip(1)
                self._add_cterm_mod(self._parse_modifications(PC.MOD_START, PC.MOD_END))
                return

            elif cur == PC.INTERVAL_START:
                if dummy_interval is not None:
                    raise ProFormaFormatError(
                        "Overlapping intervals!", self.position, self.sequence
                    )
                dummy_interval = Interval(
                    start=len(self.amino_acids),
                    end=len(self.amino_acids),
                    ambiguous=False,
                    mods=None,
                )
                self._skip(1)

            elif cur == PC.INTERVAL_END:
                if dummy_interval is None:
                    raise ProFormaFormatError(
                        "Interval ended without starting!", self.position, self.sequence
                    )
                dummy_interval.end = len(self.amino_acids)
                self._skip(1)
                if self.position < seq_len and seq[self.position] == PC.MOD_START:
                    dummy_interval.mods += self._parse_modifications(
                        PC.MOD_START, PC.MOD_END
                    )
                self._add_interval(dummy_interval)
                dummy_interval = None

            elif cur == PC.UNKNOWN:
                if dummy_interval is None:
                    raise ProFormaFormatError(
                        "Unknown mod outside of interval", self.position, self.sequence
                    )
                dummy_interval.ambiguous = True
                self._skip(1)

            else:
                raise ProFormaFormatError(
                    f"Expected either '{PC.MOD_START}', '{PC.INTERVAL_START}', '{PC.UNKNOWN}', '{PC.TERM_MOD}', '{PC.CONNECTED}', '{PC.CHARGE_SEP}' or '{PC.CHIMERIC}' but got: {cur}",
                    self.position,
                    self.sequence,
                )

    def _parse_sequence_end(self) -> None:
        """Parse the end of the sequence."""
        seq = self.sequence
        seq_len = len(seq)

        while self.position < seq_len:
            cur = seq[self.position]

            if cur == PC.CHARGE_SEP:
                self._skip(1)
                # Check for // (crosslink)
                if self.position < seq_len and seq[self.position] == PC.CONNECTED:
                    self._skip(1)
                    self.current_connection = True
                    return
                self.charge = self._parse_integer()
                # Check for charge adducts
                if self.position < seq_len and seq[self.position] == PC.MOD_START:
                    self._add_charge_adducts(
                        self._parse_modifications(PC.MOD_START, PC.MOD_END)
                    )

            elif cur == PC.CHIMERIC:
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
        """Parse single character and advance position."""
        aa = self.sequence[self.position]
        self.position += 1
        return aa

    def _parse_modifications(
        self, opening_bracket: str = PC.MOD_START, closing_bracket: str = PC.MOD_END
    ) -> list[Mod]:
        """Parse modifications from sequence."""
        mods: list[Mod] = []
        seq = self.sequence
        seq_len = len(seq)

        while self.position < seq_len:
            if seq[self.position] == opening_bracket:
                mods.append(self._parse_modification(opening_bracket, closing_bracket))
            else:
                break
        return mods

    def _parse_modification(
        self, opening_bracket: str = PC.MOD_START, closing_bracket: str = PC.MOD_END
    ) -> Mod:
        """Parse a single modification."""
        seq = self.sequence
        seq_len = len(seq)

        self.position += 1
        start = self.position
        bracket_depth = 1

        # Fast bracket matching
        while self.position < seq_len and bracket_depth > 0:
            char = seq[self.position]
            if char == opening_bracket:
                bracket_depth += 1
            elif char == closing_bracket:
                bracket_depth -= 1
            self.position += 1

        if bracket_depth != 0:
            raise ProFormaFormatError(
                f"Unmatched {opening_bracket} at position {self.position}",
                self.position,
                self.sequence,
            )

        mod = seq[start : self.position - 1]
        multiplier = 1

        # Check for multiplier
        if self.position < seq_len and seq[self.position] == PC.MULTI:
            self.position += 1
            multiplier_start = self.position
            while self.position < seq_len and seq[self.position].isdigit():
                self.position += 1
            multiplier = int(seq[multiplier_start : self.position])

        mod_value = convert_type(mod)

        return Mod(mod_value, multiplier)

    def _parse_integer(self) -> int:
        """Parse integer from sequence."""
        seq = self.sequence
        seq_len = len(seq)
        start = self.position
        digit_count = 0

        while self.position < seq_len:
            char = seq[self.position]
            if char.isdigit():
                digit_count += 1
                self.position += 1
            elif digit_count == 0 and char in (PC.PLUS, PC.MINUS):
                self.position += 1
            else:
                break

        return int(seq[start : self.position])

    def _current(self) -> str:
        return self.sequence[self.position]

    def _peek(self) -> str | None:
        return self.sequence[self.position] if not self._end_of_sequence() else None

    def _skip(self, n: int) -> None:
        self.position += n

    def _end_of_sequence(self) -> bool:
        return self.position >= self.length
