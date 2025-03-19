"""
Errors
"""


class AmbiguousAminoAcidError(ValueError):
    def __init__(self, aa: str, msg: str, *args):
        self.aa = aa
        self.msg = msg
        message = f"Ambiguous amino acid: {aa}! {msg}"
        super().__init__(message, *args)


class AmbiguousModificationError(ValueError):
    """Exception raised for errors due to ambiguous modifications."""


class AmbiguousSequenceError(ValueError):
    """Exception raised for errors due to ambiguous sequences."""


class InvalidModificationMassError(ValueError):
    """Exception raised for errors due to invalid modifications."""

    def __init__(self, modification_str: str, *args):
        self.modification_str = modification_str
        message = f'Cannot determine mass for modification: "{modification_str}"'
        super().__init__(message, *args)


class InvalidSequenceError(ValueError):
    """Exception raised for errors due to invalid sequences."""


class UnknownModificationError(ValueError):
    """Exception raised for errors due to unknown modifications."""

    def __init__(self, modification, *args):
        self.modification = modification
        message = f"Unknown modification: {modification}"
        super().__init__(message, *args)


class UnknownModificationMassError(ValueError):
    """Exception raised for errors due to unknown masses."""

    def __init__(self, modification, *args):
        self.mass = modification
        message = f"Unknown modification: {modification}"


class UnknownAminoAcidError(ValueError):
    """Exception raised for errors due to unknown amino acids."""

    def __init__(self, amino_acid, *args):
        self.amino_acid = amino_acid
        message = f"Unknown amino acid: {amino_acid}"
        super().__init__(message, *args)


class InvalidDeltaMassError(ValueError):
    """Exception raised for errors due to invalid delta mass."""

    def __init__(self, mass, *args):
        self.mass = mass
        message = f"Invalid delta mass: {mass}"
        super().__init__(message, *args)


class InvalidCompositionError(ValueError):
    """Exception raised for errors due to invalid composition."""

    def __init__(self, composition, *args):
        self.composition = composition
        message = f"Cannot retrieve composition for: {composition}"
        super().__init__(message, *args)


class DeltaMassCompositionError(ValueError):
    """Exception raised for errors due to invalid composition."""

    def __init__(self, composition, *args):
        self.composition = composition
        message = f"Cannot retrieve composition for: {composition}"
        super().__init__(message, *args)


class InvalidChemFormulaError(ValueError):
    """Exception raised for errors due to invalid formula."""

    def __init__(self, formula, msg, *args):
        self.formula = formula
        self.msg = msg
        message = f'Error parsing chem formula: "{formula}". {msg}'
        super().__init__(message, *args)


class InvalidGlycanFormulaError(ValueError):
    """Exception raised for errors due to invalid formula."""

    def __init__(self, formula, msg, *args):
        self.formula = formula
        self.msg = msg
        message = f'Error parsing glycan formula: "{formula}". {msg}'
        super().__init__(message, *args)


class ProFormaFormatError(ValueError):
    """Exception raised for errors due to invalid ProForma format."""

    def __init__(self, msg, index, sequence, *args):
        # Creating a visual marker for the error position in the sequence.
        # Ensure index is within bounds to avoid IndexError.
        if 0 <= index < len(sequence):
            highlighted_sequence = sequence[:index] + ">>>" + sequence[index] + "<<<" + sequence[index + 1:]
        else:
            highlighted_sequence = sequence

        self.msg = msg
        message = (f"Invalid ProForma format detected:\n"
                   f"  Error: {msg}\n"
                   f"  At index: {index}\n"
                   f"  In sequence: {highlighted_sequence}\n"
                   f"  Note: The character in error is indicated by >>>ERROR<<<.")
        super().__init__(message, *args)
