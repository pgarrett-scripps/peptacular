import regex as re


def convert_ip2_sequence(sequence: str) -> str:
    """
    Converts a IP2-Like sequence to a proforma2.0 compatible sequence.

    :param sequence: The sequence to be converted.
    :type sequence: str

    :return: Proforma2.0 compatable sequence.
    :rtype: str

    .. code-block:: python

        >>> convert_ip2_sequence('K.PEP(phospho)TIDE.K')
        'PEP[phospho]TIDE'

        >>> convert_ip2_sequence('K.(-1)PEP(phospho)TIDE.K')
        '[-1]-PEP[phospho]TIDE'

        >>> convert_ip2_sequence('K.PEPTIDE(2).K')
        'PEPTIDE[2]'

        >>> convert_ip2_sequence('K.PEPTIDE(2)(3).K')
        'PEPTIDE[2]-[3]'

        >>> convert_ip2_sequence('-.(1)PEP(phospho)TIDE(2)(3).-')
        '[1]-PEP[phospho]TIDE[2]-[3]'

        >>> convert_ip2_sequence('P')
        'P'

        >>> convert_ip2_sequence('')
        ''

        >>> convert_ip2_sequence('PEPTIDE')
        'PEPTIDE'
    """

    # Use regex to check if sequence starts and ends with the specified pattern
    if re.match(r"^([A-Z]|-)\..*\.([A-Z]|-)$", sequence):
        # If it matches, remove the leading and trailing characters (first and last two characters)
        sequence = sequence[2:-2]

    # Step 2: Replace () with []
    sequence = re.sub(r"\(([^)]+)\)", r"[\1]", sequence)

    # Step 3: Handle modifications at the start (can be any content, not just numbers)
    sequence = re.sub(r"^\[([^\]]+)\]", r"[\1]-", sequence)

    # Step 4: Convert consecutive modifications to use a dash
    sequence = re.sub(r"\]\[", r"]-[", sequence)

    return sequence


def convert_diann_sequence(sequence: str) -> str:
    """
    Converts a IP2-Like sequence to a proforma2.0 compatible sequence.

    :param sequence: The sequence to be converted.
    :type sequence: str

    :return: Proforma2.0 compatable sequence.
    :rtype: str

    .. code-block:: python

        >>> convert_diann_sequence('_YMGTLRGC[Carbamidomethyl]LLRLYHD_')
        'YMGTLRGC[Carbamidomethyl]LLRLYHD'

        >>> convert_diann_sequence('_[Acytel]YMGTLRGC[Carbamidomethyl]LLRLYHD_')
        '[Acytel]-YMGTLRGC[Carbamidomethyl]LLRLYHD'

        >>> convert_diann_sequence('_[Acytel]YMGTLRGC[Carbamidomethyl]LLRLYHD[1.0]_[Methyl]')
        '[Acytel]-YMGTLRGC[Carbamidomethyl]LLRLYHD[1.0]-[Methyl]'

    """

    # Check if sequence starts and ends with underscores and remove them
    if sequence.startswith("_"):
        sequence = sequence[1:]

        # Check for a modification at the start of the sequence
        if re.match(r"^\[[^\]]+\]", sequence):
            sequence = re.sub(r"^\[([^\]]+)\]", r"[\1]-", sequence)

    if sequence.endswith("_"):
        sequence = sequence[:-1]

    elif re.search(r"_\[[^\]]+\]$", sequence):
        sequence = re.sub(r"_\[([^\]]+)\]$", r"-[\1]", sequence)

    return sequence


def convert_casanovo_sequence(sequence: str) -> str:
    """
    Converts a sequence with modifications to a proforma2.0 compatible sequence.
    :param sequence: The sequence to be converted.
    :type sequence: str
    :return: Proforma2.0 compatable sequence.
    :rtype: str

    .. code-block:: python

        >>> convert_casanovo_sequence('+43.006P+100EPTIDE')
        '[+43.006]-P[+100]EPTIDE'

    """
    new_sequence_comps: list[str] = []
    in_mod = False  # Tracks if we are within a modification
    is_nterm = False  # Tracks if the current modification is at the N-terminus

    for _, char in enumerate(sequence):
        if char in {"+", "-"}:
            # Check if it's at the start (N-terminal)
            is_nterm = len(new_sequence_comps) == 0

            # Start a new modification block
            new_sequence_comps.append("[")
            new_sequence_comps.append(char)
            in_mod = True
        elif in_mod and char.isalpha():
            # End the modification block
            new_sequence_comps.append("]")

            if is_nterm:
                # Add a dash if it's an N-terminal modification
                new_sequence_comps.append("-")
                is_nterm = False

            # Add the current character and close modification
            in_mod = False
            new_sequence_comps.append(char)
        else:
            # Add regular characters
            new_sequence_comps.append(char)

    # Close any unclosed modification at the end of the sequence
    if in_mod:
        new_sequence_comps.append("]")

    return "".join(new_sequence_comps)
