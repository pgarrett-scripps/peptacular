from dataclasses import dataclass
from typing import List

from peptacular.mass import calculate_mz
from peptacular.sequence import get_fragment_sequences, get_internal_fragment_sequences, strip_modifications


@dataclass(frozen=True)
class Fragment:
    """
    Represents a fragment resulting from breaking a sequence at specific points.

    Attributes:
        sequence (str): The sequence of the fragment.
        mass (float): The mass-to-charge (m/z) ratio of the fragment.
        charge (int): The charge state of the fragment.
        ion_type (str): The type of ion formed by the fragmentation.
        number (int): The fragment number (starting from the C-terminus).
        internal (bool): True if the fragment is an internal fragment, False if it is a terminal fragment.
        parent_number (int): The number of the parent fragment from which the internal fragment is derived.
    """

    sequence: str
    mass: float
    charge: int
    ion_type: str
    number: int
    internal: bool
    parent_number: int

    @property
    def start(self):
        if self.ion_type in 'abc':
            if self.internal is False:
                return 0
            else:
                return self.parent_number - self.number
        elif self.ion_type in 'xyz':
            return -self.parent_number

    @property
    def end(self):
        if self.ion_type in 'abc':
            return self.parent_number
        elif self.ion_type in 'xyz':
            if self.internal is False:
                return None
            else:
                return self.start + self.number

    @property
    def label(self):
        return f"{'+' * self.charge}" \
               f"{self.ion_type}" \
               f"{self.parent_number}" \
               f"{'i' if self.internal else ''}" \
               f"{self.number if self.internal else ''}"


def build_fragments(sequence: str, ion_types: List[str], charges: List[int], monoisotopic: bool,
                    internal_fragments: bool) -> List[Fragment]:
    """
    Generate a list of fragments for a given peptide sequence.

    Args:
        sequence (str): The peptide sequence (can be modified).
        ion_types (List[str]): A list of ion types to consider, e.g., ['b', 'y'].
        charges (List[int]): A list of charge states for the fragments.
        monoisotopic (bool): If True, use monoisotopic masses for calculations.
        internal_fragments (bool): If True, include internal fragments; otherwise, only terminal fragments are generated.

    Returns:
        List[Fragment]: A list of Fragment objects representing the generated fragments.

    Example:
        sequence = "PEPTIDE"
        ion_types = ['b', 'y']
        charges = [1, 2]
        monoisotopic = True
        internal_fragments = True

        fragments = build_fragments(sequence, ion_types, charges, monoisotopic, internal_fragments)
        for fragment in fragments:
            print(fragment.sequence, fragment.mass, fragment.charge, fragment.ion_type, fragment.number,
                  fragment.internal, fragment.parent_number)
    """
    fragments = []
    for ion_type in ion_types:
        seq_length = len(strip_modifications(sequence))
        for number, fragment_sequence in enumerate(get_fragment_sequences(sequence=sequence, ion_type=ion_type)):

            frag_number = seq_length - number
            for charge in charges:
                mass = calculate_mz(fragment_sequence, charge=charge, ion_type=ion_type, monoisotopic=monoisotopic)
                fragments.append(Fragment(sequence=fragment_sequence, mass=mass, charge=charge, ion_type=ion_type,
                                          number=frag_number, internal=False, parent_number=frag_number))
            if internal_fragments is True:
                fragment_seq_len = len(strip_modifications(fragment_sequence))
                for internal_number, internal_fragment_sequence in enumerate(
                        get_internal_fragment_sequences(sequence=fragment_sequence, ion_type=ion_type)):

                    if internal_number == 0:  # skip first since its already created
                        continue

                    internal_frag_number = fragment_seq_len - internal_number

                    for charge in charges:
                        mass = calculate_mz(internal_fragment_sequence, charge=charge, ion_type=ion_type,
                                            monoisotopic=monoisotopic)
                        fragments.append(
                            Fragment(sequence=internal_fragment_sequence, mass=mass, charge=charge, ion_type=ion_type,
                                     number=internal_frag_number, internal=True, parent_number=frag_number))

    return fragments
