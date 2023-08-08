def calculate_mass_array(sequence: str, monoisotopic: bool = True):
    if monoisotopic is True:
        aa_masses = MONO_ISOTOPIC_AA_MASSES
    else:
        aa_masses = AVERAGE_AA_MASSES

    sequence_mods = parse_modified_sequence(sequence)
    stripped_sequence = strip_modifications(sequence)

    mass_arr = np.zeros(len(stripped_sequence), dtype=np.float32)
    for i, aa in enumerate(stripped_sequence):

        mod_mass = 0
        if i in sequence_mods:
            mod_mass += float(sequence_mods[i])

        if i == 0 and -1 in sequence_mods:
            mod_mass += float(sequence_mods[-1])

        mass_arr[i] = np.float32(aa_masses[aa] + mod_mass)

    mass_arr = np.array(mass_arr, dtype=np.float32)

    return mass_arr


def create_ion_table(sequence: str, max_len: int = 50, ion_type: str = 'y', charge: int = 0,
                     enzyme: str = 'non-specific'):
    if ion_type == 'y':
        sequence = sequence[::-1]

    sequence_cumulative_sum = np.cumsum(calculate_mass_array(sequence))
    len_sequence = len(sequence)

    # Add padding zeros to make our life easier when computing window sums
    sequence_cumulative_sum = np.pad(sequence_cumulative_sum, (1, 0), 'constant')

    # Window start and end indices
    window_start = np.arange(len_sequence)[:, None]
    window_end = window_start + np.arange(1, max_len + 1)

    # Limit window_end indices to array bounds
    window_end = np.minimum(window_end, len_sequence)
    # Compute window sums
    ion_table = sequence_cumulative_sum[window_end] - sequence_cumulative_sum[window_start]

    # Add (HYDROGEN * 2 + OXYGEN) to the computed mass, only where window_end > window_start
    if ion_type == 'y':
        ion_table[window_end > window_start] += (MONO_ISOTOPIC_ATOMIC_MASSES['HYDROGEN'] * 2 +
                                                 MONO_ISOTOPIC_ATOMIC_MASSES['OXYGEN']) + \
                                                (MONO_ISOTOPIC_ATOMIC_MASSES['PROTON'] * charge)
    elif ion_type == 'b':
        ion_table[window_end > window_start] += (MONO_ISOTOPIC_ATOMIC_MASSES['PROTON'] * charge)
    else:
        raise ValueError("Invalid ion type. Must be either 'y' or 'b'")

    if charge != 0:
        ion_table = ion_table / charge

    # set invalid indexes to 0
    for i, j in enumerate(range(len(ion_table) - 1, len(ion_table) - min(max_len, len_sequence) - 1, -1)):
        ion_table[j, i + 1:] = 0

    if len_sequence < max_len:
        ion_table[:, len_sequence:] = 0

    if enzyme != 'non-specific':
        filter_ion_table(ion_table, sequence, enzyme)

    return ion_table


def filter_ion_table(ion_table: np.ndarray, sequence: str, enzyme: str):
    enzyme_sites = set(identify_cleavage_sites(sequence, enzyme))

    for i in range(len(ion_table)):
        if i not in enzyme_sites:
            ion_table[i, :] = 0


def get_sorted_ion_table_indexes(matrix: np.ndarray) -> np.ndarray:
    # Get 2D indexes without flattening
    # Apply mass filters if provided
    indexes_2d = np.argwhere(matrix)

    # Sort by values at these indexes
    sorted_2d_indexes = indexes_2d[np.argsort(matrix[indexes_2d[:, 0], indexes_2d[:, 1]], kind='mergesort')]
    return sorted_2d_indexes



def test_convert_to_mass_array(self):
    FRAGS = np.array([97.05276384885, 129.04259308796998, 97.05276384885, 101.04767846841,
                      113.08406397713001, 115.02694302383001, 129.04259308796998], dtype=np.float32)
    mass_arr = calculate_mass_array('PEPTIDE')
    for m in FRAGS:
        self.assertTrue(m in mass_arr)

    FRAGS = np.array([197.05276384885, 129.04259308796998, 197.05276384885, 101.04767846841,
                      113.08406397713001, 115.02694302383001, 229.04259308796998], dtype=np.float32)
    mass_arr = calculate_mass_array('(100)PEP(100)TIDE(100)')
    for m in FRAGS:
        self.assertTrue(m in mass_arr)
