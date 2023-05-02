from peptacular.protein import digest_protein, peptides_to_df

protein_sequence = 'PEPKTIDEPERPTIDE'

enzyme_regexes = (
    [('([KR])([^P])', 1)],

    [])

missed_cleavages = 1
min_len = 3
max_len = 20
non_enzymatic = False
semi_enzymatic = False

peptides = digest_protein(
    protein_sequence,
    enzyme_regexes,
    missed_cleavages,
    min_len,
    max_len,
    non_enzymatic,
    semi_enzymatic
)

peptide_df = peptides_to_df(peptides)
print(peptide_df)