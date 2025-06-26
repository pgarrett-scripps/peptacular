# Testing Rules (For AI)

- Tests should be basic (only test one feature & and be relatively short)
- Limit the number of asserts in each test (ideally only assert final sequence)
- Limit the complexity of peptide sequences being tested (dont use super complex peptides with many different modifications)
- Assert equal checks should ideally be between serialize() and a peptide (including modifications) and not checking attributes. This is easier to read.
- Use minimum total tests to get adequete coverage of features