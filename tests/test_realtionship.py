import unittest
from peptacular.protein import PeptideProteinRelationship  # replace with actual module

class TestPeptideProteinRelationship(unittest.TestCase):

    def setUp(self):
        self.relation = PeptideProteinRelationship()

    def test_add(self):
        protein = 'protein1'
        peptides = ['peptide1', 'peptide2']
        locus = 'locu1'

        self.relation.add(protein, locus, peptides)
        self.assertEqual(self.relation.protein_sequence_to_protein_id[protein], 0)
        self.assertEqual(self.relation.protein_id_to_protein_sequence[0], protein)
        self.assertEqual(self.relation.peptide_sequence_to_peptide_id[peptides[0]], 0)
        self.assertEqual(self.relation.peptide_id_to_peptide_sequence[0], peptides[0])
        self.assertIn(0, self.relation.protein_id_to_peptide_ids[0])
        self.assertIn(0, self.relation.peptide_id_to_protein_ids[0])

    def test_get_proteins(self):
        protein = 'protein1'
        peptides = ['peptide1', 'peptide2']
        locus = 'locu1'

        self.relation.add(protein, locus, peptides)
        proteins = self.relation.get_proteins('peptide1')
        self.assertIn(protein, proteins)

    def test_get_peptides(self):
        protein = 'protein1'
        peptides = ['peptide1', 'peptide2']
        locus = 'locu1'

        self.relation.add(protein, locus, peptides)
        peptides = self.relation.get_peptides('protein1')
        self.assertIn('peptide1', peptides)
        self.assertIn('peptide2', peptides)

    def test_get_proteins_no_peptide(self):
        proteins = self.relation.get_proteins('peptide1')
        self.assertEqual(proteins, [])

    def test_get_peptides_no_protein(self):
        peptides = self.relation.get_peptides('protein1')
        self.assertEqual(peptides, [])


if __name__ == '__main__':
    unittest.main()
