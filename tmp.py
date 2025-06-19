import peptacular as pt
from peptacular.proforma.proforma_parser import ProFormaAnnotation


print(pt.modification_coverage("PEPTIDE[Phospho]", ["TIDE[Phospho]"]))

print(pt.mass("PEPTIDE/0"), pt.mz("PEPTIDE/0"))
print(pt.mass("PEPTIDE/1"), pt.mz("PEPTIDE/1"))
print(pt.mass("PEPTIDE/2"), pt.mz("PEPTIDE/2"))
print(pt.mass("PEPTIDE/3"), pt.mz("PEPTIDE/3"))

print(pt.mz("PEPTIDE", charge=1, ion_type="y", monoisotopic=True))

print(pt.mass('<13C><[Formula:[13C6]H20]@T>PEPTIDE', precision=3))\

mods = ProFormaAnnotation(sequence='', nterm_mods = [42.0, -20.0])
mods.pop_delta_mass_mods()

print(mods)

print(pt.condense_to_mass_mods('PEPTIDE-[Amidated]'))