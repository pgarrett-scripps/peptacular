import peptacular as pt

if __name__ == "__main__":
    # Parse a ProForma sequence
    annot = pt.parse("PEPT[Phospho]IDE-[Acetyl]")

    # Calculate mass and m/z
    mass = annot.mass()
    mz = annot.mz(charge=2)

    print(annot.set_charge(2).set_peptide_name("Peptacular").serialize())