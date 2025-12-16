from gen_elements import gen as gen_elements
from gen_monosachs import run as gen_monosachs
from gen_unimod import run as gen_unimod
from gen_psimod import run as gen_psimod
from gen_amino_acids import gen as gen_amino_acids
from gen_fragment_ions import gen_fragment_ions
from gen_proteases import gen_proteases
from gen_refmol import gen_refmol
from gen_neutral_deltas import gen_neutral_deltas

def gen():
    print("\n" + "#"*60)
    print("#" + " "*58 + "#")
    print("#" + "  PEPTACULAR DATA GENERATION".center(58) + "#")
    print("#" + " "*58 + "#")
    print("#"*60)
    
    gen_elements()
    gen_monosachs()
    gen_unimod()
    gen_psimod()
    gen_amino_acids()
    gen_fragment_ions()
    gen_proteases()
    gen_refmol()
    gen_neutral_deltas()
    
    print("\n" + "#"*60)
    print("#" + " "*58 + "#")
    print("#" + "  ✅ ALL DATA GENERATION COMPLETE".center(58) + "#")
    print("#" + " "*58 + "#")
    print("#"*60 + "\n")
    
if __name__ == "__main__":
    gen()