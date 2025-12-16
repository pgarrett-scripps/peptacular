from gen_elements import gen as gen_elements
from gen_monosachs import run as gen_monosachs
from gen_unimod import run as gen_unimod
from gen_psimod import run as gen_psimod
from gen_amino_acids import gen as gen_amino_acids
from gen_fragment_ions import gen_fragment_ions
from gen_proteases import gen_proteases
from gen_refmol import gen_refmol
from gen_neutral_deltas import gen_neutral_deltas
import sys
import io
from datetime import datetime
from contextlib import redirect_stdout, redirect_stderr


class TeeOutput:
    """Write to both a file and stdout/stderr"""
    def __init__(self, file_obj, original):
        self.file = file_obj
        self.original = original
    
    def write(self, data):
        self.file.write(data)
        self.original.write(data)
    
    def flush(self):
        self.file.flush()
        self.original.flush()


def gen():
    # Create log file with timestamp
    log_file_path = f"data_gen/generation_log.txt"
    
    with open(log_file_path, 'w', encoding='utf-8') as log_file:
        # Create tee outputs for both stdout and stderr
        tee_stdout = TeeOutput(log_file, sys.stdout)
        tee_stderr = TeeOutput(log_file, sys.stderr)
        
        # Redirect both stdout and stderr to capture everything
        original_stdout = sys.stdout
        original_stderr = sys.stderr
        
        try:
            sys.stdout = tee_stdout
            sys.stderr = tee_stderr
            
            print("\n" + "#"*60)
            print("#" + " "*58 + "#")
            print("#" + "  PEPTACULAR DATA GENERATION".center(58) + "#")
            print("#" + f"  Log file: {log_file_path}".center(58) + "#")
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
            
        except Exception as e:
            print("\n" + "#"*60)
            print("#" + " "*58 + "#")
            print("#" + "  ❌ ERROR DURING GENERATION".center(58) + "#")
            print("#" + " "*58 + "#")
            print("#"*60)
            print(f"\nError: {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()
            raise
        finally:
            # Restore original stdout/stderr
            sys.stdout = original_stdout
            sys.stderr = original_stderr
            
            print(f"\n📝 Generation log saved to: {log_file_path}")
    
if __name__ == "__main__":
    gen()