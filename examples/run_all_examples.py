"""
Run All Examples
================
Execute all example scripts and display their output.
"""

import sys
from pathlib import Path

# Add parent directory to path to import peptacular
examples_dir = Path(__file__).parent
sys.path.insert(0, str(examples_dir.parent / "src"))

def run_example(name: str, description: str):
    """Run a single example and display output."""
    print("\n" + "=" * 80)
    print(f"RUNNING: {name}")
    print(f"Description: {description}")
    print("=" * 80 + "\n")
    
    # Dynamically import and run the example module
    module = __import__(name)
    module.run()
    

def main():
    """Run all examples in sequence."""
    print("=" * 80)
    print("PEPTACULAR EXAMPLES")
    print("=" * 80)
    
    examples = [
        ("proforma_notation", "ProForma 2.0 notation examples"),
        ("annotation", "Creating and manipulating annotations"),
        ("mass_mz_comp", "Mass, m/z, and composition calculations"),
        ("physiochemical_properties", "Physicochemical property calculations"),
        ("converters", "Format conversion (IP2, DIANN, Casanovo, MS2PIP)"),
        ("digest", "Protein digestion and cleavage"),
        ("fragment", "Fragment ion generation"),
        ("parrallel", "Parallel processing examples"),
    ]
    
    for name, description in examples:
        run_example(name, description)
    print("\n" + "=" * 80)
    print("ALL EXAMPLES COMPLETED")
    print("=" * 80)

if __name__ == "__main__":
    main()
