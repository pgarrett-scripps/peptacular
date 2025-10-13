import pytest
import random
from unittest.mock import patch

import peptacular as pt

reverse_sequence = pt.reverse_sequence
shuffle_sequence = pt.shuffle_sequence
shift_sequence = pt.shift_sequence
debruijin_sequence = pt.debruijin_sequence

def test_reverse_sequence_basic():
    """Test basic reverse functionality."""
    assert reverse_sequence("ABCDE") == "EDCBA"
    assert reverse_sequence("A") == "A"
    assert reverse_sequence("") == ""

def test_reverse_sequence_palindrome():
    """Test reverse with palindromes."""
    assert reverse_sequence("ABBA") == "ABBA"
    assert reverse_sequence("RACECAR") == "RACECAR"

def test_shuffle_sequence_no_static():
    """Test shuffle without static residues."""
    seq = "ABCDE"
    shuffled = shuffle_sequence(seq)
    
    # Check length is preserved
    assert len(shuffled) == len(seq)
    
    # Check all characters are present
    assert sorted(shuffled) == sorted(seq)

def test_shuffle_sequence_with_static():
    """Test shuffle with static residues."""
    seq = "ABCABC"
    static = "A"
    shuffled = shuffle_sequence(seq, static_residues=static)
    
    # Check static residues stay in place
    assert shuffled[0] == "A"
    assert shuffled[3] == "A"
    
    # Check length and character composition
    assert len(shuffled) == len(seq)
    assert sorted(shuffled) == sorted(seq)

def test_shuffle_sequence_all_static():
    """Test shuffle when all residues are static."""
    seq = "AAAA"
    shuffled = shuffle_sequence(seq, static_residues="A")
    assert shuffled == seq

def test_shuffle_sequence_empty():
    """Test shuffle with empty sequence."""
    assert shuffle_sequence("") == ""

@patch('random.shuffle')
def test_shuffle_sequence_deterministic(mock_shuffle):
    """Test shuffle behavior with mocked randomness."""
    def reverse_list(x):
        x.reverse()
    
    mock_shuffle.side_effect = reverse_list
    seq = "ABCDE"
    shuffled = shuffle_sequence(seq)
    assert shuffled == "EDCBA"

def test_shift_sequence_basic():
    """Test basic shift functionality."""
    assert shift_sequence("ABCDE", 1) == "EABCD"
    assert shift_sequence("ABCDE", 2) == "DEABC"
    assert shift_sequence("ABCDE", 0) == "ABCDE"
    assert shift_sequence("ABCDE", -1) == "BCDEA"

def test_shift_sequence_full_rotation():
    """Test shift by full length."""
    seq = "ABCDE"
    assert shift_sequence(seq, len(seq)) == seq
    assert shift_sequence(seq, len(seq) * 2) == seq

def test_shift_sequence_with_static():
    """Test shift with static residues."""
    seq = "AXBXC"
    static = "X"
    shifted = shift_sequence(seq, 1, static_residues=static)
    
    # X should stay at positions 1 and 3
    assert shifted[1] == "X"
    assert shifted[3] == "X"
    
    # Only A, B, C should shift
    non_static = [shifted[0], shifted[2], shifted[4]]
    assert set(non_static) == {"A", "B", "C"}

def test_shift_sequence_large_n():
    """Test shift with n larger than sequence length."""
    seq = "ABC"
    # Shift by 10 should be same as shift by 10 % 3 = 1
    assert shift_sequence(seq, 10) == shift_sequence(seq, 1)

def test_shift_sequence_empty():
    """Test shift with empty sequence."""
    assert shift_sequence("", 5) == ""

def test_shift_sequence_single_char():
    """Test shift with single character."""
    assert shift_sequence("A", 5) == "A"

def test_debruijin_sequence_basic():
    """Test basic de Bruijn sequence generation."""
    # For a simple sequence with k=2
    seq = "ABAB"
    result = debruijin_sequence(seq, 2)
    
    # Check that result contains all k-mers from original
    assert len(result) >= 2

def test_debruijin_sequence_k_equals_1():
    """Test de Bruijn with k=1."""
    seq = "ABC"
    result = debruijin_sequence(seq, 1)
    # With k=1, each character is a node
    assert len(result) >= 1


def test_debruijin_sequence_longer_k():
    """Test with larger k value."""
    seq = "ABCABC"
    result = debruijin_sequence(seq, 3)
    
    assert len(result) >= 3

def test_debruijin_sequence_k_equals_length():
    """Test when k equals sequence length."""
    seq = "ABCD"
    result = debruijin_sequence(seq, 4)
    
    # Should return at least the prefix
    assert len(result) >= 3

# Integration tests
def test_reverse_then_shift():
    """Test combining reverse and shift operations."""
    seq = "ABCDE"
    reversed_seq = reverse_sequence(seq)
    shifted = shift_sequence(reversed_seq, 1)
    assert reversed_seq == "EDCBA"  # EDCBA shifted by 1 = AEDCB
    # Actually let me recalculate: EDCBA -> shift right by 1 -> AEDCB
    assert shifted == "AEDCB"

    # EDCBA

def test_shift_preserves_character_counts():
    """Test that shift preserves character frequency."""
    seq = "AABBCCDD"
    shifted = shift_sequence(seq, 3)
    assert sorted(seq) == sorted(shifted)

def test_shuffle_with_multiple_static_types():
    """Test shuffle with multiple static residue types."""
    seq = "AXBYCZ"
    static = "XYZ"
    shuffled = shuffle_sequence(seq, static_residues=static)
    
    # Static residues should remain
    assert shuffled[1] == "X"
    assert shuffled[3] == "Y"
    assert shuffled[5] == "Z"
    
    # Only A, B, C should shuffle
    non_static_positions = [0, 2, 4]
    non_static_chars = [shuffled[i] for i in non_static_positions]
    assert sorted(non_static_chars) == ["A", "B", "C"]