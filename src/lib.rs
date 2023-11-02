use pyo3::prelude::*;
use std::collections::HashMap;


fn extract_c_terminal_modifications(peptide: &str) -> (String, String) {
    // Check if peptide ends with ']'
    if !peptide.ends_with(']') {
        return (peptide.to_string(), String::new());
    }

    let mut open_bracket_index = None;
    let mut bracket_count = 0;

    for (i, c) in peptide.char_indices().rev() {
        match c {
            '[' => {
                bracket_count -= 1;
                if bracket_count == 0 {
                    open_bracket_index = Some(i);
                    break;
                }
            }
            ']' => bracket_count += 1,
            _ => {}
        }
    }

    if let Some(index) = open_bracket_index {
        let modifications = &peptide[index + 1..peptide.len() - 1];
        let remaining_sequence = &peptide[..index];
        (remaining_sequence.to_string(), modifications.to_string())
    } else {
        (peptide.to_string(), String::new())
    }
}

fn extract_n_terminal_modifications(peptide: &str) -> (String, String) {
    // Check if peptide starts with '['
    if !peptide.starts_with('[') {
        return (peptide.to_string(), String::new());
    }

    let mut close_bracket_index = None;
    let mut bracket_count = 0;

    for (i, c) in peptide.char_indices() {
        match c {
            ']' => {
                bracket_count -= 1;
                if bracket_count == 0 {
                    close_bracket_index = Some(i);
                    break;
                }
            }
            '[' => bracket_count += 1,
            _ => {}
        }
    }

    if let Some(index) = close_bracket_index {
        let modifications = &peptide[1..index];
        let remaining_sequence = &peptide[index + 1..];
        (remaining_sequence.to_string(), modifications.to_string())
    } else {
        (peptide.to_string(), String::new())
    }
}


fn extract_modifications(peptide: &str) -> HashMap<isize, String> {
    let mut modifications = HashMap::new();
    let (peptide, c_terminal_mod) = extract_c_terminal_modifications(peptide);
    let (peptide, n_terminal_mod) = extract_n_terminal_modifications(&peptide);

    // Add N-terminal and C-terminal modifications if they exist
    if !n_terminal_mod.is_empty() {
        modifications.insert(-1, n_terminal_mod);
    }

    let mut unmodified_index = 0isize;
    let mut current_modification = String::new();
    let mut parenthesis_count = 0;
    let mut unmodified_length = 0isize;

    for c in peptide.chars() {
        if c == '(' {
            parenthesis_count += 1;
            if parenthesis_count > 1 {
                current_modification.push(c);
            }
        } else if c == ')' {
            parenthesis_count -= 1;
            if parenthesis_count == 0 {
                modifications.insert(unmodified_index-1, current_modification.clone());
                current_modification.clear();
            } else {
                current_modification.push(c);
            }
        } else {
            if parenthesis_count == 0 {
                unmodified_index += 1;
                unmodified_length += 1;
            } else {
                current_modification.push(c);
            }
        }
    }

    if !c_terminal_mod.is_empty() {
        modifications.insert(unmodified_length, c_terminal_mod);
    }

    modifications
}

fn strip_peptide(modified_peptide: &str) -> String {

    let (modified_peptide, c_terminal_mod) = extract_c_terminal_modifications(modified_peptide);
    let (modified_peptide, n_terminal_mod) = extract_n_terminal_modifications(&modified_peptide);

    let mut unmodified_peptide = String::new();
    let mut parenthesis_count = 0;
    let mut square_bracket_count = 0;

    for c in modified_peptide.chars() {
        match c {
            '(' => parenthesis_count += 1,
            ')' => {
                if parenthesis_count > 0 {
                    parenthesis_count -= 1;
                }
            }
            '[' => square_bracket_count += 1,
            ']' => {
                if square_bracket_count > 0 {
                    square_bracket_count -= 1;
                }
            }
            _ => {
                if parenthesis_count == 0 && square_bracket_count == 0 {
                    unmodified_peptide.push(c);
                }
            }
        }
    }

    unmodified_peptide
}


// Python binding function
#[pyfunction]
fn extract_modifications_py(peptide: &str) -> PyResult<HashMap<isize, String>> {
    let modifications = extract_modifications(peptide);
    Ok(modifications) // Wrap the result in `Ok(...)`
}

// Python bindings for strip_peptide
#[pyfunction]
fn strip_peptide_py(modified_peptide: &str) -> PyResult<String> {
    let unmodified_peptide = strip_peptide(modified_peptide);
    Ok(unmodified_peptide) // Wrap the result in `Ok(...)`

}

#[pymodule]
fn rust_bindings(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(extract_modifications_py, m)?)?;
    m.add_function(wrap_pyfunction!(strip_peptide_py, m)?)?;
    Ok(())
}
