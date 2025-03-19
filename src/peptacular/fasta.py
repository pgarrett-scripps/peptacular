"""
Fasta.py
"""


def parse_fasta(input_data):
    """
    Parse FASTA formatted data from various input types.
    
    Parameters:
    -----------
    input_data : str, pathlib.Path, or file-like object
        The input can be:
        - A string containing FASTA formatted text
        - A path to a FASTA file (as string or Path object)
        - A file-like object (already opened file or StringIO)
    
    Returns:
    --------
    list of tuples
        Each tuple contains (header, sequence)

    """

    text = ""

    # Check if input is a file-like object with 'read' method
    if hasattr(input_data, 'read'):
        content = input_data.read()
        # Handle case where read() returns bytes (like with Streamlit's UploadedFile)
        if isinstance(content, bytes):
            text = content.decode('utf-8')
        else:
            text = content

    # Check if input is a string
    elif isinstance(input_data, str):
        # Try to open as file path if it doesn't look like FASTA content
        if not input_data.lstrip().startswith(">") and "\n" not in input_data[:100]:
            try:
                with open(input_data, 'r') as f:
                    text = f.read()
            except (FileNotFoundError, IOError):
                # If file not found, assume it's FASTA text
                text = input_data
        else:
            text = input_data

    # Handle pathlib.Path objects
    elif hasattr(input_data, 'is_file') and hasattr(input_data, 'open'):
        try:
            with input_data.open('r') as f:
                text = f.read()
        except IOError as err:
            raise ValueError(f"Could not open file: {input_data}") from err

    else:
        raise TypeError("Input must be a string, file path, or file-like object")

    # Now parse the text using the existing logic
    sequences = []
    header = None
    seq = ""
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if header is not None:
                sequences.append((header, seq))
            header = line[1:].strip()  # Remove '>' and extra whitespace
            seq = ""
        else:
            seq += line.upper()  # Convert to uppercase for consistency
    if header is not None:
        sequences.append((header, seq))
    return sequences
