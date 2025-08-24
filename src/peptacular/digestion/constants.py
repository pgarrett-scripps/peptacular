from __future__ import annotations

import regex as re

PROTEASES: dict[str, str] = {
    "arg-c": "(?<=R)",
    "asp-n": "(?=D)",
    "chymotrypsin": "(?<=[FWYL])(?!P)",
    "chymotrypsin/P": "(?<=[FWYL])",
    "promega-chymotrypsin-high-specificity": "(?<=[YFW])",
    "promega-chymotrypsin-low-specificity": "(?<=[YFWLM])",
    "glu-c": "(?<=E)",
    "lys-c": "(?<=K)",
    "lys-n": "(?=K)",
    "proteinase k": "(?<=[AEFILTVWY])",
    "trypsin": "(?<=[KR])(?=[^P])",
    "trypsin/P": "(?<=[KR])",
    "proalanase": "(?<=[PA])",
    "elastase": r"(?<=[AGSVLI])",
    "pepsin": r"(?<=[FLWY])",
    "thermolysin": r"(?<=[LFIAVM])",
    "proalanase-low-specificity": "(?<=[PASG])",
    "non-specific": "()",
    "no-cleave": "_",
}

PROTEASES_COMPILED: dict[str, re.Pattern[str]] = {}
for key, value in PROTEASES.items():
    try:
        PROTEASES_COMPILED[key] = re.compile(value)
    except Exception as e:
        print(f"Error compiling regex for {key}: {e}")
        raise e
