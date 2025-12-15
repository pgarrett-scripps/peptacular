from typing import Iterable
from ..dclass import MonosaccharideInfo
from .data import MONOSACCHARIDES
import warnings


PROFORMA_TO_OBO: dict[str, str] = {
    "aHex": "a-Hex",
    "Dec": "Dec",
    "en,aHex": "en,a-Hex",
    "dHex": "d-Hex",
    "Fuc": "Fuc",
    "Hep": "Hep",
    "Hex": "Hex",
    "HexN": "HexN",
    "HexNAc": "HexNAc",
    "HexNAcS": "HexNAc(S)",
    "HexNS": "HexNS",
    "HexP": "HexP",
    "HexS": "HexS",
    "Neu": "Neu",
    "NeuAc": "Neu5Ac",
    "NeuGc": "Neu5Gc",
    "Non": "Non",
    "Oct": "Oct",
    "Pen": "Pen",
    "Phosphate": "phosphate",
    "Sug": "Sug",
    "Sulfate": "sulfate",
    "Tet": "Tet",
    "Tri": "Tri",
}

OBO_TO_PROFORMA: dict[str, str] = {
    obo_name: proforma_name for proforma_name, obo_name in PROFORMA_TO_OBO.items()
}


class MonosaccharideLookup:
    def __init__(self, monosaccharide_data: dict[str, MonosaccharideInfo]) -> None:
        self.obo_name_to_monosaccharide = monosaccharide_data
        self.proforma_to_monosaccharide: dict[str, MonosaccharideInfo] = {}

        for obo_name, ms in monosaccharide_data.items():
            try:
                proforma_name = OBO_TO_PROFORMA[obo_name]
                self.proforma_to_monosaccharide[proforma_name] = ms.update(
                    name=proforma_name
                )
            except KeyError:
                if (
                    obo_name in ("Acetyl", "Kdn", "Me")
                ):  # Extra monosaccharides (in data_gen/data/additional_monosacharrides.obo.obo)
                    continue

                warnings.warn(
                    f"Monosaccharide OBO name '{obo_name}' (ID: {ms.id}) does not have a "
                    f"corresponding ProForma name mapping. This monosaccharide will not be "
                    f"accessible via ProForma notation.",
                    UserWarning,
                    stacklevel=2,
                )

        # convert all keys to lowercase for case-insensitive lookup
        self.obo_name_to_monosaccharide = {
            k.lower(): v for k, v in self.obo_name_to_monosaccharide.items()
        }
        self.proforma_to_monosaccharide = {
            k.lower(): v for k, v in self.proforma_to_monosaccharide.items()
        }

    def __getitem__(self, key: str) -> MonosaccharideInfo:
        info = self._query_proforma(key)
        if info is not None:
            return info

        info = self._query_obo(key)
        if info is not None:
            return info

        raise KeyError(f"Monosaccharide '{key}' not found by ProForma or OBO name.")

    def __contains__(self, key: str) -> bool:
        try:
            self[key]
            return True
        except KeyError:
            return False

    def get(self, key: str) -> MonosaccharideInfo | None:
        try:
            return self[key]
        except KeyError:
            return None

    def _query_proforma(self, name: str) -> MonosaccharideInfo | None:
        return self.proforma_to_monosaccharide.get(name.lower())

    def _query_obo(self, name: str) -> MonosaccharideInfo | None:
        return self.obo_name_to_monosaccharide.get(name.lower())

    def proforma(self, name: str) -> MonosaccharideInfo:
        val = self._query_proforma(name)
        if val is None:
            raise KeyError(f"Monosaccharide '{name}' not found by ProForma name.")
        return val

    def obo(self, name: str) -> MonosaccharideInfo:
        val = self._query_obo(name)
        if val is None:
            raise KeyError(f"Monosaccharide '{name}' not found by OBO name.")
        return val

    def __iter__(self) -> Iterable[MonosaccharideInfo]:
        """Iterator over all MonosaccharideInfo entries in the lookup."""
        return iter(self.obo_name_to_monosaccharide.values())


MONOSACCHARIDE_LOOKUP = MonosaccharideLookup(
    monosaccharide_data=MONOSACCHARIDES,
)
