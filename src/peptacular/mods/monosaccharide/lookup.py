from typing import Iterable
from ..dclass import MonosaccharideInfo
from .data import MONOSACCHARIDES, Monosaccharide


class MonosaccharideLookup:
    def __init__(self, monosaccharide_data: dict[str, MonosaccharideInfo]) -> None:
        self.proforma_to_monosaccharide: dict[str, MonosaccharideInfo] = {
            k.lower(): v for k, v in monosaccharide_data.items()
        }

    def __getitem__(self, key: str | Monosaccharide) -> MonosaccharideInfo:
        info: MonosaccharideInfo | None = self._query_proforma(key)
        if info is not None:
            return info

        raise KeyError(f"Monosaccharide '{key}' not found.")

    def __contains__(self, key: str) -> bool:
        try:
            self[key]
            return True
        except KeyError:
            return False

    def get(self, key: str | Monosaccharide) -> MonosaccharideInfo | None:
        try:
            return self[key]
        except KeyError:
            return None

    def _query_proforma(self, name: str) -> MonosaccharideInfo | None:
        return self.proforma_to_monosaccharide.get(name.lower())

    def proforma(self, name: str) -> MonosaccharideInfo:
        val: MonosaccharideInfo | None = self._query_proforma(name)
        if val is None:
            raise KeyError(f"Monosaccharide '{name}' not found by ProForma name.")
        return val

    def __iter__(self) -> Iterable[MonosaccharideInfo]:
        """Iterator over all MonosaccharideInfo entries in the lookup."""
        return iter(self.proforma_to_monosaccharide.values())


monos: dict[str, MonosaccharideInfo] = {
    str(mono): info for mono, info in MONOSACCHARIDES.items()
}

MONOSACCHARIDE_LOOKUP = MonosaccharideLookup(
    monosaccharide_data=monos,
)
