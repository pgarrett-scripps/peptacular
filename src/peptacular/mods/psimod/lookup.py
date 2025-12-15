from typing import Iterable
from .data import PSI_MODIFICATIONS
from ..dclass import PsimodInfo


class PsimodLookup:
    def __init__(self, data: dict[str, PsimodInfo]) -> None:
        self.id_to_info: dict[str, PsimodInfo] = data
        self.name_to_info: dict[str, PsimodInfo] = {
            info.name: info for info in data.values()
        }

        # make keys lowercase for case-insensitive lookup
        self.id_to_info = {k.lower(): v for k, v in self.id_to_info.items()}
        self.name_to_info = {k.lower(): v for k, v in self.name_to_info.items()}

    def query_id(self, mod_id: str | int) -> PsimodInfo | None:
        if isinstance(mod_id, int):
            mod_id = str(mod_id)

        # Strip MOD: or PSI-MOD: prefix if present (case-insensitive)
        if ":" in mod_id:
            prefix, value = mod_id.split(":", 1)
            if prefix.lower() == "mod":
                mod_id = value

        return self.id_to_info.get(mod_id.lower())

    def query_name(self, name: str) -> PsimodInfo | None:
        # Strip MOD: or PSI-MOD: prefix if present (case-insensitive)
        if ":" in name:
            prefix, value = name.split(":", 1)
            if prefix.lower() == "m":
                name = value

        return self.name_to_info.get(name.lower())

    def __getitem__(self, key: str) -> PsimodInfo:
        info = self.query_name(key)
        if info is not None:
            return info

        info = self.query_id(key)
        if info is not None:
            return info

        raise KeyError(f"PSI-MOD modification '{key}' not found by name or ID.")

    def __contains__(self, key: str) -> bool:
        try:
            self[key]
            return True
        except KeyError:
            return False

    def get(self, key: str) -> PsimodInfo | None:
        try:
            return self[key]
        except KeyError:
            return None

    def __iter__(self) -> Iterable[PsimodInfo]:
        """Iterator over all PsimodInfo entries in the lookup."""
        return iter(self.name_to_info.values())


PSIMOD_LOOKUP = PsimodLookup(PSI_MODIFICATIONS)
