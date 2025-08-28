from __future__ import annotations
import itertools
from typing import TYPE_CHECKING, Any, Generator, Iterable, Literal, Mapping, Optional
import warnings
import regex as re

from .dclasses.modlist import ModList
from ..mod import Mod

if TYPE_CHECKING:
    from .annotation import ProFormaAnnotation


def get_mod_index_from_aa(peptide: str, mod_aa: str) -> set[int]:
    # Given a peptide sequence and a modification, return the indices of the modification
    sites: set[int] = set()
    for aa in mod_aa:
        sites.update(i for i, peptide_aa in enumerate(peptide) if peptide_aa == aa)
    return sites


def get_mod_index_from_regex(peptide: str, mod: str) -> set[int]:
    # given a regex, ensure that it is just a single site
    # then return a set of all indexes

    sites: set[int] = set()
    match = re.search(mod, peptide)
    if match:
        sites.add(match.start())
        if match.end() - match.start() != 0:
            warnings.warn(
                "Using start of regex", UserWarning
            )
    return sites


def get_mod_index(peptide: str, mod: str, is_regex: bool) -> set[int]:
    if is_regex:
        return get_mod_index_from_regex(peptide, mod)
    return get_mod_index_from_aa(peptide, mod)


def get_sites(
    peptide: str, mods: Mapping[str, Iterable[str | float | int | Mod]], is_regex: bool
) -> dict[int, list[str | float | int | Mod]]:
    sites: dict[int, list[str | float | int | Mod]] = {}
    for mod_pattern, mod_values in mods.items():
        mod_sites = get_mod_index(peptide, mod_pattern, is_regex)
        for site in mod_sites:
            if site not in sites:
                sites[site] = []
            sites[site].extend(mod_values)
    return sites


def ensure_single_static_mod(mods: dict[int, list[str | float | int | Mod]]) -> None:
    # warn if multiple static mods are applied to the same site
    for site, mod in mods.items():
        if len(mod) > 1:
            warnings.warn(
                f"Multiple static modifications applied to site {site}: {mod}",
                UserWarning,
            )

def update_mod_list(
    mod_list: ModList,
    merge_strat: Literal["update", "merge", "error"],
    mods: Iterable[str | float | int | Mod],
) -> ModList:
    match merge_strat:
        case "update":
            mod_list.clear()
            mod_list.extend(mods)
        case "merge":
            mod_list.extend(mods)
        case "error":
            if mod_list.has_mods:
                raise ValueError(
                    f"ModList already has modifications: {mod_list}. Cannot apply static mods with 'error' strategy."
                )
            mod_list.extend(mods)
        case _:
            raise ValueError(f"Unknown merge strategy: {merge_strat}")
        
    return mod_list

def apply_mods(
    annotation: ProFormaAnnotation,
    nterm: Mapping[str, Iterable[str | float | int | Mod]] | None = None,
    cterm: Mapping[str, Iterable[str | float | int | Mod]] | None = None,
    internal: Mapping[str, Iterable[str | float | int | Mod]] | None = None,
    inplace: bool = True,
    merge_strat: Literal["update", "merge", "error"] = "error",
    is_regex: bool = False,
) -> ProFormaAnnotation:

    if inplace is False:
        return apply_mods(
            annotation.copy(),
            nterm=nterm,
            cterm=cterm,
            internal=internal,
            inplace=True,
            merge_strat=merge_strat,
        )

    # get sites
    nterm_sites = get_sites(annotation.sequence, nterm, is_regex) if nterm else {}
    cterm_sites = get_sites(annotation.sequence, cterm, is_regex) if cterm else {}
    internal_sites = (
        get_sites(annotation.sequence, internal, is_regex) if internal else {}
    )

    # warn if more than 1 static mod per site
    ensure_single_static_mod(nterm_sites)
    ensure_single_static_mod(cterm_sites)
    ensure_single_static_mod(internal_sites)

    update_mod_list(
        annotation.get_nterm_mod_list(),
        merge_strat=merge_strat,
        mods=nterm_sites.get(0, []),
    )
    update_mod_list(
        annotation.get_cterm_mod_list(),
        merge_strat=merge_strat,
        mods=cterm_sites.get(len(annotation.sequence)-1, []),
    )
    for site, mods in internal_sites.items():
        mod_dict = annotation.get_internal_mod_dict()
        mod_dict[site] = update_mod_list(mod_dict.get(site, ModList()), merge_strat=merge_strat, mods=mods)

    return annotation

def build_mods(
    annotation: ProFormaAnnotation,
    nterm_static: Mapping[str, Iterable[str | float | int | Mod]] | None = None,
    cterm_static: Mapping[str, Iterable[str | float | int | Mod]] | None = None,
    internal_static: Mapping[str, Iterable[str | float | int | Mod]] | None = None,
    labile_static: Mapping[str, Iterable[str | float | int | Mod]] | None = None,
    nterm_variable: Mapping[str, Iterable[str | float | int | Mod]] | None = None,
    cterm_variable: Mapping[str, Iterable[str | float | int | Mod]] | None = None,
    internal_variable: Mapping[str, Iterable[str | float | int | Mod]] | None = None,
    labile_variable: Mapping[str, Iterable[str | float | int | Mod]] | None = None,
    max_variable_mods: int = 2,
    use_regex: bool = False,
) -> Generator[ProFormaAnnotation, None, None]:
    """
    Generate all possible combinations of modifications for a peptide.
    
    Args:
        annotation: The base ProFormaAnnotation to modify
        nterm_static: Static modifications for N-terminus
        cterm_static: Static modifications for C-terminus  
        internal_static: Static modifications for internal residues
        labile_static: Static labile modifications
        nterm_variable: Variable modifications for N-terminus
        cterm_variable: Variable modifications for C-terminus
        internal_variable: Variable modifications for internal residues
        labile_variable: Variable labile modifications
        max_variable_mods: Maximum number of variable modifications to apply
        use_regex: Whether to use regex for modification site matching
        
    Yields:
        ProFormaAnnotation: Each possible modification combination
    """

    # Get all modification sites
    nterm_variable_sites = (
        get_sites(annotation.sequence, nterm_variable, use_regex)
        if nterm_variable
        else {}
    )
    cterm_variable_sites = (
        get_sites(annotation.sequence, cterm_variable, use_regex)
        if cterm_variable
        else {}
    )
    internal_variable_sites = (
        get_sites(annotation.sequence, internal_variable, use_regex)
        if internal_variable
        else {}
    )

    # Get terminal modifications
    nterm_variable_mods: list[str | float | int | Mod] = nterm_variable_sites.get(0, [])
    cterm_variable_mods: list[str | float | int | Mod] = cterm_variable_sites.get(
        len(annotation.sequence) - 1, []
    )

    # Create list of all possible variable modification sites with their mods
    variable_site_mod_pairs: list[tuple[str, int, str | float | int | Mod]] = []
    
    # Add N-term variable mods
    for mod in nterm_variable_mods:
        variable_site_mod_pairs.append(('nterm', 0, mod))
    
    # Add C-term variable mods  
    for mod in cterm_variable_mods:
        variable_site_mod_pairs.append(('cterm', len(annotation.sequence) - 1, mod))
    
    # Add internal variable mods (including terminal positions since they're treated separately)
    for site, mods in internal_variable_sites.items():
        for mod in mods:
            variable_site_mod_pairs.append(('internal', site, mod))
    
    # Add labile variable mods
    if labile_variable:
        labile_variable_sites = get_sites(annotation.sequence, labile_variable, use_regex)
        for site, mods in labile_variable_sites.items():
            for mod in mods:
                variable_site_mod_pairs.append(('labile', site, mod))


    static_modified_annotation = apply_mods(
        annotation,
        nterm=nterm_static,
        cterm=cterm_static, 
        internal=internal_static,
        inplace=False,
        merge_strat="merge",
        is_regex=use_regex
    )

    # Generate all combinations of variable modifications up to max_variable_mods
    for num_var_mods in range(max_variable_mods + 1):
        for var_mod_combination in itertools.combinations(variable_site_mod_pairs, num_var_mods):
            # Check for site conflicts (only one mod per site per type)
            site_conflicts = {}
            has_conflict = False
            
            for mod_type, site, mod in var_mod_combination:
                key = (mod_type, site)
                if key in site_conflicts:
                    has_conflict = True
                    break
                site_conflicts[key] = True
            
            if has_conflict:
                continue
                
            # Start with base annotation and apply static mods
            modified_annotation = static_modified_annotation.copy()
            
            # Apply static labile modifications if any
            if labile_static:
                labile_static_sites = get_sites(annotation.sequence, labile_static, use_regex)
                for site, mods in labile_static_sites.items():
                    for mod in mods:
                        modified_annotation.append_labile_mod(mod)
            
            # Apply variable modifications directly by site
            for mod_type, site, mod in var_mod_combination:
                if mod_type == 'nterm':
                    mod_list = modified_annotation.get_nterm_mod_list()
                    mod_list.append(mod)
                elif mod_type == 'cterm':
                    mod_list = modified_annotation.get_cterm_mod_list()
                    mod_list.append(mod)
                elif mod_type == 'internal':      
                    mod_dict = modified_annotation.get_internal_mod_dict()
                    mod_dict[site] = update_mod_list(mod_dict.get(site, ModList()), merge_strat='merge', mods=[mod])
                elif mod_type == 'labile':
                    modified_annotation.append_labile_mod(mod)

            yield modified_annotation