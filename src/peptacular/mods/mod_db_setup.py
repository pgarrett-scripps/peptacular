"""
mode_db_setup.py
"""

import dataclasses
import os
import pickle
import re
import tempfile
import warnings
from collections import Counter
from functools import cached_property
from typing import List, Dict, IO, Any, Union, Iterator, Tuple

from peptacular.constants import ISOTOPIC_ATOMIC_MASSES
from peptacular.chem.chem_util import write_chem_formula, _parse_isotope_component, parse_chem_formula, chem_mass, \
    _parse_split_chem_formula
from peptacular.errors import InvalidChemFormulaError, InvalidGlycanFormulaError
from peptacular.util import convert_type
from peptacular.types import ChemComposition


@dataclasses.dataclass
class ModEntry:
    id: str
    name: str
    mono_mass: Union[float, None]
    avg_mass: Union[float, None]
    composition: Union[str, None]
    synonyms: Union[List[str], None] = None
    parents: Union[List[str], None] = None
    entry_type: Union[Any, None] = None

    @cached_property
    def calc_mono_mass(self):
        if self.composition is not None:
            return chem_mass(self.composition, monoisotopic=True)
        return None

    @cached_property
    def calc_avg_mass(self):
        if self.composition is not None:
            return chem_mass(self.composition, monoisotopic=False)
        return None

    def dict(self):
        return {
            'id': self.id,
            'name': self.name,
            'mono_mass': self.mono_mass,
            'avg_mass': self.avg_mass,
            'composition': self.composition,
            'synonyms': self.synonyms,
            'parents': self.parents,
            'entry_type': self.entry_type
        }


def _parse_glycan_formula(formula: str, sep: str) -> ChemComposition:
    """
    Here to avoid a circular import. Real function is in peptacular.glycan
    """
    d = {}

    if formula == '':
        return d

    original_formula = formula

    if sep != '':
        return _parse_split_chem_formula(formula, sep)

    while formula != '':
        for glycan_name in MONOSACCHARIDES_DB.names_sorted:
            if formula.startswith(glycan_name):
                formula = formula[len(glycan_name):]
                count = ''

                # get the count (can have +- or digits)
                for c in formula:
                    if c.isdigit() or c in '+-.':
                        count += c
                    else:
                        break

                # remove the count from the formula
                formula = formula[len(count):]

                # add the count to the dictionary
                if count == '':
                    count = 1

                count = convert_type(count)

                if isinstance(count, str):
                    raise InvalidGlycanFormulaError(original_formula, f'Invalid count: "{count}"!')

                d[glycan_name] = count

                break
        else:
            raise InvalidGlycanFormulaError(original_formula, f'Unknown glycan: "{formula}"!')

    return d


def _glycan_comp(glycan: Union[ChemComposition, str], sep: str = '') -> ChemComposition:
    """
    Here to avoid a circular import. Real function is in peptacular.glycan
    """
    if isinstance(glycan, str):
        glycan = _parse_glycan_formula(glycan, sep)  # raises InvalidGlycanFormulaError

    counts = {}
    for component, count in glycan.items():

        if MONOSACCHARIDES_DB.contains_name(component):
            entry = MONOSACCHARIDES_DB.get_entry_by_name(component)
        elif MONOSACCHARIDES_DB.contains_synonym(component):
            entry = MONOSACCHARIDES_DB.get_entry_by_synonym(component)
        else:
            raise InvalidGlycanFormulaError(glycan, f'Unknown glycan: "{component}"!')

        chem_formula = parse_chem_formula(entry.composition)
        for element, element_count in chem_formula.items():
            counts[element] = counts.get(element, 0) + element_count * count

    return counts


class DbType:
    UNIMOD = 'unimod'
    PSI_MOD = 'psi-mod'
    GNO = 'gno'
    RESID = 'resid'
    XLMOD = 'xlmod'
    MONOSACCHARIDES = 'monosaccharide'


_DB_FILE_PATHS = {
    DbType.UNIMOD: 'unimod',
    DbType.PSI_MOD: 'psi-mod',
    DbType.GNO: 'gno',
    DbType.RESID: 'psi-mod',
    DbType.XLMOD: 'xlmod',
    DbType.MONOSACCHARIDES: 'monosaccharides'
}

_OBO_FILE_SOURCES = {
    DbType.UNIMOD: 'https://www.unimod.org/obo/unimod.obo',
    DbType.PSI_MOD: 'https://raw.githubusercontent.com/HUPO-PSI/psi-mod-CV/master/PSI-MOD.obo',
    DbType.XLMOD: 'https://raw.githubusercontent.com/HUPO-PSI/xlmod-CV/main/XLMOD.obo',
    DbType.MONOSACCHARIDES: 'https://raw.githubusercontent.com/HUPO-PSI/ProForma/master/monosaccharides/mono.obo',
    DbType.RESID: 'https://raw.githubusercontent.com/HUPO-PSI/psi-mod-CV/master/PSI-MOD.obo',
    DbType.GNO: 'http://purl.obolibrary.org/obo/gno.obo',
}


def _read_obo(file: IO[str]) -> List[Dict[str, Any]]:
    file.seek(0)

    info, elems, skip, d = {}, [], None, None

    for line in file:

        line = line.rstrip()
        if line == '':
            continue

        if line.startswith('[Typedef]'):
            skip = True
            continue

        if line.startswith('[Term]'):
            skip = False
            if d is not None:
                elems.append(d)
            d = {}
            continue

        if d is None:
            key, value = line.split(': ', 1)
            info[key] = value
            continue

        if skip:
            continue

        if line:
            key, value = line.split(': ', 1)

            if key not in d:
                d[key] = [value]
            else:
                d[key].append(value)

    if d is not None:
        elems.append(d)

    return elems


def _is_obsolete(term: Dict[str, Any]) -> bool:
    is_obsolete = term.get('is_obsolete', ['false'])[0]
    if is_obsolete in ['true', 'false']:
        is_obsolete = bool(is_obsolete == 'true')
    else:
        raise ValueError(f'Invalid is_obsolete value: {is_obsolete}')

    return is_obsolete


def _fix_entry(entry: ModEntry):
    if entry.id is None:
        raise ValueError('Entry id is None')

    if entry.name is None:
        raise ValueError('Entry name is None')

    if entry.composition is not None:
        # try to parse formula
        try:
            _ = parse_chem_formula(entry.composition)
        except Exception as err:
            warnings.warn(f'Error parsing {entry.id} {entry.name} {entry.composition}, {err}')
            entry.composition = None

    if entry.mono_mass is None and entry.composition is not None:
        try:
            entry.mono_mass = chem_mass(entry.composition, monoisotopic=True)
        except Exception as err:
            raise ValueError(f'Error parsing {entry.id} {entry.name} {entry.composition}, {err}') from err
            entry.mono_mass = None

    if entry.avg_mass is None and entry.composition is not None:
        try:
            entry.avg_mass = chem_mass(entry.composition, monoisotopic=False)
        except Exception as err:
            raise ValueError(f'Error parsing {entry.id} {entry.name} {entry.composition}, {err}') from err
            entry.avg_mass = None


def _get_id_and_name(term: Dict[str, Any]) -> (str, str):
    term_id = term.get('id', [])
    term_name = term.get('name', [])

    if len(term_id) > 1:
        warnings.warn(f'Multiple ids for {term_name} {term_id}')
        term_id = term_id[0]
    elif len(term_id) == 1:
        term_id = term_id[0]
    else:
        raise ValueError('Entry name is None')

    if len(term_name) > 1:
        warnings.warn(f'Multiple names for {term_id} {term_name}')
        term_name = term_name[0]
    elif len(term_name) == 1:
        term_name = term_name[0]
    else:
        raise ValueError('Entry id is None')

    return term_id, term_name


def _get_parent_ids(term: Dict[str, Any]) -> List[str]:
    parents = term.get('is_a', [])
    if not isinstance(parents, List):
        parents = [parents]

    parent_ids = [p.split(' ')[0] for p in parents]
    return parent_ids


def _get_unimod_entries(terms: List[Dict[str, Any]]) -> List[ModEntry]:
    for term in terms:

        term_id, term_name = _get_id_and_name(term)
        term_id = term_id.replace('UNIMOD:', '')
        parent_ids = [pid.replace('UNIMOD:', '') for pid in _get_parent_ids(term)]

        if _is_obsolete(term):
            continue

        if term_name == 'unimod root node':
            continue

        property_values = {}
        for val in term.get('xref', []):
            elems = val.split('"')
            k = elems[0].rstrip()
            v = elems[1].strip()

            property_values.setdefault(k, []).append(v)

        delta_composition = property_values.get('delta_composition', [None])
        delta_monoisotopic_mass = property_values.get('delta_mono_mass', [None])
        delta_average_mass = property_values.get('delta_avge_mass', [None])

        if len(delta_composition) > 1:
            warnings.warn(
                f'[{DbType.UNIMOD}] Multiple delta compositions for {term_id} {term_name} {delta_composition}')
            delta_composition = delta_composition[0]
        elif len(delta_composition) == 1:
            delta_composition = delta_composition[0]

        if len(delta_monoisotopic_mass) > 1:
            warnings.warn(
                f'[{DbType.UNIMOD}] Multiple delta mono masses for {term_id} {term_name} {delta_monoisotopic_mass}')
            delta_monoisotopic_mass = delta_monoisotopic_mass[0]
        elif len(delta_monoisotopic_mass) == 1:
            delta_monoisotopic_mass = delta_monoisotopic_mass[0]

        if len(delta_average_mass) > 1:
            warnings.warn(
                f'[{DbType.UNIMOD}] Multiple delta average masses for {term_id} {term_name} {delta_average_mass}')
            delta_average_mass = delta_average_mass[0]
        elif len(delta_average_mass) == 1:
            delta_average_mass = delta_average_mass[0]

        synonyms = term.get('synonym', [])
        if not isinstance(synonyms, List):
            synonyms = [synonyms]
        synonyms = [syn.split('"')[1] for syn in synonyms]

        delta_formula = None
        if delta_composition is not None and isinstance(delta_composition, str):

            delta_composition = delta_composition.replace('(', ' ').replace(')', '')
            delta_composition = parse_chem_formula(delta_composition, sep=' ')

            combined_comp = {}
            for k, v in delta_composition.items():
                if k =='Me' or k == 'Ac' or k not in ISOTOPIC_ATOMIC_MASSES:
                    glycan_composition = _glycan_comp(k)
                    for elem, count in glycan_composition.items():
                        combined_comp[elem] = combined_comp.get(elem, 0) + count * v
                else:
                    combined_comp[k] = combined_comp.get(k, 0) + v

            delta_formula = write_chem_formula(combined_comp)

        if delta_monoisotopic_mass is not None:
            delta_monoisotopic_mass = float(delta_monoisotopic_mass)

        if delta_average_mass is not None:
            delta_average_mass = float(delta_average_mass)

        try:
            _ = chem_mass(delta_formula)
        except InvalidChemFormulaError as err:
            orig_comp = property_values.get('delta_composition', [None])
            warnings.warn(f'Cannot parse Unimod: {term_id}, {delta_formula}, {orig_comp}. {err}')

        mod = ModEntry(
            id=term_id,
            name=term_name,
            mono_mass=delta_monoisotopic_mass,
            avg_mass=delta_average_mass,
            composition=delta_formula,
            synonyms=synonyms,
            parents=parent_ids,
            entry_type=DbType.UNIMOD
        )

        yield mod


def _get_psimod_entries(terms: List[Dict[str, Any]]) -> List[ModEntry]:
    for term in terms:

        term_id, term_name = _get_id_and_name(term)
        term_id = term_id.replace('MOD:', '')
        parent_ids = [pid.replace('MOD:', '') for pid in _get_parent_ids(term)]

        if _is_obsolete(term):
            continue

        property_values = {}
        for val in term.get('xref', []):
            elems = val.split('"')
            if len(elems) < 2:
                continue

            k = elems[0].rstrip().replace(':', '')
            v = elems[1].strip()

            property_values.setdefault(k, []).append(v)

        delta_composition = property_values.get('DiffFormula', [None])
        delta_monoisotopic_mass = property_values.get('DiffMono', [None])
        delta_average_mass = property_values.get('DiffAvg', [None])

        if len(delta_composition) > 1:
            warnings.warn(
                f'[{DbType.PSI_MOD}] Multiple delta compositions for {term_id} {term_name} {delta_composition}')
            delta_composition = delta_composition[0]
        elif len(delta_composition) == 1:
            delta_composition = delta_composition[0]

        if len(delta_monoisotopic_mass) > 1:
            warnings.warn(
                f'[{DbType.PSI_MOD}] Multiple delta mono masses for {term_id} {term_name} {delta_monoisotopic_mass}')
            delta_monoisotopic_mass = delta_monoisotopic_mass[0]
        elif len(delta_monoisotopic_mass) == 1:
            delta_monoisotopic_mass = delta_monoisotopic_mass[0]

        if len(delta_average_mass) > 1:
            warnings.warn(
                f'[{DbType.PSI_MOD}] Multiple delta average masses for {term_id} {term_name} {delta_average_mass}')
            delta_average_mass = delta_average_mass[0]
        elif len(delta_average_mass) == 1:
            delta_average_mass = delta_average_mass[0]

        # check if delta_monoisotopic_mass is  'none'
        if delta_monoisotopic_mass == 'none':
            delta_monoisotopic_mass = None

        if delta_average_mass == 'none':
            delta_average_mass = None

        if delta_composition == 'none':
            delta_composition = None

        synonyms = term.get('synonym', [])
        if not isinstance(synonyms, List):
            synonyms = [synonyms]
        synonyms = [syn.split('"')[1] for syn in synonyms]

        delta_formula = None
        if isinstance(delta_composition, str):
            delta_composition = delta_composition.split(' ')

            elem_counter = Counter()
            for comp, num in zip(delta_composition[::2], delta_composition[1::2]):
                elem_counter[comp.replace('(', '').replace(')', '')] = int(num)

            delta_formula = write_chem_formula(elem_counter)

            try:
                _ = parse_chem_formula(delta_formula)
            except InvalidChemFormulaError as e:  # could be a glycan composition
                delta_formula = None
                warnings.warn(
                    f'[{DbType.PSI_MOD}] Error parsing {term_id} {term_name} {property_values.get("DiffFormula")}, {e}')

        if delta_monoisotopic_mass is not None:
            delta_monoisotopic_mass = float(delta_monoisotopic_mass)

        if delta_average_mass is not None:
            delta_average_mass = float(delta_average_mass)

        yield ModEntry(
            id=term_id,
            name=term_name,
            mono_mass=delta_monoisotopic_mass,
            avg_mass=delta_average_mass,
            composition=delta_formula,
            synonyms=synonyms,
            parents=parent_ids,
            entry_type=DbType.PSI_MOD
        )


def _get_resid_entries(terms: List[Dict[str, Any]]) -> List[ModEntry]:
    for term in terms:

        term_id, term_name = _get_id_and_name(term)
        term_id = term_id.replace('RESID:', '')
        parent_ids = [pid.replace('RESID:', '') for pid in _get_parent_ids(term)]

        if _is_obsolete(term):
            continue

        # grab RESID from PSI def
        definition = term.get('def', [None])

        if len(definition) > 1:
            warnings.warn(f'[{DbType.RESID}] Multiple definitions for {term_id} {term_name} {definition}')
            definition = definition[0]
        elif len(definition) == 1:
            definition = definition[0]

        regex_str = r'RESID:(\w+)'
        res_ids = []
        if definition:
            res_ids = re.findall(regex_str, definition)

        property_values = {}
        for val in term.get('xref', []):
            elems = val.split('"')
            if len(elems) < 2:
                continue

            k = elems[0].rstrip().replace(':', '')
            v = elems[1].strip()

            property_values.setdefault(k, []).append(v)

        delta_composition = property_values.get('DiffFormula', [None])
        delta_monoisotopic_mass = property_values.get('DiffMono', [None])
        delta_average_mass = property_values.get('DiffAvg', [None])

        if len(delta_composition) > 1:
            warnings.warn(f'[{DbType.RESID}] Multiple delta compositions for {term_id} {term_name} {delta_composition}')
            delta_composition = delta_composition[0]
        elif len(delta_composition) == 1:
            delta_composition = delta_composition[0]

        if len(delta_monoisotopic_mass) > 1:
            warnings.warn(
                f'[{DbType.RESID}] Multiple delta mono masses for {term_id} {term_name} {delta_monoisotopic_mass}')
            delta_monoisotopic_mass = delta_monoisotopic_mass[0]
        elif len(delta_monoisotopic_mass) == 1:
            delta_monoisotopic_mass = delta_monoisotopic_mass[0]

        if len(delta_average_mass) > 1:
            warnings.warn(
                f'[{DbType.RESID}] Multiple delta average masses for {term_id} {term_name} {delta_average_mass}')
            delta_average_mass = delta_average_mass[0]
        elif len(delta_average_mass) == 1:
            delta_average_mass = delta_average_mass[0]

        # check if delta_monoisotopic_mass is  'none'
        if delta_monoisotopic_mass == 'none':
            delta_monoisotopic_mass = None

        if delta_average_mass == 'none':
            delta_average_mass = None

        if delta_composition == 'none':
            delta_composition = None

        if delta_monoisotopic_mass is not None:
            delta_monoisotopic_mass = float(delta_monoisotopic_mass)

        if delta_average_mass is not None:
            delta_average_mass = float(delta_average_mass)

        synonyms = term.get('synonym', [])
        if not isinstance(synonyms, List):
            synonyms = [synonyms]
        synonyms = [syn.split('"')[1] for syn in synonyms]

        delta_formula = None
        if delta_composition:
            delta_composition = delta_composition.split(' ')

            elem_counter = Counter()
            for comp, num in zip(delta_composition[::2], delta_composition[1::2]):
                elem_counter[comp.replace('(', '').replace(')', '')] = int(num)

            delta_formula = write_chem_formula(elem_counter)

            try:
                _ = parse_chem_formula(delta_formula)
            except InvalidChemFormulaError as e:  # could be a glycan composition
                delta_formula = None
                warnings.warn(
                    f'[{DbType.RESID}] Error parsing {term_id} {term_name} {property_values.get("DiffFormula")}, {e}')

        if delta_monoisotopic_mass is not None:
            delta_monoisotopic_mass = float(delta_monoisotopic_mass)

        if delta_average_mass is not None:
            delta_average_mass = float(delta_average_mass)

        for res_id in res_ids:
            yield ModEntry(
                id=res_id,
                name=term_name,
                mono_mass=delta_monoisotopic_mass,
                avg_mass=delta_average_mass,
                composition=delta_formula,
                synonyms=synonyms,
                parents=parent_ids,
                entry_type=DbType.RESID
            )


def _get_xlmod_entries(terms: List[Dict[str, Any]]) -> List[ModEntry]:
    for term in terms:

        term_id, term_name = _get_id_and_name(term)
        term_id = term_id.replace('XLMOD:', '')
        parent_ids = [pid.replace('XLMOD:', '') for pid in _get_parent_ids(term)]

        if _is_obsolete(term):
            continue

        property_values = {}
        for val in term.get('property_value', []):
            elems = val.split('"')
            if len(elems) < 2:
                continue

            k = elems[0].rstrip().replace(':', '')
            v = elems[1].strip()

            property_values.setdefault(k, []).append(v)

        delta_composition = property_values.get('bridgeFormula', property_values.get('deadEndFormula', [None]))
        delta_monoisotopic_mass = property_values.get('monoIsotopicMass', [None])
        delta_average_mass = None

        if len(delta_composition) > 1:
            warnings.warn(f'[{DbType.XLMOD}] Multiple delta compositions for {term_id} {term_name} {delta_composition}')
            delta_composition = delta_composition[0]
        elif len(delta_composition) == 1:
            delta_composition = delta_composition[0]

        if len(delta_monoisotopic_mass) > 1:
            warnings.warn(
                f'[{DbType.XLMOD}] Multiple delta mono masses for {term_id} {term_name} {delta_monoisotopic_mass}')
            delta_monoisotopic_mass = delta_monoisotopic_mass[0]
        elif len(delta_monoisotopic_mass) == 1:
            delta_monoisotopic_mass = delta_monoisotopic_mass[0]

        synonyms = term.get('synonym', [])
        if not isinstance(synonyms, List):
            synonyms = [synonyms]
        synonyms = [syn.split('"')[1] for syn in synonyms]

        delta_formula = None
        if delta_composition:
            delta_composition = delta_composition.split(' ')
            elem_counter = Counter()
            for comp in delta_composition:
                if comp.startswith('-'):
                    comp = comp[1:]
                    comps = _parse_isotope_component(comp)
                    for elem, count in comps.items():
                        elem_counter[elem] -= count
                else:
                    comps = _parse_isotope_component(comp)
                    for elem, count in comps.items():
                        elem_counter[elem] += count

            delta_formula = write_chem_formula(elem_counter)

            try:
                _ = chem_mass(delta_formula, monoisotopic=True)
                _ = chem_mass(delta_formula, monoisotopic=False)
            except InvalidChemFormulaError as e:
                warnings.warn(f'[{DbType.XLMOD}] Error parsing {term_id} {term_name} {delta_formula}, {e}')
                delta_formula = None
                delta_average_mass = None

        if delta_monoisotopic_mass is not None:
            delta_monoisotopic_mass = float(delta_monoisotopic_mass)

        if delta_average_mass is not None:
            delta_average_mass = float(delta_average_mass)

        yield ModEntry(
            id=term_id,
            name=term_name,
            mono_mass=delta_monoisotopic_mass,
            avg_mass=delta_average_mass,
            composition=delta_formula,
            synonyms=synonyms,
            parents=parent_ids,
            entry_type=DbType.XLMOD
        )


def _get_gno_entries(terms: List[Dict[str, Any]]) -> List[ModEntry]:
    for term in terms:

        term_id, term_name = _get_id_and_name(term)
        term_id = term_id.replace('GNO:', '')
        parent_ids = [pid.replace('GNO:', '') for pid in _get_parent_ids(term)]

        if _is_obsolete(term):
            continue

        property_values = {}
        for val in term.get('property_value', []):
            try:
                elems = val.split('"')
                k = elems[0].rstrip()
                v = elems[1].strip()

                property_values.setdefault(k, []).append(v)

            except IndexError:
                continue

        synonyms = term.get('synonym', [])
        if not isinstance(synonyms, List):
            synonyms = [synonyms]
        synonyms = [syn.split('"')[1] for syn in synonyms]

        delta_formula = None
        delta_average_mass = None
        delta_monoisotopic_mass = None

        val = property_values.get('GNO:00000202', [])
        if len(val) > 1:
            warnings.warn(f'[{DbType.GNO}] Multiple delta compositions for {term_id} {term_name} {val}')
            val = val[0]
        elif len(val) == 1:
            val = val[0]
        else:
            val = None

        if val:
            composition = Counter()
            tokens = re.findall(r"([A-Za-z0-9]+)\((\d+)\)", val)
            for symbol, count in tokens:
                try:
                    comp = _glycan_comp(symbol)
                except ValueError as err:
                    warnings.warn(f'[{DbType.GNO}] Error parsing {term_id} {term_name} {val}, {err}')
                    composition = None
                    break

                for elem, c in comp.items():
                    composition[elem] += c * int(count)

            if composition is not None:
                delta_formula = write_chem_formula(composition)

                try:
                    delta_monoisotopic_mass_recalc = chem_mass(delta_formula, monoisotopic=True)
                    delta_average_mass_recalc = chem_mass(delta_formula, monoisotopic=False)
                except InvalidChemFormulaError as err:
                    warnings.warn(f'[{DbType.GNO}] Error parsing {term_id} {term_name} {delta_formula}, {err}')
                    delta_formula = None

        if delta_monoisotopic_mass is not None:
            delta_monoisotopic_mass = float(delta_monoisotopic_mass)

        if delta_average_mass is not None:
            delta_average_mass = float(delta_average_mass)

        yield ModEntry(
            id=term_id,
            name=term_name,
            mono_mass=delta_monoisotopic_mass,
            avg_mass=delta_average_mass,
            composition=delta_formula,
            synonyms=synonyms,
            parents=parent_ids,
            entry_type=DbType.GNO
        )


def _get_monosaccharide_entries(terms: List[Dict[str, Any]]) -> List[ModEntry]:
    for term in terms:

        term_id, term_name = _get_id_and_name(term)
        term_id = term_id.replace('MONO:', '')
        parent_ids = [pid.replace('MONO:', '') for pid in _get_parent_ids(term)]

        if _is_obsolete(term):
            continue

        property_values: Dict[str, List[str]] = {}
        for val in term.get('property_value', []):
            elems = val.split('"')
            k = elems[0].rstrip()
            v = elems[1].strip()

            property_values.setdefault(k, []).append(v)

        delta_formula = property_values.get('has_chemical_formula', [None])
        delta_monoisotopic_mass = property_values.get('has_monoisotopic_mass', [None])
        delta_average_mass = property_values.get('has_average_mass', [None])

        if len(delta_formula) > 1:
            warnings.warn(
                f'[{DbType.MONOSACCHARIDES}] Multiple delta compositions for {term_id} {term_name} {delta_formula}')
            delta_formula = delta_formula[0]
        elif len(delta_formula) == 1:
            delta_formula = delta_formula[0]

        if len(delta_monoisotopic_mass) > 1:
            warnings.warn(
                f'[{DbType.MONOSACCHARIDES}]Multiple delta mono masses for {term_id} {term_name} '
                f'{delta_monoisotopic_mass}')
            delta_monoisotopic_mass = delta_monoisotopic_mass[0]
        elif len(delta_monoisotopic_mass) == 1:
            delta_monoisotopic_mass = delta_monoisotopic_mass[0]

        if len(delta_average_mass) > 1:
            warnings.warn(
                f'[{DbType.MONOSACCHARIDES}] Multiple delta average masses for {term_id} {term_name} '
                f'{delta_average_mass}')
            delta_average_mass = delta_average_mass[0]
        elif len(delta_average_mass) == 1:
            delta_average_mass = delta_average_mass[0]

        synonyms = term.get('synonym', [])
        if not isinstance(synonyms, List):
            synonyms = [synonyms]
        synonyms = [syn.split('"')[1] for syn in synonyms]

        # try to parse formula
        delta_monoisotopic_mass_recalc = None
        delta_average_mass_recalc = None
        try:
            delta_monoisotopic_mass_recalc = chem_mass(delta_formula, monoisotopic=True)
            delta_average_mass_recalc = chem_mass(delta_formula, monoisotopic=False)
        except InvalidChemFormulaError as e:
            warnings.warn(f'[{DbType.MONOSACCHARIDES}]Error parsing {term_id} {term_name} {delta_formula}, {e}')
            delta_formula = None

        if delta_monoisotopic_mass is not None:
            delta_monoisotopic_mass = float(delta_monoisotopic_mass)
        else:
            delta_monoisotopic_mass = delta_monoisotopic_mass_recalc

        if delta_average_mass is not None:
            delta_average_mass = float(delta_average_mass)
        else:
            delta_average_mass = delta_average_mass_recalc

        yield ModEntry(
            id=term_id,
            name=term_name,
            mono_mass=delta_monoisotopic_mass,
            avg_mass=delta_average_mass,
            composition=delta_formula,
            synonyms=synonyms,
            parents=parent_ids,
            entry_type=DbType.MONOSACCHARIDES
        )


def get_entries(db_type: DbType, file_path: str) -> List[ModEntry]:
    """
    Get the entries from a database file.

    :param db_type: The type of database to get the entries from.
    :type db_type: DbType

    :return: The entries from the database file.
    :rtype: List[ModEntry]

    """

    with open(file_path, 'r') as f:
        data = _read_obo(f)

    if db_type == DbType.UNIMOD:
        return list(_get_unimod_entries(data))

    if db_type == DbType.PSI_MOD:
        return list(_get_psimod_entries(data))

    elif db_type == DbType.GNO:
        return list(_get_gno_entries(data))

    if db_type == DbType.RESID:
        return list(_get_resid_entries(data))

    if db_type == DbType.XLMOD:
        return list(_get_xlmod_entries(data))

    if db_type == DbType.MONOSACCHARIDES:
        return list(_get_monosaccharide_entries(data))

    raise ValueError(f"Invalid type: {db_type}")


class EntryDb:
    """
    A database of entries.
    """

    def __init__(self, entry_type: str, entries: List[ModEntry], use_synonyms: bool = False):
        self.entry_type = entry_type
        self.entries = entries
        self.id_map = {}
        self.name_map = {}
        self.synonym_map = {}
        self.mono_mass_map = []
        self.avg_mass_map = []
        self.names_sorted = []

        self.use_synonyms = use_synonyms
        self.setup(entries, use_synonyms)

    def __repr__(self) -> str:
        return f'<EntryDb: {self.entry_type}, Entries: {len(self)}, Use Synonyms: {self.use_synonyms}>'

    def __str__(self) -> str:
        return f'EntryDb: {self.entry_type}'

    def __len__(self) -> int:
        return len(self.id_map)

    def __getitem__(self, key: str) -> ModEntry:
        return self.get_entry_by_id(key)

    def __iter__(self) -> Iterator[ModEntry]:
        return iter(self.id_map.values())

    def __contains__(self, key: str) -> bool:
        return key in self.id_map

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, EntryDb):
            return False
        return self.entry_type == other.entry_type and self.id_map == other.id_map

    def __ne__(self, other: Any) -> bool:
        return not self.__eq__(other)

    def __hash__(self) -> int:
        return hash((self.entry_type, frozenset(self.id_map.items())))

    def setup(self, entries: List[ModEntry], synonyms: bool = False) -> None:
        self._setup_id_map(entries)
        self._setup_name_map(entries)
        self._setup_synonym_map(entries)
        self.setup_mono_list(entries)
        self._setup_avg_list(entries)
        self.use_synonyms = synonyms

        self.names_sorted = self._get_names_sorted()

    def contains_id(self, id: str) -> bool:
        return id in self.id_map

    def contains_name(self, name: str) -> bool:
        return name in self.name_map

    def contains_synonym(self, synonym: str) -> bool:
        return synonym in self.synonym_map

    def contains_mono_mass(self, mass: float, tol: float = 0.1) -> bool:
        return bool(self.get_entries_by_mono_mass(mass, tol))

    def contains_avg_mass(self, mass: float, tol: float = 0.1) -> bool:
        return bool(self.get_entries_by_avg_mass(mass, tol))

    def _add_entry_to_id_map(self, entry: ModEntry) -> None:
        if entry.id in self.id_map and self.id_map[entry.id] != entry:
            warnings.warn(f'[{self.entry_type}] Duplicate id: {entry.id}\n{entry}\n{self.id_map[entry.id]}')
        self.id_map[entry.id] = entry

    def _setup_id_map(self, entries: List[ModEntry]) -> None:
        for entry in entries:
            self._add_entry_to_id_map(entry)

    def _add_entry_to_name_map(self, entry: ModEntry) -> None:
        if entry.name in self.name_map and self.name_map[entry.name] != entry:
            warnings.warn(f'[{self.entry_type}] Duplicate name: {entry.name}\n{entry}\n{self.name_map[entry.name]}')
        self.name_map[entry.name] = entry

    def _setup_name_map(self, entries: List[ModEntry]) -> None:
        for entry in entries:
            self._add_entry_to_name_map(entry)

    def _add_entry_to_synonym_map(self, entry: ModEntry) -> None:

        if self.use_synonyms is False:
            return

        for syn in entry.synonyms:
            if syn in self.synonym_map and self.synonym_map[syn] != entry:
                warnings.warn(f'[{self.entry_type}] Duplicate synonym: {syn}\n{entry}\n{self.synonym_map[syn]}')
            self.synonym_map[syn] = entry

    def _setup_synonym_map(self, entries: List[ModEntry]) -> None:
        for entry in entries:
            self._add_entry_to_synonym_map(entry)

    def _add_entry_to_mono_list(self, entry: ModEntry) -> None:
        if entry.mono_mass is not None:
            self.mono_mass_map.append(entry)
            self.mono_mass_map.sort(key=lambda x: x.mono_mass)

    def setup_mono_list(self, entries: List[ModEntry]) -> None:
        valid_mono_entries = [entry for entry in entries if entry.mono_mass is not None]
        self.mono_mass_map = sorted(valid_mono_entries, key=lambda x: x.mono_mass)

    def _add_entry_to_avg_list(self, entry: ModEntry) -> None:
        if entry.avg_mass is not None:
            self.avg_mass_map.append(entry)
            self.avg_mass_map.sort(key=lambda x: x.avg_mass)

    def add_entry(self, entries: Union[ModEntry, List[ModEntry]]) -> None:
        """
        Add an entry or a list of entries to the database.
        """

        if isinstance(entries, ModEntry):
            entries = [entries]

        for entry in entries:
            self._add_entry_to_id_map(entry)
            self._add_entry_to_name_map(entry)
            self._add_entry_to_synonym_map(entry)

            if entry.mono_mass is not None:
                self.mono_mass_map.append(entry)

            if entry.avg_mass is not None:
                self.avg_mass_map.append(entry)

        self.mono_mass_map.sort(key=lambda x: x.mono_mass)
        self.avg_mass_map.sort(key=lambda x: x.avg_mass)
        self.names_sorted = self._get_names_sorted()

    def _setup_avg_list(self, entries: List[ModEntry]) -> None:
        valid_avg_entries = [entry for entry in entries if entry.avg_mass is not None]
        self.avg_mass_map = sorted(valid_avg_entries, key=lambda x: x.avg_mass)

    def get_entry_by_id(self, id: str) -> ModEntry:
        return self.id_map.get(id)

    def get_entry_by_name(self, name: str) -> ModEntry:
        return self.name_map.get(name)

    def get_entry_by_synonym(self, synonym: str) -> ModEntry:
        return self.synonym_map.get(synonym)

    def get_entries_by_mono_mass(self, mass: float, tol: float = 0.1) -> List[ModEntry]:
        return [entry for entry in self.mono_mass_map if abs(entry.mono_mass - mass) < tol]

    def get_entries_by_avg_mass(self, mass: float, tol: float = 0.1) -> List[ModEntry]:
        return [entry for entry in self.avg_mass_map if abs(entry.avg_mass - mass) < tol]

    def reset(self) -> None:
        self.id_map = {}
        self.name_map = {}
        self.synonym_map = {}
        self.mono_mass_map = []
        self.avg_mass_map = []
        self.names_sorted = []

    def save(self, file: str) -> None:
        with open(file, 'wb') as f:
            pickle.dump(self, f)

    @staticmethod
    def load(file: str) -> None:
        with open(file, 'rb') as f:
            return pickle.load(f)

    def reload_from_online(self, use_synonyms: bool = False) -> None:
        """
        Download obo file from the internet and reload the database
        """
        try:
            import requests
        except ImportError as err:
            raise ImportError("The requests module is required to download the obo file from the internet") from err

        file_url = _OBO_FILE_SOURCES.get(self.entry_type)
        if file_url is None:
            raise ValueError(f"Invalid type: {self.entry_type}")

        # Download file to a temporary location
        with requests.get(file_url, stream=True) as r:
            r.raise_for_status()
            with tempfile.NamedTemporaryFile(mode='w+b', delete=False) as tmp_file:
                for chunk in r.iter_content(chunk_size=8192):
                    tmp_file.write(chunk)
                tmp_file_path = tmp_file.name

        # Reload database from the downloaded file
        self.reload_from_file(tmp_file_path, use_synonyms)

    def reload_from_file(self, file: str, use_synonyms: bool = False) -> None:
        with open(file, 'r') as f:
            data = _read_obo(f)

        if self.entry_type == DbType.UNIMOD:
            entries = list(_get_unimod_entries(data))
        elif self.entry_type == DbType.PSI_MOD:
            entries = list(_get_psimod_entries(data))
        elif self.entry_type == DbType.GNO:
            entries = list(_get_gno_entries(data))
        elif self.entry_type == DbType.RESID:
            entries = list(_get_resid_entries(data))
        elif self.entry_type == DbType.XLMOD:
            entries = list(_get_xlmod_entries(data))
        elif self.entry_type == DbType.MONOSACCHARIDES:
            entries = list(_get_monosaccharide_entries(data))
        else:
            raise ValueError(f"Invalid type: {self.entry_type}")

        old_maps = (
            self.id_map, self.name_map, self.synonym_map, self.mono_mass_map, self.avg_mass_map, self.use_synonyms,
            self.names_sorted, self.entries
        )

        try:
            self.reset()
            self.setup(entries, use_synonyms)
        except Exception as err:
            self.id_map, self.name_map, self.synonym_map, self.mono_mass_map, self.avg_mass_map, self.use_synonyms, \
                self.names_sorted, self.entries = old_maps
            raise err

        # print(f'Reloaded {self.entry_type} database from {file}')
        # print('Entries:', len(self))

    def _get_names_sorted(self) -> List[str]:
        synonyms = set(self.synonym_map.keys())
        names = set(self.name_map.keys())
        return sorted(list(synonyms.union(names)), key=lambda x: len(x), reverse=True)


def count_invalid_entries(entries: List[ModEntry]) -> Tuple[int, int, int]:
    """
    Count the number of invalid entries in a list of ModEntry objects.
    """
    none_mono, none_avg, none_comp = 0, 0, 0

    for entry in entries:
        if entry.mono_mass is None:
            none_mono += 1
        if entry.avg_mass is None:
            none_avg += 1
        if entry.composition is None:
            none_comp += 1

    return none_mono, none_avg, none_comp


# import pkgutil
# data = pkgutil.get_data(__name__, "../data/monosaccharides_updated.obo")

_dir_name = os.path.dirname(__file__)
_obo_path = os.path.join(_dir_name, "..", "data")
MONOSACCHARIDES_DB = EntryDb(DbType.MONOSACCHARIDES, [], True)
UNIMOD_DB = EntryDb(DbType.UNIMOD, [], False)
PSI_MOD_DB = EntryDb(DbType.PSI_MOD, [], False)
XLMOD_DB = EntryDb(DbType.XLMOD, [], False)
GNO_DB = EntryDb(DbType.GNO, [], False)
RESID_DB = EntryDb(DbType.RESID, [], False)

MONOSACCHARIDES_DB.reload_from_file(os.path.join(_obo_path, "monosaccharides_updated.obo"))
UNIMOD_DB.reload_from_file(os.path.join(_obo_path, "unimod.obo"))
PSI_MOD_DB.reload_from_file(os.path.join(_obo_path, "psi-mod.obo"))
XLMOD_DB.reload_from_file(os.path.join(_obo_path, "xlmod.obo"))
#GNO_DB.reload_from_file(os.path.join(_obo_path, "gno.obo"))
#RESID_DB.reload_from_file(os.path.join(_obo_path, "psi-mod.obo"))


def reload_all_databases_from_online() -> None:
    """
    Reload all databases from the online sources
    """
    MONOSACCHARIDES_DB.reload_from_online()
    UNIMOD_DB.reload_from_online()
    PSI_MOD_DB.reload_from_online()
    XLMOD_DB.reload_from_online()
    GNO_DB.reload_from_online()
    RESID_DB.reload_from_online()


def reload_all_databases() -> None:
    """
    Reload all databases from the local files
    """
    reset_all_databases()
    MONOSACCHARIDES_DB.reload_from_file(os.path.join(_obo_path, "monosaccharides_updated.obo"))
    UNIMOD_DB.reload_from_file(os.path.join(_obo_path, "unimod.obo"))
    PSI_MOD_DB.reload_from_file(os.path.join(_obo_path, "psi-mod.obo"))
    XLMOD_DB.reload_from_file(os.path.join(_obo_path, "xlmod.obo"))
    GNO_DB.reload_from_file(os.path.join(_obo_path, "gno.obo"))
    RESID_DB.reload_from_file(os.path.join(_obo_path, "psi-mod.obo"))


def reset_all_databases() -> None:
    """
    Reset all databases to empty
    """
    MONOSACCHARIDES_DB.reset()
    UNIMOD_DB.reset()
    PSI_MOD_DB.reset()
    GNO_DB.reset()
    RESID_DB.reset()
    XLMOD_DB.reset()
