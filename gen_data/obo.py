import dataclasses
from typing import Dict, Any, List

from typing.io import IO


def read_obo(file: IO[str]) -> List[Dict[str, Any]]:
    file.seek(0)

    info = {}
    l = []

    skip = False

    d = None
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
                l.append(d)
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
        l.append(d)

    return l


@dataclasses.dataclass
class OboTerm:
    id: str
    name: str
    synonyms: List[str]
    data: Dict[str, Any]
    children: List['OboTerm']


def map_obo(file: IO[str]) -> Dict[str, Dict[str, Any]]:
    file.seek(0)

    info = {}
    l = {}

    skip = False

    d = None
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
                l[d['id'][0]] = d
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
        l[d['id'][0]] = d

    return l





