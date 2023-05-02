from typing import Union


class RefString:
    """A string-like class that stores a memory view of the string and slices by reference"""
    def __init__(self, value: Union[str, memoryview]):
        """Initializes a RefString object.

        Args:
            value (Union[str, memoryview]): The value to initialize the RefString with. It must be either a memory view object or a string.

        Raises:
            ValueError: If the input data is not a memory view object or a string.
        """
        if isinstance(value, memoryview):
            self._value = value
        elif isinstance(value, str):
            self._value = memoryview(value.encode('ascii'))
        else:
            raise ValueError("Input data must be either a memory view object or a string.")

    def __repr__(self):
        return f"RefString({self._value.tobytes().decode()})"

    def __str__(self):
        return self._value.tobytes().decode()

    def __getitem__(self, index):
        if isinstance(index, int):
            return RefString(self._value[index:index + 1])
        elif isinstance(index, slice):
            return RefString(self._value[index])

    def __len__(self):
        return len(self._value)

    def __eq__(self, other):
        if isinstance(other, RefString):
            return self._value == other._value
        else:
            return str(self) == str(other)

    def __hash__(self):
        return hash(str(self))

    def startswith(self, prefix):
        return str(self).startswith(str(prefix))

    def endswith(self, suffix):
        return str(self).endswith(str(suffix))

    def find(self, sub, start=0, end=None):
        end = end or len(self._value)
        return self._value.tobytes().find(sub.encode(), start, end)

    def rfind(self, sub, start=0, end=None):
        return str(self).rfind(sub, start, end)

    def count(self, sub):
        return str(self).count(sub)

    def split(self, sep=None, maxsplit=-1):
        return [RefString(s) for s in str(self).split(sep, maxsplit)]

    def replace(self, old, new, count=-1):
        return RefString(str(self).replace(old, new, count))

    def join(self, iterable):
        return RefString(str(self).join(iterable))

    def lower(self):
        return RefString(str(self).lower())

    def upper(self):
        return RefString(str(self).upper())

    def strip(self, chars=None):
        return RefString(str(self).strip(chars))

    def lstrip(self, chars=None):
        return RefString(str(self).lstrip(chars))

    def rstrip(self, chars=None):
        return RefString(str(self).rstrip(chars))

    def ljust(self, width, fillchar=' '):
        return RefString(str(self).ljust(width, fillchar))

    def rjust(self, width, fillchar=' '):
        return RefString(str(self).rjust(width, fillchar))

    def center(self, width, fillchar=' '):
        return RefString(str(self).center(width, fillchar))
