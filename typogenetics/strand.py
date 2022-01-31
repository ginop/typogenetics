import typing
from dataclasses import dataclass
from typing import Optional

BaseType = str
PURINES = {"A", "G"}
PYRIMIDINES = {"C", "T"}
BASES = PURINES | PYRIMIDINES
COMPLIMENTS = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
}
INVERSIONS = {
    "A": "\u2200",
    "T": "\u2534",
    "G": "\u05E4",
    "C": "\u0186",
}


@dataclass
class Unit:
    base: BaseType
    right: Optional["Unit"] = None
    left: Optional["Unit"] = None
    pair: Optional["Unit"] = None

    @property
    def left_most(self):
        return self if self.left is None else self.left.left_most

    @property
    def strand(self):
        strand = []
        unit = self.left_most
        while unit is not None:
            strand.append(unit)
            unit = unit.right
        return strand

    @property
    def strand_str(self):
        return "".join(unit.base for unit in self.strand)


class Strand(list):
    @classmethod
    def make_from_code(cls, code: BaseType):
        units = [Unit(base=base) for base in code]
        for ii in range(1, len(units)):
            units[ii-1].right = units[ii]
            units[ii].left = units[ii-1]
        return cls(units)
