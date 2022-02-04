import collections
import typing

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
Unit = collections.namedtuple("Unit", ["facing", "opposite"])


class Strand(typing.List[Unit]):

    @classmethod
    def make_from_code(cls, code: str):
        return cls(Unit(base, None) for base in code)

    def __str__(self):
        if any(unit.opposite is not None for unit in self):
            return (
                "".join(unit.opposite or " " for unit in self)
                + "\n"
                + "".join(unit.facing or " " for unit in self)
            )
        return "".join(unit.facing or " " for unit in self)

    def unzip(self):
        facing_string = "".join(unit.facing or " " for unit in self)
        for code in facing_string.split():
            yield Strand.make_from_code(code)
        opposite_string = "".join(unit.opposite or " " for unit in reversed(self))
        for code in opposite_string.split():
            yield Strand.make_from_code(code)
