import typing
from typing import List

from amino_acids import AminoAcidType, ENZYME_LEFT_FOLDS
from strand import Strand, PYRIMIDINES, PURINES, COMPLIMENTS, BaseType, Unit

BINDING_PREFERENCES = {
    0: "A",
    1: "C",
    2: "G",
    3: "T",
}  # Enzyme shape expressed as net left folds


class Enzyme:
    def __init__(self, code: List[AminoAcidType]):
        self.amino_acids: List[AminoAcidType] = code
        self.next_index: int = 0
        self.copy_mode: bool = False

        self.strand: typing.Optional[Strand] = None
        self.unit_index: int = 0

    @property
    def binding_preference(self):
        net_left_folds = 0
        if len(self.amino_acids) > 2:
            inner_amino_acids = self.amino_acids[1:-2]
            for amino_acid in inner_amino_acids:
                net_left_folds += ENZYME_LEFT_FOLDS[amino_acid]
            net_left_folds %= 4
        return BINDING_PREFERENCES[net_left_folds]

    @property
    def is_active_and_attached(self):
        return self.next_index < len(self.amino_acids) \
               and self.strand is not None \
               and self.unit_index < len(self.strand) \
               and self.strand[self.unit_index].facing is not None

    def operate(self) -> typing.Iterable[Strand]:
        while self.is_active_and_attached:
            next_amino_acid = self.amino_acids[self.next_index]
            next_operation = getattr(self, "_" + next_amino_acid)
            new_strand_group = next_operation()
            if new_strand_group is not None:
                yield from new_strand_group.unzip()
            self.next_index += 1
        yield from self.strand.unzip()

    def _cut(self):
        cut_index = self.unit_index + 1
        removed_strand = Strand(self.strand[cut_index:]) if len(self.strand) > cut_index else None
        self.strand = Strand(self.strand[:cut_index])
        return removed_strand

    def _del(self):
        self.strand.pop(self.unit_index)

    def _swi(self):
        self.strand = Strand(Unit(unit.opposite, unit.facing) for unit in reversed(self.strand))
        self.unit_index = len(self.strand) - self.unit_index - 1

    def _mvr(self):
        self.unit_index += 1
        self._maybe_copy()

    def _mvl(self):
        self.unit_index -= 1
        self._maybe_copy()

    def _cop(self):
        self.copy_mode = True
        self._maybe_copy()

    def _off(self):
        self.copy_mode = False

    def _ina(self):
        self._insert_base("A")

    def _inc(self):
        self._insert_base("C")

    def _ing(self):
        self._insert_base("G")

    def _int(self):
        self._insert_base("T")

    def _rpy(self):
        self._mvr()
        while self.unit_index < len(self.strand) \
                and self.strand[self.unit_index].facing is not None \
                and self.strand[self.unit_index].facing not in PYRIMIDINES:
            self._mvr()

    def _rpu(self):
        self._mvr()
        while self.unit_index < len(self.strand) \
                and self.strand[self.unit_index].facing is not None \
                and self.strand[self.unit_index].facing not in PURINES:
            self._mvr()

    def _lpy(self):
        self._mvl()
        while self.unit_index < len(self.strand) \
                and self.strand[self.unit_index].facing is not None \
                and self.strand[self.unit_index].facing not in PYRIMIDINES:
            self._mvl()

    def _lpu(self):
        self._mvl()
        while self.unit_index < len(self.strand) \
                and self.strand[self.unit_index].facing is not None \
                and self.strand[self.unit_index].facing not in PURINES:
            self._mvl()

    def _maybe_copy(self):
        if self.copy_mode \
                and self.unit_index < len(self.strand) \
                and self.strand[self.unit_index].facing is not None \
                and self.strand[self.unit_index].opposite is None:
            self.strand[self.unit_index] = Unit(
                self.strand[self.unit_index].facing,
                COMPLIMENTS[self.strand[self.unit_index].facing],
            )

    def _insert_base(self, base: BaseType):
        insert_index = self.unit_index + 1
        self.strand.insert(insert_index, Unit(base, None))
        self._mvr()
