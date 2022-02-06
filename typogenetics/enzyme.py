import typing
from typing import List

from amino_acids import AminoAcidType, ENZYME_LEFT_FOLDS, DUPLETS
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

        self._strand: typing.Optional[Strand] = None
        self._unit_index: int = 0

    def __str__(self):
        return " - ".join(self.amino_acids)

    @classmethod
    def make_from_genes(cls, genes: BaseType):
        duplets = [genes[i-2:i] for i in range(2, len(genes)+1, 2)]
        amino_acids = []
        for duplet in duplets:
            if duplet == "AA":
                yield cls(amino_acids)
                amino_acids = []
            else:
                amino_acids.append(DUPLETS[duplet])
        yield cls(amino_acids)

    @property
    def binding_preference(self):
        net_left_folds = 0
        if len(self.amino_acids) > 2:
            inner_amino_acids = self.amino_acids[1:-1]
            for amino_acid in inner_amino_acids:
                net_left_folds += ENZYME_LEFT_FOLDS[amino_acid]
            net_left_folds %= 4
        return BINDING_PREFERENCES[net_left_folds]

    @property
    def is_active_and_attached(self):
        return self.next_index < len(self.amino_acids) \
               and self._strand is not None \
               and 0 <= self._unit_index < len(self._strand) \
               and self._strand[self._unit_index].facing is not None

    def operate(self, strand: Strand, starting_unit: int) -> typing.Iterable[Strand]:
        self._strand = strand
        self._unit_index = starting_unit

        while self.is_active_and_attached:
            next_amino_acid = self.amino_acids[self.next_index]
            next_operation = getattr(self, "_" + next_amino_acid)
            new_strand_group = next_operation()
            if new_strand_group is not None:
                yield from new_strand_group.unzip()
            self.next_index += 1
        yield from self._strand.unzip()

        self._strand = None
        self._unit_index = 0

    def _cut(self):
        cut_index = self._unit_index + 1
        removed_strand = Strand(self._strand[cut_index:]) if len(self._strand) > cut_index else None
        self._strand = Strand(self._strand[:cut_index])
        return removed_strand

    def _del(self):
        self._strand.pop(self._unit_index)

    def _swi(self):
        self._strand = Strand(Unit(unit.opposite, unit.facing) for unit in reversed(self._strand))
        self._unit_index = len(self._strand) - self._unit_index - 1

    def _mvr(self):
        self._unit_index += 1
        self._maybe_copy()

    def _mvl(self):
        self._unit_index -= 1
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
        while self._unit_index < len(self._strand) \
                and self._strand[self._unit_index].facing is not None \
                and self._strand[self._unit_index].facing not in PYRIMIDINES:
            self._mvr()

    def _rpu(self):
        self._mvr()
        while self._unit_index < len(self._strand) \
                and self._strand[self._unit_index].facing is not None \
                and self._strand[self._unit_index].facing not in PURINES:
            self._mvr()

    def _lpy(self):
        self._mvl()
        while 0 <= self._unit_index < len(self._strand) \
                and self._strand[self._unit_index].facing is not None \
                and self._strand[self._unit_index].facing not in PYRIMIDINES:
            self._mvl()

    def _lpu(self):
        self._mvl()
        while 0 <= self._unit_index < len(self._strand) \
                and self._strand[self._unit_index].facing is not None \
                and self._strand[self._unit_index].facing not in PURINES:
            self._mvl()

    def _maybe_copy(self):
        if self.copy_mode \
                and self._unit_index < len(self._strand) \
                and self._strand[self._unit_index].facing is not None \
                and self._strand[self._unit_index].opposite is None:
            self._strand[self._unit_index] = Unit(
                self._strand[self._unit_index].facing,
                COMPLIMENTS[self._strand[self._unit_index].facing],
            )

    def _insert_base(self, base: BaseType):
        insert_index = self._unit_index + 1
        self._strand.insert(insert_index, Unit(base, None))
        self._mvr()
