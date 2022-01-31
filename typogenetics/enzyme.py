import random
from typing import List

from typogenetics import amino_acids, strand

BINDING_PREFERENCES = {
    0: "A",
    1: "C",
    2: "G",
    3: "T",
}  # Enzyme shape expressed as net left folds


class Enzyme:
    def __init__(self, code: List[amino_acids.AminoAcidType]):
        self.amino_acids = code
        self.next_index = 0
        self.unit = None
        self.copy_mode = False

    @property
    def binding_preference(self):
        net_left_folds = 0
        if len(self.amino_acids) > 2:
            inner_amino_acids = self.amino_acids[1:-2]
            for amino_acid in inner_amino_acids:
                net_left_folds += amino_acids.ENZYME_LEFT_FOLDS[amino_acid]
            net_left_folds %= 4
        return BINDING_PREFERENCES[net_left_folds]

    def bind_to_strand(self, strand_member: strand.Unit):
        binding_preference = self.binding_preference
        candidate_units = []

        # Scan left for matching bases
        cursor = strand_member
        while True:
            if cursor is None:
                break
            if cursor.base == binding_preference:
                candidate_units.append(cursor)
            cursor = cursor.left

        # Scan right for matching bases
        cursor = strand_member.right
        while True:
            if cursor is None:
                break
            if cursor.base == binding_preference:
                candidate_units.append(cursor)
            cursor = cursor.right

        # Bind to a randomly-chosen candidate, if any
        if candidate_units:
            self.unit = random.choice(candidate_units)

    def operate(self):
        while self.is_active:
            self.step()

    @property
    def is_active(self):
        return self.unit is not None \
               and self.amino_acids \
               and self.next_index < len(self.amino_acids)

    def step(self):
        assert self.is_active
        next_amino_acid = self.amino_acids[self.next_index]
        next_operation = getattr(self, "_" + next_amino_acid)
        next_operation()
        self.next_index += 1

    def _cut(self):
        right_unit = self.unit.right
        pair_unit = self.unit.pair

        # For right, remove references each way between units
        if right_unit is not None:
            self.unit.right = None
            right_unit.left = None

        # If paired, cut there also (but on its left)
        if pair_unit is not None:
            pair_left = pair_unit.left
            if pair_left is not None:
                pair_left = None
                pair_unit.left = None

    def _del(self):
        right_unit = self.unit.right
        left_unit = self.unit.left
        pair_unit = self.unit.pair

        # For each, remove references each way between units
        if right_unit is not None:
            self.unit.right = None
            right_unit.left = None
        if left_unit is not None:
            self.unit.left = None
            left_unit.right = None
        if pair_unit is not None:
            self.unit.pair = None
            pair_unit.pair = None

        # At the end, bind to the next unit to the right
        self.unit = right_unit

    def _swi(self):
        self.unit = self.unit.pair

    def _mvr(self):
        self.unit = self.unit.right
        self._maybe_copy()

    def _mvl(self):
        self.unit = self.unit.left
        self._maybe_copy()

    def _cop(self):
        self.copy_mode = True

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
        while self.unit is not None \
                and self.unit.base not in strand.PYRIMIDINES:
            self._mvr()

    def _rpu(self):
        self._mvr()
        while self.unit is not None \
                and self.unit.base not in strand.PURINES:
            self._mvr()

    def _lpy(self):
        self._mvl()
        while self.unit is not None \
                and self.unit.base not in strand.PYRIMIDINES:
            self._mvl()

    def _lpu(self):
        self._mvl()
        while self.unit is not None \
                and self.unit.base not in strand.PURINES:
            self._mvl()

    def _maybe_copy(self):
        if self.copy_mode and self.unit.pair is None:
            new_pair_unit = strand.Unit(base=strand.COMPLIMENTS[self.unit.base])

            # Attach the new unit across from the current one
            self.unit.pair = new_pair_unit
            new_pair_unit.pair = self.unit.pair

            # Check if new unit has a neighbors on either side
            if self.unit.right is not None \
                    and self.unit.right.pair is not None:
                new_pair_unit.left = self.unit.right.pair
                self.unit.right.pair.right = new_pair_unit
            if self.unit.left is not None \
                    and self.unit.left.pair is not None:
                new_pair_unit.right = self.unit.left.pair
                self.unit.left.pair.left = new_pair_unit

    def _insert_base(self, base: strand.BaseType):
        # Note the unit to the right, then cut the strand
        right_unit = self.unit.right
        self._cut()

        # Attach a new base to the right of the current unit
        new_unit = strand.Unit(base=base)
        self.unit.right = new_unit
        new_unit.left = self.unit

        # If needed, attach the previous right unit after the new one
        if right_unit is not None:
            new_unit.right = right_unit
            right_unit.left = new_unit

        self._maybe_copy()
