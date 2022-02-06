import random

from typogenetics.strand import Strand, BASES
from typogenetics.enzyme import Enzyme


def population_simulation():
    strand = Strand.make_from_code("".join(random.choices(tuple(BASES), k=random.randint(50, 101))))
    strands = [strand]
    enzymes = list(Enzyme.make_from_genes(str(strand)))

    while True:
        print(f"{len(strands)} strands and {len(enzymes)} enzymes")

        # Get one random enzyme and one strand
        enzyme = enzymes.pop(random.randrange(len(enzymes)))
        strand = strands.pop(random.randrange(len(strands)))

        print(strand, enzyme, sep="\n")

        # Process the strand
        new_strands = list(enzyme.operate(strand, 0))

        # Add new strands and enzymes to population
        strands.extend(new_strands)
        for strand in new_strands:
            enzymes.extend(list(Enzyme.make_from_genes(str(strand))))


if __name__ == '__main__':
    population_simulation()
