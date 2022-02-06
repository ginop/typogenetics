import collections
import itertools
import random

from typogenetics.strand import Strand, BASES
from typogenetics.enzyme import Enzyme


class Population(collections.Counter):
    def sample(self):
        assert self.total > 0
        random_index = random.randrange(self.total)
        for species, count in self.items():
            random_index -= count
            if random_index <= 0:
                if count == 1:
                    del self[species]
                else:
                    self[species] -= 1
                return species
        else:
            raise AssertionError

    @property
    def total(self):
        return sum(self.values())


def population_simulation():
    original_strand = "".join(random.choices(tuple(BASES), k=random.randint(50, 101)))
    strands = Population([original_strand])
    enzymes = Population(str(enzyme) for enzyme in Enzyme.make_from_genes(original_strand))

    for n_gen in itertools.count():
        print(f"Gen {n_gen}: \n"
              f"{strands.total} strands in {len(strands)} species and \n"
              f"{enzymes.total} enzymes in {len(enzymes)} species")

        # If no enzymes, create more
        if enzymes.total == 0:
            print("No enzymes available")
            strand = strands.sample()
            for enzyme in Enzyme.make_from_genes(strand):
                enzymes[str(enzyme)] += 1
            strands[strand] += 1
            continue

        # Get one random enzyme and one strand
        enzyme = Enzyme(enzymes.sample().split(" - "))
        strand = Strand.make_from_code(strands.sample())
        print(strand, enzyme, sep="\n")

        # Find preferred binding site
        binding_options = [n for n, unit in enumerate(strand) if unit.facing == enzyme.binding_preference]
        if not binding_options:
            print("No binding site found")
            strands[str(strand)] += 1
            enzymes[str(enzyme)] += 1
            continue

        # Process the strand
        new_strands = list(enzyme.operate(strand, random.choice(binding_options)))

        # Add new strands and enzymes to population
        for strand in new_strands:
            strands[str(strand)] += 1
            for enzyme in Enzyme.make_from_genes(str(strand)):
                enzymes[str(enzyme)] += 1

        print(strands.most_common(5))


if __name__ == '__main__':
    population_simulation()
