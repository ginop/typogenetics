import collections
import itertools
import random

from tqdm import tqdm

from typogenetics.strand import Strand, BASES
from typogenetics.enzyme import Enzyme


class Population(collections.Counter):
    def sample(self):
        random_index = random.randrange(self.total())
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


def population_simulation(log_file=None):
    strands = Population()
    enzymes = Population()

    for _ in range(1):
        original_strand = "".join(random.choices(tuple(BASES), k=random.randint(50, 101)))
        strands[original_strand] += 1
        for enzyme in Enzyme.make_from_genes(original_strand):
            enzymes[str(enzyme)] += 1

    for n_gen in tqdm(itertools.count()):
        # print(f"Gen {n_gen}: \n"
        #       f"{strands.total} strands in {len(strands)} species and \n"
        #       f"{enzymes.total} enzymes in {len(enzymes)} species")

        # If no enzymes, create more
        if enzymes.total() == 0:
            # print("No enzymes available")
            strand = strands.sample()
            for enzyme in Enzyme.make_from_genes(strand):
                enzymes[str(enzyme)] += 1
            strands[strand] += 1
            continue

        # Get one random enzyme and one strand
        enzyme = Enzyme(enzymes.sample().split(" - "))
        strand = Strand.make_from_code(strands.sample())
        # print(strand, enzyme, sep="\n")

        # Find preferred binding site
        binding_options = [n for n, unit in enumerate(strand) if unit.facing == enzyme.binding_preference]
        if not binding_options:
            # print("No binding site found")
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

        # print(strands.most_common(5))

        if log_file is not None:
            log_file.write(
                f"(total: {strands.total()}), " +
                ", ".join(f"({species}: {count})" for species, count in strands.most_common(10)) +
                "\n"
            )


if __name__ == '__main__':
    with open("population_log.csv", "w") as file:
        population_simulation(file)
