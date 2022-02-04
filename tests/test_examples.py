from collections import Counter

from typogenetics.enzyme import Enzyme
from typogenetics.strand import Strand


def test_geb_1():
    example_strand = Strand.make_from_code("ACA")
    example_enzyme = Enzyme(["del", "mvr", "int"])

    produced_strands = Counter(str(strand) for strand in
                               example_enzyme.operate(example_strand, 0))

    expected_strands = Counter(["CAT"])
    assert produced_strands == expected_strands


def test_geb_2():
    example_strand = Strand.make_from_code("ACA")
    example_enzyme = Enzyme(["del", "mvr", "int"])

    produced_strands = Counter(str(strand) for strand in
                               example_enzyme.operate(example_strand, 2))

    expected_strands = Counter(["AC"])
    assert produced_strands == expected_strands


def test_geb_3():
    example_strand = Strand.make_from_code("CAAAGAGAATCCTCTTTGAT")
    example_enzyme = Enzyme(["rpy", "cop", "rpu", "cut"])

    produced_strands = Counter(str(strand) for strand in
                               example_enzyme.operate(example_strand, 2))

    expected_strands = Counter([
        "CAAAGAGAATCCTCTTTG",
        "AT",
        "CAAAGAGGA",
    ])
    assert produced_strands == expected_strands


def test_geb_4():
    example_genes = "TAGATCCAGTCCACATCGA"
    coded_enzymes = list(enzyme.amino_acids for enzyme in
                         Enzyme.make_from_genes(example_genes))
    expexted_enzymes = [["rpy", "ina", "rpu",
                         "mvr", "int", "mvl",
                         "cut", "swi", "cop"]]
    assert coded_enzymes == expexted_enzymes


def test_geb_5():
    example_enzyme = Enzyme(["rpy", "ina", "rpu",
                             "mvr", "int", "mvl",
                             "cut", "swi", "cop"])
    assert example_enzyme.binding_preference == "C"
