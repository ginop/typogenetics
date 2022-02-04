from collections import Counter

from typogenetics.enzyme import Enzyme
from typogenetics.strand import Strand


def test_geb_1():
    example_strand = Strand.make_from_code("ACA")
    example_enzyme = Enzyme(["del", "mvr", "int"])

    example_enzyme.strand = example_strand
    example_enzyme.unit_index = 0
    produced_strands = Counter(str(strand) for strand in example_enzyme.operate())

    expected_strands = Counter(["CAT"])
    assert produced_strands == expected_strands


def test_geb_2():
    example_strand = Strand.make_from_code("ACA")
    example_enzyme = Enzyme(["del", "mvr", "int"])

    example_enzyme.strand = example_strand
    example_enzyme.unit_index = 2
    produced_strands = Counter(str(strand) for strand in example_enzyme.operate())

    expected_strands = Counter(["AC"])
    assert produced_strands == expected_strands


def test_geb_3():
    example_strand = Strand.make_from_code("CAAAGAGAATCCTCTTTGAT")
    example_enzyme = Enzyme(["rpy", "cop", "rpu", "cut"])

    example_enzyme.strand = example_strand
    example_enzyme.unit_index = 2
    produced_strands = Counter(str(strand) for strand in example_enzyme.operate())

    expected_strands = Counter([
        "CAAAGAGAATCCTCTTTG",
        "AT",
        "CAAAGAGGA",
    ])
    assert produced_strands == expected_strands
