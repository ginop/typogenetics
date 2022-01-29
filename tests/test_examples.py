from typogenetics import strand, enzyme


def test_geb_1():
    example_strand = strand.Strand.make_from_code("ACA")
    example_enzyme = enzyme.Enzyme(["del", "mvr", "int"])

    example_enzyme.unit = example_strand[0]
    example_enzyme.operate()

    assert example_strand[1].strand_str == "CAT"


def test_geb_2():
    example_strand = strand.Strand.make_from_code("ACA")
    example_enzyme = enzyme.Enzyme(["del", "mvr", "int"])

    example_enzyme.unit = example_strand[2]
    example_enzyme.operate()

    assert example_strand[1].strand_str == "AC"


def test_geb_3():
    example_strand = strand.Strand.make_from_code("CAAAGAGAATCCTCTTTGAT")
    example_enzyme = enzyme.Enzyme(["rpy", "cop", "rpu", "cut"])

    example_enzyme.unit = example_strand[2]
    example_enzyme.operate()

    # TODO
    assert example_strand[1].strand_str == "CAAAGAGAATCCTCTTTG"
    assert example_strand[-1].strand_str == "AT"
