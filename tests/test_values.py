from typogenetics import strand, amino_acids


def test_bases():
    assert len(strand.BASES) == 4
    for base in strand.BASES:
        assert base in "ACGT"


def test_amino_acids():
    assert len(amino_acids.AMINO_ACIDS) == 15


def test_duplets():
    assert len(amino_acids.DUPLETS) == 15
    for duplet in amino_acids.DUPLETS:
        assert len(duplet) == 2
        for base in duplet:
            assert base in strand.BASES


def test_enzyme_folds():
    assert len(amino_acids.ENZYME_FOLDS) == 15
    for fold in amino_acids.ENZYME_FOLDS.values():
        assert fold in "srl"
