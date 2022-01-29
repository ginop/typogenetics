from typogenetics import strand, amino_acids


def test_bases():
    assert len(strand.bases) == 4
    for base in strand.bases:
        assert base in "ACGT"


def test_amino_acids():
    assert len(amino_acids.amino_acids) == 15


def test_duplets():
    assert len(amino_acids.duplets) == 15
    for duplet in amino_acids.duplets:
        assert len(duplet) == 2
        for base in duplet:
            assert base in strand.bases


def test_enzyme_folds():
    assert len(amino_acids.enzyme_folds) == 15
    for fold in amino_acids.enzyme_folds.values():
        assert fold in "srl"
