AminoAcidType = str
AMINO_ACIDS = {
    "cut",  # cut strand(s)
    "del",  # delete a base from strand
    "swi",  # switch enzyme to other strand
    "mvr",  # move one unit to the right
    "mvl",  # move one unit to the left
    "cop",  # turn on Copy mode
    "off",  # turn off Copy mode
    "ina",  # insert A to the right of this unit
    "inc",  # insert C to the right of this unit
    "ing",  # insert G to the right of this unit
    "int",  # insert T to the right of this unit
    "rpy",  # search for the nearest pyrimidine to the right
    "rpu",  # search for the nearest purine to the right
    "lpy",  # search for the nearest pyrimidine to the left
    "lpu",  # search for the nearest purine to the left
}

DupletType = str
DUPLETS = {
    "AC": "cut",
    "AG": "del",
    "AT": "swi",
    "CA": "mvr",
    "CC": "mvl",
    "CG": "cop",
    "CT": "off",
    "GA": "ina",
    "GC": "inc",
    "GG": "ing",
    "GT": "int",
    "TA": "rpy",
    "TC": "rpu",
    "TG": "lpy",
    "TT": "lpu",
}

ENZYME_FOLDS = {
    "cut": "s",
    "del": "s",
    "swi": "r",
    "mvr": "s",
    "mvl": "s",
    "cop": "r",
    "off": "l",
    "ina": "s",
    "inc": "r",
    "ing": "r",
    "int": "l",
    "rpy": "r",
    "rpu": "l",
    "lpy": "l",
    "lpu": "l",
}
ENZYME_LEFT_FOLDS = {
    amino_acid:
        {"s": 0, "l": 1, "r": -1}[fold]
    for amino_acid, fold in ENZYME_FOLDS.items()
}
