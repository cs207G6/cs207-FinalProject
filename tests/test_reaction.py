from chemkin.parser import DataParser
from chemkin.nasa import NASACoeffs
from os.path import join
import numpy as np

nasa = NASACoeffs()


def get_example_data_file(file):
    return join("chemkin/example_data", file)


def parse_file(file_name):
    return DataParser().parse_file(get_example_data_file(file_name), nasa)


def test_unkownNASA():
    rd = parse_file("rxns_unknownNASA.xml")
    try:
        rd.get_progress_rate(np.ones(len(rd.species)), 1000)
        assert False
    except NotImplementedError:
        assert True


def test_nunkownT():
    rd = parse_file("rxns_reversible.xml")
    try:
        rd.get_progress_rate(np.ones(len(rd.species)), 100)
        assert False
    except NotImplementedError:
        assert True


def test_reversible():
    rd = parse_file("rxns_reversible.xml")
    rd.get_progress_rate(np.ones(len(rd.species)), 3000)


def test_reversible_mixed():
    rd = parse_file("rxns_reversible_mixed.xml")
    rd.get_progress_rate(np.ones(len(rd.species)), 3000)
