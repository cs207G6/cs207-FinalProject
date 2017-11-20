from chemkin.parser import DataParser
from chemkin.nasa import NASACoeffs
from os.path import join
import numpy as np

nasa = NASACoeffs()


def get_example_data_file(file):
    return join("chemkin/example_data", file)


def parse_file(file_name):
    return DataParser().parse_file(get_example_data_file(file_name), nasa)


def test_xmlErrors():
    try:
        parse_file('test_parse.xml')  # no phase labels
        assert False
    except Exception:
        assert True


def test_negA():
    try:
        parse_file("test_negA.xml")
        assert False
    except Exception:
        assert True


def test_btype():
    try:
        parse_file("test_btype.xml")
        assert False
    except Exception:
        assert True


def test_Atype():
    try:
        parse_file("test_Atype.xml")
        assert False
    except Exception:
        assert True


def test_negk():
    try:
        parse_file("test_negk.xml")
        assert False
    except Exception:
        assert True


def test_errorChem():
    try:
        parse_file("test_errorChem.xml")
        assert False
    except Exception:
        assert True


def test_oneMoreReaction():
    try:
        parse_file("test_oneMoreReaction.xml")
    except Exception:
        assert False


def test_length():
    assert (len(parse_file('rxns.xml')) == 3)


def test_idCollision():
    try:
        parse_file("test_idCollision.xml")
        assert False
    except Exception:
        assert True


def test_progRateWrongDimension():
    cd = parse_file('rxns.xml')
    try:
        cd.get_progress_rate([1, 2, 3, 4, 5], 1000)
        assert False
    except Exception:
        assert True


def test_progRateWrongSp():
    try:
        parse_file('test_wrongSp.xml')
        assert False
    except Exception:
        assert True


def test_progRateWrongSp2():
    try:
        parse_file('test_wrongSp2.xml')
        assert False
    except Exception:
        assert True


def test_progRateNonIrr():
    cd = parse_file('test_nonirr.xml')
    try:
        cd.get_progress_rate([1, 2, 3, 4, 5, 6], 1000)
        assert False
    except Exception:
        assert True


def test_progRateNonEle():
    cd = parse_file('test_nonele.xml')
    try:
        cd.get_progress_rate([1, 2, 3, 4, 5, 6], 1000)
        assert False
    except Exception:
        assert True


def test_mixreverse():
    cd = parse_file('test_mixreverse.xml')
    try:
        cd.get_progress_rate([1, 2, 3, 4, 5, 6, 7, 8], 1000)
        assert (True)
    except Exception as err:
        assert (False)


def test_allreverse():
    try:
        cd = parse_file('test_allreverse.xml')
        cd.get_progress_rate([1, 2, 3, 4, 5, 6, 7, 8], 1000)
        assert (True)
    except Exception as err:
        assert (False)


def test_RxWtWrongSp():
    try:
        cd = parse_file('test_RxWtWrongSp.xml')
        cd.get_progress_rate([1, 2, 3, 4, 5, 6, 7, 8], 1000)
        assert (False)
    except Exception as err:
        assert (True)


def test_WrongSp_in_array():
    try:
        parse_file('test_WrongSp_in_array.xml')
        assert False
    except Exception:
        assert True


def test_wrong_num_products():
    try:
        parse_file('test_wrong_num_products.xml')
        assert False
    except Exception:
        assert True


def test_wrong_num_reactants():
    try:
        parse_file('test_wrong_num_reactants.xml')
        assert False
    except Exception:
        assert True


def test_unknown_coeff():
    try:
        parse_file('test_notimplementedCoeff.xml')
        assert False
    except NotImplementedError:
        assert True

def test_db():
    nasa = NASACoeffs()
    assert(nasa.get_coeffs("H", 'high')[0] == np.array([3.282538 ,1.483088e-03 ,-7.579667e-07, 2.094706e-10, -2.167178e-14, -1088.457720, 5.453231]))
