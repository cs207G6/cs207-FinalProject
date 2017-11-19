from chemkin.parser import DataParser
from chemkin.nasa import NASACoeffs
from os.path import join

nasa = NASACoeffs()


def get_example_data_file(file):
    return join("chemkin/example_data", file)


def parser_file(file_name):
    return DataParser().parse_file(get_example_data_file(file_name), nasa)


def test_xmlErrors():
    try:
        parser_file('test_parse.xml')  # no phase labels
        assert False
    except Exception:
        assert True


def test_negA():
    try:
        parser_file("test_negA.xml")
        assert False
    except Exception:
        assert True


def test_btype():
    try:
        parser_file("test_btype.xml")
        assert False
    except Exception:
        assert True


def test_Atype():
    try:
        parser_file("test_Atype.xml")
        assert False
    except Exception:
        assert True


def test_negk():
    try:
        parser_file("test_negk.xml")
        assert False
    except Exception:
        assert True


def test_errorChem():
    try:
        parser_file("test_errorChem.xml")
        assert False
    except Exception:
        assert True


def test_oneMoreReaction():
    try:
        parser_file("test_oneMoreReaction.xml")
    except Exception:
        assert False


def test_length():
    assert (len(parser_file('rxns.xml')) == 3)


def test_idCollision():
    try:
        parser_file("test_idCollision.xml")
        assert False
    except Exception:
        assert True


def test_progRateWrongDimension():
    cd = parser_file('rxns.xml')
    try:
        cd.get_progress_rate([1, 2, 3, 4, 5], 1000)
        assert False
    except Exception:
        assert True


def test_progRateWrongSp():
    try:
        parser_file('test_wrongSp.xml')
        assert False
    except Exception:
        assert True


def test_progRateWrongSp2():
    try:
        parser_file('test_wrongSp2.xml')
        assert False
    except Exception:
        assert True


def test_progRateNonIrr():
    cd = parser_file('test_nonirr.xml')
    try:
        cd.get_progress_rate([1, 2, 3, 4, 5, 6], 1000)
        assert False
    except Exception:
        assert True


def test_progRateNonEle():
    cd = parser_file('test_nonele.xml')
    try:
        cd.get_progress_rate([1, 2, 3, 4, 5, 6], 1000)
        assert False
    except Exception:
        assert True


def test_mixreverse():
    cd = parser_file('test_mixreverse.xml')
    try:
        cd.get_progress_rate([1, 2, 3, 4, 5, 6, 7, 8], 1000)
        assert (True)
    except Exception as err:
        assert (False)


def test_allreverse():
    try:
        cd = parser_file('test_allreverse.xml')
        cd.get_progress_rate([1, 2, 3, 4, 5, 6, 7, 8], 1000)
        assert (True)
    except Exception as err:
        assert (False)


def test_RxWtWrongSp():
    try:
        cd = parser_file('test_RxWtWrongSp.xml')
        cd.get_progress_rate([1, 2, 3, 4, 5, 6, 7, 8], 1000)
        assert (False)
    except Exception as err:
        assert (True)


def test_WrongSp_in_array():
    try:
        parser_file('test_WrongSp_in_array.xml')
        assert False
    except Exception:
        assert True
