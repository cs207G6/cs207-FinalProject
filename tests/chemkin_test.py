
from ..chemkin.parser import DataParser
from os.path import join


def get_example_data_file(file):
    return join("data", file)


def test_xmlErrors():
    try:
        DataParser().parse_file(get_example_data_file('test_parse.xml'))  # no phase labels
        assert False
    except Exception:
        assert True


def test_negA():
    try:
        DataParser().parse_file(get_example_data_file("test_negA.xml"))
        assert False
    except Exception:
        assert True


def test_btype():
    try:
        DataParser().parse_file(get_example_data_file("test_btype.xml"))
        assert False
    except Exception:
        assert True


def test_Atype():
    try:
        DataParser().parse_file(get_example_data_file("test_Atype.xml"))
        assert False
    except Exception:
        assert True


def test_negk():
    try:
        DataParser().parse_file(get_example_data_file("test_negk.xml"))
        assert False
    except Exception:
        assert True


def test_errorChem():
    try:
        DataParser().parse_file(get_example_data_file("test_errorChem.xml"))
        assert False
    except Exception:
        assert True


def test_oneMoreReaction():
    try:
        DataParser().parse_file(get_example_data_file("test_oneMoreReaction.xml"))
    except Exception:
        assert False


def test_length():
    assert (len(DataParser().parse_file(get_example_data_file('rxns.xml'))) == 3)


def test_idCollision():
    try:
        DataParser().parse_file(get_example_data_file("test_idCollision.xml"))
        assert False
    except Exception:
        assert True


def test_progRateWrongDimension():
    cd = DataParser().parse_file(get_example_data_file('rxns.xml'))
    try:
        cd.get_progress_rate([1, 2, 3, 4, 5], 100)
        assert False
    except Exception:
        assert True


def test_progRateWrongSp():
    try:
        DataParser().parse_file(get_example_data_file('test_wrongSp.xml'))
        assert False
    except Exception:
        assert True


def test_progRateWrongSp2():
    try:
        DataParser().parse_file(get_example_data_file('test_wrongSp2.xml'))
        assert False
    except Exception:
        assert True


def test_progRateNonIrr():
    cd = DataParser().parse_file(get_example_data_file('test_nonirr.xml'))
    try:
        cd.get_progress_rate([1, 2, 3, 4, 5, 6], 100)
        assert True
    except Exception:
        assert False


def test_progRateNonEle():
    cd = DataParser().parse_file(get_example_data_file('test_nonele.xml'))
    try:
        cd.get_progress_rate([1, 2, 3, 4, 5, 6], 100)
        assert False
    except Exception:
        assert True
        
def test_mixreverse():
	cd = DataParser().parse_file(get_example_data_file('test_mixreverse.xml'))
	try:
		cd.get_progress_rate([1,2,3,4,5,6,7,8],100)
		assert(True)
	except Exception as err:
		assert(False)

def test_allreverse():
	try:
		cd = DataParser().parse_file(get_example_data_file('test_allreverse.xml'))
		cd.get_progress_rate([1,2,3,4,5,6,7,8],100)
		assert(True)
	except Exception as err:
		assert(False)

def test_RxWtWrongSp():
	try:
		cd = DataParser().parse_file(get_example_data_file('test_RxWtWrongSp.xml'))
		cd.get_progress_rate([1,2,3,4,5,6,7,8],100)
		assert(False)
	except Exception as err:
		assert(True)


def test_eqtIncst():
	try:
		cd = DataParser().parse_file(get_example_data_file('test_eqtIncst.xml'))
		cd.get_progress_rate([1,2,3,4,5,6,7,8],100)
		assert(False)
	except Exception as err:
		assert(True)

def test_incstReverse():
	try:
		cd = DataParser().parse_file(get_example_data_file('test_incstReverse.xml'))
		cd.get_progress_rate([1,2,3,4,5,6,7,8],100)
		assert(False)
	except Exception as err:
		assert(True)

def test_WrongSp_in_array():
    try:
        DataParser().parse_file(get_example_data_file('test_WrongSp_in_array.xml'))
        assert False
    except Exception:
        assert True

