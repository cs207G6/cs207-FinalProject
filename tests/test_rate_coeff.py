from chemkin.rate_coeff import RateCoeff, ModifiedArrhenius, Arrhenius, Constant


def test_notimplemented():
    try:
        RateCoeff().get_K(0)
        assert False
    except NotImplementedError:
        assert True


def test_negativeModifiedArrhenius():
    try:
        ModifiedArrhenius(-10, 20, 30).get_K(50)
        assert False
    except ValueError:
        assert True
    try:
        ModifiedArrhenius(10, 20, 30, R=-1).get_K(50)
        assert False
    except ValueError:
        assert True
    try:
        ModifiedArrhenius(10, 20, 30).get_K(-50)
        assert False
    except ValueError:
        assert True


def test_negativeArrhenius():
    try:
        Arrhenius(-10, 20, 30).get_K(50)
        assert False
    except ValueError:
        assert True
    try:
        Arrhenius(10, 20, R=-1).get_K(50)
        assert False
    except ValueError:
        assert True
    try:
        Arrhenius(10, 20, 30).get_K(-50)
        assert False
    except ValueError:
        assert True


def test_negativeConstant():
    try:
        Constant(-10).get_K(50)
        assert False
    except ValueError:
        assert True
