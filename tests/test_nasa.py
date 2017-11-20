import os
from os.path import join

import numpy as np

from chemkin.nasa import NASACoeffs


def get_example_data_file(file):
    return join("chemkin/example_data", file)


def test_create_db():
    try:
        nasa = NASACoeffs("temp.sqlite")
        nasa.create_db(get_example_data_file("example_thermo.xml"))
        assert np.allclose(nasa.get_coeffs("H", 'high')[0], np.array(
            [2.500000, -2.308430e-11, 1.615619e-14, - 4.735152e-18, 4.981974e-22, 25473.659900, - 0.446683]))
    except:
        raise
    finally:
        os.remove("temp.sqlite")
