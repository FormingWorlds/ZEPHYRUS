import pytest
import numpy as np
from numpy.testing import assert_allclose
import mors

from zephyrus.constants import ergcm2stoWm2, au2cm, au2m
from zephyrus.planets_parameters import Me, Re, e_earth, a_earth
from zephyrus.escape import EL_escape

TEST_DATA = (
    (  150.0, (7.45671140e+28, 1.11900066e+29, 0.06630242, 20336.8630)),
    ( 5813.0, (8.99939840e+27, 3.08011304e+28, 0.01415193,  4340.8062)),
    (10020.0, (7.17720092e+27, 3.36660286e+28, 0.01452269,  4454.5274)),
)

@pytest.mark.parametrize("inp,expected", TEST_DATA)
def test_earth(inp, expected):

    age_star = inp

    # Load the stellar evolution tracks from MORS
    mors.DownloadEvolutionTracks('Spada')
    star = mors.Star(Mstar=1.0, Omega=1.0)

    Lx = star.Value(age_star, 'Lx')
    Leuv = star.Value(age_star, 'Leuv')
    Fxuv_star_SI = (Lx+Leuv) * ergcm2stoWm2 / (4*np.pi*(a_earth*au2cm)**2)
    escape = EL_escape(False, a_earth*au2m, e_earth, Me, 1.0, 0.15, Re, Re, Fxuv_star_SI)   # Compute EL escape     [kg s-1]

    ret = (
        Lx,
        Leuv,
        Fxuv_star_SI,
        escape,
        )

    assert_allclose(ret, expected, rtol=1e-5, atol=0)
