import pytest
import numpy as np
from numpy.testing import assert_allclose
import mors

from zephyrus.constants import ergcm2stoWm2, au2cm, au2m
from zephyrus.planets_parameters import Me, Re, e_earth, a_earth
from zephyrus.escape import EL_escape

TEST_DATA = (
    (  150.0, (7.53422377e+28, 1.12676493e+29, 0.06685410, 20506.0813)),
    ( 2480.0, (1.78184954e+28, 4.48318430e+28, 0.02227668,  6832.8987)),
    ( 5813.0, (8.98042499e+27, 3.07577786e+28, 0.01412977,  4334.0088)),
    ( 7927.0, (7.02400992e+27, 2.82190252e+28, 0.01253142,  3843.7476)),
    ( 9101.0, (6.53002620e+27, 2.88714085e+28, 0.01258774,  3861.0233)),
    (10020.0, (6.67781787e+27, 3.20900636e+28, 0.01378475,  4228.1816)),
    (11050.0, (4.33094381e+27, 3.69943595e+28, 0.01469410,  4507.1043)),
)

@pytest.mark.parametrize("inp,expected", TEST_DATA)
def test_earth(inp, expected):

    age_star = inp

    # Load the stellar evolution tracks from MORS
    mors.DownloadEvolutionTracks('Spada')
    star = mors.Star(Mstar=1.0, Omega=1.0)

    #Lx = np.interp(age_star, star.Tracks['Age'], star.Tracks['Lx'])
    #Leuv = np.interp(age_star, star.Tracks['Age'], star.Tracks['Leuv'])
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
