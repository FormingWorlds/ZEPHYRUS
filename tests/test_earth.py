import pytest
import numpy as np
from numpy.testing import assert_allclose
import mors

from zephyrus.constants import ergcm2stoWm2, au2cm, au2m
from zephyrus.planets_parameters import Me, Re, e_earth, a_earth
from zephyrus.escape import EL_escape

TEST_DATA = (
    (  150.0, (0.06685410, 20506.0813)),
    ( 2480.0, (0.02227668,  6832.8987)),
#    (10020.0, (0.01378475,  4228.1816)),
)

@pytest.mark.parametrize("inp,expected", TEST_DATA)
def test_earth(inp, expected):

    age_star = inp

    # Load the stellar evolution tracks from MORS
    mors.DownloadEvolutionTracks('Spada')
    star = mors.Star(Mstar=1.0, Omega=1.0)

    Lxuv_star = (star.Value(age_star, 'Lx')+star.Value(age_star, 'Leuv'))
    Fxuv_star_SI = Lxuv_star * ergcm2stoWm2 / (4*np.pi*(a_earth*au2cm)**2) 
    escape = EL_escape(False, a_earth*au2m, e_earth, Me, 1.0, 0.15, Re, Re, Fxuv_star_SI)   # Compute EL escape     [kg s-1]

    ret = (
        Fxuv_star_SI,
        escape,
        )

    assert_allclose(ret, expected, rtol=1e-5, atol=0)
