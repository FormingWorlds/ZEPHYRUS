import matplotlib.pyplot as plt

from zephyrus.constants import *
from zephyrus.planets_parameters import *
from zephyrus.XUV_flux import Fxuv_mors

########################### Fxuv computation ####################################

star_age_sun, star_fxuv_sun = Fxuv_mors(1.0, 1.0, 1.0)
star_age_2, star_fxuv_2 = Fxuv_mors(1.0, 0.5, 1.0)

########################### Plot ####################################

plt.figure(figsize=(8, 6))

plt.loglog(star_age_sun,star_fxuv_sun, color='gold', linestyle='-', label='1.0 AU, 1.0 Msun, 1.0 OmegaSun')
plt.loglog(star_age_2,star_fxuv_2, color='red', linestyle='-', label='1.0 AU, 0.5 Msun, 1.0 OmegaSun')

plt.xlabel('Age [Myr]', fontsize=15)
plt.ylabel(r'XUV Flux received by the planet [W $m^{-2}$]', fontsize=15)
plt.legend()
plt.grid(alpha=0.5)
plt.title('Test of the Fxuv_mors() function, using Star()', fontsize=15)


plt.savefig('output/test_Fxuv_mors_function.pdf',dpi=180)




