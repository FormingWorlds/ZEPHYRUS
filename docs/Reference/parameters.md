# ZEPHYRUS parameter reference

This is a reference page for all parameters and constants used in ZEPHYRUS. For the physical model, see the [model overview](../Explanations/model.md).

---

## Physical constants (`constants.py`)

| Name | Symbol | Value | Units |
|---|---|---|---|
| `kb` | $k_B$ | $1.38 \times 10^{-23}$ | J K⁻¹ |
| `G` | $G$ | $6.6743 \times 10^{-11}$ | m³ kg⁻¹ s⁻² |
| `G_cgs` | $G_\mathrm{cgs}$ | $6.6743 \times 10^{-8}$ | cm³ g⁻¹ s⁻² |
| `c` | $c$ | $2.99792458 \times 10^{8}$ | m s⁻¹ |

## Unit conversions (`constants.py`)

| Name | Factor | Conversion |
|---|---|---|
| `s2yr` | $1/(3600 \cdot 24 \cdot 365)$ | seconds → years |
| `erg2joule` | $10^{-7}$ | erg → J |
| `au2m` | $1.496 \times 10^{11}$ | au → m |
| `au2cm` | $1.496 \times 10^{13}$ | au → cm |
| `ergpersecondtowatt` | $10^{-7}$ | erg s⁻¹ → W |
| `ergcm2stoWm2` | $10^{-3}$ | erg s⁻¹ cm⁻² → W m⁻² |

---

## Sun–Earth reference values (`planets_parameters.py`)

### Sun

| Name | Symbol | Value | Units |
|---|---|---|---|
| `Rs` | $R_\odot$ | $6.957 \times 10^{8}$ | m |
| `Ms` | $M_\odot$ | $1.98847 \times 10^{30}$ | kg |
| `Ls` | $L_\odot$ | $3.828 \times 10^{26}$ | W |
| `age_sun` | – | $4.603 \times 10^{9}$ | yr |


### Earth

| Name | Symbol | Value | Units |
|---|---|---|---|
| `Re` | $R_\oplus$ | $6.378 \times 10^{6}$ | m |
| `Me` | $M_\oplus$ | $5.9722 \times 10^{24}$ | kg |
| `Me_atm` | – | $5.15 \times 10^{18}$ | kg |
| `Fxuv_earth_10Myr` | $F_\mathrm{XUV,\oplus}(10\,\mathrm{Myr})$ | $14.67$ | W m⁻² |
| `Fxuv_earth_today` | $F_\mathrm{XUV,\oplus}$ | $4.64 \times 10^{-3}$ | W m⁻² |
| `age_earth` | – | $4.543 \times 10^{9}$ | yr |
| `e_earth` | – | $0.017$ | – |
| `a_earth` | – | $1$ | au |

`Fxuv_earth_10Myr` is taken from Fig. 9 of Wordsworth et al. (2018).

---

## TOI-561 reference values (`planets_parameters.py`)

### TOI-561 (star, Weiss et al. 2021)

| Name | Value | Errorbar | Units |
|---|---|---|---|
| `R_TOI561` | $0.832\,R_\odot$ | $0.019\,R_\odot$ | m |
| `M_TOI561` | $0.805\,M_\odot$ | $0.030\,M_\odot$ | kg |
| `L_TOI561` | $0.522\,L_\odot$ | $0.017\,L_\odot$ | W |
| `age_TOI561` | $10 \times 10^{9}$ | $3 \times 10^{9}$ | yr |


### TOI-561 b (planet, Brinkman et al. 2023)

| Name | Value | Errorbar | Units |
|---|---|---|---|
| `R_TOI561b` | $1.37\,R_\oplus$ | $0.04\,R_\oplus$ | m |
| `M_TOI561b` | $2.24\,M_\oplus$ | $0.20\,M_\oplus$ | kg |
| `a_TOI561b` | $0.0106$ | $0.0004$ | au |
