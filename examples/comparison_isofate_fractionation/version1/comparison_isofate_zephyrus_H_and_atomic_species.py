#!/usr/bin/env python
from zephyrus_binary import isocalc_zephyrus
from isofate_binary import isocalc
from isofunks import *
from constants import *
import numpy as np
import matplotlib.pyplot as plt

# Choose species mixture: options 'H/O', 'H/N', 'H/S'
species_mixture = 'H/O'

# Map mixture to species labels and molecular weights
species_map = {
    'H/O': dict(light='H', heavy='O', mu_light=mu_H, mu_heavy=mu_O, protosolar_ratio=OtoH_protosolar, color_heavy='red'),
    'H/N': dict(light='H', heavy='N', mu_light=mu_H, mu_heavy=mu_N, protosolar_ratio=NtoH_protosolar, color_heavy='deeppink'),
    'H/S': dict(light='H', heavy='S', mu_light=mu_H, mu_heavy=mu_S, protosolar_ratio=StoH_protosolar, color_heavy='gold'),
    'H/C': dict(light='H', heavy='C', mu_light=mu_H, mu_heavy=mu_C, protosolar_ratio=CtoH_protosolar, color_heavy='coral'),
}

# Initialization parameters
f_atm = 0.0001
Mp = 1.0 * Me
Mstar = 1.0 * Ms
Rstar = 1.0 * Rs
Tstar = 5600  # [K]
T = 1000      # [K]
d = 1.0 * au2m
L = Luminosity(Rstar, Tstar)
Fp = Insolation(L, d)
F0 = Fp * 10**(-3.5) * (Mstar / Ms)

# Function to run, plot, and save for a given mixture
def process_mixture(species_mixture):
    spec = species_map[species_mixture]
    light = spec['light']; heavy = spec['heavy']
    mu_l = spec['mu_light']; mu_h = spec['mu_heavy']
    protosolar_ratio = spec['protosolar_ratio']
    color_h = spec['color_heavy']

    # Run models
    iso = isocalc(f_atm, Mp, Mstar, F0, Fp, T, d,
                  time=5e9, mechanism='XUV', species=species_mixture, rad_evol=True,
                  mu=mu_solar, eps=0.15, activity='medium', flux_model='power law', stellar_type='G',
                  Rp_override=False, t_sat=50e6, f_atm_final='null', n_TO_final='null',
                  n_steps=int(1e5), t0=1e6, rho_rcb=1.0, Johnson=False, RR=True,
                  f_pred=False, thermal=True, beta=-1.23)
    zep = isocalc_zephyrus(f_atm, Mp, Mstar, F0, Fp, T, d,
                            time=5e9, mechanism='XUV', species=species_mixture, rad_evol=True,
                            mu=mu_solar, eps=0.15, activity='medium', flux_model='power law', stellar_type='G',
                            Rp_override=False, t_sat=50e6, f_atm_final='null', n_TO_final='null',
                            n_steps=int(1e5), t0=1e6, rho_rcb=1.0, Johnson=False, RR=True,
                            f_pred=False, thermal=True, beta=-1.23)

    t_iso = iso['time']; t_zep = zep['time']

    # Plotting setup
    fig, axs = plt.subplots(4, 2, figsize=(10, 10))
    axs = axs.flatten()

    # Panel 1: Flux comparison
    axs[0].plot(t_iso * s2yr, iso['phi'][0], '-', label='IsoFATE', color='black')
    axs[0].plot(t_zep * s2yr, zep['phi'][0], ':', label='Zephyrus', color='black', linewidth=2)
    axs[0].plot(t_iso * s2yr, iso['phi'][0], '-', label='Total flux', color='darkviolet')
    axs[0].plot(t_iso * s2yr, iso['phic'][0], '--', label='Critical flux', color='orange')
    axs[0].plot(t_iso * s2yr, iso['Phi1'][0] * mu_l, '-', label=f'{light} flux', color='dodgerblue')
    axs[0].plot(t_iso * s2yr, iso['Phi2'][0] * mu_h, '-', label=f'{heavy} flux', color=color_h)
    axs[0].plot(t_zep * s2yr, zep['phi'][0], ':', color='darkviolet', linewidth=2)
    axs[0].plot(t_zep * s2yr, zep['phic'][0], ':', color='orange', linewidth=2)
    axs[0].plot(t_zep * s2yr, zep['Phi1'][0] * mu_l, ':', color='dodgerblue', linewidth=2)
    axs[0].plot(t_zep * s2yr, zep['Phi2'][0] * mu_h, ':', color=color_h, linewidth=2)
    axs[0].set_xscale('log'); axs[0].set_yscale('log')
    axs[0].set_xlabel('Time [yr]', fontsize=14); axs[0].set_ylabel(r'$\Phi$ [kg m$^{-2}$ s$^{-1}$]', fontsize=14); axs[0].legend()


    # Panel 2: Moles
    axs[1].plot(t_iso*s2yr, iso['N1'][0]/avogadro, '-', label=f'IsoFATE', color='black')
    axs[1].plot(t_zep*s2yr, zep['N1'][0]/avogadro, ':', label='Zephyrus', color='black', linewidth=2)
    axs[1].plot(t_iso*s2yr, iso['N1'][0]/avogadro, '-', label=f'$N_{light}$', color='dodgerblue')
    axs[1].plot(t_iso*s2yr, iso['N2'][0]/avogadro, '-', label=f'$N_{heavy}$', color=color_h)
    axs[1].plot(t_zep*s2yr, zep['N1'][0]/avogadro, ':', color='dodgerblue', linewidth=2)
    axs[1].plot(t_zep*s2yr, zep['N2'][0]/avogadro, ':', color=color_h, linewidth=2)
    axs[1].set_xscale('log'); axs[1].set_yscale('log')
    axs[1].set_xlabel('Time [yr]', fontsize=14); axs[1].set_ylabel(r'N$_i$ [mol]', fontsize=14); axs[1].legend()

    # Panel 8: Radius
    axs[7].plot(t_iso*s2yr, iso['rp'][0]/Re, '-', label='IsoFATE', color='green')
    axs[7].plot(t_zep*s2yr, zep['rp'][0]/Re, ':', label='Zephyrus', color='green', linewidth=2)
    axs[7].set_xscale('log'); axs[7].set_xlabel('Time [yr]', fontsize=14); axs[7].set_ylabel(r'Radius [R$_\oplus$]', fontsize=14); axs[7].legend()

    # Panel 4: Species fraction
    axs[3].plot(t_iso*s2yr, iso['x1'][0]*100, '-', label='IsoFATE', color='black')
    axs[3].plot(t_zep*s2yr, zep['x1'][0]*100, ':', label='Zephyrus', color='black', linewidth=2)
    axs[3].plot(t_iso*s2yr, iso['x1'][0]*100, '-', label=f'$x_{light}$', color='dodgerblue')
    axs[3].plot(t_iso*s2yr, iso['x2'][0]*100, '-', label=f'$x_{heavy}$', color=color_h)
    axs[3].plot(t_zep*s2yr, zep['x1'][0]*100, ':', color='dodgerblue', linewidth=2)
    axs[3].plot(t_zep*s2yr, zep['x2'][0]*100, ':', color=color_h, linewidth=2)
    axs[3].set_xscale('log'); axs[3].set_xlabel('Time [yr]', fontsize=14); axs[3].set_ylabel(r'$\chi_{i}^{ATM}$ [%]', fontsize=14); axs[3].legend()

    # Panel 5: Mole ratio relative
    axs[4].plot(t_iso*s2yr, (iso['N2'][0]/iso['N1'][0])/protosolar_ratio, '-', label='IsoFATE', color='maroon')
    axs[4].plot(t_zep*s2yr, (zep['N2'][0]/zep['N1'][0])/protosolar_ratio, ':', label='Zephyrus', color='maroon', linewidth=2)
    axs[4].set_xscale('log'); axs[4].set_yscale('log'); axs[4].set_xlabel('Time [yr]', fontsize=14); axs[4].set_ylabel(f'{heavy}/{light} [Solar]', fontsize=14); axs[4].legend()

    # Panel 6: Atmospheric fraction
    axs[5].plot(t_iso*s2yr, iso['fenv'][0]*100, '-', label='IsoFATE', color='gold')
    axs[5].plot(t_zep*s2yr, zep['fenv'][0]*100, ':', label='Zephyrus', color='gold', linewidth=2)
    axs[5].set_xscale('log'); axs[5].set_xlabel('Time [yr]', fontsize=14); axs[5].set_ylabel(r'f$_{atm}$ [%]', fontsize=14); axs[5].legend()

    # Panel 7: Atmosphere mass
    axs[6].plot(t_iso*s2yr, iso['menv'][0]/Me, '-', label='IsoFATE', color='darkblue')
    axs[6].plot(t_zep*s2yr, zep['menv'][0]/Me, ':', label='Zephyrus', color='darkblue', linewidth=2)
    axs[6].set_xscale('log'); axs[6].set_yscale('log'); axs[6].set_xlabel('Time [yr]', fontsize=14); axs[6].set_ylabel(r'M$_{atm}$ [M$_{\oplus}$]', fontsize=14); axs[6].legend()

    # Panel 3: Mass loss
    axs[2].plot(t_iso*s2yr, iso['mloss'][0]/Me, '-', label='IsoFATE', color='black')
    axs[2].plot(t_zep*s2yr, zep['mloss'][0]/Me, ':', label='Zephyrus', color='black', linewidth=2)
    axs[2].plot(t_iso*s2yr, iso['mloss'][0]/Me, '-', label='Total loss', color='green')
    axs[2].plot(t_zep*s2yr, zep['mloss'][0]/Me, ':', color='green', linewidth=2)
    axs[2].plot(t_iso*s2yr, iso['mlossl'][0]/Me, '-', label=f'{light} loss', color='dodgerblue')
    axs[2].plot(t_iso*s2yr, iso['mlossh'][0]/Me, '-', label=f'{heavy} loss', color=color_h)
    axs[2].plot(t_zep*s2yr, zep['mlossl'][0]/Me, ':', color='dodgerblue', linewidth=2)
    axs[2].plot(t_zep*s2yr, zep['mlossh'][0]/Me, ':', color=color_h, linewidth=2)
    axs[2].set_xscale('log'); axs[2].set_yscale('log'); axs[2].set_xlabel('Time [yr]', fontsize=14); axs[2].set_ylabel(r'$\Delta$ mass [M$_\oplus$]', fontsize=14); axs[2].legend()

    # Box input parameters
    param_text = (
        f"$f_{{atm}}$ = {f_atm*100} % \n"
        f"$M_p$ = {Mp/Me:.2f} $M_\\oplus$\n"
        f"$T_{{p}}$ = {T} K\n"
        f"$d$ = {d/au2m:.2f} AU\n"
        f"$M_{{star}}$ = {Mstar/Ms:.2f} $M_\\odot$\n"
        f"$R_{{star}}$ = {Rstar/Rs:.2f} $R_\\odot$\n"
        f"$T_{{star}}$ = {Tstar} K")
    fig.text(0.595, 0.302, param_text, fontsize=10, va='bottom', ha='left',
            bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="black", alpha=0.8))

    # Title and save
    fig.suptitle(f'Comparison fractionation IsoFATE vs Zephyrus : {species_mixture}', fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    o_species = species_mixture.replace('/', '')
    output_dir = 'plot_comparison_isofate_zephyrus_H_and_atomic_species'
    filename = f'{output_dir}/8panels_comparison_{o_species}.pdf'
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"Saved {species_mixture} plot as {filename}")


# Loop over all mixtures
for mixture in species_map.keys():
    process_mixture(mixture)
