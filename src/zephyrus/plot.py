"""Figures for dispatcher results.

Two figures: a flux sweep with each point marked by its regime label, and
one atmosphere dispatched along a stellar history. Both take plain arrays,
so they plot any sequence of `dispatch` results, and both return the
matplotlib figure without saving it.

The colors follow the PROTEUS visual language. The hex values are fixed
here, so the figures look the same with or without the optional
`proteus-mpl` package; `apply_theme` adds its style sheet and fonts when it
is installed.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.figure import Figure

from zephyrus.dispatcher import REGIME_LABELS

ESCAPE_COLOR = '#1B6FA8'  # PROTEUS domain color of chemistry and escape
RULE_COLOR = '#3E4A55'

# Regime colors from the shared categorical cycle; boil-off takes the one red mark.
REGIME_COLORS = {
    'boiloff': '#E23D28',
    'hydrodynamic:EL': '#1B6FA8',
    'hydrodynamic:RR': '#4FA3D9',
    'hydrodynamic:PL': '#14406B',
    'hydrostatic': '#7A8894',
    'roche_overflow': '#593E74',
}

# One marker per label, so no figure encodes a regime in color alone.
REGIME_MARKERS = {
    'boiloff': 'D',
    'hydrodynamic:EL': 'o',
    'hydrodynamic:RR': 's',
    'hydrodynamic:PL': 'v',
    'hydrostatic': '^',
    'roche_overflow': 'X',
}

ELEMENT_COLORS = {
    'H': '#1B6FA8',
    'He': '#4FA3D9',
    'C': '#E23D28',
    'N': '#57A05C',
    'O': '#C8860F',
}

MONO_FAMILY = 'Spline Sans Mono'


def apply_theme() -> bool:
    """Apply the PROTEUS matplotlib theme when `proteus-mpl` is installed.

    Returns
    -------
    bool
        True when the theme was applied, False when the package is absent.
        The figures draw with the same colors either way.
    """
    try:
        import proteus_mpl
    except ImportError:
        return False
    proteus_mpl.use()
    return True


def mono_family() -> str:
    """Font family for in-axes annotations.

    Returns
    -------
    str
        'Spline Sans Mono' when matplotlib has that font registered (it
        arrives with `proteus-mpl` 1.1 and later through `apply_theme`),
        otherwise the generic 'monospace', so matplotlib never warns about
        a missing face.
    """
    names = {entry.name for entry in font_manager.fontManager.ttflist}
    return MONO_FAMILY if MONO_FAMILY in names else 'monospace'


def _check_regimes(regimes: Sequence[str]) -> None:
    unknown = sorted(set(regimes) - set(REGIME_LABELS))
    if unknown:
        raise ValueError(f'unknown regime labels {unknown}; expected labels from REGIME_LABELS')


def _check_lengths(n: int, **arrays) -> None:
    for name, values in arrays.items():
        if len(values) != n:
            raise ValueError(f'{name} has {len(values)} entries, expected {n}')


def plot_regime_sweep(panels: Sequence[Mapping], rate_floor: float = 1.0e3) -> Figure:
    """Dispatched rate against XUV flux, one panel per planet, marked by regime.

    Parameters
    ----------
    panels : sequence of mapping
        One mapping per panel, with keys ``title`` (str), ``F_xuv`` (XUV
        fluxes, W m-2), ``mdot`` (dispatched rates, kg s-1), and ``regime``
        (labels), all of one length, and optionally ``boundaries`` (fluxes
        of label changes, W m-2, drawn as dashed lines) and ``band`` (a pair
        of fluxes, W m-2, shaded as the width of a boundary).
    rate_floor : float
        Lowest rate on the axis [kg s-1]. A rate below it is drawn at the
        floor with an open marker, since on these axes it reads as zero.

    Returns
    -------
    matplotlib.figure.Figure
        The figure, unsaved.

    Raises
    ------
    ValueError
        If ``panels`` is empty, ``rate_floor`` is not positive, a panel's
        arrays differ in length, or a label is not a regime label.
    """
    if not panels:
        raise ValueError('panels is empty')
    if not rate_floor > 0.0:
        raise ValueError(f'rate_floor must be positive, got {rate_floor}')
    mono = mono_family()
    fig, axes = plt.subplots(1, len(panels), figsize=(5.5 * len(panels), 4.6), sharey=True)
    axes = np.atleast_1d(axes)
    top = rate_floor
    for ax, panel in zip(axes, panels):
        flux = np.asarray(panel['F_xuv'], dtype=float)
        mdot = np.asarray(panel['mdot'], dtype=float)
        regimes = list(panel['regime'])
        _check_lengths(flux.size, mdot=mdot, regime=regimes)
        _check_regimes(regimes)
        shown = np.maximum(mdot, rate_floor)
        top = max(top, float(shown.max()))
        ax.plot(flux, shown, color=ESCAPE_COLOR, linewidth=1.2, alpha=0.5, zorder=0)
        for label in REGIME_LABELS:
            in_label = np.array([r == label for r in regimes])
            for on_scale in (True, False):
                mask = in_label & ((mdot >= rate_floor) == on_scale)
                if not mask.any():
                    continue
                ax.plot(
                    flux[mask],
                    shown[mask],
                    linestyle='none',
                    marker=REGIME_MARKERS[label],
                    markersize=8,
                    markerfacecolor=REGIME_COLORS[label] if on_scale else 'none',
                    markeredgecolor=REGIME_COLORS[label],
                    label=label,
                )
        band = panel.get('band')
        if band:
            ax.axvspan(
                min(band), max(band), color=REGIME_COLORS['hydrostatic'], alpha=0.18, lw=0
            )
            ax.text(
                math.sqrt(band[0] * band[1]),
                0.86,
                'boundary width',
                transform=ax.get_xaxis_transform(),
                ha='center',
                fontsize=11,
                color=RULE_COLOR,
                fontfamily=mono,
            )
        for flux_b in panel.get('boundaries', ()):
            ax.axvline(flux_b, color=RULE_COLOR, linestyle='--', linewidth=1.0)
            ax.annotate(
                f'{flux_b:.3g} W m$^{{-2}}$',
                xy=(flux_b, 0.94),
                xycoords=ax.get_xaxis_transform(),
                xytext=(3, 0),
                textcoords='offset points',
                fontsize=11,
                color=RULE_COLOR,
                fontfamily=mono,
            )
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('XUV flux [W m$^{-2}$]')
        ax.set_title(panel['title'])
    axes[0].set_ylim(rate_floor / 3.0, top * 30.0)
    axes[0].set_ylabel('Mass loss rate [kg s$^{-1}$]')
    axes[0].text(
        0.03,
        0.14,
        'open markers: rate reads as zero',
        transform=axes[0].transAxes,
        fontsize=11,
        fontfamily=mono,
    )
    # One legend entry per label across all panels, filled marker preferred.
    entries = {}
    for ax in axes:
        for handle, label in zip(*ax.get_legend_handles_labels()):
            if label not in entries or handle.get_markerfacecolor() != 'none':
                entries[label] = handle
    ordered = [label for label in REGIME_LABELS if label in entries]
    fig.legend(
        [entries[label] for label in ordered],
        ordered,
        loc='upper center',
        ncol=min(len(ordered), 3),
        bbox_to_anchor=(0.5, 0.11),
    )
    fig.subplots_adjust(bottom=0.25, top=0.9, left=0.09, right=0.98, wspace=0.08)
    return fig


def plot_track(
    age_myr: Sequence[float],
    mdot: Sequence[float],
    regimes: Sequence[str],
    per_species: Sequence[Mapping[str, float]],
    inconsistent: Sequence[bool] | None = None,
    species_floor: float = 1.0e-25,
) -> Figure:
    """Bulk and per-element rates of one atmosphere along a stellar history.

    The upper panel is the bulk rate marked by regime, the lower one the
    per-element rates. Dashed lines mark each label change, at the geometric
    mean of the two ages around it, and the ages the self-consistency screen
    rejects are shaded.

    Parameters
    ----------
    age_myr : sequence of float
        Stellar ages [Myr], increasing.
    mdot : sequence of float
        Dispatched bulk rates [kg s-1].
    regimes : sequence of str
        Regime labels.
    per_species : sequence of mapping
        Per-element rates [kg s-1], one mapping per age.
    inconsistent : sequence of bool, optional
        `diagnostics['self_consistency']['inconsistent']` per age.
    species_floor : float
        Lowest element rate on the axis [kg s-1]; smaller rates are drawn
        at it.

    Returns
    -------
    matplotlib.figure.Figure
        The figure, unsaved.

    Raises
    ------
    ValueError
        If the sequences are empty or differ in length, the ages do not
        increase, ``species_floor`` is not positive, or a label is not a
        regime label.
    """
    ages = np.asarray(age_myr, dtype=float)
    rates = np.asarray(mdot, dtype=float)
    regimes = list(regimes)
    if ages.size == 0:
        raise ValueError('age_myr is empty')
    if not species_floor > 0.0:
        raise ValueError(f'species_floor must be positive, got {species_floor}')
    flags = list(inconsistent) if inconsistent is not None else [False] * ages.size
    _check_lengths(ages.size, mdot=rates, regimes=regimes, per_species=per_species)
    _check_lengths(ages.size, inconsistent=flags)
    if np.any(np.diff(ages) <= 0.0):
        raise ValueError('age_myr must increase')
    _check_regimes(regimes)

    mono = mono_family()
    fig, (ax_rate, ax_species) = plt.subplots(2, 1, figsize=(8.5, 7.0), sharex=True)
    bad = ages[np.asarray(flags, dtype=bool)]
    if bad.size:
        for ax in (ax_rate, ax_species):
            ax.axvspan(
                bad.min(), bad.max(), color=REGIME_COLORS['hydrostatic'], alpha=0.15, lw=0
            )
        ax_rate.text(
            math.sqrt(bad.min() * bad.max()),
            0.9,
            'rejected by the consistency screen',
            transform=ax_rate.get_xaxis_transform(),
            ha='center',
            fontsize=11,
            color=RULE_COLOR,
            fontfamily=mono,
        )

    positive = rates[rates > 0.0]
    rate_floor = positive.min() / 3.0 if positive.size else species_floor
    ax_rate.plot(ages, np.maximum(rates, rate_floor), color=ESCAPE_COLOR, lw=1.2, alpha=0.5)
    for label in REGIME_LABELS:
        mask = np.array([r == label for r in regimes])
        if mask.any():
            ax_rate.plot(
                ages[mask],
                np.maximum(rates[mask], rate_floor),
                linestyle='none',
                marker=REGIME_MARKERS[label],
                markersize=7,
                color=REGIME_COLORS[label],
                label=label,
            )
    for i in range(ages.size - 1):
        if regimes[i] != regimes[i + 1]:
            crossing = math.sqrt(ages[i] * ages[i + 1])
            for ax in (ax_rate, ax_species):
                ax.axvline(crossing, color=RULE_COLOR, linestyle='--', lw=1.0)
            ax_rate.annotate(
                f'{regimes[i]} ends near {crossing:.0f} Myr',
                xy=(crossing, 0.06),
                xycoords=ax_rate.get_xaxis_transform(),
                xytext=(-6, 0),
                textcoords='offset points',
                ha='right',
                fontsize=11,
                color=RULE_COLOR,
                fontfamily=mono,
            )

    elements = sorted({el for rates_i in per_species for el in rates_i})
    for element in elements:
        values = [max(rates_i.get(element, 0.0), species_floor) for rates_i in per_species]
        ax_species.plot(
            ages, values, color=ELEMENT_COLORS.get(element, RULE_COLOR), label=element
        )

    ax_rate.set_yscale('log')
    ax_rate.set_ylabel('Bulk rate [kg s$^{-1}$]')
    ax_rate.legend(loc='lower left', title='regime')
    ax_species.set_xscale('log')
    ax_species.set_yscale('log')
    ax_species.set_ylim(bottom=species_floor / 3.0)
    ax_species.set_xlabel('Stellar age [Myr]')
    ax_species.set_ylabel('Element rate [kg s$^{-1}$]')
    ax_species.legend(loc='lower left', title='element', ncol=3)
    fig.subplots_adjust(left=0.12, right=0.98, top=0.97, bottom=0.09, hspace=0.08)
    return fig
