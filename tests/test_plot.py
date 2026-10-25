"""Tests for ``src/zephyrus/plot.py``.

``plot.py`` is a utility source (figures for dispatcher results), so it is
exempt from the physics-invariant requirement. What these tests guard is
what a reader of a figure relies on: every regime gets its own marker and
color, a rate below the axis floor is drawn open rather than dropped, a
label change is marked where it happens, the optional theme package is
optional, and malformed input raises instead of drawing something
misleading.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

from __future__ import annotations

import sys
import types

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

from zephyrus import plot  # noqa: E402
from zephyrus.dispatcher import REGIME_LABELS  # noqa: E402

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def _sweep_panel(**overrides):
    """A three-point sweep crossing from hydrostatic to EL, one rate below the floor."""
    panel = dict(
        title='CO2',
        F_xuv=[0.1, 1.0, 10.0],
        mdot=[1.0e-120, 5.0e5, 2.0e6],
        regime=['hydrostatic', 'hydrodynamic:EL', 'hydrodynamic:EL'],
        boundaries=[0.75],
        band=(0.59, 2.6),
    )
    panel.update(overrides)
    return panel


def test_every_regime_label_has_a_distinct_marker_and_color():
    """Each label of REGIME_LABELS has its own marker and its own color.

    A figure that encoded a regime by color alone, or reused a marker,
    would make two regimes indistinguishable for a reader with a color
    vision deficiency or in grayscale print.
    """
    assert set(plot.REGIME_MARKERS) == set(REGIME_LABELS)
    assert set(plot.REGIME_COLORS) == set(REGIME_LABELS)
    assert len(set(plot.REGIME_MARKERS.values())) == len(REGIME_LABELS)
    assert len(set(plot.REGIME_COLORS.values())) == len(REGIME_LABELS)


def test_apply_theme_is_optional(monkeypatch):
    """Without proteus-mpl the theme is skipped; with it, use() is called once.

    The figures must draw without the optional package, and must apply its
    theme when present rather than reading keys an older release lacks.
    """
    monkeypatch.setitem(sys.modules, 'proteus_mpl', None)
    assert plot.apply_theme() is False

    calls = []
    fake = types.SimpleNamespace(use=lambda: calls.append('use'))
    monkeypatch.setitem(sys.modules, 'proteus_mpl', fake)
    assert plot.apply_theme() is True
    assert calls == ['use']


def test_mono_family_requires_the_font_to_be_registered(monkeypatch):
    """The mono family is named only when matplotlib has the font.

    A package that imports but ships no font file (proteus-mpl 1.0) made
    matplotlib warn ``findfont: Font family 'Spline Sans Mono' not found``;
    the check is on the registered fonts, not on the package.
    """
    monkeypatch.setattr(plot.font_manager.fontManager, 'ttflist', [])
    assert plot.mono_family() == 'monospace'

    entry = types.SimpleNamespace(name=plot.MONO_FAMILY)
    monkeypatch.setattr(plot.font_manager.fontManager, 'ttflist', [entry])
    assert plot.mono_family() == 'Spline Sans Mono'


def test_regime_sweep_draws_rates_below_the_floor_open():
    """A rate below the axis floor is drawn at the floor with an open marker.

    The hydrostatic point at 1e-120 kg/s must stay visible and be marked as
    reading zero, while the wind points above the floor are filled. Each
    regime appears once in the figure legend.
    """
    fig = plot.plot_regime_sweep(
        [_sweep_panel(), _sweep_panel(title='N2-O2')], rate_floor=1.0e3
    )
    ax = fig.axes[0]
    markers = [line for line in ax.lines if line.get_linestyle() == 'None']
    hydrostatic = [line for line in markers if line.get_label() == 'hydrostatic']
    assert len(hydrostatic) == 1
    assert hydrostatic[0].get_markerfacecolor() == 'none'
    assert hydrostatic[0].get_ydata()[0] == pytest.approx(1.0e3, rel=1e-12, abs=0.0)
    wind = [line for line in markers if line.get_label() == 'hydrodynamic:EL']
    assert wind[0].get_markerfacecolor() != 'none'
    legend_labels = [text.get_text() for text in fig.legends[0].get_texts()]
    assert legend_labels == ['hydrodynamic:EL', 'hydrostatic']
    assert len(fig.axes) == 2
    plt.close(fig)


def test_regime_sweep_marks_boundaries_and_band():
    """Each boundary is a dashed line and the band is one shaded span.

    The band is the width of a boundary across its criterion, so a missing
    span would print a boundary as a line, which the docs forbid.
    """
    fig = plot.plot_regime_sweep([_sweep_panel(boundaries=[0.75, 500.0])])
    ax = fig.axes[0]
    dashed = [line for line in ax.lines if line.get_linestyle() == '--']
    assert sorted(line.get_xdata()[0] for line in dashed) == [0.75, 500.0]
    assert len(ax.patches) == 1
    # Edge case: a panel without boundaries or band draws neither.
    plt.close(fig)
    bare = plot.plot_regime_sweep([_sweep_panel(boundaries=(), band=None)])
    assert not [line for line in bare.axes[0].lines if line.get_linestyle() == '--']
    assert len(bare.axes[0].patches) == 0
    plt.close(bare)


def test_regime_sweep_rejects_malformed_panels():
    """Empty input, unequal lengths, unknown labels, and a bad floor raise.

    Drawing a mislabeled or truncated sweep would be worse than failing, so
    each case raises ValueError before a figure exists.
    """
    before = plt.get_fignums()
    with pytest.raises(ValueError, match='empty'):
        plot.plot_regime_sweep([])
    with pytest.raises(ValueError, match='mdot has 2 entries'):
        plot.plot_regime_sweep([_sweep_panel(mdot=[1.0, 2.0])])
    with pytest.raises(ValueError, match='unknown regime'):
        plot.plot_regime_sweep([_sweep_panel(regime=['hydrostatic', 'EL', 'EL'])])
    with pytest.raises(ValueError, match='rate_floor'):
        plot.plot_regime_sweep([_sweep_panel()], rate_floor=0.0)
    plt.close('all')
    assert plt.get_fignums() == [] or plt.get_fignums() == before


def _track(**overrides):
    """A five-sample track: a wind that ends, then a flat exosphere rate."""
    track = dict(
        age_myr=[1.0, 10.0, 100.0, 1000.0, 5000.0],
        mdot=[1.0e7, 1.0e6, 2.0e5, 4.7e-2, 4.7e-2],
        regimes=['hydrodynamic:EL'] * 3 + ['hydrostatic'] * 2,
        per_species=[{'H': 1.0e6, 'C': 2.0e6}] * 3 + [{'H': 4.7e-2, 'C': 1.0e-140}] * 2,
        inconsistent=[False, True, True, False, False],
    )
    track.update(overrides)
    return track


def test_track_marks_the_label_change_and_floors_tiny_rates():
    """One label change gives one dashed line per panel, at the geometric mean age.

    The heavy-element rate of 1e-140 kg/s is drawn at the species floor
    rather than vanishing from a log axis, and the screen's rejected ages
    are shaded in both panels.
    """
    fig = plot.plot_track(**_track(), species_floor=1.0e-25)
    ax_rate, ax_species = fig.axes[:2]
    dashed = [line for line in ax_rate.lines if line.get_linestyle() == '--']
    assert len(dashed) == 1
    assert dashed[0].get_xdata()[0] == pytest.approx((100.0 * 1000.0) ** 0.5, rel=1e-12)
    carbon = [line for line in ax_species.lines if line.get_label() == 'C'][0]
    assert min(carbon.get_ydata()) == pytest.approx(1.0e-25, rel=1e-12, abs=0.0)
    assert len(ax_rate.patches) == 1 and len(ax_species.patches) == 1
    plt.close(fig)


def test_track_without_label_change_or_screen_draws_neither():
    """A track with one label and no rejected ages has no dashed line and no span.

    Edge case of the change detection: a spurious crossing line would claim
    a regime change the sequence does not contain.
    """
    fig = plot.plot_track(**_track(regimes=['hydrostatic'] * 5, inconsistent=None))
    ax_rate = fig.axes[0]
    assert not [line for line in ax_rate.lines if line.get_linestyle() == '--']
    assert len(ax_rate.patches) == 0
    plt.close(fig)


def test_track_rejects_malformed_input():
    """Empty, unequal, non-increasing, or mislabeled tracks raise ValueError."""
    with pytest.raises(ValueError, match='empty'):
        plot.plot_track([], [], [], [])
    with pytest.raises(ValueError, match='per_species has 1 entries'):
        plot.plot_track(**_track(per_species=[{'H': 1.0}]))
    with pytest.raises(ValueError, match='increase'):
        plot.plot_track(**_track(age_myr=[1.0, 10.0, 10.0, 1000.0, 5000.0]))
    with pytest.raises(ValueError, match='unknown regime'):
        plot.plot_track(**_track(regimes=['wind'] * 5))
    with pytest.raises(ValueError, match='species_floor'):
        plot.plot_track(**_track(), species_floor=-1.0)
