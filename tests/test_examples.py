"""Tests for the shipped examples under examples/.

Cross-cutting by design: an example is a script plus the documentation page
that quotes its output, so these tests exercise the example's own entry
points and pin the results the pages state. What they guard is the regime
sequence of the flux sweep, the per-species closure that the coupled path
relies on, the limit behavior at zero XUV flux, the error contract on a
malformed state, and the label change along a stellar history with the
stellar lookup mocked.

Tier: smoke, because these run the real dispatch path end to end rather
than a mocked one. See docs/How-to/run_tests.md.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import math
import re
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from zephyrus.dispatcher import dispatch

pytestmark = [pytest.mark.smoke, pytest.mark.timeout(60)]

_EXAMPLE = (
    Path(__file__).resolve().parents[1] / 'examples' / 'demo_dispatcher' / 'demo_dispatcher.py'
)


def _load_example():
    """Import the example module from its path, without a package."""
    spec = importlib.util.spec_from_file_location('demo_dispatcher', _EXAMPLE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _fake_star(n: int = 12):
    """A stand-in for a MORS star with a plausible declining XUV history.

    Saturated for the first few samples then decaying as a power law, with
    a bolometric luminosity of order solar: the shape matters here, because
    a constant history could not produce a regime change and a zero one
    would hide a dilution error.
    """
    age = np.logspace(0, 4, n)  # Myr
    l_x = 1.0e30 * np.minimum(1.0, (age / 30.0) ** -1.5)  # erg/s
    l_euv = 2.0 * l_x
    l_bol = np.full(n, 3.828e33)  # erg/s, one solar luminosity

    class Star:
        Tracks = {'Age': age, 'Lx': l_x, 'Leuv': l_euv, 'Lbol': l_bol}

    return Star()


@pytest.mark.physics_invariant
def test_dispatcher_example_sweep_crosses_two_boundaries():
    """The documented flux sweep crosses both boundaries, in order.

    The tutorial states that the carbon dioxide planet is hydrostatic at
    the bottom of the swept range, energy limited in the middle, and
    recombination limited at the top. Pin that ordering, the closure of the
    per-species split, and the monotonicity of the rate inside the wind.
    """
    example = _load_example()
    rows = example.flux_sweep('CO2')
    labels = [row['regime'] for row in rows]

    # The two edges of the swept range are the two extreme regimes.
    assert labels[0] == 'hydrostatic'
    assert labels[-1] == 'hydrodynamic:RR'
    # Every label seen is one the framework declares, and the sweep visits
    # all three in order: hydrostatic, then EL, then RR, with no return.
    ordered = [label for i, label in enumerate(labels) if i == 0 or label != labels[i - 1]]
    assert ordered == ['hydrostatic', 'hydrodynamic:EL', 'hydrodynamic:RR']
    # Inside the wind the rate rises with the flux; a sign or exponent slip
    # in the energy-limited chain would break the ordering.
    wind = [row for row in rows if row['regime'].startswith('hydrodynamic')]
    assert len(wind) > 10
    assert all(b['mdot'] > a['mdot'] for a, b in zip(wind[:-1], wind[1:]))
    # The hydrostatic corner is below the reporting floor, which is the
    # tutorial's point about a computed rate that is not a physical one.
    assert rows[0]['mdot'] < example.RATE_FLOOR
    assert not rows[0]['above_floor']


@pytest.mark.physics_invariant
def test_dispatcher_example_verdict_closes_and_respects_limits():
    """The per-species split closes, and zero XUV flux drives no wind.

    Mass closure over the per-species rates is the contract a coupled run
    depends on when it debits element inventories. The zero-flux limit is
    the input edge case: with no XUV heating there is no wind to sustain,
    so the state falls to the hydrostatic branch and the rate drops far
    below the reporting floor.
    """
    example = _load_example()
    result = dispatch(example.build_state('CO2', 1.0, 1.0, 10.0))
    total = sum(result.per_species.values())
    assert total == pytest.approx(result.mdot, rel=1e-12, abs=0.0)
    # Scale guard: the wind rate on this planet is of order 1e6 kg/s, not
    # 1e0 (a dropped geometric factor) or 1e12 (a cgs slip).
    assert 1.0e5 < result.mdot < 1.0e8
    assert result.regime == 'hydrodynamic:EL'
    # A clean call reports no flags at all.
    assert result.flags == {}

    quiet = dispatch(example.build_state('CO2', 1.0, 1.0, 0.0))
    assert quiet.regime == 'hydrostatic'
    assert quiet.mdot < example.RATE_FLOOR
    assert quiet.mdot >= 0.0
    # The closure holds on the hydrostatic branch too, where the split comes
    # from the branch itself rather than from the fractionation solver.
    assert sum(quiet.per_species.values()) == pytest.approx(quiet.mdot, rel=1e-12, abs=0.0)


@pytest.mark.physics_invariant
def test_dispatcher_example_boundary_bisection_brackets_a_label_change():
    """The bisected boundary flux has different labels on either side.

    The documentation quotes boundary fluxes from this bisection, so the
    property that matters is that the returned flux is a boundary at all:
    a bisection that returned an endpoint, or converged on the wrong side,
    would put the same label on both sides.
    """
    example = _load_example()
    flux = example.boundary_flux('CO2', 0.1, 10.0)
    below = dispatch(example.build_state('CO2', 1.0, 1.0, flux * 0.9)).regime
    above = dispatch(example.build_state('CO2', 1.0, 1.0, flux * 1.1)).regime
    assert below != above
    assert (below, above) == ('hydrostatic', 'hydrodynamic:EL')
    # The boundary lies strictly inside the bracket it was given, so a
    # bisection that fell back to an endpoint fails here.
    assert 0.1 < flux < 10.0
    # The band edges come back in threshold order, 3 first then 0.1. The
    # stricter threshold (0.1) demands a denser sonic point, so it needs
    # more flux to call a wind: the band is ordered, not just wide. A
    # swapped or shared threshold would collapse this to equality.
    at_kn3, at_kn0p1 = example.boundary_band('CO2')
    assert at_kn0p1 > at_kn3
    assert at_kn0p1 / at_kn3 > 1.5
    assert at_kn3 < flux < at_kn0p1


def test_dispatcher_example_rejects_a_malformed_state():
    """A malformed state raises, and returns nothing.

    The framework reserves exceptions for input that is not physically
    posed. An eccentricity of one is outside the documented domain, so the
    call must raise rather than dispatch something.
    """
    example = _load_example()
    state = example.build_state('CO2', 1.0, 1.0, 10.0)
    state.e = 1.0
    with pytest.raises(ValueError, match='e must be in'):
        dispatch(state)
    # A negative flux is rejected on the same contract, and the valid
    # state either side of these edits still dispatches.
    state.e = 0.0
    state.F_xuv = -1.0
    with pytest.raises(ValueError, match='F_xuv'):
        dispatch(state)
    state.F_xuv = 10.0
    assert dispatch(state).regime == 'hydrodynamic:EL'


@pytest.mark.physics_invariant
def test_dispatcher_example_track_changes_regime_as_the_star_quiets():
    """Along a declining XUV history the label changes once, downward.

    With the stellar lookup mocked, the same frozen atmosphere is
    dispatched at each age. The physical content is that the label belongs
    to the state: a decaying flux drives the state out of the wind regime,
    and the rate falls with it. A dilution or unit slip in the flux
    conversion would leave the label fixed for the whole track.
    """
    example = _load_example()
    with patch('mors.Star', return_value=_fake_star()):
        rows = example.stellar_track(n_samples=12)

    labels = [row['regime'] for row in rows]
    assert labels[0] == 'hydrodynamic:EL'
    assert labels[-1] == 'hydrostatic'
    changes = sum(1 for a, b in zip(labels[:-1], labels[1:]) if a != b)
    assert changes == 1
    # The flux declines and so does the rate, by orders of magnitude across
    # the crossing rather than by a small step.
    assert rows[0]['F_xuv'] > 10.0 * rows[-1]['F_xuv']
    assert rows[0]['mdot'] > 1.0e6 * rows[-1]['mdot']
    # After the crossing the hydrostatic rate is flux independent, because
    # that branch carries no XUV physics: the last two samples agree.
    tail = [row['mdot'] for row in rows if row['regime'] == 'hydrostatic']
    assert tail[-1] == pytest.approx(tail[0], rel=1e-12, abs=0.0)
    assert all(math.isfinite(row['mdot']) for row in rows)


TUTORIAL = Path(__file__).resolve().parents[1] / 'docs' / 'Tutorials' / 'dispatch.md'


def _tutorial_blocks():
    """Every python fence in the dispatcher tutorial with the output it quotes."""
    text = TUTORIAL.read_text()
    pattern = re.compile(r'```python\n(.*?)```(.*?)(?=```python|\Z)', re.S)
    for index, (code, after) in enumerate(pattern.findall(text), start=1):
        quoted = re.search(r'```text\n(.*?)```', after, re.S)
        yield index, code, (quoted.group(1).rstrip('\n') if quoted else None)


def test_tutorial_snippets_print_what_the_page_quotes():
    """Every tutorial snippet runs in order and prints its quoted output.

    The page states that every printed number is the verbatim output of a
    snippet the reader can run, which is a claim about the documentation that
    only a test can hold. The snippets share one namespace and run in the
    order they appear, as a reader would execute them, and each quoted output
    block must match what the preceding snippet printed, character for
    character. This is the guard against the drift that has to be repaired by
    hand otherwise: a coefficient change three modules away moves a number
    here, and nothing else notices.
    """
    namespace: dict = {}
    compared = 0
    for index, code, quoted in _tutorial_blocks():
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            exec(compile(code, f'<tutorial block {index}>', 'exec'), namespace)
        printed = buffer.getvalue().rstrip('\n')
        if quoted is None:
            continue
        assert printed == quoted, (
            f'tutorial block {index} prints something other than the page quotes:\n'
            f'--- page ---\n{quoted}\n--- code ---\n{printed}'
        )
        compared += 1
    # Guard against the extraction silently finding nothing, which would make
    # the assertions above vacuous.
    assert compared >= 15, f'only {compared} tutorial output blocks were compared'
