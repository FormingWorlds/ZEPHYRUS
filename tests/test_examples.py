"""Tests for the shipped examples under examples/.

Cross-cutting by design: the dispatcher example is the script behind the
tutorial page, which includes each of its sections verbatim and quotes what
each prints. These tests run the script section by section, in page order,
and check that the page includes every section once and quotes its output.
They also pin what the page says about the results: the regime sequence of
the flux sweep, the per-species closure that the coupled path relies on,
the limit behavior at zero XUV flux, the boundary search, the boil-off and
overflow labels, the error contract on a malformed state, and the label
change along a stellar history with the stellar lookup mocked.

Tier: smoke, because these run the real dispatch path end to end rather
than a mocked one. See docs/How-to/run_tests.md.
"""

from __future__ import annotations

import contextlib
import io
import math
import re
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from zephyrus.constants import rate_floor
from zephyrus.dispatcher import dispatch

pytestmark = [pytest.mark.smoke, pytest.mark.timeout(60)]

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'examples' / 'demo_dispatcher' / 'demo_dispatcher.py'
TUTORIAL = ROOT / 'docs' / 'Tutorials' / 'dispatch.md'
SCRIPT_REF = 'examples/demo_dispatcher/demo_dispatcher.py'

# Sections that read the MORS tracks under FWL_DATA, which this tier does not
# download; test_tutorial_track.py compares their output on the real tracks.
DATA_DEPENDENT = ('step7_track', 'step7_figure')


def script_sections() -> dict[str, str]:
    """The named sections of the example script, in file order."""
    pattern = re.compile(r'^# --8<-- \[start:(\w+)\]\n(.*?)^# --8<-- \[end:\1\]\n', re.S | re.M)
    return dict(pattern.findall(SCRIPT.read_text()))


def tutorial_blocks() -> list[tuple[str, str | None]]:
    """Each python fence on the tutorial page as (section, quoted output or None).

    Every fence must be an include of one script section; anything else on
    the page would be code that no test runs.
    """
    text = TUTORIAL.read_text()
    pattern = re.compile(r'```python\n(.*?)```(.*?)(?=```python|\Z)', re.S)
    include = re.compile(r'--8<-- "' + re.escape(SCRIPT_REF) + r':(\w+)"\n')
    blocks = []
    for code, after in pattern.findall(text):
        match = include.fullmatch(code)
        if match is None:
            raise AssertionError(f'tutorial fence is not a script include:\n{code}')
        quoted = re.search(r'```text\n(.*?)```', after, re.S)
        blocks.append((match.group(1), quoted.group(1).rstrip('\n') if quoted else None))
    return blocks


def run_sections(names, namespace: dict) -> dict[str, str]:
    """Execute the named script sections in order and return what each printed."""
    sections = script_sections()
    printed = {}
    for name in names:
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            exec(compile(sections[name], f'<{SCRIPT_REF}:{name}>', 'exec'), namespace)
        printed[name] = buffer.getvalue().rstrip('\n')
    return printed


_FLOAT = re.compile(r'-?\d+\.\d+(?:[eE][-+]?\d+)?|-?\d+[eE][-+]?\d+')


def same_output(printed, quoted, rel=1e-12):
    """Whether printed output matches the quoted text up to float rounding.

    Text between numbers must match exactly; each float may differ by
    ``rel``, far above the last-place differences between platforms'
    floating-point libraries and far below any change a coefficient makes.
    """
    if _FLOAT.sub('#', printed) != _FLOAT.sub('#', quoted):
        return False
    pairs = zip(_FLOAT.findall(printed), _FLOAT.findall(quoted), strict=True)
    return all(math.isclose(float(a), float(b), rel_tol=rel, abs_tol=0.0) for a, b in pairs)


@pytest.fixture(scope='module')
def tutorial(tmp_path_factory):
    """The script run once, every section but the stellar track, in a scratch directory.

    Returns the shared namespace and what each section printed. The script
    writes its figures to ``output/`` under the working directory, so the
    run happens in a temporary one.
    """
    workdir = tmp_path_factory.mktemp('tutorial')
    names = [name for name in script_sections() if name not in DATA_DEPENDENT]
    namespace: dict = {}
    with pytest.MonkeyPatch.context() as mp:
        mp.chdir(workdir)
        printed = run_sections(names, namespace)
    return dict(namespace=namespace, printed=printed, workdir=workdir)


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


def test_tutorial_includes_every_script_section_once_in_order():
    """The page includes each script section once, in the order of the script.

    A section the page omits is code a reader of the page never sees run; a
    section included twice or out of order would run differently from the
    script. The setup section comes first because every later one needs it.
    """
    included = [name for name, _ in tutorial_blocks()]
    assert included == list(script_sections())
    assert included[0] == 'setup'
    assert len(set(included)) == len(included)


def test_tutorial_sections_print_what_the_page_quotes(tutorial):
    """Every section the page quotes prints that output, run in page order.

    Text must match exactly and floats to a relative 1e-12, so a
    coefficient change anywhere in the package that moves a quoted number
    fails here while a last-place difference between platforms does not.
    The stellar track is compared in the integration tier.
    """
    compared = 0
    for name, quoted in tutorial_blocks():
        if quoted is None or name in DATA_DEPENDENT:
            continue
        printed = tutorial['printed'][name]
        assert same_output(printed, quoted), (
            f'section {name} prints something other than the page quotes:\n'
            f'--- page ---\n{quoted}\n--- script ---\n{printed}'
        )
        compared += 1
    # Guard against the extraction silently finding nothing.
    assert compared >= 12, f'only {compared} output blocks were compared'
    # Both figures the page shows are written by the sweep step.
    assert (tutorial['workdir'] / 'output' / 'dispatcher_regime_sweep.png').is_file()


def test_output_comparison_allows_only_float_rounding():
    """The tutorial comparison accepts a last-place difference and nothing more.

    A one-unit change in the sixteenth significant digit (the macOS and Linux
    difference that motivated it) passes; a relative change of 1e-9, a changed
    label, and a dropped number each fail.
    """
    page = "hydrodynamic:EL 660971.2983440236 {'C': 1.5e-3}"
    assert same_output("hydrodynamic:EL 660971.2983440235 {'C': 1.5e-3}", page)
    assert not same_output("hydrodynamic:EL 660971.2989 {'C': 1.5e-3}", page)
    assert not same_output("hydrodynamic:RR 660971.2983440236 {'C': 1.5e-3}", page)
    assert not same_output("hydrodynamic:EL 660971.2983440236 {'C': }", page)


@pytest.mark.physics_invariant
def test_dispatcher_example_sweep_crosses_two_boundaries(tutorial):
    """The documented flux sweep crosses both boundaries, in order.

    The tutorial states that the carbon dioxide planet is hydrostatic at
    the bottom of the swept range, energy limited in the middle, and
    recombination limited at the top. Pin that ordering, and the
    monotonicity of the rate inside the wind.
    """
    sweep = tutorial['namespace']['sweeps']['CO2']
    labels = [out.regime for out in sweep]
    ordered = [label for i, label in enumerate(labels) if i == 0 or label != labels[i - 1]]
    assert ordered == ['hydrostatic', 'hydrodynamic:EL', 'hydrodynamic:RR']
    # Inside the wind the rate rises with the flux; a sign or exponent slip
    # in the energy-limited chain would break the ordering.
    wind = [out.mdot for out in sweep if out.regime.startswith('hydrodynamic')]
    assert len(wind) > 10
    assert all(b > a for a, b in zip(wind[:-1], wind[1:]))
    # The hydrostatic corner is below the reporting floor, as the page says.
    assert sweep[0].mdot < rate_floor
    assert not sweep[0].diagnostics['rate_floor']['above_floor']


@pytest.mark.physics_invariant
def test_dispatcher_example_verdict_closes_and_respects_limits(tutorial):
    """The per-species split closes, and zero XUV flux drives no wind.

    Mass closure over the per-species rates is the contract a coupled run
    depends on when it debits element inventories. The zero-flux limit is
    the input edge case: with no XUV heating there is no wind to sustain,
    so the state falls to the hydrostatic branch and the rate drops far
    below the reporting floor.
    """
    ns = tutorial['namespace']
    result = ns['result']
    total = sum(result.per_species.values())
    assert total == pytest.approx(result.mdot, rel=1e-12, abs=0.0)
    # Scale guard: the wind rate on this planet is of order 1e6 kg/s, not
    # 1e0 (a dropped geometric factor) or 1e12 (a cgs slip).
    assert 1.0e5 < result.mdot < 1.0e8
    assert result.regime == 'hydrodynamic:EL'
    assert result.flags == {}

    quiet = dispatch(ns['build_state']('CO2', 1.0, 1.0, 0.0))
    assert quiet.regime == 'hydrostatic'
    assert 0.0 <= quiet.mdot < rate_floor
    # The closure holds on the hydrostatic branch too, where the split comes
    # from the branch itself rather than from the fractionation solver.
    assert sum(quiet.per_species.values()) == pytest.approx(quiet.mdot, rel=1e-12, abs=0.0)


@pytest.mark.physics_invariant
def test_dispatcher_example_boundary_search_brackets_a_label_change(tutorial):
    """The boundary flux has different labels on either side, and the band orders.

    The page quotes boundary fluxes from this search, so the property that
    matters is that each returned flux is a boundary at all: a search that
    returned an endpoint, or converged on the wrong side, would put the same
    label on both sides.
    """
    ns = tutorial['namespace']
    flux = ns['boundaries']['CO2'][0]
    below = dispatch(ns['build_state']('CO2', 1.0, 1.0, flux * 0.9)).regime
    above = dispatch(ns['build_state']('CO2', 1.0, 1.0, flux * 1.1)).regime
    assert (below, above) == ('hydrostatic', 'hydrodynamic:EL')
    assert 0.1 < flux < 10.0
    # Band edges come back in threshold order, 3 then 0.1: the stricter 0.1
    # needs a denser sonic point and so more flux. A swapped or shared
    # threshold would collapse this to equality.
    at_kn3, at_kn0p1 = ns['bands']['CO2']
    assert at_kn0p1 / at_kn3 > 1.5
    assert at_kn3 < flux < at_kn0p1


def test_dispatcher_example_reaches_the_two_extreme_labels(tutorial):
    """The boil-off and overflow cases return those labels, and the default does not.

    The overflow case admits the post-gate residual. With the default
    settings the same state is a bound XUV wind, which the page shows next,
    so the label must differ between the two calls on one state.
    """
    ns = tutorial['namespace']
    assert ns['inflated'].regime == 'boiloff'
    assert ns['puffy'].regime == 'roche_overflow'
    assert ns['default'].regime.startswith('hydrodynamic')
    for out in (ns['inflated'], ns['puffy'], ns['default']):
        assert math.isfinite(out.mdot) and out.mdot > 0.0
    # The rename keeps the winning rate: the admitted residual's own value.
    assert ns['puffy'].mdot == pytest.approx(
        ns['default'].diagnostics['bolometric']['rate_kg_s'], rel=1e-12, abs=0.0
    )


def test_dispatcher_example_rejects_a_malformed_state(tutorial):
    """A malformed state raises, and the valid state either side still dispatches.

    The framework reserves exceptions for input that is not physically
    posed. An eccentricity of one is outside the documented domain, so the
    call must raise rather than dispatch something.
    """
    state = tutorial['namespace']['build_state']('CO2', 1.0, 1.0, 10.0)
    state.e = 1.0
    with pytest.raises(ValueError, match='e must be in'):
        dispatch(state)
    state.e = 0.0
    state.F_xuv = -1.0
    with pytest.raises(ValueError, match='F_xuv'):
        dispatch(state)
    state.F_xuv = 10.0
    assert dispatch(state).regime == 'hydrodynamic:EL'


@pytest.mark.physics_invariant
def test_dispatcher_example_track_changes_regime_as_the_star_quiets(tutorial, tmp_path):
    """Along a declining XUV history the label changes once, downward.

    With the stellar lookup mocked, the same frozen atmosphere is
    dispatched at each age. The physical content is that the label belongs
    to the state: a decaying flux drives the state out of the wind regime,
    and the rate falls with it. A dilution or unit slip in the flux
    conversion would leave the label fixed for the whole track.
    """
    namespace = dict(tutorial['namespace'])
    with patch('mors.Star', return_value=_fake_star()), pytest.MonkeyPatch.context() as mp:
        mp.chdir(tmp_path)
        (tmp_path / 'output').mkdir()
        run_sections(DATA_DEPENDENT, namespace)
    track = namespace['track']
    labels = [row['result'].regime for row in track]
    assert labels[0] == 'hydrodynamic:EL'
    assert labels[-1] == 'hydrostatic'
    assert sum(1 for a, b in zip(labels[:-1], labels[1:]) if a != b) == 1
    # The flux declines and so does the rate, by orders of magnitude across
    # the crossing rather than by a small step.
    assert track[0]['F_xuv'] > 10.0 * track[-1]['F_xuv']
    assert track[0]['result'].mdot > 1.0e6 * track[-1]['result'].mdot
    # After the crossing the hydrostatic rate is flux independent, because
    # that branch carries no XUV physics: the last two samples agree.
    tail = [row['result'].mdot for row in track if row['result'].regime == 'hydrostatic']
    assert tail[-1] == pytest.approx(tail[0], rel=1e-12, abs=0.0)
    assert (tmp_path / 'output' / 'dispatcher_track.png').is_file()
