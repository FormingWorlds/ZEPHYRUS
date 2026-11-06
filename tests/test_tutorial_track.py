"""Integration test for the tutorial's stellar-history step.

The dispatcher tutorial closes on a track dispatched along a real MORS
stellar history, which reads the Spada grid under ``FWL_DATA``. The unit and
smoke tier does not download that data, so this section's verbatim
comparison runs here, where the nightly workflow has it.

The invariant under test:

- Documentation fidelity: the page states that every printed number is the
  output of the script section it includes, and the track's numbers come
  from the real stellar lookup rather than a mock, so only a run with the
  reference data present can hold that claim for them.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

from __future__ import annotations

import pytest

# The tests directory is not a package, so pytest puts it on the path and
# a sibling module imports by its bare name.
from test_examples import (
    DATA_DEPENDENT,
    run_sections,
    same_output,
    script_sections,
    tutorial_blocks,
)

pytestmark = [pytest.mark.integration, pytest.mark.timeout(300)]


def test_tutorial_stellar_track_prints_what_the_page_quotes(monkeypatch, tmp_path):
    """The track section reproduces its quoted output on the real tracks.

    Runs the whole script in order so the track inherits the namespace the
    earlier sections build, then compares the track's output against the
    page. Mocking the lookup would compare against numbers the page does
    not quote.
    """
    monkeypatch.chdir(tmp_path)
    printed = run_sections(list(script_sections()), {})
    quoted = dict(tutorial_blocks())
    compared = 0
    for name in DATA_DEPENDENT:
        if quoted[name] is None:
            continue
        assert same_output(printed[name], quoted[name]), (
            f'section {name} prints something other than the page quotes:\n'
            f'--- page ---\n{quoted[name]}\n--- script ---\n{printed[name]}'
        )
        compared += 1
    assert compared == 1, f'expected one data-dependent output block, compared {compared}'
    assert (tmp_path / 'output' / 'dispatcher_track.png').is_file()
