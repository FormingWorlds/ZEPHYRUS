"""Integration test for the tutorial's stellar-history section.

The dispatcher tutorial closes on a track dispatched along a real MORS
stellar history, which reads the Spada grid under ``FWL_DATA``. The unit and
smoke tier does not download that data, so the snippet is executed here
instead, where the nightly workflow has it: the same character-for-character
comparison the smoke tier applies to the other eighteen output blocks, on the
one block it cannot run.

The invariant under test:

- Documentation fidelity: the page states that every printed number is the
  verbatim output of a snippet the reader can run, and the track's numbers
  come from the real stellar lookup rather than a mock, so only a run with the
  reference data present can hold that claim for them.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

from __future__ import annotations

import contextlib
import io
from pathlib import Path

import pytest

# The tests directory is not a package, so pytest puts it on the path and
# a sibling module imports by its bare name.
from test_examples import TUTORIAL_DATA_DEPENDENT, _tutorial_blocks

pytestmark = [pytest.mark.integration, pytest.mark.timeout(300)]


def test_tutorial_stellar_track_prints_what_the_page_quotes(monkeypatch):
    """The track snippet reproduces its quoted output on the real tracks.

    Runs the whole page in order so the track inherits the namespace the
    earlier snippets build, then compares the track block's output against
    the page. The smoke tier covers the other blocks; this one exists because
    the stellar lookup needs reference data that tier does not fetch, and
    mocking the lookup would compare against numbers the page does not quote.
    """
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]))
    namespace: dict = {}
    compared = 0
    for index, code, quoted in _tutorial_blocks(skip_data_dependent=False):
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            exec(compile(code, f'<tutorial block {index}>', 'exec'), namespace)
        if quoted is None or TUTORIAL_DATA_DEPENDENT not in code:
            continue
        printed = buffer.getvalue().rstrip('\n')
        assert printed == quoted, (
            f'tutorial block {index} prints something other than the page quotes:\n'
            f'--- page ---\n{quoted}\n--- code ---\n{printed}'
        )
        compared += 1
    assert compared == 1, f'expected one data-dependent block, compared {compared}'
