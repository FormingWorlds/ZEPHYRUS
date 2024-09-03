# Contributing guidelines

## Development

### Building the documentation

The documentation is written in [markdown](https://www.markdownguide.org/basic-syntax/), and uses [mkdocs](https://www.mkdocs.org/) to generate the pages.

To build the documentation for yourself:

```console
pip install -e .[docs]
mkdocs serve
```

You can find the documentation source in the [docs](https://github.com/FormingWorlds/ZEPHYRUS/tree/main/docs) directory.
If you are adding new pages, make sure to update the listing in the [`mkdocs.yml`](https://github.com/FormingWorlds/ZEPHYRUS/blob/main/mkdocs.yml) under the `nav` entry.

The documentation is hosted on [readthedocs](https://readthedocs.io/projects/fwl-zephyrus).

### Running tests

ZEPHYRUS uses [pytest](https://docs.pytest.org/en/latest/) to run the tests. You can run the tests for yourself using:

```console
pytest
```

To check coverage:

```console
coverage run -m pytest
coverage report  # to output to terminal
coverage html    # to generate html report
```


### Making a release

The versioning scheme we use is [CalVer](https://calver.org/).

0. Update requirements files:

```console
python tools/generate_requirements_txt.py
pip-compile -o requirements_full.txt pyproject.toml
```

1. Bump the version (`release`/`patch`) as needed

```console
bump-my-version bump release
# 24.06.26
```

2. Commit and push your changes.

3. Make a new [release](https://github.com/FormingWorlds/ZEPHYRUS/releases). Make sure to set the tag to the specified version, e.g. `24.06.26`.

4. The [upload to pypi](https://pypi.org/project/fwl-zephyrus) is triggered when a release is published and handled by [this workflow](https://github.com/FormingWorlds/ZEPHYRUS/actions/workflows/publish.yaml).
