# Documentation development

This page explains how to build, preview, and extend the ZEPHYRUS documentation locally.

## Overview

ZEPHYRUS' documentation is written in [Markdown](https://www.markdownguide.org/basic-syntax/) and generated using [Zensical](https://zensical.org/docs/get-started/).

Zensical is fully compatible with the project's `mkdocs.yml` configuration, which is used for the documentation build process and site structure.

The documentation source files live in the `docs/` directory of the repository.

## Local development

To build and preview the documentation locally, install the documentation dependencies and start the development server:

```console
pip install -e '.[docs]'
zensical serve
```

This command will generate the documentation and start a local server for previewing changes as you edit.

Once the server is running, open the URL shown in your terminal in a browser. In most cases, it will look similar to:

```console
http://127.0.0.1:8000
```
!!! note "Using a different port"
    If port 8000 is already busy, you can build the documentation on a different port:
    ```console
    zensical serve -a localhost:8001
    ```

!!! note "SSH connection"
    If you are connected to a remote server via SSH, for example on the Habrok, Snellius, or Kapteyn cluster, you may need to create an additional SSH tunnel. On your **local machine**, run:

    ```console
    ssh -L 4000:localhost:8000 account@server
    ```

    Then open `http://localhost:4000` in your browser.

## Documentation structure

Documentation source files are located in the [`docs`](https://github.com/FormingWorlds/ZEPHYRUS/tree/main/docs) directory, and mostly follow the [Diátaxis](https://diataxis.fr/) approach of documentation. This means they are categorised into:

- **Tutorials**: learning-oriented guides for new users, such as a first workflow tutorial.
- **How-to guides**: task-oriented instructions, such as installation, configuration, or test development.
- **Explanations**: conceptual material, such as model overviews or code architecture.
- **Reference**: lookup material, such as reference data, API-like documentation, or bibliographies.

Using this structure helps keep pages focused on a single purpose. When adding a new page, try to decide first which of these four categories it belongs to.

## Layout of the `docs/` directory

Within the `docs/` directory, documentation pages are generally organised according to the categories above. In addition, several extra directories contain assets and project-level documentation material.

The current `docs/` folder structure is shown below:

```text
.
├── assets                  # Images, diagrams and logos
├── Community               # Contact, contributing, code of conduct
├── Explanations            # Explanation pages such as a model overview
├── getting_started.md      # Getting started page with quicklinks
├── How-to                  # How-to guides such as this one
├── index.md                # Landing page (home)
├── javascripts             # Javascripts to configure features such as the clickable PROTEUS logo
├── overrides               # HTML overrides on top of the Zensical configurable style
├── Reference               # Reference pages, such as an API reference
│    └──  api
├── stylesheets             # CSS overrides on top of the Zensical configurable style
├── proteus_framework.md    # Overview of PROTEUS and links to other modules
└── Tutorials               # ZEPHYRUS tutorials 
```

## Adding a new page

When adding a new documentation page:

1. Create the Markdown file in the appropriate location under `docs/`.
2. Add the page to the navigation in `mkdocs.yml` under the `nav` section.
3. Start the local preview with `zensical serve`. If the preview is already open, saving `mkdocs.yml` should refresh the page (check this in the terminal).
4. Check that the page renders correctly and appears in the expected place in the navigation.
5. Verify that internal links, images, and code blocks display correctly.

For example, a page stored at `docs/Tutorials/model_superearth.md` might be added to `mkdocs.yml` like this:

```yaml
nav:
  - Tutorials:
    - Simulate a super-Earth: Tutorials/model_superearth.md

```

!!! info "Documentation files are relative to `docs/` directory"
    Zensical looks for documentation files relative to the `docs/` directory. In `mkdocs.yml`, add files like this: `Tutorials/model_superearth.md` and not like this: `docs/Tutorials/model_superearth.md`.