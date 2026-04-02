# Contributing to SolarSim

Thank you for helping improve SolarSim. This document describes how we expect changes to be prepared and reviewed.

## Getting started

1. **Fork** [fraware/solarsim](https://github.com/fraware/solarsim) and clone your fork.
2. Create a **branch** from `main` (or `master`, depending on the default branch name in the repo).
3. Create a **virtual environment** (recommended) and install with dev dependencies:

   ```bash
   pip install -e ".[dev]"
   ```

## What to run before opening a PR

From the repository root:

```bash
ruff check src tests
ruff format src tests          # or: ruff format --check src tests to verify only
mypy src
pytest
```

Use `MPLBACKEND=Agg` (or the equivalent on Windows) if your environment has no display or Matplotlib tries to open a GUI backend during tests.

## Code style

- **Python:** 3.11+ syntax and typing. Prefer explicit types on public APIs (`Body`, `SolarSystemSim`, ephemeris helpers).
- **Lint / format:** [Ruff](https://docs.astral.sh/ruff/) configuration lives in [`pyproject.toml`](pyproject.toml).
- **Types:** [Mypy](https://mypy-lang.org/) runs in strict mode on `src/`; third-party stubs are relaxed where needed (see `[tool.mypy.overrides]` in `pyproject.toml`).
- **Architecture:** Keep **simulation and ephemeris** free of `tkinter` and `matplotlib` imports so logic stays testable headless. New UI code belongs under `src/solarsim/gui/`.

## Tests

- Add or update **pytest** tests under [`tests/`](tests/) for behavior you change (physics, ephemeris, timestep policy, collisions, regressions).
- Avoid tests that require a live display or manual interaction unless clearly marked and skipped in CI.

## Pull requests

1. **Scope:** One coherent change per PR when possible (feature, fix, or doc update).
2. **Description:** Summarize motivation, what changed, and how you verified it (commands run).
3. **CI:** Fix any failures on the GitHub Actions workflow (Ruff, Mypy, Pytest on Ubuntu and Windows).

## Packaging and releases

- Version is defined in [`pyproject.toml`](pyproject.toml) under `[project] version`.
- PyInstaller layout is maintained in [`solarsim.spec`](solarsim.spec). If you add optional native dependencies, update hidden imports or hooks as needed and note that in the PR.

## Questions

Open a [GitHub issue](https://github.com/fraware/solarsim/issues) to discuss larger design changes before investing significant implementation time.
