# Contributing to GeoMix

Thank you for your interest in contributing to GeoMix! Contributions of all kinds are welcome — bug reports, feature requests, documentation improvements, and code patches.

## Reporting Bugs

Please open a [GitHub issue](https://github.com/bradleywakefield/geomix-paper/issues) and include:

- A clear, descriptive title.
- Steps to reproduce the issue.
- The Python version and relevant package versions (`pip freeze`).
- The full error traceback, if applicable.

## Suggesting Enhancements

Open a GitHub issue describing:

- The motivation for the enhancement and the use case it addresses.
- A proposed interface or implementation sketch (optional but helpful).

## Pull Requests

1. **Fork** the repository and create a feature branch:
   ```bash
   git checkout -b feature/my-feature
   ```

2. **Install** the development dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. **Make your changes**, keeping commits focused and atomic.

4. **Run the tests** to ensure nothing is broken:
   ```bash
   pytest tests/
   ```

5. **Open a pull request** against the `main` branch with a clear description of what was changed and why.

## Code Style

- Follow [PEP 8](https://peps.python.org/pep-0008/) for Python code.
- Use [ruff](https://docs.astral.sh/ruff/) for linting: `ruff check src/ scripts/ tests/`.
- Add docstrings to public functions and classes (NumPy docstring style preferred).

## Reproducibility

If your contribution affects numerical results or figures, please update the relevant script in `scripts/` and document the change in the pull request description.

## License

By contributing, you agree that your contributions will be licensed under the [MIT License](LICENSE).
