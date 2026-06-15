# Contributing

This repository is being reorganized to preserve the original thesis material and build a modern Python layer.

## Main rules

1. Do not directly modify the original files in `src/fortran_original/` and `src/mathematica_original/`.
2. Do not directly modify the original thesis files in `Thesis/`.
3. Every new Python feature must have a test.
4. Every Fortran → Python transcription must indicate the source Fortran file.
5. Commits should be small and descriptive.

## Branch convention

```text
chore/...   organization and maintenance
docs/...    documentation
feat/...    new features
fix/...     bug fixes
tests/...   tests
slides/...  thesis presentation
```

## Commit convention

Examples:

```text
docs: add Fortran code map
feat: implement manufactured case 12 in Python
tests: add convergence tests for finite differences
chore: reorganize original thesis sources
slides: restructure thesis defense presentation
```
