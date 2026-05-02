# Release Notes — v0.6.0

## Summary

Introduces the source-term porting layer for MMS Case 12.

This release does not fake the long algebraic expressions from Fortran. Instead,
it creates the public source-term contract, documents the mapping from Fortran
arrays to Python arrays, and adds tests that prevent accidental use of
unimplemented source terms.

## Added

- `src/mms_viscoelastic/source_terms/__init__.py`
- `src/mms_viscoelastic/source_terms/case_12.py`
- `tests/test_source_terms_case12_contract.py`
- `docs/case-12-source-term-porting-map.md`
- `docs/thesis-equations-overview.md`

## Design decision

`evaluate_case_12_source_terms()` raises `NotImplementedError` until the exact
Fortran expressions are ported and validated.

This avoids silently returning zeros or placeholders in verification workflows.

## Next release

`v0.6.1` should port `tfwz` first and validate it against selected Fortran
reference values.
