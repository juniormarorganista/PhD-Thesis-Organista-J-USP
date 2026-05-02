"""Manufactured source-term interfaces.

This subpackage will contain Python ports of the Fortran manufactured source
terms used in the thesis verification code.

The original Fortran sources remain preserved under:

``src/fortran_original/4Order_LPOG_Euler/``

The first target is MMS Case 12, whose original Fortran arrays are:

- ``tfwz``
- ``tfTxx``
- ``tfTxy``
- ``tfTyy``
"""

from mms_viscoelastic.source_terms.case_12 import (
    CASE12_FORTRAN_ARRAYS,
    Case12SourceTermParameters,
    Case12SourceTerms,
    available_case12_source_terms,
    case12_source_term_status,
    evaluate_case_12_source_terms,
)

__all__ = [
    "CASE12_FORTRAN_ARRAYS",
    "Case12SourceTermParameters",
    "Case12SourceTerms",
    "available_case12_source_terms",
    "case12_source_term_status",
    "evaluate_case_12_source_terms",
]
