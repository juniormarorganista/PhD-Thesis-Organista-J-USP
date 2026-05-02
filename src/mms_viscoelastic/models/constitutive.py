"""Constitutive model identifiers used during the Fortran-to-Python port."""

from __future__ import annotations

from enum import StrEnum


class ConstitutiveModel(StrEnum):
    """Supported model labels."""

    NEWTONIAN = "newtonian"
    UCM = "ucm"
    OLDROYD_B = "oldroyd_b"
    GIESEKUS = "giesekus"
    LPTT = "lptt"
