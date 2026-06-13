"""Constitutive model identifiers used during the Fortran-to-Python port."""

from __future__ import annotations

from enum import Enum


class ConstitutiveModel(str, Enum):
    """Supported model labels."""

    NEWTONIAN = "newtonian"
    UCM = "ucm"
    OLDROYD_B = "oldroyd_b"
    GIESEKUS = "giesekus"
    LPTT = "lptt"
