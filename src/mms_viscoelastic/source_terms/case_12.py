"""Source terms for MMS Case 12.

This module ports selected algebraic source terms from the original Fortran
routine:

``src/fortran_original/4Order_LPOG_Euler/case_12.f::source_term_case_12(t)``

Current implementation status:

- ``tfwz``: ported to NumPy;
- ``tfTxx``: ported to NumPy;
- ``tfTxy``: ported to NumPy;
- ``tfTyy``: ported to NumPy.

All four Case 12 source-term components are now available in Python and
validated against the independent Fortran reference (see
``docs/case-12-tfwz-fortran-reference-validation.md``,
``docs/case-12-tfTxx-fortran-reference-validation.md``,
``docs/case-12-tfTxy-fortran-reference-validation.md`` and
``docs/case-12-tfTyy-fortran-reference-validation.md``).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from numpy.typing import NDArray

from mms_viscoelastic.grids import Grid2D

FloatArray = NDArray[np.float64]

CASE12_FORTRAN_ARRAYS: tuple[str, ...] = ("tfwz", "tfTxx", "tfTxy", "tfTyy")
CASE12_PORTED_ARRAYS: tuple[str, ...] = ("tfwz", "tfTxx", "tfTxy", "tfTyy")
CASE12_PENDING_ARRAYS: tuple[str, ...] = ()


@dataclass(frozen=True)
class Case12SourceTermParameters:
    """Physical and constitutive parameters required by Case 12 source terms."""

    reynolds: float = 1.0
    wi: float = 1.0
    beta: float = 0.5
    epsilon: float = 0.0
    xi: float = 0.0
    alpha_g: float = 0.0
    at: float = 0.0

    a1xx: float = 1.0
    a2xx: float = 1.0
    a3xx: float = 1.0
    b1xx: float = 1.0
    b2xx: float = 1.0
    b3xx: float = 1.0

    a1xy: float = 1.0
    a2xy: float = 1.0
    a3xy: float = 1.0
    b1xy: float = 1.0
    b2xy: float = 1.0
    b3xy: float = 1.0

    a1yy: float = 1.0
    a2yy: float = 1.0
    a3yy: float = 1.0
    b1yy: float = 1.0
    b2yy: float = 1.0
    b3yy: float = 1.0

    def validate(self) -> None:
        """Validate parameters needed by the ported source terms."""

        if self.reynolds == 0.0:
            raise ValueError("reynolds must be nonzero.")
        if self.wi == 0.0:
            raise ValueError("wi must be nonzero.")
        if not np.isfinite(self.beta):
            raise ValueError("beta must be finite.")


@dataclass(frozen=True)
class Case12SourceTerms:
    """Container for MMS Case 12 source-term arrays.

    Pending components are represented by ``None`` until their full algebraic
    Fortran expressions are ported and validated.
    """

    tfwz: Optional[FloatArray] = None
    tfTxx: Optional[FloatArray] = None
    tfTxy: Optional[FloatArray] = None
    tfTyy: Optional[FloatArray] = None

    def available_components(self) -> tuple[str, ...]:
        """Return source-term components currently stored as arrays."""

        names = ("tfwz", "tfTxx", "tfTxy", "tfTyy")
        return tuple(name for name in names if getattr(self, name) is not None)

    def pending_components(self) -> tuple[str, ...]:
        """Return source-term components not yet stored as arrays."""

        names = ("tfwz", "tfTxx", "tfTxy", "tfTyy")
        return tuple(name for name in names if getattr(self, name) is None)


def available_case12_source_terms() -> tuple[str, ...]:
    """Return the Fortran source-term arrays targeted by the Python port."""

    return CASE12_FORTRAN_ARRAYS


def ported_case12_source_terms() -> tuple[str, ...]:
    """Return source-term arrays already ported to Python."""

    return CASE12_PORTED_ARRAYS


def pending_case12_source_terms() -> tuple[str, ...]:
    """Return source-term arrays still pending."""

    return CASE12_PENDING_ARRAYS


def case12_source_term_status() -> dict[str, str]:
    """Return current implementation status for each source-term component."""

    return {
        "tfwz": "ported to NumPy; regression against Fortran reference available",
        "tfTxx": "ported to NumPy; regression against Fortran reference available",
        "tfTxy": "ported to NumPy; regression against Fortran reference pending",
        "tfTyy": "ported to NumPy; regression against Fortran reference pending",
    }


def _case12_common(
    grid: Grid2D,
    t: float,
    parameters: Case12SourceTermParameters | None,
) -> dict[str, object]:
    """Collect common arrays and parameters for Case 12 source terms."""

    params = parameters or Case12SourceTermParameters()
    params.validate()

    x = grid.X
    y = grid.Y
    pi = np.pi

    cos_px = np.cos(pi * x)
    sin_px = np.sin(pi * x)
    cos_py = np.cos(pi * y)
    sin_py = np.sin(pi * y)
    x_poly = (-1.0 + x) * x

    return {
        "p": params,
        "x": x,
        "y": y,
        "pi": pi,
        "exp_at": np.exp(params.at * t),
        "exp_2at": np.exp(2.0 * params.at * t),
        "cos_px": cos_px,
        "sin_px": sin_px,
        "cos_py": cos_py,
        "sin_py": sin_py,
        "cos_2px": np.cos(2.0 * pi * x),
        "sin_2px": np.sin(2.0 * pi * x),
        "cos_2py": np.cos(2.0 * pi * y),
        "sin_2py": np.sin(2.0 * pi * y),
        "x_poly": x_poly,
        "dpsi_dy_like": pi * y * cos_py + sin_py,
        "dudx_shape": pi * x_poly * cos_px + (-1.0 + 2.0 * x) * sin_px,
        "uxx_shape": (
            2.0 * pi * (1.0 - 2.0 * x) * cos_px
            + (-2.0 + pi**2 * x_poly) * sin_px
        ),
        "mixed_xy_shape": -2.0 * cos_py + pi * y * sin_py,
    }


def _tensor_shapes(
    c: dict[str, object],
) -> dict[str, FloatArray]:
    """Return manufactured tensor shape factors used in Case 12."""

    p: Case12SourceTermParameters = c["p"]  # type: ignore[assignment]
    cos_px: FloatArray = c["cos_px"]  # type: ignore[assignment]
    sin_px: FloatArray = c["sin_px"]  # type: ignore[assignment]
    cos_py: FloatArray = c["cos_py"]  # type: ignore[assignment]
    sin_py: FloatArray = c["sin_py"]  # type: ignore[assignment]

    Axx = p.a1xx + p.a3xx * cos_px + p.a2xx * sin_px
    Bxx = p.b1xx + p.b3xx * cos_py + p.b2xx * sin_py

    Axy = p.a1xy + p.a3xy * cos_px + p.a2xy * sin_px
    Bxy = p.b1xy + p.b3xy * cos_py + p.b2xy * sin_py

    Ayy = p.a1yy + p.a3yy * cos_px + p.a2yy * sin_px
    Byy = p.b1yy + p.b3yy * cos_py + p.b2yy * sin_py

    return {
        "Axx": Axx,
        "Bxx": Bxx,
        "Axy": Axy,
        "Bxy": Bxy,
        "Ayy": Ayy,
        "Byy": Byy,
    }


def evaluate_case_12_tfwz(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case12SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfwz`` source term for MMS Case 12."""

    c = _case12_common(grid, t, parameters)
    p: Case12SourceTermParameters = c["p"]  # type: ignore[assignment]

    x: FloatArray = c["x"]  # type: ignore[assignment]
    y: FloatArray = c["y"]  # type: ignore[assignment]
    pi: float = c["pi"]  # type: ignore[assignment]
    exp_at: FloatArray = c["exp_at"]  # type: ignore[assignment]
    exp_2at: FloatArray = c["exp_2at"]  # type: ignore[assignment]
    cos_px: FloatArray = c["cos_px"]  # type: ignore[assignment]
    sin_px: FloatArray = c["sin_px"]  # type: ignore[assignment]
    cos_py: FloatArray = c["cos_py"]  # type: ignore[assignment]
    sin_py: FloatArray = c["sin_py"]  # type: ignore[assignment]
    cos_2px: FloatArray = c["cos_2px"]  # type: ignore[assignment]
    sin_2px: FloatArray = c["sin_2px"]  # type: ignore[assignment]
    cos_2py: FloatArray = c["cos_2py"]  # type: ignore[assignment]
    sin_2py: FloatArray = c["sin_2py"]  # type: ignore[assignment]
    x_poly: FloatArray = c["x_poly"]  # type: ignore[assignment]

    Rey = p.reynolds
    betann = p.beta
    at = p.at

    term_base = (
        -(Rey * y)
        + 2.0 * Rey * x * y
        + Rey * y * cos_2px
        - 2.0 * Rey * x * y * cos_2px
        + 2.0 * pi**2 * Rey * x * y * cos_2px
        - 6.0 * pi**2 * Rey * x**2 * y * cos_2px
        + 4.0 * pi**2 * Rey * x**3 * y * cos_2px
        + pi * Rey * y * sin_2px
        - 2.0 * pi * Rey * x * y * sin_2px
        + 2.0 * pi * Rey * x**2 * y * sin_2px
        - 2.0 * pi**3 * Rey * x**2 * y * sin_2px
        + 4.0 * pi**3 * Rey * x**3 * y * sin_2px
        - 2.0 * pi**3 * Rey * x**4 * y * sin_2px
    )

    term_cos2y = (
        Rey
        * y
        * cos_2py
        * (
            -(
                (-1.0 + 2.0 * x)
                * (1.0 + 2.0 * pi**2 * x_poly - cos_2px)
            )
            + pi * (-1.0 - 2.0 * x_poly) * sin_2px
        )
    )

    term_cosx = (
        2.0
        * exp_at
        * pi
        * cos_px
        * (
            -(p.a3xy * p.b1xy * (-1.0 + betann) * pi * Rey)
            + pi
            * (
                (p.a2xx * p.b2xx - p.a2yy * p.b2yy) * Rey
                + betann
                * (
                    -8.0
                    - p.a2xx * p.b2xx * Rey
                    + p.a2yy * p.b2yy * Rey
                    + 16.0 * x
                )
            )
            * cos_py
            + (
                (p.a2xx * p.b3xx - p.a2yy * p.b3yy)
                * (-1.0 + betann)
                * pi
                * Rey
                - 2.0 * (4.0 * betann * pi**2 - at * Rey) * (-1.0 + 2.0 * x) * y
            )
            * sin_py
        )
    )

    term_exp = exp_at * (
        -2.0 * p.a2xy * p.b1xy * (-1.0 + betann) * pi**2 * Rey * sin_px
        + 2.0
        * pi
        * cos_py
        * (
            p.a1xy * p.b3xy * (-1.0 + betann) * pi * Rey
            + (
                Rey
                * (
                    -(p.a3xx * p.b2xx * pi)
                    + p.a3yy * p.b2yy * pi
                    + 2.0 * at * x_poly
                )
                + betann
                * (
                    8.0
                    + p.a3xx * p.b2xx * pi * Rey
                    - pi * (p.a3yy * p.b2yy * Rey + 8.0 * pi * x_poly)
                )
            )
            * sin_px
        )
        + 2.0
        * (
            p.a1xy * p.b2xy * (-1.0 + betann) * pi**2 * Rey
            + (
                -(
                    (p.a3xx * p.b3xx - p.a3yy * p.b3yy)
                    * (-1.0 + betann)
                    * pi**2
                    * Rey
                )
                + 2.0
                * (
                    at * Rey
                    + pi**2
                    * (
                        -(at * Rey * x_poly)
                        + 2.0 * betann * (-4.0 + pi**2 * x_poly)
                    )
                )
                * y
            )
            * sin_px
        )
        * sin_py
    )

    term_sin2y = (
        pi
        * Rey
        * (
            (-1.0 + 2.0 * x)
            * (y**2 + x_poly * (-1.0 + 2.0 * pi**2 * y**2))
            + (-1.0 + 2.0 * x) * (x_poly - y**2) * cos_2px
            - pi
            * (((-1.0 + x) ** 2) * x**2 + (-1.0 - 2.0 * x_poly) * y**2)
            * sin_2px
        )
        * sin_2py
    )

    numerator = term_base + term_cos2y + term_cosx + term_exp + term_sin2y
    return numerator / (2.0 * exp_2at * Rey)


def evaluate_case_12_tfTxx(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case12SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTxx`` source term for MMS Case 12."""

    c = _case12_common(grid, t, parameters)
    s = _tensor_shapes(c)
    p: Case12SourceTermParameters = c["p"]  # type: ignore[assignment]

    x: FloatArray = c["x"]  # type: ignore[assignment]
    y: FloatArray = c["y"]  # type: ignore[assignment]
    pi: float = c["pi"]  # type: ignore[assignment]
    exp_at: FloatArray = c["exp_at"]  # type: ignore[assignment]
    exp_2at: FloatArray = c["exp_2at"]  # type: ignore[assignment]
    cos_px: FloatArray = c["cos_px"]  # type: ignore[assignment]
    sin_px: FloatArray = c["sin_px"]  # type: ignore[assignment]
    cos_py: FloatArray = c["cos_py"]  # type: ignore[assignment]
    sin_py: FloatArray = c["sin_py"]  # type: ignore[assignment]
    x_poly: FloatArray = c["x_poly"]  # type: ignore[assignment]
    dpsi_dy_like: FloatArray = c["dpsi_dy_like"]  # type: ignore[assignment]
    dudx_shape: FloatArray = c["dudx_shape"]  # type: ignore[assignment]
    mixed_xy_shape: FloatArray = c["mixed_xy_shape"]  # type: ignore[assignment]
    uxx_shape: FloatArray = c["uxx_shape"]  # type: ignore[assignment]

    Rey = p.reynolds
    Wi = p.wi
    betann = p.beta
    epsylon = p.epsilon
    xi = p.xi
    alphaG = p.alpha_g
    at = p.at

    Axx = s["Axx"]
    Bxx = s["Bxx"]
    Axy = s["Axy"]
    Bxy = s["Bxy"]
    Ayy = s["Ayy"]
    Byy = s["Byy"]

    x_derivative_shape = -(p.a2xx * cos_px) + p.a3xx * sin_px
    y_derivative_shape = p.b2xx * cos_py - p.b3xx * sin_py

    first = 2.0 * (-1.0 + betann) * exp_at * dudx_shape * dpsi_dy_like / Rey
    giesekus = alphaG * (1.0 - betann) * Rey * Wi * (Axx**2 * Bxx**2 + Axy**2 * Bxy**2)
    relaxation = (
        (1.0 - betann)
        * exp_at
        * Axx
        * Bxx
        * (1.0 + (epsylon * Rey * Wi * Axx * Bxx) / exp_at + (epsylon * Rey * Wi * Ayy * Byy) / exp_at)
    )

    wi_inner = (
        at * (-1.0 + betann) * exp_at * Axx * Bxx
        + (-1.0 + betann) * pi * y * cos_py * Axx * dudx_shape * Bxx
        + (-1.0 + betann) * Axx * dudx_shape * sin_py * Bxx
        - (-1.0 + betann) * pi * x_poly * cos_px * Axx * dpsi_dy_like * Bxx
        - (-1.0 + betann) * (-1.0 + x) * sin_px * Axx * dpsi_dy_like * Bxx
        - (-1.0 + betann) * x * sin_px * Axx * dpsi_dy_like * Bxx
        + (-1.0 + betann) * pi * x_poly * sin_px * x_derivative_shape * dpsi_dy_like * Bxx
        + 2.0 * (-1.0 + betann) * Axx * dudx_shape * dpsi_dy_like * Bxx
        + (-1.0 + betann) * pi * y * Axx * dudx_shape * sin_py * y_derivative_shape
        - 2.0 * (-1.0 + betann) * pi * x_poly * sin_px * Axy * Bxy * mixed_xy_shape
        + (1.0 - betann)
        * xi
        * (
            2.0 * Axx * dudx_shape * dpsi_dy_like * Bxx
            + y * Axy * uxx_shape * sin_py * Bxy
            - pi * x_poly * sin_px * Axy * Bxy * mixed_xy_shape
        )
    )

    return (first + giesekus + relaxation + Wi * wi_inner) / exp_2at


def evaluate_case_12_tfTxy(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case12SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTxy`` source term for MMS Case 12."""

    c = _case12_common(grid, t, parameters)
    s = _tensor_shapes(c)
    p: Case12SourceTermParameters = c["p"]  # type: ignore[assignment]

    x: FloatArray = c["x"]  # type: ignore[assignment]
    y: FloatArray = c["y"]  # type: ignore[assignment]
    pi: float = c["pi"]  # type: ignore[assignment]
    exp_at: FloatArray = c["exp_at"]  # type: ignore[assignment]
    exp_2at: FloatArray = c["exp_2at"]  # type: ignore[assignment]
    cos_px: FloatArray = c["cos_px"]  # type: ignore[assignment]
    sin_px: FloatArray = c["sin_px"]  # type: ignore[assignment]
    cos_py: FloatArray = c["cos_py"]  # type: ignore[assignment]
    sin_py: FloatArray = c["sin_py"]  # type: ignore[assignment]
    x_poly: FloatArray = c["x_poly"]  # type: ignore[assignment]
    dpsi_dy_like: FloatArray = c["dpsi_dy_like"]  # type: ignore[assignment]
    dudx_shape: FloatArray = c["dudx_shape"]  # type: ignore[assignment]
    uxx_shape: FloatArray = c["uxx_shape"]  # type: ignore[assignment]

    Rey = p.reynolds
    Wi = p.wi
    betann = p.beta
    epsylon = p.epsilon
    xi = p.xi
    alphaG = p.alpha_g
    at = p.at

    Axx = s["Axx"]
    Bxx = s["Bxx"]
    Axy = s["Axy"]
    Bxy = s["Bxy"]
    Ayy = s["Ayy"]
    Byy = s["Byy"]

    x_derivative_xy = -(p.a2xy * cos_px) + p.a3xy * sin_px
    y_derivative_xy = p.b2xy * cos_py - p.b3xy * sin_py

    shear_shape = pi * x_poly * cos_py * sin_px - y * (pi * (-1.0 + 2.0 * x) * cos_px + sin_px) * sin_py

    tensor_sum = (
        p.a1xx * p.b1xx
        + p.a1yy * p.b1yy
        + cos_py
        * (
            p.a1xx * p.b3xx
            + p.a1yy * p.b3yy
            + (p.a2xx * p.b3xx + p.a2yy * p.b3yy) * sin_px
        )
        + (p.a1xx * p.b2xx + p.a1yy * p.b2yy) * sin_py
        + sin_px
        * (
            p.a2xx * p.b1xx
            + p.a2yy * p.b1yy
            + (p.a2xx * p.b2xx + p.a2yy * p.b2yy) * sin_py
        )
        + cos_px
        * (
            p.a3xx * p.b1xx
            + p.a3yy * p.b1yy
            + (p.a3xx * p.b3xx + p.a3yy * p.b3yy) * cos_py
            + (p.a3xx * p.b2xx + p.a3yy * p.b2yy) * sin_py
        )
    )

    first = 2.0 * (-1.0 + betann) * exp_at * shear_shape / Rey

    giesekus = alphaG * (1.0 - betann) * Rey * Wi * Axy * Bxy * (Axx * Bxx + Ayy * Byy)

    relaxation = (
        (1.0 - betann)
        * exp_at
        * Axy
        * Bxy
        * (1.0 + (epsylon * Rey * Wi * Axx * Bxx) / exp_at + (epsylon * Rey * Wi * Ayy * Byy) / exp_at)
    )

    wi_inner = (
        y * Axx * uxx_shape * sin_py * Bxx
        + at * exp_at * Axy * Bxy
        + pi * y * cos_py * Axy * dudx_shape * Bxy
        + Axy * dudx_shape * sin_py * Bxy
        - pi * x_poly * cos_px * Axy * dpsi_dy_like * Bxy
        - (-1.0 + x) * sin_px * Axy * dpsi_dy_like * Bxy
        - x * sin_px * Axy * dpsi_dy_like * Bxy
        + pi * x_poly * sin_px * x_derivative_xy * dpsi_dy_like * Bxy
        + pi * y * Axy * dudx_shape * sin_py * y_derivative_xy
        + pi * x_poly * sin_px * Ayy * Byy * (2.0 * cos_py - pi * y * sin_py)
        + xi
        * (
            -(pi * x_poly * cos_py * sin_px)
            + y * (pi * (-1.0 + 2.0 * x) * cos_px + sin_px) * sin_py
        )
        * tensor_sum
    )

    return (first + giesekus + relaxation + (-1.0 + betann) * Wi * wi_inner) / exp_2at


def evaluate_case_12_tfTyy(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case12SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTyy`` source term for MMS Case 12.

    This is a NumPy translation of the original Fortran assignment
    ``tfTyy(i,j) = ...`` in ``source_term_case_12(t)``.
    """

    c = _case12_common(grid, t, parameters)
    s = _tensor_shapes(c)
    p: Case12SourceTermParameters = c["p"]  # type: ignore[assignment]

    x: FloatArray = c["x"]  # type: ignore[assignment]
    y: FloatArray = c["y"]  # type: ignore[assignment]
    pi: float = c["pi"]  # type: ignore[assignment]
    exp_at: FloatArray = c["exp_at"]  # type: ignore[assignment]
    exp_2at: FloatArray = c["exp_2at"]  # type: ignore[assignment]
    cos_px: FloatArray = c["cos_px"]  # type: ignore[assignment]
    sin_px: FloatArray = c["sin_px"]  # type: ignore[assignment]
    cos_py: FloatArray = c["cos_py"]  # type: ignore[assignment]
    sin_py: FloatArray = c["sin_py"]  # type: ignore[assignment]
    x_poly: FloatArray = c["x_poly"]  # type: ignore[assignment]
    dpsi_dy_like: FloatArray = c["dpsi_dy_like"]  # type: ignore[assignment]
    dudx_shape: FloatArray = c["dudx_shape"]  # type: ignore[assignment]
    mixed_xy_shape: FloatArray = c["mixed_xy_shape"]  # type: ignore[assignment]
    uxx_shape: FloatArray = c["uxx_shape"]  # type: ignore[assignment]

    Rey = p.reynolds
    Wi = p.wi
    betann = p.beta
    epsylon = p.epsilon
    xi = p.xi
    alphaG = p.alpha_g
    at = p.at

    Axx = s["Axx"]
    Bxx = s["Bxx"]
    Axy = s["Axy"]
    Bxy = s["Bxy"]
    Ayy = s["Ayy"]
    Byy = s["Byy"]

    x_derivative_yy = -(p.a2yy * cos_px) + p.a3yy * sin_px
    y_derivative_yy = p.b2yy * cos_py - p.b3yy * sin_py

    first = -2.0 * (-1.0 + betann) * exp_at * dudx_shape * dpsi_dy_like / Rey

    relaxation = (
        (1.0 - betann)
        * exp_at
        * Ayy
        * Byy
        * (
            1.0
            + (epsylon * Rey * Wi * Axx * Bxx) / exp_at
            + (epsylon * Rey * Wi * Ayy * Byy) / exp_at
        )
    )

    giesekus = alphaG * (1.0 - betann) * Rey * Wi * (
        Axy**2 * Bxy**2 + Ayy**2 * Byy**2
    )

    wi_inner = (
        2.0 * (-1.0 + betann) * y * Axy * uxx_shape * sin_py * Bxy
        + at * (-1.0 + betann) * exp_at * Ayy * Byy
        + (-1.0 + betann) * pi * y * cos_py * Ayy * dudx_shape * Byy
        + (-1.0 + betann) * Ayy * dudx_shape * sin_py * Byy
        - (-1.0 + betann) * pi * x_poly * cos_px * Ayy * dpsi_dy_like * Byy
        - (-1.0 + betann) * (-1.0 + x) * sin_px * Ayy * dpsi_dy_like * Byy
        - (-1.0 + betann) * x * sin_px * Ayy * dpsi_dy_like * Byy
        + (-1.0 + betann) * pi * x_poly * sin_px * x_derivative_yy * dpsi_dy_like * Byy
        - 2.0 * (-1.0 + betann) * Ayy * dudx_shape * dpsi_dy_like * Byy
        + (-1.0 + betann) * pi * y * Ayy * dudx_shape * sin_py * y_derivative_yy
        + (1.0 - betann)
        * xi
        * (
            y * Axy * uxx_shape * sin_py * Bxy
            - 2.0 * Ayy * dudx_shape * dpsi_dy_like * Byy
            - pi * x_poly * sin_px * Axy * Bxy * mixed_xy_shape
        )
    )

    return (first + relaxation + giesekus + Wi * wi_inner) / exp_2at


def evaluate_case_12_source_terms(
    grid: Grid2D | None = None,
    t: float = 0.0,
    parameters: Case12SourceTermParameters | None = None,
    *,
    allow_partial: bool = False,
) -> Case12SourceTerms:
    """Evaluate all currently ported MMS Case 12 source terms.

    All four source-term components are now available in Python. The
    ``allow_partial`` keyword is retained for backward compatibility with
    previous phases; it is no longer required.
    """

    if grid is None:
        raise ValueError("grid is required.")

    return Case12SourceTerms(
        tfwz=evaluate_case_12_tfwz(grid, t=t, parameters=parameters),
        tfTxx=evaluate_case_12_tfTxx(grid, t=t, parameters=parameters),
        tfTxy=evaluate_case_12_tfTxy(grid, t=t, parameters=parameters),
        tfTyy=evaluate_case_12_tfTyy(grid, t=t, parameters=parameters),
    )


__all__ = [
    "CASE12_FORTRAN_ARRAYS",
    "CASE12_PENDING_ARRAYS",
    "CASE12_PORTED_ARRAYS",
    "Case12SourceTermParameters",
    "Case12SourceTerms",
    "available_case12_source_terms",
    "case12_source_term_status",
    "evaluate_case_12_source_terms",
    "evaluate_case_12_tfTxx",
    "evaluate_case_12_tfTxy",
    "evaluate_case_12_tfTyy",
    "evaluate_case_12_tfwz",
    "pending_case12_source_terms",
    "ported_case12_source_terms",
]
