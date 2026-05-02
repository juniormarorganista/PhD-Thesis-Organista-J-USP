"""Physical and numerical parameters used by the Python port."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FlowParameters:
    """Dimensionless parameters for incompressible viscoelastic simulations.

    Parameters
    ----------
    re:
        Reynolds number.
    wi:
        Weissenberg number.
    beta:
        Solvent-viscosity ratio.  In many formulations, 0 < beta <= 1.
    epsilon:
        PTT-type nonlinear parameter.  Keep zero for UCM/Oldroyd-B/Giesekus cases.
    xi:
        Slip/asymmetry parameter used by some PTT-type formulations.
    alpha_g:
        Giesekus mobility parameter.
    dt:
        Time step.
    t_final:
        Final simulation time.
    """

    re: float = 1.0
    wi: float = 1.0
    beta: float = 0.5
    epsilon: float = 0.0
    xi: float = 0.0
    alpha_g: float = 0.0
    dt: float = 1.0e-6
    t_final: float = 0.05

    def validate(self) -> None:
        """Raise ``ValueError`` if parameters are outside basic admissible ranges."""
        if self.re <= 0:
            raise ValueError("re must be positive.")
        if self.wi <= 0:
            raise ValueError("wi must be positive.")
        if not (0 < self.beta <= 1):
            raise ValueError("beta must satisfy 0 < beta <= 1.")
        if self.dt <= 0:
            raise ValueError("dt must be positive.")
        if self.t_final <= 0:
            raise ValueError("t_final must be positive.")
        if self.alpha_g < 0:
            raise ValueError("alpha_g must be non-negative.")
        if self.epsilon < 0:
            raise ValueError("epsilon must be non-negative.")
