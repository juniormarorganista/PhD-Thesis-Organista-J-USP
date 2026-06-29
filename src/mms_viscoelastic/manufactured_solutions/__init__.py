"""Manufactured solutions for verification studies."""

from .case_01 import Case01Parameters, case_01_feature_names, evaluate_case_01_exact
from .case_02 import Case02Parameters, case_02_feature_names, evaluate_case_02_exact
from .case_03 import Case03Parameters, case_03_feature_names, evaluate_case_03_exact
from .case_04 import Case04Parameters, case_04_feature_names, evaluate_case_04_exact
from .case_05 import Case05Parameters, case_05_feature_names, evaluate_case_05_exact
from .case_06 import Case06Parameters, case_06_feature_names, evaluate_case_06_exact
from .case_12 import (
    Case12Coefficients,
    Case12Parameters,
    case_12_feature_names,
    evaluate_case_12_exact,
)
from .toy import evaluate_toy_solution

__all__ = [
    "Case01Parameters",
    "Case02Parameters",
    "Case03Parameters",
    "Case04Parameters",
    "Case05Parameters",
    "Case06Parameters",
    "Case12Coefficients",
    "Case12Parameters",
    "case_01_feature_names",
    "case_02_feature_names",
    "case_03_feature_names",
    "case_04_feature_names",
    "case_05_feature_names",
    "case_06_feature_names",
    "case_12_feature_names",
    "evaluate_case_01_exact",
    "evaluate_case_02_exact",
    "evaluate_case_03_exact",
    "evaluate_case_04_exact",
    "evaluate_case_05_exact",
    "evaluate_case_06_exact",
    "evaluate_case_12_exact",
    "evaluate_toy_solution",
]
