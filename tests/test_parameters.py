import pytest

from mms_viscoelastic.parameters import FlowParameters


def test_default_parameters_are_valid() -> None:
    FlowParameters().validate()


@pytest.mark.parametrize(
    "field,value",
    [
        ("re", 0.0),
        ("re", -1.0),
        ("wi", 0.0),
        ("wi", -1.0),
        ("beta", 0.0),
        ("beta", 1.5),
        ("dt", 0.0),
        ("dt", -1.0e-6),
        ("t_final", 0.0),
        ("t_final", -0.05),
        ("alpha_g", -0.1),
        ("epsilon", -0.1),
    ],
)
def test_validate_rejects_out_of_range_values(field: str, value: float) -> None:
    params = FlowParameters(**{field: value})
    with pytest.raises(ValueError):
        params.validate()


def test_beta_equal_to_one_is_valid() -> None:
    FlowParameters(beta=1.0).validate()
