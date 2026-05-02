from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import evaluate_toy_solution
from mms_viscoelastic.ml import field_to_feature_table
from mms_viscoelastic.numerics import derivative_x, l2_error


def main() -> None:
    grid = Grid2D.from_bounds(nx=32, ny=32)
    field = evaluate_toy_solution(grid)
    table = field_to_feature_table(grid, field)

    dudx = derivative_x(field.u, grid.dx)
    zero_error = l2_error(field.u, field.u)

    print(table.head())
    print(f"dudx shape: {dudx.shape}")
    print(f"zero error: {zero_error:.3e}")


if __name__ == "__main__":
    main()
