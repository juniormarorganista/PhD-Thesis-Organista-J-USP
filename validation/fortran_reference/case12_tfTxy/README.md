# Case 12 `tfTxy` Fortran reference

Isolated Fortran reference program for the MMS Case 12 source term `tfTxy`.

It does not modify the original Fortran thesis code under `src/fortran_original/`.

## Compile

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case12_tfTxy\case12_tfTxy_reference.f90 `
  -o validation\fortran_reference\case12_tfTxy\case12_tfTxy_reference.exe
```

## Generate reference CSV

```powershell
.\validation\fortran_reference\case12_tfTxy\case12_tfTxy_reference.exe `
  validation\reference_data\case12_tfTxy_3x3_default.csv 3 3

.\validation\fortran_reference\case12_tfTxy\case12_tfTxy_reference.exe `
  validation\reference_data\case12_tfTxy_5x5_default.csv 5 5
```
