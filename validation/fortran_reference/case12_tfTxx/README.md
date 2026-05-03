# Case 12 `tfTxx` Fortran reference

This directory contains an isolated Fortran reference program for the MMS Case
12 source term:

```text
tfTxx
```

The program is intentionally separated from the original thesis Fortran code.
It does not modify the original files under:

```text
src/fortran_original/
```

## Compile

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case12_tfTxx\case12_tfTxx_reference.f90 `
  -o validation\fortran_reference\case12_tfTxx\case12_tfTxx_reference.exe
```

## Generate reference CSV

3 x 3:

```powershell
.\validation\fortran_reference\case12_tfTxx\case12_tfTxx_reference.exe `
  validation\reference_data\case12_tfTxx_3x3_default.csv 3 3
```

5 x 5:

```powershell
.\validation\fortran_reference\case12_tfTxx\case12_tfTxx_reference.exe `
  validation\reference_data\case12_tfTxx_5x5_default.csv 5 5
```

## CSV schema

```text
x,y,tfTxx
```

The grid uses the same convention as the Python validation layer:

```text
x in [0, 1]
y in [0, 1]
```
