# Thesis figures used in the README

This document lists the thesis figures selected for the GitHub README.

The original thesis figures are mostly EPS files. GitHub does not render EPS reliably, so the selected figures were converted to PNG and placed under:

```text
assets/thesis/case12/
```

## Selected figures

| README asset | Description | Original source pattern |
|---|---|---|
| `assets/thesis/case12/thesis_streamlines.png` | Streamlines | `Thesis/Slides/AAUgraphics/streamslines.png` |
| `assets/thesis/case12/exact_map_psi.png` | Mapa exato da função corrente psi | `Thesis/Slides/Figures/Exact_Map_*_MMS_12_Psi.eps` |
| `assets/thesis/case12/exact_map_u.png` | Mapa exato da velocidade u | `Thesis/Slides/Figures/Exact_Map_*_MMS_12_U.eps` |
| `assets/thesis/case12/exact_map_v.png` | Mapa exato da velocidade v | `Thesis/Slides/Figures/Exact_Map_*_MMS_12_V.eps` |
| `assets/thesis/case12/exact_map_wz.png` | Mapa exato da vorticidade omega_z | `Thesis/Slides/Figures/Exact_Map_*_MMS_12_Wz.eps` |
| `assets/thesis/case12/exact_map_txx.png` | Mapa exato de Txx | `Thesis/Slides/Figures/Exact_Map_*_MMS_12_Txx.eps` |
| `assets/thesis/case12/exact_map_txy.png` | Mapa exato de Txy | `Thesis/Slides/Figures/Exact_Map_*_MMS_12_Txy.eps` |
| `assets/thesis/case12/exact_map_tyy.png` | Mapa exato de Tyy | `Thesis/Slides/Figures/Exact_Map_*_MMS_12_Tyy.eps` |
| `assets/thesis/case12/exact_surface_psi.png` | Superfície exata de psi | `Thesis/Slides/Figures/Exact_Surf_*_MMS_12_Psi.eps` |
| `assets/thesis/case12/exact_surface_u.png` | Superfície exata de u | `Thesis/Slides/Figures/Exact_Surf_*_MMS_12_U.eps` |
| `assets/thesis/case12/exact_surface_wz.png` | Superfície exata de omega_z | `Thesis/Slides/Figures/Exact_Surf_*_MMS_12_Wz.eps` |
| `assets/thesis/case12/exact_surface_txx.png` | Superfície exata de Txx | `Thesis/Slides/Figures/Exact_Surf_*_MMS_12_Txx.eps` |
| `assets/thesis/case12/convergence_re100_wi1_alphaG0_psi.png` | Convergência de psi | `Thesis/Slides/Figures/NormErr_2nd_Re_100_Wi_1_*_alphaG_0_*_MMS_12_Psi.eps` |
| `assets/thesis/case12/convergence_re100_wi1_alphaG0_u.png` | Convergência de u | `Thesis/Slides/Figures/NormErr_2nd_Re_100_Wi_1_*_alphaG_0_*_MMS_12_U.eps` |
| `assets/thesis/case12/convergence_re100_wi1_alphaG0_v.png` | Convergência de v | `Thesis/Slides/Figures/NormErr_2nd_Re_100_Wi_1_*_alphaG_0_*_MMS_12_V.eps` |
| `assets/thesis/case12/convergence_re100_wi1_alphaG0_wz.png` | Convergência de omega_z | `Thesis/Slides/Figures/NormErr_2nd_Re_100_Wi_1_*_alphaG_0_*_MMS_12_Wz.eps` |
| `assets/thesis/case12/convergence_re100_wi1_alphaG0_txx.png` | Convergência de Txx | `Thesis/Slides/Figures/NormErr_2nd_Re_100_Wi_1_*_alphaG_0_*_MMS_12_Txx.eps` |
| `assets/thesis/case12/convergence_re100_wi1_alphaG0_txy.png` | Convergência de Txy | `Thesis/Slides/Figures/NormErr_2nd_Re_100_Wi_1_*_alphaG_0_*_MMS_12_Txy.eps` |
| `assets/thesis/case12/convergence_re100_wi1_alphaG0_tyy.png` | Convergência de Tyy | `Thesis/Slides/Figures/NormErr_2nd_Re_100_Wi_1_*_alphaG_0_*_MMS_12_Tyy.eps` |
| `assets/thesis/case12/slice_x_u.png` | Corte em x para u | `Thesis/Slides/Figures/Slice_x_Tog_Numerical_*_MMS_12_*_U.eps` |
| `assets/thesis/case12/slice_y_wz.png` | Corte em y para omega_z | `Thesis/Slides/Figures/Slice_y_Tog_Numerical_*_MMS_12_*_Wz.eps` |


## Maintenance rule

Do not commit the full thesis figure tree into `assets/`. Keep only curated README figures here.

If a new thesis figure is needed:

1. copy or convert it to PNG;
2. place it under `assets/thesis/case12/` or another clear subfolder;
3. add a row to this document;
4. reference it in `README.md`.

Recommended conversion command for EPS files:

```powershell
magick -density 140 input.eps -background white -alpha remove -alpha off -trim +repage -resize "1400x1000>" output.png
```
