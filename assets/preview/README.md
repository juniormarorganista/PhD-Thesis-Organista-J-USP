# Imagens do README

Esta pasta contém imagens usadas diretamente no `README.md`.

## Imagens incluídas neste pacote

```text
mms_pipeline.svg
repository_layout.svg
```

## Imagens geradas localmente

Depois de executar:

```powershell
python examples\03_case12_visualization_and_ml_dataset.py
```

copie:

```powershell
Copy-Item outputs\case12\figures\case12_panel.png assets\preview\case12_panel.png -Force
Copy-Item outputs\case12\figures\case12_u.png assets\preview\case12_u.png -Force
Copy-Item outputs\case12\figures\case12_v.png assets\preview\case12_v.png -Force
Copy-Item outputs\case12\figures\case12_psi.png assets\preview\case12_psi.png -Force
Copy-Item outputs\case12\figures\case12_wz.png assets\preview\case12_wz.png -Force
```

Essas imagens aparecerão na galeria do README.
