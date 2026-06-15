# README images

This folder contains images used directly in `README.md`.

## Images included in this package

```text
mms_pipeline.svg
repository_layout.svg
```

## Locally generated images

After running:

```powershell
python examples\03_case12_visualization_and_ml_dataset.py
```

copy:

```powershell
Copy-Item outputs\case12\figures\case12_panel.png assets\preview\case12_panel.png -Force
Copy-Item outputs\case12\figures\case12_u.png assets\preview\case12_u.png -Force
Copy-Item outputs\case12\figures\case12_v.png assets\preview\case12_v.png -Force
Copy-Item outputs\case12\figures\case12_psi.png assets\preview\case12_psi.png -Force
Copy-Item outputs\case12\figures\case12_wz.png assets\preview\case12_wz.png -Force
```

These images will appear in the README gallery.
