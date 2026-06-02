# Fase 4 - Visualização do case 12 e datasets para Machine Learning

Esta fase adiciona uma camada prática sobre o port Python do caso manufaturado 12.
O objetivo é transformar a solução exata em dois produtos úteis:

1. figuras de inspeção visual para a tese, apresentação e README;
2. datasets tabulares para futuros experimentos de Machine Learning.

## Arquivos criados

```text
src/mms_viscoelastic/visualization/__init__.py
src/mms_viscoelastic/visualization/fields.py
src/mms_viscoelastic/ml/datasets.py
examples/03_case12_visualization_and_ml_dataset.py
tests/test_case12_visualization.py
tests/test_case12_datasets.py
```

## Convenção de eixos

A função `plot_scalar_field` valida explicitamente que o array tem forma:

```text
(ny, nx)
```

Isso evita a troca acidental de eixos entre a convenção Fortran e a convenção
NumPy/Matplotlib. Esse ponto é crítico para comparar corretamente os resultados
numéricos com as figuras da tese.

## Como gerar figuras e datasets

Na raiz do repositório, após instalar o pacote em modo desenvolvimento:

```powershell
python examples\03_case12_visualization_and_ml_dataset.py
```

Com grade mais fina:

```powershell
python examples\03_case12_visualization_and_ml_dataset.py --nx 129 --ny 129
```

Com pasta de saída personalizada:

```powershell
python examples\03_case12_visualization_and_ml_dataset.py --output-dir outputs\case12
```

## Saídas esperadas

```text
outputs/case12/figures/case12_panel.png
outputs/case12/figures/case12_u.png
outputs/case12/figures/case12_v.png
outputs/case12/figures/case12_psi.png
outputs/case12/figures/case12_wz.png
outputs/case12/figures/case12_txx.png
outputs/case12/figures/case12_txy.png
outputs/case12/figures/case12_tyy.png
outputs/case12/datasets/case12_full.csv
outputs/case12/datasets/case12_train.csv
outputs/case12/datasets/case12_validation.csv
outputs/case12/datasets/case12_metadata.txt
```

## Colunas do dataset

As colunas de alvo são:

```text
u, v, psi, wz, txx, txy, tyy
```

As colunas de entrada incluem:

```text
x, y, time, beta, x2, y2, xy, sin_pi_x, cos_pi_x, sin_pi_y, cos_pi_y
```

Esse formato já permite treinar modelos de regressão supervisionada para campos
manufaturados e, futuramente, comparar modelos clássicos de ML, PINNs e operadores
neurais.

## Limitação

Esta fase ainda não usa dados produzidos pelo Fortran. Ela trabalha com a solução
exata Python do `case_12`. A comparação Fortran vs Python deve vir em fase própria,
após a exportação padronizada dos resultados de referência.
