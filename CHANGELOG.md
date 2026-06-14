# Changelog

Todas as mudanças relevantes deste projeto serão documentadas aqui.

O formato segue [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

### Adicionado
- `numerics/poisson.py` com `solve_poisson_dirichlet`, um solver direto
  (`scipy.sparse` + `spsolve`) para `laplacian(phi) = source` com condições
  de contorno de Dirichlet, usando o operador de cinco pontos de segunda
  ordem.
- `numerics/poisson.py` com `run_case12_poisson_convergence_study`,
  verificando `solve_poisson_dirichlet` contra a MMS do Case 12
  (`laplacian(psi) = wz`, BC = `psi_exact`), com ordem observada ≈ 2 para
  `L2` e `Linf`.
- `examples/12_case12_poisson_solver.py`.
- `tests/test_case12_poisson_solver.py`.
- `docs/case-12-poisson-solver.md`, documentando a discretização, a
  convenção de sinal `wz = laplacian(psi)` e a simplificação de não portar o
  solver multigrid do Fortran (`poisson.f`).
- `scipy>=1.11` adicionado como dependência principal (antes apenas em `ml`).

- Programa Fortran isolado de referência para os campos exatos do Case 12
  (`u`, `v`, `wz`, `psi`, `txx`, `txy`, `tyy`) no instante `t=0`, porta
  verbatim de `exact_solutions_case_12(t)`.
- CSVs de referência 3x3 e 5x5 para os campos exatos do Case 12.
- `validation/solver_comparison.py` com `CASE12_FIELD_COMPONENTS` e
  `compare_case12_field_reference`, comparando `evaluate_case_12_exact`
  contra a referência Fortran.
- `examples/11_case12_solver_comparison.py`.
- `tests/test_case12_solver_comparison.py`, cobrindo os sete campos para as
  referências 3x3 e 5x5 (`max_abs_error < 1e-10`, `max_rel_error < 1e-10`).
- `docs/case-12-solver-comparison.md`, documentando o escopo reduzido da
  FASE VII (comparação dos campos exatos no instante t=0) e o trabalho futuro
  de comparação completa do solver.

## [0.2.0] - 2026-06-13

### Corrigido
- `models/constitutive.py`: `ConstitutiveModel` agora usa `str, Enum` em vez de
  `enum.StrEnum` (incompatível com Python 3.10).
- `ml/features.py`: `field_to_feature_table` agora estreita corretamente o tipo
  de dataclass antes de chamar `asdict`.
- `visualization/fields.py`: `extent` do `imshow` agora é uma tupla, conforme
  esperado pelos stubs do Matplotlib.
- `source_terms/case_12.py`: removida variável local não utilizada em
  `evaluate_case_12_tfTxy`.

### Adicionado
- `models/ucm.py` com `VelocityGradient2D`, `StressTensor2D`,
  `rate_of_deformation_tensor`, `upper_convected_derivative`, `ucm_residual`
  e `ucm_steady_simple_shear`, implementando o modelo constitutivo
  Upper-Convected Maxwell (UCM) adimensional e sua solução analítica de
  cisalhamento simples (Couette).
- `models/oldroyd_b.py` com `solvent_stress`, `total_stress`,
  `oldroyd_b_residual` e `oldroyd_b_steady_simple_shear`, implementando o
  modelo Oldroyd-B como decomposição `tau = tau_p + tau_s` sobre o UCM.
- `docs/constitutive-models.md` com a derivação rigorosa da derivada
  upper-convected em forma de componentes e das soluções analíticas de UCM e
  Oldroyd-B em cisalhamento simples, usadas na validação dos novos modelos.
- `tests/test_model_ucm.py` e `tests/test_model_oldroyd_b.py`, validando os
  modelos UCM e Oldroyd-B contra as soluções analíticas de Couette e
  verificando que o resíduo das equações constitutivas é nulo (`atol=1e-12`)
  para uma grade de `(Wi, beta, shear_rate)`, incluindo os limites
  newtonianos `beta=0` e `beta=1`.
- `ml/splits.py` com `stratified_train_validation_test_split`, dividindo um
  dataframe em train/validation/test estratificado por uma coluna (ex.:
  `case_id`), garantindo que cada split contenha amostras de todas as
  combinações de parâmetros.
- `ml/datasets.py`: `build_case12_parametric_dataset`, que varre combinações
  de `(Re, Wi, beta, t)` via `itertools.product`, combina os campos exatos do
  Case 12 com os quatro termos fonte (`tfwz`, `tfTxx`, `tfTxy`, `tfTyy`) em
  uma única tabela de features por combinação, e retorna um
  `ParametricDatasetBundle` com splits estratificados por `case_id`.
- `ml/datasets.py`: `save_parametric_dataset_bundle`, `ParametricDatasetBundle`,
  `PARAMETRIC_TARGET_COLUMNS` e `SOURCE_TERM_COLUMNS`.
- `infer_feature_columns` agora aceita `exclude_columns`, usado para remover
  identificadores como `case_id` das features sem removê-los do dataframe.
- `examples/10_case12_ml_dataset_parametric.py`, exemplo que gera o dataset
  paramétrico e exporta os CSVs de full/train/validation/test.
- `tests/test_ml_splits.py` e `tests/test_ml_parametric.py`, cobrindo a
  validação das frações do split estratificado, a cobertura de todas as
  combinações de parâmetros nos três splits e a reprodutibilidade via
  `random_state`.
- Job `lint` no CI (`ruff check` + `mypy`) executado antes da matriz de testes.
- `numerics/convergence.py` com `run_convergence_study`, verificando a ordem
  de convergência (L2 e Linf) das relações `u = d(psi)/dy`, `v = -d(psi)/dx`
  e `wz = d(u)/dy - d(v)/dx` para o Case 12.
- `examples/09_case12_convergence_verification.py`.
- `tests/test_case12_convergence.py` (ordem L2 >= 1.8 e Linf >= 1.5 para
  grids 8->16->32->64).
- `docs/case-12-convergence-verification.md` com a verificação simbólica das
  identidades e a análise de erro de truncamento dos operadores.
- `docs/case-12-complete-validation-report.md`, consolidando a solução exata,
  a validação dos quatro termos fonte (tfwz, tfTxx, tfTxy, tfTyy) contra a
  referência Fortran e o estudo de convergência, incluindo a discussão da
  convenção de sinal de `wz` em relação à formulação geral do README.
- `examples/notebooks/Case12_MMS_Pipeline.ipynb`, notebook reprodutível que
  executa a solução exata, os termos fonte, a validação Fortran e o estudo
  de convergência do Case 12.
- Seção "Convergence results (Case 12)" no README, com a tabela de
  convergência e a nota sobre a convenção de sinal de `wz`.

### Alterado
- `source_terms/case_12.py`: docstring do módulo atualizada para refletir que
  os quatro termos fonte do Case 12 (tfwz, tfTxx, tfTxy, tfTyy) já estão
  validados contra a referência Fortran.
- README: status do `tfTyy` atualizado de "pending validation" para
  "validated"; seção "Pending work" atualizada.

---

## [0.6.7.5] - 2026-05-03

### Adicionado
- `assets/thesis/case12/` com figuras selecionadas da tese (campos exatos, superfícies de erro, curvas de convergência, slices numéricos).
- Seção "Thesis figures in this README" no README.
- `docs/thesis-figures-used-in-readme.md` com mapa de origem das figuras.

---

## [0.6.7.4] - 2026-05-03

### Adicionado
- README científico completo com: visão geral, formulação matemática, explicação do MMS, equações de verificação numérica, status do Case 12, validação dos termos fonte, figuras, comandos de instalação/teste/exemplos, fluxo Git local e roadmap.
- SVGs: `assets/preview/mms_validation_pipeline.svg`, `assets/preview/repository_architecture.svg`, `assets/preview/case12_source_terms_status.svg`.
- `docs/readme-maintenance-policy.md` - política de manutenção do README.

### Regra
O README não deve ser reduzido. Fases futuras podem expandir ou melhorar, nunca substituir por versão mínima.

---

## [0.6.7.3] - 2026-05-03

### Corrigido
- Restaurados `mean_abs_error` e `mean_rel_error` em `SourceTermReferenceComparison`.

---

## [0.6.7.2] - 2026-05-03

### Corrigido
- Restaurados `nx` e `ny` em `SourceTermReferenceComparison`.
- Reconstrução do grid de referência agora retorna `(grid, nx, ny)`.

---

## [0.6.7.1] - 2026-05-03

### Corrigido
- Removidos argumentos `x0`, `x1`, `y0`, `y1` da reconstrução do grid de validação (API `Grid2D.from_bounds` não os aceita).
- Adicionada validação explícita de que as coordenadas do CSV Fortran formam grids uniformes em `[0, 1]`.
- Testes de regressão para reconstrução do grid de referência.

---

## [0.6.7] - 2026-05-02

### Adicionado
- `evaluate_case_12_tfTyy` - quarto e último termo fonte do Case 12.
- `Case12SourceTerms` completo com `tfwz`, `tfTxx`, `tfTxy`, `tfTyy`.
- Registro atualizado de termos fonte.
- Testes de shape, valores finitos e referência padrão para `tfTyy`.
- `docs/case-12-tfTyy-port.md`.

### Status
Todos os quatro termos fonte do Case 12 estão implementados em Python.

---

## [0.6.6] - 2026-05-02

### Adicionado
- Programa Fortran isolado de referência para `tfTxy`.
- CSVs de referência 3×3 e 5×5 para `tfTxy`.
- Utilitário de validação Python para `tfTxy`.
- `examples/07_validate_case12_tfTxy_reference.py`.
- Testes de regressão pytest contra dados Fortran para `tfTxy`.
- `docs/case-12-tfTxy-fortran-reference-validation.md`.

---

## [0.6.5] - 2026-05-02

### Adicionado
- `evaluate_case_12_tfTxy`.
- Objeto parcial atualizado com `tfwz`, `tfTxx` e `tfTxy`.

### Corrigido
- Referência padrão de teste para `tfTxx` na origem corrigida para `2.0`.

---

## [0.6.4] - 2026-05-02

### Adicionado
- Programa Fortran isolado de referência para `tfTxx`.
- CSVs de referência 3×3 e 5×5 para `tfTxx`.
- Utilitários de validação Python para comparação CSV de termos fonte.
- `examples/06_validate_case12_tfTxx_reference.py`.
- Testes de regressão pytest contra dados Fortran para `tfTxx`.
- `docs/case-12-tfTxx-fortran-reference-validation.md`.

---

## [0.6.3] - 2026-05-02

### Adicionado
- `evaluate_case_12_tfTxx`.
- Objeto parcial atualizado com `tfwz` e `tfTxx`.

---

## [0.6.2] - 2026-05-02

### Adicionado
- Exportador Fortran isolado para `tfwz`.
- CSVs de referência 3×3 e 5×5 para `tfwz`.
- Utilitário de validação Python para comparar `tfwz` contra CSVs de referência.
- `examples/05_validate_case12_tfwz_reference.py`.
- Testes de regressão pytest contra dados Fortran para `tfwz`.
- `docs/case-12-tfwz-fortran-reference-validation.md`.

---

## [0.6.1] - 2026-05-02

### Adicionado
- `evaluate_case_12_tfwz` - primeiro termo fonte portado do Fortran.
- Suporte parcial de termos fonte via `allow_partial=True`.
- Registro explícito de componentes portados e pendentes.
- `examples/04_case12_tfwz_source_term.py`.
- Testes de shape, valores finitos, comportamento parcial e referência padrão para `tfwz`.
- `docs/case-12-tfwz-port.md`.

---

## [0.6.0] - 2026-05-02

### Adicionado
- `src/mms_viscoelastic/source_terms/__init__.py`
- `src/mms_viscoelastic/source_terms/case_12.py` - contrato público de termos fonte (levanta `NotImplementedError` até que as expressões sejam portadas).
- `tests/test_source_terms_case12_contract.py`.
- `docs/case-12-source-term-porting-map.md`.
- `docs/thesis-equations-overview.md`.

### Decisão de design
`evaluate_case_12_source_terms()` levanta `NotImplementedError` para evitar que zeros ou placeholders entrem silenciosamente em workflows de verificação.

---

## [0.5.4] - 2026-05-02

### Corrigido
- `Grid2D` agora instanciado via `Grid2D.from_bounds(nx=..., ny=...)`.
- `field_to_feature_table` preserva o padrão conservador: coordenadas + campos físicos.
- Features ML engineered são adicionadas apenas quando explicitamente solicitadas.
- Exportação com nomes retrocompatíveis: `plot_field`, `plot_field_panel`.
- Exemplo de visualização fecha figuras Matplotlib após salvar.

---

## [0.5.3] - 2026-05-02

### Adicionado
- `src/mms_viscoelastic/ml/features.py`.
- Separação entre extração de features de baixo nível e construtores de dataset.
- `tests/test_ml_imports.py`.
- `docs/ml-package-structure.md`.

### Corrigido
- Import circular em `mms_viscoelastic.ml`.

---

## [0.5.0] - 2026-05-02

### Adicionado
- README profissional com estrutura visual para o GitHub.
- Seção de preview do Case 12.
- `docs/github-about.md`.
- `docs/local-git-workflow.md`.
- `docs/case-12-preliminary-results.md`.
- Pasta `assets/preview/` com imagens para o README.

---

## [0.1.0] - 2026-05-02

### Adicionado
- Reorganização inicial do repositório.
- Preservação dos códigos Fortran originais em `src/fortran_original/`.
- Preservação dos notebooks Mathematica originais em `src/mathematica_original/`.
- Fundação da biblioteca Python `mms_viscoelastic` (parâmetros, malhas, campos, erros, exportação).
- Testes mínimos da camada Python.
- Documentação inicial do projeto.
