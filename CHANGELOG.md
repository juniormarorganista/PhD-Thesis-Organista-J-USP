# Changelog

Todas as mudanças relevantes deste projeto serão documentadas aqui.

O formato segue [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

As versões `0.1.0` a `0.6.7.5` cobrem a fase inicial de reorganização do
repositório e da documentação, antes da existência de um pacote Python
formal. A partir da criação do `pyproject.toml`, o versionamento semântico
do pacote `mms_viscoelastic` recomeçou em `0.1.0`/`0.2.0`; é por isso que a
numeração não é monotonicamente crescente ao longo do histórico abaixo.

---

## [Unreleased]

### Adicionado
- MMS Case 02 portado de `case_02.f` (segundo dos 13 casos restantes,
  Etapa 14): `manufactured_solutions/case_02.py`
  (`evaluate_case_02_exact`, `Case02Parameters`) e `source_terms/case_02.py`
  (`evaluate_case_02_tfwz`/`tftxx`/`tftxy`/`tftyy`,
  `Case02SourceTermParameters`). `u`,`v`,`wz`,`psi` são idênticos ao Case
  01; o campo de tensões é polinomial quadrático em x,y, reaproveitando
  `Case12Coefficients`. Simplificação algébrica do termo de relaxação
  (`Axx*Bxx+Ayy*Byy`) verificada com SymPy. Os quatro termos fonte da
  referência Fortran foram gerados por um script de tradução mecânica
  (preserva parênteses do Fortran original caractere a caractere), que
  pegou um parêntese extra introduzido na primeira transcrição manual de
  `tfTxy`; a contrapartida Python tinha um `Wi*` faltando em `tfTyy`,
  também corrigido. Validado contra novo programa Fortran de referência
  isolado (`validation/fortran_reference/case02_fields/`) nas 11
  quantidades exportadas, a precisão de máquina. Documentado em
  `docs/case-02-port.md`.
- MMS Case 01 portado de `case_01.f` (primeiro dos 13 casos restantes,
  Etapa 14): `manufactured_solutions/case_01.py`
  (`evaluate_case_01_exact`, `Case01Parameters`) e `source_terms/case_01.py`
  (`evaluate_case_01_tfwz`/`tftxx`/`tftxy`/`tftyy`,
  `Case01SourceTermParameters`). Diferente do Case 12, o campo de tensões
  manufaturado do Case 01 é espacialmente uniforme (`c1=1`, `c2=2`, `c3=3`
  de `par.nn`, escalados por `(1-beta)`); só `u`,`v`,`wz`,`psi` variam com a
  posição. Identidades de função de corrente (`u=dpsi/dy`, `v=-dpsi/dx`,
  `wz=du/dy-dv/dx`) verificadas simbolicamente com SymPy. Validado contra
  novo programa Fortran de referência isolado
  (`validation/fortran_reference/case01_fields/`) nas 11 quantidades
  exportadas (campos + termos fonte), a `max_abs_error < 1e-10` (precisão
  de máquina observada). Documentado em `docs/case-01-port.md`.
- Esquema compacto de 4ª ordem (Padé) e filtro de Lele, portados de
  `derivs.f` (`coef_4a`, `coeff_4a`, `rhsx_4a`, `rhsxx_4a`, `rhsy_4a`,
  `rhsyy_4a`), `solvers.f` (`tridx`, `tridy`) e `loop.f`
  (`rhs_tridfx_wz`/`rhs_tridfy_wz`, `lhs_tridfx`/`lhs_tridfy`), substituindo
  o placeholder de 2ª ordem usado até aqui:
  - `numerics/derivatives.py`: `solve_tridiagonal` (algoritmo de Thomas),
    `compact_derivative_x`/`compact_derivative_y`,
    `compact_second_derivative_x`/`compact_second_derivative_y`;
  - `numerics/filters.py` (novo): `lele_filter_x`/`lele_filter_y`,
    `filter_coefficients`, com `alphaf=0.48` padrão de `par.nn`;
  - verificação simbólica (SymPy) de que o fechamento de contorno da 2ª
    derivada compacta é de 3ª ordem (uma ordem abaixo do interior de 4ª
    ordem), enquanto o fechamento da 1ª derivada é de 4ª ordem completo;
    documentado em `validation/fortran_reference/case_compact_derivatives/README.md`
    e nos docstrings dos operadores afetados;
  - dois novos programas Fortran de referência isolados
    (`validation/fortran_reference/case_compact_derivatives/`,
    `validation/fortran_reference/case_lele_filter/`) e testes de regressão
    comparando o port Python contra eles (`max_abs_error < 1e-10`);
  - `numerics/convergence.py` reapontado para os operadores compactos; o
    estudo de convergência do Case 12 passa a recuperar ordem ~4.0-4.5
    (antes ~2.0), consistente com a tese e o artigo JNNFM.
  - `numerics/poisson.py`, `numerics/time_stepping.py` e
    `numerics/time_marching.py` continuam usando os operadores de 2ª ordem
    deliberadamente (documentado nos próprios programas Fortran de
    referência existentes); portá-los para o esquema compacto exige também
    reescrever esses programas de referência e o solver de Poisson
    pentadiagonal (`band5_poi`/`banbk5_poi`), deixado para uma etapa
    futura.
- `LICENSE` (MIT) na raiz do repositório, cobrindo a biblioteca Python,
  exemplos, testes e documentação novos; o material Fortran/Mathematica
  original em `src/fortran_original/`, `src/mathematica_original/` e
  `Thesis/` mantém os termos de `Thesis/LICENSE` (BSD-3-Clause,
  Copyright (c) 2015, Humberto Lidio Antonelli). Nova seção "License" no
  `README.md` esclarece essa distinção; necessário para o GitHub detectar
  a licença do repositório. O `LICENSE` contém apenas o texto padrão da MIT
  (sem anexos), para que o detector de licenças do GitHub o reconheça como
  "MIT" em vez de "Other".
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

- `numerics/derivatives.py` com `second_derivative_x` e `second_derivative_y`
  (derivada segunda de segunda ordem, com fronteiras zeradas).
- `numerics/time_stepping.py` com `EulerStepParameters`, `wz_time_derivative`,
  `stress_time_derivatives`, `euler_step_wz` e `euler_step_stress`, porta
  verbatim de `subroutine drv` (`drvs.f`) e do passo de Euler explícito de
  `subroutine solver_euler` (`euler_rk.f`) para `wz`, `Txx`, `Txy` e `Tyy`.
  Apenas os pontos internos da malha são atualizados; a reimposição das
  condições de contorno e o recálculo de `psi`/`ux`/`uy` em `t=dt` ficam para
  a Etapa 10.
- Programa Fortran isolado de referência para um passo de Euler do Case 12
  (`validation/fortran_reference/case12_euler_step/`), com gerador de CSV de
  entrada em Python (`generate_input_csv.py`) a partir dos campos exatos e
  termos fonte já validados, evitando re-transcrever `tfTxx`/`tfTxy`/`tfTyy`.
- CSV de referência 9x9 para o passo de Euler do Case 12
  (`dt = 1e-5`).
- `validation/solver_comparison.py` com `CASE12_EULER_STEP_COMPONENTS` e
  `compare_case12_euler_step_reference`.
- `examples/13_case12_euler_step_comparison.py`.
- `tests/test_case12_euler_step.py`, cobrindo `wz`, `Txx`, `Txy`, `Tyy`
  (`max_abs_error < 1e-12`, `max_rel_error < 1e-12`).
- `docs/case-12-euler-step-comparison.md`, documentando a formulação do passo
  de Euler, a escolha de parâmetros (`Re=1`, `Wi=1`, `beta=0.5`, `dt=1e-5`) e
  o escopo (discretização de segunda ordem, atualização apenas dos pontos
  internos).

- `numerics/time_marching.py` com `TimeMarchingState` e `case12_time_step`,
  combinando o passo de Euler (Etapa 9) com a reimposição das condições de
  contorno de Dirichlet em `t+dt`, a nova resolução de
  `laplacian(psi) = wz` (Etapa 8) e o recálculo de `u`/`v` a partir do `psi`
  atualizado. A condição de contorno de `wz` é reimposta a partir da solução
  exata em `t+dt` (em vez do stencil unilateral de 4ª ordem de
  `boundarie_vorticity_cavity_flow.f`), simplificação documentada e
  consistente com as Etapas 8 e 9.
- Programa Fortran isolado de referência para um passo completo de marcha no
  tempo do Case 12 (`validation/fortran_reference/case12_time_step/`):
  passo de Euler + reimposição de condições de contorno em `t+dt` + resolução
  densa do sistema de Poisson (eliminação gaussiana) + recálculo de `u`/`v`.
- Gerador de CSV de entrada com a coluna `psi`
  (`validation/fortran_reference/case12_time_step/generate_input_csv.py`) e
  CSV de referência 9x9 para o passo completo do Case 12 (`dt = 1e-5`).
- `validation/solver_comparison.py` com `CASE12_TIME_STEP_COMPONENTS` e
  `compare_case12_time_step_reference`.
- `examples/14_case12_time_step_comparison.py`.
- `tests/test_case12_time_step.py`, cobrindo `wz`, `Txx`, `Txy`, `Tyy`, `psi`,
  `u`, `v` (`max_abs_error < 1e-12`, `max_rel_error < 1e-12`).
- `docs/case-12-time-marching-comparison.md`, documentando a sequência
  portada (Euler + condições de contorno + Poisson + `u`/`v`), a decisão de
  reimpor `wz` na fronteira a partir da solução exata e o trabalho futuro de
  marcha multi-passo.

- `numerics/time_marching.py` com `case12_multistep_march`, encadeando
  `case12_time_step` (Etapa 10) por `n_steps` incrementos de `dt`.
- Programa Fortran isolado de referência para a marcha multi-passo do Case 12
  (`validation/fortran_reference/case12_multistep/`), reaproveitando o CSV de
  entrada `t=0` da Etapa 10 a cada passo (válido para o parâmetro padrão
  `at=0`, sob o qual os campos exatos e os termos fonte manufaturados são
  independentes do tempo).
- CSV de referência 9x9 para 10 passos (`dt = 1e-5`).
- `validation/solver_comparison.py` com `compare_case12_multistep_reference`.
- `examples/15_case12_multistep_comparison.py`.
- `tests/test_case12_multistep_march.py`, cobrindo `wz`, `Txx`, `Txy`, `Tyy`,
  `psi`, `u`, `v` após 10 passos (`max_abs_error < 1e-12`,
  `max_rel_error < 1e-12`).
- `docs/case-12-multistep-marching-comparison.md`, documentando
  `case12_multistep_march`, a independência temporal dos campos/termos fonte
  do Case 12 para `at=0` e o trabalho futuro de comparação da evolução
  temporal completa.
- Programa Fortran isolado de referência para a marcha completa do Case 12
  até o critério de convergência de `solver_euler`
  (`validation/fortran_reference/case12_full_march/`), com detecção de
  divergência (`NaN`/overflow) e histórico opcional por passo
  (`errowz`/`erroTxx`/`erroTxy`/`erroTyy`/`emv` + máximos absolutos dos
  campos).
- CSVs de referência de entrada e histórico para as malhas `9x9` e `17x17`
  (`validation/reference_data/case12_full_march_*`).
- `docs/case-12-full-march-reference.md` (Etapa 12 da FASE IX), documentando
  o achado de que a referência de 2ª ordem diverge em `t~4.5506` (malha
  `9x9`) em vez de convergir, e o re-escopo das Etapas 13-15.

### Corrigido
- `CITATION.cff`: `repository-code` apontava para um placeholder
  (`SEU_USUARIO/SEU_REPOSITORIO`); corrigido para a URL real do repositório.
- `CONTRIBUTING.md`: regras de preservação referenciavam `thesis/` (minúsculo,
  inexistente); corrigido para `Thesis/` e para incluir
  `src/mathematica_original/`.

### Alterado
- `README.md`: seção "Repository structure" corrigida para refletir os
  diretórios reais (`Thesis/`, `Papers/<...>/`, e os subpacotes `models/` e
  `io/` de `src/mms_viscoelastic`).
- `README.md`: "Current Python status" e "Pending work" atualizados para
  refletir os modelos constitutivos (UCM, Oldroyd-B) e o solver de marcha no
  tempo do Case 12 (Etapas 8-11).
- `README.md`: nova seção "Time-marching solver (Case 12)" com uma tabela
  resumindo o Poisson solver, o passo de Euler, o passo completo e a marcha
  multi-passo, com links para a documentação de cada etapa.
- `README.md`: seção "Local Git workflow" reescrita como "GitHub workflow",
  refletindo o fluxo real baseado em pull requests
  (`main → branch → commits → push → PR → merge via GitHub → main`).
- `README.md`: seção "GitHub About" atualizada para refletir a descrição e os
  tópicos já configurados no repositório, em vez de apresentá-los como
  sugestão.
- `README.md`: seção "Citation" atualizada para referenciar o `CITATION.cff`
  existente e o recurso "Cite this repository" do GitHub.
- `docs/local-git-workflow.md` e `docs/github-about.md` atualizados para o
  mesmo fluxo baseado em PR e para os tópicos/descrição configurados.

### Removido
- 22 artefatos de build (`__pycache__/*.pyc`) versionados por engano,
  já cobertos por `.gitignore`.

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
