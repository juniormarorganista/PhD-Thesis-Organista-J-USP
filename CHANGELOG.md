# Changelog

Todas as mudanças relevantes deste projeto serão documentadas aqui.

O formato segue [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

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
- Job `lint` no CI (`ruff check` + `mypy`) executado antes da matriz de testes.
- `numerics/convergence.py` com `run_convergence_study`, verificando a ordem
  de convergência (L2 e Linf) das relações `u = d(psi)/dy`, `v = -d(psi)/dx`
  e `wz = d(u)/dy - d(v)/dx` para o Case 12.
- `examples/09_case12_convergence_verification.py`.
- `tests/test_case12_convergence.py` (ordem L2 >= 1.8 e Linf >= 1.5 para
  grids 8->16->32->64).
- `docs/case-12-convergence-verification.md` com a verificação simbólica das
  identidades e a análise de erro de truncamento dos operadores.

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
