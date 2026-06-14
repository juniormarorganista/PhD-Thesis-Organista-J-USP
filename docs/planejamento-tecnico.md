# Planejamento Técnico - mms_viscoelastic

> Documento de referência estratégica para o desenvolvimento da biblioteca.
> Atualizado em: 2026-06-13.

---

## Estado atual do projeto

| Item | Estado |
|---|---|
| Testes | ✅ 363/363 passam |
| Local ↔ GitHub | ✅ Sincronizado |
| PRs abertos | - nenhum |
| Fundação Python | ✅ `grids`, `fields`, `parameters`, `numerics`, `ml`, `io`, `visualization` |
| Case 12 - solução manufaturada | ✅ portada e testada |
| Case 12 - termos fonte (tfwz, tfTxx, tfTxy, tfTyy) | ✅ portados e validados vs Fortran |
| Módulo de convergência | ✅ `numerics/convergence.py` com `run_convergence_study`, ordem L2 ≈ 2 verificada |
| Dataset ML paramétrico | ✅ `ml/datasets.py` com `build_case12_parametric_dataset` e `ml/splits.py` com split estratificado por `case_id` |
| Modelos constitutivos | ✅ `models/ucm.py` (UCM) e `models/oldroyd_b.py` (Oldroyd-B), validados contra solução analítica de cisalhamento simples (Couette) |
| Relatório científico do Case 12 | ✅ `docs/case-12-complete-validation-report.md` |
| Notebook científico | ✅ `examples/notebooks/Case12_MMS_Pipeline.ipynb` |
| CI/CD (GitHub Actions) | ✅ pytest em 3.10/3.11/3.12 + lint (ruff + mypy) |
| Versão do pacote | [0.2.0](https://github.com/juniormarorganista/Thesis_USP_Codes_Organista_J/tree/v0.2.0) |
| Comparação Fortran x Python - campos exatos do Case 12 (t=0) | ✅ `validation/solver_comparison.py`, ver `docs/case-12-solver-comparison.md` |

---

## Habilidades e seus contextos de uso

| Skill | Quando usar |
|---|---|
| `dev` | Toda tarefa que cria/edita/deleta arquivos Python, roda testes, executa comandos de sistema. Solicitar todas as permissões no início da sessão. |
| `doutor-em-matematica` | Derivações, demonstrações, análises de convergência, documentação de equações. Todo resultado deve ser respaldado por referência ou demonstração. Nenhuma afirmação matemática sem fundamento. |
| `escritor` | Notebooks científicos, relatórios, seções do README, documentação técnica de alto nível. Escrita natural sem padrões detectáveis de IA. |
| `git-github` | Criação de branches, commits, PRs, merges, tags. Seguir prefixos: `feat/`, `fix/`, `docs/`, `chore/`, `test/`. PRs e merges sempre pelo GitHub. |

---

## Modelo por tipo de tarefa

| Modelo | Uso típico |
|---|---|
| **haiku** | Operações mecânicas e bem definidas: git, linting, criação de arquivos de teste a partir de padrões existentes, expansão de datasets ML seguindo padrões, atualização de configurações. |
| **sonnet** | Raciocínio e complexidade: derivações matemáticas, ports de novos termos fonte, análise de convergência, implementação de física nova (modelos constitutivos), escrita científica de qualidade, refatoração multi-arquivo. |

---

## Fases de desenvolvimento

### ▸ FASE I - Organização e baseline de qualidade ✅ CONCLUÍDA
**Modelo: haiku** | **Skills: `dev` + `git-github`**

**Objetivo:** garantir que o projeto está limpo, lint-free e com CI sólido antes de qualquer expansão.

| Tarefa | Arquivo / Comando | Branch |
|---|---|---|
| Rodar ruff e corrigir avisos | `ruff check src/ tests/` | `chore/lint-and-type-check` |
| Rodar mypy e corrigir type errors | `mypy src/` | `chore/lint-and-type-check` |
| Verificar se GitHub Actions roda corretamente em push | `.github/workflows/tests.yml` | `chore/lint-and-type-check` |
| Adicionar ruff ao CI | `.github/workflows/tests.yml` | `chore/lint-and-type-check` |
| PR + merge via GitHub | - | - |

**Status:** [PR #7](https://github.com/juniormarorganista/Thesis_USP_Codes_Organista_J/pull/7) mergeada. `ruff check` e `mypy src/` limpos; job `lint` adicionado ao CI.

---

### ▸ FASE II - Verificação de convergência numérica (núcleo do MMS) ✅ CONCLUÍDA
**Modelo: sonnet** | **Skills: `doutor-em-matematica` + `dev`**

**Objetivo:** implementar e demonstrar que os operadores de diferenças finitas convergem com a ordem correta quando aplicados às soluções manufaturadas do Case 12.

**Contexto matemático:** A verificação MMS exige que o erro numérico $\|u_h - u_{exata}\|_{L_2}$ decresça como $\mathcal{O}(h^p)$ ao refinar a malha, onde $p$ é a ordem nominal do esquema. Os operadores de 2ª ordem em `derivatives.py` devem exibir $p \approx 2$.

| Tarefa | Arquivo | Branch |
|---|---|---|
| Criar `numerics/convergence.py` com `run_convergence_study(case, grids)` | `src/mms_viscoelastic/numerics/convergence.py` | `feat/case12-convergence` |
| Implementar tabela de convergência (h, L2, Linf, ordem estimada) | idem | idem |
| Criar `examples/09_case12_convergence_verification.py` | `examples/09_case12_convergence_verification.py` | idem |
| Testes: verificar ordem >= 1.8 para grids 8→16→32→64 | `tests/test_case12_convergence.py` | idem |
| Documentar equações com rigor (skill `doutor-em-matematica`) | `docs/case-12-convergence-verification.md` | idem |
| PR + merge via GitHub | - | - |

**Status:** [PR #8](https://github.com/juniormarorganista/Thesis_USP_Codes_Organista_J/pull/8) mergeada. Identidades `u = d(psi)/dy`, `v = -d(psi)/dx`, `wz = d(u)/dy - d(v)/dx` verificadas simbolicamente (SymPy) e numericamente: ordem L2 ≈ 2.00-2.12 para grids 8→16→32→64.

---

### ▸ FASE III - Relatório científico e notebook do Case 12 ✅ CONCLUÍDA
**Modelo: sonnet** | **Skills: `escritor` + `doutor-em-matematica` + `git-github`**

**Objetivo:** consolidar toda a evidência científica do Case 12 em um documento e notebook reprodutíveis.

| Tarefa | Arquivo | Branch |
|---|---|---|
| Criar relatório completo de validação | `docs/case-12-complete-validation-report.md` | `docs/case12-scientific-report` |
| Criar notebook científico reprodutível | `examples/notebooks/Case12_MMS_Pipeline.ipynb` | idem |
| Atualizar README com seção "Convergence results" | `README.md` | idem |
| PR + merge via GitHub | - | - |

> **Regra do `escritor`:** redigir na voz do pesquisador. Evitar fraseado genérico de IA. Citar a tese e os papers originais quando mencionar equações.

**Status:** [PR #10](https://github.com/juniormarorganista/Thesis_USP_Codes_Organista_J/pull/10) mergeada. Relatório consolidado em `docs/case-12-complete-validation-report.md` (solução exata, validação dos 4 termos fonte vs. Fortran, estudo de convergência e nota sobre a convenção de sinal de `wz`); notebook `examples/notebooks/Case12_MMS_Pipeline.ipynb` executado com saídas reais; README atualizado com a seção "Convergence results (Case 12)" e status do `tfTyy` corrigido para validado.

---

### ▸ FASE IV - Expansão dos datasets ML ✅ CONCLUÍDA
**Modelo: haiku** | **Skills: `dev` + `git-github`**

**Objetivo:** ampliar o pipeline de datasets para cobrir múltiplos regimes de parâmetros, habilitando treinamento de modelos surrogates.

| Tarefa | Arquivo | Branch |
|---|---|---|
| Expandir `ml/datasets.py` para varredura paramétrica (Re, Wi, beta, t) | `src/mms_viscoelastic/ml/datasets.py` | `feat/ml-parametric-dataset` |
| Criar `ml/splits.py` com train/val/test stratificado | `src/mms_viscoelastic/ml/splits.py` | idem |
| Criar `examples/10_case12_ml_dataset_parametric.py` | `examples/10_case12_ml_dataset_parametric.py` | idem |
| Testes para datasets e splits | `tests/test_ml_parametric.py`, `tests/test_ml_splits.py` | idem |
| PR + merge via GitHub | - | - |

**Status:** [PR #12](https://github.com/juniormarorganista/Thesis_USP_Codes_Organista_J/pull/12) mergeada. `build_case12_parametric_dataset` varre combinações de `(Re, Wi, beta, t)`, combina os campos exatos do Case 12 com os quatro termos fonte em uma única tabela por combinação e retorna um `ParametricDatasetBundle`; `ml/splits.py` adiciona `stratified_train_validation_test_split`, garantindo que train/validation/test contenham amostras de cada `case_id`.

---

### ▸ FASE V - Modelos constitutivos (UCM e Oldroyd-B) ✅ CONCLUÍDA
**Modelo: sonnet** | **Skills: `doutor-em-matematica` + `dev`**

**Objetivo:** implementar os primeiros modelos constitutivos com física real, validados contra soluções analíticas de cisalhamento simples (Couette).

**Contexto matemático:**
- **UCM** (Upper-Convected Maxwell): $\lambda \overset{\nabla}{\mathbf{T}} + \mathbf{T} = 2\eta_p \mathbf{D}$
- **Oldroyd-B**: $\mathbf{T} = \mathbf{T}_p + \mathbf{T}_s$ com $\mathbf{T}_s = 2\eta_s \mathbf{D}$

| Tarefa | Arquivo | Branch |
|---|---|---|
| Implementar `models/ucm.py` com extra-stress UCM | `src/mms_viscoelastic/models/ucm.py` | `feat/model-ucm-oldroyd-b` |
| Implementar `models/oldroyd_b.py` | `src/mms_viscoelastic/models/oldroyd_b.py` | idem |
| Validar contra solução analítica de Couette | `tests/test_model_ucm.py`, `tests/test_model_oldroyd_b.py` | idem |
| Documentar equações com rigor (skill `doutor-em-matematica`) | `docs/constitutive-models.md` | idem |
| PR + merge via GitHub | - | - |

**Status:** [PR #14](https://github.com/juniormarorganista/Thesis_USP_Codes_Organista_J/pull/14) mergeada. `models/ucm.py` implementa a derivada upper-convected em forma de componentes, o resíduo `ucm_residual` e a solução analítica `ucm_steady_simple_shear`; `models/oldroyd_b.py` decompõe `tau = tau_p + tau_s` com `solvent_stress`, `total_stress`, `oldroyd_b_residual` e `oldroyd_b_steady_simple_shear`. `docs/constitutive-models.md` documenta a derivação completa, e os testes verificam as soluções de Couette e o resíduo nulo (`atol=1e-12`) para uma grade de `(Wi, beta, shear_rate)`, incluindo `beta=0` e `beta=1`.

---

### ▸ FASE VI - Versioning e CI ✅ CONCLUÍDA
**Modelo: haiku** | **Skills: `dev` + `git-github`**

**Objetivo:** versionar o primeiro ciclo completo de desenvolvimento e atualizar o CI.

| Tarefa | Arquivo | Branch |
|---|---|---|
| Atualizar versão para `0.2.0` em `_version.py` e `pyproject.toml` | ambos | `chore/release-0.2.0` |
| Adicionar CHANGELOG entry para 0.2.0 | `CHANGELOG.md` | idem |
| Adicionar `ruff` e `mypy` ao GitHub Actions | `.github/workflows/tests.yml` | idem |
| Criar tag `v0.2.0` no GitHub | - | - |
| PR + merge via GitHub | - | - |

**Status:** [PR #16](https://github.com/juniormarorganista/Thesis_USP_Codes_Organista_J/pull/16) mergeada. Versão `0.2.0` em `pyproject.toml`, `_version.py` e `CITATION.cff`; entradas de `[Unreleased]` movidas para `[0.2.0] - 2026-06-13` no `CHANGELOG.md`; `ruff` e `mypy` já estavam habilitados no job `lint` do CI desde a Fase I. Tag [`v0.2.0`](https://github.com/juniormarorganista/Thesis_USP_Codes_Organista_J/tree/v0.2.0) criada.

---

### ▸ FASE VII - Comparação solver-level Fortran × Python ✅ CONCLUÍDA (escopo reduzido)
**Modelo: sonnet** | **Skills: `doutor-em-matematica` + `dev` + `escritor`**

**Objetivo:** demonstrar equivalência numérica entre o solver Fortran original e a implementação Python para o Case 12.

| Tarefa | Arquivo | Branch |
|---|---|---|
| Criar pipeline de comparação campo a campo | `src/mms_viscoelastic/validation/solver_comparison.py` | `feat/solver-comparison` |
| Executar Fortran reference e capturar outputs | `validation/fortran_reference/case12_fields/` | idem |
| Comparar Python vs Fortran nos campos u, v, wz, psi, T | `examples/11_case12_solver_comparison.py` | idem |
| Documentar divergências e justificar tolerâncias | `docs/case-12-solver-comparison.md` | idem |
| PR + merge via GitHub | - | - |

**Status:** [PR #18](https://github.com/juniormarorganista/Thesis_USP_Codes_Organista_J/pull/18) mergeada, com **escopo reduzido**: o pacote `mms_viscoelastic` ainda não possui um solver time-marching (sem solver de Poisson para `psi`, sem integração temporal da vorticidade e do tensor extra-tensão), então a comparação solver-level completa não era viável nesta fase. Em vez disso, foi implementada a comparação dos **campos exatos do Case 12 no instante inicial t=0**: um programa Fortran isolado (`validation/fortran_reference/case12_fields/case12_fields_reference.f90`) porta verbatim `exact_solutions_case_12(t)` e é comparado contra `evaluate_case_12_exact` para os sete campos (`u, v, wz, psi, txx, txy, tyy`), com erro de arredondamento de ponto flutuante (~1e-16 a 1e-21). A comparação completa do solver time-marching fica para a Fase VIII.

---

### ▸ FASE VIII - Solver time-marching Python e comparação completa com Fortran (EM ANDAMENTO)
**Modelo: sonnet** | **Skills: `doutor-em-matematica` + `dev` + `escritor`**

**Objetivo:** portar o laço de integração temporal do Case 12 (`src/fortran_original/4Order_LPOG_Euler/`) para Python e comparar a evolução temporal dos campos `u, v, wz, psi, T` contra o solver Fortran original.

**Contexto:** esta fase foi adiada a partir da FASE VII (ver status acima). Requer, no mínimo, um solver de Poisson para `psi`, a integração temporal da equação de transporte da vorticidade e do tensor extra-tensão, e a injeção dos termos fonte manufaturados (`tfwz`, `tfTxx`, `tfTxy`, `tfTyy`) já portados. É um esforço substancialmente maior que as fases anteriores e foi dividida em sub-etapas.

| Tarefa | Arquivo | Branch | Status |
|---|---|---|---|
| Etapa 8 - solver de Poisson para `psi` (`solve_poisson_dirichlet`) | `src/mms_viscoelastic/numerics/poisson.py`, `docs/case-12-poisson-solver.md` | `feat/case12-poisson-solver` | ✅ PR #20 |
| Etapa 9 - passo de Euler explícito (`wz`, `Txx`, `Txy`, `Tyy`) vs Fortran | `src/mms_viscoelastic/numerics/time_stepping.py`, `docs/case-12-euler-step-comparison.md` | `feat/case12-euler-step` | ✅ PR #21 |
| Etapa 10 - um passo completo de marcha no tempo (Euler + BC + Poisson + `u`,`v`) vs Fortran | `src/mms_viscoelastic/numerics/time_marching.py`, `docs/case-12-time-marching-comparison.md` | `feat/case12-time-marching` | em andamento |
| Etapa 11 (futuro) - marcha multi-passo e comparação da evolução temporal completa | `examples/`, `docs/case-12-solver-comparison.md` | a definir | pendente |

---

## Fluxo de trabalho padrão

```
main → branch → commits → push → PR no GitHub → merge via GitHub → main
```

### Convenções de commit

```
feat: nova funcionalidade
fix: correção de bug
docs: documentação apenas
chore: manutenção, infraestrutura
test: testes apenas
refactor: refatoração sem mudança de comportamento
```

### Regras de preservação

- `src/fortran_original/` - **nunca alterar**. Referência imutável.
- `src/mathematica_original/` - **nunca alterar**. Referência imutável.
- `Thesis/` - **nunca alterar**. Material original da tese.

---

## Checklist de início de sessão

- [ ] `git pull origin main` - verificar sincronização
- [ ] `python -m pytest -q` - confirmar que todos os testes passam
- [ ] Identificar a fase atual neste documento
- [ ] Solicitar permissões necessárias (skill `dev`) antes de começar
- [ ] Usar modelo indicado na fase

## Checklist de fim de sessão

- [ ] Todos os testes passam (`pytest -q`)
- [ ] `ruff check src/ tests/` sem erros
- [ ] Branch pushed, PR criado ou merged
- [ ] `.claude/PLANEJAMENTO.md` atualizado com log da sessão
- [ ] Este documento atualizado se houver mudança de fase
