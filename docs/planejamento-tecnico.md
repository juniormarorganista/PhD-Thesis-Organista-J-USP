# Planejamento Técnico - mms_viscoelastic

> Documento de referência estratégica para o desenvolvimento da biblioteca.
> Atualizado em: 2026-06-03.

---

## Estado atual do projeto

| Item | Estado |
|---|---|
| Testes | ✅ 36/36 passam |
| Local ↔ GitHub | ✅ Sincronizado |
| PRs abertos | - nenhum |
| Fundação Python | ✅ `grids`, `fields`, `parameters`, `numerics`, `ml`, `io`, `visualization` |
| Case 12 - solução manufaturada | ✅ portada e testada |
| Case 12 - termos fonte (tfwz, tfTxx, tfTxy, tfTyy) | ✅ portados e validados vs Fortran |
| Módulo de convergência | ⏳ `numerics/errors.py` existe, mas pipeline de convergência está ausente |
| Modelos constitutivos | ⏳ enum `ConstitutiveModel` existe; física ainda não implementada |
| Relatório científico do Case 12 | ⏳ pendente |
| Notebook científico | ⏳ pendente |
| CI/CD (GitHub Actions) | ✅ pytest em 3.10/3.11/3.12 |
| Versão do pacote | 0.1.0 |

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

### ▸ FASE I - Organização e baseline de qualidade
**Modelo: haiku** | **Skills: `dev` + `git-github`**

**Objetivo:** garantir que o projeto está limpo, lint-free e com CI sólido antes de qualquer expansão.

| Tarefa | Arquivo / Comando | Branch |
|---|---|---|
| Rodar ruff e corrigir avisos | `ruff check src/ tests/` | `chore/lint-and-type-check` |
| Rodar mypy e corrigir type errors | `mypy src/` | `chore/lint-and-type-check` |
| Verificar se GitHub Actions roda corretamente em push | `.github/workflows/tests.yml` | `chore/lint-and-type-check` |
| Adicionar ruff ao CI | `.github/workflows/tests.yml` | `chore/lint-and-type-check` |
| PR + merge via GitHub | - | - |

---

### ▸ FASE II - Verificação de convergência numérica (núcleo do MMS)
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

---

### ▸ FASE III - Relatório científico e notebook do Case 12
**Modelo: sonnet** | **Skills: `escritor` + `doutor-em-matematica` + `git-github`**

**Objetivo:** consolidar toda a evidência científica do Case 12 em um documento e notebook reprodutíveis.

| Tarefa | Arquivo | Branch |
|---|---|---|
| Criar relatório completo de validação | `docs/case-12-complete-validation-report.md` | `docs/case12-scientific-report` |
| Criar notebook científico reprodutível | `examples/notebooks/Case12_MMS_Pipeline.ipynb` | idem |
| Atualizar README com seção "Convergence results" | `README.md` | idem |
| PR + merge via GitHub | - | - |

> **Regra do `escritor`:** redigir na voz do pesquisador. Evitar fraseado genérico de IA. Citar a tese e os papers originais quando mencionar equações.

---

### ▸ FASE IV - Expansão dos datasets ML
**Modelo: haiku** | **Skills: `dev` + `git-github`**

**Objetivo:** ampliar o pipeline de datasets para cobrir múltiplos regimes de parâmetros, habilitando treinamento de modelos surrogates.

| Tarefa | Arquivo | Branch |
|---|---|---|
| Expandir `ml/datasets.py` para varredura paramétrica (Re, Wi, beta, t) | `src/mms_viscoelastic/ml/datasets.py` | `feat/ml-parametric-dataset` |
| Criar `ml/splits.py` com train/val/test stratificado | `src/mms_viscoelastic/ml/splits.py` | idem |
| Criar `examples/10_case12_ml_dataset_parametric.py` | `examples/10_case12_ml_dataset_parametric.py` | idem |
| Testes para datasets e splits | `tests/test_ml_parametric.py`, `tests/test_ml_splits.py` | idem |
| PR + merge via GitHub | - | - |

---

### ▸ FASE V - Modelos constitutivos (UCM e Oldroyd-B)
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

---

### ▸ FASE VI - Versioning e CI
**Modelo: haiku** | **Skills: `dev` + `git-github`**

**Objetivo:** versionar o primeiro ciclo completo de desenvolvimento e atualizar o CI.

| Tarefa | Arquivo | Branch |
|---|---|---|
| Atualizar versão para `0.2.0` em `_version.py` e `pyproject.toml` | ambos | `chore/release-0.2.0` |
| Adicionar CHANGELOG entry para 0.2.0 | `CHANGELOG.md` | idem |
| Adicionar `ruff` e `mypy` ao GitHub Actions | `.github/workflows/tests.yml` | idem |
| Criar tag `v0.2.0` no GitHub | - | - |
| PR + merge via GitHub | - | - |

---

### ▸ FASE VII - Comparação solver-level Fortran × Python
**Modelo: sonnet** | **Skills: `doutor-em-matematica` + `dev` + `escritor`**

**Objetivo:** demonstrar equivalência numérica entre o solver Fortran original e a implementação Python para o Case 12.

| Tarefa | Arquivo | Branch |
|---|---|---|
| Criar pipeline de comparação campo a campo | `src/mms_viscoelastic/validation/solver_comparison.py` | `feat/solver-comparison` |
| Executar Fortran reference e capturar outputs | `validation/fortran_reference/case12_solver/` | idem |
| Comparar Python vs Fortran nos campos u, v, wz, psi, T | `examples/11_case12_solver_comparison.py` | idem |
| Documentar divergências e justificar tolerâncias | `docs/case-12-solver-comparison.md` | idem |
| PR + merge via GitHub | - | - |

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
