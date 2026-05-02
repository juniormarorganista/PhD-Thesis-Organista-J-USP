from pathlib import Path
from datetime import date

ROOT = Path.cwd()

def write(path: str, content: str, overwrite: bool = True) -> None:
    p = ROOT / path
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.exists() and not overwrite:
        print(f"[SKIP] {path} ja existe")
        return
    if p.exists() and path == "README.md":
        backup = ROOT / "README.before-phase2.md"
        if not backup.exists():
            backup.write_text(p.read_text(encoding="utf-8", errors="replace"), encoding="utf-8", newline="\n")
            print("[BACKUP] README.md -> README.before-phase2.md")
    p.write_text(content.strip() + "\n", encoding="utf-8", newline="\n")
    print(f"[OK] {path}")

readme = r'''
# Verification of High-Order Viscoelastic Flow Solvers Using Manufactured Solutions

Repositório de organização, documentação e modernização dos códigos associados à tese de doutorado de **Juniormar Organista**.

O projeto reúne códigos originais em Fortran, notebooks Mathematica, arquivos LaTeX da tese, apresentação da defesa, material de artigo científico e uma fundação em Python para futura transcrição, validação e uso em fluxos de *Scientific Machine Learning*.

> Este repositório está em fase de reorganização. Os arquivos originais em Fortran, LaTeX e Mathematica são preservados como fontes históricas e científicas. A nova camada Python está sendo construída separadamente para validação, reprodutibilidade e futuras aplicações com aprendizado de máquina.

## Palavras-chave

`computational-fluid-dynamics`, `viscoelastic-flows`, `manufactured-solutions`, `MMS`, `Oldroyd-B`, `Giesekus`, `LPTT`, `UCM`, `log-conformation`, `finite-differences`, `compact-schemes`, `Fortran`, `Python`, `scientific-machine-learning`, `machine-learning`, `numerical-verification`, `CFD`, `non-newtonian-fluids`

## Sobre o trabalho

A tese trata da **verificação numérica de um solucionador incompressível para escoamentos viscoelásticos**, utilizando o **Método das Soluções Manufaturadas** (*Method of Manufactured Solutions*, MMS). O foco está na avaliação da consistência computacional de esquemas numéricos de alta ordem aplicados a modelos constitutivos não newtonianos.

O problema científico central pode ser resumido assim:

> Como verificar, de maneira sistemática, se a implementação numérica de um código para escoamentos viscoelásticos reproduz corretamente a ordem de convergência esperada para diferentes modelos constitutivos e diferentes variáveis do escoamento?

A resposta construída na tese passa por quatro elementos:

1. formulação matemática dos escoamentos incompressíveis viscoelásticos;
2. construção de soluções manufaturadas e respectivos termos fonte;
3. implementação de esquemas numéricos de alta ordem;
4. análise de erro e estimativa da ordem de convergência.

## Modelos e variáveis de interesse

O material original contempla formulações e experimentos associados a modelos como:

- fluido Newtoniano como caso de referência;
- UCM (*Upper-Convected Maxwell*);
- Oldroyd-B;
- Giesekus;
- LPTT (*Linear Phan-Thien--Tanner*);
- formulações com transformação log-conformação.

As principais variáveis envolvidas incluem:

- velocidade horizontal `u`;
- velocidade vertical `v`;
- função corrente `psi`;
- vorticidade `wz`;
- tensões poliméricas `Txx`, `Txy`, `Tyy`;
- variáveis associadas à formulação log-conformação.

## Estrutura atual do repositório

```text
.
├── docs/                         # documentação técnica e planos de reorganização
├── examples/                     # exemplos mínimos da biblioteca Python
├── papers/                       # material original de artigo científico
├── src/
│   ├── fortran_original/          # código Fortran original preservado
│   ├── mathematica_original/      # notebooks Mathematica originais preservados
│   └── mms_viscoelastic/          # fundação da futura biblioteca Python
├── tests/                        # testes da camada Python
├── thesis/                       # tese, apresentação e arquivos LaTeX originais
├── pyproject.toml                # configuração da biblioteca Python
├── requirements.txt              # dependências mínimas
└── requirements-dev.txt          # dependências de desenvolvimento
```

## Código Fortran original

O código Fortran original foi preservado em:

```text
src/fortran_original/4Order_LPOG_Euler/
```

Arquivos centrais:

| Arquivo | Função esperada |
|---|---|
| `main.f` | fluxo principal da simulação |
| `case_01.f`--`case_14.f` | casos manufaturados e termos fonte |
| `derivs.f` | operadores diferenciais e derivadas |
| `poisson.f` | solução da equação de Poisson |
| `solvers.f` | rotinas auxiliares de solução numérica |
| `euler_rk.f` | integração temporal |
| `conform.f` | formulação associada à log-conformação |
| `normerror.f` | cálculo de normas de erro |
| `write_solutions.f` | escrita das soluções/resultados |

A reorganização do repositório não altera esses arquivos. Eles são a referência para a transcrição progressiva em Python.

## Fundação da biblioteca Python

A nova camada Python está em:

```text
src/mms_viscoelastic/
```

Objetivos da camada Python:

1. transcrever progressivamente a lógica dos códigos Fortran;
2. permitir testes unitários e validação cruzada Fortran/Python;
3. gerar datasets estruturados para aprendizado de máquina;
4. facilitar experimentos futuros com regressão, classificação, operadores neurais e modelos híbridos física-dados;
5. transformar o projeto em uma biblioteca Python científica reutilizável.

Instalação em modo desenvolvimento:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip setuptools wheel
python -m pip install -r requirements-dev.txt
python -m pip install -e .
```

Execução dos testes:

```powershell
python -m pytest -q
```

Execução do exemplo mínimo:

```powershell
python examples\01_smoke_python_package.py
```

## Plano de transcrição Fortran → Python

A transcrição será incremental. A ordem recomendada é:

1. mapear cada arquivo Fortran;
2. transcrever estruturas de parâmetros e malhas;
3. transcrever soluções manufaturadas;
4. transcrever operadores diferenciais;
5. transcrever cálculo de normas de erro;
6. implementar validação de ordem de convergência;
7. transcrever rotinas de Poisson e integração temporal;
8. criar comparação Fortran vs Python;
9. gerar datasets para ML;
10. documentar exemplos e notebooks.

A regra é simples: **cada módulo Python novo deve nascer com teste**. Código científico sem teste é só superstição com `for`.

## Resultados e evidências numéricas

Os resultados originais da tese e do artigo incluem famílias de figuras e tabelas relacionadas a:

- mapas de erro;
- superfícies de erro;
- cortes em linhas específicas do domínio;
- comparação entre solução exata/manufaturada e solução numérica;
- normas de erro;
- estimativas de ordem de convergência;
- variações de parâmetros como `Re`, `Wi`, `beta`, `alphaG`, `epsilon`, `xi` e `dt`.

As figuras originais permanecem preservadas nos diretórios da tese e do artigo. Uma etapa posterior deverá selecionar figuras representativas e convertê-las para `PNG` em:

```text
assets/preview/
```

Essas figuras serão usadas no README e na documentação, evitando depender diretamente de arquivos EPS pesados.

## Como citar

Use o arquivo `CITATION.cff` deste repositório quando disponível. Em trabalhos acadêmicos, cite também a tese de doutorado original e os artigos derivados quando aplicável.

## Status do projeto

- [x] Clone corrigido no Windows com caminhos longos.
- [x] Reorganização inicial de diretórios.
- [x] Fundação Python criada.
- [x] Testes mínimos funcionando.
- [ ] README com figuras representativas.
- [ ] Mapa completo do Fortran.
- [ ] Transcrição do primeiro caso manufaturado.
- [ ] Validação Fortran/Python.
- [ ] Geração de datasets para Machine Learning.
- [ ] Notebooks explicativos.
- [ ] GitHub Actions.

## Autor

**Juniormar Organista**  
Doutorado em Matemática Aplicada e Computacional  
Ênfase: Mecânica dos Fluidos Computacional, Métodos Numéricos, Escoamentos Viscoelásticos e Verificação Numérica.
'''

project_map = r'''
# Mapa do Projeto

Este documento descreve a organização atual do repositório após a reorganização inicial.

## Princípio de organização

A estrutura separa quatro camadas:

1. **fontes originais preservadas**: Fortran, Mathematica, LaTeX da tese e material do artigo;
2. **documentação técnica**: mapas, planos e instruções de uso;
3. **camada Python nova**: biblioteca em desenvolvimento;
4. **testes e exemplos**: validação da nova camada Python.

## Diretórios principais

```text
src/fortran_original/
```

Contém o código Fortran original. Não deve ser alterado durante a reorganização. Qualquer melhoria deve ser feita em camada separada ou em branch específica de manutenção.

```text
src/mathematica_original/
```

Contém os notebooks Mathematica originais usados na geração de soluções manufaturadas e/ou termos simbólicos.

```text
src/mms_viscoelastic/
```

Contém a nova biblioteca Python. Este diretório será expandido para comportar:

- parâmetros físicos;
- malhas;
- campos numéricos;
- soluções manufaturadas;
- operadores diferenciais;
- solucionadores;
- rotinas de erro;
- exportação de dados;
- integração com Machine Learning.

```text
thesis/
```

Contém a tese e sua apresentação original. Os arquivos LaTeX originais devem ser preservados.

```text
papers/
```

Contém o material original de artigo científico e figuras associadas.

```text
docs/
```

Contém documentação técnica e planos de desenvolvimento.

```text
tests/
```

Contém testes automatizados da camada Python.

```text
examples/
```

Contém exemplos mínimos executáveis.

## Convenção de evolução

Novas mudanças devem ser feitas em branches temáticas:

```text
chore/...   reorganização, infraestrutura e manutenção
docs/...    documentação
feat/...    novas funcionalidades
fix/...     correções
tests/...   testes
slides/...  apresentação da tese
```

## Regra de preservação

Arquivos originais de tese, Fortran e Mathematica devem ser considerados referência. A transcrição em Python deve ser validada contra esses arquivos, não substituí-los.
'''

fortran_map = r'''
# Mapa Inicial do Código Fortran

Diretório original reorganizado:

```text
src/fortran_original/4Order_LPOG_Euler/
```

Este mapa orienta a futura transcrição para Python.

## Arquivos principais

| Arquivo | Papel no código original | Destino Python planejado |
|---|---|---|
| `main.f` | programa principal e orquestração da simulação | `mms_viscoelastic.runners` |
| `comm.var`, `comm.vare` | variáveis comuns/parâmetros globais | `mms_viscoelastic.parameters` |
| `par.nn`, `par.nn.template` | parâmetros de execução | `mms_viscoelastic.config` |
| `case_01.f`--`case_14.f` | casos MMS e termos fonte | `mms_viscoelastic.manufactured.cases` |
| `derivs.f` | derivadas e operadores diferenciais | `mms_viscoelastic.numerics.derivatives` |
| `drvs.f` | possível variante/auxiliar de derivadas | avaliar antes da transcrição |
| `poisson.f` | solução da equação de Poisson | `mms_viscoelastic.solvers.poisson` |
| `solvers.f` | rotinas auxiliares de solução | `mms_viscoelastic.solvers.linear` |
| `euler_rk.f` | integração temporal Euler/RK | `mms_viscoelastic.solvers.time` |
| `conform.f` | transformação log-conformação | `mms_viscoelastic.models.log_conformation` |
| `condictions.f` | condições iniciais/de contorno | `mms_viscoelastic.boundary` |
| `boundarie_vorticity_cavity_flow.f` | condições de vorticidade no problema da cavidade | `mms_viscoelastic.boundary.vorticity` |
| `normerror.f` | cálculo de normas de erro | `mms_viscoelastic.numerics.errors` |
| `write_solutions.f` | escrita de resultados | `mms_viscoelastic.io` |
| `loop.f` | laços principais da simulação | `mms_viscoelastic.runners` |
| `formatted.f` | formatação/saída auxiliar | `mms_viscoelastic.io` |
| `makefile` | compilação Fortran | manter em `src/fortran_original` |

## Ordem recomendada de transcrição

1. `comm.var`, `comm.vare`, `par.nn` → parâmetros Python;
2. `case_12.f` → primeiro caso MMS de referência;
3. `normerror.f` → erros e ordem de convergência;
4. `derivs.f` → derivadas de segunda ordem inicialmente;
5. `poisson.f` → solver de Poisson simplificado;
6. `euler_rk.f` → integração temporal;
7. `conform.f` → log-conformação;
8. demais casos `case_01.f`--`case_14.f`.

## Critério de aceite por módulo

Um módulo transcrito só deve ser considerado aceito quando:

- possuir teste unitário;
- possuir exemplo mínimo;
- tiver documentação curta;
- reproduzir uma saída conhecida ou uma propriedade matemática esperada;
- não alterar os arquivos Fortran originais.
'''

python_plan = r'''
# Plano de Transcrição Fortran para Python

## Objetivo

Construir uma biblioteca Python chamada `mms_viscoelastic`, baseada nos códigos Fortran originais da tese, preparada para:

- verificação numérica;
- geração de dados;
- análise de erro;
- comparação Fortran/Python;
- integração futura com Machine Learning e Scientific Machine Learning.

## Filosofia da transcrição

A transcrição não deve ser uma cópia cega linha a linha. O Fortran original é procedural e orientado à execução numérica direta. A camada Python deve ser modular, testável e extensível.

Portanto, a tradução será conceitual:

```text
Fortran original -> interpretação matemática -> API Python limpa -> testes -> validação
```

## Módulos planejados

```text
mms_viscoelastic/
├── parameters.py
├── grids.py
├── fields.py
├── manufactured/
│   ├── __init__.py
│   ├── base.py
│   ├── case_12.py
│   └── cases.py
├── numerics/
│   ├── derivatives.py
│   ├── errors.py
│   ├── convergence.py
│   └── filters.py
├── models/
│   ├── newtonian.py
│   ├── ucm.py
│   ├── oldroyd_b.py
│   ├── giesekus.py
│   ├── lptt.py
│   └── log_conformation.py
├── solvers/
│   ├── poisson.py
│   ├── linear.py
│   └── time.py
├── ml/
│   ├── datasets.py
│   ├── features.py
│   └── splits.py
└── io/
    ├── readers.py
    └── writers.py
```

## Fases

### Fase 1 — Fundação

Já iniciada:

- pacote instalável;
- objetos de malha;
- objetos de campo;
- erros básicos;
- conversão campo → tabela;
- testes mínimos.

### Fase 2 — Documentação e mapa do Fortran

Atual:

- README profissional;
- mapa do projeto;
- mapa dos arquivos Fortran;
- plano de transcrição;
- instruções GitHub.

### Fase 3 — Primeiro caso MMS

Implementar `case_12.f` como primeiro caso Python:

```text
src/mms_viscoelastic/manufactured/case_12.py
```

Deve conter funções para:

- `u_exact(x, y, t, params)`;
- `v_exact(x, y, t, params)`;
- `psi_exact(x, y, t, params)`;
- `wz_exact(x, y, t, params)`;
- `stress_exact(x, y, t, params)`;
- termos fonte quando transcritos.

### Fase 4 — Derivadas e verificação

Implementar derivadas com diferenças finitas e testar usando funções analíticas suaves.

### Fase 5 — Dados para Machine Learning

Criar um pipeline que gere tabelas com colunas como:

```text
x, y, t, Re, Wi, beta, alphaG, epsilon, u, v, psi, wz, Txx, Txy, Tyy
```

Esses dados poderão alimentar:

- regressão supervisionada;
- modelos de erro;
- surrogates;
- operadores neurais;
- classificadores de regimes;
- modelos híbridos físico-estatísticos.

## Regra de qualidade

Toda transcrição deve seguir o ciclo:

```text
ler Fortran -> documentar função -> implementar Python -> testar -> comparar -> versionar
```

Nada de transcrição heroica de madrugada sem teste. Isso é como resolver Navier-Stokes no guardanapo: romântico, mas perigoso.
'''

github_about = r'''
# GitHub About

## Descrição curta sugerida

```text
Verification and Python modernization of high-order viscoelastic flow solvers using the Method of Manufactured Solutions.
```

## Descrição curta em português

```text
Organização, verificação e modernização em Python de códigos de escoamentos viscoelásticos usando o Método das Soluções Manufaturadas.
```

## Website

Use o link do artigo, da tese, do currículo, do laboratório ou deixe em branco por enquanto.

## Topics sugeridos

```text
computational-fluid-dynamics
viscoelastic-flows
manufactured-solutions
method-of-manufactured-solutions
mms
oldroyd-b
giesekus
lptt
ucm
non-newtonian-fluids
finite-differences
compact-schemes
log-conformation
fortran
python
scientific-computing
scientific-machine-learning
machine-learning
numerical-methods
numerical-verification
cfd
applied-mathematics
```

## Instruções para configurar no GitHub

1. Abra o repositório no GitHub.
2. Clique na engrenagem ao lado de **About**.
3. Cole a descrição curta.
4. Adicione os topics sugeridos.
5. Marque, se quiser:
   - Releases;
   - Packages, se futuramente publicar a biblioteca;
   - Environments, se usar GitHub Actions.

## Sugestão de nome futuro do pacote

Nome atual no `pyproject.toml`:

```text
mms-viscoelastic
```

Nome de importação:

```python
import mms_viscoelastic
```
'''

windows_setup = r'''
# Configuração no Windows

Este projeto contém muitos arquivos LaTeX, imagens e caminhos longos. No Windows, isso pode causar falhas no clone ou na extração.

## 1. Ativar caminhos longos no Git

```powershell
git config --global core.longpaths true
```

## 2. Ativar caminhos longos no Windows

Abra o PowerShell como administrador:

```powershell
New-ItemProperty -Path "HKLM:\SYSTEM\CurrentControlSet\Control\FileSystem" `
  -Name "LongPathsEnabled" `
  -Value 1 `
  -PropertyType DWORD `
  -Force
```

Reinicie o computador.

## 3. Clonar em caminho curto

Evite OneDrive, Área de Trabalho e pastas com nomes longos.

Use:

```powershell
cd C:\
mkdir dev
cd C:\dev
git clone URL_DO_REPOSITORIO tese
cd tese
```

## 4. Instalar camada Python

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip setuptools wheel
python -m pip install -r requirements-dev.txt
python -m pip install -e .
python -m pytest -q
```

Se a ativação do ambiente virtual for bloqueada:

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass -Force
.\.venv\Scripts\Activate.ps1
```
'''

presentation_plan = r'''
# Plano de Melhoria da Apresentação da Tese

## Diagnóstico

A apresentação original contém muito material técnico e segue de perto a estrutura da tese. Para defesa ou apresentação pública, a narrativa deve ser mais seletiva.

A apresentação melhorada deve enfatizar:

1. problema científico;
2. dificuldade numérica;
3. estratégia MMS;
4. formulação mínima necessária;
5. método numérico;
6. resultados centrais;
7. contribuição científica.

## Estrutura recomendada

```text
1. Título
2. Mensagem central da tese
3. Problema físico e numérico
4. Por que escoamentos viscoelásticos são difíceis
5. Objetivos
6. Formulação incompressível
7. Modelos constitutivos
8. Método numérico
9. Método das Soluções Manufaturadas
10. Pipeline de verificação
11. Resultados principais
12. Comparativo de modelos
13. Contribuições
14. Limitações
15. Trabalhos futuros
16. Conclusão
17. Agradecimentos
```

## Regra de ouro

Cada slide deve responder a uma pergunta. Se um slide exige cinco minutos de explicação antes de fazer sentido, ele provavelmente precisa virar três slides ou ir para backup.

## Próximas ações

- criar `presentation/improved/`;
- criar versão limpa da defesa;
- mover slides longos para backup;
- converter figuras importantes para PNG;
- criar uma narrativa visual com menos equações e mais fluxos.
'''

contributing = r'''
# Contribuindo

Este repositório está sendo reorganizado para preservar o material original da tese e construir uma camada Python moderna.

## Regras principais

1. Não alterar diretamente os arquivos originais em `src/fortran_original/`.
2. Não alterar diretamente os arquivos originais da tese em `thesis/` sem branch específica.
3. Toda nova funcionalidade Python deve ter teste.
4. Toda transcrição Fortran → Python deve indicar o arquivo Fortran de origem.
5. Commits devem ser pequenos e descritivos.

## Convenção de branches

```text
chore/...   organização e manutenção
docs/...    documentação
feat/...    novas funcionalidades
fix/...     correções
tests/...   testes
slides/...  apresentação da tese
```

## Convenção de commits

Exemplos:

```text
docs: add Fortran code map
feat: implement manufactured case 12 in Python
tests: add convergence tests for finite differences
chore: reorganize original thesis sources
slides: restructure thesis defense presentation
```
'''

changelog = r'''
# Changelog

Todas as mudanças relevantes deste projeto serão documentadas aqui.

## [0.1.0] - Em desenvolvimento

### Adicionado

- Reorganização inicial do repositório.
- Preservação dos códigos Fortran originais em `src/fortran_original/`.
- Preservação dos notebooks Mathematica originais em `src/mathematica_original/`.
- Fundação da biblioteca Python `mms_viscoelastic`.
- Testes mínimos da camada Python.
- Documentação inicial do projeto.

### Planejado

- Mapa completo do código Fortran.
- Transcrição do primeiro caso MMS.
- Validação Fortran/Python.
- Geração de datasets para Machine Learning.
- Melhoria da apresentação da tese.
'''

citation = r'''
cff-version: 1.2.0
message: "If you use this repository, please cite the thesis and this software repository."
title: "Verification of High-Order Viscoelastic Flow Solvers Using Manufactured Solutions"
authors:
  - family-names: "Organista"
    given-names: "Juniormar"
version: "0.1.0"
date-released: "2026-05-02"
keywords:
  - "computational fluid dynamics"
  - "viscoelastic flows"
  - "manufactured solutions"
  - "Fortran"
  - "Python"
  - "scientific machine learning"
license: "MIT"
repository-code: "https://github.com/SEU_USUARIO/SEU_REPOSITORIO"
'''

gitignore = r'''
# Python
__pycache__/
*.py[cod]
*.pyo
*.pyd
*.egg-info/
.eggs/
build/
dist/
.venv/
venv/
.env
.pytest_cache/
.mypy_cache/
.ruff_cache/
.coverage
htmlcov/

# Jupyter
.ipynb_checkpoints/

# LaTeX build artifacts
*.aux
*.bbl
*.bcf
*.blg
*.fdb_latexmk
*.fls
*.idx
*.ilg
*.ind
*.lof
*.log
*.lot
*.nav
*.out
*.run.xml
*.snm
*.synctex.gz
*.toc
*.vrb
*.xdv

# Fortran build artifacts
*.o
*.obj
*.mod
*.smod
*.exe
*.out
*.a
*.so
*.dll

# OS/editor
.DS_Store
Thumbs.db
.vscode/
.idea/

# Project local outputs
outputs/
results/
tmp/
temp/
'''

gitattributes = r'''
# Normalize text files
* text=auto

# Keep scripts with LF where useful
*.sh text eol=lf
*.py text eol=lf
*.md text eol=lf
*.tex text eol=lf
*.bib text eol=lf
*.f text eol=lf
*.f90 text eol=lf

# Binary/scientific assets
*.pdf binary
*.png binary
*.jpg binary
*.jpeg binary
*.eps binary
*.nb binary

# Suggested Git LFS patterns for future large assets.
# Existing tracked files are not migrated automatically by these lines.
*.eps filter=lfs diff=lfs merge=lfs -text
*.pdf filter=lfs diff=lfs merge=lfs -text
*.png filter=lfs diff=lfs merge=lfs -text
*.jpg filter=lfs diff=lfs merge=lfs -text
*.jpeg filter=lfs diff=lfs merge=lfs -text
*.tif filter=lfs diff=lfs merge=lfs -text
*.tiff filter=lfs diff=lfs merge=lfs -text
'''

workflow_tests = r'''
name: Python tests

on:
  push:
    branches: [main]
  pull_request:
    branches: [main]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ["3.10", "3.11", "3.12"]

    steps:
      - name: Checkout repository
        uses: actions/checkout@v4

      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}

      - name: Install dependencies
        run: |
          python -m pip install --upgrade pip setuptools wheel
          python -m pip install -r requirements-dev.txt
          python -m pip install -e .

      - name: Run tests
        run: |
          python -m pytest -q
'''

write("README.md", readme)
write("docs/project-map.md", project_map)
write("docs/fortran-code-map.md", fortran_map)
write("docs/python-porting-plan.md", python_plan)
write("docs/github-about.md", github_about)
write("docs/windows-setup.md", windows_setup)
write("docs/presentation-improvement-plan.md", presentation_plan)
write("CONTRIBUTING.md", contributing)
write("CHANGELOG.md", changelog)
write("CITATION.cff", citation)
write(".gitignore", gitignore)
write(".gitattributes", gitattributes)
write(".github/workflows/python-tests.yml", workflow_tests)

print("\nFase 2 concluida: documentacao e arquivos GitHub criados/atualizados.")
print("Revise os arquivos antes do commit.")
