# Mapa inicial do repositório

Este repositório reúne os arquivos originais da tese, códigos Fortran, notebooks
Mathematica, artigo científico e a nova fundação Python.

## Regra de preservação

Os códigos Fortran e os arquivos LaTeX originais não devem ser alterados nesta
fase. Eles foram apenas reorganizados com `git mv`, preservando o histórico do Git.

## Estrutura após a reorganização

```text
src/
  fortran_original/        # código Fortran original da tese
  mathematica_original/    # notebooks Mathematica originais
  mms_viscoelastic/        # nova biblioteca Python em construção

thesis/
  original/                # tese LaTeX original

papers/
  original/                # artigo e figuras originais

docs/                      # documentação técnica do repositório
examples/                  # exemplos executáveis da biblioteca Python
tests/                     # testes automatizados da biblioteca Python
```

## Objetivo da fundação Python

A pasta `src/mms_viscoelastic` foi criada para evoluir para uma biblioteca Python.
A transcrição Fortran -> Python deve ser incremental, testável e orientada por
comparação contra a implementação original.
