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
