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
