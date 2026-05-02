# Port Python do caso manufaturado 12

Este documento registra a primeira transcrição efetiva do código Fortran para a
biblioteca Python `mms_viscoelastic`.

## Arquivo de origem

O módulo Python foi baseado na subrotina Fortran:

```text
src/fortran_original/4Order_LPOG_Euler/case_12.f
```

A rotina transcrita nesta fase é:

```fortran
exact_solutions_case_12(t)
```

## Arquivo Python criado

```text
src/mms_viscoelastic/manufactured_solutions/case_12.py
```

## Campos implementados

A função `evaluate_case_12_exact` retorna um objeto `Field2D` com:

- `u`: velocidade na direção x;
- `v`: velocidade na direção y;
- `wz`: vorticidade;
- `psi`: função de corrente;
- `txx`: componente normal xx do tensor de tensões;
- `txy`: componente cisalhante xy do tensor de tensões;
- `tyy`: componente normal yy do tensor de tensões.

## Convenção de indexação

No Fortran, os campos são acessados como:

```fortran
campo(i,j)
```

No Python/NumPy, a convenção adotada é:

```python
campo[j, i]
```

ou seja, os arrays têm forma:

```text
(ny, nx)
```

Isso é essencial para evitar troca de eixos nos gráficos e nos datasets.

## Preparação para Machine Learning

O exemplo:

```text
examples/02_case12_manufactured_solution.py
```

permite gerar uma tabela com colunas:

```text
x, y, u, v, psi, wz, txx, txy, tyy
```

Essa tabela é a primeira base limpa para futuros estudos com Machine Learning,
como regressão de campos, operadores neurais, surrogates ou modelos híbridos.

## Limitação desta fase

A subrotina `source_term_case_12(t)` ainda não foi transcrita. Ela contém uma
expressão algébrica muito extensa e deve ser portada em fase própria, com testes
separados e comparação cuidadosa com os dados gerados pelo Fortran.

## Testes criados

```text
tests/test_case_12_manufactured.py
```

Os testes verificam:

1. formato dos campos;
2. valores conhecidos em pontos específicos;
3. conversão do campo para `pandas.DataFrame`;
4. compatibilidade inicial com pipeline de Machine Learning.
