# Plano de transcrição Fortran -> Python

## Princípio

A versão Python deve nascer como biblioteca científica, não como um conjunto solto
de scripts. O código Fortran original permanece como referência.

## Ordem sugerida de transcrição

1. Parâmetros físicos e numéricos.
2. Malhas e convenções de indexação.
3. Campos bidimensionais e saída de dados.
4. Soluções manufaturadas dos arquivos `case_01.f` até `case_14.f`.
5. Operadores diferenciais de segunda ordem.
6. Operadores compactos de alta ordem.
7. Normas de erro e estimativa de ordem de convergência.
8. Solver de Poisson.
9. Integração temporal Euler/RK.
10. Modelos UCM, Oldroyd-B, Giesekus e LPTT.
11. Comparação automatizada Fortran vs Python.
12. Geração de datasets para Machine Learning.

## Preparação para Machine Learning

A função `field_to_feature_table` já converte malhas e campos em `pandas.DataFrame`,
criando uma base natural para regressão, classificação, surrogates, PINNs/FNOs ou
outros métodos de Scientific Machine Learning.

## Cuidados críticos

- Manter a convenção NumPy: campos com shape `(ny, nx)`.
- Documentar qualquer diferença em relação aos índices Fortran.
- Testar cada módulo isoladamente.
- Não validar o Python apenas visualmente; comparar normas e ordens de convergência.
