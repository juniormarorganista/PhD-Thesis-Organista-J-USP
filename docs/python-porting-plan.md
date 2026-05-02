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
