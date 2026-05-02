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
