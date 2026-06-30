# Changelog

Todas as mudanças relevantes deste projeto serão documentadas aqui.

O formato segue [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

As versões `0.1.0` a `0.6.7.5` cobrem a fase inicial de reorganização do
repositório e da documentação, antes da existência de um pacote Python
formal. A partir da criação do `pyproject.toml`, o versionamento semântico
do pacote `mms_viscoelastic` recomeçou em `0.1.0`/`0.2.0`; é por isso que a
numeração não é monotonicamente crescente ao longo do histórico abaixo.

---

## [Unreleased]

## [1.0.0] - 2026-06-29

Primeira versão estável do pacote `mms_viscoelastic`. Consolida todas as 17
etapas de portagem e verificação: os 14 casos MMS completos, o esquema compacto
de 4ª ordem + filtro de Lele, o solver multigrid V-cycle FAS (`poi_psi4a`) e o
caminho log-conformation (variante Euler). Toda a suite de testes (508 casos)
passa com limpeza em ruff e mypy. Esta release marca o repositório como pronto
para publicação pública.

## [0.3.1] - 2026-06-29

Consolida a Etapa 17: marcha temporal de 4ª ordem real com multigrid V-cycle FAS
(`poi_psi4a`), filtro de Lele e derivadas compactas de 4ª ordem como pipeline
completo. Inclui achado de código morto (`band5_poi`/`banbk5_poi` nunca
chamados), cinco referências Fortran reescritas e nova documentação técnica.

### Adicionado
- Etapa 17 (marcha temporal de 4ª ordem real): achado de que
  `band5_poi`/`banbk5_poi` (`solvers.f`) é código morto — nunca chamado em
  nenhum lugar do Fortran original; o solver de Poisson real é `poi_psi4a`
  (multigrid V-cycle FAS sobre stencil compacto de 9 pontos tipo
  Mehrstellen). Portado em `numerics/poisson_multigrid.py`
  (`solve_poisson_multigrid_4th_order`), validado bit-a-bit contra
  referência Fortran isolada (`validation/fortran_reference/case12_poisson_multigrid/`)
  e à ordem 4 de convergência. `numerics/time_stepping.py`,
  `numerics/time_stepping_log.py` e `numerics/time_marching.py`
  reapontados para as derivadas compactas de 4ª ordem (Etapa 13) + filtro
  de Lele (aplicado uma vez por passo, entre a imposição de contorno e o
  resolve de Poisson, replicando `subroutine loop` do `main.f`) + o
  multigrid real, substituindo o placeholder de 2ª ordem usado desde as
  Etapas 8-11. As cinco referências Fortran isoladas associadas
  (`case12_euler_step`, `case12_log_euler_step`, `case12_time_step`,
  `case12_multistep`, `case12_full_march`) foram reescritas e os CSVs
  regenerados; toda comparação continua batendo a precisão de máquina
  (`<1e-12`). Achado adicional: a instabilidade documentada em
  `docs/case-12-full-march-reference.md` (Etapa 12, divergência do Euler
  explícito em malha 9x9) persiste mesmo com o esquema real e dissipativo
  (filtro de Lele), só que mais cedo — refutando a hipótese original de que
  o esquema dissipativo não compartilharia essa instabilidade. Documentado
  em `docs/case-12-4th-order-time-marching.md`.

## [0.3.0] - 2026-06-29

Consolida as Etapas 13-16: esquema compacto de 4ª ordem + filtro de Lele,
os 14 casos MMS portados, o caminho log-conformation (Euler) e a
finalização de documentação/README. Ver entradas individuais abaixo.

### Adicionado
- Etapa 16 (finalização): `README.md` atualizado (status dos 14 casos MMS,
  esquema de 4ª ordem, log-conformation; nota explícita classificando
  `ml/` como extensão exploratória fora do núcleo da tese);
  `docs/fortran-code-map.md`, `docs/python-porting-plan.md`,
  `docs/project-map.md`, `docs/repository-map.md` atualizados para o estado
  real do repositório; docstring de `src/mms_viscoelastic/ml/__init__.py`
  com a mesma nota de escopo.
- Caminho log-conformation portado de `conform.f` (Etapa 15, não discutido
  no texto da tese — portado por fidelidade ao código-fonte):
  `models/log_conformation.py` (`decompose_symmetric_2x2`, forma fechada da
  decomposição espectral de matriz simétrica 2x2 que reproduz exatamente a
  única rotação de Jacobi executada por `decomp_matrix` para `dimen=2`;
  `stress_to_log_conformation`/`log_conformation_to_stress`, port de
  `T_to_Psi`/`Psi_to_T`; `source_term_to_log_conformation`, port de
  `tf_conformation`; `kernel_conformation`, port de `kernel_conformation`) e
  `numerics/time_stepping_log.py` (`log_conformation_time_derivatives`,
  port de `drv_Psi`/`nlterms_logPsi`; `euler_step_log_conformation`, port do
  ramo `ty_sim=22`/`solver_log_euler`). Validado por precisão de máquina
  contra referência Fortran isolada que mantém o algoritmo de Jacobi
  genérico (não a forma fechada), em
  `validation/fortran_reference/case12_log_euler_step/`
  (`compare_case12_log_euler_step_reference`,
  `tests/test_case12_log_conformation_fortran_reference.py`), mais testes
  unitários cruzados contra `numpy.linalg.eigh`
  (`tests/test_model_log_conformation.py`,
  `tests/test_time_stepping_log.py`). A variante RK4 (`ty_sim=2`,
  `solver_log_rk4`) e a marcha temporal multi-passo completa ficam
  deferidas. Documentado em `docs/case-12-log-conformation.md`.
- Cases 13 e 14 portados como módulos triviais (décimo segundo e
  décimo terceiro/último da Etapa 14): `manufactured_solutions/case_13.py`,
  `case_14.py`, `source_terms/case_13.py`, `case_14.py`. Investigação
  revelou que, diferente dos Casos 01-12, esses dois **não são casos
  MMS** — a solução exata e o termo fonte são identicamente zero em
  ambos no Fortran original, e por isso não têm notebook Mathematica
  correspondente. São, na verdade, configurações de benchmark físico
  da cavidade com tampa móvel (Case 13: velocidade constante da tampa;
  Case 14: mesma config. com rampa suave via `tanh`). Após confirmação
  com o autor, portados como módulos triviais (retornam campos/termos
  zerados) apenas por paridade estrutural com os outros 12 casos; sem
  validação cruzada Fortran (ambos os lados seriam trivialmente zero).
  A condição de contorno da tampa (único conteúdo não-trivial) não foi
  portada, seguindo a mesma convenção dos demais casos. Com isso, a
  Etapa 14 está completa: 11 casos MMS genuínos (01-11) totalmente
  validados contra referência Fortran isolada, mais 2 casos triviais
  (13-14) por paridade estrutural. Documentado em
  `docs/case-13-14-port.md`.
- MMS Case 11 portado de `case_11.f` (décimo primeiro dos 13 casos
  restantes, Etapa 14): `manufactured_solutions/case_11.py`
  (`evaluate_case_11_exact`, `Case11Parameters`) e `source_terms/case_11.py`
  (`evaluate_case_11_tfwz`/`tftxx`/`tftxy`/`tftyy`,
  `Case11SourceTermParameters`). `u`,`v`,`wz`,`psi` idênticos ao Case 10
  (e ao Case 12); tensões polinomiais quadráticas, igual ao Case 02.
  Bateu precisão de máquina contra a referência Fortran isolada
  (`validation/fortran_reference/case11_fields/`) de primeira, nas 11
  quantidades exportadas, em duas malhas (5x5, 9x9). Documentado em
  `docs/case-11-port.md`.
- MMS Case 10 portado de `case_10.f` (décimo dos 13 casos restantes,
  Etapa 14): `manufactured_solutions/case_10.py`
  (`evaluate_case_10_exact`, `Case10Parameters`) e `source_terms/case_10.py`
  (`evaluate_case_10_tfwz`/`tftxx`/`tftxy`/`tftyy`,
  `Case10SourceTermParameters`). `u`,`v`,`wz`,`psi` idênticos ao Case 12
  (ansatz trigonométrico-vezes-polinomial); tensões uniformes no espaço,
  idênticas ao Case 01. Bateu precisão de máquina contra a referência
  Fortran isolada (`validation/fortran_reference/case10_fields/`) de
  primeira, nas 11 quantidades exportadas, em duas malhas (5x5, 9x9).
  Documentado em `docs/case-10-port.md`.
- MMS Case 09 portado de `case_09.f` (nono dos 13 casos restantes,
  Etapa 14, o maior e mais denso de todo o projeto — 790 linhas):
  `manufactured_solutions/case_09.py` (`evaluate_case_09_exact`,
  `Case09Parameters`) e `source_terms/case_09.py`
  (`evaluate_case_09_tfwz`/`tftxx`/`tftxy`/`tftyy`,
  `Case09SourceTermParameters`). `u`,`v`,`wz`,`psi` idênticos ao Case 08;
  tensões trigonométricas, igual ao Case 03/06/12. A combinação do
  jacobiano convectivo explícito (Case08) com o ansatz trigonométrico
  (Case06) torna `tfTxy` a expressão mais densa de todo o projeto (296
  pares de parênteses). Nenhuma parte foi retranscrita manualmente —
  módulo Python montado por script a partir dos termos mecânicos,
  normalizado via `ruff format`. Bateu precisão de máquina contra a
  referência Fortran isolada
  (`validation/fortran_reference/case09_fields/`) de primeira, nas 11
  quantidades exportadas, em duas malhas (5x5, 9x9), mesmo nessa escala
  inédita. Documentado em `docs/case-09-port.md`.
- MMS Case 08 portado de `case_08.f` (oitavo dos 13 casos restantes,
  Etapa 14, um dos dois maiores e mais densos — 678 linhas):
  `manufactured_solutions/case_08.py` (`evaluate_case_08_exact`,
  `Case08Parameters`) e `source_terms/case_08.py`
  (`evaluate_case_08_tfwz`/`tftxx`/`tftxy`/`tftyy`,
  `Case08SourceTermParameters`). `u`,`v`,`wz`,`psi` decaem no tempo via
  `tanh(4-8*at*t)` como no Case 07, mas com estrutura espacial puramente
  polinomial; tensões polinomiais quadráticas, igual ao Case 02. O
  jacobiano convectivo não-linear da solução manufaturada aparece
  explícito nos termos fonte (em vez de se anular), tornando `tfTxy` a
  expressão mais densa portada até agora (212 pares de parênteses).
  Nenhuma parte foi retranscrita manualmente — módulo Python montado
  por script a partir dos termos mecânicos, depois normalizado via
  `ruff format` (só formatação, nunca semântica). Bateu precisão de
  máquina contra a referência Fortran isolada
  (`validation/fortran_reference/case08_fields/`) de primeira, nas 11
  quantidades exportadas, em duas malhas (5x5, 9x9). Documentado em
  `docs/case-08-port.md`.
- MMS Case 07 portado de `case_07.f` (sétimo dos 13 casos restantes,
  Etapa 14): `manufactured_solutions/case_07.py`
  (`evaluate_case_07_exact`, `Case07Parameters`) e `source_terms/case_07.py`
  (`evaluate_case_07_tfwz`/`tftxx`/`tftxy`/`tftyy`,
  `Case07SourceTermParameters`). Único caso cujos campos de
  velocidade/vorticidade/função de corrente decaem no tempo via
  `tanh(4-8*at*t)` em vez de `exp(-at*t)`; tensões uniformes no espaço,
  idênticas ao Case 01. Introduziu dois casos novos para os scripts
  mecânicos: funções hiperbólicas (`Tanh`/`Cosh` → `np.tanh`/`np.cosh`)
  e literais Fortran com sufixo de expoente `d0` (`8.0d0` → `8.0e0`).
  A validação cruzada revelou um artefato de ordem de avaliação em
  ponto flutuante do gfortran (não um erro de transcrição — verificado
  via comparação byte a byte e teste com `-O0`/`-O2`/`-ffp-contract=off`)
  em `tfTyy`: a expressão Fortran de referência foi decomposta em
  variáveis intermediárias nomeadas (mesma técnica já usada no port
  Python) para eliminar o artefato. Bateu precisão de máquina contra a
  referência Fortran isolada (`validation/fortran_reference/case07_fields/`)
  nas 11 quantidades exportadas, em duas malhas (5x5, 9x9). Documentado
  em `docs/case-07-port.md`.
- MMS Case 06 portado de `case_06.f` (sexto dos 13 casos restantes,
  Etapa 14): `manufactured_solutions/case_06.py`
  (`evaluate_case_06_exact`, `Case06Parameters`) e `source_terms/case_06.py`
  (`evaluate_case_06_tfwz`/`tftxx`/`tftxy`/`tftyy`,
  `Case06SourceTermParameters`). `u`,`v`,`wz`,`psi` seguem o ansatz
  Taylor-Green dos Casos 04/05; o campo de tensões usa o ansatz
  trigonométrico do Case 03/Case 12 — primeira combinação
  Taylor-Green+trigonométrico entre os casos portados até agora. Todo
  termo (mesmo o `tfwz` compacto) foi colado literalmente do script
  mecânico Fortran→Python. Bateu precisão de máquina contra a referência
  Fortran isolada (`validation/fortran_reference/case06_fields/`) de
  primeira, nas 11 quantidades exportadas, em duas malhas (5x5, 9x9).
  Documentado em `docs/case-06-port.md`.
- MMS Case 05 portado de `case_05.f` (quinto dos 13 casos restantes,
  Etapa 14): `manufactured_solutions/case_05.py`
  (`evaluate_case_05_exact`, `Case05Parameters`) e `source_terms/case_05.py`
  (`evaluate_case_05_tfwz`/`tftxx`/`tftxy`/`tftyy`,
  `Case05SourceTermParameters`). `u`,`v`,`wz`,`psi` seguem o mesmo ansatz
  Taylor-Green do Case 04; o campo de tensões é polinomial quadrático,
  igual ao Case 02 (reaproveita `Case12Coefficients`). Aplicando a lição
  do Case 04, todo termo (mesmo o compacto `tfwz`) foi colado literalmente
  do script mecânico Fortran→Python, sem retranscrição manual em nenhum
  ponto — bateu precisão de máquina contra a referência Fortran isolada
  (`validation/fortran_reference/case05_fields/`) de primeira, nas 11
  quantidades exportadas, em duas malhas (5x5, 9x9). Documentado em
  `docs/case-05-port.md`.
- MMS Case 04 portado de `case_04.f` (quarto dos 13 casos restantes,
  Etapa 14): `manufactured_solutions/case_04.py`
  (`evaluate_case_04_exact`, `Case04Parameters`) e `source_terms/case_04.py`
  (`evaluate_case_04_tfwz`/`tftxx`/`tftxy`/`tftyy`,
  `Case04SourceTermParameters`). `u`,`v`,`wz`,`psi` seguem um ansatz tipo
  Taylor-Green (trigonométrico, sem fator polinomial); o campo de tensões
  é uniforme no espaço, idêntico ao do Case 01. Por serem termos fonte
  bem mais compactos que os Casos 02/03, foram transcritos manualmente em
  vez de copiados literalmente do script mecânico — e isso reintroduziu
  exatamente o tipo de erro que a cópia literal evita: um erro de sinal
  (`(1 - betann)` em vez de `(-1 + betann)`, sinais opostos) nos três
  termos de tensão. A validação cruzada contra a referência Fortran
  isolada pegou o erro de imediato (erro absoluto constante ~100-200 em
  todo ponto da malha, independente de x/y — assinatura de erro de sinal
  em termo que não depende de posição); corrigido restaurando
  `(-1.0 + betann)`. Documentado em `docs/case-04-port.md`. Validado
  contra novo programa Fortran de referência isolado
  (`validation/fortran_reference/case04_fields/`) nas 11 quantidades
  exportadas, a precisão de máquina.
- MMS Case 03 portado de `case_03.f` (terceiro dos 13 casos restantes,
  Etapa 14): `manufactured_solutions/case_03.py`
  (`evaluate_case_03_exact`, `Case03Parameters`) e `source_terms/case_03.py`
  (`evaluate_case_03_tfwz`/`tftxx`/`tftxy`/`tftyy`,
  `Case03SourceTermParameters`). `u`,`v`,`wz`,`psi` são idênticos aos Casos
  01/02 (verificado com SymPy); o campo de tensões usa o mesmo ansatz
  trigonométrico do Case 12, reaproveitando `Case12Coefficients`. Dada a
  densidade das expressões de `case_03.f` (até 225 pares de parênteses por
  termo), tanto a referência Fortran isolada quanto o módulo Python foram
  gerados por scripts mecânicos de tradução (preservam parênteses do
  Fortran original caractere a caractere; o script Fortran→Python divide
  o resultado em termos aditivos de topo, copiados sem retranscrição
  manual), em vez de transcritos à mão termo a termo. As duas
  transcrições mecânicas independentes concordam entre si e com a
  referência a precisão de máquina nas 11 quantidades exportadas, em duas
  malhas (5x5, 9x9). Documentado em `docs/case-03-port.md`.
- MMS Case 02 portado de `case_02.f` (segundo dos 13 casos restantes,
  Etapa 14): `manufactured_solutions/case_02.py`
  (`evaluate_case_02_exact`, `Case02Parameters`) e `source_terms/case_02.py`
  (`evaluate_case_02_tfwz`/`tftxx`/`tftxy`/`tftyy`,
  `Case02SourceTermParameters`). `u`,`v`,`wz`,`psi` são idênticos ao Case
  01; o campo de tensões é polinomial quadrático em x,y, reaproveitando
  `Case12Coefficients`. Simplificação algébrica do termo de relaxação
  (`Axx*Bxx+Ayy*Byy`) verificada com SymPy. Os quatro termos fonte da
  referência Fortran foram gerados por um script de tradução mecânica
  (preserva parênteses do Fortran original caractere a caractere), que
  pegou um parêntese extra introduzido na primeira transcrição manual de
  `tfTxy`; a contrapartida Python tinha um `Wi*` faltando em `tfTyy`,
  também corrigido. Validado contra novo programa Fortran de referência
  isolado (`validation/fortran_reference/case02_fields/`) nas 11
  quantidades exportadas, a precisão de máquina. Documentado em
  `docs/case-02-port.md`.
- MMS Case 01 portado de `case_01.f` (primeiro dos 13 casos restantes,
  Etapa 14): `manufactured_solutions/case_01.py`
  (`evaluate_case_01_exact`, `Case01Parameters`) e `source_terms/case_01.py`
  (`evaluate_case_01_tfwz`/`tftxx`/`tftxy`/`tftyy`,
  `Case01SourceTermParameters`). Diferente do Case 12, o campo de tensões
  manufaturado do Case 01 é espacialmente uniforme (`c1=1`, `c2=2`, `c3=3`
  de `par.nn`, escalados por `(1-beta)`); só `u`,`v`,`wz`,`psi` variam com a
  posição. Identidades de função de corrente (`u=dpsi/dy`, `v=-dpsi/dx`,
  `wz=du/dy-dv/dx`) verificadas simbolicamente com SymPy. Validado contra
  novo programa Fortran de referência isolado
  (`validation/fortran_reference/case01_fields/`) nas 11 quantidades
  exportadas (campos + termos fonte), a `max_abs_error < 1e-10` (precisão
  de máquina observada). Documentado em `docs/case-01-port.md`.
- Esquema compacto de 4ª ordem (Padé) e filtro de Lele, portados de
  `derivs.f` (`coef_4a`, `coeff_4a`, `rhsx_4a`, `rhsxx_4a`, `rhsy_4a`,
  `rhsyy_4a`), `solvers.f` (`tridx`, `tridy`) e `loop.f`
  (`rhs_tridfx_wz`/`rhs_tridfy_wz`, `lhs_tridfx`/`lhs_tridfy`), substituindo
  o placeholder de 2ª ordem usado até aqui:
  - `numerics/derivatives.py`: `solve_tridiagonal` (algoritmo de Thomas),
    `compact_derivative_x`/`compact_derivative_y`,
    `compact_second_derivative_x`/`compact_second_derivative_y`;
  - `numerics/filters.py` (novo): `lele_filter_x`/`lele_filter_y`,
    `filter_coefficients`, com `alphaf=0.48` padrão de `par.nn`;
  - verificação simbólica (SymPy) de que o fechamento de contorno da 2ª
    derivada compacta é de 3ª ordem (uma ordem abaixo do interior de 4ª
    ordem), enquanto o fechamento da 1ª derivada é de 4ª ordem completo;
    documentado em `validation/fortran_reference/case_compact_derivatives/README.md`
    e nos docstrings dos operadores afetados;
  - dois novos programas Fortran de referência isolados
    (`validation/fortran_reference/case_compact_derivatives/`,
    `validation/fortran_reference/case_lele_filter/`) e testes de regressão
    comparando o port Python contra eles (`max_abs_error < 1e-10`);
  - `numerics/convergence.py` reapontado para os operadores compactos; o
    estudo de convergência do Case 12 passa a recuperar ordem ~4.0-4.5
    (antes ~2.0), consistente com a tese e o artigo JNNFM.
  - `numerics/poisson.py`, `numerics/time_stepping.py` e
    `numerics/time_marching.py` continuam usando os operadores de 2ª ordem
    deliberadamente (documentado nos próprios programas Fortran de
    referência existentes); portá-los para o esquema compacto exige também
    reescrever esses programas de referência e o solver de Poisson
    pentadiagonal (`band5_poi`/`banbk5_poi`), deixado para uma etapa
    futura.
- `LICENSE` (MIT) na raiz do repositório, cobrindo a biblioteca Python,
  exemplos, testes e documentação novos; o material Fortran/Mathematica
  original em `src/fortran_original/`, `src/mathematica_original/` e
  `Thesis/` mantém os termos de `Thesis/LICENSE` (BSD-3-Clause,
  Copyright (c) 2015, Humberto Lidio Antonelli). Nova seção "License" no
  `README.md` esclarece essa distinção; necessário para o GitHub detectar
  a licença do repositório. O `LICENSE` contém apenas o texto padrão da MIT
  (sem anexos), para que o detector de licenças do GitHub o reconheça como
  "MIT" em vez de "Other".
- `numerics/poisson.py` com `solve_poisson_dirichlet`, um solver direto
  (`scipy.sparse` + `spsolve`) para `laplacian(phi) = source` com condições
  de contorno de Dirichlet, usando o operador de cinco pontos de segunda
  ordem.
- `numerics/poisson.py` com `run_case12_poisson_convergence_study`,
  verificando `solve_poisson_dirichlet` contra a MMS do Case 12
  (`laplacian(psi) = wz`, BC = `psi_exact`), com ordem observada ≈ 2 para
  `L2` e `Linf`.
- `examples/12_case12_poisson_solver.py`.
- `tests/test_case12_poisson_solver.py`.
- `docs/case-12-poisson-solver.md`, documentando a discretização, a
  convenção de sinal `wz = laplacian(psi)` e a simplificação de não portar o
  solver multigrid do Fortran (`poisson.f`).
- `scipy>=1.11` adicionado como dependência principal (antes apenas em `ml`).

- Programa Fortran isolado de referência para os campos exatos do Case 12
  (`u`, `v`, `wz`, `psi`, `txx`, `txy`, `tyy`) no instante `t=0`, porta
  verbatim de `exact_solutions_case_12(t)`.
- CSVs de referência 3x3 e 5x5 para os campos exatos do Case 12.
- `validation/solver_comparison.py` com `CASE12_FIELD_COMPONENTS` e
  `compare_case12_field_reference`, comparando `evaluate_case_12_exact`
  contra a referência Fortran.
- `examples/11_case12_solver_comparison.py`.
- `tests/test_case12_solver_comparison.py`, cobrindo os sete campos para as
  referências 3x3 e 5x5 (`max_abs_error < 1e-10`, `max_rel_error < 1e-10`).
- `docs/case-12-solver-comparison.md`, documentando o escopo reduzido da
  FASE VII (comparação dos campos exatos no instante t=0) e o trabalho futuro
  de comparação completa do solver.

- `numerics/derivatives.py` com `second_derivative_x` e `second_derivative_y`
  (derivada segunda de segunda ordem, com fronteiras zeradas).
- `numerics/time_stepping.py` com `EulerStepParameters`, `wz_time_derivative`,
  `stress_time_derivatives`, `euler_step_wz` e `euler_step_stress`, porta
  verbatim de `subroutine drv` (`drvs.f`) e do passo de Euler explícito de
  `subroutine solver_euler` (`euler_rk.f`) para `wz`, `Txx`, `Txy` e `Tyy`.
  Apenas os pontos internos da malha são atualizados; a reimposição das
  condições de contorno e o recálculo de `psi`/`ux`/`uy` em `t=dt` ficam para
  a Etapa 10.
- Programa Fortran isolado de referência para um passo de Euler do Case 12
  (`validation/fortran_reference/case12_euler_step/`), com gerador de CSV de
  entrada em Python (`generate_input_csv.py`) a partir dos campos exatos e
  termos fonte já validados, evitando re-transcrever `tfTxx`/`tfTxy`/`tfTyy`.
- CSV de referência 9x9 para o passo de Euler do Case 12
  (`dt = 1e-5`).
- `validation/solver_comparison.py` com `CASE12_EULER_STEP_COMPONENTS` e
  `compare_case12_euler_step_reference`.
- `examples/13_case12_euler_step_comparison.py`.
- `tests/test_case12_euler_step.py`, cobrindo `wz`, `Txx`, `Txy`, `Tyy`
  (`max_abs_error < 1e-12`, `max_rel_error < 1e-12`).
- `docs/case-12-euler-step-comparison.md`, documentando a formulação do passo
  de Euler, a escolha de parâmetros (`Re=1`, `Wi=1`, `beta=0.5`, `dt=1e-5`) e
  o escopo (discretização de segunda ordem, atualização apenas dos pontos
  internos).

- `numerics/time_marching.py` com `TimeMarchingState` e `case12_time_step`,
  combinando o passo de Euler (Etapa 9) com a reimposição das condições de
  contorno de Dirichlet em `t+dt`, a nova resolução de
  `laplacian(psi) = wz` (Etapa 8) e o recálculo de `u`/`v` a partir do `psi`
  atualizado. A condição de contorno de `wz` é reimposta a partir da solução
  exata em `t+dt` (em vez do stencil unilateral de 4ª ordem de
  `boundarie_vorticity_cavity_flow.f`), simplificação documentada e
  consistente com as Etapas 8 e 9.
- Programa Fortran isolado de referência para um passo completo de marcha no
  tempo do Case 12 (`validation/fortran_reference/case12_time_step/`):
  passo de Euler + reimposição de condições de contorno em `t+dt` + resolução
  densa do sistema de Poisson (eliminação gaussiana) + recálculo de `u`/`v`.
- Gerador de CSV de entrada com a coluna `psi`
  (`validation/fortran_reference/case12_time_step/generate_input_csv.py`) e
  CSV de referência 9x9 para o passo completo do Case 12 (`dt = 1e-5`).
- `validation/solver_comparison.py` com `CASE12_TIME_STEP_COMPONENTS` e
  `compare_case12_time_step_reference`.
- `examples/14_case12_time_step_comparison.py`.
- `tests/test_case12_time_step.py`, cobrindo `wz`, `Txx`, `Txy`, `Tyy`, `psi`,
  `u`, `v` (`max_abs_error < 1e-12`, `max_rel_error < 1e-12`).
- `docs/case-12-time-marching-comparison.md`, documentando a sequência
  portada (Euler + condições de contorno + Poisson + `u`/`v`), a decisão de
  reimpor `wz` na fronteira a partir da solução exata e o trabalho futuro de
  marcha multi-passo.

- `numerics/time_marching.py` com `case12_multistep_march`, encadeando
  `case12_time_step` (Etapa 10) por `n_steps` incrementos de `dt`.
- Programa Fortran isolado de referência para a marcha multi-passo do Case 12
  (`validation/fortran_reference/case12_multistep/`), reaproveitando o CSV de
  entrada `t=0` da Etapa 10 a cada passo (válido para o parâmetro padrão
  `at=0`, sob o qual os campos exatos e os termos fonte manufaturados são
  independentes do tempo).
- CSV de referência 9x9 para 10 passos (`dt = 1e-5`).
- `validation/solver_comparison.py` com `compare_case12_multistep_reference`.
- `examples/15_case12_multistep_comparison.py`.
- `tests/test_case12_multistep_march.py`, cobrindo `wz`, `Txx`, `Txy`, `Tyy`,
  `psi`, `u`, `v` após 10 passos (`max_abs_error < 1e-12`,
  `max_rel_error < 1e-12`).
- `docs/case-12-multistep-marching-comparison.md`, documentando
  `case12_multistep_march`, a independência temporal dos campos/termos fonte
  do Case 12 para `at=0` e o trabalho futuro de comparação da evolução
  temporal completa.
- Programa Fortran isolado de referência para a marcha completa do Case 12
  até o critério de convergência de `solver_euler`
  (`validation/fortran_reference/case12_full_march/`), com detecção de
  divergência (`NaN`/overflow) e histórico opcional por passo
  (`errowz`/`erroTxx`/`erroTxy`/`erroTyy`/`emv` + máximos absolutos dos
  campos).
- CSVs de referência de entrada e histórico para as malhas `9x9` e `17x17`
  (`validation/reference_data/case12_full_march_*`).
- `docs/case-12-full-march-reference.md` (Etapa 12 da FASE IX), documentando
  o achado de que a referência de 2ª ordem diverge em `t~4.5506` (malha
  `9x9`) em vez de convergir, e o re-escopo das Etapas 13-15.

### Corrigido
- `CITATION.cff`: `repository-code` apontava para um placeholder
  (`SEU_USUARIO/SEU_REPOSITORIO`); corrigido para a URL real do repositório.
- `CONTRIBUTING.md`: regras de preservação referenciavam `thesis/` (minúsculo,
  inexistente); corrigido para `Thesis/` e para incluir
  `src/mathematica_original/`.

### Alterado
- `README.md`: seção "Repository structure" corrigida para refletir os
  diretórios reais (`Thesis/`, `Papers/<...>/`, e os subpacotes `models/` e
  `io/` de `src/mms_viscoelastic`).
- `README.md`: "Current Python status" e "Pending work" atualizados para
  refletir os modelos constitutivos (UCM, Oldroyd-B) e o solver de marcha no
  tempo do Case 12 (Etapas 8-11).
- `README.md`: nova seção "Time-marching solver (Case 12)" com uma tabela
  resumindo o Poisson solver, o passo de Euler, o passo completo e a marcha
  multi-passo, com links para a documentação de cada etapa.
- `README.md`: seção "Local Git workflow" reescrita como "GitHub workflow",
  refletindo o fluxo real baseado em pull requests
  (`main → branch → commits → push → PR → merge via GitHub → main`).
- `README.md`: seção "GitHub About" atualizada para refletir a descrição e os
  tópicos já configurados no repositório, em vez de apresentá-los como
  sugestão.
- `README.md`: seção "Citation" atualizada para referenciar o `CITATION.cff`
  existente e o recurso "Cite this repository" do GitHub.
- `docs/local-git-workflow.md` e `docs/github-about.md` atualizados para o
  mesmo fluxo baseado em PR e para os tópicos/descrição configurados.

### Removido
- 22 artefatos de build (`__pycache__/*.pyc`) versionados por engano,
  já cobertos por `.gitignore`.

## [0.2.0] - 2026-06-13

### Corrigido
- `models/constitutive.py`: `ConstitutiveModel` agora usa `str, Enum` em vez de
  `enum.StrEnum` (incompatível com Python 3.10).
- `ml/features.py`: `field_to_feature_table` agora estreita corretamente o tipo
  de dataclass antes de chamar `asdict`.
- `visualization/fields.py`: `extent` do `imshow` agora é uma tupla, conforme
  esperado pelos stubs do Matplotlib.
- `source_terms/case_12.py`: removida variável local não utilizada em
  `evaluate_case_12_tfTxy`.

### Adicionado
- `models/ucm.py` com `VelocityGradient2D`, `StressTensor2D`,
  `rate_of_deformation_tensor`, `upper_convected_derivative`, `ucm_residual`
  e `ucm_steady_simple_shear`, implementando o modelo constitutivo
  Upper-Convected Maxwell (UCM) adimensional e sua solução analítica de
  cisalhamento simples (Couette).
- `models/oldroyd_b.py` com `solvent_stress`, `total_stress`,
  `oldroyd_b_residual` e `oldroyd_b_steady_simple_shear`, implementando o
  modelo Oldroyd-B como decomposição `tau = tau_p + tau_s` sobre o UCM.
- `docs/constitutive-models.md` com a derivação rigorosa da derivada
  upper-convected em forma de componentes e das soluções analíticas de UCM e
  Oldroyd-B em cisalhamento simples, usadas na validação dos novos modelos.
- `tests/test_model_ucm.py` e `tests/test_model_oldroyd_b.py`, validando os
  modelos UCM e Oldroyd-B contra as soluções analíticas de Couette e
  verificando que o resíduo das equações constitutivas é nulo (`atol=1e-12`)
  para uma grade de `(Wi, beta, shear_rate)`, incluindo os limites
  newtonianos `beta=0` e `beta=1`.
- `ml/splits.py` com `stratified_train_validation_test_split`, dividindo um
  dataframe em train/validation/test estratificado por uma coluna (ex.:
  `case_id`), garantindo que cada split contenha amostras de todas as
  combinações de parâmetros.
- `ml/datasets.py`: `build_case12_parametric_dataset`, que varre combinações
  de `(Re, Wi, beta, t)` via `itertools.product`, combina os campos exatos do
  Case 12 com os quatro termos fonte (`tfwz`, `tfTxx`, `tfTxy`, `tfTyy`) em
  uma única tabela de features por combinação, e retorna um
  `ParametricDatasetBundle` com splits estratificados por `case_id`.
- `ml/datasets.py`: `save_parametric_dataset_bundle`, `ParametricDatasetBundle`,
  `PARAMETRIC_TARGET_COLUMNS` e `SOURCE_TERM_COLUMNS`.
- `infer_feature_columns` agora aceita `exclude_columns`, usado para remover
  identificadores como `case_id` das features sem removê-los do dataframe.
- `examples/10_case12_ml_dataset_parametric.py`, exemplo que gera o dataset
  paramétrico e exporta os CSVs de full/train/validation/test.
- `tests/test_ml_splits.py` e `tests/test_ml_parametric.py`, cobrindo a
  validação das frações do split estratificado, a cobertura de todas as
  combinações de parâmetros nos três splits e a reprodutibilidade via
  `random_state`.
- Job `lint` no CI (`ruff check` + `mypy`) executado antes da matriz de testes.
- `numerics/convergence.py` com `run_convergence_study`, verificando a ordem
  de convergência (L2 e Linf) das relações `u = d(psi)/dy`, `v = -d(psi)/dx`
  e `wz = d(u)/dy - d(v)/dx` para o Case 12.
- `examples/09_case12_convergence_verification.py`.
- `tests/test_case12_convergence.py` (ordem L2 >= 1.8 e Linf >= 1.5 para
  grids 8->16->32->64).
- `docs/case-12-convergence-verification.md` com a verificação simbólica das
  identidades e a análise de erro de truncamento dos operadores.
- `docs/case-12-complete-validation-report.md`, consolidando a solução exata,
  a validação dos quatro termos fonte (tfwz, tfTxx, tfTxy, tfTyy) contra a
  referência Fortran e o estudo de convergência, incluindo a discussão da
  convenção de sinal de `wz` em relação à formulação geral do README.
- `examples/notebooks/Case12_MMS_Pipeline.ipynb`, notebook reprodutível que
  executa a solução exata, os termos fonte, a validação Fortran e o estudo
  de convergência do Case 12.
- Seção "Convergence results (Case 12)" no README, com a tabela de
  convergência e a nota sobre a convenção de sinal de `wz`.

### Alterado
- `source_terms/case_12.py`: docstring do módulo atualizada para refletir que
  os quatro termos fonte do Case 12 (tfwz, tfTxx, tfTxy, tfTyy) já estão
  validados contra a referência Fortran.
- README: status do `tfTyy` atualizado de "pending validation" para
  "validated"; seção "Pending work" atualizada.

---

## [0.6.7.5] - 2026-05-03

### Adicionado
- `assets/thesis/case12/` com figuras selecionadas da tese (campos exatos, superfícies de erro, curvas de convergência, slices numéricos).
- Seção "Thesis figures in this README" no README.
- `docs/thesis-figures-used-in-readme.md` com mapa de origem das figuras.

---

## [0.6.7.4] - 2026-05-03

### Adicionado
- README científico completo com: visão geral, formulação matemática, explicação do MMS, equações de verificação numérica, status do Case 12, validação dos termos fonte, figuras, comandos de instalação/teste/exemplos, fluxo Git local e roadmap.
- SVGs: `assets/preview/mms_validation_pipeline.svg`, `assets/preview/repository_architecture.svg`, `assets/preview/case12_source_terms_status.svg`.
- `docs/readme-maintenance-policy.md` - política de manutenção do README.

### Regra
O README não deve ser reduzido. Fases futuras podem expandir ou melhorar, nunca substituir por versão mínima.

---

## [0.6.7.3] - 2026-05-03

### Corrigido
- Restaurados `mean_abs_error` e `mean_rel_error` em `SourceTermReferenceComparison`.

---

## [0.6.7.2] - 2026-05-03

### Corrigido
- Restaurados `nx` e `ny` em `SourceTermReferenceComparison`.
- Reconstrução do grid de referência agora retorna `(grid, nx, ny)`.

---

## [0.6.7.1] - 2026-05-03

### Corrigido
- Removidos argumentos `x0`, `x1`, `y0`, `y1` da reconstrução do grid de validação (API `Grid2D.from_bounds` não os aceita).
- Adicionada validação explícita de que as coordenadas do CSV Fortran formam grids uniformes em `[0, 1]`.
- Testes de regressão para reconstrução do grid de referência.

---

## [0.6.7] - 2026-05-02

### Adicionado
- `evaluate_case_12_tfTyy` - quarto e último termo fonte do Case 12.
- `Case12SourceTerms` completo com `tfwz`, `tfTxx`, `tfTxy`, `tfTyy`.
- Registro atualizado de termos fonte.
- Testes de shape, valores finitos e referência padrão para `tfTyy`.
- `docs/case-12-tfTyy-port.md`.

### Status
Todos os quatro termos fonte do Case 12 estão implementados em Python.

---

## [0.6.6] - 2026-05-02

### Adicionado
- Programa Fortran isolado de referência para `tfTxy`.
- CSVs de referência 3×3 e 5×5 para `tfTxy`.
- Utilitário de validação Python para `tfTxy`.
- `examples/07_validate_case12_tfTxy_reference.py`.
- Testes de regressão pytest contra dados Fortran para `tfTxy`.
- `docs/case-12-tfTxy-fortran-reference-validation.md`.

---

## [0.6.5] - 2026-05-02

### Adicionado
- `evaluate_case_12_tfTxy`.
- Objeto parcial atualizado com `tfwz`, `tfTxx` e `tfTxy`.

### Corrigido
- Referência padrão de teste para `tfTxx` na origem corrigida para `2.0`.

---

## [0.6.4] - 2026-05-02

### Adicionado
- Programa Fortran isolado de referência para `tfTxx`.
- CSVs de referência 3×3 e 5×5 para `tfTxx`.
- Utilitários de validação Python para comparação CSV de termos fonte.
- `examples/06_validate_case12_tfTxx_reference.py`.
- Testes de regressão pytest contra dados Fortran para `tfTxx`.
- `docs/case-12-tfTxx-fortran-reference-validation.md`.

---

## [0.6.3] - 2026-05-02

### Adicionado
- `evaluate_case_12_tfTxx`.
- Objeto parcial atualizado com `tfwz` e `tfTxx`.

---

## [0.6.2] - 2026-05-02

### Adicionado
- Exportador Fortran isolado para `tfwz`.
- CSVs de referência 3×3 e 5×5 para `tfwz`.
- Utilitário de validação Python para comparar `tfwz` contra CSVs de referência.
- `examples/05_validate_case12_tfwz_reference.py`.
- Testes de regressão pytest contra dados Fortran para `tfwz`.
- `docs/case-12-tfwz-fortran-reference-validation.md`.

---

## [0.6.1] - 2026-05-02

### Adicionado
- `evaluate_case_12_tfwz` - primeiro termo fonte portado do Fortran.
- Suporte parcial de termos fonte via `allow_partial=True`.
- Registro explícito de componentes portados e pendentes.
- `examples/04_case12_tfwz_source_term.py`.
- Testes de shape, valores finitos, comportamento parcial e referência padrão para `tfwz`.
- `docs/case-12-tfwz-port.md`.

---

## [0.6.0] - 2026-05-02

### Adicionado
- `src/mms_viscoelastic/source_terms/__init__.py`
- `src/mms_viscoelastic/source_terms/case_12.py` - contrato público de termos fonte (levanta `NotImplementedError` até que as expressões sejam portadas).
- `tests/test_source_terms_case12_contract.py`.
- `docs/case-12-source-term-porting-map.md`.
- `docs/thesis-equations-overview.md`.

### Decisão de design
`evaluate_case_12_source_terms()` levanta `NotImplementedError` para evitar que zeros ou placeholders entrem silenciosamente em workflows de verificação.

---

## [0.5.4] - 2026-05-02

### Corrigido
- `Grid2D` agora instanciado via `Grid2D.from_bounds(nx=..., ny=...)`.
- `field_to_feature_table` preserva o padrão conservador: coordenadas + campos físicos.
- Features ML engineered são adicionadas apenas quando explicitamente solicitadas.
- Exportação com nomes retrocompatíveis: `plot_field`, `plot_field_panel`.
- Exemplo de visualização fecha figuras Matplotlib após salvar.

---

## [0.5.3] - 2026-05-02

### Adicionado
- `src/mms_viscoelastic/ml/features.py`.
- Separação entre extração de features de baixo nível e construtores de dataset.
- `tests/test_ml_imports.py`.
- `docs/ml-package-structure.md`.

### Corrigido
- Import circular em `mms_viscoelastic.ml`.

---

## [0.5.0] - 2026-05-02

### Adicionado
- README profissional com estrutura visual para o GitHub.
- Seção de preview do Case 12.
- `docs/github-about.md`.
- `docs/local-git-workflow.md`.
- `docs/case-12-preliminary-results.md`.
- Pasta `assets/preview/` com imagens para o README.

---

## [0.1.0] - 2026-05-02

### Adicionado
- Reorganização inicial do repositório.
- Preservação dos códigos Fortran originais em `src/fortran_original/`.
- Preservação dos notebooks Mathematica originais em `src/mathematica_original/`.
- Fundação da biblioteca Python `mms_viscoelastic` (parâmetros, malhas, campos, erros, exportação).
- Testes mínimos da camada Python.
- Documentação inicial do projeto.
