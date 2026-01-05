# Code Stastics (SCC CST 2026-01-05)

This is the code stastics report for the RSSN project.

This version is writen after the v0.2.2 release.

## Core Code

```text
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Language                              Files     Lines   Blanks  Comments     Code Complexity Complexity/Lines
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Rust                                    164    141765    26142     20293    95330       7235          1047.15
MaxLine / MeanLine                      225        22
(ULOC)                                          49187
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Total                                   164    141765    26142     20293    95330       7235          1047.15
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Unique Lines of Code (ULOC)                     49187
DRYness %                                        0.35
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Estimated Cost to Develop (semi-detached) $11,124,294
Estimated Schedule Effort (semi-detached) 27.94 months
Estimated People Required (semi-detached) 35.38
Processed 3287467 bytes, 3.287 megabytes (SI)
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
```

```text
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Language                              Files     Lines   Blanks  Comments     Code Complexity Complexity/Lines
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Rust                                    164    141765    26142     20293    95330       7235          1047.15
MaxLine / MeanLine                      225        22
(ULOC)                                          49187
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
src/symbolic/calculus.rs                         5272      834       283     4155        407             9.80
src/input/parser.rs                              4823       85        50     4688         25             0.53
src/symbolic/simplify_dag.rs                     3953      690       402     2861        380            13.28
src/symbolic/pde.rs                              3768      739       370     2659        261             9.82
src/symbolic/simplify.rs                         3471      508       255     2708        228             8.42
src/symbolic/graph_algorithms.rs                 3228      711       368     2149        301            14.01
src/symbolic/ode.rs                              2814      500       287     2027        200             9.87
src/numerical/matrix.rs                          2565      590       280     1695        226            13.33
src/symbolic/core/api.rs                         2546      458       333     1755         69             3.93
src/nightly/matrix.rs                            2546      578       277     1691        224            13.25
src/symbolic/core/expr_impl.rs                   2540      430        80     2030        114             5.62
src/symbolic/polynomial.rs                       2490      449       440     1601        151             9.43
src/numerical/computer_graphics.rs               2360      429       242     1689         48             2.84
src/symbolic/solve.rs                            2258      385       190     1683        141             8.38
src/numerical/physics.rs                         2183      493       414     1276        100             7.84
src/symbolic/matrix.rs                           2160      466       289     1405        193            13.74
src/numerical/error_correction.rs                2064      478       556     1030        153            14.85
~rc/numerical/fractal_geometry_and_chaos.rs      1921      409       485     1027         83             8.08
src/symbolic/special_functions.rs                1778      278       555      945        108            11.43
src/numerical/testing.rs                         1755      305       136     1314        123             9.36
src/symbolic/core/to_expr.rs                     1608      418        42     1148        137            11.93
src/symbolic/computer_graphics.rs                1599      172       264     1163         16             1.38
src/numerical/physics_fea.rs                     1499      304       205      990         33             3.33
src/symbolic/coordinates.rs                      1454      249       238      967         63             6.51
src/symbolic/special.rs                          1453      229       674      550         54             9.82
src/output/plotting.rs                           1422      275        68     1079         34             3.15
src/symbolic/number_theory.rs                    1392      277       160      955        120            12.57
src/numerical/physics_cfd.rs                     1378      335       213      830         59             7.11
src/symbolic/core/ast_impl.rs                    1374      147        21     1206         31             2.57
src/symbolic/error_correction_helper.rs          1331      283       238      810         99            12.22
src/numerical/physics_md.rs                      1315      296       241      778         58             7.46
src/symbolic/integration.rs                      1313      246       118      949         40             4.21
src/symbolic/topology.rs                         1299      265       239      795         78             9.81
src/physics/physics_rkm.rs                       1245      261        80      904         47             5.20
src/symbolic/poly_factorization.rs               1197      242       145      810         59             7.28
src/physics/physics_fdm.rs                       1175      277        84      814         78             9.58
src/symbolic/stats_probability.rs                1172      200        66      906          9             0.99
src/symbolic/combinatorics.rs                    1167      219       200      748         57             7.62
src/symbolic/logic.rs                            1145      202       109      834         84            10.07
src/symbolic/cas_foundations.rs                  1142      189        70      883         65             7.36
src/symbolic/transforms.rs                       1131      134       229      768         38             4.95
src/symbolic/grobner.rs                          1113      164       412      537         45             8.38
src/numerical/special.rs                         1112      334        75      703        101            14.37
src/physics/physics_fvm.rs                       1091      258        74      759         56             7.38
src/symbolic/tensor.rs                           1085      202       204      679         66             9.72
src/numerical/stats.rs                           1076      292       166      618         45             7.28
src/symbolic/finite_field.rs                     1020      185       122      713         51             7.15
src/output/pretty_print.rs                        943      202        10      731         65             8.89
src/output/io.rs                                  938      208       118      612         42             6.86
src/numerical/optimize.rs                         911      172        72      667         32             4.80
src/symbolic/unit_unification.rs                  898       94       119      685         28             4.09
src/symbolic/graph_operations.rs                  886      149       246      491         63            12.83
src/symbolic/core/dag_mgr.rs                      867      116       240      511         34             6.65
src/symbolic/core/expr.rs                         857       49       371      437          0             0.00
src/symbolic/error_correction.rs                  854      183       249      422         52            12.32
src/symbolic/complex_analysis.rs                  852      162       115      575         20             3.48
src/symbolic/series.rs                            849      126       164      559         32             5.72
src/compute/engine.rs                             832      103       327      402         11             2.74
src/plugins/manager.rs                            822      151       103      568         23             4.05
src/physics/physics_fem.rs                        804      190        66      548         56            10.22
src/physics/physics_mtm.rs                        797      184        52      561         46             8.20
src/symbolic/cryptography.rs                      797      136       175      486         41             8.44
src/numerical/differential_geometry.rs            780      166        77      537         67            12.48
src/symbolic/elementary.rs                        777      161        56      560         21             3.75
src/jit/engine.rs                                 755      130        43      582         51             8.76
src/symbolic/integral_equations.rs                739      100       149      490          7             1.43
src/physics/physics_sm.rs                         729      179        33      517         10             1.93
src/symbolic/geometric_algebra.rs                 687      138       115      434         42             9.68
src/symbolic/cad.rs                               685      156        49      480         45             9.38
src/numerical/sparse.rs                           660      169        90      401         38             9.48
~physics/physics_sim/navier_stokes_fluid.rs       633      140        82      411         46            11.19
src/physics/physics_em.rs                         632      157       102      373          9             2.41
src/numerical/geometric_algebra.rs                630       48        51      531          6             1.13
src/constant.rs                                   617       88        81      448          0             0.00
src/symbolic/quantum_mechanics.rs                 616      125        53      438          0             0.00
src/symbolic/vector.rs                            602       86       145      371          4             1.08
src/symbolic/rewriting.rs                         601      111        45      445         46            10.34
src/physics/physics_bem.rs                        579      129        75      375         24             6.40
src/physics/physics_cnm.rs                        573      166        40      367         33             8.99
src/symbolic/solid_state_physics.rs               571       77       100      394          0             0.00
src/numerical/polynomial.rs                       569      128        62      379         33             8.71
src/symbolic/fractal_geometry_and_chaos.rs        568      116        99      353         14             3.97
src/numerical/graph.rs                            563      148        77      338         45            13.31
src/numerical/integrate.rs                        556      110       167      279         28            10.04
src/symbolic/lie_groups_and_algebras.rs           555       97       139      319         13             4.08
src/numerical/elementary.rs                       555      152       113      290         23             7.93
src/numerical/finite_field.rs                     550      122        77      351         41            11.68
src/numerical/complex_analysis.rs                 550       78        93      379          3             0.79
src/symbolic/proof.rs                             549      117        49      383         38             9.92
src/physics/physics_mm.rs                         541      128        43      370         21             5.68
src/numerical/interpolate.rs                      530      114       145      271         31            11.44
src/numerical/vector.rs                           511      107       129      275         39            14.18
src/symbolic/differential_geometry.rs             508       77        67      364         20             5.49
src/numerical/tensor.rs                           496      103        63      330         25             7.58
src/symbolic/graph.rs                             490       93       108      289         14             4.84
src/symbolic/discrete_groups.rs                   490       98        58      334         26             7.78
src/numerical/vector_calculus.rs                  472      102       120      250         17             6.80
src/symbolic/real_roots.rs                        458       97        73      288         29            10.07
src/output/latex.rs                               451       85         9      357         17             4.76
src/symbolic/convergence.rs                       445       77       121      247         42            17.00
src/symbolic/multi_valued.rs                      444       74       121      249          0             0.00
src/verification/symbolic_core.rs                 439      162        19      258          6             2.33
src/symbolic/radicals.rs                          437       85        38      314         34            10.83
src/numerical/calculus.rs                         421       77       158      186         12             6.45
src/symbolic/handles.rs                           415       41       274      100          1             1.00
src/symbolic/group_theory.rs                      414       96        27      291         34            11.68
src/symbolic/stats_information_theory.rs          412       66       107      239         12             5.02
src/symbolic/classical_mechanics.rs               406       73        58      275          0             0.00
src/numerical/real_roots.rs                       406      106        65      235         54            22.98
src/numerical/topology.rs                         403      102        40      261         30            11.49
~c/physics/physics_sim/linear_elasticity.rs       401       86        50      265         12             4.53
src/lib.rs                                        399       10       295       94          1             1.06
src/symbolic/functional_analysis.rs               389       65        60      264          7             2.65
src/numerical/ode.rs                              386       72        60      254          7             2.76
src/numerical/convergence.rs                      384       91       103      190         25            13.16
~physics/physics_sim/schrodinger_quantum.rs       375       75        57      243         13             5.35
src/symbolic/optimize.rs                          369       63        66      240         19             7.92
src/symbolic/numeric.rs                           359       40        31      288         18             6.25
src/symbolic/relativity.rs                        358       66        39      253          0             0.00
src/numerical/number_theory.rs                    349      105        47      197         49            24.87
src/symbolic/stats_inference.rs                   346       44        34      268          0             0.00
src/symbolic/quantum_field_theory.rs              346       66        23      257          2             0.78
~symbolic/graph_isomorphism_and_coloring.rs       340       72        42      226         29            12.83
src/symbolic/thermodynamics.rs                    333       57        37      239          1             0.42
src/symbolic/electromagnetism.rs                  324       56        45      223          0             0.00
src/output/typst.rs                               323       62         4      257         12             4.67
src/symbolic/vector_calculus.rs                   321       51        62      208          0             0.00
~c/physics/physics_sim/ising_statistical.rs       310       76        35      199         18             9.05
src/numerical/combinatorics.rs                    297       78        58      161         35            21.74
src/numerical/transforms.rs                       297       68       121      108         19            17.59
src/numerical/coordinates.rs                      294       68        75      151         16            10.60
src/numerical/solve.rs                            290       63        46      181         23            12.71
~physics/physics_sim/geodesic_relativity.rs       280       53        52      175          3             1.71
~c/physics/physics_sim/gpe_superfluidity.rs       274       55        24      195          3             1.54
src/numerical/functional_analysis.rs              273       59        64      150          9             6.00
src/numerical/multi_valued.rs                     253       68        27      158         10             6.33
~hysics/physics_sim/fdtd_electrodynamics.rs       252       67        26      159         16            10.06
src/symbolic/stats_regression.rs                  248       46        52      150          5             3.33
src/symbolic/stats.rs                             217       42        59      116          7             6.03
src/symbolic/core/mod.rs                          217        4       201       12          0             0.00
src/compute/cache.rs                              212       39        54      119          2             1.68
src/numerical/signal.rs                           210       40        90       80         14            17.50
src/symbolic/calculus_of_variations.rs            204       22        97       85          0             0.00
src/numerical/calculus_of_variations.rs           193       20        92       81          0             0.00
src/plugins/plugin_c.rs                           181       33        53       95          2             2.11
src/compute/mod.rs                                156        1       150        5          0             0.00
src/numerical/series.rs                           150       35        50       65          3             4.62
src/jit/instructions.rs                           114       14        45       55          0             0.00
src/compute/computation.rs                        100       14        20       66          0             0.00
src/symbolic/mod.rs                                90        1        21       68          0             0.00
src/plugins/stable_abi.rs                          82       18        13       51          0             0.00
src/numerical/mod.rs                               56        1         8       47          0             0.00
src/physics/mod.rs                                 52        2        38       12          0             0.00
src/compute/state.rs                               36        7         7       22          1             4.55
src/output/mod.rs                                  27        1        21        5          0             0.00
src/compute/computable.rs                          27        3        15        9          0             0.00
src/numerical/pde.rs                               26        4        10       12          0             0.00
src/physics/physics_sim/mod.rs                     18        1        10        7          0             0.00
src/jit/mod.rs                                     16        3         2       11          0             0.00
src/plugins/mod.rs                                 15        1        11        3          0             0.00
src/verification/mod.rs                             8        1         6        1          0             0.00
src/ffi_blindings/mod.rs                            4        1         2        1          0             0.00
src/input/mod.rs                                    3        0         2        1          0             0.00
src/nightly/mod.rs                                  1        0         0        1          0             0.00
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Total                                   164    141765    26142     20293    95330       7235          1047.15
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Unique Lines of Code (ULOC)                     49187
DRYness %                                        0.35
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Estimated Cost to Develop (semi-detached) $11,124,294
Estimated Schedule Effort (semi-detached) 27.94 months
Estimated People Required (semi-detached) 35.38
Processed 3287467 bytes, 3.287 megabytes (SI)
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
```

## Main repository

```text
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Language                              Files     Lines   Blanks  Comments     Code Complexity Complexity/Lines
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Rust                                   1204    342810    67668     56233   218909      12844          4092.99
MaxLine / MeanLine                      278        23
(ULOC)                                          93534
-------------------------------------------------------------------------------------------------------------
Markdown                                 20      5554      837         0     4717          0             0.00
MaxLine / MeanLine                      655        51
(ULOC)                                           3374
-------------------------------------------------------------------------------------------------------------
TOML                                      8       385       18        15      352          2             1.13
MaxLine / MeanLine                      126        22
(ULOC)                                            326
-------------------------------------------------------------------------------------------------------------
Python                                    6       630       64        54      512         36            54.45
MaxLine / MeanLine                      128        37
(ULOC)                                            411
-------------------------------------------------------------------------------------------------------------
SVG                                       6      4589      108         2     4479        790           145.40
MaxLine / MeanLine                   324860       218
(ULOC)                                           1553
-------------------------------------------------------------------------------------------------------------
Shell                                     5       327       50        59      218         22            40.97
MaxLine / MeanLine                      114        32
(ULOC)                                            254
-------------------------------------------------------------------------------------------------------------
YAML                                      5       351       57        44      250          0             0.00
MaxLine / MeanLine                      293        37
(ULOC)                                            257
-------------------------------------------------------------------------------------------------------------
Batch                                     3       169       24        12      133         28            40.22
MaxLine / MeanLine                      101        27
(ULOC)                                            109
-------------------------------------------------------------------------------------------------------------
XML                                       3        27        0         0       27          0             0.00
MaxLine / MeanLine                      113        37
(ULOC)                                             18
-------------------------------------------------------------------------------------------------------------
Fortran Modern                            2       256       50        44      162          9            10.84
MaxLine / MeanLine                      105        35
(ULOC)                                            171
-------------------------------------------------------------------------------------------------------------
C Header                                  1     40700     2578     28171     9951          0             0.00
MaxLine / MeanLine                      142        33
(ULOC)                                           9446
-------------------------------------------------------------------------------------------------------------
C++                                       1        68       12         6       50          7            14.00
MaxLine / MeanLine                       81        25
(ULOC)                                             52
-------------------------------------------------------------------------------------------------------------
C++ Header                                1     40786     2587     28184    10015         22             0.22
MaxLine / MeanLine                      142        32
(ULOC)                                           9471
-------------------------------------------------------------------------------------------------------------
HTML                                      1         7        0         0        7          0             0.00
MaxLine / MeanLine                      113        37
(ULOC)                                              7
-------------------------------------------------------------------------------------------------------------
Handlebars                                1        16        3         0       13          0             0.00
MaxLine / MeanLine                      103        23
(ULOC)                                             14
-------------------------------------------------------------------------------------------------------------
License                                   1        73       32         0       41          0             0.00
MaxLine / MeanLine                      952       139
(ULOC)                                             42
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Total                                  1268    436748    74088    112824   249836      13760          4400.22
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Unique Lines of Code (ULOC)                    111913
DRYness %                                        0.26
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Estimated Cost to Develop (semi-detached) $32,727,207
Estimated Schedule Effort (semi-detached) 40.75 months
Estimated People Required (semi-detached) 71.34
Processed 12627225 bytes, 12.627 megabytes (SI)
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
```

```text
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Language                              Files     Lines   Blanks  Comments     Code Complexity Complexity/Lines
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Rust                                   1204    342810    67668     56233   218909      12844          4092.99
MaxLine / MeanLine                      278        23
(ULOC)                                          93534
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
src/ffi_apis/ffi_api.rs                          7963     1386      1154     5423        396             7.30
src/prelude.rs                                   6120      588       553     4979          0             0.00
src/symbolic/calculus.rs                         5272      834       283     4155        407             9.80
src/input/parser.rs                              4823       85        50     4688         25             0.53
src/symbolic/simplify_dag.rs                     3953      690       402     2861        380            13.28
src/symbolic/pde.rs                              3768      739       370     2659        261             9.82
src/symbolic/simplify.rs                         3471      508       255     2708        228             8.42
src/symbolic/graph_algorithms.rs                 3228      711       368     2149        301            14.01
tests/input_parser_test.rs                       2903      436        12     2455          0             0.00
src/symbolic/ode.rs                              2814      500       287     2027        200             9.87
src/numerical/matrix.rs                          2565      590       280     1695        226            13.33
src/symbolic/core/api.rs                         2546      458       333     1755         69             3.93
src/nightly/matrix.rs                            2546      578       277     1691        224            13.25
src/symbolic/core/expr_impl.rs                   2540      430        80     2030        114             5.62
src/symbolic/polynomial.rs                       2490      449       440     1601        151             9.43
src/numerical/computer_graphics.rs               2360      429       242     1689         48             2.84
src/symbolic/solve.rs                            2258      385       190     1683        141             8.38
src/numerical/physics.rs                         2183      493       414     1276        100             7.84
src/symbolic/matrix.rs                           2160      466       289     1405        193            13.74
src/numerical/error_correction.rs                2064      478       556     1030        153            14.85
~rc/numerical/fractal_geometry_and_chaos.rs      1921      409       485     1027         83             8.08
src/symbolic/special_functions.rs                1778      278       555      945        108            11.43
src/numerical/testing.rs                         1754      305       136     1313        123             9.37
src/symbolic/core/to_expr.rs                     1608      418        42     1148        137            11.93
src/symbolic/computer_graphics.rs                1599      172       264     1163         16             1.38
src/numerical/physics_fea.rs                     1499      304       205      990         33             3.33
src/symbolic/coordinates.rs                      1454      249       238      967         63             6.51
src/symbolic/special.rs                          1453      229       674      550         54             9.82
src/output/plotting.rs                           1422      275        68     1079         34             3.15
src/symbolic/number_theory.rs                    1392      277       160      955        120            12.57
src/numerical/physics_cfd.rs                     1378      335       213      830         59             7.11
src/symbolic/core/ast_impl.rs                    1374      147        21     1206         31             2.57
src/symbolic/error_correction_helper.rs          1331      283       238      810         99            12.22
src/ffi_apis/numerical_physics_ffi/json.rs       1321      162       529      630         17             2.70
src/numerical/physics_md.rs                      1315      296       241      778         58             7.46
src/symbolic/integration.rs                      1313      246       118      949         40             4.21
~is/numerical_computer_graphics_ffi/json.rs      1301      205       223      873         18             2.06
src/symbolic/topology.rs                         1299      265       239      795         78             9.81
src/ffi_apis/numerical_vector_ffi/json.rs        1276      162       193      921         45             4.89
src/physics/physics_rkm.rs                       1245      261        80      904         47             5.20
~i_apis/symbolic_cryptography_ffi/handle.rs      1227      233       475      519         50             9.63
src/symbolic/poly_factorization.rs               1197      242       145      810         59             7.28
tests/numerical_error_correction_test.rs         1193      331        82      780         11             1.41
src/physics/physics_fdm.rs                       1175      277        84      814         78             9.58
src/symbolic/stats_probability.rs                1172      200        66      906          9             0.99
src/symbolic/combinatorics.rs                    1167      219       200      748         57             7.62
src/symbolic/logic.rs                            1145      202       109      834         84            10.07
src/symbolic/cas_foundations.rs                  1142      189        70      883         65             7.36
src/symbolic/transforms.rs                       1131      134       229      768         38             4.95
src/symbolic/grobner.rs                          1113      164       412      537         45             8.38
src/numerical/special.rs                         1112      334        75      703        101            14.37
src/ffi_apis/numerical_special_ffi/json.rs       1095      111       476      508         16             3.15
src/physics/physics_fvm.rs                       1091      258        74      759         56             7.38
tests/numerical_computer_graphics_test.rs        1090      352        46      692         22             3.18
src/symbolic/tensor.rs                           1085      202       204      679         66             9.72
src/numerical/stats.rs                           1076      292       166      618         45             7.28
~is/symbolic_complex_analysis_ffi/handle.rs      1039      294       253      492         36             7.32
~c/ffi_apis/symbolic_transforms_ffi/json.rs      1038      283       250      505         22             4.36
~ffi_apis/symbolic_transforms_ffi/handle.rs      1035      220       290      525         38             7.24
src/symbolic/finite_field.rs                     1020      185       122      713         51             7.15
~apis/symbolic_complex_analysis_ffi/json.rs       991      259       251      481         35             7.28
~sts/numerical_error_correction_ffi_test.rs       969      273        30      666          0             0.00
~cal_fractal_geometry_and_chaos_ffi/json.rs       960      141       191      628         15             2.39
~_apis/numerical_special_ffi/bincode_api.rs       954       94       450      410         16             3.90
~pis/numerical_error_correction_ffi/json.rs       954      142       198      614         21             3.42
src/ffi_apis/numerical_stats_ffi/json.rs          954      105       391      458         13             2.84
~mbolic_complex_analysis_ffi/bincode_api.rs       949      232       277      440         35             7.95
src/output/pretty_print.rs                        943      202        10      731         65             8.89
tests/symbolic_special_functions_test.rs          941      224        99      618         22             3.56
src/output/io.rs                                  938      208       118      612         42             6.86
~pis/symbolic_transforms_ffi/bincode_api.rs       935      261       103      571         22             3.85
src/numerical/optimize.rs                         911      172        72      667         32             4.80
src/ffi_apis/symbolic_special_ffi/json.rs         902      212       243      447         56            12.53
src/symbolic/unit_unification.rs                  898       94       119      685         28             4.09
~s/symbolic_special_functions_ffi/handle.rs       897      190       233      474         51            10.76
src/symbolic/graph_operations.rs                  886      149       246      491         63            12.83
~pis/symbolic_special_functions_ffi/json.rs       882      214       233      435         33             7.59
~rc/ffi_apis/numerical_matrix_ffi/handle.rs       876      214       212      450         48            10.67
src/symbolic/core/dag_mgr.rs                      867      116       240      511         34             6.65
src/symbolic/core/expr.rs                         857       49       371      437          0             0.00
tests/symbolic_graph_algorithms_test.rs           857      214        20      623         11             1.77
src/symbolic/error_correction.rs                  854      183       249      422         52            12.32
src/symbolic/complex_analysis.rs                  852      162       115      575         20             3.48
~merical_fractal_geometry_and_chaos_test.rs       849      239        78      532         16             3.01
src/symbolic/series.rs                            849      126       164      559         32             5.72
~s/symbolic_computer_graphics_ffi/handle.rs       846      174       201      471         50            10.62
src/ffi_apis/nightly_ffi/handle.rs                845      223       172      450         48            10.67
src/compute/engine.rs                             832      103       327      402         11             2.74
~i_apis/numerical_vector_ffi/bincode_api.rs       830      119       129      582         30             5.15
~bolic_special_functions_ffi/bincode_api.rs       828      283       104      441         33             7.48
src/plugins/manager.rs                            822      151       103      568         23             4.05
~ffi_apis/numerical_physics_fea_ffi/json.rs       810       96       313      401          8             2.00
src/physics/physics_fem.rs                        804      190        66      548         56            10.22
~c/ffi_apis/numerical_physics_ffi/handle.rs       800      112       403      285          0             0.00
src/physics/physics_mtm.rs                        797      184        52      561         46             8.20
src/symbolic/cryptography.rs                      797      136       175      486         41             8.44
tests/numerical_physics_fea_test.rs               794      222        61      511          8             1.57
tests/numerical_physics_md_test.rs                790      179        53      558          2             0.36
src/numerical/differential_geometry.rs            780      166        77      537         67            12.48
~/ffi_apis/numerical_physics_md_ffi/json.rs       779       94       238      447         10             2.24
src/symbolic/elementary.rs                        777      161        56      560         21             3.75
tests/numerical_physics_test.rs                   770      203        79      488         24             4.92
tests/symbolic_special_test.rs                    768      197        84      487          3             0.62
~ffi_apis/numerical_physics_cfd_ffi/json.rs       763       87       309      367          7             1.91
src/jit/engine.rs                                 755      130        43      582         51             8.76
~c/ffi_apis/symbolic_calculus_ffi/handle.rs       753      147       185      421         49            11.64
tests/symbolic_error_correction_test.rs           745      224        46      475          1             0.21
~erical_error_correction_ffi/bincode_api.rs       740      125       134      481         21             4.37
src/symbolic/integral_equations.rs                739      100       149      490          7             1.43
src/ffi_apis/symbolic_pde_ffi/handle.rs           734      135       108      491         64            13.03
src/physics/physics_sm.rs                         729      179        33      517         10             1.93
~ffi_apis/symbolic_lie_groups_ffi/handle.rs       723      128       346      249         14             5.62
tests/symbolic_computer_graphics_test.rs          709      182        26      501         20             3.99
tests/numerical_special_test.rs                   695      152        86      457          1             0.22
tests/numerical_physics_cfd_test.rs               688      228        47      413         17             4.12
~i_apis/symbolic_group_theory_ffi/handle.rs       688      117       285      286          8             2.80
src/symbolic/geometric_algebra.rs                 687      138       115      434         42             9.68
src/symbolic/cad.rs                               685      156        49      480         45             9.38
~fi_apis/numerical_stats_ffi/bincode_api.rs       683       79       298      306         11             3.59
~i_apis/symbolic_special_ffi/bincode_api.rs       674      183        48      443         56            12.64
src/numerical/sparse.rs                           660      169        90      401         38             9.48
~s/numerical_error_correction_ffi/handle.rs       645      136       163      346         26             7.51
~ffi_apis/symbolic_cryptography_ffi/json.rs       645      136       118      391         32             8.18
~s/symbolic_quantum_mechanics_ffi/handle.rs       638      141       160      337         34            10.09
~c/ffi_apis/numerical_special_ffi/handle.rs       635       91       347      197          0             0.00
~s/symbolic_stats_probability_ffi/handle.rs       635      169       165      301         12             3.99
~physics/physics_sim/navier_stokes_fluid.rs       633      140        82      411         46            11.19
src/physics/physics_em.rs                         632      157       102      373          9             2.41
src/numerical/geometric_algebra.rs                630       48        51      531          6             1.13
~rc/ffi_apis/numerical_vector_ffi/handle.rs       625      152       123      350         34             9.71
tests/symbolic_finite_field_test.rs               622      106        22      494          0             0.00
tests/symbolic_pde_test.rs                        621      140        29      452          0             0.00
src/constant.rs                                   617       88        81      448          0             0.00
src/symbolic/quantum_mechanics.rs                 616      125        53      438          0             0.00
~ts/numerical_computer_graphics_ffi_test.rs       615      176        13      426          0             0.00
src/ffi_apis/numerical_stats_ffi/handle.rs        607      141       125      341         55            16.13
~is/symbolic_error_correction_ffi/handle.rs       607      125       165      317         29             9.15
tests/symbolic_logic_test.rs                      607      151        33      423          2             0.47
src/symbolic/vector.rs                            602       86       145      371          4             1.08
src/symbolic/rewriting.rs                         601      111        45      445         46            10.34
~cal_fractal_geometry_and_chaos_ffi_test.rs       600      167         8      425          1             0.24
src/physics/physics_bem.rs                        579      129        75      375         24             6.40
src/ffi_apis/nightly_ffi/json.rs                  576       99        61      416         26             6.25
src/ffi_apis/numerical_matrix_ffi/json.rs         576       99        61      416         26             6.25
src/physics/physics_cnm.rs                        573      166        40      367         33             8.99
tests/symbolic_numeric_test.rs                    572      157        27      388          0             0.00
src/symbolic/solid_state_physics.rs               571       77       100      394          0             0.00
~pis/symbolic_vector_calculus_ffi/handle.rs       571      121        16      434         52            11.98
tests/numerical_stats_ffi_test.rs                 570      123        25      422          0             0.00
src/numerical/polynomial.rs                       569      128        62      379         33             8.71
src/symbolic/fractal_geometry_and_chaos.rs        568      116        99      353         14             3.97
src/ffi_apis/symbolic_ode_ffi/handle.rs           567       90        96      381         51            13.39
~pis/symbolic_stats_probability_ffi/json.rs       567      163       162      242         10             4.13
~c/ffi_apis/symbolic_ode_ffi/bincode_api.rs       564      112         9      443         47            10.61
tests/symbolic_geometric_algebra_test.rs          564      130        22      412          0             0.00
src/numerical/graph.rs                            563      148        77      338         45            13.31
~is/numerical_geometric_algebra_ffi/json.rs       563       79        97      387         17             4.39
src/numerical/integrate.rs                        556      110       167      279         28            10.04
src/numerical/elementary.rs                       555      152       113      290         23             7.93
src/symbolic/lie_groups_and_algebras.rs           555       97       139      319         13             4.08
~symbolic_functional_analysis_ffi/handle.rs       551      104       225      222          4             1.80
src/numerical/finite_field.rs                     550      122        77      351         41            11.68
src/numerical/complex_analysis.rs                 550       78        93      379          3             0.79
src/symbolic/proof.rs                             549      117        49      383         38             9.92
~olic_error_correction_helper_ffi/handle.rs       549      112       116      321         25             7.79
src/ffi_apis/symbolic_ode_ffi/json.rs             543      112         9      422         47            11.14
src/physics/physics_mm.rs                         541      128        43      370         21             5.68
tests/symbolic_handles_test.rs                    538      155         8      375          9             2.40
~_apis/symbolic_calculus_ffi/bincode_api.rs       537      123        13      401         35             8.73
src/ffi_apis/symbolic_calculus_ffi/json.rs        534      123        13      398         35             8.79
src/numerical/interpolate.rs                      530      114       145      271         31            11.44
tests/symbolic_polynomial_test.rs                 529      122        34      373          2             0.54
tests/symbolic_elementary_test.rs                 519      144        44      331          0             0.00
tests/numerical_special_ffi_test.rs               519      126        27      366          0             0.00
tests/numerical_physics_ffi_test.rs               514      154         7      353         14             3.97
~i_apis/numerical_combinatorics_ffi/json.rs       514       82       109      323         10             3.10
~apis/symbolic_error_correction_ffi/json.rs       512      113       135      264         31            11.74
~s/symbolic_cryptography_ffi/bincode_api.rs       511      117        17      377         32             8.49
~c/ffi_apis/symbolic_lie_groups_ffi/json.rs       511       80       237      194         15             7.73
src/numerical/vector.rs                           511      107       129      275         39            14.18
src/symbolic/differential_geometry.rs             508       77        67      364         20             5.49
tests/symbolic_calculus_test.rs                   508      127        31      350          7             2.00
~is/symbolic_graph_algorithms_ffi/handle.rs       501      132        57      312         20             6.41
tests/numerical_graph_ffi_test.rs                 499      114        25      360          0             0.00
src/numerical/tensor.rs                           496      103        63      330         25             7.58
tests/symbolic_tensor_test.rs                     495      111        22      362          2             0.55
tests/symbolic_cryptography_test.rs               490      132        38      320          5             1.56
src/symbolic/discrete_groups.rs                   490       98        58      334         26             7.78
src/symbolic/graph.rs                             490       93       108      289         14             4.84
src/ffi_apis/physics_rkm_ffi/json.rs              488       68       142      278          8             2.88
~pis/symbolic_lie_groups_ffi/bincode_api.rs       479       72       249      158         15             9.49
~s/symbolic_group_theory_ffi/bincode_api.rs       477       70       252      155         15             9.68
src/numerical/vector_calculus.rs                  472      102       120      250         17             6.80
~ffi_apis/symbolic_group_theory_ffi/json.rs       469       71       246      152         15             9.87
~_apis/numerical_physics_ffi/bincode_api.rs       469       58       200      211          7             3.32
src/ffi_apis/symbolic_graph_ffi/handle.rs         464      121        38      305         23             7.54
src/ffi_apis/symbolic_proof_ffi/handle.rs         459       87        63      309         39            12.62
src/symbolic/real_roots.rs                        458       97        73      288         29            10.07
~/ffi_apis/numerical_optimize_ffi/handle.rs       458       83       128      247         12             4.86
~apis/symbolic_graph_algorithms_ffi/json.rs       457       96       130      231         15             6.49
examples/fem_structural_analysis.rs               456      101        27      328         17             5.18
~/numerical_computer_graphics_ffi/handle.rs       455      103        76      276         19             6.88
tests/symbolic_solve_test.rs                      455       88        19      348          4             1.15
src/output/latex.rs                               451       85         9      357         17             4.76
~bolic_stats_probability_ffi/bincode_api.rs       447      148        58      241         10             4.15
src/symbolic/convergence.rs                       445       77       121      247         42            17.00
src/ffi_apis/numerical_graph_ffi/json.rs          444       84        73      287          9             3.14
~ctal_geometry_and_chaos_ffi/bincode_api.rs       444       73        81      290          9             3.10
src/symbolic/multi_valued.rs                      444       74       121      249          0             0.00
tests/numerical_combinatorics_ffi_test.rs         443      113        15      315          0             0.00
src/verification/symbolic_core.rs                 439      162        19      258          6             2.33
~s/symbolic_error_correction_helper_test.rs       439      126        43      270          5             1.85
src/symbolic/radicals.rs                          437       85        38      314         34            10.83
~pis/symbolic_computer_graphics_ffi/json.rs       436       86       100      250         17             6.80
~fi_apis/symbolic_series_ffi/bincode_api.rs       435      110        39      286          8             2.80
~mbolic_graph_algorithms_ffi/bincode_api.rs       432       96       104      232         15             6.47
src/ffi_apis/physics_fdm_ffi/json.rs              430       49       145      236          9             3.81
tests/symbolic_lie_groups_test.rs                 429      119        23      287         15             5.23
~/ffi_apis/symbolic_rewriting_ffi/handle.rs       428       92       123      213         23            10.80
src/ffi_apis/symbolic_series_ffi/json.rs          428      111        39      278          8             2.88
src/ffi_apis/symbolic_series_ffi/handle.rs        428      139        39      250          8             3.20
tests/numerical_physics_fea_ffi_test.rs           424      128         7      289          4             1.38
~c_fractal_geometry_and_chaos_ffi/handle.rs       421      111        14      296         25             8.45
src/numerical/calculus.rs                         421       77       158      186         12             6.45
tests/symbolic_poly_factorization_test.rs         420      104        34      282          2             0.71
~/numerical_geometric_algebra_ffi/handle.rs       419      116        89      214         24            11.21
examples/polynomial_demo.rs                       415       74        23      318          0             0.00
src/symbolic/handles.rs                           415       41       274      100          1             1.00
src/symbolic/group_theory.rs                      414       96        27      291         34            11.68
tests/compute_engine_test.rs                      413      126        31      256          7             2.73
src/symbolic/stats_information_theory.rs          412       66       107      239         12             5.02
~mbolic_error_correction_ffi/bincode_api.rs       411       98        46      267         30            11.24
tests/symbolic_multi_valued_test.rs               409      128        28      253          0             0.00
src/numerical/real_roots.rs                       406      106        65      235         54            22.98
~/symbolic_stats_information_theory_test.rs       406       75        32      299         19             6.35
src/symbolic/classical_mechanics.rs               406       73        58      275          0             0.00
src/numerical/topology.rs                         403      102        40      261         30            11.49
src/ffi_apis/symbolic_pde_ffi/json.rs             403       81         8      314         26             8.28
~/ffi_apis/numerical_calculus_ffi/handle.rs       402       92        36      274         31            11.31
~rical_geometric_algebra_ffi/bincode_api.rs       402       63        65      274         10             3.65
src/ffi_apis/symbolic_vector_ffi/json.rs          402      115        32      255         27            10.59
~c/physics/physics_sim/linear_elasticity.rs       401       86        50      265         12             4.53
~is/symbolic_coordinates_ffi/bincode_api.rs       400      103        37      260         17             6.54
tests/numerical_multi_valued_ffi_test.rs          399       97         6      296          0             0.00
src/lib.rs                                        399       10       295       94          1             1.06
tests/symbolic_vector_calculus_test.rs            396       81        13      302          4             1.32
tests/numerical_physics_cfd_ffi_test.rs           396      122        11      263          7             2.66
~numerical_combinatorics_ffi/bincode_api.rs       395       72        73      250         10             4.00
src/symbolic/functional_analysis.rs               389       65        60      264          7             2.65
~fi_apis/symbolic_coordinates_ffi/handle.rs       387       40       153      194          8             4.12
src/numerical/ode.rs                              386       72        60      254          7             2.76
src/ffi_apis/numerical_graph_ffi/handle.rs        385       96        75      214         16             7.48
src/numerical/convergence.rs                      384       91       103      190         25            13.16
tests/numerical_physics_md_ffi_test.rs            380      105        11      264          2             0.76
~ffi_apis/symbolic_graph_ffi/bincode_api.rs       380       82        96      202         12             5.94
~physics/physics_sim/schrodinger_quantum.rs       375       75        57      243         13             5.35
~rc/ffi_apis/numerical_sparse_ffi/handle.rs       370       89        78      203         18             8.87
src/symbolic/optimize.rs                          369       63        66      240         19             7.92
benches/symbolic_polynomial.rs                    368       69         3      296          3             1.01
~/ffi_apis/symbolic_coordinates_ffi/json.rs       366       47       136      183         21            11.48
~/symbolic_integral_equations_ffi/handle.rs       365       79        10      276         38            13.77
tests/symbolic_cad_test.rs                        364       95        23      246          1             0.41
~is/numerical_vector_calculus_ffi/handle.rs       364       83        33      248         30            12.10
~fi_apis/numerical_graph_ffi/bincode_api.rs       363       77        49      237          9             3.80
tests/symbolic_stats_regression_test.rs           360       69        26      265         23             8.68
src/symbolic/numeric.rs                           359       40        31      288         18             6.25
src/ffi_apis/common.rs                            359       89        65      205         13             6.34
~c/ffi_apis/symbolic_pde_ffi/bincode_api.rs       359       70         7      282         22             7.80
tests/numerical_stats_test.rs                     358       99        20      239          0             0.00
src/symbolic/relativity.rs                        358       66        39      253          0             0.00
tests/symbolic_solid_state_physics_test.rs        357       99         2      256          0             0.00
~apis/numerical_combinatorics_ffi/handle.rs       352       88        73      191         12             6.28
~rical_computer_graphics_ffi/bincode_api.rs       352       67        49      236          6             2.54
~lic_graph_isomorphism_and_coloring_test.rs       350       77        11      262          4             1.53
src/numerical/number_theory.rs                    349      105        47      197         49            24.87
tests/symbolic_complex_analysis_test.rs           347       84        20      243          5             2.06
~i_apis/symbolic_multi_valued_ffi/handle.rs       347      112         9      226         19             8.41
src/symbolic/quantum_field_theory.rs              346       66        23      257          2             0.78
tests/symbolic_number_theory_test.rs              346       70        20      256          6             2.34
~bolic_computer_graphics_ffi/bincode_api.rs       346       73        16      257         17             6.61
src/symbolic/stats_inference.rs                   346       44        34      268          0             0.00
~ffi_apis/symbolic_elementary_ffi/handle.rs       345       83        85      177         10             5.65
~i_apis/numerical_interpolate_ffi/handle.rs       344       80        51      213         16             7.51
src/ffi_apis/numerical_tensor_ffi/json.rs         341       53        25      263         12             4.56
~fi_apis/numerical_multi_valued_ffi/json.rs       341       48        73      220          5             2.27
~symbolic/graph_isomorphism_and_coloring.rs       340       72        42      226         29            12.83
tests/numerical_graph_test.rs                     335      100        32      203         11             5.42
~mbolic_error_correction_helper_ffi/json.rs       335       67        68      200         17             8.50
src/symbolic/thermodynamics.rs                    333       57        37      239          1             0.42
tests/symbolic_stats_probability_test.rs          332       64        15      253         15             5.93
~lic_stats_information_theory_ffi/handle.rs       331       92        84      155         16            10.32
tests/symbolic_radicals_test.rs                   330       70         6      254         32            12.60
tests/symbolic_vector_test.rs                     328       71        28      229          9             3.93
~l_fractal_geometry_and_chaos_ffi/handle.rs       327       70        60      197         14             7.11
~rc/ffi_apis/numerical_tensor_ffi/handle.rs       327       79        57      191         16             8.38
src/symbolic/electromagnetism.rs                  324       56        45      223          0             0.00
src/output/typst.rs                               323       62         4      257         12             4.67
src/symbolic/vector_calculus.rs                   321       51        62      208          0             0.00
~is/symbolic_graph_operations_ffi/handle.rs       321      102         8      211         19             9.00
~rc/ffi_apis/symbolic_special_ffi/handle.rs       320       91        48      181          0             0.00
~_apis/numerical_multi_valued_ffi/handle.rs       319       75        65      179          5             2.79
~c/ffi_apis/symbolic_topology_ffi/handle.rs       319       93        19      207         19             9.18
tests/symbolic_stats_inference_test.rs            316       51        18      247         15             6.07
tests/symbolic_matrix_test.rs                     315       60         9      246         10             4.07
~mbolic_differential_geometry_ffi/handle.rs       313       87        15      211         20             9.48
tests/symbolic_functional_analysis_test.rs        310       56        23      231          3             1.30
~c/physics/physics_sim/ising_statistical.rs       310       76        35      199         18             9.05
~apis/numerical_vector_calculus_ffi/json.rs       307       36        97      174          6             3.45
src/ffi_apis/mod.rs                               306        3       139      164          0             0.00
~s/symbolic_finite_field_ffi/bincode_api.rs       305       74         9      222         20             9.01
~is/symbolic_electromagnetism_ffi/handle.rs       304       57        56      191         29            15.18
~fi_apis/numerical_elementary_ffi/handle.rs       300       90        35      175          9             5.14
~symbolic_classical_mechanics_ffi/handle.rs       300       60        64      176         21            11.93
~ctal_geometry_and_chaos_ffi/bincode_api.rs       299       73        21      205         16             7.80
src/ffi_apis/nightly_ffi/bincode_api.rs           299       46        33      220         16             7.27
~i_apis/numerical_matrix_ffi/bincode_api.rs       299       46        33      220         16             7.27
src/numerical/combinatorics.rs                    297       78        58      161         35            21.74
src/numerical/transforms.rs                       297       68       121      108         19            17.59
src/numerical/coordinates.rs                      294       68        75      151         16            10.60
~lic_fractal_geometry_and_chaos_ffi/json.rs       293       73        21      199         16             8.04
~ffi_apis/symbolic_multi_valued_ffi/json.rs       292       67         9      216         11             5.09
src/ffi_apis/symbolic_graph_ffi/json.rs           291       73        20      198         12             6.06
~s/symbolic_multi_valued_ffi/bincode_api.rs       291       66         9      216         11             5.09
src/numerical/solve.rs                            290       63        46      181         23            12.71
src/ffi_apis/symbolic_stats_ffi/handle.rs         288       71        56      161         21            13.04
benches/symbolic_elementary.rs                    287       61         0      226          0             0.00
~/numerical_multi_valued_ffi/bincode_api.rs       287       43        57      187          5             2.67
tests/symbolic_stats_test.rs                      285       55         5      225         15             6.67
~symbolic_solid_state_physics_ffi/handle.rs       284       59        65      160         20            12.50
~ffi_apis/symbolic_finite_field_ffi/json.rs       284       73         9      202         20             9.90
tests/symbolic_combinatorics_test.rs              284       78        25      181         18             9.94
~rc/ffi_apis/numerical_optimize_ffi/json.rs       284       35        42      207          7             3.38
src/ffi_apis/constant_ffi/handle.rs               282       79        40      163          7             4.29
~ymbolic_quantum_field_theory_ffi/handle.rs       282       56        64      162         22            13.58
~physics/physics_sim/geodesic_relativity.rs       280       53        52      175          3             1.71
~bolic_stats_information_theory_ffi/json.rs       279       82        77      120         17            14.17
~sts/symbolic_differential_geometry_test.rs       279       62        18      199         11             5.53
examples/quantum_tunneling_demo.rs                278       66        15      197         13             6.60
~pis/symbolic_stats_inference_ffi/handle.rs       278       70        70      138         13             9.42
~mbolic_graph_operations_ffi/bincode_api.rs       278       66        64      148          8             5.41
~error_correction_helper_ffi/bincode_api.rs       277       58        12      207         17             8.21
~merical_vector_calculus_ffi/bincode_api.rs       276       32        88      156          6             3.85
~c/physics/physics_sim/gpe_superfluidity.rs       274       55        24      195          3             1.54
tests/numerical_matrix_test.rs                    273       66         9      198          0             0.00
src/numerical/functional_analysis.rs              273       59        64      150          9             6.00
tests/numerical_solve_test.rs                     273       52        30      191         26            13.61
src/ffi_apis/compute_cache_ffi/handle.rs          272       65        36      171         20            11.70
src/ffi_apis/symbolic_handles_ffi/json.rs         272       70        27      175          9             5.14
~fi_apis/numerical_polynomial_ffi/handle.rs       271       72        73      126         12             9.52
tests/symbolic_transforms_test.rs                 270       67        20      183          3             1.64
examples/fdm_wave_3d_ripple.rs                    270       58         9      203         14             6.90
tests/symbolic_graph_operations_test.rs           269       68         7      194          0             0.00
src/ffi_apis/constant_ffi/json.rs                 269       71        40      158          4             2.53
src/ffi_apis/numerical_sparse_ffi/json.rs         266       45        48      173          8             4.62
tests/ffi_constant_three_version_test.rs          265       86         4      175          0             0.00
tests/symbolic_rewriting_test.rs                  265       55        18      192          2             1.04
~ests/numerical_functional_analysis_test.rs       264       65        20      179          1             0.56
~fi_apis/numerical_finite_field_ffi/json.rs       264       52        37      175          9             5.14
~is/symbolic_stats_regression_ffi/handle.rs       263       63        44      156         16            10.26
~umerical_differential_geometry_ffi/json.rs       262       49        49      164          8             4.88
tests/numerical_combinatorics_test.rs             262       61        24      177          6             3.39
tests/symbolic_integral_equations_test.rs         262       58        43      161          3             1.86
~ts/symbolic_calculus_of_variations_test.rs       261       61        17      183          0             0.00
~c/ffi_apis/symbolic_elementary_ffi/json.rs       258       53        78      127          9             7.09
tests/symbolic/calculus.rs                        258       54        25      179          2             1.12
~pis/symbolic_polynomial_ffi/bincode_api.rs       255       58         7      190         18             9.47
tests/numerical_sparse_test.rs                    254       54         1      199          1             0.50
src/numerical/multi_valued.rs                     253       68        27      158         10             6.33
src/ffi_apis/numerical_solve_ffi/handle.rs        253       54        71      128         12             9.38
tests/physics_rkm_test.rs                         253       67         9      177          5             2.82
~c/ffi_apis/symbolic_polynomial_ffi/json.rs       252       58         7      187         18             9.63
~hysics/physics_sim/fdtd_electrodynamics.rs       252       67        26      159         16            10.06
~ffi_apis/numerical_interpolate_ffi/json.rs       251       50        49      152          7             4.61
~/symbolic_quantum_field_theory_ffi/json.rs       250       53        48      149          8             5.37
~i_apis/numerical_coordinates_ffi/handle.rs       249       62        44      143         11             7.69
tests/numerical_geometric_algebra_test.rs         248       70        13      165          0             0.00
src/symbolic/stats_regression.rs                  248       46        52      150          5             3.33
~rc/ffi_apis/numerical_calculus_ffi/json.rs       248       39        37      172          6             3.49
~fi_apis/symbolic_combinatorics_ffi/json.rs       247       37       124       86          4             4.65
tests/symbolic_dag_migration_test.rs              244       80        17      147          0             0.00
~c_differential_geometry_ffi/bincode_api.rs       243       53         7      183          8             4.37
~erical_differential_geometry_ffi/handle.rs       242       58        34      150         13             8.67
src/ffi_apis/physics_mm_ffi/handle.rs             240       53        43      144         11             7.64
~c/ffi_apis/symbolic_optimize_ffi/handle.rs       240       55        12      173         15             8.67
tests/symbolic_cas_foundations_test.rs            237       50        13      174          4             2.30
~pis/symbolic_cas_foundations_ffi/handle.rs       235       48        56      131         15            11.45
tests/numerical_complex_analysis_test.rs          234       54         6      174          2             1.15
tests/numerical_integrate_ffi_test.rs             234       54         5      175          2             1.14
~symbolic_differential_geometry_ffi/json.rs       233       53         7      173          8             4.62
~s/numerical_complex_analysis_ffi/handle.rs       232       55        25      152         17            11.18
tests/symbolic_ode_test.rs                        232       63        23      146          0             0.00
tests/numerical_calculus_test.rs                  232       38        14      180          0             0.00
~i_apis/symbolic_handles_ffi/bincode_api.rs       232       61         8      163          7             4.29
src/ffi_apis/numerical_signal_ffi/json.rs         232       31        93      108          3             2.78
tests/numerical_optimize_test.rs                  231       47         2      182          3             1.65
tests/physics_bem_test.rs                         231       49        10      172         12             6.98
src/ffi_apis/constant_ffi/bincode_api.rs          231       65        28      138          4             2.90
~olic_integral_equations_ffi/bincode_api.rs       231       48         6      177          7             3.95
tests/symbolic_operators_test.rs                  231       57        21      153          0             0.00
~tats_information_theory_ffi/bincode_api.rs       230       74        28      128         17            13.28
examples/rkm_adaptive_solver_comparison.rs        229       45         8      176          1             0.57
~l_differential_geometry_ffi/bincode_api.rs       229       44        33      152          8             5.26
tests/symbolic_real_roots_test.rs                 229       55        17      157          5             3.18
tests/physics_fdm_test.rs                         228       52        12      164         15             9.15
tests/symbolic_coordinates_test.rs                227       56        14      157          4             2.55
~s/numerical_physics_fea_ffi/bincode_api.rs       226       27        91      108          3             2.78
~s/symbolic_geometric_algebra_ffi/handle.rs       226       73         8      145         11             7.59
~s/numerical_physics_cfd_ffi/bincode_api.rs       225       26        94      105          3             2.86
~_apis/symbolic_number_theory_ffi/handle.rs       225       53        17      155          9             5.81
~apis/numerical_optimize_ffi/bincode_api.rs       225       22        40      163          4             2.45
tests/numerical_interpolate_ffi_test.rs           225       57         5      163          0             0.00
~rc/ffi_apis/numerical_signal_ffi/handle.rs       224       52        41      131          7             5.34
~_apis/numerical_finite_field_ffi/handle.rs       224       56        46      122          9             7.38
examples/fdm_wave_2d_ripple.rs                    223       52        10      161         12             7.45
~fi_apis/symbolic_matrix_ffi/bincode_api.rs       223       65        24      134         10             7.46
~ests/numerical_vector_calculus_ffi_test.rs       223       52         2      169          0             0.00
~apis/symbolic_graph_operations_ffi/json.rs       221       59        13      149          8             5.37
tests/numerical_vector_test.rs                    221       64         1      156          8             5.13
tests/numerical_integrate_test.rs                 221       51        16      154          0             0.00
examples/ode_lorenz_attractor.rs                  221       38        18      165          6             3.64
~s/symbolic_functional_analysis_ffi/json.rs       221       31       106       84          8             9.52
examples/ode_lorenz_attractor_1.rs                221       38        18      165          6             3.64
examples/fdm_wave_2d_surface.rs                   221       54         4      163         12             7.36
~i_apis/numerical_convergence_ffi/handle.rs       220       51        71       98          7             7.14
~bolic_calculus_of_variations_ffi/handle.rs       219       38        41      140         17            12.14
~rical_calculus_of_variations_ffi/handle.rs       219       39        18      162         17            10.49
~/symbolic_combinatorics_ffi/bincode_api.rs       219       31       122       66          4             6.06
src/symbolic/core/mod.rs                          217        4       201       12          0             0.00
~ffi_apis/numerical_coordinates_ffi/json.rs       217       39        25      153          8             5.23
src/symbolic/stats.rs                             217       42        59      116          7             6.03
src/ffi_apis/symbolic_matrix_ffi/json.rs          217       66        24      127         10             7.87
examples/relativity_black_hole_orbits.rs          216       27         5      184          1             0.54
tests/symbolic_graph_test.rs                      215       56         5      154          0             0.00
tests/symbolic_group_theory_test.rs               215       51        13      151          4             2.65
~merical_calculus_of_variations_ffi_test.rs       214       46         0      168          0             0.00
~apis/numerical_calculus_ffi/bincode_api.rs       214       35        25      154          6             3.90
~_apis/symbolic_topology_ffi/bincode_api.rs       213       49         8      156         13             8.33
src/compute/cache.rs                              212       39        54      119          2             1.68
~ffi_apis/numerical_convergence_ffi/json.rs       212       32        37      143          6             4.20
~ymbolic_fractal_geometry_and_chaos_test.rs       210       49        29      132          1             0.76
~i_apis/numerical_number_theory_ffi/json.rs       210       39        25      146          7             4.79
tests/symbolic_integration_test.rs                210       50        28      132          0             0.00
~s/numerical_interpolate_ffi/bincode_api.rs       210       45        33      132          7             5.30
~_apis/symbolic_combinatorics_ffi/handle.rs       210       31       122       57          0             0.00
examples/pde_heat_equation.rs                     210       42        18      150          9             6.00
src/ffi_apis/physics_fvm_ffi/json.rs              210       25        71      114          5             4.39
src/numerical/signal.rs                           210       40        90       80         14            17.50
tests/symbolic_grobner_test.rs                    209       55         4      150          2             1.33
tests/numerical_transforms_ffi_test.rs            208       47         2      159          0             0.00
src/symbolic/calculus_of_variations.rs            204       22        97       85          0             0.00
~lic_functional_analysis_ffi/bincode_api.rs       204       30       106       68          8            11.76
~sts/numerical_complex_analysis_ffi_test.rs       203       52         0      151          3             1.99
tests/numerical_signal_ffi_test.rs                201       48         2      151          0             0.00
~/numerical_functional_analysis_ffi/json.rs       201       32        37      132          5             3.79
src/ffi_apis/plugins_ffi/handle.rs                201       34        43      124         13            10.48
examples/fvm_advection_diffusion.rs               200       44        21      135         15            11.11
examples/fluid_dynamics_karman_vortex.rs          199       43        10      146         12             8.22
~ffi_apis/symbolic_polynomial_ffi/handle.rs       198       53        27      118          1             0.85
~_calculus_of_variations_ffi/bincode_api.rs       198       44         4      150         15            10.00
~i_apis/numerical_signal_ffi/bincode_api.rs       198       27        81       90          3             3.33
~ymbolic_calculus_of_variations_ffi/json.rs       197       44         4      149         15            10.07
tests/numerical_calculus_ffi_test.rs              197       48         0      149          0             0.00
src/ffi_apis/numerical_series_ffi/json.rs         197       23        66      108          4             3.70
tests/numerical_topology_test.rs                  196       43        19      134          7             5.22
src/ffi_apis/physics_mtm_ffi/json.rs              196       19        65      112          4             3.57
~umerical_functional_analysis_ffi/handle.rs       196       46        41      109         12            11.01
src/ffi_apis/macros.rs                            195       54        13      128          8             6.25
~i_apis/numerical_physics_fea_ffi/handle.rs       194       46        37      111          3             2.70
tests/numerical_ode_ffi_test.rs                   193       45         1      147          0             0.00
src/numerical/calculus_of_variations.rs           193       20        92       81          0             0.00
~is/symbolic_integral_equations_ffi/json.rs       193       44         6      143          7             4.90
src/ffi_apis/symbolic_topology_ffi/json.rs        193       49         8      136         13             9.56
~/symbolic_number_theory_ffi/bincode_api.rs       192       48        17      127          8             6.30
~ffi_apis/symbolic_real_roots_ffi/handle.rs       191       44        11      136         13             9.56
~ymbolic_vector_calculus_ffi/bincode_api.rs       189       42         5      142          4             2.82
~rc/ffi_apis/symbolic_handles_ffi/handle.rs       188       36        65       87          4             4.60
tests/numerical_series_ffi_test.rs                187       47         1      139          0             0.00
benches/compute_engine.rs                         186       41         0      145          0             0.00
~fi_apis/symbolic_number_theory_ffi/json.rs       186       49        17      120          8             6.67
~is/numerical_physics_md_ffi/bincode_api.rs       185       23        61      101          3             2.97
src/ffi_apis/physics_sm_ffi/json.rs               185       21        73       91          2             2.20
tests/symbolic_series_test.rs                     184       42        16      126          0             0.00
examples/symbolic_operators_demo.rs               184       58         8      118          0             0.00
~apis/symbolic_thermodynamics_ffi/handle.rs       184       36        41      107         16            14.95
~rc/ffi_apis/numerical_topology_ffi/json.rs       183       24        66       93          2             2.15
src/ffi_apis/symbolic_matrix_ffi/handle.rs        183       65        24       94          1             1.06
src/plugins/plugin_c.rs                           181       33        53       95          2             2.11
src/ffi_apis/physics_em_ffi/json.rs               181       18        36      127          5             3.94
~rc/ffi_apis/numerical_series_ffi/handle.rs       180       37        26      117         10             8.55
tests/numerical/interpolate.rs                    179       32        13      134          7             5.22
~pis/symbolic_geometric_algebra_ffi/json.rs       179       50         7      122         11             9.02
tests/numerical/calculus.rs                       178       40        11      127          3             2.36
tests/symbolic_matrix_symbolic_test.rs            178       45        25      108         13            12.04
tests/symbolic_unit_unification_test.rs           178       44         5      129         15            11.63
~bolic_geometric_algebra_ffi/bincode_api.rs       177       49         7      121         11             9.09
~pis/symbolic_elementary_ffi/bincode_api.rs       177       40        21      116          9             7.76
~i_apis/numerical_series_ffi/bincode_api.rs       176       20        60       96          4             4.17
tests/symbolic_optimize_test.rs                   176       43         8      125          2             1.60
~s/numerical_convergence_ffi/bincode_api.rs       176       30        25      121          4             3.31
~/ffi_apis/numerical_polynomial_ffi/json.rs       175       35        25      115          6             5.22
~/symbolic_poly_factorization_ffi/handle.rs       171       43         6      122         10             8.20
tests/symbolic/simplify.rs                        170       45        25      100          0             0.00
~pis/symbolic_real_roots_ffi/bincode_api.rs       170       39         3      128         14            10.94
src/ffi_apis/physics_em_ffi/bincode_api.rs        169       17        33      119          5             4.20
tests/numerical_elementary_test.rs                169       35         5      129          1             0.78
src/ffi_apis/symbolic_tensor_ffi/json.rs          169       40        16      113          9             7.96
tests/symbolic_electromagnetism_test.rs           168       55         2      111          0             0.00
~/ffi_apis/compute_cache_ffi/bincode_api.rs       168       36         7      125         12             9.60
~c/ffi_apis/symbolic_real_roots_ffi/json.rs       168       39         3      126         14            11.11
tests/symbolic_dag_serialization_test.rs          168       42        20      106          0             0.00
tests/numerical_convergence_test.rs               166       46        20      100          8             8.00
~s/numerical_coordinates_ffi/bincode_api.rs       165       22        17      126          6             4.76
tests/numerical_topology_ffi_test.rs              165       44         1      120          0             0.00
~apis/symbolic_stats_regression_ffi/json.rs       165       41        38       86          7             8.14
~umerical_differential_geometry_ffi_test.rs       165       42         0      123          0             0.00
src/ffi_apis/compute_cache_ffi/json.rs            164       36         9      119         18            15.13
~fi_apis/numerical_transforms_ffi/handle.rs       163       37        25      101          8             7.92
~apis/symbolic_rewriting_ffi/bincode_api.rs       162       40         5      117          5             4.27
~apis/numerical_topology_ffi/bincode_api.rs       162       21        60       81          2             2.47
~i_apis/symbolic_finite_field_ffi/handle.rs       162       50         6      106          9             8.49
~cal_functional_analysis_ffi/bincode_api.rs       161       28        25      108          5             4.63
tests/symbolic_topology_test.rs                   161       43        14      104          6             5.77
tests/numerical_vector_calculus_test.rs           160       35         7      118          0             0.00
~_apis/symbolic_stats_inference_ffi/json.rs       158       40        34       84          3             3.57
tests/physics_em_test.rs                          158       35         6      117          4             3.42
~_apis/symbolic_vector_calculus_ffi/json.rs       158       39         5      114          4             3.51
tests/numerical/polynomial.rs                     157       27        34       96          4             4.17
~rc/ffi_apis/symbolic_rewriting_ffi/json.rs       157       37        17      103          5             4.85
src/compute/mod.rs                                156        1       150        5          0             0.00
~pis/numerical_complex_analysis_ffi/json.rs       155       21        25      109          4             3.67
src/ffi_apis/plugins_ffi/json.rs                  153       24        13      116          5             4.31
~i_apis/numerical_tensor_ffi/bincode_api.rs       151       21         9      121          6             4.96
src/numerical/series.rs                           150       35        50       65          3             4.62
~fi_apis/numerical_real_roots_ffi/handle.rs       150       32        41       77          8            10.39
~s/numerical_calculus_of_variations_test.rs       150       26        19      105          0             0.00
src/ffi_apis/physics_mtm_ffi/handle.rs            150       38        29       83          5             6.02
~/ffi_apis/numerical_transforms_ffi/json.rs       149       21        61       67          2             2.99
~ffi_apis/symbolic_relativity_ffi/handle.rs       148       32        36       80          6             7.50
examples/numerical_vector_demo.rs                 147       36         7      104          1             0.96
tests/symbolic_classical_mechanics_test.rs        147       60         6       81          0             0.00
~i_apis/numerical_physics_cfd_ffi/handle.rs       147       33        17       97          0             0.00
src/ffi_apis/physics_rkm_ffi/handle.rs            147       48         5       94         12            12.77
tests/symbolic_simplify_dag_test.rs               146       40        12       94          0             0.00
src/ffi_apis/symbolic_tensor_ffi/handle.rs        145       42        16       87          4             4.60
tests/numerical_multi_valued_test.rs              145       38         6      101          1             0.99
examples/grobner_basis_demo_2.rs                  145       30        15      100          2             2.00
~i_apis/physics_sim_schrodinger_ffi/json.rs       143       15        38       90          4             4.44
examples/grobner_basis_demo_1.rs                  143       23        13      107          4             3.74
tests/numerical_interpolate_test.rs               142       32         2      108          1             0.93
src/ffi_apis/numerical_ode_ffi/handle.rs          142       31         9      102          9             8.82
~/ffi_apis/numerical_topology_ffi/handle.rs       141       33        25       83          6             7.23
tests/numerical_finite_field_test.rs              141       39        12       90          2             2.22
src/ffi_apis/symbolic_proof_ffi/json.rs           140       25        24       91          3             3.30
~mbolic_stats_regression_ffi/bincode_api.rs       139       37        14       88          7             7.95
src/ffi_apis/symbolic_cad_ffi/handle.rs           139       38         6       95         13            13.68
tests/symbolic_convergence_test.rs                138       26         7      105          0             0.00
~_apis/symbolic_cas_foundations_ffi/json.rs       138       34         5       99          4             4.04
~_apis/symbolic_optimize_ffi/bincode_api.rs       137       23         3      111          5             4.50
~ts/numerical_geometric_algebra_ffi_test.rs       137       33         6       98          0             0.00
tests/symbolic_quantum_mechanics_test.rs          137       36         1      100          0             0.00
src/ffi_apis/symbolic_optimize_ffi/json.rs        137       24         3      110          5             4.55
examples/grobner_basis_demo_3.rs                  137       22        14      101          2             1.98
tests/symbolic_proof_test.rs                      136       28         5      103          0             0.00
tests/numerical_coordinates_test.rs               136       31        13       92          0             0.00
src/ffi_apis/physics_mm_ffi/json.rs               136       18        57       61          1             1.64
~fi_apis/symbolic_tensor_ffi/bincode_api.rs       136       29        12       95          6             6.32
~ymbolic_cas_foundations_ffi/bincode_api.rs       136       33         5       98          4             4.08
~ests/symbolic_quantum_field_theory_test.rs       134       43         0       91          0             0.00
src/ffi_apis/physics_fdm_ffi/handle.rs            134       35        28       71          7             9.86
~/numerical_finite_field_ffi/bincode_api.rs       134       25        17       92          4             4.35
src/ffi_apis/plugins_ffi/bincode_api.rs           133       18         9      106          5             4.72
tests/numerical_optimize_ffi_test.rs              132       33         6       93          0             0.00
~ymbolic_stats_inference_ffi/bincode_api.rs       132       35        12       85          3             3.53
~/numerical_functional_analysis_ffi_test.rs       132       33         0       99          0             0.00
~ffi_apis/symbolic_solve_ffi/bincode_api.rs       132       31        12       89          5             5.62
src/ffi_apis/physics_bem_ffi/json.rs              131       17        33       81          3             3.70
~ffi_apis/numerical_integrate_ffi/handle.rs       131       22        27       82          6             7.32
src/ffi_apis/jit_ffi/json.rs                      130       22        14       94          5             5.32
~i_apis/symbolic_grobner_ffi/bincode_api.rs       130       15        38       77          4             5.19
examples/finite_field_debug.rs                    130       14         6      110          0             0.00
src/ffi_apis/symbolic_solve_ffi/json.rs           129       32        12       85          5             5.88
src/ffi_apis/symbolic_stats_ffi/json.rs           128       37         6       85          8             9.41
~physics_sim_schrodinger_ffi/bincode_api.rs       128       13        36       79          4             5.06
~is/numerical_transforms_ffi/bincode_api.rs       128       18        55       55          2             3.64
src/ffi_apis/symbolic_solve_ffi/handle.rs         128       37        12       79          3             3.80
~c/ffi_apis/numerical_integrate_ffi/json.rs       127       19        23       85          4             4.71
~rc/ffi_apis/physics_rkm_ffi/bincode_api.rs       127       19        34       74          2             2.70
tests/symbolic_thermodynamics_test.rs             127       51         5       71          0             0.00
~erical_complex_analysis_ffi/bincode_api.rs       127       19        17       91          4             4.40
tests/numerical_ode_test.rs                       126       26         6       94          0             0.00
~fi_apis/symbolic_vector_ffi/bincode_api.rs       126       44        16       66          6             9.09
~ests/physics_sim_linear_elasticity_test.rs       125       34         3       88          4             4.55
tests/numerical_number_theory_test.rs             125       28         0       97          1             1.03
~is/numerical_real_roots_ffi/bincode_api.rs       124       18        28       78          3             3.85
build.rs                                          123       26         3       94          2             2.13
benches/symbolic_matrix.rs                        123       16         3      104          0             0.00
~i_apis/numerical_sparse_ffi/bincode_api.rs       122       18         9       95          4             4.21
~ffi_apis/symbolic_stats_ffi/bincode_api.rs       122       29         6       87          8             9.20
src/ffi_apis/symbolic_grobner_ffi/json.rs         121       16        35       70          4             5.71
tests/numerical_transforms_test.rs                121       28         1       92          2             2.17
tests/numerical_tensor_test.rs                    121       31         0       90          0             0.00
tests/numerical_real_roots_test.rs                120       27        10       83          8             9.64
tests/compute_computation_test.rs                 120       17         8       95          0             0.00
~apis/physics_sim_schrodinger_ffi/handle.rs       120       21         9       90          6             6.67
~/ffi_apis/numerical_real_roots_ffi/json.rs       120       16        31       73          3             4.11
~lic_solid_state_physics_ffi/bincode_api.rs       120       25         4       91          3             3.30
tests/numerical_polynomial_test.rs                119       26         0       93          0             0.00
examples/matrix_demo.rs                           119       32         9       78          1             1.28
src/ffi_apis/jit_ffi/handle.rs                    119       26        41       52          4             7.69
~rc/ffi_apis/physics_bem_ffi/bincode_api.rs       119       15        32       72          3             4.17
src/ffi_apis/physics_em_ffi/handle.rs             118       37         4       77         11            14.29
~is/numerical_elementary_ffi/bincode_api.rs       118       22         9       87          4             4.60
~is/symbolic_poly_factorization_ffi/json.rs       117       27         4       86         11            12.79
~olic_poly_factorization_ffi/bincode_api.rs       117       26         4       87         11            12.64
benches/symbolic_simplify_dag.rs                  117       24         4       89          0             0.00
~s/symbolic_solid_state_physics_ffi/json.rs       117       26         4       87          3             3.45
~is/physics_sim_navier_stokes_ffi/handle.rs       117       21        17       79          4             5.06
src/ffi_apis/symbolic_logic_ffi/handle.rs         116       31        21       64          5             7.81
examples/solid_state_physics_demo.rs              116       25         4       87          0             0.00
tests/physics_fvm_test.rs                         115       33         3       79          9            11.39
src/ffi_apis/physics_sim_fdtd_ffi/json.rs         115       14        36       65          4             6.15
tests/symbolic_relativity_test.rs                 114       45         0       69          0             0.00
src/jit/instructions.rs                           114       14        45       55          0             0.00
~apis/physics_sim_navier_stokes_ffi/json.rs       114       13        37       64          2             3.12
src/ffi_apis/compute_state_ffi/json.rs            114       26         8       80         10            12.50
examples/simplify_with_relations_demo.rs          113       19        17       77          2             2.60
~_apis/symbolic_discrete_groups_ffi/json.rs       112       20        53       39          1             2.56
~pis/symbolic_discrete_groups_ffi/handle.rs       112       19        53       40          1             2.50
~ffi_apis/symbolic_proof_ffi/bincode_api.rs       112       19        16       77          2             2.60
~lic_classical_mechanics_ffi/bincode_api.rs       112       27         3       82          8             9.76
examples/statistical_ising_model.rs               111       26         3       82          3             3.66
~ymbolic_discrete_groups_ffi/bincode_api.rs       111       19        53       39          1             2.56
tests/numerical_signal_test.rs                    110       29         7       74          0             0.00
~/ffi_apis/numerical_elementary_ffi/json.rs       110       20        15       75          4             5.33
tests/physics_sim_gpe_test.rs                     110       24         6       80          1             1.25
~s/symbolic_classical_mechanics_ffi/json.rs       109       27         3       79          8            10.13
src/ffi_apis/numerical_ode_ffi/json.rs            109       13        34       62          2             3.23
src/ffi_apis/symbolic_vector_ffi/handle.rs        109       44        16       49          0             0.00
tests/physics_mm_test.rs                          108       16         3       89          1             1.12
~rc/ffi_apis/physics_fdm_ffi/bincode_api.rs       107       15        34       58          1             1.72
src/ffi_apis/symbolic_logic_ffi/json.rs           107       29         9       69          9            13.04
tests/physics_sim_schrodinger_test.rs             107       21         4       82          2             2.44
src/ffi_apis/physics_bem_ffi/handle.rs            107       20        12       75         10            13.33
tests/numerical_solve_ffi_test.rs                 107       30         8       69          0             0.00
tests/physics_rkm_ffi_test.rs                     106       29         3       74          0             0.00
~fi_apis/numerical_solve_ffi/bincode_api.rs       105       17         9       79          3             3.80
src/ffi_apis/physics_fem_ffi/json.rs              105       11        32       62          2             3.23
~apis/numerical_number_theory_ffi/handle.rs       105       28        17       60          1             1.67
tests/physics_cnm_test.rs                         104       26         6       72          4             5.56
~rc/ffi_apis/symbolic_grobner_ffi/handle.rs       103       14        34       55          2             3.64
tests/numerical_convergence_ffi_test.rs           102       24         5       73          0             0.00
~ysics_sim_navier_stokes_ffi/bincode_api.rs       102       11        36       55          2             3.64
tests/physics_sim_schrodinger_ffi_test.rs         101       33         1       67          0             0.00
src/ffi_apis/numerical_solve_ffi/json.rs          100       15        13       72          3             4.17
examples/latex_output_demo.rs                     100       18         9       73          1             1.37
src/compute/computation.rs                        100       14        20       66          0             0.00
tests/physics_fdm_ffi_test.rs                     100       30         3       67          0             0.00
tests/numerical_real_roots_ffi_test.rs             99       28         4       67          0             0.00
~ts/numerical_differential_geometry_test.rs        99       17        18       64          0             0.00
~graph_isomorphism_and_coloring_ffi/json.rs        97       20        30       47          3             6.38
~ffi_apis/symbolic_logic_ffi/bincode_api.rs        97       28         9       60          8            13.33
~/ffi_apis/numerical_ode_ffi/bincode_api.rs        96       11        29       56          2             3.57
tests/plugins/mod.rs                               96       21         0       75          6             8.00
src/ffi_apis/compute_state_ffi/handle.rs           95       26         7       62          6             9.68
src/ffi_apis/physics_sim_gpe_ffi/json.rs           95        8        36       51          2             3.92
src/ffi_apis/physics_sim_ising_ffi/json.rs         94       13        34       47          1             2.13
~is/numerical_polynomial_ffi/bincode_api.rs        94       18         9       67          3             4.48
~aph_isomorphism_and_coloring_ffi/handle.rs        94       26         8       60          4             6.67
~/physics_sim_linear_elasticity_ffi/json.rs        94        8        35       51          2             3.92
src/ffi_apis/physics_fvm_ffi/handle.rs             94       25        19       50          2             4.00
~fi_apis/numerical_physics_md_ffi/handle.rs        93       29        16       48          4             8.33
~numerical_number_theory_ffi/bincode_api.rs        93       18         9       66          3             4.55
examples/repro_test_origin.rs                      93       24         3       66          7            10.61
examples/repro_test.rs                             93       24         3       66          7            10.61
~rc/ffi_apis/physics_fem_ffi/bincode_api.rs        91        9        29       53          2             3.77
src/ffi_apis/physics_cnm_ffi/json.rs               91       12        33       46          1             2.17
~hysics_sim_linear_elasticity_ffi/handle.rs        91       25         4       62          5             8.06
~somorphism_and_coloring_ffi/bincode_api.rs        91       20        24       47          3             6.38
benches/compute_computation.rs                     90       11         0       79          1             1.27
src/symbolic/mod.rs                                90        1        21       68          0             0.00
tests/compute_cache_test.rs                        89       27         2       60          2             3.33
~rc/ffi_apis/physics_mtm_ffi/bincode_api.rs        89        9        30       50          2             4.00
~rc/ffi_apis/physics_fvm_ffi/bincode_api.rs        89       10        33       46          1             2.17
tests/symbolic_discrete_groups_test.rs             88       29         5       54          6            11.11
tests/numerical_series_test.rs                     87       18         1       68          0             0.00
~merical_calculus_of_variations_ffi/json.rs        87       16        13       58          2             3.45
~c/ffi_apis/physics_sim_ising_ffi/handle.rs        87       16        16       55          1             1.82
src/ffi_apis/physics_cnm_ffi/handle.rs             87       21        18       48          2             4.17
tests/physics_sim_ising_test.rs                    86       21         4       61          5             8.20
tests/physics_mtm_test.rs                          86       25         3       58          4             6.90
examples/logic_debug.rs                            86       14         3       69          0             0.00
~i_apis/physics_sim_fdtd_ffi/bincode_api.rs        86       10        33       43          3             6.98
tests/physics_sim_geodesic_test.rs                 85       20         7       58          1             1.72
~mbolic_electromagnetism_ffi/bincode_api.rs        85       17         3       65          2             3.08
~pis/symbolic_quantum_mechanics_ffi/json.rs        85       22         3       60          3             5.00
~fi_apis/symbolic_integration_ffi/handle.rs        84       24         2       58          5             8.62
~bolic_quantum_mechanics_ffi/bincode_api.rs        83       21         3       59          3             5.08
~_apis/physics_sim_ising_ffi/bincode_api.rs        83       11        33       39          1             2.56
tests/regression_test.rs                           83       21         2       60          3             5.00
~fi_apis/physics_sim_gpe_ffi/bincode_api.rs        82        6        34       42          2             4.76
src/ffi_apis/physics_sm_ffi/bincode_api.rs         82       10        35       37          1             2.70
src/plugins/stable_abi.rs                          82       18        13       51          0             0.00
~/ffi_apis/physics_sim_geodesic_ffi/json.rs        81        9        35       37          1             2.70
~pis/numerical_integrate_ffi/bincode_api.rs        81       13         9       59          2             3.39
~s_sim_linear_elasticity_ffi/bincode_api.rs        81        6        33       42          2             4.76
tests/physics_mm_ffi_test.rs                       80       24         1       55          0             0.00
tests/physics_bem_ffi_test.rs                      79       24         1       54          0             0.00
~rc/ffi_apis/physics_cnm_ffi/bincode_api.rs        78       10        30       38          1             2.63
tests/debug_limit_test.rs                          78       16         0       62          1             1.61
~c/ffi_apis/symbolic_simplify_ffi/handle.rs        78       17        23       38          2             5.26
~_calculus_of_variations_ffi/bincode_api.rs        78       14         9       55          2             3.64
~i_apis/symbolic_thermodynamics_ffi/json.rs        78       20         3       55          2             3.64
tests/physics_sim_fdtd_test.rs                     77       16         4       57          2             3.51
~symbolic_thermodynamics_ffi/bincode_api.rs        77       19         3       55          2             3.64
tests/physics_cnm_ffi_test.rs                      76       21         1       54          0             0.00
src/ffi_apis/physics_mm_ffi/bincode_api.rs         75       10        29       36          1             2.78
benches/compute_cache.rs                           75       18         0       57          0             0.00
tests/physics_sim_navier_stokes_test.rs            75       18         4       53          3             5.66
tests/constant_test.rs                             75       25         0       50          0             0.00
~/ffi_apis/compute_state_ffi/bincode_api.rs        74       16         8       50          2             4.00
tests/physics_sim_ising_ffi_test.rs                74       20         1       53          1             1.89
~apis/symbolic_electromagnetism_ffi/json.rs        74       18         3       53          2             3.77
src/ffi_apis/physics_fem_ffi/handle.rs             73       20        11       42          2             4.76
tests/symbolic/unit_unification.rs                 71       12         1       58          5             8.62
tests/symbolic/extra_simplify_tests.rs             70       26         0       44          0             0.00
src/ffi_apis/symbolic_cad_ffi/json.rs              70       16        10       44          4             9.09
~pis/symbolic_relativity_ffi/bincode_api.rs        70       19         4       47          6            12.77
~ic_quantum_field_theory_ffi/bincode_api.rs        69       16        16       37          2             5.41
~c/ffi_apis/symbolic_relativity_ffi/json.rs        69       20         4       45          6            13.33
examples/debug_proof.rs                            69       12         0       57          3             5.26
~ests/physics_sim_navier_stokes_ffi_test.rs        68       21         1       46          0             0.00
~is/physics_sim_geodesic_ffi/bincode_api.rs        68        7        32       29          1             3.45
tests/physics_mtm_ffi_test.rs                      67       20         1       46          0             0.00
tests/dag_consistency_test.rs                      67       14         3       50          2             4.00
tests/physics_sm_ffi_test.rs                       66       20         1       45          0             0.00
tests/physics_sim_gpe_ffi_test.rs                  64       20         1       43          0             0.00
examples/debug_singularity.rs                      63       14         1       48          3             6.25
tests/physics_fvm_ffi_test.rs                      62       20         1       41          0             0.00
~/ffi_apis/symbolic_integration_ffi/json.rs        62       15         2       45          3             6.67
tests/physics_sim_geodesic_ffi_test.rs             62       20         1       41          0             0.00
~c/ffi_apis/symbolic_cad_ffi/bincode_api.rs        61       14         3       44          4             9.09
tests/physics_em_ffi_test.rs                       61       18         0       43          0             0.00
~apis/symbolic_special_functions_ffi/mod.rs        61        2        56        3          0             0.00
~is/symbolic_integration_ffi/bincode_api.rs        61       14         2       45          3             6.67
tests/physics_sim_fdtd_ffi_test.rs                 61       20         1       40          0             0.00
examples/lib_examples.rs                           60       12         6       42          0             0.00
tests/physics_fem_ffi_test.rs                      59       19         1       39          0             0.00
benches/numerical_vector.rs                        59       13         0       46          0             0.00
~/physics_sim_linear_elasticity_ffi_test.rs        59       20         2       37          0             0.00
benches/compute_computable.rs                      58       12         0       46          1             2.17
benches/prelude.rs                                 57       12        21       24          0             0.00
benches/lib.rs                                     57       12        22       23          0             0.00
src/numerical/mod.rs                               56        1         8       47          0             0.00
tests/physics_sm_test.rs                           55       16         3       36          1             2.78
~rc/ffi_apis/physics_sim_fdtd_ffi/handle.rs        55       11         2       42          2             4.76
tests/compute_computable_test.rs                   54       14         0       40          1             2.50
src/ffi_apis/physics_sm_ffi/handle.rs              53       14         5       34          3             8.82
~ffi_apis/symbolic_complex_analysis_test.rs        53       17         0       36          0             0.00
src/ffi_apis/symbolic_simplify_ffi/json.rs         53       16         3       34          4            11.76
src/physics/mod.rs                                 52        2        38       12          0             0.00
~_apis/symbolic_simplify_ffi/bincode_api.rs        52       15         3       34          4            11.76
~fi_apis/physics_sim_geodesic_ffi/handle.rs        51       11         2       38          1             2.63
tests/mod.rs                                       51        6        42        3          0             0.00
src/ffi_apis/physics_sim_gpe_ffi/handle.rs         51        8         2       41          1             2.44
~c/ffi_apis/symbolic_radicals_ffi/handle.rs        50       16         2       32          2             6.25
src/ffi_apis/symbolic_radicals_ffi/json.rs         50       15         2       33          4            12.12
tests/physics_fem_test.rs                          50       13         2       35          1             2.86
tests/ffi_apis/mod.rs                              49        6        42        1          0             0.00
~_apis/symbolic_radicals_ffi/bincode_api.rs        49       14         2       33          4            12.12
src/ffi_apis/symbolic_special_ffi/mod.rs           49        1        45        3          0             0.00
tests/numerical/geometric_algebra.rs               48        6        25       17          0             0.00
tests/symbolic/cryptography.rs                     48        6        25       17          0             0.00
tests/symbolic/rewriting.rs                        48        6        25       17          0             0.00
tests/symbolic/series.rs                           48        6        25       17          0             0.00
tests/symbolic/numeric.rs                          48        6        25       17          0             0.00
tests/symbolic/solid_state_physics.rs              48        6        25       17          0             0.00
tests/symbolic/real_roots.rs                       48        6        25       17          0             0.00
tests/symbolic/multi_valued.rs                     48        6        25       17          0             0.00
tests/symbolic/solve.rs                            48        6        25       17          0             0.00
tests/symbolic/special.rs                          48        6        25       17          0             0.00
tests/symbolic/special_functions.rs                48        6        25       17          0             0.00
tests/symbolic/group_theory.rs                     48        6        25       17          0             0.00
tests/symbolic/proof.rs                            48        6        25       17          0             0.00
benches/lib_bench.rs                               48       12         0       36          0             0.00
tests/symbolic/radicals.rs                         48        6        25       17          0             0.00
tests/symbolic/quantum_mechanics.rs                48        6        25       17          0             0.00
tests/symbolic/stats_probability.rs                48        6        25       17          0             0.00
tests/symbolic/stats_regression.rs                 48        6        25       17          0             0.00
tests/symbolic/quantum_field_theory.rs             48        6        25       17          0             0.00
tests/symbolic/topology.rs                         48        6        25       17          0             0.00
tests/symbolic/tensor.rs                           48        6        25       17          0             0.00
tests/symbolic/transforms.rs                       48        6        25       17          0             0.00
tests/symbolic/vector.rs                           48        6        25       17          0             0.00
tests/symbolic/matrix.rs                           48        6        25       17          0             0.00
tests/symbolic/vector_calculus.rs                  48        6        25       17          0             0.00
tests/output/plotting.rs                           48        6        25       17          0             0.00
tests/output/io.rs                                 48        6        25       17          0             0.00
tests/physics/mod.rs                               48        6        25       17          0             0.00
tests/output/latex.rs                              48        6        25       17          0             0.00
tests/symbolic/polynomial.rs                       48        6        25       17          0             0.00
tests/physics/physics_em.rs                        48        6        25       17          0             0.00
tests/symbolic/poly_factorization.rs               48        6        25       17          0             0.00
tests/symbolic/pde.rs                              48        6        25       17          0             0.00
tests/symbolic/optimize.rs                         48        6        25       17          0             0.00
tests/symbolic/number_theory.rs                    48        6        25       17          0             0.00
tests/physics/physics_cnm.rs                       48        6        25       17          0             0.00
tests/symbolic/relativity.rs                       48        6        25       17          0             0.00
tests/symbolic/grobner.rs                          48        6        25       17          0             0.00
tests/symbolic/stats_inference.rs                  48        6        25       17          0             0.00
tests/symbolic/stats_information_theory.rs         48        6        25       17          0             0.00
tests/output/typst.rs                              48        6        25       17          0             0.00
tests/symbolic/thermodynamics.rs                   48        6        25       17          0             0.00
tests/symbolic/stats.rs                            48        6        25       17          0             0.00
tests/symbolic/logic.rs                            48        6        25       17          0             0.00
tests/symbolic/lie_groups_and_algebras.rs          48        6        25       17          0             0.00
tests/numerical/complex_analysis.rs                48        6        25       17          0             0.00
tests/symbolic/integration.rs                      48        6        25       17          0             0.00
tests/symbolic/integral_equations.rs               48        6        25       17          0             0.00
tests/numerical/convergence.rs                     48        6        25       17          0             0.00
tests/numerical/coordinates.rs                     48        6        25       17          0             0.00
tests/numerical/computer_graphics.rs               48        6        25       17          0             0.00
tests/numerical/elementary.rs                      48        6        25       17          0             0.00
~symbolic/graph_isomorphism_and_coloring.rs        48        6        25       17          0             0.00
tests/symbolic/graph_algorithms.rs                 48        6        25       17          0             0.00
tests/numerical/calculus_of_variations.rs          48        6        25       17          0             0.00
tests/symbolic/graph.rs                            48        6        25       17          0             0.00
tests/symbolic/graph_operations.rs                 48        6        25       17          0             0.00
tests/numerical/combinatorics.rs                   48        6        25       17          0             0.00
tests/symbolic/complex_analysis.rs                 48        6        25       17          0             0.00
tests/numerical/differential_geometry.rs           48        6        25       17          0             0.00
~s/physics/physics_sim/linear_elasticity.rs        48        6        25       17          0             0.00
tests/symbolic/geometric_algebra.rs                48        6        25       17          0             0.00
tests/numerical/testing.rs                         48        6        25       17          0             0.00
tests/symbolic/functional_analysis.rs              48        6        25       17          0             0.00
~sts/symbolic/fractal_geometry_and_chaos.rs        48        6        25       17          0             0.00
tests/output/mod.rs                                48        6        25       17          0             0.00
tests/numerical/error_correction.rs                48        6        25       17          0             0.00
tests/numerical/finite_field.rs                    48        6        25       17          0             0.00
tests/symbolic/finite_field.rs                     48        6        25       17          0             0.00
tests/symbolic/error_correction_helper.rs          48        6        25       17          0             0.00
tests/symbolic/error_correction.rs                 48        6        25       17          0             0.00
tests/symbolic/elementary.rs                       48        6        25       17          0             0.00
tests/symbolic/convergence.rs                      48        6        25       17          0             0.00
~ts/numerical/fractal_geometry_and_chaos.rs        48        6        25       17          0             0.00
tests/symbolic/electromagnetism.rs                 48        6        25       17          0             0.00
tests/numerical/graph.rs                           48        6        25       17          0             0.00
tests/numerical/functional_analysis.rs             48        6        25       17          0             0.00
tests/symbolic/discrete_groups.rs                  48        6        25       17          0             0.00
tests/symbolic/differential_geometry.rs            48        6        25       17          0             0.00
tests/numerical/optimize.rs                        48        6        25       17          0             0.00
tests/symbolic/coordinates.rs                      48        6        25       17          0             0.00
tests/numerical/pde.rs                             48        6        25       17          0             0.00
tests/symbolic/classical_mechanics.rs              48        6        25       17          0             0.00
tests/symbolic/computer_graphics.rs                48        6        25       17          0             0.00
tests/symbolic/combinatorics.rs                    48        6        25       17          0             0.00
tests/symbolic/cas_foundations.rs                  48        6        25       17          0             0.00
tests/numerical/physics.rs                         48        6        25       17          0             0.00
tests/numerical/number_theory.rs                   48        6        25       17          0             0.00
tests/numerical/matrix.rs                          48        6        25       17          0             0.00
tests/symbolic/calculus_of_variations.rs           48        6        25       17          0             0.00
tests/symbolic/cad.rs                              48        6        25       17          0             0.00
tests/numerical/multi_valued.rs                    48        6        25       17          0             0.00
~physics/physics_sim/schrodinger_quantum.rs        48        6        25       17          0             0.00
tests/numerical/physics_md.rs                      48        6        25       17          0             0.00
~physics/physics_sim/navier_stokes_fluid.rs        48        6        25       17          0             0.00
tests/physics/physics_sim/mod.rs                   48        6        25       17          0             0.00
tests/numerical/physics_cfd.rs                     48        6        25       17          0             0.00
tests/numerical/signal.rs                          48        6        25       17          0             0.00
tests/numerical/solve.rs                           48        6        25       17          0             0.00
tests/numerical/special.rs                         48        6        25       17          0             0.00
~s/physics/physics_sim/ising_statistical.rs        48        6        25       17          0             0.00
tests/numerical/stats.rs                           48        6        25       17          0             0.00
~s/physics/physics_sim/gpe_superfluidity.rs        48        6        25       17          0             0.00
tests/numerical/real_roots.rs                      48        6        25       17          0             0.00
tests/numerical/sparse.rs                          48        6        25       17          0             0.00
~physics/physics_sim/geodesic_relativity.rs        48        6        25       17          0             0.00
tests/output/pretty_print.rs                       48        6        25       17          0             0.00
tests/numerical/series.rs                          48        6        25       17          0             0.00
tests/physics/physics_rkm.rs                       48        6        25       17          0             0.00
tests/numerical/transforms.rs                      48        6        25       17          0             0.00
~hysics/physics_sim/fdtd_electrodynamics.rs        48        6        25       17          0             0.00
tests/physics/physics_fdm.rs                       48        6        25       17          0             0.00
tests/physics/physics_sm.rs                        48        6        25       17          0             0.00
tests/numerical/vector.rs                          48        6        25       17          0             0.00
tests/numerical/physics_fea.rs                     48        6        25       17          0             0.00
tests/numerical/vector_calculus.rs                 48        6        25       17          0             0.00
tests/physics/physics_fem.rs                       48        6        25       17          0             0.00
tests/physics/physics_mm.rs                        48        6        25       17          0             0.00
tests/physics/physics_mtm.rs                       48        6        25       17          0             0.00
tests/physics/physics_fvm.rs                       48        6        25       17          0             0.00
tests/numerical/topology.rs                        48        6        25       17          0             0.00
tests/physics/physics_bem.rs                       48        6        25       17          0             0.00
tests/numerical/tensor.rs                          48        6        25       17          0             0.00
tests/ffi_blindings/mod.rs                         47        5        42        0          0             0.00
tests/ffi_apis/ffi_api.rs                          47        5        42        0          0             0.00
examples/plugins/example_plugin/src/lib.rs         47        7         9       31          1             3.23
tests/prelude.rs                                   47        5        42        0          0             0.00
tests/lib.rs                                       47        5        42        0          0             0.00
examples/calculus_demo.rs                          47       15         3       29          0             0.00
~/ffi_apis/symbolic_convergence_ffi/json.rs        44       14         5       25          1             4.00
~fi_apis/symbolic_convergence_ffi/handle.rs        44       15         5       24          1             4.17
benches/symbolic/vector_calculus.rs                43        8        23       12          0             0.00
benches/numerical/vector_calculus.rs               43        8        23       12          0             0.00
examples/number_theory_debug.rs                    43       10         1       32          0             0.00
benches/ffi_apis/ffi_api.rs                        43        8        23       12          0             0.00
benches/ffi_blindings/mod.rs                       43        8        23       12          0             0.00
benches/ffi_apis/mod.rs                            43        8        23       12          0             0.00
benches/numerical/calculus.rs                      43        8        23       12          0             0.00
benches/numerical/combinatorics.rs                 43        8        23       12          0             0.00
benches/numerical/complex_analysis.rs              43        8        23       12          0             0.00
benches/numerical/convergence.rs                   43        8        23       12          0             0.00
benches/numerical/coordinates.rs                   43        8        23       12          0             0.00
~enches/numerical/calculus_of_variations.rs        43        8        23       12          0             0.00
benches/numerical/differential_geometry.rs         43        8        23       12          0             0.00
benches/numerical/elementary.rs                    43        8        23       12          0             0.00
benches/numerical/computer_graphics.rs             43        8        23       12          0             0.00
benches/numerical/error_correction.rs              43        8        23       12          0             0.00
~es/numerical/fractal_geometry_and_chaos.rs        43        8        23       12          0             0.00
benches/numerical/finite_field.rs                  43        8        23       12          0             0.00
benches/numerical/geometric_algebra.rs             43        8        23       12          0             0.00
benches/numerical/functional_analysis.rs           43        8        23       12          0             0.00
benches/numerical/graph.rs                         43        8        23       12          0             0.00
benches/numerical/matrix.rs                        43        8        23       12          0             0.00
benches/numerical/number_theory.rs                 43        8        23       12          0             0.00
benches/numerical/integrate.rs                     43        8        23       12          0             0.00
benches/numerical/interpolate.rs                   43        8        23       12          0             0.00
benches/numerical/mod.rs                           43        8        23       12          0             0.00
benches/numerical/optimize.rs                      43        8        23       12          0             0.00
benches/numerical/ode.rs                           43        8        23       12          0             0.00
benches/numerical/multi_valued.rs                  43        8        23       12          0             0.00
benches/numerical/physics.rs                       43        8        23       12          0             0.00
benches/numerical/physics_fea.rs                   43        8        23       12          0             0.00
benches/numerical/physics_md.rs                    43        8        23       12          0             0.00
benches/numerical/physics_cfd.rs                   43        8        23       12          0             0.00
benches/numerical/polynomial.rs                    43        8        23       12          0             0.00
benches/numerical/pde.rs                           43        8        23       12          0             0.00
benches/numerical/real_roots.rs                    43        8        23       12          0             0.00
benches/numerical/signal.rs                        43        8        23       12          0             0.00
benches/numerical/solve.rs                         43        8        23       12          0             0.00
benches/numerical/sparse.rs                        43        8        23       12          0             0.00
benches/numerical/series.rs                        43        8        23       12          0             0.00
benches/numerical/special.rs                       43        8        23       12          0             0.00
benches/numerical/tensor.rs                        43        8        23       12          0             0.00
benches/numerical/testing.rs                       43        8        23       12          0             0.00
benches/numerical/transforms.rs                    43        8        23       12          0             0.00
benches/numerical/topology.rs                      43        8        23       12          0             0.00
benches/numerical/stats.rs                         43        8        23       12          0             0.00
benches/numerical/vector.rs                        43        8        23       12          0             0.00
~is/symbolic_convergence_ffi/bincode_api.rs        43       13         5       25          1             4.00
benches/output/io.rs                               43        8        23       12          0             0.00
benches/output/mod.rs                              43        8        23       12          0             0.00
benches/output/plotting.rs                         43        8        23       12          0             0.00
benches/output/pretty_print.rs                     43        8        23       12          0             0.00
benches/output/typst.rs                            43        8        23       12          0             0.00
benches/physics/mod.rs                             43        8        23       12          0             0.00
benches/physics/physics_bem.rs                     43        8        23       12          0             0.00
benches/physics/physics_cnm.rs                     43        8        23       12          0             0.00
benches/output/latex.rs                            43        8        23       12          0             0.00
benches/physics/physics_em.rs                      43        8        23       12          0             0.00
benches/physics/physics_fem.rs                     43        8        23       12          0             0.00
benches/physics/physics_fdm.rs                     43        8        23       12          0             0.00
benches/physics/physics_fvm.rs                     43        8        23       12          0             0.00
benches/physics/physics_mtm.rs                     43        8        23       12          0             0.00
~s/physics/physics_sim/ising_statistical.rs        43        8        23       12          0             0.00
~s/physics/physics_sim/linear_elasticity.rs        43        8        23       12          0             0.00
~s/physics/physics_sim/gpe_superfluidity.rs        43        8        23       12          0             0.00
benches/physics/physics_sim/mod.rs                 43        8        23       12          0             0.00
~physics/physics_sim/schrodinger_quantum.rs        43        8        23       12          0             0.00
benches/plugins/mod.rs                             43        8        23       12          0             0.00
benches/symbolic/cad.rs                            43        8        23       12          0             0.00
~physics/physics_sim/geodesic_relativity.rs        43        8        23       12          0             0.00
benches/symbolic/calculus.rs                       43        8        23       12          0             0.00
benches/physics/physics_mm.rs                      43        8        23       12          0             0.00
~physics/physics_sim/navier_stokes_fluid.rs        43        8        23       12          0             0.00
benches/symbolic/cas_foundations.rs                43        8        23       12          0             0.00
~hysics/physics_sim/fdtd_electrodynamics.rs        43        8        23       12          0             0.00
benches/symbolic/calculus_of_variations.rs         43        8        23       12          0             0.00
benches/physics/physics_sm.rs                      43        8        23       12          0             0.00
benches/physics/physics_rkm.rs                     43        8        23       12          0             0.00
benches/symbolic/classical_mechanics.rs            43        8        23       12          0             0.00
benches/symbolic/combinatorics.rs                  43        8        23       12          0             0.00
benches/symbolic/complex_analysis.rs               43        8        23       12          0             0.00
benches/symbolic/convergence.rs                    43        8        23       12          0             0.00
benches/symbolic/computer_graphics.rs              43        8        23       12          0             0.00
benches/symbolic/core.rs                           43        8        23       12          0             0.00
benches/symbolic/cryptography.rs                   43        8        23       12          0             0.00
benches/symbolic/vector.rs                         43        8        23       12          0             0.00
benches/symbolic/differential_geometry.rs          43        8        23       12          0             0.00
benches/symbolic/topology.rs                       43        8        23       12          0             0.00
benches/symbolic/transforms.rs                     43        8        23       12          0             0.00
benches/symbolic/thermodynamics.rs                 43        8        23       12          0             0.00
benches/symbolic/tensor.rs                         43        8        23       12          0             0.00
benches/symbolic/discrete_groups.rs                43        8        23       12          0             0.00
benches/symbolic/coordinates.rs                    43        8        23       12          0             0.00
benches/symbolic/stats_probability.rs              43        8        23       12          0             0.00
~enches/symbolic/error_correction_helper.rs        43        8        23       12          0             0.00
benches/symbolic/stats_regression.rs               43        8        23       12          0             0.00
benches/symbolic/finite_field.rs                   43        8        23       12          0             0.00
~nches/symbolic/stats_information_theory.rs        43        8        23       12          0             0.00
benches/symbolic/special_functions.rs              43        8        23       12          0             0.00
benches/symbolic/solve.rs                          43        8        23       12          0             0.00
benches/symbolic/functional_analysis.rs            43        8        23       12          0             0.00
benches/symbolic/stats_inference.rs                43        8        23       12          0             0.00
~hes/symbolic/fractal_geometry_and_chaos.rs        43        8        23       12          0             0.00
benches/symbolic/stats.rs                          43        8        23       12          0             0.00
benches/symbolic/electromagnetism.rs               43        8        23       12          0             0.00
benches/symbolic/special.rs                        43        8        23       12          0             0.00
benches/symbolic/geometric_algebra.rs              43        8        23       12          0             0.00
benches/symbolic/solid_state_physics.rs            43        8        23       12          0             0.00
benches/symbolic/elementary.rs                     43        8        23       12          0             0.00
benches/symbolic/simplify.rs                       43        8        23       12          0             0.00
benches/symbolic/relativity.rs                     43        8        23       12          0             0.00
benches/symbolic/series.rs                         43        8        23       12          0             0.00
benches/symbolic/rewriting.rs                      43        8        23       12          0             0.00
benches/symbolic/graph.rs                          43        8        23       12          0             0.00
benches/symbolic/proof.rs                          43        8        23       12          0             0.00
benches/symbolic/ode.rs                            43        8        23       12          0             0.00
benches/symbolic/graph_algorithms.rs               43        8        23       12          0             0.00
benches/symbolic/real_roots.rs                     43        8        23       12          0             0.00
benches/symbolic/radicals.rs                       43        8        23       12          0             0.00
benches/symbolic/polynomial.rs                     43        8        23       12          0             0.00
benches/symbolic/graph_operations.rs               43        8        23       12          0             0.00
benches/symbolic/quantum_mechanics.rs              43        8        23       12          0             0.00
benches/symbolic/quantum_field_theory.rs           43        8        23       12          0             0.00
benches/symbolic/group_theory.rs                   43        8        23       12          0             0.00
~symbolic/graph_isomorphism_and_coloring.rs        43        8        23       12          0             0.00
benches/symbolic/number_theory.rs                  43        8        23       12          0             0.00
benches/symbolic/grobner.rs                        43        8        23       12          0             0.00
benches/symbolic/poly_factorization.rs             43        8        23       12          0             0.00
benches/symbolic/pde.rs                            43        8        23       12          0             0.00
benches/symbolic/integration.rs                    43        8        23       12          0             0.00
benches/symbolic/optimize.rs                       43        8        23       12          0             0.00
benches/symbolic/numeric.rs                        43        8        23       12          0             0.00
benches/symbolic/multi_valued.rs                   43        8        23       12          0             0.00
benches/symbolic/mod.rs                            43        8        23       12          0             0.00
~enches/symbolic/lie_groups_and_algebras.rs        43        8        23       12          0             0.00
benches/symbolic/matrix.rs                         43        8        23       12          0             0.00
benches/symbolic/logic.rs                          43        8        23       12          0             0.00
benches/symbolic/integral_equations.rs             43        8        23       12          0             0.00
benches/symbolic/error_correction.rs               43        8        23       12          0             0.00
src/ffi_apis/constant_ffi/mod.rs                   41        2        33        6          0             0.00
~mbolic_unit_unification_ffi/bincode_api.rs        40       10         4       26          3            11.54
~i_apis/symbolic_simplify_dag_ffi/handle.rs        40        9        12       19          1             5.26
src/ffi_apis/symbolic_numeric_ffi/json.rs          40       13         4       23          4            17.39
~apis/symbolic_unit_unification_ffi/json.rs        39       11         4       24          3            12.50
tests/debug_json_test.rs                           39        9         0       30          0             0.00
~i_apis/symbolic_numeric_ffi/bincode_api.rs        39       12         4       23          4            17.39
tests/symbolic/core.rs                             37        6        12       19          1             5.26
src/compute/state.rs                               36        7         7       22          1             4.55
examples/debug_integration.rs                      35        8         1       26          0             0.00
~is/symbolic_unit_unification_ffi/handle.rs        33        9         4       20          1             5.00
tests/prelude_test.rs                              33       11         6       16          0             0.00
benches/rssn_benches.rs                            31        2         0       29          0             0.00
~ffi_apis/symbolic_simplify_dag_ffi/json.rs        30        9         2       19          2            10.53
~s/symbolic_simplify_dag_ffi/bincode_api.rs        29        8         2       19          2            10.53
benches/constant.rs                                29        5         0       24          0             0.00
examples/symbolic_differentiation_poly.rs          29        6         5       18          0             0.00
~_apis/symbolic_error_correction_ffi/mod.rs        27        1        23        3          0             0.00
src/compute/computable.rs                          27        3        15        9          0             0.00
src/output/mod.rs                                  27        1        21        5          0             0.00
examples/test_simplify_bug.rs                      26        7         0       19          0             0.00
src/numerical/pde.rs                               26        4        10       12          0             0.00
tests/compute_state_test.rs                        25        8         0       17          0             0.00
~rc/ffi_apis/symbolic_numeric_ffi/handle.rs        25        9         4       12          0             0.00
tests/lib_test.rs                                  22       10         0       12          0             0.00
examples/readme_demo.rs                            21        5         4       12          0             0.00
~rc/ffi_apis/symbolic_transforms_ffi/mod.rs        20        1        16        3          0             0.00
~/ffi_apis/symbolic_cryptography_ffi/mod.rs        20        1        16        3          0             0.00
benches/prelude_bench.rs                           20        5         0       15          0             0.00
~xamples/symbolic_differentiation_modern.rs        20        5         3       12          0             0.00
src/ffi_apis/compute_cache_ffi/mod.rs              18        2        10        6          0             0.00
src/ffi_apis/compute_state_ffi/mod.rs              18        2        10        6          0             0.00
src/physics/physics_sim/mod.rs                     18        1        10        7          0             0.00
~rc/ffi_apis/symbolic_elementary_ffi/mod.rs        18        2        10        6          0             0.00
tests/numerical/ode.rs                             17        1         0       16          0             0.00
src/ffi_apis/symbolic_handles_ffi/mod.rs           17        2         9        6          0             0.00
src/ffi_apis/symbolic_simplify_ffi/mod.rs          17        2         9        6          0             0.00
src/ffi_apis/symbolic_rewriting_ffi/mod.rs         17        2         9        6          0             0.00
src/ffi_apis/symbolic_ode_ffi/mod.rs               17        2         9        6          0             0.00
~pis/symbolic_integral_equations_ffi/mod.rs        17        2         9        6          0             0.00
~/ffi_apis/symbolic_simplify_dag_ffi/mod.rs        17        2         9        6          0             0.00
src/ffi_apis/symbolic_pde_ffi/mod.rs               17        2         9        6          0             0.00
src/ffi_apis/symbolic_calculus_ffi/mod.rs          17        2         9        6          0             0.00
~i_apis/symbolic_vector_calculus_ffi/mod.rs        17        2         9        6          0             0.00
~i_apis/symbolic_cas_foundations_ffi/mod.rs        17        2         9        6          0             0.00
~ymbolic_error_correction_helper_ffi/mod.rs        16        1        12        3          0             0.00
benches/compute_state.rs                           16        3         0       13          0             0.00
~c/ffi_apis/symbolic_integration_ffi/mod.rs        16        1        12        3          0             0.00
src/jit/mod.rs                                     16        3         2       11          0             0.00
~apis/symbolic_computer_graphics_ffi/mod.rs        16        1        12        3          0             0.00
src/plugins/mod.rs                                 15        1        11        3          0             0.00
~c/ffi_apis/numerical_transforms_ffi/mod.rs        13        1         9        3          0             0.00
~_apis/symbolic_unit_unification_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/symbolic_radicals_ffi/mod.rs          13        1         9        3          0             0.00
src/ffi_apis/plugins_ffi/mod.rs                    13        1         9        3          0             0.00
src/ffi_apis/symbolic_proof_ffi/mod.rs             13        1         9        3          0             0.00
~rc/ffi_apis/symbolic_real_roots_ffi/mod.rs        13        1         9        3          0             0.00
~fi_apis/numerical_combinatorics_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/physics_bem_ffi/mod.rs                13        1         9        3          0             0.00
~umerical_calculus_of_variations_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/nightly_ffi/mod.rs                    13        1         9        3          0             0.00
~apis/numerical_complex_analysis_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/numerical_calculus_ffi/mod.rs         13        1         9        3          0             0.00
~rc/ffi_apis/symbolic_polynomial_ffi/mod.rs        13        1         9        3          0             0.00
~fi_apis/physics_sim_schrodinger_ffi/mod.rs        13        1         9        3          0             0.00
~_apis/symbolic_electromagnetism_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/physics_em_ffi/mod.rs                 13        1         9        3          0             0.00
~pis/numerical_computer_graphics_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/numerical_vector_ffi/mod.rs           13        1         9        3          0             0.00
~_apis/physics_sim_navier_stokes_ffi/mod.rs        13        1         9        3          0             0.00
~apis/symbolic_quantum_mechanics_ffi/mod.rs        13        1         9        3          0             0.00
~_apis/numerical_vector_calculus_ffi/mod.rs        13        1         9        3          0             0.00
~s/symbolic_quantum_field_theory_ffi/mod.rs        13        1         9        3          0             0.00
~rc/ffi_apis/symbolic_relativity_ffi/mod.rs        13        1         9        3          0             0.00
~/ffi_apis/numerical_convergence_ffi/mod.rs        13        1         9        3          0             0.00
~pis/symbolic_poly_factorization_ffi/mod.rs        13        1         9        3          0             0.00
~s/physics_sim_linear_elasticity_ffi/mod.rs        13        1         9        3          0             0.00
~symbolic_calculus_of_variations_ffi/mod.rs        13        1         9        3          0             0.00
~is/symbolic_solid_state_physics_ffi/mod.rs        13        1         9        3          0             0.00
~/ffi_apis/numerical_coordinates_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/symbolic_vector_ffi/mod.rs            13        1         9        3          0             0.00
~c/ffi_apis/symbolic_coordinates_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/physics_sm_ffi/mod.rs                 13        1         9        3          0             0.00
~apis/numerical_error_correction_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/symbolic_series_ffi/mod.rs            13        1         9        3          0             0.00
~fi_apis/symbolic_thermodynamics_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/physics_fdm_ffi/mod.rs                13        1         9        3          0             0.00
~/ffi_apis/numerical_physics_fea_ffi/mod.rs        13        1         9        3          0             0.00
~/ffi_apis/numerical_physics_cfd_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/symbolic_optimize_ffi/mod.rs          13        1         9        3          0             0.00
~c/ffi_apis/symbolic_convergence_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/physics_sim_ising_ffi/mod.rs          13        1         9        3          0             0.00
src/ffi_apis/numerical_topology_ffi/mod.rs         13        1         9        3          0             0.00
src/ffi_apis/numerical_stats_ffi/mod.rs            13        1         9        3          0             0.00
src/ffi_apis/numerical_optimize_ffi/mod.rs         13        1         9        3          0             0.00
~ffi_apis/numerical_finite_field_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/numerical_special_ffi/mod.rs          13        1         9        3          0             0.00
src/ffi_apis/physics_fem_ffi/mod.rs                13        1         9        3          0             0.00
src/ffi_apis/physics_fvm_ffi/mod.rs                13        1         9        3          0             0.00
src/ffi_apis/numerical_tensor_ffi/mod.rs           13        1         9        3          0             0.00
~apis/symbolic_stats_probability_ffi/mod.rs        13        1         9        3          0             0.00
~numerical_differential_geometry_ffi/mod.rs        13        1         9        3          0             0.00
~ical_fractal_geometry_and_chaos_ffi/mod.rs        13        1         9        3          0             0.00
~i_apis/symbolic_stats_inference_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/symbolic_numeric_ffi/mod.rs           13        1         9        3          0             0.00
src/ffi_apis/numerical_ode_ffi/mod.rs              13        1         9        3          0             0.00
~ffi_apis/symbolic_number_theory_ffi/mod.rs        13        1         9        3          0             0.00
~is/symbolic_classical_mechanics_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/physics_cnm_ffi/mod.rs                13        1         9        3          0             0.00
src/ffi_apis/numerical_sparse_ffi/mod.rs           13        1         9        3          0             0.00
src/ffi_apis/numerical_signal_ffi/mod.rs           13        1         9        3          0             0.00
src/ffi_apis/physics_sim_gpe_ffi/mod.rs            13        1         9        3          0             0.00
~pis/numerical_geometric_algebra_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/numerical_solve_ffi/mod.rs            13        1         9        3          0             0.00
src/ffi_apis/symbolic_topology_ffi/mod.rs          13        1         9        3          0             0.00
src/ffi_apis/symbolic_solve_ffi/mod.rs             13        1         9        3          0             0.00
~mbolic_stats_information_theory_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/physics_mm_ffi/mod.rs                 13        1         9        3          0             0.00
~c/ffi_apis/physics_sim_geodesic_ffi/mod.rs        13        1         9        3          0             0.00
~c/ffi_apis/numerical_real_roots_ffi/mod.rs        13        1         9        3          0             0.00
~/ffi_apis/numerical_interpolate_ffi/mod.rs        13        1         9        3          0             0.00
~/ffi_apis/symbolic_multi_valued_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/symbolic_tensor_ffi/mod.rs            13        1         9        3          0             0.00
src/ffi_apis/numerical_series_ffi/mod.rs           13        1         9        3          0             0.00
src/ffi_apis/physics_rkm_ffi/mod.rs                13        1         9        3          0             0.00
src/ffi_apis/numerical_matrix_ffi/mod.rs           13        1         9        3          0             0.00
~_apis/symbolic_stats_regression_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/numerical_physics_ffi/mod.rs          13        1         9        3          0             0.00
src/ffi_apis/numerical_graph_ffi/mod.rs            13        1         9        3          0             0.00
~ffi_apis/numerical_multi_valued_ffi/mod.rs        13        1         9        3          0             0.00
~fi_apis/numerical_number_theory_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/physics_mtm_ffi/mod.rs                13        1         9        3          0             0.00
src/ffi_apis/physics_sim_fdtd_ffi/mod.rs           13        2         8        3          0             0.00
src/ffi_apis/symbolic_stats_ffi/mod.rs             13        1         9        3          0             0.00
~c/ffi_apis/numerical_physics_md_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/symbolic_matrix_ffi/mod.rs            13        1         9        3          0             0.00
~c/ffi_apis/numerical_polynomial_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/symbolic_graph_ffi/mod.rs             12        1         8        3          0             0.00
~i_apis/symbolic_discrete_groups_ffi/mod.rs        12        1         8        3          0             0.00
src/ffi_apis/symbolic_logic_ffi/mod.rs             12        1         8        3          0             0.00
~/symbolic_differential_geometry_ffi/mod.rs        12        1         8        3          0             0.00
~/ffi_apis/symbolic_finite_field_ffi/mod.rs        12        1         8        3          0             0.00
~is/symbolic_functional_analysis_ffi/mod.rs        12        1         8        3          0             0.00
~olic_fractal_geometry_and_chaos_ffi/mod.rs        12        1         8        3          0             0.00
~_apis/symbolic_complex_analysis_ffi/mod.rs        12        1         8        3          0             0.00
~_apis/symbolic_graph_operations_ffi/mod.rs        12        1         8        3          0             0.00
~_graph_isomorphism_and_coloring_ffi/mod.rs        12        1         8        3          0             0.00
src/ffi_apis/symbolic_cad_ffi/mod.rs               12        1         8        3          0             0.00
src/ffi_apis/symbolic_grobner_ffi/mod.rs           11        1         7        3          0             0.00
~rc/ffi_apis/numerical_integrate_ffi/mod.rs        11        1         7        3          0             0.00
src/ffi_apis/jit_ffi/mod.rs                        11        1         8        2          0             0.00
~_apis/symbolic_graph_algorithms_ffi/mod.rs        11        1         7        3          0             0.00
~s/numerical_functional_analysis_ffi/mod.rs        11        1         7        3          0             0.00
~ffi_apis/symbolic_combinatorics_ffi/mod.rs        11        1         7        3          0             0.00
~rc/ffi_apis/symbolic_lie_groups_ffi/mod.rs        11        1         7        3          0             0.00
~apis/symbolic_geometric_algebra_ffi/mod.rs        11        1         7        3          0             0.00
~/ffi_apis/symbolic_group_theory_ffi/mod.rs        11        1         7        3          0             0.00
~c/ffi_apis/numerical_elementary_ffi/mod.rs        10        1         6        3          0             0.00
tests/symbolic/mod.rs                               8        2         2        4          0             0.00
src/verification/mod.rs                             8        1         6        1          0             0.00
src/ffi_blindings/mod.rs                            4        1         2        1          0             0.00
src/input/mod.rs                                    3        0         2        1          0             0.00
tests/numerical/mod.rs                              2        0         0        2          0             0.00
src/nightly/mod.rs                                  1        0         0        1          0             0.00
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Markdown                                 20      5554      837         0     4717          0             0.00
MaxLine / MeanLine                      655        51
(ULOC)                                           3374
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
CODE_STASTICS.md                                 1977        8         0     1969          0             0.00
ATTRIBUTIONS.md                                   874      220         0      654          0             0.00
docs/ffi_three_version_complete.md                336       70         0      266          0             0.00
docs/ffi_strategy.md                              293       53         0      240          0             0.00
docs/ast_to_dag_migration.md                      285       49         0      236          0             0.00
docs/ast_to_dag_phase2_complete.md                261       61         0      200          0             0.00
docs/ast_to_dag_phase1_complete.md                195       44         0      151          0             0.00
docs/tasks_20251119-20260129.md                   191       16         0      175          0             0.00
CONTRIBUTING.md                                   190       53         0      137          0             0.00
docs/progress_report.md                           186       41         0      145          0             0.00
~documentation_and_simplify_enhancements.md       185       49         0      136          0             0.00
README.md                                         140       51         0       89          0             0.00
docs/list_variants_refactoring.md                 120       33         0       87          0             0.00
docs/TECHNICAL_DEBT.md                             80       15         0       65          0             0.00
ARCHITECTURE.md                                    66       22         0       44          0             0.00
.github/ISSUE_TEMPLATE/bug_report.md               49       15         0       34          0             0.00
.github/ISSUE_TEMPLATE/feature_request.md          45       17         0       28          0             0.00
docs/TECHNICAL_DEBT_FIX_SUMMARY.md                 38        7         0       31          0             0.00
CODE_OF_CONDUCT.md                                 32        9         0       23          0             0.00
SECURITY.md                                        11        4         0        7          0             0.00
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
TOML                                      8       385       18        15      352          2             1.13
MaxLine / MeanLine                      126        22
(ULOC)                                            326
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Cargo.toml                                        196       12         7      177          2             1.13
cbindgen.toml                                     104        0         1      103          0             0.00
rustfmt.toml                                       38        0         3       35          0             0.00
kani.toml                                          15        2         0       13          0             0.00
about.toml                                         12        1         0       11          0             0.00
examples/plugins/example_plugin/Cargo.toml         10        2         0        8          0             0.00
clippy.toml                                         5        0         4        1          0             0.00
.cargo/config.toml                                  5        1         0        4          0             0.00
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Python                                    6       630       64        54      512         36            54.45
MaxLine / MeanLine                      128        37
(ULOC)                                            411
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
scripts/again.py                                  203       23        15      165         15             9.09
scripts/rename_prelude.py                         117        1         1      115          0             0.00
scripts/dedupe.py                                  99        1         1       97          1             1.03
scripts/update_vis.py                              96       14        11       71         11            15.49
scripts/count_rs.py                                64       15        18       31          8            25.81
examples/ode_solve.py                              51       10         8       33          1             3.03
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
SVG                                       6      4589      108         2     4479        790           145.40
MaxLine / MeanLine                   324860       218
(ULOC)                                           1553
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
doc/ode_lorenz_attractor_1_fg.svg                2134       23         2     2109        130             6.16
doc/fem_bridge.svg                                491       17         0      474        132            27.85
doc/fvm_advection.svg                             491       17         0      474        132            27.85
doc/ode_lorenz_attractor_3_fg.svg                 491       17         0      474        132            27.85
doc/ode_lorenz_attractor_4_fg.svg                 491       17         0      474        132            27.85
doc/ode_lorenz_attractor_2_fg.svg                 491       17         0      474        132            27.85
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Shell                                     5       327       50        59      218         22            40.97
MaxLine / MeanLine                      114        32
(ULOC)                                            254
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
release_cpu_targets_linux.sh                       98       21        25       52          5             9.62
test_features.sh                                   93       15        17       61          5             8.20
build_all.sh                                       72       10        10       52         11            21.15
collect_binaries.sh                                59        3         6       50          1             2.00
pre-commit.sh                                       5        1         1        3          0             0.00
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
YAML                                      5       351       57        44      250          0             0.00
MaxLine / MeanLine                      293        37
(ULOC)                                            257
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
.gitea/workflows/ci.yml                           164       30         1      133          0             0.00
.github/workflows/codeql.yml                      101       10        43       48          0             0.00
.github/workflows/ci.yml                           50       11         0       39          0             0.00
.github/workflows/dependency-review.yml            23        5         0       18          0             0.00
.github/dependabot.yml                             13        1         0       12          0             0.00
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Batch                                     3       169       24        12      133         28            40.22
MaxLine / MeanLine                      101        27
(ULOC)                                            109
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
build_all.bat                                     101       13         4       84         21            25.00
collect_binaries.bat                               65       11         8       46          7            15.22
pre-commit.bat                                      3        0         0        3          0             0.00
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
XML                                       3        27        0         0       27          0             0.00
MaxLine / MeanLine                      113        37
(ULOC)                                             18
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
.idea/vcs.xml                                      13        0         0       13          0             0.00
.idea/modules.xml                                   8        0         0        8          0             0.00
.idea/rust.xml                                      6        0         0        6          0             0.00
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Fortran Modern                            2       256       50        44      162          9            10.84
MaxLine / MeanLine                      105        35
(ULOC)                                            171
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
examples/fortran/test.f90                         167       34        33      100          6             6.00
examples/numerical_vector_ffi_demo.f90             89       16        11       62          3             4.84
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
C Header                                  1     40700     2578     28171     9951          0             0.00
MaxLine / MeanLine                      142        33
(ULOC)                                           9446
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
rssn.h                                          40700     2578     28171     9951          0             0.00
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
C++                                       1        68       12         6       50          7            14.00
MaxLine / MeanLine                       81        25
(ULOC)                                             52
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
examples/numerical_vector_ffi_demo.cpp             68       12         6       50          7            14.00
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
C++ Header                                1     40786     2587     28184    10015         22             0.22
MaxLine / MeanLine                      142        32
(ULOC)                                           9471
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
rssn.hpp                                        40786     2587     28184    10015         22             0.22
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
HTML                                      1         7        0         0        7          0             0.00
MaxLine / MeanLine                      113        37
(ULOC)                                              7
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
doc-header.html                                     7        0         0        7          0             0.00
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Handlebars                                1        16        3         0       13          0             0.00
MaxLine / MeanLine                      103        23
(ULOC)                                             14
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
about.hbs                                          16        3         0       13          0             0.00
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
License                                   1        73       32         0       41          0             0.00
MaxLine / MeanLine                      952       139
(ULOC)                                             42
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
LICENSE                                            73       32         0       41          0             0.00
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Total                                  1268    436748    74088    112824   249836      13760          4400.22
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Unique Lines of Code (ULOC)                    111913
DRYness %                                        0.26
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Estimated Cost to Develop (semi-detached) $32,727,207
Estimated Schedule Effort (semi-detached) 40.75 months
Estimated People Required (semi-detached) 71.34
Processed 12627225 bytes, 12.627 megabytes (SI)
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
```

## The src directory

```text
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Language                              Files     Lines   Blanks  Comments     Code Complexity Complexity/Lines
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Rust                                    686    975987    76690     53791   845506      36303          3578.78
MaxLine / MeanLine                      278        37
(ULOC)                                         158734
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Total                                   686    975987    76690     53791   845506      36303          3578.78
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Unique Lines of Code (ULOC)                    158734
DRYness %                                        0.16
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Estimated Cost to Develop (semi-detached) $128,205,315
Estimated Schedule Effort (semi-detached) 65.72 months
Estimated People Required (semi-detached) 173.30
Processed 37382330 bytes, 37.382 megabytes (SI)
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
```

```text
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Language                              Files     Lines   Blanks  Comments     Code Complexity Complexity/Lines
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Rust                                    686    975987    76690     53791   845506      36303          3578.78
MaxLine / MeanLine                      278        37
(ULOC)                                         158734
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
src/full.rs                                    406339    25984      3623   376732      12166             3.23
src/full_bak.rs                                301281        0      3623   297658      12167             4.09
src/ffi_apis/ffi_api.rs                          7963     1386      1154     5423        396             7.30
src/prelude.rs                                   6120      588       553     4979          0             0.00
src/symbolic/calculus.rs                         5272      834       283     4155        407             9.80
src/input/parser.rs                              4823       85        50     4688         25             0.53
src/symbolic/simplify_dag.rs                     3953      690       402     2861        380            13.28
src/symbolic/pde.rs                              3768      739       370     2659        261             9.82
src/symbolic/simplify.rs                         3471      508       255     2708        228             8.42
src/symbolic/graph_algorithms.rs                 3228      711       368     2149        301            14.01
src/symbolic/ode.rs                              2814      500       287     2027        200             9.87
src/numerical/matrix.rs                          2565      590       280     1695        226            13.33
src/symbolic/core/api.rs                         2546      458       333     1755         69             3.93
src/nightly/matrix.rs                            2546      578       277     1691        224            13.25
src/symbolic/core/expr_impl.rs                   2540      430        80     2030        114             5.62
src/symbolic/polynomial.rs                       2490      449       440     1601        151             9.43
src/numerical/computer_graphics.rs               2360      429       242     1689         48             2.84
src/symbolic/solve.rs                            2258      385       190     1683        141             8.38
src/numerical/physics.rs                         2183      493       414     1276        100             7.84
src/symbolic/matrix.rs                           2160      466       289     1405        193            13.74
src/numerical/error_correction.rs                2064      478       556     1030        153            14.85
~rc/numerical/fractal_geometry_and_chaos.rs      1921      409       485     1027         83             8.08
src/symbolic/special_functions.rs                1778      278       555      945        108            11.43
src/numerical/testing.rs                         1754      305       136     1313        123             9.37
src/symbolic/core/to_expr.rs                     1608      418        42     1148        137            11.93
src/symbolic/computer_graphics.rs                1599      172       264     1163         16             1.38
src/numerical/physics_fea.rs                     1499      304       205      990         33             3.33
src/symbolic/coordinates.rs                      1454      249       238      967         63             6.51
src/symbolic/special.rs                          1453      229       674      550         54             9.82
src/output/plotting.rs                           1422      275        68     1079         34             3.15
src/symbolic/number_theory.rs                    1392      277       160      955        120            12.57
src/numerical/physics_cfd.rs                     1378      335       213      830         59             7.11
src/symbolic/core/ast_impl.rs                    1374      147        21     1206         31             2.57
src/symbolic/error_correction_helper.rs          1331      283       238      810         99            12.22
src/ffi_apis/numerical_physics_ffi/json.rs       1321      162       529      630         17             2.70
src/numerical/physics_md.rs                      1315      296       241      778         58             7.46
src/symbolic/integration.rs                      1313      246       118      949         40             4.21
~is/numerical_computer_graphics_ffi/json.rs      1301      205       223      873         18             2.06
src/symbolic/topology.rs                         1299      265       239      795         78             9.81
src/ffi_apis/numerical_vector_ffi/json.rs        1276      162       193      921         45             4.89
src/physics/physics_rkm.rs                       1245      261        80      904         47             5.20
~i_apis/symbolic_cryptography_ffi/handle.rs      1227      233       475      519         50             9.63
src/symbolic/poly_factorization.rs               1197      242       145      810         59             7.28
src/physics/physics_fdm.rs                       1175      277        84      814         78             9.58
src/symbolic/stats_probability.rs                1172      200        66      906          9             0.99
src/symbolic/combinatorics.rs                    1167      219       200      748         57             7.62
src/symbolic/logic.rs                            1145      202       109      834         84            10.07
src/symbolic/cas_foundations.rs                  1142      189        70      883         65             7.36
src/symbolic/transforms.rs                       1131      134       229      768         38             4.95
src/symbolic/grobner.rs                          1113      164       412      537         45             8.38
src/numerical/special.rs                         1112      334        75      703        101            14.37
src/ffi_apis/numerical_special_ffi/json.rs       1095      111       476      508         16             3.15
src/physics/physics_fvm.rs                       1091      258        74      759         56             7.38
src/symbolic/tensor.rs                           1085      202       204      679         66             9.72
src/numerical/stats.rs                           1076      292       166      618         45             7.28
~is/symbolic_complex_analysis_ffi/handle.rs      1039      294       253      492         36             7.32
~c/ffi_apis/symbolic_transforms_ffi/json.rs      1038      283       250      505         22             4.36
~ffi_apis/symbolic_transforms_ffi/handle.rs      1035      220       290      525         38             7.24
src/symbolic/finite_field.rs                     1020      185       122      713         51             7.15
~apis/symbolic_complex_analysis_ffi/json.rs       991      259       251      481         35             7.28
~cal_fractal_geometry_and_chaos_ffi/json.rs       960      141       191      628         15             2.39
~_apis/numerical_special_ffi/bincode_api.rs       954       94       450      410         16             3.90
src/ffi_apis/numerical_stats_ffi/json.rs          954      105       391      458         13             2.84
~pis/numerical_error_correction_ffi/json.rs       954      142       198      614         21             3.42
~mbolic_complex_analysis_ffi/bincode_api.rs       949      232       277      440         35             7.95
src/output/pretty_print.rs                        943      202        10      731         65             8.89
src/output/io.rs                                  938      208       118      612         42             6.86
~pis/symbolic_transforms_ffi/bincode_api.rs       935      261       103      571         22             3.85
src/numerical/optimize.rs                         911      172        72      667         32             4.80
src/ffi_apis/symbolic_special_ffi/json.rs         902      212       243      447         56            12.53
src/symbolic/unit_unification.rs                  898       94       119      685         28             4.09
~s/symbolic_special_functions_ffi/handle.rs       897      190       233      474         51            10.76
src/symbolic/graph_operations.rs                  886      149       246      491         63            12.83
~pis/symbolic_special_functions_ffi/json.rs       882      214       233      435         33             7.59
~rc/ffi_apis/numerical_matrix_ffi/handle.rs       876      214       212      450         48            10.67
src/symbolic/core/dag_mgr.rs                      867      116       240      511         34             6.65
src/symbolic/core/expr.rs                         857       49       371      437          0             0.00
src/symbolic/error_correction.rs                  854      183       249      422         52            12.32
src/symbolic/complex_analysis.rs                  852      162       115      575         20             3.48
src/symbolic/series.rs                            849      126       164      559         32             5.72
~s/symbolic_computer_graphics_ffi/handle.rs       846      174       201      471         50            10.62
src/ffi_apis/nightly_ffi/handle.rs                845      223       172      450         48            10.67
src/compute/engine.rs                             832      103       327      402         11             2.74
~i_apis/numerical_vector_ffi/bincode_api.rs       830      119       129      582         30             5.15
~bolic_special_functions_ffi/bincode_api.rs       828      283       104      441         33             7.48
src/plugins/manager.rs                            822      151       103      568         23             4.05
~ffi_apis/numerical_physics_fea_ffi/json.rs       810       96       313      401          8             2.00
src/physics/physics_fem.rs                        804      190        66      548         56            10.22
~c/ffi_apis/numerical_physics_ffi/handle.rs       800      112       403      285          0             0.00
src/symbolic/cryptography.rs                      797      136       175      486         41             8.44
src/physics/physics_mtm.rs                        797      184        52      561         46             8.20
src/numerical/differential_geometry.rs            780      166        77      537         67            12.48
~/ffi_apis/numerical_physics_md_ffi/json.rs       779       94       238      447         10             2.24
src/symbolic/elementary.rs                        777      161        56      560         21             3.75
~ffi_apis/numerical_physics_cfd_ffi/json.rs       763       87       309      367          7             1.91
src/jit/engine.rs                                 755      130        43      582         51             8.76
~c/ffi_apis/symbolic_calculus_ffi/handle.rs       753      147       185      421         49            11.64
~erical_error_correction_ffi/bincode_api.rs       740      125       134      481         21             4.37
src/symbolic/integral_equations.rs                739      100       149      490          7             1.43
src/ffi_apis/symbolic_pde_ffi/handle.rs           734      135       108      491         64            13.03
src/physics/physics_sm.rs                         729      179        33      517         10             1.93
~ffi_apis/symbolic_lie_groups_ffi/handle.rs       723      128       346      249         14             5.62
~i_apis/symbolic_group_theory_ffi/handle.rs       688      117       285      286          8             2.80
src/symbolic/geometric_algebra.rs                 687      138       115      434         42             9.68
src/symbolic/cad.rs                               685      156        49      480         45             9.38
~fi_apis/numerical_stats_ffi/bincode_api.rs       683       79       298      306         11             3.59
~i_apis/symbolic_special_ffi/bincode_api.rs       674      183        48      443         56            12.64
src/numerical/sparse.rs                           660      169        90      401         38             9.48
~ffi_apis/symbolic_cryptography_ffi/json.rs       645      136       118      391         32             8.18
~s/numerical_error_correction_ffi/handle.rs       645      136       163      346         26             7.51
~s/symbolic_quantum_mechanics_ffi/handle.rs       638      141       160      337         34            10.09
~c/ffi_apis/numerical_special_ffi/handle.rs       635       91       347      197          0             0.00
~s/symbolic_stats_probability_ffi/handle.rs       635      169       165      301         12             3.99
~physics/physics_sim/navier_stokes_fluid.rs       633      140        82      411         46            11.19
src/physics/physics_em.rs                         632      157       102      373          9             2.41
src/numerical/geometric_algebra.rs                630       48        51      531          6             1.13
~rc/ffi_apis/numerical_vector_ffi/handle.rs       625      152       123      350         34             9.71
src/constant.rs                                   617       88        81      448          0             0.00
src/symbolic/quantum_mechanics.rs                 616      125        53      438          0             0.00
src/ffi_apis/numerical_stats_ffi/handle.rs        607      141       125      341         55            16.13
~is/symbolic_error_correction_ffi/handle.rs       607      125       165      317         29             9.15
src/symbolic/vector.rs                            602       86       145      371          4             1.08
src/symbolic/rewriting.rs                         601      111        45      445         46            10.34
src/physics/physics_bem.rs                        579      129        75      375         24             6.40
src/ffi_apis/nightly_ffi/json.rs                  576       99        61      416         26             6.25
src/ffi_apis/numerical_matrix_ffi/json.rs         576       99        61      416         26             6.25
src/physics/physics_cnm.rs                        573      166        40      367         33             8.99
~pis/symbolic_vector_calculus_ffi/handle.rs       571      121        16      434         52            11.98
src/symbolic/solid_state_physics.rs               571       77       100      394          0             0.00
src/numerical/polynomial.rs                       569      128        62      379         33             8.71
src/symbolic/fractal_geometry_and_chaos.rs        568      116        99      353         14             3.97
src/ffi_apis/symbolic_ode_ffi/handle.rs           567       90        96      381         51            13.39
~pis/symbolic_stats_probability_ffi/json.rs       567      163       162      242         10             4.13
~c/ffi_apis/symbolic_ode_ffi/bincode_api.rs       564      112         9      443         47            10.61
~is/numerical_geometric_algebra_ffi/json.rs       563       79        97      387         17             4.39
src/numerical/graph.rs                            563      148        77      338         45            13.31
src/numerical/integrate.rs                        556      110       167      279         28            10.04
src/numerical/elementary.rs                       555      152       113      290         23             7.93
src/symbolic/lie_groups_and_algebras.rs           555       97       139      319         13             4.08
~symbolic_functional_analysis_ffi/handle.rs       551      104       225      222          4             1.80
src/numerical/finite_field.rs                     550      122        77      351         41            11.68
src/numerical/complex_analysis.rs                 550       78        93      379          3             0.79
src/symbolic/proof.rs                             549      117        49      383         38             9.92
~olic_error_correction_helper_ffi/handle.rs       549      112       116      321         25             7.79
src/ffi_apis/symbolic_ode_ffi/json.rs             543      112         9      422         47            11.14
src/physics/physics_mm.rs                         541      128        43      370         21             5.68
~_apis/symbolic_calculus_ffi/bincode_api.rs       537      123        13      401         35             8.73
src/ffi_apis/symbolic_calculus_ffi/json.rs        534      123        13      398         35             8.79
src/numerical/interpolate.rs                      530      114       145      271         31            11.44
~i_apis/numerical_combinatorics_ffi/json.rs       514       82       109      323         10             3.10
~apis/symbolic_error_correction_ffi/json.rs       512      113       135      264         31            11.74
~s/symbolic_cryptography_ffi/bincode_api.rs       511      117        17      377         32             8.49
src/numerical/vector.rs                           511      107       129      275         39            14.18
~c/ffi_apis/symbolic_lie_groups_ffi/json.rs       511       80       237      194         15             7.73
src/symbolic/differential_geometry.rs             508       77        67      364         20             5.49
~is/symbolic_graph_algorithms_ffi/handle.rs       501      132        57      312         20             6.41
src/numerical/tensor.rs                           496      103        63      330         25             7.58
src/symbolic/discrete_groups.rs                   490       98        58      334         26             7.78
src/symbolic/graph.rs                             490       93       108      289         14             4.84
src/ffi_apis/physics_rkm_ffi/json.rs              488       68       142      278          8             2.88
~pis/symbolic_lie_groups_ffi/bincode_api.rs       479       72       249      158         15             9.49
~s/symbolic_group_theory_ffi/bincode_api.rs       477       70       252      155         15             9.68
src/numerical/vector_calculus.rs                  472      102       120      250         17             6.80
~ffi_apis/symbolic_group_theory_ffi/json.rs       469       71       246      152         15             9.87
~_apis/numerical_physics_ffi/bincode_api.rs       469       58       200      211          7             3.32
src/ffi_apis/symbolic_graph_ffi/handle.rs         464      121        38      305         23             7.54
src/ffi_apis/symbolic_proof_ffi/handle.rs         459       87        63      309         39            12.62
src/symbolic/real_roots.rs                        458       97        73      288         29            10.07
~/ffi_apis/numerical_optimize_ffi/handle.rs       458       83       128      247         12             4.86
~apis/symbolic_graph_algorithms_ffi/json.rs       457       96       130      231         15             6.49
~/numerical_computer_graphics_ffi/handle.rs       455      103        76      276         19             6.88
src/output/latex.rs                               451       85         9      357         17             4.76
~bolic_stats_probability_ffi/bincode_api.rs       447      148        58      241         10             4.15
src/symbolic/convergence.rs                       445       77       121      247         42            17.00
src/ffi_apis/numerical_graph_ffi/json.rs          444       84        73      287          9             3.14
~ctal_geometry_and_chaos_ffi/bincode_api.rs       444       73        81      290          9             3.10
src/symbolic/multi_valued.rs                      444       74       121      249          0             0.00
src/verification/symbolic_core.rs                 439      162        19      258          6             2.33
src/symbolic/radicals.rs                          437       85        38      314         34            10.83
~pis/symbolic_computer_graphics_ffi/json.rs       436       86       100      250         17             6.80
~fi_apis/symbolic_series_ffi/bincode_api.rs       435      110        39      286          8             2.80
~mbolic_graph_algorithms_ffi/bincode_api.rs       432       96       104      232         15             6.47
src/ffi_apis/physics_fdm_ffi/json.rs              430       49       145      236          9             3.81
src/ffi_apis/symbolic_series_ffi/handle.rs        428      139        39      250          8             3.20
~/ffi_apis/symbolic_rewriting_ffi/handle.rs       428       92       123      213         23            10.80
src/ffi_apis/symbolic_series_ffi/json.rs          428      111        39      278          8             2.88
src/numerical/calculus.rs                         421       77       158      186         12             6.45
~c_fractal_geometry_and_chaos_ffi/handle.rs       421      111        14      296         25             8.45
~/numerical_geometric_algebra_ffi/handle.rs       419      116        89      214         24            11.21
src/symbolic/handles.rs                           415       41       274      100          1             1.00
src/symbolic/group_theory.rs                      414       96        27      291         34            11.68
src/symbolic/stats_information_theory.rs          412       66       107      239         12             5.02
~mbolic_error_correction_ffi/bincode_api.rs       411       98        46      267         30            11.24
src/numerical/real_roots.rs                       406      106        65      235         54            22.98
src/symbolic/classical_mechanics.rs               406       73        58      275          0             0.00
src/numerical/topology.rs                         403      102        40      261         30            11.49
src/ffi_apis/symbolic_pde_ffi/json.rs             403       81         8      314         26             8.28
src/ffi_apis/symbolic_vector_ffi/json.rs          402      115        32      255         27            10.59
~/ffi_apis/numerical_calculus_ffi/handle.rs       402       92        36      274         31            11.31
~rical_geometric_algebra_ffi/bincode_api.rs       402       63        65      274         10             3.65
~c/physics/physics_sim/linear_elasticity.rs       401       86        50      265         12             4.53
~is/symbolic_coordinates_ffi/bincode_api.rs       400      103        37      260         17             6.54
src/lib.rs                                        399       10       295       94          1             1.06
~numerical_combinatorics_ffi/bincode_api.rs       395       72        73      250         10             4.00
src/symbolic/functional_analysis.rs               389       65        60      264          7             2.65
~fi_apis/symbolic_coordinates_ffi/handle.rs       387       40       153      194          8             4.12
src/numerical/ode.rs                              386       72        60      254          7             2.76
src/ffi_apis/numerical_graph_ffi/handle.rs        385       96        75      214         16             7.48
src/numerical/convergence.rs                      384       91       103      190         25            13.16
~ffi_apis/symbolic_graph_ffi/bincode_api.rs       380       82        96      202         12             5.94
~physics/physics_sim/schrodinger_quantum.rs       375       75        57      243         13             5.35
~rc/ffi_apis/numerical_sparse_ffi/handle.rs       370       89        78      203         18             8.87
src/symbolic/optimize.rs                          369       63        66      240         19             7.92
~/ffi_apis/symbolic_coordinates_ffi/json.rs       366       47       136      183         21            11.48
~/symbolic_integral_equations_ffi/handle.rs       365       79        10      276         38            13.77
~is/numerical_vector_calculus_ffi/handle.rs       364       83        33      248         30            12.10
~fi_apis/numerical_graph_ffi/bincode_api.rs       363       77        49      237          9             3.80
~c/ffi_apis/symbolic_pde_ffi/bincode_api.rs       359       70         7      282         22             7.80
src/symbolic/numeric.rs                           359       40        31      288         18             6.25
src/ffi_apis/common.rs                            359       89        65      205         13             6.34
src/symbolic/relativity.rs                        358       66        39      253          0             0.00
~apis/numerical_combinatorics_ffi/handle.rs       352       88        73      191         12             6.28
~rical_computer_graphics_ffi/bincode_api.rs       352       67        49      236          6             2.54
src/numerical/number_theory.rs                    349      105        47      197         49            24.87
~i_apis/symbolic_multi_valued_ffi/handle.rs       347      112         9      226         19             8.41
src/symbolic/stats_inference.rs                   346       44        34      268          0             0.00
~bolic_computer_graphics_ffi/bincode_api.rs       346       73        16      257         17             6.61
src/symbolic/quantum_field_theory.rs              346       66        23      257          2             0.78
~ffi_apis/symbolic_elementary_ffi/handle.rs       345       83        85      177         10             5.65
~i_apis/numerical_interpolate_ffi/handle.rs       344       80        51      213         16             7.51
src/ffi_apis/numerical_tensor_ffi/json.rs         341       53        25      263         12             4.56
~fi_apis/numerical_multi_valued_ffi/json.rs       341       48        73      220          5             2.27
~symbolic/graph_isomorphism_and_coloring.rs       340       72        42      226         29            12.83
~mbolic_error_correction_helper_ffi/json.rs       335       67        68      200         17             8.50
src/symbolic/thermodynamics.rs                    333       57        37      239          1             0.42
~lic_stats_information_theory_ffi/handle.rs       331       92        84      155         16            10.32
~rc/ffi_apis/numerical_tensor_ffi/handle.rs       327       79        57      191         16             8.38
~l_fractal_geometry_and_chaos_ffi/handle.rs       327       70        60      197         14             7.11
src/symbolic/electromagnetism.rs                  324       56        45      223          0             0.00
src/output/typst.rs                               323       62         4      257         12             4.67
src/symbolic/vector_calculus.rs                   321       51        62      208          0             0.00
~is/symbolic_graph_operations_ffi/handle.rs       321      102         8      211         19             9.00
~rc/ffi_apis/symbolic_special_ffi/handle.rs       320       91        48      181          0             0.00
~c/ffi_apis/symbolic_topology_ffi/handle.rs       319       93        19      207         19             9.18
~_apis/numerical_multi_valued_ffi/handle.rs       319       75        65      179          5             2.79
~mbolic_differential_geometry_ffi/handle.rs       313       87        15      211         20             9.48
~c/physics/physics_sim/ising_statistical.rs       310       76        35      199         18             9.05
~apis/numerical_vector_calculus_ffi/json.rs       307       36        97      174          6             3.45
src/ffi_apis/mod.rs                               306        3       139      164          0             0.00
~s/symbolic_finite_field_ffi/bincode_api.rs       305       74         9      222         20             9.01
~is/symbolic_electromagnetism_ffi/handle.rs       304       57        56      191         29            15.18
~symbolic_classical_mechanics_ffi/handle.rs       300       60        64      176         21            11.93
~fi_apis/numerical_elementary_ffi/handle.rs       300       90        35      175          9             5.14
src/ffi_apis/nightly_ffi/bincode_api.rs           299       46        33      220         16             7.27
~i_apis/numerical_matrix_ffi/bincode_api.rs       299       46        33      220         16             7.27
~ctal_geometry_and_chaos_ffi/bincode_api.rs       299       73        21      205         16             7.80
src/numerical/transforms.rs                       297       68       121      108         19            17.59
src/numerical/combinatorics.rs                    297       78        58      161         35            21.74
src/numerical/coordinates.rs                      294       68        75      151         16            10.60
~lic_fractal_geometry_and_chaos_ffi/json.rs       293       73        21      199         16             8.04
~ffi_apis/symbolic_multi_valued_ffi/json.rs       292       67         9      216         11             5.09
~s/symbolic_multi_valued_ffi/bincode_api.rs       291       66         9      216         11             5.09
src/ffi_apis/symbolic_graph_ffi/json.rs           291       73        20      198         12             6.06
src/numerical/solve.rs                            290       63        46      181         23            12.71
src/ffi_apis/symbolic_stats_ffi/handle.rs         288       71        56      161         21            13.04
~/numerical_multi_valued_ffi/bincode_api.rs       287       43        57      187          5             2.67
~ffi_apis/symbolic_finite_field_ffi/json.rs       284       73         9      202         20             9.90
~rc/ffi_apis/numerical_optimize_ffi/json.rs       284       35        42      207          7             3.38
~symbolic_solid_state_physics_ffi/handle.rs       284       59        65      160         20            12.50
src/ffi_apis/constant_ffi/handle.rs               282       79        40      163          7             4.29
~ymbolic_quantum_field_theory_ffi/handle.rs       282       56        64      162         22            13.58
~physics/physics_sim/geodesic_relativity.rs       280       53        52      175          3             1.71
~bolic_stats_information_theory_ffi/json.rs       279       82        77      120         17            14.17
~mbolic_graph_operations_ffi/bincode_api.rs       278       66        64      148          8             5.41
~pis/symbolic_stats_inference_ffi/handle.rs       278       70        70      138         13             9.42
~error_correction_helper_ffi/bincode_api.rs       277       58        12      207         17             8.21
~merical_vector_calculus_ffi/bincode_api.rs       276       32        88      156          6             3.85
~c/physics/physics_sim/gpe_superfluidity.rs       274       55        24      195          3             1.54
src/numerical/functional_analysis.rs              273       59        64      150          9             6.00
src/ffi_apis/compute_cache_ffi/handle.rs          272       65        36      171         20            11.70
src/ffi_apis/symbolic_handles_ffi/json.rs         272       70        27      175          9             5.14
~fi_apis/numerical_polynomial_ffi/handle.rs       271       72        73      126         12             9.52
src/ffi_apis/constant_ffi/json.rs                 269       71        40      158          4             2.53
src/ffi_apis/numerical_sparse_ffi/json.rs         266       45        48      173          8             4.62
~fi_apis/numerical_finite_field_ffi/json.rs       264       52        37      175          9             5.14
~is/symbolic_stats_regression_ffi/handle.rs       263       63        44      156         16            10.26
~umerical_differential_geometry_ffi/json.rs       262       49        49      164          8             4.88
~c/ffi_apis/symbolic_elementary_ffi/json.rs       258       53        78      127          9             7.09
~pis/symbolic_polynomial_ffi/bincode_api.rs       255       58         7      190         18             9.47
src/ffi_apis/numerical_solve_ffi/handle.rs        253       54        71      128         12             9.38
src/numerical/multi_valued.rs                     253       68        27      158         10             6.33
~hysics/physics_sim/fdtd_electrodynamics.rs       252       67        26      159         16            10.06
~c/ffi_apis/symbolic_polynomial_ffi/json.rs       252       58         7      187         18             9.63
~ffi_apis/numerical_interpolate_ffi/json.rs       251       50        49      152          7             4.61
~/symbolic_quantum_field_theory_ffi/json.rs       250       53        48      149          8             5.37
~i_apis/numerical_coordinates_ffi/handle.rs       249       62        44      143         11             7.69
src/symbolic/stats_regression.rs                  248       46        52      150          5             3.33
~rc/ffi_apis/numerical_calculus_ffi/json.rs       248       39        37      172          6             3.49
~fi_apis/symbolic_combinatorics_ffi/json.rs       247       37       124       86          4             4.65
~c_differential_geometry_ffi/bincode_api.rs       243       53         7      183          8             4.37
~erical_differential_geometry_ffi/handle.rs       242       58        34      150         13             8.67
src/ffi_apis/physics_mm_ffi/handle.rs             240       53        43      144         11             7.64
~c/ffi_apis/symbolic_optimize_ffi/handle.rs       240       55        12      173         15             8.67
~pis/symbolic_cas_foundations_ffi/handle.rs       235       48        56      131         15            11.45
~symbolic_differential_geometry_ffi/json.rs       233       53         7      173          8             4.62
src/ffi_apis/numerical_signal_ffi/json.rs         232       31        93      108          3             2.78
~i_apis/symbolic_handles_ffi/bincode_api.rs       232       61         8      163          7             4.29
~s/numerical_complex_analysis_ffi/handle.rs       232       55        25      152         17            11.18
src/ffi_apis/constant_ffi/bincode_api.rs          231       65        28      138          4             2.90
~olic_integral_equations_ffi/bincode_api.rs       231       48         6      177          7             3.95
~tats_information_theory_ffi/bincode_api.rs       230       74        28      128         17            13.28
~l_differential_geometry_ffi/bincode_api.rs       229       44        33      152          8             5.26
~s/symbolic_geometric_algebra_ffi/handle.rs       226       73         8      145         11             7.59
~s/numerical_physics_fea_ffi/bincode_api.rs       226       27        91      108          3             2.78
~_apis/symbolic_number_theory_ffi/handle.rs       225       53        17      155          9             5.81
~s/numerical_physics_cfd_ffi/bincode_api.rs       225       26        94      105          3             2.86
~apis/numerical_optimize_ffi/bincode_api.rs       225       22        40      163          4             2.45
~rc/ffi_apis/numerical_signal_ffi/handle.rs       224       52        41      131          7             5.34
~_apis/numerical_finite_field_ffi/handle.rs       224       56        46      122          9             7.38
~fi_apis/symbolic_matrix_ffi/bincode_api.rs       223       65        24      134         10             7.46
~apis/symbolic_graph_operations_ffi/json.rs       221       59        13      149          8             5.37
~s/symbolic_functional_analysis_ffi/json.rs       221       31       106       84          8             9.52
~i_apis/numerical_convergence_ffi/handle.rs       220       51        71       98          7             7.14
~bolic_calculus_of_variations_ffi/handle.rs       219       38        41      140         17            12.14
~/symbolic_combinatorics_ffi/bincode_api.rs       219       31       122       66          4             6.06
~rical_calculus_of_variations_ffi/handle.rs       219       39        18      162         17            10.49
~ffi_apis/numerical_coordinates_ffi/json.rs       217       39        25      153          8             5.23
src/symbolic/stats.rs                             217       42        59      116          7             6.03
src/ffi_apis/symbolic_matrix_ffi/json.rs          217       66        24      127         10             7.87
src/symbolic/core/mod.rs                          217        4       201       12          0             0.00
~apis/numerical_calculus_ffi/bincode_api.rs       214       35        25      154          6             3.90
~_apis/symbolic_topology_ffi/bincode_api.rs       213       49         8      156         13             8.33
~ffi_apis/numerical_convergence_ffi/json.rs       212       32        37      143          6             4.20
src/compute/cache.rs                              212       39        54      119          2             1.68
~_apis/symbolic_combinatorics_ffi/handle.rs       210       31       122       57          0             0.00
~s/numerical_interpolate_ffi/bincode_api.rs       210       45        33      132          7             5.30
~i_apis/numerical_number_theory_ffi/json.rs       210       39        25      146          7             4.79
src/ffi_apis/physics_fvm_ffi/json.rs              210       25        71      114          5             4.39
src/numerical/signal.rs                           210       40        90       80         14            17.50
src/symbolic/calculus_of_variations.rs            204       22        97       85          0             0.00
~lic_functional_analysis_ffi/bincode_api.rs       204       30       106       68          8            11.76
~/numerical_functional_analysis_ffi/json.rs       201       32        37      132          5             3.79
src/ffi_apis/plugins_ffi/handle.rs                201       34        43      124         13            10.48
~_calculus_of_variations_ffi/bincode_api.rs       198       44         4      150         15            10.00
~i_apis/numerical_signal_ffi/bincode_api.rs       198       27        81       90          3             3.33
~ffi_apis/symbolic_polynomial_ffi/handle.rs       198       53        27      118          1             0.85
src/ffi_apis/numerical_series_ffi/json.rs         197       23        66      108          4             3.70
~ymbolic_calculus_of_variations_ffi/json.rs       197       44         4      149         15            10.07
src/ffi_apis/physics_mtm_ffi/json.rs              196       19        65      112          4             3.57
~umerical_functional_analysis_ffi/handle.rs       196       46        41      109         12            11.01
src/ffi_apis/macros.rs                            195       54        13      128          8             6.25
~i_apis/numerical_physics_fea_ffi/handle.rs       194       46        37      111          3             2.70
src/ffi_apis/symbolic_topology_ffi/json.rs        193       49         8      136         13             9.56
~is/symbolic_integral_equations_ffi/json.rs       193       44         6      143          7             4.90
src/numerical/calculus_of_variations.rs           193       20        92       81          0             0.00
~/symbolic_number_theory_ffi/bincode_api.rs       192       48        17      127          8             6.30
~ffi_apis/symbolic_real_roots_ffi/handle.rs       191       44        11      136         13             9.56
~ymbolic_vector_calculus_ffi/bincode_api.rs       189       42         5      142          4             2.82
~rc/ffi_apis/symbolic_handles_ffi/handle.rs       188       36        65       87          4             4.60
~fi_apis/symbolic_number_theory_ffi/json.rs       186       49        17      120          8             6.67
src/ffi_apis/physics_sm_ffi/json.rs               185       21        73       91          2             2.20
~is/numerical_physics_md_ffi/bincode_api.rs       185       23        61      101          3             2.97
~apis/symbolic_thermodynamics_ffi/handle.rs       184       36        41      107         16            14.95
src/ffi_apis/symbolic_matrix_ffi/handle.rs        183       65        24       94          1             1.06
~rc/ffi_apis/numerical_topology_ffi/json.rs       183       24        66       93          2             2.15
src/ffi_apis/physics_em_ffi/json.rs               181       18        36      127          5             3.94
src/plugins/plugin_c.rs                           181       33        53       95          2             2.11
~rc/ffi_apis/numerical_series_ffi/handle.rs       180       37        26      117         10             8.55
~pis/symbolic_geometric_algebra_ffi/json.rs       179       50         7      122         11             9.02
~pis/symbolic_elementary_ffi/bincode_api.rs       177       40        21      116          9             7.76
~bolic_geometric_algebra_ffi/bincode_api.rs       177       49         7      121         11             9.09
~s/numerical_convergence_ffi/bincode_api.rs       176       30        25      121          4             3.31
~i_apis/numerical_series_ffi/bincode_api.rs       176       20        60       96          4             4.17
~/ffi_apis/numerical_polynomial_ffi/json.rs       175       35        25      115          6             5.22
~/symbolic_poly_factorization_ffi/handle.rs       171       43         6      122         10             8.20
~pis/symbolic_real_roots_ffi/bincode_api.rs       170       39         3      128         14            10.94
src/ffi_apis/symbolic_tensor_ffi/json.rs          169       40        16      113          9             7.96
src/ffi_apis/physics_em_ffi/bincode_api.rs        169       17        33      119          5             4.20
~c/ffi_apis/symbolic_real_roots_ffi/json.rs       168       39         3      126         14            11.11
~/ffi_apis/compute_cache_ffi/bincode_api.rs       168       36         7      125         12             9.60
~apis/symbolic_stats_regression_ffi/json.rs       165       41        38       86          7             8.14
~s/numerical_coordinates_ffi/bincode_api.rs       165       22        17      126          6             4.76
src/ffi_apis/compute_cache_ffi/json.rs            164       36         9      119         18            15.13
~fi_apis/numerical_transforms_ffi/handle.rs       163       37        25      101          8             7.92
~apis/numerical_topology_ffi/bincode_api.rs       162       21        60       81          2             2.47
~i_apis/symbolic_finite_field_ffi/handle.rs       162       50         6      106          9             8.49
~apis/symbolic_rewriting_ffi/bincode_api.rs       162       40         5      117          5             4.27
~cal_functional_analysis_ffi/bincode_api.rs       161       28        25      108          5             4.63
~_apis/symbolic_vector_calculus_ffi/json.rs       158       39         5      114          4             3.51
~_apis/symbolic_stats_inference_ffi/json.rs       158       40        34       84          3             3.57
~rc/ffi_apis/symbolic_rewriting_ffi/json.rs       157       37        17      103          5             4.85
src/compute/mod.rs                                156        1       150        5          0             0.00
~pis/numerical_complex_analysis_ffi/json.rs       155       21        25      109          4             3.67
src/ffi_apis/plugins_ffi/json.rs                  153       24        13      116          5             4.31
~i_apis/numerical_tensor_ffi/bincode_api.rs       151       21         9      121          6             4.96
src/ffi_apis/physics_mtm_ffi/handle.rs            150       38        29       83          5             6.02
src/numerical/series.rs                           150       35        50       65          3             4.62
~fi_apis/numerical_real_roots_ffi/handle.rs       150       32        41       77          8            10.39
~/ffi_apis/numerical_transforms_ffi/json.rs       149       21        61       67          2             2.99
~ffi_apis/symbolic_relativity_ffi/handle.rs       148       32        36       80          6             7.50
src/ffi_apis/physics_rkm_ffi/handle.rs            147       48         5       94         12            12.77
~i_apis/numerical_physics_cfd_ffi/handle.rs       147       33        17       97          0             0.00
src/ffi_apis/symbolic_tensor_ffi/handle.rs        145       42        16       87          4             4.60
~i_apis/physics_sim_schrodinger_ffi/json.rs       143       15        38       90          4             4.44
src/ffi_apis/numerical_ode_ffi/handle.rs          142       31         9      102          9             8.82
~/ffi_apis/numerical_topology_ffi/handle.rs       141       33        25       83          6             7.23
src/ffi_apis/symbolic_proof_ffi/json.rs           140       25        24       91          3             3.30
src/ffi_apis/symbolic_cad_ffi/handle.rs           139       38         6       95         13            13.68
~mbolic_stats_regression_ffi/bincode_api.rs       139       37        14       88          7             7.95
~_apis/symbolic_cas_foundations_ffi/json.rs       138       34         5       99          4             4.04
src/ffi_apis/symbolic_optimize_ffi/json.rs        137       24         3      110          5             4.55
~_apis/symbolic_optimize_ffi/bincode_api.rs       137       23         3      111          5             4.50
~fi_apis/symbolic_tensor_ffi/bincode_api.rs       136       29        12       95          6             6.32
~ymbolic_cas_foundations_ffi/bincode_api.rs       136       33         5       98          4             4.08
src/ffi_apis/physics_mm_ffi/json.rs               136       18        57       61          1             1.64
src/ffi_apis/physics_fdm_ffi/handle.rs            134       35        28       71          7             9.86
~/numerical_finite_field_ffi/bincode_api.rs       134       25        17       92          4             4.35
src/ffi_apis/plugins_ffi/bincode_api.rs           133       18         9      106          5             4.72
~ymbolic_stats_inference_ffi/bincode_api.rs       132       35        12       85          3             3.53
~ffi_apis/symbolic_solve_ffi/bincode_api.rs       132       31        12       89          5             5.62
~ffi_apis/numerical_integrate_ffi/handle.rs       131       22        27       82          6             7.32
src/ffi_apis/physics_bem_ffi/json.rs              131       17        33       81          3             3.70
src/ffi_apis/jit_ffi/json.rs                      130       22        14       94          5             5.32
~i_apis/symbolic_grobner_ffi/bincode_api.rs       130       15        38       77          4             5.19
src/ffi_apis/symbolic_solve_ffi/json.rs           129       32        12       85          5             5.88
~physics_sim_schrodinger_ffi/bincode_api.rs       128       13        36       79          4             5.06
src/ffi_apis/symbolic_stats_ffi/json.rs           128       37         6       85          8             9.41
src/ffi_apis/symbolic_solve_ffi/handle.rs         128       37        12       79          3             3.80
~is/numerical_transforms_ffi/bincode_api.rs       128       18        55       55          2             3.64
~erical_complex_analysis_ffi/bincode_api.rs       127       19        17       91          4             4.40
~rc/ffi_apis/physics_rkm_ffi/bincode_api.rs       127       19        34       74          2             2.70
~c/ffi_apis/numerical_integrate_ffi/json.rs       127       19        23       85          4             4.71
~fi_apis/symbolic_vector_ffi/bincode_api.rs       126       44        16       66          6             9.09
~is/numerical_real_roots_ffi/bincode_api.rs       124       18        28       78          3             3.85
~ffi_apis/symbolic_stats_ffi/bincode_api.rs       122       29         6       87          8             9.20
~i_apis/numerical_sparse_ffi/bincode_api.rs       122       18         9       95          4             4.21
src/ffi_apis/symbolic_grobner_ffi/json.rs         121       16        35       70          4             5.71
~/ffi_apis/numerical_real_roots_ffi/json.rs       120       16        31       73          3             4.11
~apis/physics_sim_schrodinger_ffi/handle.rs       120       21         9       90          6             6.67
~lic_solid_state_physics_ffi/bincode_api.rs       120       25         4       91          3             3.30
~rc/ffi_apis/physics_bem_ffi/bincode_api.rs       119       15        32       72          3             4.17
src/ffi_apis/jit_ffi/handle.rs                    119       26        41       52          4             7.69
src/ffi_apis/physics_em_ffi/handle.rs             118       37         4       77         11            14.29
~is/numerical_elementary_ffi/bincode_api.rs       118       22         9       87          4             4.60
~s/symbolic_solid_state_physics_ffi/json.rs       117       26         4       87          3             3.45
~olic_poly_factorization_ffi/bincode_api.rs       117       26         4       87         11            12.64
~is/symbolic_poly_factorization_ffi/json.rs       117       27         4       86         11            12.79
~is/physics_sim_navier_stokes_ffi/handle.rs       117       21        17       79          4             5.06
src/ffi_apis/symbolic_logic_ffi/handle.rs         116       31        21       64          5             7.81
src/ffi_apis/physics_sim_fdtd_ffi/json.rs         115       14        36       65          4             6.15
src/jit/instructions.rs                           114       14        45       55          0             0.00
~apis/physics_sim_navier_stokes_ffi/json.rs       114       13        37       64          2             3.12
src/ffi_apis/compute_state_ffi/json.rs            114       26         8       80         10            12.50
~_apis/symbolic_discrete_groups_ffi/json.rs       112       20        53       39          1             2.56
~pis/symbolic_discrete_groups_ffi/handle.rs       112       19        53       40          1             2.50
~ffi_apis/symbolic_proof_ffi/bincode_api.rs       112       19        16       77          2             2.60
~lic_classical_mechanics_ffi/bincode_api.rs       112       27         3       82          8             9.76
~ymbolic_discrete_groups_ffi/bincode_api.rs       111       19        53       39          1             2.56
~/ffi_apis/numerical_elementary_ffi/json.rs       110       20        15       75          4             5.33
~s/symbolic_classical_mechanics_ffi/json.rs       109       27         3       79          8            10.13
src/ffi_apis/symbolic_vector_ffi/handle.rs        109       44        16       49          0             0.00
src/ffi_apis/numerical_ode_ffi/json.rs            109       13        34       62          2             3.23
src/ffi_apis/symbolic_logic_ffi/json.rs           107       29         9       69          9            13.04
~rc/ffi_apis/physics_fdm_ffi/bincode_api.rs       107       15        34       58          1             1.72
src/ffi_apis/physics_bem_ffi/handle.rs            107       20        12       75         10            13.33
~fi_apis/numerical_solve_ffi/bincode_api.rs       105       17         9       79          3             3.80
src/ffi_apis/physics_fem_ffi/json.rs              105       11        32       62          2             3.23
~apis/numerical_number_theory_ffi/handle.rs       105       28        17       60          1             1.67
~rc/ffi_apis/symbolic_grobner_ffi/handle.rs       103       14        34       55          2             3.64
~ysics_sim_navier_stokes_ffi/bincode_api.rs       102       11        36       55          2             3.64
src/ffi_apis/numerical_solve_ffi/json.rs          100       15        13       72          3             4.17
src/compute/computation.rs                        100       14        20       66          0             0.00
~ffi_apis/symbolic_logic_ffi/bincode_api.rs        97       28         9       60          8            13.33
~graph_isomorphism_and_coloring_ffi/json.rs        97       20        30       47          3             6.38
~/ffi_apis/numerical_ode_ffi/bincode_api.rs        96       11        29       56          2             3.57
src/ffi_apis/compute_state_ffi/handle.rs           95       26         7       62          6             9.68
src/ffi_apis/physics_sim_gpe_ffi/json.rs           95        8        36       51          2             3.92
~is/numerical_polynomial_ffi/bincode_api.rs        94       18         9       67          3             4.48
~/physics_sim_linear_elasticity_ffi/json.rs        94        8        35       51          2             3.92
~aph_isomorphism_and_coloring_ffi/handle.rs        94       26         8       60          4             6.67
src/ffi_apis/physics_sim_ising_ffi/json.rs         94       13        34       47          1             2.13
src/ffi_apis/physics_fvm_ffi/handle.rs             94       25        19       50          2             4.00
~numerical_number_theory_ffi/bincode_api.rs        93       18         9       66          3             4.55
~fi_apis/numerical_physics_md_ffi/handle.rs        93       29        16       48          4             8.33
src/ffi_apis/physics_cnm_ffi/json.rs               91       12        33       46          1             2.17
~hysics_sim_linear_elasticity_ffi/handle.rs        91       25         4       62          5             8.06
~somorphism_and_coloring_ffi/bincode_api.rs        91       20        24       47          3             6.38
~rc/ffi_apis/physics_fem_ffi/bincode_api.rs        91        9        29       53          2             3.77
src/symbolic/mod.rs                                90        1        21       68          0             0.00
~rc/ffi_apis/physics_fvm_ffi/bincode_api.rs        89       10        33       46          1             2.17
~rc/ffi_apis/physics_mtm_ffi/bincode_api.rs        89        9        30       50          2             4.00
src/ffi_apis/physics_cnm_ffi/handle.rs             87       21        18       48          2             4.17
~c/ffi_apis/physics_sim_ising_ffi/handle.rs        87       16        16       55          1             1.82
~merical_calculus_of_variations_ffi/json.rs        87       16        13       58          2             3.45
~i_apis/physics_sim_fdtd_ffi/bincode_api.rs        86       10        33       43          3             6.98
~mbolic_electromagnetism_ffi/bincode_api.rs        85       17         3       65          2             3.08
~pis/symbolic_quantum_mechanics_ffi/json.rs        85       22         3       60          3             5.00
~fi_apis/symbolic_integration_ffi/handle.rs        84       24         2       58          5             8.62
~_apis/physics_sim_ising_ffi/bincode_api.rs        83       11        33       39          1             2.56
~bolic_quantum_mechanics_ffi/bincode_api.rs        83       21         3       59          3             5.08
src/ffi_apis/physics_sm_ffi/bincode_api.rs         82       10        35       37          1             2.70
~fi_apis/physics_sim_gpe_ffi/bincode_api.rs        82        6        34       42          2             4.76
src/plugins/stable_abi.rs                          82       18        13       51          0             0.00
~pis/numerical_integrate_ffi/bincode_api.rs        81       13         9       59          2             3.39
~/ffi_apis/physics_sim_geodesic_ffi/json.rs        81        9        35       37          1             2.70
~s_sim_linear_elasticity_ffi/bincode_api.rs        81        6        33       42          2             4.76
~_calculus_of_variations_ffi/bincode_api.rs        78       14         9       55          2             3.64
~c/ffi_apis/symbolic_simplify_ffi/handle.rs        78       17        23       38          2             5.26
~i_apis/symbolic_thermodynamics_ffi/json.rs        78       20         3       55          2             3.64
~rc/ffi_apis/physics_cnm_ffi/bincode_api.rs        78       10        30       38          1             2.63
~symbolic_thermodynamics_ffi/bincode_api.rs        77       19         3       55          2             3.64
src/ffi_apis/physics_mm_ffi/bincode_api.rs         75       10        29       36          1             2.78
~/ffi_apis/compute_state_ffi/bincode_api.rs        74       16         8       50          2             4.00
~apis/symbolic_electromagnetism_ffi/json.rs        74       18         3       53          2             3.77
src/ffi_apis/physics_fem_ffi/handle.rs             73       20        11       42          2             4.76
~pis/symbolic_relativity_ffi/bincode_api.rs        70       19         4       47          6            12.77
src/ffi_apis/symbolic_cad_ffi/json.rs              70       16        10       44          4             9.09
~ic_quantum_field_theory_ffi/bincode_api.rs        69       16        16       37          2             5.41
~c/ffi_apis/symbolic_relativity_ffi/json.rs        69       20         4       45          6            13.33
~is/physics_sim_geodesic_ffi/bincode_api.rs        68        7        32       29          1             3.45
~/ffi_apis/symbolic_integration_ffi/json.rs        62       15         2       45          3             6.67
~c/ffi_apis/symbolic_cad_ffi/bincode_api.rs        61       14         3       44          4             9.09
~apis/symbolic_special_functions_ffi/mod.rs        61        2        56        3          0             0.00
~is/symbolic_integration_ffi/bincode_api.rs        61       14         2       45          3             6.67
src/numerical/mod.rs                               56        1         8       47          0             0.00
~rc/ffi_apis/physics_sim_fdtd_ffi/handle.rs        55       11         2       42          2             4.76
src/ffi_apis/physics_sm_ffi/handle.rs              53       14         5       34          3             8.82
src/ffi_apis/symbolic_simplify_ffi/json.rs         53       16         3       34          4            11.76
~_apis/symbolic_simplify_ffi/bincode_api.rs        52       15         3       34          4            11.76
src/physics/mod.rs                                 52        2        38       12          0             0.00
src/ffi_apis/physics_sim_gpe_ffi/handle.rs         51        8         2       41          1             2.44
~fi_apis/physics_sim_geodesic_ffi/handle.rs        51       11         2       38          1             2.63
~c/ffi_apis/symbolic_radicals_ffi/handle.rs        50       16         2       32          2             6.25
src/ffi_apis/symbolic_radicals_ffi/json.rs         50       15         2       33          4            12.12
src/ffi_apis/symbolic_special_ffi/mod.rs           49        1        45        3          0             0.00
~_apis/symbolic_radicals_ffi/bincode_api.rs        49       14         2       33          4            12.12
~/ffi_apis/symbolic_convergence_ffi/json.rs        44       14         5       25          1             4.00
~fi_apis/symbolic_convergence_ffi/handle.rs        44       15         5       24          1             4.17
~is/symbolic_convergence_ffi/bincode_api.rs        43       13         5       25          1             4.00
src/ffi_apis/constant_ffi/mod.rs                   41        2        33        6          0             0.00
~mbolic_unit_unification_ffi/bincode_api.rs        40       10         4       26          3            11.54
~i_apis/symbolic_simplify_dag_ffi/handle.rs        40        9        12       19          1             5.26
src/ffi_apis/symbolic_numeric_ffi/json.rs          40       13         4       23          4            17.39
~i_apis/symbolic_numeric_ffi/bincode_api.rs        39       12         4       23          4            17.39
~apis/symbolic_unit_unification_ffi/json.rs        39       11         4       24          3            12.50
src/compute/state.rs                               36        7         7       22          1             4.55
~is/symbolic_unit_unification_ffi/handle.rs        33        9         4       20          1             5.00
~ffi_apis/symbolic_simplify_dag_ffi/json.rs        30        9         2       19          2            10.53
~s/symbolic_simplify_dag_ffi/bincode_api.rs        29        8         2       19          2            10.53
~_apis/symbolic_error_correction_ffi/mod.rs        27        1        23        3          0             0.00
src/compute/computable.rs                          27        3        15        9          0             0.00
src/output/mod.rs                                  27        1        21        5          0             0.00
src/numerical/pde.rs                               26        4        10       12          0             0.00
~rc/ffi_apis/symbolic_numeric_ffi/handle.rs        25        9         4       12          0             0.00
~rc/ffi_apis/symbolic_transforms_ffi/mod.rs        20        1        16        3          0             0.00
~/ffi_apis/symbolic_cryptography_ffi/mod.rs        20        1        16        3          0             0.00
src/physics/physics_sim/mod.rs                     18        1        10        7          0             0.00
src/ffi_apis/compute_state_ffi/mod.rs              18        2        10        6          0             0.00
~rc/ffi_apis/symbolic_elementary_ffi/mod.rs        18        2        10        6          0             0.00
src/ffi_apis/compute_cache_ffi/mod.rs              18        2        10        6          0             0.00
~/ffi_apis/symbolic_simplify_dag_ffi/mod.rs        17        2         9        6          0             0.00
src/ffi_apis/symbolic_calculus_ffi/mod.rs          17        2         9        6          0             0.00
src/ffi_apis/symbolic_simplify_ffi/mod.rs          17        2         9        6          0             0.00
src/ffi_apis/symbolic_rewriting_ffi/mod.rs         17        2         9        6          0             0.00
~i_apis/symbolic_cas_foundations_ffi/mod.rs        17        2         9        6          0             0.00
src/ffi_apis/symbolic_handles_ffi/mod.rs           17        2         9        6          0             0.00
~pis/symbolic_integral_equations_ffi/mod.rs        17        2         9        6          0             0.00
src/ffi_apis/symbolic_pde_ffi/mod.rs               17        2         9        6          0             0.00
src/ffi_apis/symbolic_ode_ffi/mod.rs               17        2         9        6          0             0.00
~i_apis/symbolic_vector_calculus_ffi/mod.rs        17        2         9        6          0             0.00
~c/ffi_apis/symbolic_integration_ffi/mod.rs        16        1        12        3          0             0.00
src/jit/mod.rs                                     16        3         2       11          0             0.00
~apis/symbolic_computer_graphics_ffi/mod.rs        16        1        12        3          0             0.00
~ymbolic_error_correction_helper_ffi/mod.rs        16        1        12        3          0             0.00
src/plugins/mod.rs                                 15        1        11        3          0             0.00
src/ffi_apis/symbolic_optimize_ffi/mod.rs          13        1         9        3          0             0.00
src/ffi_apis/physics_sm_ffi/mod.rs                 13        1         9        3          0             0.00
~apis/symbolic_quantum_mechanics_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/symbolic_topology_ffi/mod.rs          13        1         9        3          0             0.00
src/ffi_apis/numerical_stats_ffi/mod.rs            13        1         9        3          0             0.00
src/ffi_apis/symbolic_tensor_ffi/mod.rs            13        1         9        3          0             0.00
~_apis/symbolic_stats_regression_ffi/mod.rs        13        1         9        3          0             0.00
~apis/symbolic_stats_probability_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/numerical_tensor_ffi/mod.rs           13        1         9        3          0             0.00
~mbolic_stats_information_theory_ffi/mod.rs        13        1         9        3          0             0.00
~i_apis/symbolic_stats_inference_ffi/mod.rs        13        1         9        3          0             0.00
~c/ffi_apis/numerical_real_roots_ffi/mod.rs        13        1         9        3          0             0.00
~symbolic_calculus_of_variations_ffi/mod.rs        13        1         9        3          0             0.00
~is/symbolic_solid_state_physics_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/symbolic_stats_ffi/mod.rs             13        1         9        3          0             0.00
src/ffi_apis/numerical_calculus_ffi/mod.rs         13        1         9        3          0             0.00
~fi_apis/symbolic_thermodynamics_ffi/mod.rs        13        1         9        3          0             0.00
~c/ffi_apis/numerical_transforms_ffi/mod.rs        13        1         9        3          0             0.00
~c/ffi_apis/numerical_physics_md_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/numerical_physics_ffi/mod.rs          13        1         9        3          0             0.00
src/ffi_apis/numerical_topology_ffi/mod.rs         13        1         9        3          0             0.00
src/ffi_apis/numerical_vector_ffi/mod.rs           13        1         9        3          0             0.00
~_apis/numerical_vector_calculus_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/symbolic_solve_ffi/mod.rs             13        1         9        3          0             0.00
src/ffi_apis/numerical_matrix_ffi/mod.rs           13        1         9        3          0             0.00
~umerical_calculus_of_variations_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/numerical_special_ffi/mod.rs          13        1         9        3          0             0.00
src/ffi_apis/numerical_solve_ffi/mod.rs            13        1         9        3          0             0.00
~apis/numerical_complex_analysis_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/symbolic_series_ffi/mod.rs            13        1         9        3          0             0.00
src/ffi_apis/physics_cnm_ffi/mod.rs                13        1         9        3          0             0.00
~/ffi_apis/symbolic_multi_valued_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/numerical_series_ffi/mod.rs           13        1         9        3          0             0.00
~pis/symbolic_poly_factorization_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/physics_bem_ffi/mod.rs                13        1         9        3          0             0.00
~rc/ffi_apis/symbolic_real_roots_ffi/mod.rs        13        1         9        3          0             0.00
~pis/numerical_computer_graphics_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/physics_fdm_ffi/mod.rs                13        1         9        3          0             0.00
src/ffi_apis/symbolic_radicals_ffi/mod.rs          13        1         9        3          0             0.00
~/ffi_apis/numerical_physics_fea_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/physics_em_ffi/mod.rs                 13        1         9        3          0             0.00
~ffi_apis/numerical_multi_valued_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/physics_fem_ffi/mod.rs                13        1         9        3          0             0.00
~rc/ffi_apis/symbolic_polynomial_ffi/mod.rs        13        1         9        3          0             0.00
~rc/ffi_apis/symbolic_relativity_ffi/mod.rs        13        1         9        3          0             0.00
~/ffi_apis/numerical_convergence_ffi/mod.rs        13        1         9        3          0             0.00
~_apis/symbolic_unit_unification_ffi/mod.rs        13        1         9        3          0             0.00
~_apis/symbolic_electromagnetism_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/numerical_sparse_ffi/mod.rs           13        1         9        3          0             0.00
src/ffi_apis/symbolic_numeric_ffi/mod.rs           13        1         9        3          0             0.00
src/ffi_apis/numerical_optimize_ffi/mod.rs         13        1         9        3          0             0.00
~fi_apis/numerical_combinatorics_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/physics_mm_ffi/mod.rs                 13        1         9        3          0             0.00
src/ffi_apis/physics_fvm_ffi/mod.rs                13        1         9        3          0             0.00
~c/ffi_apis/numerical_polynomial_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/physics_mtm_ffi/mod.rs                13        1         9        3          0             0.00
~/ffi_apis/numerical_interpolate_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/symbolic_matrix_ffi/mod.rs            13        1         9        3          0             0.00
src/ffi_apis/physics_rkm_ffi/mod.rs                13        1         9        3          0             0.00
~pis/numerical_geometric_algebra_ffi/mod.rs        13        1         9        3          0             0.00
~/ffi_apis/numerical_physics_cfd_ffi/mod.rs        13        1         9        3          0             0.00
~/ffi_apis/numerical_coordinates_ffi/mod.rs        13        1         9        3          0             0.00
~numerical_differential_geometry_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/nightly_ffi/mod.rs                    13        1         9        3          0             0.00
~s/symbolic_quantum_field_theory_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/symbolic_vector_ffi/mod.rs            13        1         9        3          0             0.00
~c/ffi_apis/symbolic_coordinates_ffi/mod.rs        13        1         9        3          0             0.00
~ffi_apis/numerical_finite_field_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/symbolic_proof_ffi/mod.rs             13        1         9        3          0             0.00
src/ffi_apis/numerical_graph_ffi/mod.rs            13        1         9        3          0             0.00
~fi_apis/numerical_number_theory_ffi/mod.rs        13        1         9        3          0             0.00
~apis/numerical_error_correction_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/physics_sim_gpe_ffi/mod.rs            13        1         9        3          0             0.00
~ical_fractal_geometry_and_chaos_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/physics_sim_fdtd_ffi/mod.rs           13        2         8        3          0             0.00
src/ffi_apis/numerical_ode_ffi/mod.rs              13        1         9        3          0             0.00
~c/ffi_apis/symbolic_convergence_ffi/mod.rs        13        1         9        3          0             0.00
~c/ffi_apis/physics_sim_geodesic_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/plugins_ffi/mod.rs                    13        1         9        3          0             0.00
~fi_apis/physics_sim_schrodinger_ffi/mod.rs        13        1         9        3          0             0.00
~s/physics_sim_linear_elasticity_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/numerical_signal_ffi/mod.rs           13        1         9        3          0             0.00
~is/symbolic_classical_mechanics_ffi/mod.rs        13        1         9        3          0             0.00
~_apis/physics_sim_navier_stokes_ffi/mod.rs        13        1         9        3          0             0.00
~ffi_apis/symbolic_number_theory_ffi/mod.rs        13        1         9        3          0             0.00
src/ffi_apis/physics_sim_ising_ffi/mod.rs          13        1         9        3          0             0.00
src/ffi_apis/symbolic_logic_ffi/mod.rs             12        1         8        3          0             0.00
~olic_fractal_geometry_and_chaos_ffi/mod.rs        12        1         8        3          0             0.00
~_apis/symbolic_complex_analysis_ffi/mod.rs        12        1         8        3          0             0.00
~/ffi_apis/symbolic_finite_field_ffi/mod.rs        12        1         8        3          0             0.00
~is/symbolic_functional_analysis_ffi/mod.rs        12        1         8        3          0             0.00
~/symbolic_differential_geometry_ffi/mod.rs        12        1         8        3          0             0.00
~_graph_isomorphism_and_coloring_ffi/mod.rs        12        1         8        3          0             0.00
src/ffi_apis/symbolic_graph_ffi/mod.rs             12        1         8        3          0             0.00
~_apis/symbolic_graph_operations_ffi/mod.rs        12        1         8        3          0             0.00
src/ffi_apis/symbolic_cad_ffi/mod.rs               12        1         8        3          0             0.00
~i_apis/symbolic_discrete_groups_ffi/mod.rs        12        1         8        3          0             0.00
~_apis/symbolic_graph_algorithms_ffi/mod.rs        11        1         7        3          0             0.00
src/ffi_apis/symbolic_grobner_ffi/mod.rs           11        1         7        3          0             0.00
~ffi_apis/symbolic_combinatorics_ffi/mod.rs        11        1         7        3          0             0.00
~apis/symbolic_geometric_algebra_ffi/mod.rs        11        1         7        3          0             0.00
~rc/ffi_apis/symbolic_lie_groups_ffi/mod.rs        11        1         7        3          0             0.00
src/ffi_apis/jit_ffi/mod.rs                        11        1         8        2          0             0.00
~/ffi_apis/symbolic_group_theory_ffi/mod.rs        11        1         7        3          0             0.00
~s/numerical_functional_analysis_ffi/mod.rs        11        1         7        3          0             0.00
~rc/ffi_apis/numerical_integrate_ffi/mod.rs        11        1         7        3          0             0.00
~c/ffi_apis/numerical_elementary_ffi/mod.rs        10        1         6        3          0             0.00
src/verification/mod.rs                             8        1         6        1          0             0.00
src/ffi_blindings/mod.rs                            4        1         2        1          0             0.00
src/input/mod.rs                                    3        0         2        1          0             0.00
src/nightly/mod.rs                                  1        0         0        1          0             0.00
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Total                                   686    975987    76690     53791   845506      36303          3578.78
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Unique Lines of Code (ULOC)                    158734
DRYness %                                        0.16
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
Estimated Cost to Develop (semi-detached) $128,205,315
Estimated Schedule Effort (semi-detached) 65.72 months
Estimated People Required (semi-detached) 173.30
Processed 37382330 bytes, 37.382 megabytes (SI)
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
```
