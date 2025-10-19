# Code Stastics (SCC CST 2025-10-19)

```
───────────────────────────────────────────────────────────────────────────────
Language                 Files     Lines   Blanks  Comments     Code Complexity
───────────────────────────────────────────────────────────────────────────────
Rust                       427     56379     2268     15824    38287       4796
Markdown                     4      1233      182         0     1051          0
TOML                         3       128       14         1      113          2
YAML                         2       306       53        17      236          0
Batch                        1         3        0         0        3          0
Fortran Modern               1       167       34        33      100          6
License                      1        73       32         0       41          0
Shell                        1         5        1         1        3          0
───────────────────────────────────────────────────────────────────────────────
Total                      440     58294     2584     15876    39834       4804
───────────────────────────────────────────────────────────────────────────────
Estimated Cost to Develop (organic) $1,293,784
Estimated Schedule Effort (organic) 15.17 months
Estimated People Required (organic) 7.58
───────────────────────────────────────────────────────────────────────────────
Processed 2291136 bytes, 2.291 megabytes (SI)
───────────────────────────────────────────────────────────────────────────────

───────────────────────────────────────────────────────────────────────────────
Language                 Files     Lines   Blanks  Comments     Code Complexity
───────────────────────────────────────────────────────────────────────────────
Rust                       427     56379     2268     15824    38287       4796
───────────────────────────────────────────────────────────────────────────────
src\symbolic\core.rs                3490      147       267     3076        137
src\ffi_apis\ffi_api.rs             2940        3       130     2807        375
src\symbolic\calculus.rs            1920       22       194     1704        365
~\symbolic\graph_algorithms.rs      1446        0       297     1149        289
src\symbolic\simplify.rs            1256       12        58     1186        176
src\symbolic\simplify_dag.rs        1011       58       155      798        135
src\prelude.rs                       928        0        16      912          0
src\symbolic\polynomial.rs           910        0       208      702        130
src\symbolic\matrix.rs               888        0       213      675        187
src\symbolic\solve.rs                794        0        89      705        135
src\symbolic\ode.rs                  784        0       143      641        112
src\numerical\testing.rs             769        0       131      638        123
src\symbolic\pde.rs                  735        0       123      612        101
src\symbolic\coordinates.rs          658        0       168      490         64
src\numerical\matrix.rs              643        0       139      504        110
~symbolic\computer_graphics.rs       634        0       164      470          7
src\symbolic\number_theory.rs        541        0       148      393        108
src\symbolic\logic.rs                535        0        74      461         88
~symbolic\special_functions.rs       524        0       209      315         64
src\symbolic\transforms.rs           499        0       181      318         39
src\symbolic\tensor.rs               495        0       173      322         63
src\symbolic\topology.rs             482        0       184      298         55
~ic\error_correction_helper.rs       481        0       143      338         60
~ymbolic\poly_factorization.rs       476        0        94      382         58
src\symbolic\combinatorics.rs        457        0       166      291         47
src\numerical\optimize.rs            452        2        33      417         35
src\numerical\physics.rs             444        0       113      331         64
src\symbolic\integration.rs          402        0        44      358         35
src\symbolic\finite_field.rs         401        0        75      326         51
~symbolic\stats_probability.rs       387        0       163      224          0
src\physics\physics_fvm.rs           386        0        66      320         38
~c\symbolic\cas_foundations.rs       386        4        40      342         57
~ymbolic\integral_equations.rs       357        0       120      237          7
src\symbolic\series.rs               355        0       124      231         31
src\numerical\stats.rs               346        0       114      232         12
src\output\pretty_print.rs           330       11         8      311         36
src\physics\physics_mtm.rs           301        0        23      278         49
src\output\plotting.rs               283       10         8      265         15
src\physics\physics_fem.rs           280        0        15      265         54
src\physics\physics_rkm.rs           279        0        49      230         14
src\physics\physics_sm.rs            272        0        12      260         12
src\output\latex.rs                  268       16         8      244         14
~\symbolic\unit_unification.rs       255        0         4      251         26
~symbolic\geometric_algebra.rs       253        0        85      168         38
src\symbolic\vector.rs               252        0       122      130          3
src\symbolic\elementary.rs           247        0        52      195         11
src\lib.rs                           238        5       148       85          1
~ic\lie_groups_and_algebras.rs       238        0        76      162          7
src\symbolic\proof.rs                236        0        81      155         33
src\symbolic\rewriting.rs            231        0        39      192         46
src\symbolic\graph.rs                230        0       104      126         15
src\physics\physics_em.rs            230        0        85      145         10
src\numerical\interpolate.rs         227        0        56      171         41
src\symbolic\cad.rs                  224        0        10      214         35
src\numerical\finite_field.rs        215        0        71      144         27
~olic\differential_geometry.rs       213        0        67      146         20
src\symbolic\grobner.rs              212        0        38      174         19
~\symbolic\error_correction.rs       208        0        59      149         27
src\physics\physics_bem.rs           205        0        34      171         17
~numerical\error_correction.rs       204        0        36      168         27
src\numerical\sparse.rs              202        0        69      133         17
src\physics\physics_mm.rs            200        0         8      192         21
~mbolic\functional_analysis.rs       198        0       106       92          2
src\numerical\polynomial.rs          197        0        44      153         30
~c\stats_information_theory.rs       197        0        90      107         12
src\symbolic\real_roots.rs           194        0        49      145         28
src\physics\physics_fdm.rs           188        0        36      152         26
~c\symbolic\discrete_groups.rs       182        0        42      140         25
~umerical\computer_graphics.rs       178        0        72      106          7
src\numerical\coordinates.rs         173        0        63      110         16
src\plugins\manager.rs               172        0        22      150         15
src\physics\physics_cnm.rs           164        0        17      147         17
~rc\symbolic\thermodynamics.rs       160        0       100       60          0
~h_isomorphism_and_coloring.rs       159        0        42      117         29
~c\symbolic\vector_calculus.rs       157        0        62       95          0
src\symbolic\convergence.rs          157        0        47      110         34
~sics_sim\linear_elasticity.rs       156        0        22      134         20
~\symbolic\complex_analysis.rs       151        0        59       92          8
~numerical\complex_analysis.rs       148        0        61       87          3
src\output\typst.rs                  146       12         3      131         11
src\symbolic\group_theory.rs         146        0        61       85         14
src\symbolic\optimize.rs             145        0        45      100         18
src\numerical\integrate.rs           145        0        51       94         11
~mbolic\classical_mechanics.rs       145        0        93       52          0
~symbolic\quantum_mechanics.rs       143        0        75       68          0
~\numerical\vector_calculus.rs       134        0        44       90          7
~umerical\geometric_algebra.rs       133        0        22      111          4
src\symbolic\cryptography.rs         133        0        56       77         13
src\numerical\solve.rs               131        0        34       97         20
src\numerical\real_roots.rs          128        0        26      102         43
src\numerical\physics_cfd.rs         128        0        44       84         11
~cs_sim\schrodinger_quantum.rs       127        0        10      117         14
src\numerical\transforms.rs          125        0        41       84         19
src\symbolic\relativity.rs           125        0        67       58          0
~\symbolic\graph_operations.rs       124        0         5      119         41
src\numerical\tensor.rs              122        0        37       85         15
src\symbolic\numeric.rs              119        0        24       95         32
src\symbolic\stats.rs                117        0        59       58          7
~\symbolic\stats_regression.rs       116        0        45       71          5
~cs_sim\geodesic_relativity.rs       116        0        26       90          3
src\numerical\vector.rs              115        0        68       47         14
~rc\numerical\number_theory.rs       114        0        43       71         31
src\output\io.rs                     112        0        43       69          9
src\numerical\convergence.rs         111        0        41       70          9
src\numerical\graph.rs               111        0        44       67          7
src\symbolic\multi_valued.rs         108        0        65       43          0
src\numerical\topology.rs            108        0        25       83         19
~cs_sim\navier_stokes_fluid.rs       106        0         4      102         17
~\symbolic\electromagnetism.rs       102        0        52       50          0
~bolic\quantum_field_theory.rs       100        0        60       40          2
tests\symbolic\simplify.rs            98        9        25       64          0
src\numerical\physics_fea.rs          98        0        41       57         10
~lic\calculus_of_variations.rs        98        0        69       29          0
~rc\numerical\combinatorics.rs        97        0        40       57         15
tests\symbolic\calculus.rs            96       10        25       61          0
src\numerical\elementary.rs           93        0        29       64          1
~erical\functional_analysis.rs        91        0        47       44          4
src\numerical\physics_md.rs           90        0        35       55          4
src\numerical\ode.rs                  89        0        20       69          4
src\symbolic\mod.rs                   88        0        20       68          0
~sics_sim\ising_statistical.rs        86        0        10       76         11
~sics_sim\gpe_superfluidity.rs        85        0         6       79          5
src\symbolic\special.rs               82        0        61       21          0
src\numerical\signal.rs               82        0        26       56          4
~s_sim\fdtd_electrodynamics.rs        81        0        10       71         16
~c\symbolic\stats_inference.rs        81        0        19       62          0
~ical\differential_geometry.rs        79        0        20       59         12
tests\numerical\calculus.rs           76       14        11       51          3
~fractal_geometry_and_chaos.rs        75        0        31       44          5
~ests\numerical\interpolate.rs        72       11        13       48          7
~mbolic\solid_state_physics.rs        65        0        40       25          0
src\plugins\mod.rs                    64        0        27       37          2
src\symbolic\radicals.rs              59        0        16       43         10
src\numerical\calculus.rs             57        0        29       28          4
src\numerical\multi_valued.rs         53        0        19       34          5
src\numerical\series.rs               51        0        15       36          2
~fractal_geometry_and_chaos.rs        51        0        30       21          1
tests\mod.rs                          50        6        42        2          0
~sts\symbolic\combinatorics.rs        48        6        25       17          0
tests\numerical\solve.rs              48        6        25       17          0
tests\output\pretty_print.rs          48        6        25       17          0
~numerical\complex_analysis.rs        48        6        25       17          0
tests\physics\mod.rs                  48        6        25       17          0
~cal\calculus_of_variations.rs        48        6        25       17          0
tests\numerical\elementary.rs         48        6        25       17          0
tests\output\typst.rs                 48        6        25       17          0
~ical\differential_geometry.rs        48        6        25       17          0
~numerical\error_correction.rs        48        6        25       17          0
~sts\numerical\finite_field.rs        48        6        25       17          0
tests\physics\physics_fem.rs          48        6        25       17          0
tests\physics\physics_bem.rs          48        6        25       17          0
tests\physics\physics_em.rs           48        6        25       17          0
~umerical\geometric_algebra.rs        48        6        25       17          0
~erical\functional_analysis.rs        48        6        25       17          0
tests\physics\physics_cnm.rs          48        6        25       17          0
tests\physics\physics_fvm.rs          48        6        25       17          0
tests\physics\physics_mm.rs           48        6        25       17          0
tests\physics\physics_fdm.rs          48        6        25       17          0
tests\physics\physics_mtm.rs          48        6        25       17          0
~fractal_geometry_and_chaos.rs        48        6        25       17          0
tests\physics\physics_rkm.rs          48        6        25       17          0
tests\symbolic\logic.rs               48        6        25       17          0
~sts\symbolic\number_theory.rs        48        6        25       17          0
tests\numerical\topology.rs           48        6        25       17          0
~ic\lie_groups_and_algebras.rs        48        6        25       17          0
~ests\symbolic\multi_valued.rs        48        6        25       17          0
tests\numerical\integrate.rs          48        6        25       17          0
tests\physics\physics_sm.rs           48        6        25       17          0
tests\numerical\graph.rs              48        6        25       17          0
tests\numerical\matrix.rs             48        6        25       17          0
~sts\numerical\multi_valued.rs        48        6        25       17          0
~ts\numerical\number_theory.rs        48        6        25       17          0
tests\numerical\optimize.rs           48        6        25       17          0
tests\numerical\pde.rs                48        6        25       17          0
tests\numerical\physics.rs            48        6        25       17          0
~s_sim\fdtd_electrodynamics.rs        48        6        25       17          0
~ests\numerical\convergence.rs        48        6        25       17          0
~umerical\computer_graphics.rs        48        6        25       17          0
~cs_sim\geodesic_relativity.rs        48        6        25       17          0
~ts\numerical\combinatorics.rs        48        6        25       17          0
tests\output\mod.rs                   48        6        25       17          0
~sics_sim\gpe_superfluidity.rs        48        6        25       17          0
~sics_sim\ising_statistical.rs        48        6        25       17          0
~sics_sim\linear_elasticity.rs        48        6        25       17          0
~ests\numerical\coordinates.rs        48        6        25       17          0
tests\output\plotting.rs              48        6        25       17          0
~cs_sim\schrodinger_quantum.rs        48        6        25       17          0
~ests\numerical\physics_cfd.rs        48        6        25       17          0
~ests\numerical\physics_fea.rs        48        6        25       17          0
tests\symbolic\integration.rs         48        6        25       17          0
tests\numerical\polynomial.rs         48        6        25       17          0
~cs_sim\navier_stokes_fluid.rs        48        6        25       17          0
tests\numerical\real_roots.rs         48        6        25       17          0
tests\symbolic\vector.rs              48        6        25       17          0
~ts\physics\physics_sim\mod.rs        48        6        25       17          0
tests\symbolic\stats.rs               48        6        25       17          0
tests\numerical\series.rs             48        6        25       17          0
tests\output\latex.rs                 48        6        25       17          0
tests\symbolic\ode.rs                 48        6        25       17          0
tests\symbolic\transforms.rs          48        6        25       17          0
src\physics\mod.rs                    48        0        37       11          0
tests\output\io.rs                    48        6        25       17          0
tests\numerical\physics_md.rs         48        6        25       17          0
~c\stats_information_theory.rs        48        6        25       17          0
~s\symbolic\vector_calculus.rs        48        6        25       17          0
tests\symbolic\matrix.rs              48        6        25       17          0
~lic\calculus_of_variations.rs        48        6        25       17          0
~symbolic\stats_probability.rs        48        6        25       17          0
~s\symbolic\cas_foundations.rs        48        6        25       17          0
~mbolic\classical_mechanics.rs        48        6        25       17          0
tests\numerical\signal.rs             48        6        25       17          0
~ts\symbolic\thermodynamics.rs        48        6        25       17          0
~ymbolic\integral_equations.rs        48        6        25       17          0
tests\symbolic\topology.rs            48        6        25       17          0
~\numerical\vector_calculus.rs        48        6        25       17          0
~\symbolic\complex_analysis.rs        48        6        25       17          0
tests\symbolic\tensor.rs              48        6        25       17          0
tests\numerical\sparse.rs             48        6        25       17          0
~\symbolic\stats_regression.rs        48        6        25       17          0
tests\numerical\vector.rs             48        6        25       17          0
~symbolic\computer_graphics.rs        48        6        25       17          0
~s\symbolic\stats_inference.rs        48        6        25       17          0
tests\symbolic\convergence.rs         48        6        25       17          0
~symbolic\special_functions.rs        48        6        25       17          0
tests\symbolic\special.rs             48        6        25       17          0
tests\symbolic\solve.rs               48        6        25       17          0
~mbolic\solid_state_physics.rs        48        6        25       17          0
tests\symbolic\coordinates.rs         48        6        25       17          0
tests\numerical\special.rs            48        6        25       17          0
tests\symbolic\series.rs              48        6        25       17          0
~ests\symbolic\cryptography.rs        48        6        25       17          0
tests\numerical\transforms.rs         48        6        25       17          0
~olic\differential_geometry.rs        48        6        25       17          0
~s\symbolic\discrete_groups.rs        48        6        25       17          0
tests\symbolic\rewriting.rs           48        6        25       17          0
tests\symbolic\elementary.rs          48        6        25       17          0
~ic\error_correction_helper.rs        48        6        25       17          0
tests\numerical\tensor.rs             48        6        25       17          0
~\symbolic\error_correction.rs        48        6        25       17          0
~ests\symbolic\finite_field.rs        48        6        25       17          0
tests\numerical\stats.rs              48        6        25       17          0
~\symbolic\electromagnetism.rs        48        6        25       17          0
tests\symbolic\relativity.rs          48        6        25       17          0
~fractal_geometry_and_chaos.rs        48        6        25       17          0
~mbolic\functional_analysis.rs        48        6        25       17          0
tests\symbolic\cad.rs                 48        6        25       17          0
tests\symbolic\polynomial.rs          48        6        25       17          0
~symbolic\geometric_algebra.rs        48        6        25       17          0
tests\symbolic\graph.rs               48        6        25       17          0
tests\symbolic\real_roots.rs          48        6        25       17          0
tests\symbolic\radicals.rs            48        6        25       17          0
~bolic\quantum_field_theory.rs        48        6        25       17          0
~\symbolic\graph_algorithms.rs        48        6        25       17          0
~symbolic\quantum_mechanics.rs        48        6        25       17          0
tests\symbolic\proof.rs               48        6        25       17          0
~\symbolic\graph_operations.rs        48        6        25       17          0
~ymbolic\poly_factorization.rs        48        6        25       17          0
tests\numerical\testing.rs            48        6        25       17          0
~h_isomorphism_and_coloring.rs        48        6        25       17          0
tests\symbolic\grobner.rs             48        6        25       17          0
tests\symbolic\pde.rs                 48        6        25       17          0
tests\symbolic\numeric.rs             48        6        25       17          0
~ests\symbolic\group_theory.rs        48        6        25       17          0
tests\symbolic\optimize.rs            48        6        25       17          0
tests\lib.rs                          47        5        42        0          0
src\symbolic\handles.rs               47        0        14       33          0
~ins\example_plugin\src\lib.rs        47        7         9       31          1
tests\ffi_blindings\mod.rs            47        5        42        0          0
tests\ffi_apis\ffi_api.rs             47        5        42        0          0
tests\ffi_apis\mod.rs                 47        5        42        0          0
tests\prelude.rs                      47        5        42        0          0
src\numerical\mod.rs                  45        0         4       41          0
~nches\numerical\real_roots.rs        43        8        23       12          0
~sics_sim\gpe_superfluidity.rs        43        8        23       12          0
benches\numerical\matrix.rs           43        8        23       12          0
benches\numerical\tensor.rs           43        8        23       12          0
~fractal_geometry_and_chaos.rs        43        8        23       12          0
~mbolic\functional_analysis.rs        43        8        23       12          0
~ic\error_correction_helper.rs        43        8        23       12          0
~\symbolic\error_correction.rs        43        8        23       12          0
~symbolic\geometric_algebra.rs        43        8        23       12          0
benches\symbolic\graph.rs             43        8        23       12          0
~\symbolic\graph_algorithms.rs        43        8        23       12          0
~enches\symbolic\elementary.rs        43        8        23       12          0
~\symbolic\electromagnetism.rs        43        8        23       12          0
~s\symbolic\discrete_groups.rs        43        8        23       12          0
~h_isomorphism_and_coloring.rs        43        8        23       12          0
~olic\differential_geometry.rs        43        8        23       12          0
benches\symbolic\grobner.rs           43        8        23       12          0
benches\numerical\mod.rs              43        8        23       12          0
~symbolic\computer_graphics.rs        43        8        23       12          0
~nches\symbolic\convergence.rs        43        8        23       12          0
~es\numerical\number_theory.rs        43        8        23       12          0
~\symbolic\graph_operations.rs        43        8        23       12          0
~ymbolic\integral_equations.rs        43        8        23       12          0
~nches\symbolic\integration.rs        43        8        23       12          0
~ches\symbolic\cryptography.rs        43        8        23       12          0
~hes\symbolic\combinatorics.rs        43        8        23       12          0
benches\symbolic\core.rs              43        8        23       12          0
~ches\symbolic\group_theory.rs        43        8        23       12          0
~nches\symbolic\coordinates.rs        43        8        23       12          0
benches\numerical\optimize.rs         43        8        23       12          0
~ic\lie_groups_and_algebras.rs        43        8        23       12          0
benches\symbolic\logic.rs             43        8        23       12          0
~mbolic\classical_mechanics.rs        43        8        23       12          0
benches\symbolic\matrix.rs            43        8        23       12          0
~\symbolic\complex_analysis.rs        43        8        23       12          0
benches\symbolic\mod.rs               43        8        23       12          0
benches\numerical\solve.rs            43        8        23       12          0
~ches\symbolic\multi_valued.rs        43        8        23       12          0
benches\symbolic\cad.rs               43        8        23       12          0
benches\numerical\pde.rs              43        8        23       12          0
~hes\symbolic\number_theory.rs        43        8        23       12          0
benches\numerical\topology.rs         43        8        23       12          0
~s\symbolic\cas_foundations.rs        43        8        23       12          0
benches\symbolic\numeric.rs           43        8        23       12          0
~hes\numerical\multi_valued.rs        43        8        23       12          0
benches\lib.rs                        43        8        23       12          0
benches\symbolic\ode.rs               43        8        23       12          0
~ches\numerical\interpolate.rs        43        8        23       12          0
~lic\calculus_of_variations.rs        43        8        23       12          0
benches\symbolic\pde.rs               43        8        23       12          0
benches\numerical\graph.rs            43        8        23       12          0
benches\symbolic\calculus.rs          43        8        23       12          0
benches\numerical\physics.rs          43        8        23       12          0
benches\symbolic\optimize.rs          43        8        23       12          0
~ches\numerical\physics_fea.rs        43        8        23       12          0
~nches\numerical\physics_md.rs        43        8        23       12          0
~ches\numerical\physics_cfd.rs        43        8        23       12          0
~nches\numerical\polynomial.rs        43        8        23       12          0
benches\symbolic\proof.rs             43        8        23       12          0
~enches\numerical\integrate.rs        43        8        23       12          0
~cs_sim\schrodinger_quantum.rs        43        8        23       12          0
~enches\symbolic\polynomial.rs        43        8        23       12          0
~\numerical\vector_calculus.rs        43        8        23       12          0
benches\numerical\ode.rs              43        8        23       12          0
benches\plugins\mod.rs                43        8        23       12          0
~fractal_geometry_and_chaos.rs        43        8        23       12          0
~bolic\quantum_field_theory.rs        43        8        23       12          0
~numerical\error_correction.rs        43        8        23       12          0
benches\prelude.rs                    43        8        23       12          0
benches\numerical\stats.rs            43        8        23       12          0
~cs_sim\navier_stokes_fluid.rs        43        8        23       12          0
~symbolic\quantum_mechanics.rs        43        8        23       12          0
~sics_sim\linear_elasticity.rs        43        8        23       12          0
benches\symbolic\rewriting.rs         43        8        23       12          0
~enches\symbolic\relativity.rs        43        8        23       12          0
benches\output\latex.rs               43        8        23       12          0
~ymbolic\poly_factorization.rs        43        8        23       12          0
~ches\symbolic\finite_field.rs        43        8        23       12          0
benches\symbolic\series.rs            43        8        23       12          0
benches\symbolic\simplify.rs          43        8        23       12          0
~mbolic\solid_state_physics.rs        43        8        23       12          0
benches\symbolic\solve.rs             43        8        23       12          0
benches\numerical\signal.rs           43        8        23       12          0
benches\symbolic\radicals.rs          43        8        23       12          0
~symbolic\special_functions.rs        43        8        23       12          0
~umerical\geometric_algebra.rs        43        8        23       12          0
benches\symbolic\special.rs           43        8        23       12          0
~s_sim\fdtd_electrodynamics.rs        43        8        23       12          0
~es\physics\physics_sim\mod.rs        43        8        23       12          0
benches\numerical\special.rs          43        8        23       12          0
benches\symbolic\stats.rs             43        8        23       12          0
~s\symbolic\stats_inference.rs        43        8        23       12          0
~enches\symbolic\real_roots.rs        43        8        23       12          0
~c\stats_information_theory.rs        43        8        23       12          0
~symbolic\stats_probability.rs        43        8        23       12          0
~hes\numerical\finite_field.rs        43        8        23       12          0
~\symbolic\stats_regression.rs        43        8        23       12          0
~es\symbolic\thermodynamics.rs        43        8        23       12          0
~erical\functional_analysis.rs        43        8        23       12          0
~s\symbolic\vector_calculus.rs        43        8        23       12          0
benches\symbolic\tensor.rs            43        8        23       12          0
benches\symbolic\topology.rs          43        8        23       12          0
~enches\symbolic\transforms.rs        43        8        23       12          0
~cs_sim\geodesic_relativity.rs        43        8        23       12          0
~sics_sim\ising_statistical.rs        43        8        23       12          0
benches\symbolic\vector.rs            43        8        23       12          0
benches\numerical\testing.rs          43        8        23       12          0
benches\ffi_apis\mod.rs               43        8        23       12          0
benches\output\io.rs                  43        8        23       12          0
benches\physics\physics_sm.rs         43        8        23       12          0
benches\ffi_apis\ffi_api.rs           43        8        23       12          0
~enches\physics\physics_rkm.rs        43        8        23       12          0
benches\physics\physics_mm.rs         43        8        23       12          0
~enches\physics\physics_mtm.rs        43        8        23       12          0
benches\ffi_blindings\mod.rs          43        8        23       12          0
~nches\numerical\elementary.rs        43        8        23       12          0
benches\numerical\series.rs           43        8        23       12          0
~ches\numerical\coordinates.rs        43        8        23       12          0
~enches\physics\physics_fem.rs        43        8        23       12          0
~enches\physics\physics_fdm.rs        43        8        23       12          0
~ical\differential_geometry.rs        43        8        23       12          0
~enches\physics\physics_fvm.rs        43        8        23       12          0
~umerical\computer_graphics.rs        43        8        23       12          0
benches\physics\physics_em.rs         43        8        23       12          0
~enches\physics\physics_bem.rs        43        8        23       12          0
benches\numerical\vector.rs           43        8        23       12          0
~ches\numerical\convergence.rs        43        8        23       12          0
~numerical\complex_analysis.rs        43        8        23       12          0
~nches\numerical\transforms.rs        43        8        23       12          0
~es\numerical\combinatorics.rs        43        8        23       12          0
~enches\physics\physics_cnm.rs        43        8        23       12          0
benches\physics\mod.rs                43        8        23       12          0
~cal\calculus_of_variations.rs        43        8        23       12          0
benches\output\typst.rs               43        8        23       12          0
benches\numerical\calculus.rs         43        8        23       12          0
~enches\output\pretty_print.rs        43        8        23       12          0
benches\output\plotting.rs            43        8        23       12          0
benches\numerical\sparse.rs           43        8        23       12          0
benches\output\mod.rs                 43        8        23       12          0
~cal\calculus_of_variations.rs        40        0        23       17          0
tests\symbolic\core.rs                37        6        12       19          1
tests\regression_test.rs              36        8         2       26          3
~bolic\extra_simplify_tests.rs        35        3         0       32          0
~\symbolic\unit_unification.rs        29        4         1       24          5
src\numerical\special.rs              27        0         6       21          0
src\output\mod.rs                     21        0        16        5          0
tests\numerical\ode.rs                17        1         0       16          0
src\numerical\pde.rs                  16        0        10        6          0
~rc\physics\physics_sim\mod.rs         7        0         0        7          0
tests\symbolic\mod.rs                  6        0         2        4          0
src\ffi_apis\mod.rs                    5        0         4        1          0
src\ffi_blindings\mod.rs               2        0         2        0          0
benches\rssn_benches.rs                2        0         1        1          0
tests\numerical\mod.rs                 2        0         0        2          0
tests\plugins\mod.rs                   0        0         0        0          0
───────────────────────────────────────────────────────────────────────────────
Markdown                     4      1233      182         0     1051          0
───────────────────────────────────────────────────────────────────────────────
README.md                            936       96         0      840          0
CONTRIBUTING.md                      203       54         0      149          0
~\ISSUE_TEMPLATE\bug_report.md        49       15         0       34          0
~E_TEMPLATE\feature_request.md        45       17         0       28          0
───────────────────────────────────────────────────────────────────────────────
TOML                         3       128       14         1      113          2
───────────────────────────────────────────────────────────────────────────────
Cargo.toml                           116       12         1      103          2
~ins\example_plugin\Cargo.toml        10        2         0        8          0
.cargo\config.toml                     2        0         0        2          0
───────────────────────────────────────────────────────────────────────────────
YAML                         2       306       53        17      236          0
───────────────────────────────────────────────────────────────────────────────
.gitea\workflows\ci.yml              160       29         1      130          0
.github\workflows\ci.yml             146       24        16      106          0
───────────────────────────────────────────────────────────────────────────────
Batch                        1         3        0         0        3          0
───────────────────────────────────────────────────────────────────────────────
pre-commit.bat                         3        0         0        3          0
───────────────────────────────────────────────────────────────────────────────
Fortran Modern               1       167       34        33      100          6
───────────────────────────────────────────────────────────────────────────────
examples\fortran\test.f90            167       34        33      100          6
───────────────────────────────────────────────────────────────────────────────
License                      1        73       32         0       41          0
───────────────────────────────────────────────────────────────────────────────
LICENSE                               73       32         0       41          0
───────────────────────────────────────────────────────────────────────────────
Shell                        1         5        1         1        3          0
───────────────────────────────────────────────────────────────────────────────
pre-commit.sh                          5        1         1        3          0
───────────────────────────────────────────────────────────────────────────────
Total                      440     58294     2584     15876    39834       4804
───────────────────────────────────────────────────────────────────────────────
Estimated Cost to Develop (organic) $1,293,784
Estimated Schedule Effort (organic) 15.17 months
Estimated People Required (organic) 7.58
───────────────────────────────────────────────────────────────────────────────
Processed 2291136 bytes, 2.291 megabytes (SI)
───────────────────────────────────────────────────────────────────────────────
```

```
───────────────────────────────────────────────────────────────────────────────
Language                 Files     Lines   Blanks  Comments     Code Complexity
───────────────────────────────────────────────────────────────────────────────
Rust                       143     43564      302      9149    34113       4776
───────────────────────────────────────────────────────────────────────────────
Total                      143     43564      302      9149    34113       4776
───────────────────────────────────────────────────────────────────────────────
Estimated Cost to Develop (organic) $1,099,413
Estimated Schedule Effort (organic) 14.26 months
Estimated People Required (organic) 6.85
───────────────────────────────────────────────────────────────────────────────
Processed 1685045 bytes, 1.685 megabytes (SI)
───────────────────────────────────────────────────────────────────────────────

───────────────────────────────────────────────────────────────────────────────
Language                 Files     Lines   Blanks  Comments     Code Complexity
───────────────────────────────────────────────────────────────────────────────
Rust                       143     43564      302      9149    34113       4776
───────────────────────────────────────────────────────────────────────────────
symbolic\core.rs                    3490      147       267     3076        137
ffi_apis\ffi_api.rs                 2940        3       130     2807        375
symbolic\calculus.rs                1920       22       194     1704        365
symbolic\graph_algorithms.rs        1446        0       297     1149        289
symbolic\simplify.rs                1256       12        58     1186        176
symbolic\simplify_dag.rs            1011       58       155      798        135
prelude.rs                           928        0        16      912          0
symbolic\polynomial.rs               910        0       208      702        130
symbolic\matrix.rs                   888        0       213      675        187
symbolic\solve.rs                    794        0        89      705        135
symbolic\ode.rs                      784        0       143      641        112
numerical\testing.rs                 769        0       131      638        123
symbolic\pde.rs                      735        0       123      612        101
symbolic\coordinates.rs              658        0       168      490         64
numerical\matrix.rs                  643        0       139      504        110
symbolic\computer_graphics.rs        634        0       164      470          7
symbolic\number_theory.rs            541        0       148      393        108
symbolic\logic.rs                    535        0        74      461         88
symbolic\special_functions.rs        524        0       209      315         64
symbolic\transforms.rs               499        0       181      318         39
symbolic\tensor.rs                   495        0       173      322         63
symbolic\topology.rs                 482        0       184      298         55
~ic\error_correction_helper.rs       481        0       143      338         60
~ymbolic\poly_factorization.rs       476        0        94      382         58
symbolic\combinatorics.rs            457        0       166      291         47
numerical\optimize.rs                452        2        33      417         35
numerical\physics.rs                 444        0       113      331         64
symbolic\integration.rs              402        0        44      358         35
symbolic\finite_field.rs             401        0        75      326         51
symbolic\stats_probability.rs        387        0       163      224          0
physics\physics_fvm.rs               386        0        66      320         38
symbolic\cas_foundations.rs          386        4        40      342         57
~ymbolic\integral_equations.rs       357        0       120      237          7
symbolic\series.rs                   355        0       124      231         31
numerical\stats.rs                   346        0       114      232         12
output\pretty_print.rs               330       11         8      311         36
physics\physics_mtm.rs               301        0        23      278         49
output\plotting.rs                   283       10         8      265         15
physics\physics_fem.rs               280        0        15      265         54
physics\physics_rkm.rs               279        0        49      230         14
physics\physics_sm.rs                272        0        12      260         12
output\latex.rs                      268       16         8      244         14
symbolic\unit_unification.rs         255        0         4      251         26
symbolic\geometric_algebra.rs        253        0        85      168         38
symbolic\vector.rs                   252        0       122      130          3
symbolic\elementary.rs               247        0        52      195         11
~ic\lie_groups_and_algebras.rs       238        0        76      162          7
lib.rs                               238        5       148       85          1
symbolic\proof.rs                    236        0        81      155         33
symbolic\rewriting.rs                231        0        39      192         46
physics\physics_em.rs                230        0        85      145         10
symbolic\graph.rs                    230        0       104      126         15
numerical\interpolate.rs             227        0        56      171         41
symbolic\cad.rs                      224        0        10      214         35
numerical\finite_field.rs            215        0        71      144         27
~olic\differential_geometry.rs       213        0        67      146         20
symbolic\grobner.rs                  212        0        38      174         19
symbolic\error_correction.rs         208        0        59      149         27
physics\physics_bem.rs               205        0        34      171         17
numerical\error_correction.rs        204        0        36      168         27
numerical\sparse.rs                  202        0        69      133         17
physics\physics_mm.rs                200        0         8      192         21
~mbolic\functional_analysis.rs       198        0       106       92          2
numerical\polynomial.rs              197        0        44      153         30
~c\stats_information_theory.rs       197        0        90      107         12
symbolic\real_roots.rs               194        0        49      145         28
physics\physics_fdm.rs               188        0        36      152         26
symbolic\discrete_groups.rs          182        0        42      140         25
~umerical\computer_graphics.rs       178        0        72      106          7
numerical\coordinates.rs             173        0        63      110         16
plugins\manager.rs                   172        0        22      150         15
physics\physics_cnm.rs               164        0        17      147         17
symbolic\thermodynamics.rs           160        0       100       60          0
~h_isomorphism_and_coloring.rs       159        0        42      117         29
symbolic\vector_calculus.rs          157        0        62       95          0
symbolic\convergence.rs              157        0        47      110         34
~sics_sim\linear_elasticity.rs       156        0        22      134         20
symbolic\complex_analysis.rs         151        0        59       92          8
numerical\complex_analysis.rs        148        0        61       87          3
symbolic\group_theory.rs             146        0        61       85         14
output\typst.rs                      146       12         3      131         11
symbolic\optimize.rs                 145        0        45      100         18
numerical\integrate.rs               145        0        51       94         11
~mbolic\classical_mechanics.rs       145        0        93       52          0
symbolic\quantum_mechanics.rs        143        0        75       68          0
numerical\vector_calculus.rs         134        0        44       90          7
symbolic\cryptography.rs             133        0        56       77         13
~umerical\geometric_algebra.rs       133        0        22      111          4
numerical\solve.rs                   131        0        34       97         20
numerical\physics_cfd.rs             128        0        44       84         11
numerical\real_roots.rs              128        0        26      102         43
~cs_sim\schrodinger_quantum.rs       127        0        10      117         14
numerical\transforms.rs              125        0        41       84         19
symbolic\relativity.rs               125        0        67       58          0
symbolic\graph_operations.rs         124        0         5      119         41
numerical\tensor.rs                  122        0        37       85         15
symbolic\numeric.rs                  119        0        24       95         32
symbolic\stats.rs                    117        0        59       58          7
symbolic\stats_regression.rs         116        0        45       71          5
~cs_sim\geodesic_relativity.rs       116        0        26       90          3
numerical\vector.rs                  115        0        68       47         14
numerical\number_theory.rs           114        0        43       71         31
output\io.rs                         112        0        43       69          9
numerical\graph.rs                   111        0        44       67          7
numerical\convergence.rs             111        0        41       70          9
numerical\topology.rs                108        0        25       83         19
symbolic\multi_valued.rs             108        0        65       43          0
~cs_sim\navier_stokes_fluid.rs       106        0         4      102         17
symbolic\electromagnetism.rs         102        0        52       50          0
~bolic\quantum_field_theory.rs       100        0        60       40          2
numerical\physics_fea.rs              98        0        41       57         10
~lic\calculus_of_variations.rs        98        0        69       29          0
numerical\combinatorics.rs            97        0        40       57         15
numerical\elementary.rs               93        0        29       64          1
~erical\functional_analysis.rs        91        0        47       44          4
numerical\physics_md.rs               90        0        35       55          4
numerical\ode.rs                      89        0        20       69          4
symbolic\mod.rs                       88        0        20       68          0
~sics_sim\ising_statistical.rs        86        0        10       76         11
~sics_sim\gpe_superfluidity.rs        85        0         6       79          5
symbolic\special.rs                   82        0        61       21          0
numerical\signal.rs                   82        0        26       56          4
~s_sim\fdtd_electrodynamics.rs        81        0        10       71         16
symbolic\stats_inference.rs           81        0        19       62          0
~ical\differential_geometry.rs        79        0        20       59         12
~fractal_geometry_and_chaos.rs        75        0        31       44          5
~mbolic\solid_state_physics.rs        65        0        40       25          0
plugins\mod.rs                        64        0        27       37          2
symbolic\radicals.rs                  59        0        16       43         10
numerical\calculus.rs                 57        0        29       28          4
numerical\multi_valued.rs             53        0        19       34          5
~fractal_geometry_and_chaos.rs        51        0        30       21          1
numerical\series.rs                   51        0        15       36          2
physics\mod.rs                        48        0        37       11          0
symbolic\handles.rs                   47        0        14       33          0
numerical\mod.rs                      45        0         4       41          0
~cal\calculus_of_variations.rs        40        0        23       17          0
numerical\special.rs                  27        0         6       21          0
output\mod.rs                         21        0        16        5          0
numerical\pde.rs                      16        0        10        6          0
physics\physics_sim\mod.rs             7        0         0        7          0
ffi_apis\mod.rs                        5        0         4        1          0
ffi_blindings\mod.rs                   2        0         2        0          0
───────────────────────────────────────────────────────────────────────────────
Total                      143     43564      302      9149    34113       4776
───────────────────────────────────────────────────────────────────────────────
Estimated Cost to Develop (organic) $1,099,413
Estimated Schedule Effort (organic) 14.26 months
Estimated People Required (organic) 6.85
───────────────────────────────────────────────────────────────────────────────
Processed 1685045 bytes, 1.685 megabytes (SI)
───────────────────────────────────────────────────────────────────────────────
```