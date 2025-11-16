# Code Stastics (SCC CST 2025-11-16)

```
───────────────────────────────────────────────────────────────────────────────
Language            Files       Lines    Blanks  Comments       Code Complexity
───────────────────────────────────────────────────────────────────────────────
Rust                  439      60,921     2,610    16,104     42,207      5,187
Markdown               11       2,202       418         0      1,784          0
YAML                    5         466        87        50        329          0
TOML                    4         158        14         5        139          2
Batch                   3         169        24        12        133         28
Shell                   3         136        14        17        105         13
XML                     3          27         0         0         27          0
Fortran Modern          1         167        34        33        100          6
Handlebars              1          16         3         0         13          0
License                 1          73        32         0         41          0
───────────────────────────────────────────────────────────────────────────────
Total                 471      64,335     3,236    16,221     44,878      5,236
───────────────────────────────────────────────────────────────────────────────
Estimated Cost to Develop (semi-detached) $17,001,389
Estimated Schedule Effort (semi-detached) 16.31 months
Estimated People Required (semi-detached) 13.03
───────────────────────────────────────────────────────────────────────────────
Processed 2501974 bytes, 2.502 megabytes (SI)
───────────────────────────────────────────────────────────────────────────────

───────────────────────────────────────────────────────────────────────────────
Language            Files       Lines    Blanks  Comments       Code Complexity
───────────────────────────────────────────────────────────────────────────────
Rust                  439      60,921     2,610    16,104     42,207      5,187
───────────────────────────────────────────────────────────────────────────────
src\symbolic\core.rs            3,970       156       351      3,463        182
src\ffi_apis\ffi_api.rs         2,940         3       130      2,807        375
src\input\parser.rs             2,582       207        92      2,283         36
src\symbolic\calculus.rs        1,920        22       194      1,704        365
~c\symbolic\simplify_dag.rs     1,507        63       182      1,262        319
~mbolic\graph_algorithms.rs     1,446         0       297      1,149        289
src\symbolic\simplify.rs        1,256        12        58      1,186        179
src\prelude.rs                    946         0        16        930          0
src\symbolic\polynomial.rs        910         0       208        702        130
src\symbolic\matrix.rs            889         0       214        675        187
src\symbolic\solve.rs             794         0        89        705        135
src\symbolic\ode.rs               784         0       143        641        112
src\numerical\testing.rs          769         0       131        638        123
src\numerical\matrix.rs           758        25       162        571        153
src\symbolic\pde.rs               735         0       123        612        101
~rc\symbolic\coordinates.rs       658         0       168        490         64
~bolic\computer_graphics.rs       634         0       164        470          7
~\symbolic\number_theory.rs       541         0       148        393        108
src\symbolic\logic.rs             535         0        74        461         90
~bolic\special_functions.rs       524         0       209        315         69
src\symbolic\transforms.rs        499         0       181        318         39
src\symbolic\tensor.rs            495         0       173        322         63
src\symbolic\topology.rs          485         2       185        298         55
~error_correction_helper.rs       481         0       143        338         60
~olic\poly_factorization.rs       476         0        94        382         58
~\symbolic\combinatorics.rs       457         0       166        291         47
src\numerical\optimize.rs         452         2        33        417         35
src\numerical\physics.rs          444         0       113        331         64
~rc\symbolic\integration.rs       403         0        45        358         35
~c\symbolic\finite_field.rs       401         0        75        326         51
~bolic\stats_probability.rs       387         0       163        224          0
src\physics\physics_fvm.rs        386         0        66        320         38
~ymbolic\cas_foundations.rs       386         4        40        342         57
~olic\integral_equations.rs       357         0       120        237          7
src\symbolic\series.rs            355         0       124        231         31
src\numerical\stats.rs            346         0       114        232         12
src\output\pretty_print.rs        330        11         8        311         36
src\physics\physics_mtm.rs        301         0        23        278         49
src\output\plotting.rs            283        10         8        265         15
src\physics\physics_fem.rs        280         0        15        265         54
src\physics\physics_rkm.rs        279         0        49        230         14
src\physics\physics_sm.rs         272         0        12        260         12
src\plugins\manager.rs            271        17        22        232         21
src\output\latex.rs               268        16         8        244         14
~mbolic\unit_unification.rs       255         0         4        251         26
~bolic\geometric_algebra.rs       254         0        86        168         38
src\symbolic\vector.rs            252         0       122        130          3
src\lib.rs                        252         6       153         93          1
~rc\numerical\elementary.rs       250         2        34        214         46
src\symbolic\elementary.rs        247         0        52        195         11
~lie_groups_and_algebras.rs       238         0        76        162          7
src\symbolic\proof.rs             236         0        81        155         33
src\symbolic\rewriting.rs         231         0        39        192         46
src\physics\physics_em.rs         230         0        85        145         10
src\symbolic\graph.rs             230         0       104        126         15
~c\numerical\interpolate.rs       227         0        56        171         41
src\symbolic\cad.rs               224         0        10        214         35
~\numerical\finite_field.rs       215         0        71        144         27
~c\differential_geometry.rs       213         0        67        146         20
src\symbolic\grobner.rs           212         0        38        174         19
~mbolic\error_correction.rs       208         0        59        149         27
src\compute\engine.rs             206        18         2        186         11
src\physics\physics_bem.rs        205         0        34        171         17
~erical\error_correction.rs       204         0        36        168         27
src\numerical\sparse.rs           202         0        69        133         17
src\physics\physics_mm.rs         200         0         8        192         21
~lic\functional_analysis.rs       199         0       107         92          2
~tats_information_theory.rs       197         0        90        107         12
~rc\numerical\polynomial.rs       197         0        44        153         30
src\symbolic\real_roots.rs        194         0        49        145         28
src\physics\physics_fdm.rs        188         0        36        152         26
src\numerical\integrate.rs        184         2        59        123         16
~ymbolic\discrete_groups.rs       182         0        42        140         25
~rical\computer_graphics.rs       178         0        72        106          7
~c\numerical\coordinates.rs       173         0        63        110         16
~s_sim\linear_elasticity.rs       167         0        33        134         20
src\physics\physics_cnm.rs        164         0        17        147         17
~symbolic\thermodynamics.rs       160         0       100         60          0
~somorphism_and_coloring.rs       159         0        42        117         29
~ymbolic\vector_calculus.rs       157         0        62         95          0
~ts\numerical\polynomial.rs       157        27        34         96          4
~rc\symbolic\convergence.rs       157         0        47        110         34
~mbolic\complex_analysis.rs       151         0        59         92          8
~erical\complex_analysis.rs       148         0        61         87          3
src\output\typst.rs               146        12         3        131         11
~c\symbolic\group_theory.rs       146         0        61         85         14
~lic\classical_mechanics.rs       145         0        93         52          0
src\symbolic\optimize.rs          145         0        45        100         18
~bolic\quantum_mechanics.rs       143         0        75         68          0
~merical\vector_calculus.rs       134         0        44         90          7
~c\symbolic\cryptography.rs       133         0        56         77         13
~rical\geometric_algebra.rs       133         0        22        111          4
src\numerical\solve.rs            131         0        34         97         20
~c\numerical\physics_cfd.rs       128         0        44         84         11
~rc\numerical\real_roots.rs       128         0        26        102         43
~sim\schrodinger_quantum.rs       127         0        10        117         14
src\symbolic\relativity.rs        125         0        67         58          0
~rc\numerical\transforms.rs       125         0        41         84         19
~mbolic\graph_operations.rs       124         0         5        119         41
src\numerical\tensor.rs           122         0        37         85         15
src\symbolic\numeric.rs           119         0        24         95         32
src\symbolic\stats.rs             117         0        59         58          7
~sim\geodesic_relativity.rs       116         0        26         90          3
~mbolic\stats_regression.rs       116         0        45         71          5
src\numerical\vector.rs           115         0        68         47         14
~numerical\number_theory.rs       114         0        43         71         31
src\output\io.rs                  112         0        43         69          9
~c\numerical\convergence.rs       111         0        41         70          9
src\numerical\graph.rs            111         0        44         67          7
~c\symbolic\multi_valued.rs       108         0        65         43          0
src\numerical\topology.rs         108         0        25         83         19
~sim\navier_stokes_fluid.rs       106         0         4        102         17
~mbolic\electromagnetism.rs       102         0        52         50          0
~ic\quantum_field_theory.rs       100         0        60         40          2
tests\symbolic\simplify.rs         98         9        25         64          0
~c\numerical\physics_fea.rs        98         0        41         57         10
~\calculus_of_variations.rs        98         0        69         29          0
~numerical\combinatorics.rs        97         0        40         57         15
tests\symbolic\calculus.rs         96        10        25         61          0
~cal\functional_analysis.rs        91         0        47         44          4
~rc\numerical\physics_md.rs        90         0        35         55          4
src\numerical\ode.rs               89         0        20         69          4
src\symbolic\mod.rs                88         0        20         68          0
~s_sim\ising_statistical.rs        86         0        10         76         11
~s_sim\gpe_superfluidity.rs        85         0         6         79          5
src\symbolic\special.rs            82         0        61         21          0
src\numerical\signal.rs            82         0        26         56          4
~ymbolic\stats_inference.rs        81         0        19         62          0
~im\fdtd_electrodynamics.rs        81         0        10         71         16
~l\differential_geometry.rs        79         0        20         59         12
~ests\numerical\calculus.rs        76        14        11         51          3
~ctal_geometry_and_chaos.rs        75         0        31         44          5
~s\numerical\interpolate.rs        72        11        13         48          7
src\compute\cache.rs               66        10         0         56          2
src\plugins\plugin_c.rs            66         1        27         38          2
~lic\solid_state_physics.rs        65         0        40         25          0
src\symbolic\radicals.rs           59         0        16         43         10
src\numerical\calculus.rs          57         0        29         28          4
~\numerical\multi_valued.rs        53         0        19         34          5
src\numerical\series.rs            51         0        15         36          2
~ctal_geometry_and_chaos.rs        51         0        30         21          1
tests\mod.rs                       50         6        42          2          0
~s\numerical\coordinates.rs        48         6        25         17          0
~sts\physics\physics_mtm.rs        48         6        25         17          0
~ests\numerical\optimize.rs        48         6        25         17          0
~\numerical\multi_valued.rs        48         6        25         17          0
~numerical\number_theory.rs        48         6        25         17          0
tests\numerical\pde.rs             48         6        25         17          0
tests\symbolic\graph.rs            48         6        25         17          0
tests\numerical\physics.rs         48         6        25         17          0
tests\symbolic\special.rs          48         6        25         17          0
~s\numerical\physics_fea.rs        48         6        25         17          0
tests\numerical\graph.rs           48         6        25         17          0
tests\symbolic\ode.rs              48         6        25         17          0
~s\numerical\physics_cfd.rs        48         6        25         17          0
~ts\numerical\physics_md.rs        48         6        25         17          0
tests\numerical\matrix.rs          48         6        25         17          0
~ts\numerical\real_roots.rs        48         6        25         17          0
~sts\numerical\integrate.rs        48         6        25         17          0
tests\numerical\series.rs          48         6        25         17          0
tests\numerical\signal.rs          48         6        25         17          0
~somorphism_and_coloring.rs        48         6        25         17          0
tests\numerical\solve.rs           48         6        25         17          0
tests\numerical\sparse.rs          48         6        25         17          0
~mbolic\graph_operations.rs        48         6        25         17          0
~bolic\geometric_algebra.rs        48         6        25         17          0
tests\numerical\stats.rs           48         6        25         17          0
tests\numerical\special.rs         48         6        25         17          0
~lic\solid_state_physics.rs        48         6        25         17          0
tests\symbolic\grobner.rs          48         6        25         17          0
tests\numerical\tensor.rs          48         6        25         17          0
~rical\geometric_algebra.rs        48         6        25         17          0
tests\numerical\testing.rs         48         6        25         17          0
~ts\numerical\transforms.rs        48         6        25         17          0
tests\symbolic\optimize.rs         48         6        25         17          0
~ests\numerical\topology.rs        48         6        25         17          0
tests\symbolic\solve.rs            48         6        25         17          0
~bolic\special_functions.rs        48         6        25         17          0
tests\numerical\vector.rs          48         6        25         17          0
~error_correction_helper.rs        48         6        25         17          0
~ctal_geometry_and_chaos.rs        48         6        25         17          0
~lic\functional_analysis.rs        48         6        25         17          0
~merical\vector_calculus.rs        48         6        25         17          0
~ymbolic\stats_inference.rs        48         6        25         17          0
src\physics\mod.rs                 48         0        37         11          0
tests\symbolic\stats.rs            48         6        25         17          0
tests\symbolic\series.rs           48         6        25         17          0
~s\symbolic\finite_field.rs        48         6        25         17          0
~erical\complex_analysis.rs        48         6        25         17          0
tests\output\io.rs                 48         6        25         17          0
~tats_information_theory.rs        48         6        25         17          0
tests\output\latex.rs              48         6        25         17          0
~cal\functional_analysis.rs        48         6        25         17          0
tests\output\mod.rs                48         6        25         17          0
tests\output\plotting.rs           48         6        25         17          0
~ests\symbolic\rewriting.rs        48         6        25         17          0
~sts\output\pretty_print.rs        48         6        25         17          0
~bolic\stats_probability.rs        48         6        25         17          0
~ctal_geometry_and_chaos.rs        48         6        25         17          0
~s\symbolic\group_theory.rs        48         6        25         17          0
tests\output\typst.rs              48         6        25         17          0
tests\symbolic\tensor.rs           48         6        25         17          0
~mbolic\error_correction.rs        48         6        25         17          0
tests\symbolic\pde.rs              48         6        25         17          0
~erical\error_correction.rs        48         6        25         17          0
~sim\schrodinger_quantum.rs        48         6        25         17          0
~c\differential_geometry.rs        48         6        25         17          0
tests\symbolic\topology.rs         48         6        25         17          0
tests\symbolic\numeric.rs          48         6        25         17          0
~sts\symbolic\polynomial.rs        48         6        25         17          0
~mbolic\electromagnetism.rs        48         6        25         17          0
~sts\symbolic\transforms.rs        48         6        25         17          0
~sim\navier_stokes_fluid.rs        48         6        25         17          0
~mbolic\stats_regression.rs        48         6        25         17          0
~sts\symbolic\elementary.rs        48         6        25         17          0
~ymbolic\discrete_groups.rs        48         6        25         17          0
~s\symbolic\cryptography.rs        48         6        25         17          0
~\numerical\finite_field.rs        48         6        25         17          0
~\calculus_of_variations.rs        48         6        25         17          0
~ts\numerical\elementary.rs        48         6        25         17          0
~ymbolic\cas_foundations.rs        48         6        25         17          0
tests\physics\mod.rs               48         6        25         17          0
~sts\symbolic\relativity.rs        48         6        25         17          0
~physics\physics_sim\mod.rs        48         6        25         17          0
~sts\physics\physics_bem.rs        48         6        25         17          0
~l\differential_geometry.rs        48         6        25         17          0
~sts\physics\physics_cnm.rs        48         6        25         17          0
~ts\symbolic\coordinates.rs        48         6        25         17          0
~ymbolic\vector_calculus.rs        48         6        25         17          0
~rical\computer_graphics.rs        48         6        25         17          0
~s\numerical\convergence.rs        48         6        25         17          0
~ests\physics\physics_em.rs        48         6        25         17          0
~sts\physics\physics_fdm.rs        48         6        25         17          0
~ts\symbolic\convergence.rs        48         6        25         17          0
~sts\physics\physics_fvm.rs        48         6        25         17          0
~numerical\combinatorics.rs        48         6        25         17          0
~\symbolic\combinatorics.rs        48         6        25         17          0
~sts\physics\physics_fem.rs        48         6        25         17          0
tests\symbolic\proof.rs            48         6        25         17          0
~sts\symbolic\real_roots.rs        48         6        25         17          0
~olic\integral_equations.rs        48         6        25         17          0
~bolic\computer_graphics.rs        48         6        25         17          0
~mbolic\graph_algorithms.rs        48         6        25         17          0
~ts\symbolic\integration.rs        48         6        25         17          0
~sts\physics\physics_rkm.rs        48         6        25         17          0
~symbolic\thermodynamics.rs        48         6        25         17          0
~mbolic\complex_analysis.rs        48         6        25         17          0
~ests\physics\physics_mm.rs        48         6        25         17          0
~ests\physics\physics_sm.rs        48         6        25         17          0
~lie_groups_and_algebras.rs        48         6        25         17          0
tests\symbolic\matrix.rs           48         6        25         17          0
tests\symbolic\vector.rs           48         6        25         17          0
~\calculus_of_variations.rs        48         6        25         17          0
~lic\classical_mechanics.rs        48         6        25         17          0
~s_sim\linear_elasticity.rs        48         6        25         17          0
tests\symbolic\logic.rs            48         6        25         17          0
tests\symbolic\radicals.rs         48         6        25         17          0
~bolic\quantum_mechanics.rs        48         6        25         17          0
~ic\quantum_field_theory.rs        48         6        25         17          0
~sim\geodesic_relativity.rs        48         6        25         17          0
~im\fdtd_electrodynamics.rs        48         6        25         17          0
tests\symbolic\cad.rs              48         6        25         17          0
~s\symbolic\multi_valued.rs        48         6        25         17          0
~\symbolic\number_theory.rs        48         6        25         17          0
~s_sim\ising_statistical.rs        48         6        25         17          0
~s_sim\gpe_superfluidity.rs        48         6        25         17          0
~olic\poly_factorization.rs        48         6        25         17          0
tests\prelude.rs                   47         5        42          0          0
tests\lib.rs                       47         5        42          0          0
tests\ffi_blindings\mod.rs         47         5        42          0          0
tests\ffi_apis\mod.rs              47         5        42          0          0
tests\ffi_apis\ffi_api.rs          47         5        42          0          0
src\symbolic\handles.rs            47         0        14         33          0
~\example_plugin\src\lib.rs        47         7         9         31          1
src\numerical\mod.rs               45         0         4         41          0
~ymbolic\discrete_groups.rs        43         8        23         12          0
~hes\physics\physics_rkm.rs        43         8        23         12          0
~hes\symbolic\polynomial.rs        43         8        23         12          0
benches\symbolic\proof.rs          43         8        23         12          0
~ic\quantum_field_theory.rs        43         8        23         12          0
benches\symbolic\pde.rs            43         8        23         12          0
~bolic\quantum_mechanics.rs        43         8        23         12          0
~nches\symbolic\radicals.rs        43         8        23         12          0
~hes\symbolic\real_roots.rs        43         8        23         12          0
~nches\symbolic\optimize.rs        43         8        23         12          0
~ches\symbolic\rewriting.rs        43         8        23         12          0
~hes\symbolic\relativity.rs        43         8        23         12          0
benches\symbolic\ode.rs            43         8        23         12          0
~nches\symbolic\simplify.rs        43         8        23         12          0
benches\symbolic\series.rs         43         8        23         12          0
~enches\symbolic\numeric.rs        43         8        23         12          0
~\symbolic\number_theory.rs        43         8        23         12          0
~lic\solid_state_physics.rs        43         8        23         12          0
benches\symbolic\solve.rs          43         8        23         12          0
~enches\symbolic\special.rs        43         8        23         12          0
~s\symbolic\multi_valued.rs        43         8        23         12          0
~bolic\special_functions.rs        43         8        23         12          0
benches\symbolic\stats.rs          43         8        23         12          0
benches\symbolic\mod.rs            43         8        23         12          0
benches\symbolic\matrix.rs         43         8        23         12          0
benches\symbolic\logic.rs          43         8        23         12          0
~lie_groups_and_algebras.rs        43         8        23         12          0
~tats_information_theory.rs        43         8        23         12          0
~olic\integral_equations.rs        43         8        23         12          0
~ymbolic\stats_inference.rs        43         8        23         12          0
~es\symbolic\integration.rs        43         8        23         12          0
~s\symbolic\group_theory.rs        43         8        23         12          0
~mbolic\stats_regression.rs        43         8        23         12          0
~enches\symbolic\grobner.rs        43         8        23         12          0
~bolic\stats_probability.rs        43         8        23         12          0
~mbolic\graph_operations.rs        43         8        23         12          0
benches\symbolic\tensor.rs         43         8        23         12          0
~somorphism_and_coloring.rs        43         8        23         12          0
~symbolic\thermodynamics.rs        43         8        23         12          0
~mbolic\graph_algorithms.rs        43         8        23         12          0
benches\symbolic\graph.rs          43         8        23         12          0
~bolic\geometric_algebra.rs        43         8        23         12          0
benches\lib.rs                     43         8        23         12          0
~ctal_geometry_and_chaos.rs        43         8        23         12          0
~lic\functional_analysis.rs        43         8        23         12          0
~s\symbolic\finite_field.rs        43         8        23         12          0
~error_correction_helper.rs        43         8        23         12          0
~mbolic\error_correction.rs        43         8        23         12          0
~nches\symbolic\topology.rs        43         8        23         12          0
~mbolic\electromagnetism.rs        43         8        23         12          0
benches\symbolic\vector.rs         43         8        23         12          0
~ymbolic\vector_calculus.rs        43         8        23         12          0
~hes\symbolic\elementary.rs        43         8        23         12          0
benches\prelude.rs                 43         8        23         12          0
~c\differential_geometry.rs        43         8        23         12          0
~enches\ffi_apis\ffi_api.rs        43         8        23         12          0
~s\symbolic\cryptography.rs        43         8        23         12          0
~es\symbolic\convergence.rs        43         8        23         12          0
~es\symbolic\coordinates.rs        43         8        23         12          0
benches\symbolic\core.rs           43         8        23         12          0
~\calculus_of_variations.rs        43         8        23         12          0
~hes\symbolic\transforms.rs        43         8        23         12          0
~bolic\computer_graphics.rs        43         8        23         12          0
~mbolic\complex_analysis.rs        43         8        23         12          0
~\symbolic\combinatorics.rs        43         8        23         12          0
~lic\classical_mechanics.rs        43         8        23         12          0
~ymbolic\cas_foundations.rs        43         8        23         12          0
~nches\symbolic\calculus.rs        43         8        23         12          0
benches\symbolic\cad.rs            43         8        23         12          0
benches\ffi_apis\mod.rs            43         8        23         12          0
~sim\schrodinger_quantum.rs        43         8        23         12          0
~sim\navier_stokes_fluid.rs        43         8        23         12          0
benches\plugins\mod.rs             43         8        23         12          0
~s_sim\ising_statistical.rs        43         8        23         12          0
~s_sim\linear_elasticity.rs        43         8        23         12          0
~physics\physics_sim\mod.rs        43         8        23         12          0
~s_sim\gpe_superfluidity.rs        43         8        23         12          0
~sim\geodesic_relativity.rs        43         8        23         12          0
~im\fdtd_electrodynamics.rs        43         8        23         12          0
~nches\ffi_blindings\mod.rs        43         8        23         12          0
~\calculus_of_variations.rs        43         8        23         12          0
~olic\poly_factorization.rs        43         8        23         12          0
~ches\physics\physics_sm.rs        43         8        23         12          0
~hes\physics\physics_mtm.rs        43         8        23         12          0
~numerical\combinatorics.rs        43         8        23         12          0
~ches\physics\physics_mm.rs        43         8        23         12          0
~ches\numerical\calculus.rs        43         8        23         12          0
~rical\computer_graphics.rs        43         8        23         12          0
~s\numerical\convergence.rs        43         8        23         12          0
~hes\physics\physics_fvm.rs        43         8        23         12          0
~hes\physics\physics_fdm.rs        43         8        23         12          0
benches\physics\mod.rs             43         8        23         12          0
~hes\physics\physics_fem.rs        43         8        23         12          0
~hes\physics\physics_cnm.rs        43         8        23         12          0
~ches\physics\physics_em.rs        43         8        23         12          0
~hes\physics\physics_bem.rs        43         8        23         12          0
~numerical\number_theory.rs        43         8        23         12          0
~es\numerical\elementary.rs        43         8        23         12          0
benches\output\typst.rs            43         8        23         12          0
~hes\output\pretty_print.rs        43         8        23         12          0
benches\output\latex.rs            43         8        23         12          0
benches\output\plotting.rs         43         8        23         12          0
benches\output\mod.rs              43         8        23         12          0
benches\output\io.rs               43         8        23         12          0
~\numerical\multi_valued.rs        43         8        23         12          0
~s\numerical\coordinates.rs        43         8        23         12          0
~merical\vector_calculus.rs        43         8        23         12          0
~l\differential_geometry.rs        43         8        23         12          0
~enches\numerical\vector.rs        43         8        23         12          0
~erical\error_correction.rs        43         8        23         12          0
~ches\numerical\topology.rs        43         8        23         12          0
~es\numerical\transforms.rs        43         8        23         12          0
benches\numerical\mod.rs           43         8        23         12          0
~erical\complex_analysis.rs        43         8        23         12          0
~nches\numerical\testing.rs        43         8        23         12          0
benches\numerical\stats.rs         43         8        23         12          0
~enches\numerical\sparse.rs        43         8        23         12          0
~enches\numerical\tensor.rs        43         8        23         12          0
~s\numerical\interpolate.rs        43         8        23         12          0
benches\numerical\solve.rs         43         8        23         12          0
~nches\numerical\special.rs        43         8        23         12          0
~enches\numerical\signal.rs        43         8        23         12          0
~enches\numerical\series.rs        43         8        23         12          0
~\numerical\finite_field.rs        43         8        23         12          0
~es\numerical\physics_md.rs        43         8        23         12          0
~es\numerical\real_roots.rs        43         8        23         12          0
~es\numerical\polynomial.rs        43         8        23         12          0
~ctal_geometry_and_chaos.rs        43         8        23         12          0
~s\numerical\physics_fea.rs        43         8        23         12          0
~s\numerical\physics_cfd.rs        43         8        23         12          0
benches\numerical\pde.rs           43         8        23         12          0
benches\numerical\ode.rs           43         8        23         12          0
~nches\numerical\physics.rs        43         8        23         12          0
~ches\numerical\optimize.rs        43         8        23         12          0
~enches\numerical\matrix.rs        43         8        23         12          0
~cal\functional_analysis.rs        43         8        23         12          0
~hes\numerical\integrate.rs        43         8        23         12          0
~rical\geometric_algebra.rs        43         8        23         12          0
benches\numerical\graph.rs         43         8        23         12          0
~\calculus_of_variations.rs        40         0        23         17          0
tests\symbolic\core.rs             37         6        12         19          1
src\compute\computation.rs         37         4         1         32          0
tests\regression_test.rs           36         8         2         26          3
~ic\extra_simplify_tests.rs        35         3         0         32          0
src\constant.rs                    32         6         6         20          0
src\plugins\stable_abi.rs          32         6         0         26          0
~mbolic\unit_unification.rs        29         4         1         24          5
src\numerical\special.rs           27         0         6         21          0
src\output\mod.rs                  21         0        16          5          0
build.rs                           17         4         0         13          0
tests\numerical\ode.rs             17         1         0         16          0
src\numerical\pde.rs               16         0        10          6          0
src\compute\state.rs                7         1         1          5          0
~physics\physics_sim\mod.rs         7         0         0          7          0
src\compute\computable.rs           6         1         0          5          0
tests\symbolic\mod.rs               6         0         2          4          0
src\compute\mod.rs                  5         0         0          5          0
src\ffi_apis\mod.rs                 5         0         4          1          0
src\plugins\mod.rs                  3         0         0          3          0
tests\numerical\mod.rs              2         0         0          2          0
benches\rssn_benches.rs             2         0         1          1          0
src\ffi_blindings\mod.rs            2         0         2          0          0
src\input\mod.rs                    2         0         1          1          0
tests\plugins\mod.rs                0         0         0          0          0
───────────────────────────────────────────────────────────────────────────────
Markdown               11       2,202       418         0      1,784          0
───────────────────────────────────────────────────────────────────────────────
ATTRIBUTIONS.md                   874       220         0        654          0
CODESTASTICS.md                   674         4         0        670          0
CONTRIBUTING.md                   190        53         0        137          0
README.md                         143        52         0         91          0
TECHNICAL_DEBT.md                  80        15         0         65          0
ARCHITECTURE.md                    66        22         0         44          0
~SUE_TEMPLATE\bug_report.md        49        15         0         34          0
~EMPLATE\feature_request.md        45        17         0         28          0
~HNICAL_DEBT_FIX_SUMMARY.md        38         7         0         31          0
CODE_OF_CONDUCT.md                 32         9         0         23          0
SECURITY.md                        11         4         0          7          0
───────────────────────────────────────────────────────────────────────────────
YAML                    5         466        87        50        329          0
───────────────────────────────────────────────────────────────────────────────
.github\workflows\ci.yml          165        41         6        118          0
.gitea\workflows\ci.yml           164        30         1        133          0
~ithub\workflows\codeql.yml       101        10        43         48          0
~lows\dependency-review.yml        23         5         0         18          0
.github\dependabot.yml             13         1         0         12          0
───────────────────────────────────────────────────────────────────────────────
TOML                    4         158        14         5        139          2
───────────────────────────────────────────────────────────────────────────────
Cargo.toml                        134        11         5        118          2
about.toml                         12         1         0         11          0
~\example_plugin\Cargo.toml        10         2         0          8          0
.cargo\config.toml                  2         0         0          2          0
───────────────────────────────────────────────────────────────────────────────
Batch                   3         169        24        12        133         28
───────────────────────────────────────────────────────────────────────────────
build_all.bat                     101        13         4         84         21
collect_binaries.bat               65        11         8         46          7
pre-commit.bat                      3         0         0          3          0
───────────────────────────────────────────────────────────────────────────────
Shell                   3         136        14        17        105         13
───────────────────────────────────────────────────────────────────────────────
build_all.sh                       72        10        10         52         12
collect_binaries.sh                59         3         6         50          1
pre-commit.sh                       5         1         1          3          0
───────────────────────────────────────────────────────────────────────────────
XML                     3          27         0         0         27          0
───────────────────────────────────────────────────────────────────────────────
.idea\vcs.xml                      13         0         0         13          0
.idea\modules.xml                   8         0         0          8          0
.idea\rust.xml                      6         0         0          6          0
───────────────────────────────────────────────────────────────────────────────
Fortran Modern          1         167        34        33        100          6
───────────────────────────────────────────────────────────────────────────────
examples\fortran\test.f90         167        34        33        100          6
───────────────────────────────────────────────────────────────────────────────
Handlebars              1          16         3         0         13          0
───────────────────────────────────────────────────────────────────────────────
about.hbs                          16         3         0         13          0
───────────────────────────────────────────────────────────────────────────────
License                 1          73        32         0         41          0
───────────────────────────────────────────────────────────────────────────────
LICENSE                            73        32         0         41          0
───────────────────────────────────────────────────────────────────────────────
Total                 471      64,335     3,236    16,221     44,878      5,236
───────────────────────────────────────────────────────────────────────────────
Estimated Cost to Develop (semi-detached) $17,001,389
Estimated Schedule Effort (semi-detached) 16.31 months
Estimated People Required (semi-detached) 13.03
───────────────────────────────────────────────────────────────────────────────
Processed 2501974 bytes, 2.502 megabytes (SI)
───────────────────────────────────────────────────────────────────────────────
```

```
───────────────────────────────────────────────────────────────────────────────
Language            Files       Lines    Blanks  Comments       Code Complexity
───────────────────────────────────────────────────────────────────────────────
Rust                  154      47,980       619     9,420     37,941      5,163
───────────────────────────────────────────────────────────────────────────────
Total                 154      47,980       619     9,420     37,941      5,163
───────────────────────────────────────────────────────────────────────────────
Estimated Cost to Develop (semi-detached) $14,086,682
Estimated Schedule Effort (semi-detached) 15.27 months
Estimated People Required (semi-detached) 11.53
───────────────────────────────────────────────────────────────────────────────
Processed 1835568 bytes, 1.836 megabytes (SI)
───────────────────────────────────────────────────────────────────────────────

───────────────────────────────────────────────────────────────────────────────
Language            Files       Lines    Blanks  Comments       Code Complexity
───────────────────────────────────────────────────────────────────────────────
Rust                  154      47,980       619     9,420     37,941      5,163
───────────────────────────────────────────────────────────────────────────────
symbolic\core.rs                3,970       156       351      3,463        182
ffi_apis\ffi_api.rs             2,940         3       130      2,807        375
input\parser.rs                 2,582       207        92      2,283         36
symbolic\calculus.rs            1,920        22       194      1,704        365
symbolic\simplify_dag.rs        1,507        63       182      1,262        319
~mbolic\graph_algorithms.rs     1,446         0       297      1,149        289
symbolic\simplify.rs            1,256        12        58      1,186        179
prelude.rs                        946         0        16        930          0
symbolic\polynomial.rs            910         0       208        702        130
symbolic\matrix.rs                889         0       214        675        187
symbolic\solve.rs                 794         0        89        705        135
symbolic\ode.rs                   784         0       143        641        112
numerical\testing.rs              769         0       131        638        123
numerical\matrix.rs               758        25       162        571        153
symbolic\pde.rs                   735         0       123        612        101
symbolic\coordinates.rs           658         0       168        490         64
~bolic\computer_graphics.rs       634         0       164        470          7
symbolic\number_theory.rs         541         0       148        393        108
symbolic\logic.rs                 535         0        74        461         90
~bolic\special_functions.rs       524         0       209        315         69
symbolic\transforms.rs            499         0       181        318         39
symbolic\tensor.rs                495         0       173        322         63
symbolic\topology.rs              485         2       185        298         55
~error_correction_helper.rs       481         0       143        338         60
~olic\poly_factorization.rs       476         0        94        382         58
symbolic\combinatorics.rs         457         0       166        291         47
numerical\optimize.rs             452         2        33        417         35
numerical\physics.rs              444         0       113        331         64
symbolic\integration.rs           403         0        45        358         35
symbolic\finite_field.rs          401         0        75        326         51
~bolic\stats_probability.rs       387         0       163        224          0
~ymbolic\cas_foundations.rs       386         4        40        342         57
physics\physics_fvm.rs            386         0        66        320         38
~olic\integral_equations.rs       357         0       120        237          7
symbolic\series.rs                355         0       124        231         31
numerical\stats.rs                346         0       114        232         12
output\pretty_print.rs            330        11         8        311         36
physics\physics_mtm.rs            301         0        23        278         49
output\plotting.rs                283        10         8        265         15
physics\physics_fem.rs            280         0        15        265         54
physics\physics_rkm.rs            279         0        49        230         14
physics\physics_sm.rs             272         0        12        260         12
plugins\manager.rs                271        17        22        232         21
output\latex.rs                   268        16         8        244         14
~mbolic\unit_unification.rs       255         0         4        251         26
~bolic\geometric_algebra.rs       254         0        86        168         38
lib.rs                            252         6       153         93          1
symbolic\vector.rs                252         0       122        130          3
numerical\elementary.rs           250         2        34        214         46
symbolic\elementary.rs            247         0        52        195         11
~lie_groups_and_algebras.rs       238         0        76        162          7
symbolic\proof.rs                 236         0        81        155         33
symbolic\rewriting.rs             231         0        39        192         46
symbolic\graph.rs                 230         0       104        126         15
physics\physics_em.rs             230         0        85        145         10
numerical\interpolate.rs          227         0        56        171         41
symbolic\cad.rs                   224         0        10        214         35
numerical\finite_field.rs         215         0        71        144         27
~c\differential_geometry.rs       213         0        67        146         20
symbolic\grobner.rs               212         0        38        174         19
~mbolic\error_correction.rs       208         0        59        149         27
compute\engine.rs                 206        18         2        186         11
physics\physics_bem.rs            205         0        34        171         17
~erical\error_correction.rs       204         0        36        168         27
numerical\sparse.rs               202         0        69        133         17
physics\physics_mm.rs             200         0         8        192         21
~lic\functional_analysis.rs       199         0       107         92          2
numerical\polynomial.rs           197         0        44        153         30
~tats_information_theory.rs       197         0        90        107         12
symbolic\real_roots.rs            194         0        49        145         28
physics\physics_fdm.rs            188         0        36        152         26
numerical\integrate.rs            184         2        59        123         16
~ymbolic\discrete_groups.rs       182         0        42        140         25
~rical\computer_graphics.rs       178         0        72        106          7
numerical\coordinates.rs          173         0        63        110         16
~s_sim\linear_elasticity.rs       167         0        33        134         20
physics\physics_cnm.rs            164         0        17        147         17
symbolic\thermodynamics.rs        160         0       100         60          0
~somorphism_and_coloring.rs       159         0        42        117         29
~ymbolic\vector_calculus.rs       157         0        62         95          0
symbolic\convergence.rs           157         0        47        110         34
~mbolic\complex_analysis.rs       151         0        59         92          8
~erical\complex_analysis.rs       148         0        61         87          3
symbolic\group_theory.rs          146         0        61         85         14
output\typst.rs                   146        12         3        131         11
~lic\classical_mechanics.rs       145         0        93         52          0
symbolic\optimize.rs              145         0        45        100         18
~bolic\quantum_mechanics.rs       143         0        75         68          0
~merical\vector_calculus.rs       134         0        44         90          7
symbolic\cryptography.rs          133         0        56         77         13
~rical\geometric_algebra.rs       133         0        22        111          4
numerical\solve.rs                131         0        34         97         20
numerical\physics_cfd.rs          128         0        44         84         11
numerical\real_roots.rs           128         0        26        102         43
~sim\schrodinger_quantum.rs       127         0        10        117         14
numerical\transforms.rs           125         0        41         84         19
symbolic\relativity.rs            125         0        67         58          0
~mbolic\graph_operations.rs       124         0         5        119         41
numerical\tensor.rs               122         0        37         85         15
symbolic\numeric.rs               119         0        24         95         32
symbolic\stats.rs                 117         0        59         58          7
~sim\geodesic_relativity.rs       116         0        26         90          3
~mbolic\stats_regression.rs       116         0        45         71          5
numerical\vector.rs               115         0        68         47         14
numerical\number_theory.rs        114         0        43         71         31
output\io.rs                      112         0        43         69          9
numerical\graph.rs                111         0        44         67          7
numerical\convergence.rs          111         0        41         70          9
numerical\topology.rs             108         0        25         83         19
symbolic\multi_valued.rs          108         0        65         43          0
~sim\navier_stokes_fluid.rs       106         0         4        102         17
~mbolic\electromagnetism.rs       102         0        52         50          0
~ic\quantum_field_theory.rs       100         0        60         40          2
numerical\physics_fea.rs           98         0        41         57         10
~\calculus_of_variations.rs        98         0        69         29          0
numerical\combinatorics.rs         97         0        40         57         15
~cal\functional_analysis.rs        91         0        47         44          4
numerical\physics_md.rs            90         0        35         55          4
numerical\ode.rs                   89         0        20         69          4
symbolic\mod.rs                    88         0        20         68          0
~s_sim\ising_statistical.rs        86         0        10         76         11
~s_sim\gpe_superfluidity.rs        85         0         6         79          5
symbolic\special.rs                82         0        61         21          0
numerical\signal.rs                82         0        26         56          4
~im\fdtd_electrodynamics.rs        81         0        10         71         16
~ymbolic\stats_inference.rs        81         0        19         62          0
~l\differential_geometry.rs        79         0        20         59         12
~ctal_geometry_and_chaos.rs        75         0        31         44          5
compute\cache.rs                   66        10         0         56          2
plugins\plugin_c.rs                66         1        27         38          2
~lic\solid_state_physics.rs        65         0        40         25          0
symbolic\radicals.rs               59         0        16         43         10
numerical\calculus.rs              57         0        29         28          4
numerical\multi_valued.rs          53         0        19         34          5
~ctal_geometry_and_chaos.rs        51         0        30         21          1
numerical\series.rs                51         0        15         36          2
physics\mod.rs                     48         0        37         11          0
symbolic\handles.rs                47         0        14         33          0
numerical\mod.rs                   45         0         4         41          0
~\calculus_of_variations.rs        40         0        23         17          0
compute\computation.rs             37         4         1         32          0
constant.rs                        32         6         6         20          0
plugins\stable_abi.rs              32         6         0         26          0
numerical\special.rs               27         0         6         21          0
output\mod.rs                      21         0        16          5          0
numerical\pde.rs                   16         0        10          6          0
compute\state.rs                    7         1         1          5          0
physics\physics_sim\mod.rs          7         0         0          7          0
compute\computable.rs               6         1         0          5          0
ffi_apis\mod.rs                     5         0         4          1          0
compute\mod.rs                      5         0         0          5          0
plugins\mod.rs                      3         0         0          3          0
ffi_blindings\mod.rs                2         0         2          0          0
input\mod.rs                        2         0         1          1          0
───────────────────────────────────────────────────────────────────────────────
Total                 154      47,980       619     9,420     37,941      5,163
───────────────────────────────────────────────────────────────────────────────
Estimated Cost to Develop (semi-detached) $14,086,682
Estimated Schedule Effort (semi-detached) 15.27 months
Estimated People Required (semi-detached) 11.53
───────────────────────────────────────────────────────────────────────────────
Processed 1835568 bytes, 1.836 megabytes (SI)
───────────────────────────────────────────────────────────────────────────────

```
