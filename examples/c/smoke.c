/* Minimal C client of rssn: differentiate, evaluate, extract a matrix and
 * run a simulation. Build (from the repository root):
 *
 *   DEV=1 cargo build --release
 *   cc -I. examples/c/smoke.c target/release/librssn.a -lm -lpthread -ldl -o smoke
 */
#include <stdio.h>
#include "rssn.h"

#define CHECK(call)                                                        \
    do {                                                                   \
        RssnStatus status_ = (call);                                       \
        if (status_ != RSSN_STATUS_OK) {                                   \
            fprintf(stderr, "%s failed (%d): %s\n", #call, (int)status_,   \
                    rssn_last_error());                                    \
            return 1;                                                      \
        }                                                                  \
    } while (0)

int main(void) {
    RssnSession *s = rssn_session_new();
    printf("rssn %s\n", rssn_version());

    /* d/dx sin(x)*cos(x), built from handles. */
    uint32_t x = rssn_term_sym(s, "x");
    uint32_t sx = rssn_term_apply(s, "sin", &x, 1);
    uint32_t cx = rssn_term_apply(s, "cos", &x, 1);
    uint32_t prod[2] = {sx, cx};
    uint32_t product = rssn_term_apply(s, "mul", prod, 2);
    uint32_t args[2] = {product, x};
    uint32_t request = rssn_term_apply(s, "diff", args, 2);

    RssnAnswer answer;
    CHECK(rssn_compute(s, request, NULL, &answer));
    char *text = rssn_term_to_string(s, answer.term);
    char *latex = rssn_term_to_latex(s, answer.term);
    printf("d/dx sin(x)cos(x) = %s   [%s]\n", text, latex);
    rssn_string_free(text);
    rssn_string_free(latex);

    /* The same request, numerically at x = 0.3. */
    RssnConfig *cfg = rssn_config_new();
    CHECK(rssn_config_numeric(cfg, 1e-12));
    CHECK(rssn_config_bind(cfg, "x", 0.3));
    CHECK(rssn_compute(s, request, cfg, &answer));
    printf("at x = 0.3: %.15f (+- %g)\n", answer.value, answer.error);
    rssn_config_free(cfg);

    /* Parse, compute, extract a dense matrix. */
    uint32_t m = rssn_simplify(s, "matmul(list(list(1, 2), list(3, 4)), list(list(0, 1), list(1, 0)))");
    RssnTensor t;
    CHECK(rssn_term_tensor_data(s, m, &t));
    printf("matrix %zux%zu:", t.shape[0], t.shape[1]);
    for (size_t i = 0; i < t.len; i++) printf(" %g", t.data[i]);
    printf("\n");
    rssn_tensor_free(&t);

    /* Errors carry a message. */
    if (rssn_term_parse(s, "sin(") == RSSN_NO_TERM)
        printf("parse error reported: %s\n", rssn_last_error());

    char *json = NULL;
    CHECK(rssn_sim_run("ising",
                       "{\"width\":8,\"height\":8,\"temperature\":2.0,\"mc_steps\":10,\"seed\":3}",
                       &json));
    printf("ising: %.60s...\n", json);
    rssn_string_free(json);

    rssn_session_free(s);
    return 0;
}
