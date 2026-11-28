#include "rssn.h"

#include <stdio.h>

int main(void) {
    struct RssnSession *session = rssn_session_new();
    if (session == NULL) {
        fprintf(stderr, "%s\n", rssn_last_error());
        return 1;
    }

    const char *operator_name = "num_number_theory_gcd_u128";
    char *rule_set = NULL;
    char *name = NULL;
    uint8_t hash[16];
    uint32_t inputs[2];
    uint32_t result = 0;
    char *result_text = NULL;

    enum RssnStatus status = rssn_op_key_at(operator_name, 0, &rule_set, &name, hash);
    if (status != RSSN_STATUS_OK) goto done;

    status = rssn_value_new_int_decimal(session, "79228162514264337593543950336", &inputs[0]);
    if (status != RSSN_STATUS_OK) goto done;
    status = rssn_value_new_int_decimal(session, "237684487542793012780631851008", &inputs[1]);
    if (status != RSSN_STATUS_OK) goto done;

    status = rssn_value_call_key(session, rule_set, name, hash, inputs, 2, &result);
    if (status != RSSN_STATUS_OK) goto done;
    status = rssn_value_to_text(session, result, &result_text);

done:
    if (status == RSSN_STATUS_OK) {
        puts(result_text);
    } else {
        fprintf(stderr, "%s\n", rssn_last_error());
    }
    rssn_string_free(result_text);
    rssn_string_free(name);
    rssn_string_free(rule_set);
    rssn_session_free(session);
    return status == RSSN_STATUS_OK ? 0 : 1;
}
