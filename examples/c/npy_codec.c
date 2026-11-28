#include "rssn.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

int main(void) {
    struct RssnSession *session = rssn_session_new();
    if (session == NULL) {
        fprintf(stderr, "%s\n", rssn_last_error());
        return 1;
    }

    const size_t extent = 2;
    const double samples[2] = {-0.0, 1.5};
    uint32_t axis = 0, input[2] = {0, 0}, wire = 0, pair = 0, output = 0;
    char *hex = NULL;
    RssnTensor array = {0};
    int result = 1;

    RssnStatus status = rssn_value_new_int(session, 2, &axis);
    if (status != RSSN_STATUS_OK) goto done;
    status = rssn_value_new_list(session, &axis, 1, &input[0]);
    if (status != RSSN_STATUS_OK) goto done;
    status = rssn_value_new_tensor(session, samples, 2, &extent, 1, &input[1]);
    if (status != RSSN_STATUS_OK) goto done;

    status = rssn_value_call(session, "num_npy_encode", input, 2, &wire);
    if (status != RSSN_STATUS_OK) goto done;
    status = rssn_value_to_text(session, wire, &hex);
    if (status != RSSN_STATUS_OK) goto done;
    if (strncmp(hex, "934e554d5059", 12) != 0) {
        fprintf(stderr, "invalid NumPy magic\n");
        goto done;
    }

    // The same string can be supplied to decode without touching the filesystem.
    status = rssn_value_call(session, "num_npy_decode", &wire, 1, &pair);
    if (status != RSSN_STATUS_OK) goto done;
    status = rssn_value_at(session, pair, 1, &output);
    if (status != RSSN_STATUS_OK) goto done;
    status = rssn_value_to_tensor(session, output, &array);
    if (status != RSSN_STATUS_OK) goto done;
    if (array.rank != 1 || array.len != 2 || array.shape[0] != 2 ||
        !signbit(array.data[0]) || array.data[1] != samples[1]) {
        fprintf(stderr, "NumPy samples did not round trip\n");
        goto done;
    }

    puts("NumPy payload preserved rank, signed zero and sample order");
    result = 0;

done:
    if (status != RSSN_STATUS_OK) fprintf(stderr, "%s\n", rssn_last_error());
    rssn_tensor_free(&array);
    rssn_string_free(hex);
    rssn_session_free(session);
    return result;
}
