#include "rssn.h"

#include <stdio.h>

static RssnStatus print_text(const struct RssnSession *session, uint32_t cell) {
    char *text = NULL;
    RssnStatus status = rssn_value_to_text(session, cell, &text);
    if (status == RSSN_STATUS_OK) {
        puts(text);
    }
    rssn_string_free(text);
    return status;
}

static RssnStatus list_backends(const struct RssnSession *session) {
    uint32_t list = 0;
    size_t count = 0;
    RssnStatus status = rssn_value_call(session, "num_backend_available", NULL, 0, &list);
    if (status != RSSN_STATUS_OK) return status;

    status = rssn_value_len(session, list, &count);
    puts("Usable backend selectors");
    for (size_t index = 0; status == RSSN_STATUS_OK && index < count; ++index) {
        uint32_t name = 0;
        status = rssn_value_at(session, list, index, &name);
        if (status == RSSN_STATUS_OK) {
            status = print_text(session, name);
            rssn_value_free(session, name);
        }
    }
    rssn_value_free(session, list);
    return status;
}

int main(int argc, char **argv) {
    if (argc > 2) {
        fprintf(stderr, "Usage: %s [backend selector]\n", argv[0]);
        return 2;
    }

    struct RssnSession *session = rssn_session_new();
    if (session == NULL) {
        fprintf(stderr, "%s\n", rssn_last_error());
        return 1;
    }

    uint32_t argument = 0, response = 0;
    RssnStatus status = list_backends(session);
    if (status != RSSN_STATUS_OK) goto done;

    if (argc == 2) {
        status = rssn_value_new_text(session, argv[1], &argument);
        if (status != RSSN_STATUS_OK) goto done;
        status = rssn_value_call(session, "num_backend_select_default", &argument, 1, &response);
    } else {
        status = rssn_value_call(session, "num_backend_current_name", NULL, 0, &response);
    }
    if (status == RSSN_STATUS_OK) {
        puts("Effective backend");
        status = print_text(session, response);
    }

done:
    if (status != RSSN_STATUS_OK) {
        fprintf(stderr, "%s\n", rssn_last_error());
    }
    rssn_session_free(session);
    return status == RSSN_STATUS_OK ? 0 : 1;
}
