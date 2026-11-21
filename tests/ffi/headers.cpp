#include "rssn.hpp"

int main() {
    const auto* version = rssn::rssn_version();
    if (version == nullptr || *version == '\0') {
        return 1;
    }

    auto* session = rssn::rssn_session_new();
    if (session == nullptr) {
        return 2;
    }

    const auto term = rssn::rssn_term_sym(session, "temperature");
    rssn::rssn_session_free(session);
    return term == RSSN_NO_TERM ? 3 : 0;
}
