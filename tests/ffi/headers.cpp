#include "rssn.hpp"

#include <array>
#include <cstring>

int main() {
  const auto *version = rssn::rssn_version();
  if (version == nullptr || *version == '\0') {
    return 1;
  }

  auto *session = rssn::rssn_session_new();
  if (session == nullptr) {
    return 2;
  }

  const auto term = rssn::rssn_term_sym(session, "temperature");
  if (term == RSSN_NO_TERM) {
    rssn::rssn_session_free(session);
    return 3;
  }

  char *rules = nullptr;
  char *name = nullptr;
  std::array<std::uint8_t, 16> digest{};
  auto status = rssn::rssn_op_key_at("num_number_theory_gcd_u128", 0, &rules,
                                     &name, digest.data());
  std::uint32_t left = 0;
  std::uint32_t right = 0;
  std::uint32_t result = 0;
  using rssn::RssnStatus;
  if (status == RssnStatus::RSSN_STATUS_OK)
    status = rssn::rssn_value_new_int(session, 12, &left);
  if (status == RssnStatus::RSSN_STATUS_OK)
    status = rssn::rssn_value_new_int(session, 18, &right);
  const std::uint32_t cells[] = {left, right};
  if (status == RssnStatus::RSSN_STATUS_OK) {
    status = rssn::rssn_value_call_key(session, rules, name, digest.data(),
                                       cells, 2, &result);
  }
  char *text = nullptr;
  if (status == RssnStatus::RSSN_STATUS_OK)
    status = rssn::rssn_value_to_text(session, result, &text);
  const bool valid = status == RssnStatus::RSSN_STATUS_OK && text != nullptr &&
                     std::strcmp(text, "6") == 0;

  rssn::rssn_string_free(text);
  rssn::rssn_string_free(name);
  rssn::rssn_string_free(rules);
  rssn::rssn_session_free(session);
  return valid ? 0 : 4;
}
