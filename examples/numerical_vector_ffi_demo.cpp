#include <iostream>
#include <vector>
#include <iomanip>

extern "C" {
    // Basic handle functions
    void* rssn_num_vec_create(const double* data, size_t len);
    void rssn_num_vec_free(void* v);
    size_t rssn_num_vec_len(const void* v);
    const double* rssn_num_vec_data(const void* v);

    // Operations
    void* rssn_num_vec_add(const void* v1, const void* v2);
    void* rssn_num_vec_sub(const void* v1, const void* v2);
    void* rssn_num_vec_scalar_mul(const void* v, double s);
    int rssn_num_vec_dot_product(const void* v1, const void* v2, double* result);
    int rssn_num_vec_norm(const void* v, double* result);
}

void print_vec(void* v, const std::string& label) {
    if (!v) {
        std::cout << label << ": NULL" << std::endl;
        return;
    }
    size_t len = rssn_num_vec_len(v);
    const double* data = rssn_num_vec_data(v);
    std::cout << label << ": [";
    for (size_t i = 0; i < len; ++i) {
        std::cout << data[i] << (i == len - 1 ? "" : ", ");
    }
    std::cout << "]" << std::endl;
}

int main() {
    std::cout << "--- C++ FFI Demo: Numerical Vectors ---" << std::endl;

    double d1[] = {1.0, 2.0, 3.0};
    double d2[] = {4.0, 5.0, 6.0};

    void* v1 = rssn_num_vec_create(d1, 3);
    void* v2 = rssn_num_vec_create(d2, 3);

    print_vec(v1, "v1");
    print_vec(v2, "v2");

    // Add
    void* v_sum = rssn_num_vec_add(v1, v2);
    print_vec(v_sum, "v1 + v2");

    // Dot Product
    double dot;
    if (rssn_num_vec_dot_product(v1, v2, &dot) == 0) {
        std::cout << "v1 . v2 = " << dot << std::endl;
    }

    // Norm
    double n1;
    if (rssn_num_vec_norm(v1, &n1) == 0) {
        std::cout << "||v1|| = " << n1 << std::endl;
    }

    // Cleanup
    rssn_num_vec_free(v1);
    rssn_num_vec_free(v2);
    rssn_num_vec_free(v_sum);

    return 0;
}
