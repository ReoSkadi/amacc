/**
 * advanced_calculator.h - 高级计算器库头文件
 */

#ifndef ADVANCED_CALCULATOR_H
#define ADVANCED_CALCULATOR_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * 矩阵结构体定义
 */
typedef struct {
    double** data;
    size_t rows;
    size_t cols;
} Matrix;

/**
 * 复数结构体定义
 */
typedef struct {
    double real;
    double imag;
} Complex;

/**
 * 张量结构体定义
 */
typedef struct {
    double* data;
    size_t* dims;
    size_t ndims;
    size_t total_size;
} Tensor;

/**
 * 稀疏矩阵结构体定义
 */
typedef struct {
    double* values;       // 非零元素值
    size_t* row_indices;  // 行索引
    size_t* col_indices;  // 列索引
    size_t nnz;           // 非零元素数量
    size_t rows;          // 矩阵行数
    size_t cols;          // 矩阵列数
} SparseMatrix;

/**
 * 矩阵操作函数
 */
Matrix* matrix_create(size_t rows, size_t cols);
void matrix_free(Matrix* matrix);
Matrix* matrix_add(const Matrix* a, const Matrix* b);
Matrix* matrix_subtract(const Matrix* a, const Matrix* b);
Matrix* matrix_multiply(const Matrix* a, const Matrix* b);
Matrix* matrix_transpose(const Matrix* matrix);
double matrix_determinant(const Matrix* matrix);
Matrix* matrix_inverse(const Matrix* matrix);

/**
 * 新增矩阵高级操作函数
 */
Matrix* matrix_power(const Matrix* matrix, int power);
double matrix_trace(const Matrix* matrix);
double matrix_norm(const Matrix* matrix, const char* norm_type);
Matrix* matrix_eigenvalues(const Matrix* matrix);
Matrix* matrix_eigenvectors(const Matrix* matrix);
Matrix* matrix_svd(const Matrix* matrix, Matrix** U, Matrix** V);
Matrix* matrix_lu_decomposition(const Matrix* matrix, Matrix** L, Matrix** U, int* pivot);
Matrix* matrix_qr_decomposition(const Matrix* matrix, Matrix** Q, Matrix** R);
Matrix* matrix_cholesky(const Matrix* matrix);
Matrix* matrix_solve(const Matrix* A, const Matrix* b);
Matrix* matrix_pseudo_inverse(const Matrix* matrix);
Matrix* matrix_hadamard_product(const Matrix* a, const Matrix* b);

/**
 * 张量操作函数
 */
Tensor* tensor_create(const size_t* dims, size_t ndims);
void tensor_free(Tensor* tensor);
Tensor* tensor_add(const Tensor* a, const Tensor* b);
Tensor* tensor_subtract(const Tensor* a, const Tensor* b);
Tensor* tensor_element_multiply(const Tensor* a, const Tensor* b);
size_t tensor_get_index(const Tensor* tensor, const size_t* indices);
size_t tensor_get_total_size(const size_t* dims, size_t ndims);

/**
 * 稀疏矩阵操作函数
 */
SparseMatrix* sparse_matrix_create(size_t rows, size_t cols, size_t nnz);
void sparse_matrix_free(SparseMatrix* matrix);
SparseMatrix* sparse_matrix_from_dense(const Matrix* dense);
Matrix* sparse_matrix_to_dense(const SparseMatrix* sparse);
SparseMatrix* sparse_matrix_add(const SparseMatrix* a, const SparseMatrix* b);
SparseMatrix* sparse_matrix_multiply(const SparseMatrix* a, const SparseMatrix* b);
SparseMatrix* sparse_matrix_transpose(const SparseMatrix* matrix);

/**
 * 复数运算函数
 */
Complex complex_add(Complex a, Complex b);
Complex complex_subtract(Complex a, Complex b);
Complex complex_multiply(Complex a, Complex b);
Complex complex_divide(Complex a, Complex b);
double complex_modulus(Complex z);
double complex_argument(Complex z);
Complex complex_conjugate(Complex z);
Complex complex_exp(Complex z);
Complex complex_log(Complex z);
Complex complex_power(Complex base, Complex exponent);

/**
 * 多项式结构定义
 */
typedef struct {
    double* coefficients;
    int degree;
} Polynomial;

/**
 * 多项式操作函数
 */
Polynomial* polynomial_create(double* coefficients, int degree);
void polynomial_free(Polynomial* poly);
double polynomial_evaluate(const Polynomial* poly, double x);
Polynomial* polynomial_add(const Polynomial* a, const Polynomial* b);
Polynomial* polynomial_multiply(const Polynomial* a, const Polynomial* b);
Polynomial* polynomial_derivative(const Polynomial* poly);
Polynomial* polynomial_integral(const Polynomial* poly);
double* polynomial_find_roots(const Polynomial* poly, int* num_roots);

/**
 * 统计分析函数
 */
double* linear_regression(const double x[], const double y[], int size, double* slope, double* intercept);
double correlation_coefficient(const double x[], const double y[], int size);
double covariance(const double x[], const double y[], int size);
double* multiple_regression(double** x, const double y[], int n_samples, int n_features, double* coefficients);

/**
 * 数值分析函数
 */
double numerical_integration(double (*f)(double), double a, double b, int n);
double root_finding_bisection(double (*f)(double), double a, double b, double tol, int max_iter);
double root_finding_newton(double (*f)(double), double (*df)(double), double x0, double tol, int max_iter);

/**
 * 概率统计函数
 */
double normal_pdf(double x, double mean, double stddev);
double normal_cdf(double x, double mean, double stddev);
double poisson_pmf(int k, double lambda);
double binomial_pmf(int k, int n, double p);

/**
 * 金融高级计算函数
 */
double black_scholes_call(double S, double K, double r, double sigma, double T);
double black_scholes_put(double S, double K, double r, double sigma, double T);
double option_implied_volatility(double option_price, double S, double K, double r, double T, int is_call);

#ifdef __cplusplus
}
#endif

#endif /* ADVANCED_CALCULATOR_H */ 