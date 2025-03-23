/**
 * calculator.h - 计算器库头文件
 */

#ifndef CALCULATOR_H
#define CALCULATOR_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * 基本算术操作
 */
double add(double a, double b);
double subtract(double a, double b);
double multiply(double a, double b);
double divide(double a, double b);

/**
 * 高级数学函数
 */
double power(double base, int exponent);
double square_root(double x);
int factorial(int n);
double absolute(double x);

/**
 * 三角函数
 */
double sine(double angle_in_degrees);
double cosine(double angle_in_degrees);
double tangent(double angle_in_degrees);

/**
 * 对数函数
 */
double natural_log(double x);
double log_base10(double x);

/**
 * 统计函数
 */
double calculate_mean(const double values[], int size);
double calculate_variance(const double values[], int size);
double calculate_std_deviation(const double values[], int size);

/**
 * 转换函数
 */
double celsius_to_fahrenheit(double celsius);
double fahrenheit_to_celsius(double fahrenheit);
double degrees_to_radians(double degrees);
double radians_to_degrees(double radians);

/**
 * 金融计算函数
 */
double calculate_simple_interest(double principal, double rate, double time);
double calculate_compound_interest(double principal, double rate, double time, int frequency);
double calculate_loan_payment(double principal, double rate, int term_months);

/**
 * 位操作函数
 */
unsigned int bit_count(unsigned int value);
unsigned int bit_reverse(unsigned int value);
int is_power_of_two(unsigned int value);

#ifdef __cplusplus
}
#endif

#endif /* CALCULATOR_H */ 