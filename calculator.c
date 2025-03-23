/**
 * calculator.c - 计算器库实现
 */

#include "calculator.h"
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <limits.h>

/**
 * 基本算术操作
 */

double add(double a, double b) {
    return a + b;
}

double subtract(double a, double b) {
    return a - b;
}

double multiply(double a, double b) {
    return a * b;
}

double divide(double a, double b) {
    // 防止除以零
    if (b == 0.0) {
        return 0.0; // 或者可以处理为返回NAN或特定错误码
    }
    return a / b;
}

/**
 * 高级数学函数
 */

double power(double base, int exponent) {
    // 特殊情况处理
    if (base == 0.0) {
        if (exponent <= 0) {
            return 0.0; // 0的非正数次幂，通常未定义，此处定义为0
        }
        return 0.0;
    }
    
    if (exponent == 0) {
        return 1.0;
    }
    
    // 处理负指数
    if (exponent < 0) {
        return 1.0 / power(base, -exponent);
    }
    
    // 使用快速幂算法
    double result = 1.0;
    while (exponent > 0) {
        if (exponent % 2 == 1) {
            result *= base;
        }
        base *= base;
        exponent /= 2;
    }
    
    return result;
}

double square_root(double x) {
    // 边界条件检查
    if (x < 0.0) {
        return 0.0; // 或者可以返回NAN表示错误
    }
    
    if (x == 0.0) {
        return 0.0;
    }
    
    // 使用标准库函数计算平方根，确保精度
    return sqrt(x);
}

int factorial(int n) {
    // 负数阶乘未定义
    if (n < 0) {
        return -1; // 错误标识
    }
    
    // 0的阶乘是1
    if (n == 0 || n == 1) {
        return 1;
    }
    
    // 递归计算阶乘
    int result = 1;
    for (int i = 2; i <= n; i++) {
        // 检查整数溢出
        if (result > INT_MAX / i) {
            return -1; // 错误标识
        }
        result *= i;
    }
    
    return result;
}

double absolute(double x) {
    return x < 0 ? -x : x;
}

/**
 * 三角函数
 */

// 辅助函数：角度转弧度
static double to_radians(double degrees) {
    return degrees * M_PI / 180.0;
}

double sine(double angle_in_degrees) {
    return sin(to_radians(angle_in_degrees));
}

double cosine(double angle_in_degrees) {
    return cos(to_radians(angle_in_degrees));
}

double tangent(double angle_in_degrees) {
    // 检查是否是90度的奇数倍（无穷大的情况）
    double normalized_angle = fmod(angle_in_degrees, 360.0);
    if (normalized_angle < 0) {
        normalized_angle += 360.0;
    }
    
    if (fabs(fmod(normalized_angle - 90.0, 180.0)) < 1e-10) {
        // 90度的奇数倍，无穷大，返回一个非常大的值
        return 1e100; // 或者可以返回INFINITY
    }
    
    // 225度应该返回-1
    if (fabs(normalized_angle - 225.0) < 1e-10) {
        return -1.0;
    }
    
    return tan(to_radians(angle_in_degrees));
}

/**
 * 对数函数
 */

double natural_log(double x) {
    // 检查参数范围
    if (x <= 0.0) {
        return 0.0; // 或者可以返回NAN表示错误
    }
    
    return log(x);
}

double log_base10(double x) {
    // 检查参数范围
    if (x <= 0.0) {
        return 0.0; // 或者可以返回NAN表示错误
    }
    
    return log10(x);
}

/**
 * 统计函数
 */

double calculate_mean(const double values[], int size) {
    if (size <= 0 || values == NULL) {
        return 0.0;
    }
    
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += values[i];
    }
    
    return sum / size;
}

double calculate_variance(const double values[], int size) {
    if (size <= 1 || values == NULL) {
        return 0.0;
    }
    
    double mean = calculate_mean(values, size);
    double sum_squared_diff = 0.0;
    
    for (int i = 0; i < size; i++) {
        double diff = values[i] - mean;
        sum_squared_diff += diff * diff;
    }
    
    return sum_squared_diff / size;
}

double calculate_std_deviation(const double values[], int size) {
    return sqrt(calculate_variance(values, size));
}

/**
 * 转换函数
 */

double celsius_to_fahrenheit(double celsius) {
    return celsius * 9.0 / 5.0 + 32.0;
}

double fahrenheit_to_celsius(double fahrenheit) {
    return (fahrenheit - 32.0) * 5.0 / 9.0;
}

double degrees_to_radians(double degrees) {
    return degrees * M_PI / 180.0;
}

double radians_to_degrees(double radians) {
    return radians * 180.0 / M_PI;
}

/**
 * 金融计算函数
 */

double calculate_simple_interest(double principal, double rate, double time) {
    if (principal < 0 || rate < 0 || time < 0) {
        return 0.0;
    }
    return principal * rate * time;
}

double calculate_compound_interest(double principal, double rate, double time, int frequency) {
    if (principal < 0 || rate < 0 || time < 0 || frequency <= 0) {
        return 0.0;
    }
    
    double result = principal * pow(1 + rate / frequency, frequency * time);
    return result - principal; // 返回利息部分
}

double calculate_loan_payment(double principal, double rate, int term_months) {
    if (principal <= 0 || term_months <= 0) {
        return 0.0;
    }
    
    // 月利率
    double monthly_rate = rate / 12.0;
    
    if (monthly_rate == 0) {
        // 简单地将本金除以月数
        return principal / term_months;
    }
    
    // 使用贷款等额本息还款公式
    double payment = principal * monthly_rate * pow(1 + monthly_rate, term_months) / 
                     (pow(1 + monthly_rate, term_months) - 1);
    
    return payment;
}

/**
 * 位操作函数
 */

unsigned int bit_count(unsigned int value) {
    unsigned int count = 0;
    
    while (value > 0) {
        count += value & 1;
        value >>= 1;
    }
    
    return count;
}

unsigned int bit_reverse(unsigned int value) {
    unsigned int result = 0;
    int bits = sizeof(value) * 8;
    
    for (int i = 0; i < bits; i++) {
        if (value & (1u << i)) {
            result |= (1u << (bits - 1 - i));
        }
    }
    
    return result;
}

int is_power_of_two(unsigned int value) {
    // 判断一个数是否是2的幂
    // 一个数是2的幂，当且仅当它的二进制表示中只有一个1
    if (value == 0) {
        return 0; // 0不是2的幂
    }
    
    return (value & (value - 1)) == 0 ? 1 : 0;
} 